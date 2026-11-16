//! Whether two molecules are the same, however they were written (#58).
use crate::descriptors::fingerprints::ecfp::fnv1a;
use crate::model::atom::Atom;
use crate::model::bond::{Bond, BondOrder};
use crate::model::matcher::{
    find_embedding, find_embedding_unpruned, mix64, mix_word, pack_bond, plan_query_order,
    refine_hash_signature, CsrGraph, EmbedOutcome, EmbeddingSpec, SymmetryGroup,
};
use crate::model::molecule::Molecule;
use crate::model::stereo::{Configuration, CoordinationGeometry, Ligand};
use petgraph::stable_graph::{EdgeReference, NodeIndex, StableGraph};
use petgraph::visit::{EdgeRef, IntoEdgeReferences, NodeIndexable};
use petgraph::Undirected;
use std::sync::{Arc, OnceLock};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExactMatchError {
    /// Square-planar, trigonal-bipyramidal or octahedral stereo, which isn't compared yet.
    UnsupportedStereo,
    /// More than the search limit of atom mappings to try.
    TooSymmetric,
}

/// Core atom attributes packed without per-atom heap allocations; field order matches the
/// canonical ordering used when folding and hashing twin end-atoms.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct AtomBase {
    atomic_number: u8,
    isotope: u16,
    charge: i8,
    hydrogens: u32,
    aromatic: bool,
}

impl AtomBase {
    #[inline]
    fn packed(self) -> u64 {
        u64::from(self.atomic_number)
            | (u64::from(self.aromatic) << 8)
            | (u64::from(self.charge.cast_unsigned()) << 9)
            | (u64::from(self.isotope) << 17)
            | (u64::from(self.hydrogens) << 33)
    }
}

/// What became of a source atom.
#[derive(Clone, Copy)]
enum Fate {
    Kept(usize),
    /// Folded into its neighbour's hydrogen count.
    Hydrogen,
    /// One of two or more alike end atoms on one neighbour, folded into its key.
    Twin,
}

/// A stereo ligand, by its atom in the folded graph.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Site {
    Atom(usize),
    /// A hydrogen, by the atom it is on, so an allene's two terminal hydrogens differ.
    Hydrogen(usize),
    LonePair,
    Twin,
    Vacancy,
}

impl Site {
    fn mapped(self, mapping: &[usize]) -> Site {
        match self {
            Site::Atom(atom) => Site::Atom(mapping[atom]),
            Site::Hydrogen(atom) => Site::Hydrogen(mapping[atom]),
            other => other,
        }
    }
}

enum Stereo {
    /// A tetrahedral centre, or an allene's middle atom with its four terminal ligands.
    Parity {
        center: usize,
        ligands: [Site; 4],
        clockwise: bool,
    },
    Coordination {
        center: usize,
        ligands: [Site; 6],
        len: usize,
        configuration: Configuration,
    },
    DoubleBond {
        ends: [usize; 2],
        references: [Site; 2],
        same_side: bool,
    },
}

impl Stereo {
    fn mapped(&self, mapping: &[usize]) -> Stereo {
        match *self {
            Stereo::Parity {
                center,
                ligands,
                clockwise,
            } => Stereo::Parity {
                center: mapping[center],
                ligands: ligands.map(|site| site.mapped(mapping)),
                clockwise,
            },
            Stereo::Coordination {
                center,
                ligands,
                len,
                configuration,
            } => {
                let mut mapped = [Site::Vacancy; 6];
                for (slot, site) in mapped[..len].iter_mut().zip(&ligands[..len]) {
                    *slot = site.mapped(mapping);
                }
                Stereo::Coordination {
                    center: mapping[center],
                    ligands: mapped,
                    len,
                    configuration,
                }
            }
            Stereo::DoubleBond {
                ends,
                references,
                same_side,
            } => Stereo::DoubleBond {
                ends: ends.map(|end| mapping[end]),
                references: references.map(|site| site.mapped(mapping)),
                same_side,
            },
        }
    }

    /// The atoms it belongs to, sorted, so a mapped element can be looked up.
    fn anchor(&self, mapping: impl Fn(usize) -> usize) -> [usize; 2] {
        match self {
            Stereo::Parity { center, .. } | Stereo::Coordination { center, .. } => {
                [mapping(*center); 2]
            }
            Stereo::DoubleBond { ends, .. } => {
                let [a, b] = ends.map(mapping);
                [a.min(b), a.max(b)]
            }
        }
    }

    /// Whether this element, carried over by `mapping`, has `target`'s configuration.
    fn agrees(&self, mapping: &[usize], target: &Stereo) -> bool {
        match (self, target) {
            (
                Stereo::Parity {
                    ligands, clockwise, ..
                },
                Stereo::Parity {
                    ligands: target_ligands,
                    clockwise: target_clockwise,
                    ..
                },
            ) => {
                let mut positions = [0; 4];
                for (position, site) in positions.iter_mut().zip(ligands) {
                    let site = site.mapped(mapping);
                    let Some(found) = target_ligands.iter().position(|&other| other == site) else {
                        return false;
                    };
                    *position = found;
                }
                let inversions = (0..4)
                    .flat_map(|i| (i + 1..4).map(move |j| (i, j)))
                    .filter(|&(i, j)| positions[i] > positions[j])
                    .count();
                (inversions % 2 == 0) == (clockwise == target_clockwise)
            }
            (
                Stereo::Coordination {
                    ligands,
                    len,
                    configuration,
                    ..
                },
                Stereo::Coordination {
                    ligands: target_ligands,
                    len: target_len,
                    configuration: target_configuration,
                    ..
                },
            ) => {
                let mut mapped = [Site::Vacancy; 6];
                for (slot, site) in mapped[..*len].iter_mut().zip(&ligands[..*len]) {
                    *slot = site.mapped(mapping);
                }
                let lhs = configuration.coordination_geometry(&mapped[..*len]);
                let rhs =
                    target_configuration.coordination_geometry(&target_ligands[..*target_len]);
                lhs.is_some() && lhs == rhs
            }
            (
                Stereo::DoubleBond {
                    ends,
                    references,
                    same_side,
                },
                Stereo::DoubleBond {
                    ends: target_ends,
                    references: target_references,
                    same_side: target_same_side,
                },
            ) => {
                let target_references = if *target_ends == ends.map(|end| mapping[end]) {
                    *target_references
                } else {
                    [target_references[1], target_references[0]]
                };
                // each end's other substituent is on the opposite side
                let flips = (0..2)
                    .filter(|&i| references[i].mapped(mapping) != target_references[i])
                    .count();
                (*same_side == (flips % 2 == 0)) == *target_same_side
            }
            _ => false,
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum RankedSite {
    Atom(usize),
    Hydrogen(usize),
    LonePair,
    Twin,
    Vacancy,
}

fn site_code(site: RankedSite) -> u64 {
    match site {
        RankedSite::Atom(r) => (r as u64) << 4,
        RankedSite::Hydrogen(r) => ((r as u64) << 4) | 1,
        RankedSite::LonePair => 2,
        RankedSite::Twin => 3,
        RankedSite::Vacancy => 4,
    }
}

fn encode_coordination_leaf(
    center_rank: usize,
    ligands: &[Site],
    configuration: Configuration,
    map_site: impl Fn(Site) -> RankedSite,
) -> [u64; 8] {
    let mut sites = [RankedSite::Vacancy; 6];
    for (slot, &s) in sites[..ligands.len()].iter_mut().zip(ligands) {
        *slot = map_site(s);
    }
    let geo = configuration.coordination_geometry(&sites[..ligands.len()]);
    let rc = center_rank as u64;
    match geo {
        Some(CoordinationGeometry::SquarePlanar(pairs)) => [
            2,
            rc,
            site_code(pairs[0][0]),
            site_code(pairs[0][1]),
            site_code(pairs[1][0]),
            site_code(pairs[1][1]),
            0,
            0,
        ],
        Some(CoordinationGeometry::TrigonalBipyramidal((a, eq, b))) => [
            3,
            rc,
            site_code(a),
            site_code(eq[0]),
            site_code(eq[1]),
            site_code(eq[2]),
            site_code(b),
            0,
        ],
        Some(CoordinationGeometry::Octahedral((a, eq, b))) => [
            4,
            rc,
            site_code(a),
            site_code(eq[0]),
            site_code(eq[1]),
            site_code(eq[2]),
            site_code(eq[3]),
            site_code(b),
        ],
        None => [5, rc, 0, 0, 0, 0, 0, 0],
    }
}

impl Stereo {
    fn anchor_atoms(&self) -> (usize, Option<usize>) {
        match *self {
            Stereo::Parity { center, .. } | Stereo::Coordination { center, .. } => (center, None),
            Stereo::DoubleBond { ends, .. } => (ends[0], Some(ends[1])),
        }
    }

    fn collect_vertices(&self, out: &mut Vec<usize>) {
        out.clear();
        let push_site = |out: &mut Vec<usize>, site: Site| match site {
            Site::Atom(a) | Site::Hydrogen(a) => out.push(a),
            Site::LonePair | Site::Twin | Site::Vacancy => {}
        };
        match *self {
            Stereo::Parity {
                center, ligands, ..
            } => {
                out.push(center);
                for s in ligands {
                    push_site(out, s);
                }
            }
            Stereo::Coordination {
                center,
                ligands,
                len,
                ..
            } => {
                out.push(center);
                for &s in &ligands[..len] {
                    push_site(out, s);
                }
            }
            Stereo::DoubleBond {
                ends, references, ..
            } => {
                out.extend_from_slice(&ends);
                for s in references {
                    push_site(out, s);
                }
            }
        }
        out.sort_unstable();
        out.dedup();
    }

    fn encode_for_leaf(&self, rank_of: &[usize], csr: &CsrGraph, keys: &[AtomBase]) -> [u64; 8] {
        let map_site = |site: Site| match site {
            Site::Atom(a) => RankedSite::Atom(rank_of[a]),
            Site::Hydrogen(a) => RankedSite::Hydrogen(rank_of[a]),
            Site::LonePair => RankedSite::LonePair,
            Site::Twin => RankedSite::Twin,
            Site::Vacancy => RankedSite::Vacancy,
        };
        let site_u64 = site_code;
        match *self {
            Stereo::Parity {
                center,
                ligands,
                clockwise,
            } => {
                let mut sites = ligands.map(map_site);
                let mut odd = false;
                for i in 0..4 {
                    for j in i + 1..4 {
                        if sites[i] > sites[j] {
                            sites.swap(i, j);
                            odd = !odd;
                        }
                    }
                }
                [
                    1,
                    rank_of[center] as u64,
                    u64::from(clockwise ^ odd),
                    site_u64(sites[0]),
                    site_u64(sites[1]),
                    site_u64(sites[2]),
                    site_u64(sites[3]),
                    0,
                ]
            }
            Stereo::Coordination {
                center,
                ligands,
                len,
                configuration,
            } => {
                encode_coordination_leaf(rank_of[center], &ligands[..len], configuration, map_site)
            }
            Stereo::DoubleBond {
                ends,
                references,
                same_side,
            } => {
                let mut flipped = false;
                for side in [0, 1] {
                    let end = ends[side];
                    let partner = ends[1 - side];
                    let ref_ranked = map_site(references[side]);
                    let mut subs: Vec<RankedSite> = csr
                        .neighbors(end)
                        .iter()
                        .filter(|&&(_, v)| v != partner)
                        .map(|&(_, v)| RankedSite::Atom(rank_of[v]))
                        .collect();
                    if keys[end].hydrogens == 1 {
                        subs.push(RankedSite::Hydrogen(rank_of[end]));
                    }
                    if let Some(&min_sub) = subs.iter().min() {
                        if ref_ranked != min_sub {
                            flipped = !flipped;
                        }
                    }
                }
                let r0 = rank_of[ends[0]] as u64;
                let r1 = rank_of[ends[1]] as u64;
                [
                    6,
                    r0.min(r1),
                    r0.max(r1),
                    u64::from(same_side ^ flipped),
                    0,
                    0,
                    0,
                    0,
                ]
            }
        }
    }
}

/// A molecule as compared. Hydrogens written as atoms are folded into their neighbour's
/// count, so `[H]C([H])([H])[H]` is the same as `C`. Alike end atoms on one neighbour,
/// such as the methyls of `N(C)C` or the oxygens of `S(=O)(=O)`, are folded into its key:
/// swapping them changes nothing, but unpruned VF2 would try both orders of every such pair.
pub(crate) struct Folded {
    csr: CsrGraph,
    keys: Vec<AtomBase>,
    twin_spans: Vec<(u32, u32)>,
    twin_pool: Vec<(AtomBase, u8)>,
    refined_colors: Vec<usize>,
    refined_hashes: Vec<u64>,
    color_hist: Vec<(u64, usize)>,
    dative_count: usize,
    stereo: Vec<Stereo>,
    /// Sorted by anchor `[usize; 2]` for fast lookup without `HashMap` allocation.
    anchors: Vec<([usize; 2], usize)>,
    fake: OnceLock<Vec<bool>>,
    real_sym: OnceLock<SymmetryGroup>,
}

struct FoldedRaw {
    n: usize,
    offsets: Vec<usize>,
    edges: Vec<(u8, usize)>,
    keys: Vec<AtomBase>,
    twin_spans: Vec<(u32, u32)>,
    twin_pool: Vec<(AtomBase, u8)>,
    dative_bonds: Vec<usize>,
    fates: Vec<Fate>,
}

impl FoldedRaw {
    #[inline]
    fn twins_of(&self, u: usize) -> &[(AtomBase, u8)] {
        if self.twin_pool.is_empty() {
            return &[];
        }
        let (start, len) = self.twin_spans[u];
        &self.twin_pool[start as usize..(start + len) as usize]
    }
}

#[inline]
fn is_fake_at(fake: &[bool], idx: usize) -> bool {
    fake.get(idx).copied().unwrap_or(false)
}

#[inline]
fn edge_other(edge: EdgeReference<'_, Bond>, node: NodeIndex) -> NodeIndex {
    if edge.source() == node {
        edge.target()
    } else {
        edge.source()
    }
}

fn is_plain_hydrogen(source: &StableGraph<Atom, Bond, Undirected>, node: NodeIndex) -> bool {
    let atom = &source[node];
    let mut edges = source.edges(node);
    atom.element.atomic_number == 1
        && atom.isotope == 0
        && atom.charge == 0
        && atom.hs + atom.implicit_hs == 0
        && matches!((edges.next(), edges.next()), (Some(edge), None)
            if edge.weight().order == BondOrder::Single
                && source[edge_other(edge, node)].element.atomic_number != 1)
}

struct TwinCollapse {
    kept_deg: Vec<usize>,
    is_twin: Vec<bool>,
    raw_twins: Vec<(usize, AtomBase, u8)>,
}

fn collapse_twins(
    source: &StableGraph<Atom, Bond, Undirected>,
    kept_idx: &[usize],
    kept_nodes: &[NodeIndex],
    kept_bases: &[AtomBase],
) -> TwinCollapse {
    let n_kept = kept_nodes.len();
    let mut kept_deg = vec![0usize; n_kept];
    let mut end_cands: Vec<(usize, u8, AtomBase, usize)> = Vec::new();
    for (u, &node) in kept_nodes.iter().enumerate() {
        let mut deg = 0usize;
        let mut last_parent = 0usize;
        let mut last_order = BondOrder::Single;
        for edge in source.edges(node) {
            let other = kept_idx[edge_other(edge, node).index()];
            if other != usize::MAX {
                deg += 1;
                last_parent = other;
                last_order = edge.weight().order;
            }
        }
        kept_deg[u] = deg;
        if deg == 1 && last_order != BondOrder::Dative {
            end_cands.push((last_parent, pack_bond(last_order, false), kept_bases[u], u));
        }
    }

    let mut is_twin = vec![false; n_kept];
    let mut raw_twins: Vec<(usize, AtomBase, u8)> = Vec::new();
    if end_cands.len() > 1 {
        end_cands.sort_unstable();
        let mut i = 0;
        while i < end_cands.len() {
            let (parent, bond_code, base, _) = end_cands[i];
            let mut j = i + 1;
            while j < end_cands.len()
                && end_cands[j].0 == parent
                && end_cands[j].1 == bond_code
                && end_cands[j].2 == base
            {
                j += 1;
            }
            if j - i > 1 {
                for &(_, _, _, u) in &end_cands[i..j] {
                    is_twin[u] = true;
                    kept_deg[u] = 0;
                    kept_deg[parent] -= 1;
                    raw_twins.push((parent, base, bond_code));
                }
            }
            i = j;
        }
    }
    TwinCollapse {
        kept_deg,
        is_twin,
        raw_twins,
    }
}

type TwinPool = (Vec<(u32, u32)>, Vec<(AtomBase, u8)>);

fn build_twin_pool(
    n: usize,
    folded_of_kept: &[usize],
    mut raw_twins: Vec<(usize, AtomBase, u8)>,
) -> TwinPool {
    if raw_twins.is_empty() {
        return (Vec::new(), Vec::new());
    }
    for entry in &mut raw_twins {
        entry.0 = folded_of_kept[entry.0];
    }
    raw_twins.sort_unstable();
    let mut spans = vec![(0u32, 0u32); n];
    let mut pool = Vec::with_capacity(raw_twins.len());
    let mut i = 0;
    while i < raw_twins.len() {
        let owner = raw_twins[i].0;
        let start = pool.len();
        while i < raw_twins.len() && raw_twins[i].0 == owner {
            pool.push((raw_twins[i].1, raw_twins[i].2));
            i += 1;
        }
        let len = pool.len() - start;
        spans[owner] = (
            u32::try_from(start).unwrap_or(u32::MAX),
            u32::try_from(len).unwrap_or(u32::MAX),
        );
    }
    (spans, pool)
}

fn fold_raw(molecule: &Molecule) -> FoldedRaw {
    let source = molecule.graph();
    let nb = source.node_bound();
    let mut kept_idx = vec![usize::MAX; nb];
    let mut kept_nodes = Vec::with_capacity(source.node_count());
    let mut kept_bases = Vec::with_capacity(source.node_count());

    for node in source.node_indices() {
        if is_plain_hydrogen(source, node) {
            continue;
        }
        let atom = &source[node];
        let folded_h = source
            .neighbors(node)
            .filter(|&nbr| is_plain_hydrogen(source, nbr))
            .count();
        let hydrogens = u32::from(atom.hs)
            + u32::from(atom.implicit_hs)
            + u32::try_from(folded_h).unwrap_or(u32::MAX);
        kept_idx[node.index()] = kept_nodes.len();
        kept_nodes.push(node);
        kept_bases.push(AtomBase {
            atomic_number: atom.element.atomic_number,
            isotope: atom.isotope,
            charge: atom.charge,
            hydrogens,
            aromatic: atom.aromatic,
        });
    }

    let n_kept = kept_nodes.len();
    let collapse = collapse_twins(source, &kept_idx, &kept_nodes, &kept_bases);

    let mut fates = vec![Fate::Hydrogen; nb];
    let mut folded_of_kept = vec![usize::MAX; n_kept];
    let mut keys = Vec::with_capacity(n_kept);
    let mut offsets = Vec::with_capacity(n_kept + 1);
    offsets.push(0);
    for u in 0..n_kept {
        if collapse.is_twin[u] {
            fates[kept_nodes[u].index()] = Fate::Twin;
            continue;
        }
        let f = keys.len();
        folded_of_kept[u] = f;
        fates[kept_nodes[u].index()] = Fate::Kept(f);
        keys.push(kept_bases[u]);
        offsets.push(offsets[f] + collapse.kept_deg[u]);
    }
    let n = keys.len();
    let (twin_spans, twin_pool) = build_twin_pool(n, &folded_of_kept, collapse.raw_twins);

    let mut cursor = offsets[..n].to_vec();
    let mut edges = vec![(0u8, 0usize); offsets[n]];
    let mut dative_bonds = Vec::new();
    for edge in source.edge_references() {
        let (Fate::Kept(a), Fate::Kept(b)) =
            (fates[edge.source().index()], fates[edge.target().index()])
        else {
            continue;
        };
        let weight = edge.weight();
        if weight.donor.is_some() {
            dative_bonds.push(a);
        }
        let a_donor = weight.donor == Some(edge.source());
        let b_donor = weight.donor == Some(edge.target());
        edges[cursor[a]] = (pack_bond(weight.order, a_donor), b);
        cursor[a] += 1;
        edges[cursor[b]] = (pack_bond(weight.order, b_donor), a);
        cursor[b] += 1;
    }

    FoldedRaw {
        n,
        offsets,
        edges,
        keys,
        twin_spans,
        twin_pool,
        dative_bonds,
        fates,
    }
}

impl Folded {
    fn from_raw(raw: FoldedRaw, dative_count: usize, stereo: Vec<Stereo>) -> Folded {
        let n = raw.n;
        let mut base_hashes = Vec::with_capacity(n);
        if raw.twin_pool.is_empty() {
            for &base in &raw.keys {
                base_hashes.push(fast_atom_code(base, &[]));
            }
        } else {
            for (u, &base) in raw.keys.iter().enumerate() {
                base_hashes.push(fast_atom_code(base, raw.twins_of(u)));
            }
        }
        let csr = CsrGraph {
            n,
            offsets: raw.offsets,
            edges: raw.edges,
            base_hashes,
        };
        let (refined_colors, refined_hashes, color_hist) = refine_hash_signature(&csr);
        let mut anchors: Vec<([usize; 2], usize)> = stereo
            .iter()
            .enumerate()
            .map(|(i, element)| (element.anchor(|atom| atom), i))
            .collect();
        if anchors.len() > 16 {
            anchors.sort_unstable_by_key(|&(anchor, _)| anchor);
        }
        Folded {
            csr,
            keys: raw.keys,
            twin_spans: raw.twin_spans,
            twin_pool: raw.twin_pool,
            refined_colors,
            refined_hashes,
            color_hist,
            dative_count,
            stereo,
            anchors,
            fake: OnceLock::new(),
            real_sym: OnceLock::new(),
        }
    }

    #[inline]
    fn twins_of(&self, u: usize) -> &[(AtomBase, u8)] {
        if self.twin_pool.is_empty() {
            return &[];
        }
        let (start, len) = self.twin_spans[u];
        &self.twin_pool[start as usize..(start + len) as usize]
    }

    #[inline]
    fn atom_equal(&self, q: usize, target: &Folded, t: usize) -> bool {
        if self.keys[q] != target.keys[t] {
            return false;
        }
        if self.twin_pool.is_empty() && target.twin_pool.is_empty() {
            return true;
        }
        self.twins_of(q) == target.twins_of(t)
    }

    #[inline]
    fn anchor_get(&self, anchor: [usize; 2]) -> Option<usize> {
        if self.anchors.len() <= 16 {
            self.anchors
                .iter()
                .find_map(|&(key, idx)| (key == anchor).then_some(idx))
        } else {
            self.anchors
                .binary_search_by_key(&anchor, |&(key, _)| key)
                .ok()
                .map(|pos| self.anchors[pos].1)
        }
    }

    /// One `Folded` per connected component, so k alike salt or water parts cost k
    /// matches instead of k! mappings. No stereo element spans two components.
    pub(crate) fn components(molecule: &Molecule) -> Result<Vec<Folded>, ExactMatchError> {
        let raw = fold_raw(molecule);
        if raw.n == 0 {
            return Ok(Vec::new());
        }
        let stereo = Folded::read_stereo(molecule, &raw)?;
        let n = raw.n;

        // Connected component scan directly on the flat CSR; fast-paths single-component
        // molecules without allocating per-component graphs or remapping vertices.
        let mut component = vec![usize::MAX; n];
        let mut queue = Vec::with_capacity(n);
        component[0] = 0;
        queue.push(0usize);
        let mut head = 0usize;
        while head < queue.len() {
            let u = queue[head];
            head += 1;
            for &(_, v) in &raw.edges[raw.offsets[u]..raw.offsets[u + 1]] {
                if component[v] == usize::MAX {
                    component[v] = 0;
                    queue.push(v);
                }
            }
        }
        if queue.len() == n {
            let dative_count = raw.dative_bonds.len();
            return Ok(vec![Folded::from_raw(raw, dative_count, stereo)]);
        }

        Ok(Self::split_components(
            &raw,
            &mut component,
            &mut queue,
            stereo,
        ))
    }

    fn split_components(
        raw: &FoldedRaw,
        component: &mut [usize],
        queue: &mut Vec<usize>,
        stereo: Vec<Stereo>,
    ) -> Vec<Folded> {
        let n = raw.n;
        let mut comp_verts: Vec<Vec<usize>> = vec![Vec::new()];
        let mut local = vec![0usize; n];
        for start in 0..n {
            if component[start] != usize::MAX {
                continue;
            }
            let comp_id = comp_verts.len();
            comp_verts.push(Vec::new());
            queue.clear();
            component[start] = comp_id;
            queue.push(start);
            let mut qh = 0usize;
            while qh < queue.len() {
                let u = queue[qh];
                qh += 1;
                for &(_, v) in &raw.edges[raw.offsets[u]..raw.offsets[u + 1]] {
                    if component[v] == usize::MAX {
                        component[v] = comp_id;
                        queue.push(v);
                    }
                }
            }
        }
        for (u, &c) in component.iter().enumerate() {
            local[u] = comp_verts[c].len();
            comp_verts[c].push(u);
        }

        let num_comps = comp_verts.len();
        let mut comp_datives = vec![0usize; num_comps];
        for &atom in &raw.dative_bonds {
            comp_datives[component[atom]] += 1;
        }
        let mut comp_stereo: Vec<Vec<Stereo>> = (0..num_comps).map(|_| Vec::new()).collect();
        for element in stereo {
            let [atom, _] = element.anchor(|atom| atom);
            comp_stereo[component[atom]].push(element.mapped(&local));
        }

        let mut out = Vec::with_capacity(num_comps);
        for (c, verts) in comp_verts.iter().enumerate() {
            let cn = verts.len();
            let mut offsets = Vec::with_capacity(cn + 1);
            let mut keys = Vec::with_capacity(cn);
            let mut twin_spans = if raw.twin_pool.is_empty() {
                Vec::new()
            } else {
                Vec::with_capacity(cn)
            };
            let mut twin_pool = Vec::new();
            offsets.push(0);
            let mut total_edges = 0usize;
            for &u in verts {
                let deg = raw.offsets[u + 1] - raw.offsets[u];
                total_edges += deg;
                offsets.push(total_edges);
                keys.push(raw.keys[u]);
                if !raw.twin_pool.is_empty() {
                    let t = raw.twins_of(u);
                    let s = u32::try_from(twin_pool.len()).unwrap_or(u32::MAX);
                    let l = u32::try_from(t.len()).unwrap_or(u32::MAX);
                    twin_spans.push((s, l));
                    twin_pool.extend_from_slice(t);
                }
            }
            let mut edges = Vec::with_capacity(total_edges);
            for &u in verts {
                for &(code, v) in &raw.edges[raw.offsets[u]..raw.offsets[u + 1]] {
                    edges.push((code, local[v]));
                }
            }
            let sub_raw = FoldedRaw {
                n: cn,
                offsets,
                edges,
                keys,
                twin_spans,
                twin_pool,
                dative_bonds: Vec::new(),
                fates: Vec::new(),
            };
            out.push(Folded::from_raw(
                sub_raw,
                comp_datives[c],
                std::mem::take(&mut comp_stereo[c]),
            ));
        }
        out
    }

    /// Whether `target` is the same component using refined colors, incremental stereo checks,
    /// and lazy Indigo-style target stabilizer orbit pruning.
    fn matches(&self, target: &Folded) -> bool {
        if self.csr.n != target.csr.n
            || self.csr.edges.len() != target.csr.edges.len()
            || self.dative_count != target.dative_count
            || self.color_hist != target.color_hist
        {
            return false;
        }
        let n = self.csr.n;
        if n == 0 {
            return true;
        }

        // Fast O(n + m) path when refined colors already force a unique 1-to-1 bijection.
        // Distinct vertex colors forbid any graph automorphism, so no stereo element can be fake.
        if self.refined_hashes.len() == n {
            if self.stereo.len() != target.stereo.len() {
                return false;
            }
            let mut target_by_color = vec![0usize; n];
            for (t, &c) in target.refined_colors.iter().enumerate() {
                target_by_color[c] = t;
            }
            let mapping: Vec<usize> = self
                .refined_colors
                .iter()
                .map(|&c| target_by_color[c])
                .collect();
            let atoms_and_edges_ok = (0..n).all(|q| {
                let t = mapping[q];
                self.atom_equal(q, target, t)
                    && self.csr.neighbors(q).iter().all(|&(q_code, qn)| {
                        target.csr.edge_between(t, mapping[qn]) == Some(q_code)
                    })
            });
            if !atoms_and_edges_ok {
                return false;
            }
            return self.stereo.iter().all(|el| {
                target
                    .anchor_get(el.anchor(|a| mapping[a]))
                    .is_some_and(|j| el.agrees(&mapping, &target.stereo[j]))
            });
        }

        let q_fake = self.fake();
        let t_fake = target.fake();
        let q_real_count = self.stereo.len() - q_fake.iter().filter(|&&f| f).count();
        let t_real_count = target.stereo.len() - t_fake.iter().filter(|&&f| f).count();
        if q_real_count != t_real_count {
            return false;
        }

        let spec = ExactMatchSpec::new(self, target, q_fake, t_fake);
        let domain_sizes: Vec<usize> = self
            .refined_colors
            .iter()
            .map(|&c| target.color_hist[c].1)
            .collect();
        let plan = plan_query_order(&self.csr, &domain_sizes);

        // Try fast unpruned search first; only compute full SymmetryGroup when symmetry-induced
        // backtracking exceeds a tiny budget.
        match find_embedding_unpruned(&self.csr, &target.csr, &plan, &spec, 64) {
            EmbedOutcome::Matched(_) => return true,
            EmbedOutcome::NoMatch => return false,
            EmbedOutcome::BudgetExceeded => {}
        }

        let t_sym = target.real_sym.get_or_init(|| {
            let t_real: Vec<usize> = (0..target.stereo.len())
                .filter(|&i| !is_fake_at(t_fake, i))
                .collect();
            target.compute_sym(&t_real)
        });
        find_embedding(&self.csr, &target.csr, t_sym, &plan, &spec).is_some()
    }

    fn read_coordination(
        fates: &[Fate],
        center: NodeIndex,
        ligands: &[Ligand],
        configuration: Configuration,
    ) -> Result<Stereo, ExactMatchError> {
        let Fate::Kept(owner) = fates[center.index()] else {
            return Err(ExactMatchError::UnsupportedStereo);
        };
        let len = ligands.len();
        if len > 6 {
            return Err(ExactMatchError::UnsupportedStereo);
        }
        let mut sites = [Site::Vacancy; 6];
        for (slot, &ligand) in sites[..len].iter_mut().zip(ligands) {
            *slot = match ligand {
                Ligand::Atom(node) => match fates[node.index()] {
                    Fate::Kept(atom) => Site::Atom(atom),
                    Fate::Hydrogen => Site::Hydrogen(owner),
                    Fate::Twin => Site::Twin,
                },
                Ligand::Hydrogen => Site::Hydrogen(owner),
                Ligand::LonePair => Site::LonePair,
                Ligand::Vacancy => Site::Vacancy,
            };
        }
        if configuration.coordination_geometry(&sites[..len]).is_none() {
            return Err(ExactMatchError::UnsupportedStereo);
        }
        Ok(Stereo::Coordination {
            center: owner,
            ligands: sites,
            len,
            configuration,
        })
    }

    /// The stereo elements in folded-graph atoms, without those that two alike ligands
    /// rule out.
    fn read_stereo(molecule: &Molecule, raw: &FoldedRaw) -> Result<Vec<Stereo>, ExactMatchError> {
        let fates = &raw.fates;
        // a hydrogen or twin is an end atom, so it is never a centre, an allene terminal or a
        // double-bond end
        let atom = |node: NodeIndex| match fates[node.index()] {
            Fate::Kept(atom) => Ok(atom),
            _ => Err(ExactMatchError::UnsupportedStereo),
        };
        // a ligand on the folded atom `owner`
        let site = |ligand: Ligand, owner: usize| match ligand {
            Ligand::Atom(node) => Ok(match fates[node.index()] {
                Fate::Kept(atom) => Site::Atom(atom),
                Fate::Hydrogen => Site::Hydrogen(owner),
                Fate::Twin => Site::Twin,
            }),
            Ligand::Hydrogen => Ok(Site::Hydrogen(owner)),
            Ligand::LonePair => Ok(Site::LonePair),
            Ligand::Vacancy => Ok(Site::Vacancy),
        };
        // tag 1 is @, 2 is @@; `owners` are the atoms the ligands are on
        let parity = |center: NodeIndex, ligands: &[Ligand], owners: [usize; 4], tag: u8| {
            if ligands.len() != 4 {
                return Err(ExactMatchError::UnsupportedStereo);
            }
            let mut sites = [Site::LonePair; 4];
            for ((slot, &ligand), owner) in sites.iter_mut().zip(ligands).zip(owners) {
                *slot = site(ligand, owner)?;
            }
            Ok(Stereo::Parity {
                center: atom(center)?,
                ligands: sites,
                clockwise: match tag {
                    1 => false,
                    2 => true,
                    _ => return Err(ExactMatchError::UnsupportedStereo),
                },
            })
        };
        let mut stereo = Vec::new();
        for center in molecule.atom_stereo() {
            match center.configuration {
                Configuration::Tetrahedral(tag) => {
                    stereo.push(parity(
                        center.center,
                        &center.ligands,
                        [atom(center.center)?; 4],
                        tag,
                    )?);
                }
                Configuration::SquarePlanar(_)
                | Configuration::TrigonalBipyramidal(_)
                | Configuration::Octahedral(_) => {
                    stereo.push(Self::read_coordination(
                        fates,
                        center.center,
                        &center.ligands,
                        center.configuration,
                    )?);
                }
            }
        }
        for allene in molecule.allene_stereo() {
            let terminal = |end: Option<&NodeIndex>| end.map(|node| fates[node.index()]);
            let (Some(Fate::Kept(first)), Some(Fate::Kept(last))) =
                (terminal(allene.axis.first()), terminal(allene.axis.last()))
            else {
                // a folded terminal carries only hydrogens, so there is no stereo
                continue;
            };
            stereo.push(parity(
                allene.center,
                &allene.ligands,
                [first, first, last, last],
                allene.permutation,
            )?);
        }
        for bond in molecule.double_bond_stereo() {
            let ends = [atom(bond.ends[0])?, atom(bond.ends[1])?];
            let references = [
                site(Ligand::Atom(bond.references[0]), ends[0])?,
                site(Ligand::Atom(bond.references[1]), ends[1])?,
            ];
            stereo.push(Stereo::DoubleBond {
                ends,
                references,
                same_side: bond.same_side,
            });
        }
        // two alike ligands on one atom, two hydrogens or folded twins, make it no stereo at
        // all; twins always come at least two to an atom
        stereo.retain(|element| match element {
            Stereo::Parity { ligands, .. } => {
                (0..4).all(|i| (i + 1..4).all(|j| ligands[i] != ligands[j]))
            }
            Stereo::Coordination {
                ligands,
                len,
                configuration,
                ..
            } => !configuration.all_coordination_trivial(&ligands[..*len]),
            Stereo::DoubleBond { ends, .. } => ends
                .iter()
                .all(|&end| raw.keys[end].hydrogens < 2 && raw.twins_of(end).is_empty()),
        });
        Ok(stereo)
    }

    fn compute_sym(&self, active_indices: &[usize]) -> SymmetryGroup {
        let mut stereo_free = vec![true; self.csr.n];
        let mut verts = Vec::new();
        for &i in active_indices {
            self.stereo[i].collect_vertices(&mut verts);
            for &v in &verts {
                stereo_free[v] = false;
            }
        }
        let mut leaf_items: Vec<[u64; 8]> = Vec::with_capacity(active_indices.len());
        SymmetryGroup::compute(
            &self.csr,
            &stereo_free,
            |colors, buf| {
                self.push_stereo_refinement_keys(active_indices, colors, buf);
            },
            |rank_of, sig| {
                leaf_items.clear();
                for &i in active_indices {
                    leaf_items.push(self.stereo[i].encode_for_leaf(rank_of, &self.csr, &self.keys));
                }
                leaf_items.sort_unstable();
                for item in &leaf_items {
                    sig.extend_from_slice(item);
                }
            },
        )
    }

    fn push_stereo_refinement_keys(
        &self,
        active_indices: &[usize],
        colors: &[usize],
        buf: &mut [Vec<u64>],
    ) {
        let map_site = |site: Site| match site {
            Site::Atom(a) => RankedSite::Atom(colors[a]),
            Site::Hydrogen(a) => RankedSite::Hydrogen(colors[a]),
            Site::LonePair => RankedSite::LonePair,
            Site::Twin => RankedSite::Twin,
            Site::Vacancy => RankedSite::Vacancy,
        };
        for &i in active_indices {
            match self.stereo[i] {
                Stereo::Parity {
                    center,
                    ligands,
                    clockwise,
                } => {
                    let sites = ligands.map(map_site);
                    let mut distinct = true;
                    let mut odd = false;
                    for a in 0..4 {
                        for b in a + 1..4 {
                            if sites[a] == sites[b] {
                                distinct = false;
                                break;
                            }
                            odd ^= sites[a] > sites[b];
                        }
                    }
                    if distinct {
                        buf[center].push(10 + u64::from(clockwise ^ odd));
                    }
                }
                Stereo::Coordination {
                    center,
                    ligands,
                    len,
                    configuration,
                } => {
                    let mut sites = [RankedSite::Vacancy; 6];
                    for (slot, &s) in sites[..len].iter_mut().zip(&ligands[..len]) {
                        *slot = map_site(s);
                    }
                    if !configuration.all_coordination_trivial(&sites[..len]) {
                        let code = encode_coordination_leaf(
                            colors[center],
                            &ligands[..len],
                            configuration,
                            map_site,
                        );
                        let mut h = code[0];
                        for &word in &code[1..] {
                            h = mix_word(h, word);
                        }
                        buf[center].push(mix64(h));
                    }
                }
                Stereo::DoubleBond {
                    ends,
                    references,
                    same_side,
                } => self.push_double_bond_key(ends, references, same_side, colors, buf),
            }
        }
        for row in buf.iter_mut() {
            if row.len() > 1 {
                row.sort_unstable();
            }
        }
    }

    fn push_double_bond_key(
        &self,
        ends: [usize; 2],
        references: [Site; 2],
        same_side: bool,
        colors: &[usize],
        buf: &mut [Vec<u64>],
    ) {
        let map_site = |site: Site| match site {
            Site::Atom(a) => RankedSite::Atom(colors[a]),
            Site::Hydrogen(a) => RankedSite::Hydrogen(colors[a]),
            Site::LonePair => RankedSite::LonePair,
            Site::Twin => RankedSite::Twin,
            Site::Vacancy => RankedSite::Vacancy,
        };
        let mut parity = Some(same_side);
        for side in [0, 1] {
            let end = ends[side];
            let partner = ends[1 - side];
            let ref_c = map_site(references[side]);
            let mut other: Option<RankedSite> = None;
            for &(_, v) in self.csr.neighbors(end) {
                if v != partner {
                    let sc = RankedSite::Atom(colors[v]);
                    if sc != ref_c {
                        if other.is_some() {
                            parity = None;
                            break;
                        }
                        other = Some(sc);
                    }
                }
            }
            if self.keys[end].hydrogens == 1 {
                let hc = RankedSite::Hydrogen(colors[end]);
                if hc != ref_c {
                    if other.is_some() {
                        parity = None;
                    } else {
                        other = Some(hc);
                    }
                }
            }
            if let (Some(p), Some(o)) = (parity, other) {
                if o == ref_c {
                    parity = None;
                } else if ref_c > o {
                    parity = Some(!p);
                }
            }
        }
        if let Some(p) = parity {
            let p_word = u64::from(p);
            buf[ends[0]].push(mix64(mix_word(
                mix_word(20, colors[ends[1]] as u64),
                p_word,
            )));
            buf[ends[1]].push(mix64(mix_word(
                mix_word(20, colors[ends[0]] as u64),
                p_word,
            )));
        }
    }

    fn fake(&self) -> &[bool] {
        self.fake.get_or_init(|| self.find_fake())
    }

    fn element_has_tied_colors(&self, index: usize) -> bool {
        let map_const_site = |site: Site| match site {
            Site::Atom(a) => RankedSite::Atom(self.refined_colors[a]),
            Site::Hydrogen(h) => RankedSite::Hydrogen(self.refined_colors[h]),
            Site::LonePair => RankedSite::LonePair,
            Site::Twin => RankedSite::Twin,
            Site::Vacancy => RankedSite::Vacancy,
        };
        match self.stereo[index] {
            Stereo::Parity { ligands, .. } => {
                let s = ligands.map(map_const_site);
                (0..4).any(|a| (a + 1..4).any(|b| s[a] == s[b]))
            }
            Stereo::Coordination { ligands, len, .. } => {
                let mut s = [RankedSite::Vacancy; 6];
                for (slot, &l) in s[..len].iter_mut().zip(&ligands[..len]) {
                    *slot = map_const_site(l);
                }
                !distinct_sites(&s[..len])
            }
            Stereo::DoubleBond { ends, .. } => ends.iter().any(|&end| {
                let partner = if end == ends[0] { ends[1] } else { ends[0] };
                let mut first: Option<usize> = None;
                for &(_, v) in self.csr.neighbors(end) {
                    if v != partner {
                        let c = self.refined_colors[v];
                        if let Some(prev) = first {
                            if prev == c {
                                return true;
                            }
                        } else {
                            first = Some(c);
                        }
                    }
                }
                false
            }),
        }
    }

    /// Marks each stereo element that a symmetry of the molecule inverts while keeping all
    /// other active stereo fixed, iterating until convergence for cascading fake elements.
    /// Returns an empty `Vec` when no element is fake to avoid allocating on typical molecules.
    fn find_fake(&self) -> Vec<bool> {
        let m = self.stereo.len();
        if m == 0 || !(0..m).any(|i| self.element_has_tied_colors(i)) {
            return Vec::new();
        }
        let mut fake = vec![false; m];
        let mut any_fake = false;

        loop {
            let active: Vec<usize> = (0..m).filter(|&i| !fake[i]).collect();
            let tied: Vec<usize> = active
                .iter()
                .copied()
                .filter(|&i| self.element_has_tied_colors(i))
                .collect();
            if tied.is_empty() {
                break;
            }
            let sym = self.compute_sym(&active);
            let mut newly_fake = Vec::new();
            for &i in &tied {
                if self.is_element_fake(i, &active, &sym) {
                    newly_fake.push(i);
                }
            }
            if newly_fake.is_empty() {
                break;
            }
            any_fake = true;
            for i in newly_fake {
                fake[i] = true;
            }
        }
        if any_fake {
            fake
        } else {
            Vec::new()
        }
    }

    fn is_element_fake(&self, index: usize, active: &[usize], sym: &SymmetryGroup) -> bool {
        let target = &self.stereo[index];
        let (a0, a1) = target.anchor_atoms();
        let mut domain_sizes: Vec<usize> = self
            .refined_colors
            .iter()
            .map(|&c| self.color_hist[c].1)
            .collect();
        domain_sizes[a0] = if a1.is_some() { 2 } else { 1 };
        let plan = plan_query_order(&self.csr, &domain_sizes);

        let check_variant = |variant: Stereo| -> bool {
            let spec = FakeCheckSpec::new(self, active, index, variant);
            find_embedding(&self.csr, &self.csr, sym, &plan, &spec).is_some()
        };

        match *target {
            Stereo::Parity {
                center,
                ligands,
                clockwise,
            } => check_variant(Stereo::Parity {
                center,
                ligands,
                clockwise: !clockwise,
            }),
            Stereo::DoubleBond {
                ends,
                references,
                same_side,
            } => check_variant(Stereo::DoubleBond {
                ends,
                references,
                same_side: !same_side,
            }),
            Stereo::Coordination {
                center,
                ligands,
                len,
                configuration,
            } => {
                for num in 1..=configuration.coordination_variant_count() {
                    let Some(next) = configuration.with_coordination_number(num) else {
                        continue;
                    };
                    if next == configuration {
                        continue;
                    }
                    if !check_variant(Stereo::Coordination {
                        center,
                        ligands,
                        len,
                        configuration: next,
                    }) {
                        return false;
                    }
                }
                true
            }
        }
    }
}

fn distinct_sites(sites: &[RankedSite]) -> bool {
    (0..sites.len()).all(|i| (i + 1..sites.len()).all(|j| sites[i] != sites[j]))
}

struct StereoTrigger {
    req_verts: Vec<usize>,
    stereo_idx: usize,
}

fn build_triggers(
    n: usize,
    stereo: &[Stereo],
    indices: &[usize],
) -> (Vec<StereoTrigger>, Vec<Vec<usize>>) {
    if indices.is_empty() {
        return (Vec::new(), Vec::new());
    }
    let mut triggers = Vec::with_capacity(indices.len());
    let mut by_vertex = vec![Vec::new(); n];
    for &idx in indices {
        let mut req = Vec::new();
        stereo[idx].collect_vertices(&mut req);
        let pos = triggers.len();
        for &v in &req {
            by_vertex[v].push(pos);
        }
        triggers.push(StereoTrigger {
            req_verts: req,
            stereo_idx: idx,
        });
    }
    (triggers, by_vertex)
}

struct ExactMatchSpec<'a> {
    q: &'a Folded,
    t: &'a Folded,
    t_fake: &'a [bool],
    triggers: Vec<StereoTrigger>,
    by_vertex: Vec<Vec<usize>>,
}

impl<'a> ExactMatchSpec<'a> {
    fn new(q: &'a Folded, t: &'a Folded, q_fake: &[bool], t_fake: &'a [bool]) -> Self {
        if q.stereo.is_empty() {
            return Self {
                q,
                t,
                t_fake,
                triggers: Vec::new(),
                by_vertex: Vec::new(),
            };
        }
        let q_real: Vec<usize> = (0..q.stereo.len())
            .filter(|&i| !is_fake_at(q_fake, i))
            .collect();
        let (triggers, by_vertex) = build_triggers(q.csr.n, &q.stereo, &q_real);
        Self {
            q,
            t,
            t_fake,
            triggers,
            by_vertex,
        }
    }
}

impl EmbeddingSpec for ExactMatchSpec<'_> {
    fn atom_compatible(&self, q: usize, t: usize) -> bool {
        self.q.refined_colors[q] == self.t.refined_colors[t] && self.q.atom_equal(q, self.t, t)
    }

    fn bond_compatible(&self, q_code: u8, t_code: u8) -> bool {
        q_code == t_code
    }

    fn step_ok(&self, q: usize, mapping: &[usize]) -> bool {
        if self.by_vertex.is_empty() {
            return true;
        }
        for &trig_idx in &self.by_vertex[q] {
            let trig = &self.triggers[trig_idx];
            if trig.req_verts.iter().any(|&v| mapping[v] == usize::MAX) {
                continue;
            }
            let el = &self.q.stereo[trig.stereo_idx];
            let Some(j) = self.t.anchor_get(el.anchor(|a| mapping[a])) else {
                return false;
            };
            if is_fake_at(self.t_fake, j) || !el.agrees(mapping, &self.t.stereo[j]) {
                return false;
            }
        }
        true
    }

    fn leaf_ok(&self, _mapping: &[usize]) -> bool {
        true
    }
}

struct FakeCheckSpec<'a> {
    mol: &'a Folded,
    is_active: Vec<bool>,
    skip_idx: usize,
    variant: Stereo,
    anchor_a: usize,
    anchor_b: Option<usize>,
    triggers: Vec<StereoTrigger>,
    by_vertex: Vec<Vec<usize>>,
}

impl<'a> FakeCheckSpec<'a> {
    fn new(mol: &'a Folded, active: &[usize], skip_idx: usize, variant: Stereo) -> Self {
        let mut is_active = vec![false; mol.stereo.len()];
        for &i in active {
            is_active[i] = true;
        }
        let (anchor_a, anchor_b) = mol.stereo[skip_idx].anchor_atoms();
        let (triggers, by_vertex) = build_triggers(mol.csr.n, &mol.stereo, active);
        Self {
            mol,
            is_active,
            skip_idx,
            variant,
            anchor_a,
            anchor_b,
            triggers,
            by_vertex,
        }
    }
}

impl EmbeddingSpec for FakeCheckSpec<'_> {
    fn atom_compatible(&self, q: usize, t: usize) -> bool {
        if q == self.anchor_a && t != self.anchor_a && self.anchor_b != Some(t) {
            return false;
        }
        if self.anchor_b == Some(q) && t != self.anchor_a && self.anchor_b != Some(t) {
            return false;
        }
        self.mol.refined_colors[q] == self.mol.refined_colors[t]
            && self.mol.atom_equal(q, self.mol, t)
    }

    fn bond_compatible(&self, q_code: u8, t_code: u8) -> bool {
        q_code == t_code
    }

    fn step_ok(&self, q: usize, mapping: &[usize]) -> bool {
        if self.by_vertex.is_empty() {
            return true;
        }
        for &trig_idx in &self.by_vertex[q] {
            let trig = &self.triggers[trig_idx];
            if trig.req_verts.iter().any(|&v| mapping[v] == usize::MAX) {
                continue;
            }
            let el = &self.mol.stereo[trig.stereo_idx];
            let Some(j) = self.mol.anchor_get(el.anchor(|a| mapping[a])) else {
                return false;
            };
            if !self.is_active[j] {
                return false;
            }
            let target_el = if j == self.skip_idx {
                &self.variant
            } else {
                &self.mol.stereo[j]
            };
            if !el.agrees(mapping, target_el) {
                return false;
            }
        }
        true
    }

    fn leaf_ok(&self, _mapping: &[usize]) -> bool {
        true
    }
}

/// Fast ephemeral atom hash used by `Folded` (never persisted to disk).
#[inline]
fn fast_atom_code(base: AtomBase, twins: &[(AtomBase, u8)]) -> u64 {
    let mut h = mix_word(base.packed(), twins.len() as u64);
    for &(twin_base, twin_bond) in twins {
        h = mix_word(h, twin_base.packed());
        h = mix_word(h, u64::from(twin_bond));
    }
    mix64(h)
}

// explicit codes, so reordering BondOrder cannot change stored keys
#[inline]
fn order_code_from_pack(packed: u8) -> u64 {
    match packed >> 1 {
        0 => 1,
        1 => 2,
        2 => 3,
        3 => 4,
        4 if (packed & 1) != 0 => 5,
        _ => 6,
    }
}

/// Stable FNV-1a atom hash matching stored `.key` index files (`VERSION = 1`).
fn stable_atom_code(base: AtomBase, twins: &[(AtomBase, u8)]) -> u64 {
    let twins_iter = twins.iter().flat_map(|&(twin_base, twin_bond)| {
        [
            stable_atom_code(twin_base, &[]),
            order_code_from_pack(twin_bond),
        ]
    });
    fnv1a(
        [
            u64::from(base.atomic_number),
            u64::from(base.isotope),
            // sign-extend, like the fingerprints
            i64::from(base.charge).cast_unsigned(),
            u64::from(base.hydrogens),
            u64::from(base.aromatic),
            twins.len() as u64,
        ]
        .into_iter()
        .chain(twins_iter),
    )
}

impl Molecule {
    fn folded_components(&self) -> Result<&[Folded], ExactMatchError> {
        let cached = self
            .exact_cache
            .get_or_init(|| Arc::new(Folded::components(self)));
        cached.as_ref().as_deref().map_err(|error| *error)
    }

    /// A 64-bit key that molecules `exact_match` calls the same always share, whatever
    /// their stereo. Atom codes are refined from their neighbours' codes until that splits
    /// no more atoms, and the key hashes the sorted codes. It depends neither on atom order
    /// nor on the Rust version, so it can be stored.
    pub fn exact_key(&self) -> u64 {
        let raw = fold_raw(self);
        let n = raw.n;
        let mut codes: Vec<u64> = (0..n)
            .map(|u| stable_atom_code(raw.keys[u], raw.twins_of(u)))
            .collect();
        let mut scratch = Vec::with_capacity(n);
        let distinct_len = |scratch: &mut Vec<u64>, values: &[u64]| {
            scratch.clear();
            scratch.extend_from_slice(values);
            scratch.sort_unstable();
            scratch.dedup();
            scratch.len()
        };
        let mut classes = distinct_len(&mut scratch, &codes);
        let mut next = vec![0u64; n];
        let mut neighbours: Vec<(u64, u64)> = Vec::with_capacity(8);
        loop {
            for u in 0..n {
                neighbours.clear();
                for &(packed, v) in &raw.edges[raw.offsets[u]..raw.offsets[u + 1]] {
                    neighbours.push((order_code_from_pack(packed), codes[v]));
                }
                neighbours.sort_unstable();
                next[u] = fnv1a(
                    std::iter::once(codes[u])
                        .chain(neighbours.iter().flat_map(|&(bond, code)| [bond, code])),
                );
            }
            std::mem::swap(&mut codes, &mut next);
            let next_classes = distinct_len(&mut scratch, &codes);
            // refining never merges classes, short of a hash collision
            if next_classes <= classes {
                break;
            }
            classes = next_classes;
        }
        codes.sort_unstable();
        fnv1a(codes)
    }

    /// Whether `other` is the same molecule, however either was written: atom order,
    /// ring-closure digits, Kekulé or aromatic form, hydrogens as atoms or counts. Atoms
    /// must agree on element, isotope, charge, hydrogens and aromaticity, bonds on order
    /// and dative direction, and tetrahedral, allene and double-bond stereo must agree.
    /// Atom classes are labels, so they are ignored.
    pub fn exact_match(&self, other: &Molecule) -> Result<bool, ExactMatchError> {
        let query = self.folded_components()?;
        let target = other.folded_components()?;
        if query.len() != target.len() {
            return Ok(false);
        }
        if query.len() == 1 {
            return Ok(query[0].matches(&target[0]));
        }
        let mut used = vec![false; target.len()];
        // the first alike target component will do: being the same is an equivalence
        'parts: for part in query {
            for (slot_idx, other_part) in target.iter().enumerate() {
                if !used[slot_idx] && part.matches(other_part) {
                    used[slot_idx] = true;
                    continue 'parts;
                }
            }
            return Ok(false);
        }
        Ok(true)
    }
}
