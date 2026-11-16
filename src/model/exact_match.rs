//! Whether two molecules are the same, however they were written (#58).
use crate::descriptors::fingerprints::ecfp::fnv1a;
use crate::model::bond::{Bond, BondOrder};
use crate::model::matcher::{
    find_embedding, pack_bond, plan_query_order, refine_hash_signature, CsrGraph, EmbeddingSpec,
    SymmetryGroup,
};
use crate::model::molecule::Molecule;
use crate::model::stereo::{Configuration, CoordinationGeometry, Ligand};
use petgraph::graph::{Graph, NodeIndex as GraphIndex};
use petgraph::stable_graph::{EdgeReference, NodeIndex};
use petgraph::visit::{EdgeRef, IntoEdgeReferences, NodeIndexable};
use petgraph::Undirected;
use std::cell::OnceCell;
use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExactMatchError {
    /// Square-planar, trigonal-bipyramidal or octahedral stereo, which isn't compared yet.
    UnsupportedStereo,
    /// More than the search limit of atom mappings to try.
    TooSymmetric,
}

#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct AtomKey {
    atomic_number: u8,
    isotope: u16,
    charge: i8,
    hydrogens: usize,
    aromatic: bool,
    /// Folded twin end atoms, with their bonds, sorted.
    twins: Vec<(AtomKey, BondOrder)>,
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

    fn encode_for_leaf(&self, rank_of: &[usize], csr: &CsrGraph, keys: &[AtomKey]) -> [u64; 8] {
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
struct Folded {
    csr: CsrGraph,
    keys: Vec<AtomKey>,
    refined_colors: Vec<usize>,
    refined_hashes: Vec<u64>,
    color_hist: Vec<(u64, usize)>,
    /// (donor, acceptor) of each dative bond.
    datives: Vec<(usize, usize)>,
    stereo: Vec<Stereo>,
    /// Index into `stereo` by anchor.
    anchors: HashMap<[usize; 2], usize>,
    fake: OnceCell<Result<Vec<bool>, ExactMatchError>>,
    real_sym: OnceCell<SymmetryGroup>,
}

type FoldedGraph = Graph<AtomKey, BondOrder, Undirected>;

impl Folded {
    fn new(graph: &FoldedGraph, datives: Vec<(usize, usize)>, stereo: Vec<Stereo>) -> Folded {
        let n = graph.node_count();
        let keys: Vec<AtomKey> = (0..n).map(|i| graph[GraphIndex::new(i)].clone()).collect();
        let base_hashes: Vec<u64> = keys.iter().map(atom_code).collect();
        let mut adj: Vec<Vec<(u8, usize)>> = vec![Vec::new(); n];
        for edge in graph.edge_references() {
            let (u, v) = (edge.source().index(), edge.target().index());
            let order = *edge.weight();
            adj[u].push((pack_bond(order, datives.contains(&(u, v))), v));
            adj[v].push((pack_bond(order, datives.contains(&(v, u))), u));
        }
        let csr = CsrGraph::from_adj(n, &adj, base_hashes);
        let (refined_colors, refined_hashes) = refine_hash_signature(&csr);
        let mut counts = vec![0usize; refined_hashes.len()];
        for &c in &refined_colors {
            counts[c] += 1;
        }
        let color_hist: Vec<(u64, usize)> = refined_hashes.iter().copied().zip(counts).collect();
        let anchors = stereo
            .iter()
            .enumerate()
            .map(|(i, element)| (element.anchor(|atom| atom), i))
            .collect();
        Folded {
            csr,
            keys,
            refined_colors,
            refined_hashes,
            color_hist,
            datives,
            stereo,
            anchors,
            fake: OnceCell::new(),
            real_sym: OnceCell::new(),
        }
    }

    /// One `Folded` per connected component, so k alike salt or water parts cost k
    /// matches instead of k! mappings. No stereo element spans two components.
    fn components(molecule: &Molecule) -> Result<Vec<Folded>, ExactMatchError> {
        let (graph, fates, datives) = Folded::fold(molecule);
        let stereo = Folded::read_stereo(molecule, &graph, &fates)?;
        let mut component = vec![usize::MAX; graph.node_count()];
        // each atom's index within its component
        let mut local = vec![0; graph.node_count()];
        let mut parts: Vec<FoldedGraph> = Vec::new();
        for start in graph.node_indices() {
            if component[start.index()] != usize::MAX {
                continue;
            }
            let mut part = Graph::new_undirected();
            let mut stack = vec![start];
            component[start.index()] = parts.len();
            while let Some(node) = stack.pop() {
                local[node.index()] = part.add_node(graph[node].clone()).index();
                for next in graph.neighbors(node) {
                    if component[next.index()] == usize::MAX {
                        component[next.index()] = parts.len();
                        stack.push(next);
                    }
                }
            }
            parts.push(part);
        }
        for edge in graph.edge_references() {
            let (a, b) = (edge.source().index(), edge.target().index());
            parts[component[a]].add_edge(
                GraphIndex::new(local[a]),
                GraphIndex::new(local[b]),
                *edge.weight(),
            );
        }
        let mut part_datives = vec![Vec::new(); parts.len()];
        for (donor, acceptor) in datives {
            part_datives[component[donor]].push((local[donor], local[acceptor]));
        }
        let mut part_stereo: Vec<Vec<Stereo>> = parts.iter().map(|_| Vec::new()).collect();
        for element in stereo {
            let [atom, _] = element.anchor(|atom| atom);
            part_stereo[component[atom]].push(element.mapped(&local));
        }
        Ok(parts
            .into_iter()
            .zip(part_datives)
            .zip(part_stereo)
            .map(|((graph, datives), stereo)| Folded::new(&graph, datives, stereo))
            .collect())
    }

    /// Whether `target` is the same component using refined colors, incremental stereo checks,
    /// and Indigo-style target stabilizer orbit pruning.
    fn matches(&self, target: &Folded) -> Result<bool, ExactMatchError> {
        if self.csr.n != target.csr.n
            || self.csr.edges.len() != target.csr.edges.len()
            || self.datives.len() != target.datives.len()
            || self.color_hist != target.color_hist
        {
            return Ok(false);
        }
        let q_fake = self.fake()?;
        let t_fake = target.fake()?;
        let q_real: Vec<usize> = q_fake
            .iter()
            .enumerate()
            .filter_map(|(i, &f)| (!f).then_some(i))
            .collect();
        let t_real: Vec<usize> = t_fake
            .iter()
            .enumerate()
            .filter_map(|(i, &f)| (!f).then_some(i))
            .collect();
        if q_real.len() != t_real.len() {
            return Ok(false);
        }
        let n = self.csr.n;
        if n == 0 {
            return Ok(true);
        }

        // Fast O(n + m) path when refined colors already force a unique 1-to-1 bijection.
        if self.refined_hashes.len() == n {
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
                self.keys[q] == target.keys[t]
                    && self.csr.neighbors(q).iter().all(|&(q_code, qn)| {
                        target.csr.edge_between(t, mapping[qn]) == Some(q_code)
                    })
            });
            if !atoms_and_edges_ok {
                return Ok(false);
            }
            return Ok(q_real.iter().all(|&i| {
                let el = &self.stereo[i];
                target
                    .anchors
                    .get(&el.anchor(|a| mapping[a]))
                    .is_some_and(|&j| !t_fake[j] && el.agrees(&mapping, &target.stereo[j]))
            }));
        }

        let spec = ExactMatchSpec::new(self, target, &q_real, t_fake);
        let t_sym = target.real_sym.get_or_init(|| target.compute_sym(&t_real));
        let mut counts = vec![0usize; target.refined_hashes.len()];
        for &c in &target.refined_colors {
            counts[c] += 1;
        }
        let domain_sizes: Vec<usize> = self.refined_colors.iter().map(|&c| counts[c]).collect();
        let plan = plan_query_order(&self.csr, &domain_sizes);
        Ok(find_embedding(&self.csr, &target.csr, t_sym, &plan, &spec).is_some())
    }

    /// The folded graph, what became of each source atom, and the dative bonds.
    fn fold(molecule: &Molecule) -> (FoldedGraph, Vec<Fate>, Vec<(usize, usize)>) {
        let source = molecule.graph();
        let other = |edge: EdgeReference<'_, Bond>, node: NodeIndex| {
            if edge.source() == node {
                edge.target()
            } else {
                edge.source()
            }
        };
        // a plain hydrogen on one heavy atom; [2H], [H+] and H2 stay atoms
        let hydrogen = |node: NodeIndex| {
            let atom = &source[node];
            let mut edges = source.edges(node);
            atom.element.atomic_number == 1
                && atom.isotope == 0
                && atom.charge == 0
                && atom.hs + atom.implicit_hs == 0
                && matches!((edges.next(), edges.next()), (Some(edge), None)
                    if edge.weight().order == BondOrder::Single
                        && source[other(edge, node)].element.atomic_number != 1)
        };
        let mut keys: Vec<Option<AtomKey>> = vec![None; source.node_bound()];
        for node in source.node_indices().filter(|&node| !hydrogen(node)) {
            let atom = &source[node];
            let folded = source
                .neighbors(node)
                .filter(|&other| hydrogen(other))
                .count();
            keys[node.index()] = Some(AtomKey {
                atomic_number: atom.element.atomic_number,
                isotope: atom.isotope,
                charge: atom.charge,
                hydrogens: usize::from(atom.hs) + usize::from(atom.implicit_hs) + folded,
                aromatic: atom.aromatic,
                twins: Vec::new(),
            });
        }

        // end atoms: one covalent bond to a kept atom, and no other
        let mut ends: HashMap<(NodeIndex, BondOrder, &AtomKey), Vec<NodeIndex>> = HashMap::new();
        for node in source.node_indices() {
            let Some(key) = &keys[node.index()] else {
                continue;
            };
            let mut bonds = source
                .edges(node)
                .filter(|&edge| keys[other(edge, node).index()].is_some());
            if let (Some(edge), None) = (bonds.next(), bonds.next()) {
                if edge.weight().order != BondOrder::Dative {
                    let group = (other(edge, node), edge.weight().order, key);
                    ends.entry(group).or_default().push(node);
                }
            }
        }
        let twins: Vec<(NodeIndex, BondOrder, Vec<NodeIndex>)> = ends
            .into_iter()
            .filter(|(_, group)| group.len() > 1)
            .map(|((parent, order, _), group)| (parent, order, group))
            .collect();
        let mut fates = vec![Fate::Hydrogen; source.node_bound()];
        for (parent, order, group) in twins {
            for twin in group {
                fates[twin.index()] = Fate::Twin;
                if let (Some(key), Some(parent)) =
                    (keys[twin.index()].take(), &mut keys[parent.index()])
                {
                    parent.twins.push((key, order));
                }
            }
        }

        let mut graph = Graph::new_undirected();
        for node in source.node_indices() {
            if let Some(mut key) = keys[node.index()].take() {
                key.twins.sort();
                fates[node.index()] = Fate::Kept(graph.add_node(key).index());
            }
        }
        let mut datives = Vec::new();
        for edge in source.edge_references() {
            let (Fate::Kept(a), Fate::Kept(b)) =
                (fates[edge.source().index()], fates[edge.target().index()])
            else {
                continue;
            };
            graph.add_edge(GraphIndex::new(a), GraphIndex::new(b), edge.weight().order);
            if let Some(donor) = edge.weight().donor {
                datives.push(if donor == edge.source() {
                    (a, b)
                } else {
                    (b, a)
                });
            }
        }
        (graph, fates, datives)
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
    fn read_stereo(
        molecule: &Molecule,
        graph: &FoldedGraph,
        fates: &[Fate],
    ) -> Result<Vec<Stereo>, ExactMatchError> {
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
            Stereo::DoubleBond { ends, .. } => ends.iter().all(|&end| {
                let key = &graph[GraphIndex::new(end)];
                key.hydrogens < 2 && key.twins.is_empty()
            }),
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
                        buf[center].push(fnv1a(code));
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
            buf[ends[0]].push(fnv1a([20, colors[ends[1]] as u64, u64::from(p)]));
            buf[ends[1]].push(fnv1a([20, colors[ends[0]] as u64, u64::from(p)]));
        }
    }

    fn fake(&self) -> Result<&[bool], ExactMatchError> {
        self.fake
            .get_or_init(|| Ok(self.find_fake()))
            .as_deref()
            .map_err(|error| *error)
    }

    /// Marks each stereo element that a symmetry of the molecule inverts while keeping all
    /// other active stereo fixed, iterating until convergence for cascading fake elements.
    fn find_fake(&self) -> Vec<bool> {
        let m = self.stereo.len();
        let mut fake = vec![false; m];
        if m == 0 {
            return fake;
        }
        let n = self.csr.n;
        let mut stereo_refs = vec![0usize; n];
        let mut verts = Vec::new();
        for element in &self.stereo {
            element.collect_vertices(&mut verts);
            for &v in &verts {
                stereo_refs[v] += 1;
            }
        }
        let map_const_site = |site: Site| match site {
            Site::Atom(a) => RankedSite::Atom(self.refined_colors[a]),
            Site::Hydrogen(h) => RankedSite::Hydrogen(self.refined_colors[h]),
            Site::LonePair => RankedSite::LonePair,
            Site::Twin => RankedSite::Twin,
            Site::Vacancy => RankedSite::Vacancy,
        };

        loop {
            let active: Vec<usize> = (0..m).filter(|&i| !fake[i]).collect();
            let tied: Vec<usize> = active
                .iter()
                .copied()
                .filter(|&i| match self.stereo[i] {
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
                        let nbrs: Vec<usize> = self
                            .csr
                            .neighbors(end)
                            .iter()
                            .filter(|&&(_, v)| v != partner)
                            .map(|&(_, v)| self.refined_colors[v])
                            .collect();
                        nbrs.len() >= 2 && nbrs[0] == nbrs[1]
                    }),
                })
                .collect();
            if tied.is_empty() {
                break;
            }
            let sym = self.compute_sym(&active);
            let mut counts = vec![0usize; self.refined_hashes.len()];
            for &c in &self.refined_colors {
                counts[c] += 1;
            }
            let mut newly_fake = Vec::new();
            for &i in &tied {
                if self.is_element_fake(i, &active, &sym, &counts) {
                    newly_fake.push(i);
                }
            }
            if newly_fake.is_empty() {
                break;
            }
            for i in newly_fake {
                fake[i] = true;
            }
        }
        fake
    }

    fn is_element_fake(
        &self,
        index: usize,
        active: &[usize],
        sym: &SymmetryGroup,
        counts: &[usize],
    ) -> bool {
        let target = &self.stereo[index];
        let (a0, a1) = target.anchor_atoms();
        let mut domain_sizes: Vec<usize> = self.refined_colors.iter().map(|&c| counts[c]).collect();
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
    fn new(q: &'a Folded, t: &'a Folded, q_real: &[usize], t_fake: &'a [bool]) -> Self {
        let (triggers, by_vertex) = build_triggers(q.csr.n, &q.stereo, q_real);
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
        self.q.refined_colors[q] == self.t.refined_colors[t] && self.q.keys[q] == self.t.keys[t]
    }

    fn bond_compatible(&self, q_code: u8, t_code: u8) -> bool {
        q_code == t_code
    }

    fn step_ok(&self, q: usize, mapping: &[usize]) -> bool {
        for &trig_idx in &self.by_vertex[q] {
            let trig = &self.triggers[trig_idx];
            if trig.req_verts.iter().any(|&v| mapping[v] == usize::MAX) {
                continue;
            }
            let el = &self.q.stereo[trig.stereo_idx];
            let Some(&j) = self.t.anchors.get(&el.anchor(|a| mapping[a])) else {
                return false;
            };
            if self.t_fake[j] || !el.agrees(mapping, &self.t.stereo[j]) {
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
            && self.mol.keys[q] == self.mol.keys[t]
    }

    fn bond_compatible(&self, q_code: u8, t_code: u8) -> bool {
        q_code == t_code
    }

    fn step_ok(&self, q: usize, mapping: &[usize]) -> bool {
        for &trig_idx in &self.by_vertex[q] {
            let trig = &self.triggers[trig_idx];
            if trig.req_verts.iter().any(|&v| mapping[v] == usize::MAX) {
                continue;
            }
            let el = &self.mol.stereo[trig.stereo_idx];
            let Some(&j) = self.mol.anchors.get(&el.anchor(|a| mapping[a])) else {
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

/// Stable code of an atom: its own properties, then its folded twins in sorted order.
fn atom_code(key: &AtomKey) -> u64 {
    let twins = key
        .twins
        .iter()
        .flat_map(|(twin, order)| [atom_code(twin), order_code(*order, true)]);
    fnv1a(
        [
            u64::from(key.atomic_number),
            u64::from(key.isotope),
            // sign-extend, like the fingerprints
            i64::from(key.charge).cast_unsigned(),
            key.hydrogens as u64,
            u64::from(key.aromatic),
            key.twins.len() as u64,
        ]
        .into_iter()
        .chain(twins),
    )
}

// explicit codes, so reordering BondOrder cannot change stored keys
fn order_code(order: BondOrder, donor: bool) -> u64 {
    match order {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => 4,
        BondOrder::Dative if donor => 5,
        BondOrder::Dative => 6,
    }
}

impl Molecule {
    /// A 64-bit key that molecules `exact_match` calls the same always share, whatever
    /// their stereo. Atom codes are refined from their neighbours' codes until that splits
    /// no more atoms, and the key hashes the sorted codes. It depends neither on atom order
    /// nor on the Rust version, so it can be stored.
    pub fn exact_key(&self) -> u64 {
        let (graph, _, datives) = Folded::fold(self);
        let mut codes: Vec<u64> = graph.node_weights().map(atom_code).collect();
        let distinct = |codes: &[u64]| {
            let mut sorted = codes.to_vec();
            sorted.sort_unstable();
            sorted.dedup();
            sorted.len()
        };
        let mut classes = distinct(&codes);
        loop {
            let next: Vec<u64> = graph
                .node_indices()
                .map(|node| {
                    let mut neighbours: Vec<(u64, u64)> = graph
                        .edges(node)
                        .map(|edge| {
                            let other = if edge.source() == node {
                                edge.target()
                            } else {
                                edge.source()
                            };
                            let donor = datives.contains(&(node.index(), other.index()));
                            (order_code(*edge.weight(), donor), codes[other.index()])
                        })
                        .collect();
                    neighbours.sort_unstable();
                    fnv1a(
                        std::iter::once(codes[node.index()])
                            .chain(neighbours.into_iter().flat_map(|(bond, code)| [bond, code])),
                    )
                })
                .collect();
            let next_classes = distinct(&next);
            codes = next;
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
        let query = Folded::components(self)?;
        let mut target: Vec<Option<Folded>> =
            Folded::components(other)?.into_iter().map(Some).collect();
        if query.len() != target.len() {
            return Ok(false);
        }
        // the first alike target component will do: being the same is an equivalence
        'parts: for part in &query {
            for slot in &mut target {
                if slot.as_ref().map(|other| part.matches(other)).transpose()? == Some(true) {
                    *slot = None;
                    continue 'parts;
                }
            }
            return Ok(false);
        }
        Ok(true)
    }
}
