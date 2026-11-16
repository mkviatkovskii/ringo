//! Deterministic individualization/refinement for canonical SMILES.
//!
//! Stereo participates in refinement once the current colors distinguish its ligands.
use super::writer::{write, write_smiles, write_smiles_ranked, DoubleBondStereoPolicy, WriteError};
use crate::model::bond::{Bond, BondOrder};
use crate::model::molecule::Molecule;
use crate::model::stereo::{
    AlleneStereo, AtomStereo, Configuration, CoordinationGeometry, DoubleBondStereo, Ligand,
};
use petgraph::algo::kosaraju_scc;
use petgraph::stable_graph::NodeIndex;
use petgraph::visit::{EdgeIndexable, EdgeRef, IntoEdgeReferences, NodeIndexable};
use std::borrow::Cow;
use std::collections::BTreeMap;

/// Maximum search nodes used by the default canonicalizer.
pub const DEFAULT_SEARCH_LIMIT: usize = 100_000;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CanonicalError {
    /// An exact answer was not established within the search budget.
    SearchLimit,
    Write(WriteError),
}

/// Canonicalizes the represented graph, folding plain explicit `[H]` atoms into H counts
/// and dropping non-stereogenic stereo annotations. The result is Ringo's canonical spelling,
/// not necessarily another toolkit's spelling.
pub fn canonical_smiles(molecule: &Molecule) -> Result<String, CanonicalError> {
    canonical_smiles_with_limit(molecule, DEFAULT_SEARCH_LIMIT)
}

/// Like [`canonical_smiles`], with an explicit budget of search nodes.
/// Exhausting the budget returns an error, never a provisional canonical string.
pub fn canonical_smiles_with_limit(
    molecule: &Molecule,
    limit: usize,
) -> Result<String, CanonicalError> {
    canonical_smiles_with(molecule, limit, DoubleBondStereoPolicy::Reject)
}

/// Like [`canonical_smiles_with_limit`], with an explicit [`DoubleBondStereoPolicy`].
pub fn canonical_smiles_with(
    molecule: &Molecule,
    limit: usize,
    double_bond_stereo: DoubleBondStereoPolicy,
) -> Result<String, CanonicalError> {
    let mut remaining = limit;
    let groups = kosaraju_scc(&molecule.graph);
    if groups.len() <= 1 {
        return canonical_component_with(molecule, &mut remaining, double_bond_stereo);
    }

    let mut output = components(molecule, &groups)
        .iter()
        .map(|component| canonical_component_with(component, &mut remaining, double_bond_stereo))
        .collect::<Result<Vec<_>, _>>()?;
    output.sort_unstable();
    Ok(output.join("."))
}

fn components(molecule: &Molecule, groups: &[Vec<NodeIndex>]) -> Vec<Molecule> {
    let mut parts: Vec<_> = groups.iter().map(|_| Molecule::new()).collect();
    let mut nodes = vec![(0, NodeIndex::end()); molecule.graph.node_bound()];
    let mut edges = vec![None; molecule.graph.edge_bound()];

    for (group, members) in groups.iter().enumerate() {
        let part = &mut parts[group];
        for &old in members {
            nodes[old.index()] = (group, part.add_atom(molecule.graph[old].clone()));
            part.input_aromatic
                .push(molecule.input_aromatic[old.index()]);
        }
    }
    let node = |old: NodeIndex| nodes[old.index()].1;
    let ligand = |old: Ligand| match old {
        Ligand::Atom(atom) => Ligand::Atom(node(atom)),
        other => other,
    };
    for edge in molecule.graph.edge_references() {
        let a = edge.source();
        let b = edge.target();
        let part = &mut parts[nodes[a.index()].0];
        let bond = edge.weight();
        edges[edge.id().index()] = Some(part.add_bond(
            node(a),
            node(b),
            Bond {
                order: bond.order,
                donor: bond.donor.map(node),
            },
        ));
        part.input_orders
            .push(molecule.input_orders[edge.id().index()]);
    }
    for stereo in &molecule.stereo {
        let part = &mut parts[nodes[stereo.center.index()].0];
        let mut stereo = stereo.clone();
        stereo.center = node(stereo.center);
        stereo
            .ligands
            .iter_mut()
            .for_each(|site| *site = ligand(*site));
        part.stereo.push(stereo);
    }
    for stereo in &molecule.allene_stereo {
        let part = &mut parts[nodes[stereo.center.index()].0];
        let mut stereo = stereo.clone();
        stereo.center = node(stereo.center);
        stereo.axis.iter_mut().for_each(|atom| *atom = node(*atom));
        stereo.ligands = stereo.ligands.map(ligand);
        part.allene_stereo.push(stereo);
    }
    for stereo in &molecule.bond_stereo {
        let part = &mut parts[nodes[stereo.ends[0].index()].0];
        let mut stereo = stereo.clone();
        stereo.bond = edges[stereo.bond.index()]
            .unwrap_or_else(|| unreachable!("every bond belongs to a component"));
        stereo.ends = stereo.ends.map(node);
        stereo.references = stereo.references.map(node);
        part.bond_stereo.push(stereo);
    }
    parts
}

fn bond_code(order: BondOrder, is_donor: bool) -> u8 {
    let base = match order {
        BondOrder::Single => 0,
        BondOrder::Double => 1,
        BondOrder::Triple => 2,
        BondOrder::Aromatic => 3,
        BondOrder::Dative => 4,
    };
    (base << 1) | u8::from(is_donor)
}

fn canonical_component(
    molecule: &Molecule,
    remaining: &mut usize,
) -> Result<String, CanonicalError> {
    let mut search = Search::new(molecule, *remaining);
    let (mut colors, mut distinct) = search.initial_colors();
    distinct = search.refine(&mut colors, distinct);
    search.seed_twin_symmetries(&colors);
    let result = search.visit(colors, distinct);
    *remaining = search.remaining;
    result
}

fn canonical_component_with(
    molecule: &Molecule,
    remaining: &mut usize,
    double_bond_stereo: DoubleBondStereoPolicy,
) -> Result<String, CanonicalError> {
    let mut current = drop_trivial_stereo(fold_explicit_hydrogens(molecule));
    loop {
        let mut search = Search::new(&current, *remaining);
        let (mut colors, mut distinct) = search.initial_colors();
        distinct = search.refine(&mut colors, distinct);

        if let Some((drop_atom, drop_allene, drop_bond)) =
            degree1_twin_stereo_drops(&current, &search, &colors)
        {
            drop(search);
            let owned = current.to_mut();
            retain_marked(&mut owned.stereo, &drop_atom);
            retain_marked(&mut owned.allene_stereo, &drop_allene);
            retain_marked(&mut owned.bond_stereo, &drop_bond);
            continue;
        }

        let (tied_atom, tied_allene, tied_bond) = tied_stereo_indices(&search, &colors);
        search.seed_twin_symmetries(&colors);
        let baseline = match search.visit(colors, distinct) {
            Ok(smiles) => {
                *remaining = search.remaining;
                smiles
            }
            Err(CanonicalError::Write(WriteError::DoubleBondStereo))
                if double_bond_stereo == DoubleBondStereoPolicy::Tolerate =>
            {
                *remaining = search.remaining;
                canonical_component_tolerating_double_bonds(&current, remaining)?
            }
            Err(error) => {
                *remaining = search.remaining;
                return Err(error);
            }
        };
        if tied_atom.is_empty() && tied_allene.is_empty() && tied_bond.is_empty() {
            return Ok(baseline);
        }

        let mut fake_atom = vec![false; current.stereo.len()];
        let mut fake_allene = vec![false; current.allene_stereo.len()];
        let mut fake_bond = vec![false; current.bond_stereo.len()];
        let mut dropped_any = false;

        for i in tied_atom {
            if atom_stereo_is_fake(&current, i, &baseline, remaining, double_bond_stereo)? {
                fake_atom[i] = true;
                dropped_any = true;
            }
        }
        for i in tied_allene {
            let mut flipped = current.as_ref().clone();
            flipped.allene_stereo[i].permutation = 3 - flipped.allene_stereo[i].permutation;
            if canonical_component_solving(&flipped, remaining, double_bond_stereo)? == baseline {
                fake_allene[i] = true;
                dropped_any = true;
            }
        }
        for i in tied_bond {
            let mut flipped = current.as_ref().clone();
            flipped.bond_stereo[i].same_side = !flipped.bond_stereo[i].same_side;
            if canonical_component_solving(&flipped, remaining, double_bond_stereo)? == baseline {
                fake_bond[i] = true;
                dropped_any = true;
            }
        }
        if !dropped_any {
            return Ok(baseline);
        }
        let owned = current.to_mut();
        retain_marked(&mut owned.stereo, &fake_atom);
        retain_marked(&mut owned.allene_stereo, &fake_allene);
        retain_marked(&mut owned.bond_stereo, &fake_bond);
    }
}

fn retain_marked<T>(items: &mut Vec<T>, drop: &[bool]) {
    let mut i = 0;
    items.retain(|_| {
        let keep = !drop[i];
        i += 1;
        keep
    });
}

fn canonical_component_solving(
    molecule: &Molecule,
    remaining: &mut usize,
    double_bond_stereo: DoubleBondStereoPolicy,
) -> Result<String, CanonicalError> {
    match canonical_component(molecule, remaining) {
        Ok(smiles) => Ok(smiles),
        Err(CanonicalError::Write(WriteError::DoubleBondStereo))
            if double_bond_stereo == DoubleBondStereoPolicy::Tolerate =>
        {
            canonical_component_tolerating_double_bonds(molecule, remaining)
        }
        Err(error) => Err(error),
    }
}

fn foldable_hydrogen_parent(molecule: &Molecule, node: NodeIndex) -> Option<NodeIndex> {
    let atom = &molecule.graph[node];
    if atom.element.atomic_number != 1
        || atom.isotope != 0
        || atom.charge != 0
        || atom.atom_class != 0
        || atom.hs + atom.implicit_hs != 0
    {
        return None;
    }
    let mut edges = molecule.graph.edges(node);
    let (Some(edge), None) = (edges.next(), edges.next()) else {
        return None;
    };
    if edge.weight().order != BondOrder::Single {
        return None;
    }
    let parent = if edge.source() == node {
        edge.target()
    } else {
        edge.source()
    };
    (molecule.graph[parent].element.atomic_number != 1).then_some(parent)
}

fn fold_explicit_hydrogens(molecule: &Molecule) -> Cow<'_, Molecule> {
    let bound = molecule.graph.node_bound();
    let mut parent_of_h = vec![None; bound];
    let mut added_hs = vec![0u8; bound];
    let mut any_folded = false;
    for node in molecule.graph.node_indices() {
        if let Some(parent) = foldable_hydrogen_parent(molecule, node) {
            let existing = molecule.graph[parent]
                .hs
                .saturating_add(added_hs[parent.index()]);
            if existing < u8::MAX {
                parent_of_h[node.index()] = Some(parent);
                added_hs[parent.index()] += 1;
                any_folded = true;
            }
        }
    }
    if !any_folded {
        return Cow::Borrowed(molecule);
    }

    let mut folded = Molecule::new();
    let mut node_map = vec![None; bound];
    for old in molecule.graph.node_indices() {
        if parent_of_h[old.index()].is_some() {
            continue;
        }
        let mut atom = molecule.graph[old].clone();
        atom.hs = atom.hs.saturating_add(added_hs[old.index()]);
        node_map[old.index()] = Some(folded.add_atom(atom));
        folded
            .input_aromatic
            .push(molecule.input_aromatic[old.index()]);
    }
    let mut edge_map = vec![None; molecule.graph.edge_bound()];
    for edge in molecule.graph.edge_references() {
        let (Some(a), Some(b)) = (
            node_map[edge.source().index()],
            node_map[edge.target().index()],
        ) else {
            continue;
        };
        let weight = edge.weight();
        edge_map[edge.id().index()] = Some(folded.add_bond(
            a,
            b,
            Bond {
                order: weight.order,
                donor: weight.donor.and_then(|d| node_map[d.index()]),
            },
        ));
        folded
            .input_orders
            .push(molecule.input_orders[edge.id().index()]);
    }
    fold_stereo(molecule, &mut folded, &node_map, &edge_map);
    Cow::Owned(folded)
}

fn fold_stereo(
    molecule: &Molecule,
    folded: &mut Molecule,
    node_map: &[Option<NodeIndex>],
    edge_map: &[Option<petgraph::stable_graph::EdgeIndex>],
) {
    let map_ligand = |ligand: Ligand| match ligand {
        Ligand::Atom(atom) => node_map[atom.index()].map_or(Ligand::Hydrogen, Ligand::Atom),
        other => other,
    };
    for stereo in &molecule.stereo {
        let Some(center) = node_map[stereo.center.index()] else {
            continue;
        };
        folded.stereo.push(AtomStereo {
            center,
            ligands: stereo.ligands.iter().copied().map(map_ligand).collect(),
            configuration: stereo.configuration,
        });
    }
    for stereo in &molecule.allene_stereo {
        let Some(center) = node_map[stereo.center.index()] else {
            continue;
        };
        let Some(axis) = stereo
            .axis
            .iter()
            .map(|atom| node_map[atom.index()])
            .collect::<Option<Vec<_>>>()
        else {
            continue;
        };
        folded.allene_stereo.push(AlleneStereo {
            center,
            axis,
            ligands: stereo.ligands.map(map_ligand),
            permutation: stereo.permutation,
        });
    }
    for stereo in &molecule.bond_stereo {
        let Some(bond) = edge_map[stereo.bond.index()] else {
            continue;
        };
        let (Some(e0), Some(e1)) = (
            node_map[stereo.ends[0].index()],
            node_map[stereo.ends[1].index()],
        ) else {
            continue;
        };
        let ends = [e0, e1];
        let mut references = [NodeIndex::end(); 2];
        let mut same_side = stereo.same_side;
        let mut keep = true;
        for side in [0, 1] {
            let end = ends[side];
            let atom = &folded.graph[end];
            if u16::from(atom.hs) + u16::from(atom.implicit_hs) >= 2 {
                keep = false;
                break;
            }
            if let Some(kept) = node_map[stereo.references[side].index()] {
                references[side] = kept;
            } else {
                let partner = ends[1 - side];
                let mut others = folded.graph.neighbors(end).filter(|&n| n != partner);
                if let (Some(alt), None) = (others.next(), others.next()) {
                    references[side] = alt;
                    same_side = !same_side;
                } else {
                    keep = false;
                    break;
                }
            }
        }
        if keep {
            folded.bond_stereo.push(DoubleBondStereo {
                bond,
                ends,
                references,
                same_side,
            });
        }
    }
}

fn has_duplicate_ligands(ligands: &[Ligand]) -> bool {
    (0..ligands.len()).any(|i| (i + 1..ligands.len()).any(|j| ligands[i] == ligands[j]))
}

fn drop_trivial_stereo(mut molecule: Cow<'_, Molecule>) -> Cow<'_, Molecule> {
    let drop_atom: Vec<bool> = molecule
        .stereo
        .iter()
        .map(|s| match s.configuration {
            Configuration::Tetrahedral(_) => has_duplicate_ligands(&s.ligands),
            Configuration::SquarePlanar(_)
            | Configuration::TrigonalBipyramidal(_)
            | Configuration::Octahedral(_) => s.configuration.all_coordination_trivial(&s.ligands),
        })
        .collect();
    let drop_allene: Vec<bool> = molecule
        .allene_stereo
        .iter()
        .map(|s| s.ligands[0] == s.ligands[1] || s.ligands[2] == s.ligands[3])
        .collect();
    let drop_bond: Vec<bool> = molecule
        .bond_stereo
        .iter()
        .map(|s| {
            s.ends.iter().any(|&end| {
                let atom = &molecule.graph[end];
                u16::from(atom.hs) + u16::from(atom.implicit_hs) >= 2
            })
        })
        .collect();
    if drop_atom.iter().any(|&d| d)
        || drop_allene.iter().any(|&d| d)
        || drop_bond.iter().any(|&d| d)
    {
        let owned = molecule.to_mut();
        retain_marked(&mut owned.stereo, &drop_atom);
        retain_marked(&mut owned.allene_stereo, &drop_allene);
        retain_marked(&mut owned.bond_stereo, &drop_bond);
    }
    molecule
}

fn degree1_parent_and_bond(search: &Search<'_>, node: NodeIndex) -> Option<(usize, u8)> {
    let pos = search.node_pos[node.index()];
    let start = search.adj_offsets[pos];
    let end = search.adj_offsets[pos + 1];
    if end - start != 1 {
        return None;
    }
    let (code, parent) = search.adj_edges[start];
    Some((parent, code))
}

fn stereo_reference_counts(molecule: &Molecule) -> Vec<u16> {
    let mut counts = vec![0u16; molecule.graph.node_bound()];
    for s in &molecule.stereo {
        counts[s.center.index()] += 1;
        for &l in &s.ligands {
            if let Ligand::Atom(n) = l {
                counts[n.index()] += 1;
            }
        }
    }
    for s in &molecule.allene_stereo {
        counts[s.center.index()] += 1;
        for &a in &s.axis {
            counts[a.index()] += 1;
        }
        for l in s.ligands {
            if let Ligand::Atom(n) = l {
                counts[n.index()] += 1;
            }
        }
    }
    for s in &molecule.bond_stereo {
        for n in s.ends.into_iter().chain(s.references) {
            counts[n.index()] += 1;
        }
    }
    counts
}

fn degree1_twin_stereo_drops(
    molecule: &Molecule,
    search: &Search<'_>,
    refined: &[usize],
) -> Option<(Vec<bool>, Vec<bool>, Vec<bool>)> {
    if molecule.stereo.is_empty()
        && molecule.allene_stereo.is_empty()
        && molecule.bond_stereo.is_empty()
    {
        return None;
    }
    let refs = stereo_reference_counts(molecule);
    let twin_key = |ligand: Ligand, owner: NodeIndex| -> Option<(usize, u8)> {
        let Ligand::Atom(node) = ligand else {
            return None;
        };
        if refs[node.index()] > 1 {
            return None;
        }
        let (parent, code) = degree1_parent_and_bond(search, node)?;
        (parent == owner.index()).then_some((refined[node.index()], code))
    };

    let drop_atom: Vec<bool> = molecule
        .stereo
        .iter()
        .map(|s| match s.configuration {
            Configuration::Tetrahedral(_) => (0..s.ligands.len()).any(|i| {
                twin_key(s.ligands[i], s.center).is_some_and(|ki| {
                    (i + 1..s.ligands.len()).any(|j| twin_key(s.ligands[j], s.center) == Some(ki))
                })
            }),
            Configuration::SquarePlanar(_)
            | Configuration::TrigonalBipyramidal(_)
            | Configuration::Octahedral(_) => {
                #[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
                enum TwinOrSite {
                    Ligand(Ligand),
                    Twin(usize, u8),
                }
                let sites: Vec<TwinOrSite> = s
                    .ligands
                    .iter()
                    .copied()
                    .map(|l| {
                        twin_key(l, s.center)
                            .map_or(TwinOrSite::Ligand(l), |(c, b)| TwinOrSite::Twin(c, b))
                    })
                    .collect();
                s.configuration.all_coordination_trivial(&sites)
            }
        })
        .collect();

    let drop_allene: Vec<bool> = molecule
        .allene_stereo
        .iter()
        .map(|s| {
            let pair_fake = |a: Ligand, b: Ligand, term: NodeIndex| {
                matches!(
                    (twin_key(a, term), twin_key(b, term)),
                    (Some(ka), Some(kb)) if ka == kb
                )
            };
            pair_fake(s.ligands[0], s.ligands[1], s.axis[0])
                || pair_fake(
                    s.ligands[2],
                    s.ligands[3],
                    s.axis[s.axis.len().saturating_sub(1)],
                )
        })
        .collect();

    let drop_bond: Vec<bool> = molecule
        .bond_stereo
        .iter()
        .map(|s| {
            [0, 1].into_iter().any(|side| {
                let end = s.ends[side];
                let partner = s.ends[1 - side];
                let Some(rk) = twin_key(Ligand::Atom(s.references[side]), end) else {
                    return false;
                };
                let pos = search.node_pos[end.index()];
                let slice = &search.adj_edges[search.adj_offsets[pos]..search.adj_offsets[pos + 1]];
                slice.iter().any(|&(_, other)| {
                    other != partner.index()
                        && other != s.references[side].index()
                        && twin_key(Ligand::Atom(NodeIndex::new(other)), end) == Some(rk)
                })
            })
        })
        .collect();

    let changed = drop_atom.iter().any(|&d| d)
        || drop_allene.iter().any(|&d| d)
        || drop_bond.iter().any(|&d| d);
    changed.then_some((drop_atom, drop_allene, drop_bond))
}

fn tied_stereo_indices(
    search: &Search<'_>,
    refined: &[usize],
) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    let molecule = search.molecule;
    if molecule.stereo.is_empty()
        && molecule.allene_stereo.is_empty()
        && molecule.bond_stereo.is_empty()
    {
        return (Vec::new(), Vec::new(), Vec::new());
    }
    let ligand_color = |ligand| match ligand {
        Ligand::Hydrogen => LigandColor::Hydrogen,
        Ligand::LonePair => LigandColor::LonePair,
        Ligand::Vacancy => LigandColor::Vacancy,
        Ligand::Atom(node) => LigandColor::Atom(refined[node.index()]),
    };
    let atom_tied: Vec<usize> = molecule
        .stereo
        .iter()
        .enumerate()
        .filter(|(_, s)| {
            let sites: Vec<_> = s.ligands.iter().copied().map(ligand_color).collect();
            (0..sites.len()).any(|i| (i + 1..sites.len()).any(|j| sites[i] == sites[j]))
        })
        .map(|(i, _)| i)
        .collect();
    let allene_tied: Vec<usize> = molecule
        .allene_stereo
        .iter()
        .enumerate()
        .filter(|(_, s)| {
            let sites = s.ligands.map(ligand_color);
            color_parity(&sites[..2]).is_none() || color_parity(&sites[2..]).is_none()
        })
        .map(|(i, _)| i)
        .collect();
    let bond_tied: Vec<usize> = molecule
        .bond_stereo
        .iter()
        .enumerate()
        .filter(|(_, s)| {
            search.double_bond_flip(s, 0, refined).is_none()
                || search.double_bond_flip(s, 1, refined).is_none()
        })
        .map(|(i, _)| i)
        .collect();
    (atom_tied, allene_tied, bond_tied)
}

fn atom_stereo_is_fake(
    molecule: &Molecule,
    index: usize,
    baseline: &str,
    remaining: &mut usize,
    double_bond_stereo: DoubleBondStereoPolicy,
) -> Result<bool, CanonicalError> {
    let config = molecule.stereo[index].configuration;
    match config {
        Configuration::Tetrahedral(tag) => {
            let mut flipped = molecule.clone();
            flipped.stereo[index].configuration = Configuration::Tetrahedral(3 - tag);
            Ok(canonical_component_solving(&flipped, remaining, double_bond_stereo)? == baseline)
        }
        Configuration::SquarePlanar(_)
        | Configuration::TrigonalBipyramidal(_)
        | Configuration::Octahedral(_) => {
            for variant in 1..=config.coordination_variant_count() {
                let Some(next) = config.with_coordination_number(variant) else {
                    continue;
                };
                if next == config {
                    continue;
                }
                let mut flipped = molecule.clone();
                flipped.stereo[index].configuration = next;
                if canonical_component_solving(&flipped, remaining, double_bond_stereo)? != baseline
                {
                    return Ok(false);
                }
            }
            Ok(true)
        }
    }
}

fn canonical_component_tolerating_double_bonds(
    molecule: &Molecule,
    remaining: &mut usize,
) -> Result<String, CanonicalError> {
    const MAX_EXACT_COMBINATIONS: usize = 128;

    let viable: Vec<DoubleBondStereo> = molecule
        .bond_stereo
        .iter()
        .filter(|stereo| {
            let mut solo = molecule.clone();
            solo.bond_stereo = vec![(*stereo).clone()];
            write_smiles(&solo).is_ok()
        })
        .cloned()
        .collect();
    let mut reduced = molecule.clone();
    if viable.len() < molecule.bond_stereo.len() {
        reduced.bond_stereo = viable;
        match canonical_component(&reduced, remaining) {
            Ok(smiles) => return Ok(smiles),
            Err(CanonicalError::Write(WriteError::DoubleBondStereo)) => {}
            Err(error) => return Err(error),
        }
    }

    let m = reduced.bond_stereo.len();
    for drops in 1..=m {
        if !combinations_within(m, drops, MAX_EXACT_COMBINATIONS) {
            break;
        }
        if let Some(best) = best_dropping(&reduced, drops, remaining)? {
            return Ok(best);
        }
    }
    greedy_double_bond_subset(&reduced, remaining)
}

fn combinations_within(n: usize, k: usize, cap: usize) -> bool {
    let k = k.min(n.saturating_sub(k));
    let mut count = 1usize;
    for i in 0..k {
        count = count.saturating_mul(n - i) / (i + 1);
        if count > cap {
            return false;
        }
    }
    count <= cap
}

fn best_dropping(
    molecule: &Molecule,
    drops: usize,
    remaining: &mut usize,
) -> Result<Option<String>, CanonicalError> {
    let mut dropped = Vec::with_capacity(drops);
    let mut best: Option<String> = None;
    search_drops(molecule, drops, 0, &mut dropped, &mut best, remaining)?;
    Ok(best)
}

fn search_drops(
    molecule: &Molecule,
    drops: usize,
    start: usize,
    dropped: &mut Vec<usize>,
    best: &mut Option<String>,
    remaining: &mut usize,
) -> Result<(), CanonicalError> {
    let m = molecule.bond_stereo.len();
    if dropped.len() == drops {
        let mut candidate = molecule.clone();
        candidate.bond_stereo = molecule
            .bond_stereo
            .iter()
            .enumerate()
            .filter(|(i, _)| !dropped.contains(i))
            .map(|(_, s)| s.clone())
            .collect();
        match canonical_component(&candidate, remaining) {
            Ok(smiles) => {
                if best.as_ref().is_none_or(|prev| smiles < *prev) {
                    *best = Some(smiles);
                }
            }
            Err(CanonicalError::Write(WriteError::DoubleBondStereo)) => {}
            Err(error) => return Err(error),
        }
        return Ok(());
    }
    let needed = drops - dropped.len();
    for i in start..=m - needed {
        dropped.push(i);
        search_drops(molecule, drops, i + 1, dropped, best, remaining)?;
        dropped.pop();
    }
    Ok(())
}

fn greedy_double_bond_subset(
    molecule: &Molecule,
    remaining: &mut usize,
) -> Result<String, CanonicalError> {
    let mut search = Search::new(molecule, *remaining);
    let (mut colors, distinct) = search.initial_colors();
    search.refine(&mut colors, distinct);
    let refined = colors;
    let mut order: Vec<usize> = (0..molecule.bond_stereo.len()).collect();
    order.sort_by_key(|&i| {
        let s = &molecule.bond_stereo[i];
        let c0 = refined[s.ends[0].index()];
        let c1 = refined[s.ends[1].index()];
        let r0 = refined[s.references[0].index()];
        let r1 = refined[s.references[1].index()];
        if c0 <= c1 {
            (c0, c1, s.same_side, r0, r1)
        } else {
            (c1, c0, s.same_side, r1, r0)
        }
    });

    let mut kept = Vec::with_capacity(order.len());
    for i in order {
        kept.push(molecule.bond_stereo[i].clone());
        let mut test = molecule.clone();
        test.bond_stereo.clone_from(&kept);
        if write_smiles_ranked(&test, &refined).is_err() {
            kept.pop();
        }
    }
    loop {
        let mut test = molecule.clone();
        test.bond_stereo.clone_from(&kept);
        match canonical_component(&test, remaining) {
            Ok(smiles) => return Ok(smiles),
            Err(CanonicalError::Write(WriteError::DoubleBondStereo)) if !kept.is_empty() => {
                kept.pop();
            }
            Err(error) => return Err(error),
        }
    }
}

struct Symmetry {
    map: Vec<usize>,
    support: Vec<usize>,
}

impl Symmetry {
    fn from_map(map: Vec<usize>) -> Self {
        let support = map
            .iter()
            .enumerate()
            .filter_map(|(i, &target)| (target != i).then_some(i))
            .collect();
        Self { map, support }
    }

    fn preserves(&self, colors: &[usize]) -> bool {
        self.support
            .iter()
            .all(|&idx| colors[idx] == colors[self.map[idx]])
    }
}

struct Search<'a> {
    molecule: &'a Molecule,
    nodes: Vec<NodeIndex>,
    node_pos: Vec<usize>,
    adj_offsets: Vec<usize>,
    adj_edges: Vec<(u8, usize)>,
    has_stereo: bool,
    remaining: usize,
    leaves: BTreeMap<String, Vec<NodeIndex>>,
    symmetries: Vec<Symmetry>,
    // Reusable scratch buffers for refine() and covered()
    color_counts: Vec<usize>,
    cell_starts: Vec<usize>,
    order: Vec<usize>,
    next_colors: Vec<usize>,
    nbr_buf: Vec<u64>,
    nbr_ranges: Vec<(usize, usize)>,
    stereo_buf: Vec<Vec<StereoColor>>,
    stereo_dirty: Vec<usize>,
    covered_mark: Vec<u32>,
    covered_epoch: u32,
    covered_stack: Vec<usize>,
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum LigandColor {
    Hydrogen,
    LonePair,
    Vacancy,
    Atom(usize),
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
enum StereoColor {
    Tetrahedral(Option<bool>),
    Coordination(Option<CoordinationGeometry<LigandColor>>),
    Allene(Option<bool>),
    DoubleBond(usize, Option<bool>),
}

// A tie cannot supply a parity until later refinement distinguishes the ligands.
fn color_parity(sites: &[LigandColor]) -> Option<bool> {
    let mut odd = false;
    for (i, site) in sites.iter().enumerate() {
        for other in &sites[..i] {
            if other == site {
                return None;
            }
            odd ^= other > site;
        }
    }
    Some(odd)
}

fn double_bond_flip(
    node_pos: &[usize],
    adj_offsets: &[usize],
    adj_edges: &[(u8, usize)],
    stereo: &DoubleBondStereo,
    side: usize,
    colors: &[usize],
) -> Option<bool> {
    let pos = node_pos[stereo.ends[side].index()];
    let reference = stereo.references[side].index();
    let partner = stereo.ends[1 - side].index();
    let slice = &adj_edges[adj_offsets[pos]..adj_offsets[pos + 1]];
    let mut other = None;
    for &(_, nbr) in slice {
        if nbr != reference && nbr != partner {
            if other.is_some() {
                return None;
            }
            other = Some(nbr);
        }
    }
    let Some(other) = other else {
        return Some(false);
    };
    let (rc, oc) = (colors[reference], colors[other]);
    (rc != oc).then_some(rc > oc)
}

impl<'a> Search<'a> {
    fn new(molecule: &'a Molecule, remaining: usize) -> Self {
        let nodes: Vec<NodeIndex> = molecule.graph.node_indices().collect();
        let n = nodes.len();
        let bound = molecule.graph.node_bound();
        let mut node_pos = vec![0usize; bound];
        for (pos, &node) in nodes.iter().enumerate() {
            node_pos[node.index()] = pos;
        }
        let mut adj_offsets = Vec::with_capacity(n + 1);
        let mut adj_edges = Vec::with_capacity(molecule.graph.edge_count() * 2);
        adj_offsets.push(0);
        for &node in &nodes {
            for edge in molecule.graph.edges(node) {
                let other = if edge.source() == node {
                    edge.target()
                } else {
                    edge.source()
                };
                let code = bond_code(edge.weight().order, edge.weight().donor == Some(node));
                adj_edges.push((code, other.index()));
            }
            adj_offsets.push(adj_edges.len());
        }
        let has_stereo = !molecule.stereo.is_empty()
            || !molecule.allene_stereo.is_empty()
            || !molecule.bond_stereo.is_empty();
        Self {
            molecule,
            nodes,
            node_pos,
            adj_offsets,
            adj_edges,
            has_stereo,
            remaining,
            leaves: BTreeMap::new(),
            symmetries: Vec::new(),
            color_counts: vec![0; n + 1],
            cell_starts: vec![0; n + 1],
            order: (0..n).collect(),
            next_colors: vec![0; bound],
            nbr_buf: Vec::new(),
            nbr_ranges: vec![(0, 0); n],
            stereo_buf: if has_stereo {
                vec![Vec::new(); bound]
            } else {
                Vec::new()
            },
            stereo_dirty: Vec::new(),
            covered_mark: vec![0; bound],
            covered_epoch: 0,
            covered_stack: Vec::new(),
        }
    }

    fn initial_colors(&mut self) -> (Vec<usize>, usize) {
        let n = self.nodes.len();
        let bound = self.molecule.graph.node_bound();
        let mut bond_codes = Vec::with_capacity(self.adj_edges.len());
        let mut ranges = Vec::with_capacity(n);
        for pos in 0..n {
            let start = bond_codes.len();
            for &(code, _) in &self.adj_edges[self.adj_offsets[pos]..self.adj_offsets[pos + 1]] {
                bond_codes.push(code);
            }
            bond_codes[start..].sort_unstable();
            ranges.push((start, bond_codes.len()));
        }
        let atom_key = |pos: usize| {
            let atom = &self.molecule.graph[self.nodes[pos]];
            let (s, e) = ranges[pos];
            (
                atom.element.atomic_number,
                atom.isotope,
                atom.charge,
                u16::from(atom.hs) + u16::from(atom.implicit_hs),
                atom.aromatic,
                atom.atom_class,
                &bond_codes[s..e],
            )
        };
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_unstable_by(|&a, &b| atom_key(a).cmp(&atom_key(b)));
        let mut colors = vec![0usize; bound];
        let mut color = 0usize;
        for idx in 0..n {
            let pos = order[idx];
            if idx > 0 && atom_key(pos) != atom_key(order[idx - 1]) {
                color += 1;
            }
            colors[self.nodes[pos].index()] = color;
        }
        let distinct = if n == 0 { 0 } else { color + 1 };
        (colors, distinct)
    }

    fn seed_twin_symmetries(&mut self, colors: &[usize]) {
        let n = self.nodes.len();
        let bound = self.molecule.graph.node_bound();
        let stereo_counts = if self.has_stereo {
            stereo_reference_counts(self.molecule)
        } else {
            Vec::new()
        };
        let mut candidates: Vec<((usize, usize, u8), usize)> = Vec::new();
        for pos in 0..n {
            let start = self.adj_offsets[pos];
            let end = self.adj_offsets[pos + 1];
            if end - start != 1 {
                continue;
            }
            let node_idx = self.nodes[pos].index();
            if self.has_stereo && stereo_counts[node_idx] > 0 {
                continue;
            }
            let (code, parent) = self.adj_edges[start];
            candidates.push(((colors[node_idx], parent, code), node_idx));
        }
        if candidates.len() < 2 {
            return;
        }
        candidates.sort_unstable();
        for pair in candidates.windows(2) {
            if pair[0].0 == pair[1].0 {
                let (a, b) = (pair[0].1, pair[1].1);
                let mut map: Vec<usize> = (0..bound).collect();
                map.swap(a, b);
                self.symmetries.push(Symmetry::from_map(map));
            }
        }
    }

    fn fill_stereo_colors(&mut self, colors: &[usize]) {
        for idx in self.stereo_dirty.drain(..) {
            self.stereo_buf[idx].clear();
        }
        let molecule = self.molecule;
        if !self.has_stereo {
            return;
        }
        let ligand_color = |ligand| match ligand {
            Ligand::Hydrogen => LigandColor::Hydrogen,
            Ligand::LonePair => LigandColor::LonePair,
            Ligand::Vacancy => LigandColor::Vacancy,
            Ligand::Atom(node) => LigandColor::Atom(colors[node.index()]),
        };
        let stereo_buf = &mut self.stereo_buf;
        let stereo_dirty = &mut self.stereo_dirty;
        let mut push = |idx: usize, entry: StereoColor| {
            let slot = &mut stereo_buf[idx];
            if slot.is_empty() {
                stereo_dirty.push(idx);
            }
            slot.push(entry);
        };
        for stereo in &molecule.stereo {
            let sites: Vec<_> = stereo.ligands.iter().copied().map(ligand_color).collect();
            let entry = match stereo.configuration {
                Configuration::Tetrahedral(n) => {
                    let parity = color_parity(&sites).map(|odd| odd ^ (n == 1));
                    StereoColor::Tetrahedral(parity)
                }
                Configuration::SquarePlanar(_)
                | Configuration::TrigonalBipyramidal(_)
                | Configuration::Octahedral(_) => {
                    StereoColor::Coordination(stereo.configuration.coordination_geometry(&sites))
                }
            };
            push(stereo.center.index(), entry);
        }
        for stereo in &molecule.allene_stereo {
            let sites = stereo.ligands.map(ligand_color);
            // Reversing the axis exchanges two pairs, which leaves parity unchanged.
            let parity = color_parity(&sites[..2])
                .zip(color_parity(&sites[2..]))
                .map(|(a, b)| a ^ b ^ (stereo.permutation == 1));
            push(stereo.center.index(), StereoColor::Allene(parity));
        }
        for stereo in &molecule.bond_stereo {
            let parity = double_bond_flip(
                &self.node_pos,
                &self.adj_offsets,
                &self.adj_edges,
                stereo,
                0,
                colors,
            )
            .zip(double_bond_flip(
                &self.node_pos,
                &self.adj_offsets,
                &self.adj_edges,
                stereo,
                1,
                colors,
            ))
            .map(|(a, b)| a ^ b ^ stereo.same_side);
            for side in [0, 1] {
                push(
                    stereo.ends[side].index(),
                    StereoColor::DoubleBond(colors[stereo.ends[1 - side].index()], parity),
                );
            }
        }
        for &idx in &self.stereo_dirty {
            if self.stereo_buf[idx].len() > 1 {
                self.stereo_buf[idx].sort_unstable();
            }
        }
    }

    fn double_bond_flip(
        &self,
        stereo: &DoubleBondStereo,
        side: usize,
        colors: &[usize],
    ) -> Option<bool> {
        double_bond_flip(
            &self.node_pos,
            &self.adj_offsets,
            &self.adj_edges,
            stereo,
            side,
            colors,
        )
    }

    fn refine(&mut self, colors: &mut [usize], mut distinct: usize) -> usize {
        let n = self.nodes.len();
        loop {
            if distinct >= n {
                return distinct;
            }
            self.color_counts[..distinct].fill(0);
            for &node in &self.nodes {
                self.color_counts[colors[node.index()]] += 1;
            }

            if self.has_stereo {
                self.fill_stereo_colors(colors);
            }
            self.nbr_buf.clear();
            for pos in 0..n {
                let c = colors[self.nodes[pos].index()];
                if self.color_counts[c] <= 1 {
                    self.nbr_ranges[pos] = (0, 0);
                    continue;
                }
                let start = self.nbr_buf.len();
                for &(code, other) in
                    &self.adj_edges[self.adj_offsets[pos]..self.adj_offsets[pos + 1]]
                {
                    let key = (u64::from(code) << 56) | (colors[other] as u64);
                    self.nbr_buf.push(key);
                }
                self.nbr_buf[start..].sort_unstable();
                self.nbr_ranges[pos] = (start, self.nbr_buf.len());
            }

            // Bucket vertices by current color in O(n) and sort only non-singleton cells.
            let mut offset = 0usize;
            for c in 0..distinct {
                self.cell_starts[c] = offset;
                offset += self.color_counts[c];
            }
            for pos in 0..n {
                let c = colors[self.nodes[pos].index()];
                let slot = self.cell_starts[c];
                self.cell_starts[c] = slot + 1;
                self.order[slot] = pos;
            }

            let mut start = 0usize;
            for c in 0..distinct {
                let count = self.color_counts[c];
                let end = start + count;
                if count > 1 {
                    let nbr_buf = &self.nbr_buf;
                    let nbr_ranges = &self.nbr_ranges;
                    let nodes = &self.nodes;
                    let stereo_buf = &self.stereo_buf;
                    let has_stereo = self.has_stereo;
                    self.order[start..end].sort_unstable_by(|&a, &b| {
                        let (sa, ea) = nbr_ranges[a];
                        let (sb, eb) = nbr_ranges[b];
                        let ord = nbr_buf[sa..ea].cmp(&nbr_buf[sb..eb]);
                        if ord != std::cmp::Ordering::Equal || !has_stereo {
                            return ord;
                        }
                        stereo_buf[nodes[a].index()].cmp(&stereo_buf[nodes[b].index()])
                    });
                }
                start = end;
            }

            let mut next_color = 0usize;
            let first_pos = self.order[0];
            self.next_colors[self.nodes[first_pos].index()] = 0;
            for idx in 1..n {
                let prev = self.order[idx - 1];
                let cur = self.order[idx];
                let prev_idx = self.nodes[prev].index();
                let cur_idx = self.nodes[cur].index();
                let differs = colors[prev_idx] != colors[cur_idx]
                    || self.nbr_buf[self.nbr_ranges[prev].0..self.nbr_ranges[prev].1]
                        != self.nbr_buf[self.nbr_ranges[cur].0..self.nbr_ranges[cur].1]
                    || (self.has_stereo && self.stereo_buf[prev_idx] != self.stereo_buf[cur_idx]);
                if differs {
                    next_color += 1;
                }
                self.next_colors[cur_idx] = next_color;
            }
            let next_distinct = next_color + 1;
            for &node in &self.nodes {
                colors[node.index()] = self.next_colors[node.index()];
            }
            if next_distinct == distinct {
                return next_distinct;
            }
            distinct = next_distinct;
        }
    }

    // Only automorphisms preserving this partition can prune its children.
    fn covered(&mut self, node: NodeIndex, explored: &[NodeIndex], preserving: &[usize]) -> bool {
        self.covered_epoch = self.covered_epoch.wrapping_add(1);
        if self.covered_epoch == 0 {
            self.covered_mark.fill(0);
            self.covered_epoch = 1;
        }
        let epoch = self.covered_epoch;
        self.covered_stack.clear();
        self.covered_stack.push(node.index());
        self.covered_mark[node.index()] = epoch;
        while let Some(n) = self.covered_stack.pop() {
            if explored.iter().any(|e| e.index() == n) {
                return true;
            }
            for &gen_idx in preserving {
                let next = self.symmetries[gen_idx].map[n];
                if self.covered_mark[next] != epoch {
                    self.covered_mark[next] = epoch;
                    self.covered_stack.push(next);
                }
            }
        }
        false
    }

    fn visit(&mut self, mut colors: Vec<usize>, distinct: usize) -> Result<String, CanonicalError> {
        self.remaining = self
            .remaining
            .checked_sub(1)
            .ok_or(CanonicalError::SearchLimit)?;
        let distinct = self.refine(&mut colors, distinct);
        let n = self.nodes.len();
        if distinct >= n {
            let (smiles, order) = write(self.molecule, &colors).map_err(CanonicalError::Write)?;
            if let Some(previous) = self.leaves.get(&smiles) {
                let mut map: Vec<usize> = (0..colors.len()).collect();
                for (a, b) in previous.iter().zip(&order) {
                    map[a.index()] = b.index();
                }
                self.symmetries.push(Symmetry::from_map(map));
            } else {
                self.leaves.insert(smiles.clone(), order);
            }
            return Ok(smiles);
        }

        self.color_counts[..distinct].fill(0);
        for &node in &self.nodes {
            self.color_counts[colors[node.index()]] += 1;
        }
        let Some((target_color, &cell_len)) = self.color_counts[..distinct]
            .iter()
            .enumerate()
            .filter(|&(_, &count)| count > 1)
            .min_by_key(|&(_, &count)| count)
        else {
            unreachable!("distinct < n implies at least one cell with len > 1");
        };
        let mut cell = Vec::with_capacity(cell_len);
        for &node in &self.nodes {
            if colors[node.index()] == target_color {
                cell.push(node);
            }
        }

        let mut best: Option<String> = None;
        let mut write_error = WriteError::UnsupportedStereo;
        let mut explored = Vec::new();
        let mut checked_symmetries = 0usize;
        let mut preserving = Vec::new();
        for node in cell {
            for (index, symmetry) in self.symmetries.iter().enumerate().skip(checked_symmetries) {
                if symmetry.preserves(&colors) {
                    preserving.push(index);
                }
            }
            checked_symmetries = self.symmetries.len();
            if self.covered(node, &explored, &preserving) {
                continue;
            }
            explored.push(node);
            let mut child = colors.clone();
            for &other in &self.nodes {
                let c = &mut child[other.index()];
                if *c > target_color {
                    *c += 1;
                }
            }
            child[node.index()] += 1;
            let result = match self.visit(child, distinct + 1) {
                Ok(result) => result,
                Err(CanonicalError::Write(error)) => {
                    write_error = error;
                    continue;
                }
                Err(error) => return Err(error),
            };
            if best.as_ref().is_none_or(|previous| result < *previous) {
                best = Some(result);
            }
        }
        best.ok_or(CanonicalError::Write(write_error))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::smiles::reader::molecule::parse_molecule;
    fn canonical(s: &str) -> String {
        canonical_smiles(&parse_molecule(s).unwrap().1).unwrap()
    }
    #[test]
    fn spelling_and_idempotence() {
        for group in [
            vec!["CCO", "OCC", "C(O)C"],
            vec!["CC/C=C/C", "C/C=C/CC"],
            vec!["CC(C)=CCC/C(C)=C/CO", "OC/C=C(\\C)CCC=C(C)C"],
            vec!["c1ccccc1", "C1=CC=CC=C1", "c2ccccc2"],
            vec!["CC(C)O", "OC(C)C", "C(C)(O)C"],
            vec!["N->[Cu]", "[Cu]<-N"],
            vec!["CC.O", "O.CC"],
            vec![
                "[Cu]<-N.F/C=C/F.FC=[C@]=CCl.N[C@H](F)C",
                "N[C@H](F)C.FC=[C@]=CCl.F/C=C/F.N->[Cu]",
            ],
            vec!["C1CC2CCC1C2", "C1C2CCC(C1)C2"],
        ] {
            let expected = canonical(group[0]);
            for input in group {
                assert_eq!(canonical(input), expected, "{input}");
            }
            assert_eq!(canonical(&expected), expected);
        }
    }
    #[test]
    fn stereo_is_preserved() {
        for (a, b) in [
            ("N[C@H](F)C", "N[C@@H](F)C"),
            ("F/C=C/F", "F/C=C\\F"),
            ("FC=[C@]=CCl", "FC=[C@@]=CCl"),
        ] {
            assert_ne!(canonical(a), canonical(b));
            assert_eq!(canonical(&canonical(a)), canonical(a));
        }
    }

    #[test]
    fn canonicalizing_preserves_unspecified_double_bond_stereo() {
        // SureChEMBL SCHEMBL12328173 has a required slash next to an unspecified double bond.
        let input = r"CC1=C/C(=C(/C#N)C(=O)O)C=C(/C=C/c2ccc(N(c3ccccc3)c3ccc(-c4ccc(N(c5ccccc5)c5ccc(/C=C/C6=C/C(=C(\C#N)C(=O)O)C=C(C)O6)cc5)cc4)cc3)cc2)O1";
        let molecule = parse_molecule(input).unwrap().1;
        let output = canonical_smiles(&molecule).unwrap();
        let reread = parse_molecule(&output).unwrap().1;
        assert_eq!(reread.double_bond_stereo().len(), 5);
        assert_eq!(molecule.exact_match(&reread), Ok(true));
        assert_eq!(canonical_smiles(&reread), Ok(output));
    }

    #[test]
    fn canonical_round_trips_preserve_the_molecule() {
        for input in include_str!("../../../tests/data/surechembl_canonical_round_trip.smi").lines()
        {
            let molecule = parse_molecule(input)
                .unwrap_or_else(|e| panic!("{input}: {e:?}"))
                .1;
            let output = canonical_smiles(&molecule).unwrap_or_else(|e| panic!("{input}: {e:?}"));
            let reread = parse_molecule(&output)
                .unwrap_or_else(|e| panic!("{input} -> {output}: {e:?}"))
                .1;
            assert_eq!(
                molecule.exact_match(&reread),
                Ok(true),
                "{input} -> {output}"
            );
            assert_eq!(canonical_smiles(&reread), Ok(output), "{input}");
        }
    }

    #[test]
    fn canonicalizes_surechembl_double_bond_stereo() {
        let records = include_str!("../../../tests/data/surechembl_double_bond_stereo.smi");
        for record in records.lines() {
            let molecule = parse_molecule(record).unwrap().1;
            let output =
                canonical_smiles(&molecule).unwrap_or_else(|error| panic!("{record}: {error:?}"));
            let reread = parse_molecule(&output).unwrap().1;
            assert_eq!(molecule.exact_match(&reread), Ok(true), "{record}");
            assert_eq!(canonical(&output), output, "{record}");

            let ranks: Vec<_> = (0..molecule.count_atoms()).rev().collect();
            let reversed = super::super::writer::write_smiles_ranked(&molecule, &ranks).unwrap();
            assert_eq!(canonical(&reversed), output, "{record}");
        }
    }

    #[test]
    fn opposite_configurations_split_otherwise_identical_branches() {
        for input in [
            "F[C@H](Cl)[C@H](Br)[C@@H](F)Cl",
            "FC=[C@]=C(Cl)OC(Cl)=[C@@]=CF",
            "F/C=C/O/C=C\\F",
        ] {
            let molecule = parse_molecule(input).unwrap().1;
            let expected = canonical_smiles_with_limit(&molecule, 1)
                .unwrap_or_else(|error| panic!("{input}: {error:?}"));
            let ranks: Vec<_> = (0..molecule.count_atoms()).rev().collect();
            let reversed = super::super::writer::write_smiles_ranked(&molecule, &ranks).unwrap();
            assert_eq!(canonical(&reversed), expected, "{input} -> {reversed}");
            assert_eq!(canonical(&expected), expected, "{input}");
        }
    }

    #[test]
    fn shuffled_writer_orders() {
        for input in [
            "CC(C)(C)C",
            "c1ccccc1",
            "C12C3C4C1C5C2C3C45",
            "N[C@H](F)C",
            "CC/C=C/C",
            "CC(C)=CCC/C(C)=C/CO",
            "F/C=C\\F",
            "FC=[C@]=CCl",
            "CC.O.CC",
            "Cl[13CH2]C(=O)[O-]",
            "Cc1cc2cc(=O)cc3oc(C)cc(o1)c23",
        ] {
            let molecule = parse_molecule(input).unwrap().1;
            let expected = canonical_smiles(&molecule).unwrap();
            let mut state = 12345_u64;
            for _ in 0..10 {
                let mut ranks: Vec<_> = (0..molecule.graph.node_bound()).collect();
                for i in (1..ranks.len()).rev() {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    let j = usize::try_from(state % (i as u64 + 1)).unwrap();
                    ranks.swap(i, j);
                }
                let written = super::super::writer::write_smiles_ranked(&molecule, &ranks).unwrap();
                assert_eq!(canonical(&written), expected, "{input} -> {written}");
            }
        }
    }

    #[test]
    fn symmetric_arms_fit_a_small_budget() {
        let molecule = parse_molecule("C(C(C)(C)C)(C(C)(C)C)(C(C)(C)C)C(C)(C)C")
            .unwrap()
            .1;
        assert!(canonical_smiles_with_limit(&molecule, 1000).is_ok());
    }

    #[test]
    fn repository_stereo_records_have_consistent_support() {
        for line in include_str!("../../../tests/data/stereo.smi").lines() {
            let (_, molecule) = parse_molecule(line).unwrap();
            let output = canonical_smiles(&molecule).unwrap_or_else(|e| panic!("{line}: {e:?}"));
            let (_, reread) = parse_molecule(&output).unwrap();
            assert_eq!(
                molecule.exact_match(&reread),
                Ok(true),
                "{line} -> {output}"
            );
            assert_eq!(canonical(&output), output, "{line}");
        }
    }

    #[test]
    fn coordination_stereo_is_canonicalized_in_every_order() {
        for group in [
            vec![
                "S[As@TB1](F)(Cl)(Br)N",
                "S[As@TB2](Br)(Cl)(F)N",
                "S[As@TB5](F)(N)(Cl)Br",
                "F[As@TB10](S)(Cl)(N)Br",
                "F[As@TB15](Cl)(S)(Br)N",
                "Br[As@TB20](Cl)(S)(F)N",
            ],
            vec![
                "C[Co@](F)(Cl)(Br)(I)S",
                "F[Co@@](S)(I)(C)(Cl)Br",
                "S[Co@OH5](F)(I)(Cl)(C)Br",
                "Br[Co@OH9](C)(S)(Cl)(F)I",
                "Br[Co@OH12](Cl)(I)(F)(S)C",
                "Cl[Co@OH15](C)(Br)(F)(I)S",
                "Cl[Co@OH19](C)(I)(F)(S)Br",
                "I[Co@OH27](Cl)(Br)(F)(S)C",
            ],
            vec![
                "Cl[Pt@SP1](Cl)(<-[NH3])<-[NH3]",
                "[NH3]->[Pt@SP1](Cl)(Cl)<-[NH3]",
            ],
            vec![
                "Cl[Pt@SP2](Cl)(<-[NH3])<-[NH3]",
                "[NH3]->[Pt@SP3](Cl)(Cl)<-[NH3]",
            ],
        ] {
            let expected = canonical(group[0]);
            for input in group {
                assert_eq!(canonical(input), expected, "{input}");
            }
            assert_eq!(canonical(&expected), expected);
        }
        for (a, b) in [
            ("F[Pt@SP1](Cl)(Br)I", "F[Pt@SP2](Cl)(Br)I"),
            ("S[As@TB1](F)(Cl)(Br)N", "S[As@TB2](F)(Cl)(Br)N"),
            ("C[Co@](F)(Cl)(Br)(I)S", "C[Co@@](F)(Cl)(Br)(I)S"),
        ] {
            assert_ne!(canonical(a), canonical(b), "{a} vs {b}");
        }
    }

    #[test]
    fn search_budget_is_an_error() {
        for input in ["C1CCCCC1", "O.O"] {
            let molecule = parse_molecule(input).unwrap().1;
            assert_eq!(
                canonical_smiles_with_limit(&molecule, 1),
                Err(CanonicalError::SearchLimit)
            );
        }
    }

    #[test]
    fn tolerate_double_bond_stereo_drops_minimally_and_is_idempotent() {
        use crate::io::smiles::reader::molecule::{parse_molecule_with, DirectionalStereo};

        let (_, mol) = parse_molecule_with(
            "C/C=C\\C(/C=C\\C)=C(/C=C\\C)/C(O)=C\\C",
            DirectionalStereo::Tolerate,
        )
        .unwrap();
        assert_eq!(
            canonical_smiles_with(&mol, DEFAULT_SEARCH_LIMIT, DoubleBondStereoPolicy::Reject),
            Err(CanonicalError::Write(WriteError::DoubleBondStereo))
        );
        let tolerated =
            canonical_smiles_with(&mol, DEFAULT_SEARCH_LIMIT, DoubleBondStereoPolicy::Tolerate)
                .unwrap();
        assert_eq!(
            mol.double_bond_stereo().len() - 1,
            parse_molecule(&tolerated)
                .unwrap()
                .1
                .double_bond_stereo()
                .len(),
            "expected minimal drop of 1 bond: {tolerated}"
        );
        assert_eq!(canonical(&tolerated), tolerated);

        for seed in [12345_u64, 67890_u64, 13579_u64] {
            let mut state = seed;
            let mut perm: Vec<NodeIndex> = mol.graph.node_indices().collect();
            for i in (1..perm.len()).rev() {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                let j = usize::try_from(state % (i as u64 + 1)).unwrap();
                perm.swap(i, j);
            }
            let reordered = components(&mol, &[perm]).pop().unwrap();
            assert_eq!(
                canonical_smiles_with(
                    &reordered,
                    DEFAULT_SEARCH_LIMIT,
                    DoubleBondStereoPolicy::Tolerate,
                )
                .unwrap(),
                tolerated,
                "seed {seed}"
            );
        }

        let (_, disconnected) = parse_molecule_with(
            "O.C/C=C\\C(/C=C\\C)=C(/C=C\\C)/C(O)=C\\C",
            DirectionalStereo::Tolerate,
        )
        .unwrap();
        assert_eq!(
            canonical_smiles_with(
                &disconnected,
                DEFAULT_SEARCH_LIMIT,
                DoubleBondStereoPolicy::Tolerate,
            )
            .unwrap(),
            format!("{tolerated}.O")
        );
    }

    #[test]
    fn explicit_hydrogens_canonicalize_to_h_counts() {
        for (a, b) in [
            ("C", "[H]C([H])([H])[H]"),
            ("CO", "[H]OC"),
            ("[C@H](F)(Cl)Br", "[C@]([H])(F)(Cl)Br"),
            ("F[C@H](Cl)Br", "F[C@]([H])(Cl)Br"),
            ("[nH]1cccc1", "[n]([H])1cccc1"),
            ("F/C=C\\F", "[H]/C(F)=C/F"),
            ("OC=[C@]=CF", "OC([H])=[C@AL1]=C([H])F"),
            ("C[Pt@SP1H](Cl)F", "C[Pt@SP1]([H])(Cl)F"),
        ] {
            let out_a = canonical(a);
            let out_b = canonical(b);
            assert_eq!(out_a, out_b, "{a} vs {b}");
            assert_eq!(canonical(&out_a), out_a, "{a}");
        }
        assert_ne!(
            canonical("[2H][C@](F)(Cl)Br"),
            canonical("[H][C@](F)(Cl)Br")
        );
        assert_ne!(canonical("[2H]C"), canonical("[H]C"));
        assert_eq!(canonical("[H][H]"), "[H][H]");
    }

    #[test]
    fn fake_and_interacting_stereo_are_normalized() {
        for (a, b) in [
            ("C[C@H](C)O", "CC(C)O"),
            ("C[C@H](C)O", "C[C@@H](C)O"),
            ("[H][C@]([H])(F)Cl", "FCCl"),
            ("C/C=C(/C)C", "CC=C(C)C"),
            ("[H]/C([H])=C/F", "C=CF"),
            ("C[C@@H](O)[C@H](O)[C@@H](C)O", "C[C@@H](O)C(O)[C@@H](C)O"),
            (
                "C[C@@H](O)[C@H](O)[C@@H](C)O",
                "C[C@@H](O)[C@@H](O)[C@@H](C)O",
            ),
            ("C[C@H]1CCC(C)CC1", "CC1CCC(C)CC1"),
            (
                "F[C@H](C1CCC(C)CC1)C2CC[C@H](C)CC2",
                "FC(C1CCC(C)CC1)C2CCC(C)CC2",
            ),
            (
                "F[C@H]([NH2]->[Cu])[NH2]->[Cu]",
                "FC([NH2]->[Cu])[NH2]->[Cu]",
            ),
            ("CC(C)=[C@]=CF", "CC(C)=C=CF"),
            ("C=[C@]=C", "C=C=C"),
            ("[H]C([H])=[C@]=C([H])F", "C=C=CF"),
            ("F[Pt@SP1](C)(C)C", "F[Pt](C)(C)C"),
            ("F[Pt@SP1H3]", "F[PtH3]"),
        ] {
            let out_a = canonical(a);
            let out_b = canonical(b);
            assert_eq!(out_a, out_b, "{a} vs {b}");
            assert_eq!(canonical(&out_a), out_a, "{a}");
        }

        for (a, b) in [
            ("F[Pt@SP1H2]Cl", "F[Pt@SP3H2]Cl"),
            ("C[C@H]1CC[C@H](C)CC1", "C[C@H]1CC[C@@H](C)CC1"),
            ("C[C@H]1CC[C@H](C)CC1", "C[C@H]1CCC(C)CC1"),
            (
                "C[C@@H](O)[C@H](O)[C@H](C)O",
                "C[C@@H](O)[C@@H](O)[C@H](C)O",
            ),
            ("C[C@@H](O)[C@H](O)[C@H](C)O", "C[C@@H](O)C(O)[C@H](C)O"),
            (
                "F[C@H]([NH2]->[Cu])[NH2]<-[Cu]",
                "F[C@@H]([NH2]->[Cu])[NH2]<-[Cu]",
            ),
            ("F[C@H](Cl)Br", "F[C@@]([H])(Cl)Br"),
            ("C[Pt@SP1](Cl)(Cl)N", "C[Pt@SP3](Cl)(Cl)N"),
        ] {
            assert_ne!(canonical(a), canonical(b), "{a} vs {b}");
        }
    }
}
