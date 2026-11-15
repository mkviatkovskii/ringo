//! Deterministic individualization/refinement for canonical SMILES.
//!
//! Stereo participates in refinement once the current colors distinguish its ligands.
use super::writer::{write, WriteError};
use crate::model::bond::Bond;
use crate::model::molecule::Molecule;
use crate::model::stereo::{Configuration, DoubleBondStereo, Ligand};
use petgraph::algo::kosaraju_scc;
use petgraph::stable_graph::NodeIndex;
use petgraph::visit::{EdgeIndexable, EdgeRef, IntoEdgeReferences, NodeIndexable};
use std::collections::BTreeMap;

/// Maximum search nodes used by the default canonicalizer.
pub const DEFAULT_SEARCH_LIMIT: usize = 100_000;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CanonicalError {
    /// An exact answer was not established within the search budget.
    SearchLimit,
    Write(WriteError),
}

/// Canonicalizes the represented graph, preserving its stereo annotations.
///
/// This does not yet remove chemically impossible stereo annotations. The result
/// is Ringo's canonical spelling, not necessarily another toolkit's spelling.
pub fn canonical_smiles(molecule: &Molecule) -> Result<String, CanonicalError> {
    canonical_smiles_with_limit(molecule, DEFAULT_SEARCH_LIMIT)
}

/// Like [`canonical_smiles`], with an explicit budget of search nodes.
/// Exhausting the budget returns an error, never a provisional canonical string.
pub fn canonical_smiles_with_limit(
    molecule: &Molecule,
    limit: usize,
) -> Result<String, CanonicalError> {
    if molecule
        .stereo
        .iter()
        .any(|s| !matches!(s.configuration, Configuration::Tetrahedral(_)))
    {
        return Err(CanonicalError::Write(WriteError::UnsupportedStereo));
    }
    let mut remaining = limit;
    let groups = kosaraju_scc(&molecule.graph);
    if groups.len() <= 1 {
        return canonical_component(molecule, &mut remaining);
    }

    let mut output = components(molecule, &groups)
        .iter()
        .map(|component| canonical_component(component, &mut remaining))
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

fn canonical_component(
    molecule: &Molecule,
    remaining: &mut usize,
) -> Result<String, CanonicalError> {
    let nodes: Vec<_> = molecule.graph.node_indices().collect();
    let keys: Vec<_> = nodes
        .iter()
        .map(|&node| {
            let atom = &molecule.graph[node];
            let mut bonds: Vec<_> = molecule
                .graph
                .edges(node)
                .map(|e| (e.weight().order, e.weight().donor == Some(node)))
                .collect();
            bonds.sort_unstable();
            (
                atom.element.atomic_number,
                atom.isotope,
                atom.charge,
                u16::from(atom.hs) + u16::from(atom.implicit_hs),
                atom.aromatic,
                atom.atom_class,
                bonds,
            )
        })
        .collect();
    let colors = compress(&nodes, &keys, molecule.graph.node_bound());
    let mut search = Search {
        molecule,
        nodes,
        remaining: *remaining,
        leaves: BTreeMap::new(),
        symmetries: Vec::new(),
    };
    let result = search.visit(colors);
    *remaining = search.remaining;
    result
}

fn compress<K: Ord>(nodes: &[NodeIndex], keys: &[K], bound: usize) -> Vec<usize> {
    let mut order: Vec<_> = (0..nodes.len()).collect();
    order.sort_by_key(|&i| &keys[i]);
    let mut colors = vec![0; bound];
    let mut color = 0;
    for (position, &i) in order.iter().enumerate() {
        if position > 0 && keys[i] != keys[order[position - 1]] {
            color += 1;
        }
        colors[nodes[i].index()] = color;
    }
    colors
}

struct Search<'a> {
    molecule: &'a Molecule,
    nodes: Vec<NodeIndex>,
    remaining: usize,
    leaves: BTreeMap<String, Vec<NodeIndex>>,
    symmetries: Vec<Vec<usize>>,
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
enum LigandColor {
    Hydrogen,
    LonePair,
    Vacancy,
    Atom(usize),
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
enum StereoColor {
    Tetrahedral(Option<bool>),
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

impl Search<'_> {
    fn stereo_colors(&self, colors: &[usize]) -> Vec<Vec<StereoColor>> {
        let molecule = self.molecule;
        if molecule.stereo.is_empty()
            && molecule.allene_stereo.is_empty()
            && molecule.bond_stereo.is_empty()
        {
            return Vec::new();
        }
        let ligand_color = |ligand| match ligand {
            Ligand::Hydrogen => LigandColor::Hydrogen,
            Ligand::LonePair => LigandColor::LonePair,
            Ligand::Vacancy => LigandColor::Vacancy,
            Ligand::Atom(node) => LigandColor::Atom(colors[node.index()]),
        };
        let mut result: Vec<Vec<StereoColor>> = (0..colors.len()).map(|_| Vec::new()).collect();
        for stereo in &molecule.stereo {
            if let Configuration::Tetrahedral(n) = stereo.configuration {
                let sites: Vec<_> = stereo.ligands.iter().copied().map(ligand_color).collect();
                let parity = color_parity(&sites).map(|odd| odd ^ (n == 1));
                result[stereo.center.index()].push(StereoColor::Tetrahedral(parity));
            }
        }
        for stereo in &molecule.allene_stereo {
            let sites = stereo.ligands.map(ligand_color);
            // Reversing the axis exchanges two pairs, which leaves parity unchanged.
            let parity = color_parity(&sites[..2])
                .zip(color_parity(&sites[2..]))
                .map(|(a, b)| a ^ b ^ (stereo.permutation == 1));
            result[stereo.center.index()].push(StereoColor::Allene(parity));
        }
        for stereo in &molecule.bond_stereo {
            let parity = self
                .double_bond_flip(stereo, 0, colors)
                .zip(self.double_bond_flip(stereo, 1, colors))
                .map(|(a, b)| a ^ b ^ stereo.same_side);
            for side in [0, 1] {
                result[stereo.ends[side].index()].push(StereoColor::DoubleBond(
                    colors[stereo.ends[1 - side].index()],
                    parity,
                ));
            }
        }
        for stereo in &mut result {
            stereo.sort_unstable();
        }
        result
    }

    fn double_bond_flip(
        &self,
        stereo: &DoubleBondStereo,
        side: usize,
        colors: &[usize],
    ) -> Option<bool> {
        let reference = stereo.references[side];
        let mut alternatives = self
            .molecule
            .graph
            .neighbors(stereo.ends[side])
            .filter(|&node| node != reference && node != stereo.ends[1 - side]);
        let Some(other) = alternatives.next() else {
            return Some(false);
        };
        let (reference, other) = (colors[reference.index()], colors[other.index()]);
        (alternatives.next().is_none() && reference != other).then_some(reference > other)
    }

    fn refine(&self, mut colors: Vec<usize>) -> Vec<usize> {
        loop {
            let stereo = self.stereo_colors(&colors);
            let keys: Vec<_> = self
                .nodes
                .iter()
                .map(|&node| {
                    let mut neighbours: Vec<_> = self
                        .molecule
                        .graph
                        .edges(node)
                        .map(|e| {
                            let other = if e.source() == node {
                                e.target()
                            } else {
                                e.source()
                            };
                            (
                                e.weight().order,
                                e.weight().donor == Some(node),
                                colors[other.index()],
                            )
                        })
                        .collect();
                    neighbours.sort_unstable();
                    let stereo = stereo.get(node.index()).map_or(&[][..], Vec::as_slice);
                    (colors[node.index()], neighbours, stereo)
                })
                .collect();
            let next = compress(&self.nodes, &keys, colors.len());
            if next == colors {
                return colors;
            }
            colors = next;
        }
    }

    // Only automorphisms preserving this partition can prune its children.
    fn covered(
        node: NodeIndex,
        explored: &[NodeIndex],
        colors: &[usize],
        generators: &[&Vec<usize>],
    ) -> bool {
        let mut seen = vec![false; colors.len()];
        let mut stack = vec![node.index()];
        seen[node.index()] = true;
        while let Some(n) = stack.pop() {
            if explored.iter().any(|e| e.index() == n) {
                return true;
            }
            for generator in generators {
                let next = generator[n];
                if !std::mem::replace(&mut seen[next], true) {
                    stack.push(next);
                }
            }
        }
        false
    }

    fn visit(&mut self, colors: Vec<usize>) -> Result<String, CanonicalError> {
        self.remaining = self
            .remaining
            .checked_sub(1)
            .ok_or(CanonicalError::SearchLimit)?;
        let colors = self.refine(colors);
        let mut cells: BTreeMap<usize, Vec<NodeIndex>> = BTreeMap::new();
        for &node in &self.nodes {
            cells.entry(colors[node.index()]).or_default().push(node);
        }
        let Some(cell) = cells
            .values()
            .filter(|c| c.len() > 1)
            .min_by_key(|c| c.len())
        else {
            let (smiles, order) = write(self.molecule, &colors).map_err(CanonicalError::Write)?;
            if let Some(previous) = self.leaves.get(&smiles) {
                let mut symmetry: Vec<_> = (0..colors.len()).collect();
                for (a, b) in previous.iter().zip(&order) {
                    symmetry[a.index()] = b.index();
                }
                self.symmetries.push(symmetry);
            } else {
                self.leaves.insert(smiles.clone(), order);
            }
            return Ok(smiles);
        };
        let mut best: Option<String> = None;
        let mut write_error = WriteError::UnsupportedStereo;
        let mut explored = Vec::new();
        let mut checked_symmetries = 0;
        let mut preserving = Vec::new();
        for &node in cell {
            for index in checked_symmetries..self.symmetries.len() {
                if self
                    .nodes
                    .iter()
                    .all(|n| colors[n.index()] == colors[self.symmetries[index][n.index()]])
                {
                    preserving.push(index);
                }
            }
            checked_symmetries = self.symmetries.len();
            let generators: Vec<_> = preserving.iter().map(|&i| &self.symmetries[i]).collect();
            if Self::covered(node, &explored, &colors, &generators) {
                continue;
            }
            explored.push(node);
            let mut child = colors.clone();
            for &other in &self.nodes {
                child[other.index()] *= 2;
            }
            child[node.index()] += 1;
            let result = match self.visit(child) {
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
            if molecule
                .stereo
                .iter()
                .any(|s| !matches!(s.configuration, Configuration::Tetrahedral(_)))
            {
                assert_eq!(
                    canonical_smiles(&molecule),
                    Err(CanonicalError::Write(WriteError::UnsupportedStereo)),
                    "{line}"
                );
            } else {
                let output =
                    canonical_smiles(&molecule).unwrap_or_else(|e| panic!("{line}: {e:?}"));
                assert_eq!(canonical(&output), output, "{line}");
            }
        }
    }

    #[test]
    fn unsupported_coordination_is_rejected_in_every_order() {
        for input in [
            "F[Pt@SP1](Cl)(Br)I",
            "Cl[Pt@SP1](Br)(I)F",
            "S[As@TB5](F)(N)(Cl)Br",
            "O[Mn@OH28](Cl)(C)(N)F",
        ] {
            let molecule = parse_molecule(input).unwrap().1;
            assert_eq!(
                canonical_smiles(&molecule),
                Err(CanonicalError::Write(WriteError::UnsupportedStereo))
            );
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
}
