//! Deterministic individualization/refinement for canonical SMILES.
//!
//! This is the graph-ranking stage of #63. Stereo marks are preserved by the writer;
//! perception of impossible stereo and stereo-aware partition refinement are pending.
use super::writer::{write, WriteError};
use crate::model::molecule::Molecule;
use crate::model::stereo::Configuration;
use petgraph::stable_graph::NodeIndex;
use petgraph::visit::{EdgeRef, NodeIndexable};
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
        remaining: limit,
        leaves: BTreeMap::new(),
        symmetries: Vec::new(),
    };
    search.visit(colors)
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

impl Search<'_> {
    fn refine(&self, mut colors: Vec<usize>) -> Vec<usize> {
        loop {
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
                    (colors[node.index()], neighbours)
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
        let molecule = parse_molecule("C1CCCCC1").unwrap().1;
        assert_eq!(
            canonical_smiles_with_limit(&molecule, 1),
            Err(CanonicalError::SearchLimit)
        );
    }
}
