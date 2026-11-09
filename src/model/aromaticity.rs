use crate::model::bond::BondOrder;
use crate::model::molecule::Molecule;
use fixedbitset::FixedBitSet;
use petgraph::algo::bridges;
use petgraph::stable_graph::{EdgeIndex, NodeIndex};
use petgraph::visit::{EdgeFiltered, EdgeIndexable, EdgeRef, IntoEdgeReferences, NodeIndexable};
use std::collections::HashSet;

impl Molecule {
    /// Marks every ring that obeys Hückel's rule (4n+2 π electrons) as aromatic: its
    /// atoms aromatic, its bonds `Aromatic`. A ring written entirely with aromatic
    /// bonds that fails the rule is still taken as aromatic, as written. Any other
    /// `Aromatic` bond becomes single. Afterwards a molecule has the same atoms and bonds
    /// whether it was written in Kekulé form ("C1=CC=CC=C1") or aromatic form ("c1ccccc1").
    pub fn perceive_aromaticity(&mut self) {
        let ring_bonds = self.covalent_ring_bonds();
        let rings = self.rings(&ring_bonds);

        // aromatic by itself, so all its bonds are aromatic
        let mut alone: Vec<bool> = rings
            .iter()
            .map(|ring| self.huckel(ring.iter().copied(), &ring_bonds))
            .collect();
        let mut aromatic = alone.clone();
        // fused pairs such as azulene (5 + 7 atoms) are aromatic only around the
        // outside; the bond they share stays single
        let mut pairs = Vec::new();
        for i in 0..rings.len() {
            for j in i + 1..rings.len() {
                let shared: Vec<NodeIndex> = rings[i]
                    .iter()
                    .filter(|node| rings[j].contains(node))
                    .copied()
                    .collect();
                if let [a, b] = shared[..] {
                    if let Some((edge, _)) = self.graph.find_edge_undirected(a, b) {
                        pairs.push((i, j, edge));
                    }
                }
            }
        }
        // when one ring could pair with either of two others, as in phenalene, a shared bond
        // written `-` picks the pair, so a written molecule reads back the same in any order
        pairs.sort_by_key(|&(_, _, edge)| self.graph[edge].order == BondOrder::Aromatic);
        let mut shared_bonds: HashSet<EdgeIndex> = HashSet::new();
        for (i, j, edge) in pairs {
            if aromatic[i] || aromatic[j] {
                continue;
            }
            let union: HashSet<NodeIndex> = rings[i].iter().chain(&rings[j]).copied().collect();
            if self.huckel(union.into_iter(), &ring_bonds) {
                aromatic[i] = true;
                aromatic[j] = true;
                shared_bonds.insert(edge);
            }
        }
        for (i, ring) in rings.iter().enumerate() {
            if !aromatic[i] && self.written_aromatic(ring) {
                aromatic[i] = true;
                alone[i] = true;
            }
        }

        let mut aromatic_bonds: HashSet<EdgeIndex> = HashSet::new();
        for (i, ring) in rings.iter().enumerate().filter(|&(i, _)| aromatic[i]) {
            aromatic_bonds.extend(
                self.ring_edges(ring)
                    .filter(|edge| alone[i] || !shared_bonds.contains(edge)),
            );
            for &node in ring {
                if let Some(atom) = self.graph.node_weight_mut(node) {
                    atom.aromatic = true;
                }
            }
        }
        for edge in self.graph.edge_indices().collect::<Vec<_>>() {
            if let Some(bond) = self.graph.edge_weight_mut(edge) {
                if aromatic_bonds.contains(&edge) {
                    bond.order = BondOrder::Aromatic;
                } else if bond.order == BondOrder::Aromatic {
                    // like the bond joining biphenyl's rings, or fluorene's in its 5-ring
                    bond.order = BondOrder::Single;
                }
            }
        }
    }

    /// Whether every bond of the ring was written aromatic: with no symbol between lowercase
    /// atoms, or with `:`. A `-` or `=` stands as written, so the 4-ring of biphenylene,
    /// `c1ccc2c(c1)-c1ccccc1-2`, isn't aromatic although its atoms are.
    fn written_aromatic(&self, ring: &[NodeIndex]) -> bool {
        self.ring_edges(ring)
            .all(|edge| self.graph[edge].order == BondOrder::Aromatic)
    }

    /// Whether these atoms hold 4n+2 π electrons, each of them giving some.
    fn huckel(&self, atoms: impl Iterator<Item = NodeIndex>, ring_bonds: &FixedBitSet) -> bool {
        let mut electrons = 0;
        for node in atoms {
            match self.pi_electrons(node, ring_bonds) {
                Some(count) => electrons += count,
                None => return false,
            }
        }
        electrons % 4 == 2
    }

    /// π electrons an atom gives to an aromatic ring, or None when it cannot be in one
    /// (an sp3 carbon such as the CH2 of cyclopentadiene).
    fn pi_electrons(&self, node: NodeIndex, ring_bonds: &FixedBitSet) -> Option<u32> {
        let atom = &self.graph[node];
        let mut ring_double = false;
        let mut ring_aromatic = false;
        let mut exocyclic_double = None;
        for edge in self.graph.edges(node) {
            let in_ring = ring_bonds.contains(edge.id().index());
            match edge.weight().order {
                BondOrder::Triple => return None,
                BondOrder::Double if in_ring => ring_double = true,
                BondOrder::Aromatic if in_ring => ring_aromatic = true,
                // an aromatic bond outside rings, like biphenyl's link, ends up single
                BondOrder::Single | BondOrder::Aromatic | BondOrder::Dative => {}
                BondOrder::Double => {
                    let other = if edge.source() == node {
                        edge.target()
                    } else {
                        edge.source()
                    };
                    exocyclic_double = Some(self.graph[other].element.atomic_number);
                }
            }
        }
        match exocyclic_double {
            // the C of a C=O, as in 2-pyridone, gives none
            Some(7 | 8 | 16) => return Some(0),
            // an exocyclic C=C, as in fulvene, takes the electrons out of the ring
            Some(_) => return None,
            None => {}
        }
        if ring_double {
            return Some(1);
        }
        let neighbours = self
            .graph
            .edges(node)
            .filter(|edge| edge.weight().order != BondOrder::Dative)
            .count()
            + usize::from(atom.hs)
            + usize::from(atom.implicit_hs);
        match (atom.element.atomic_number, atom.charge) {
            // a lone pair: pyrrole's N, furan's O, thiophene's S
            (7 | 15, 0) if neighbours == 3 => Some(2),
            (7, -1) | (8 | 16 | 34 | 52, 0) if neighbours == 2 => Some(2),
            // the cyclopentadienyl anion, and an empty p orbital: tropylium, boron
            (6, -1) => Some(2),
            (6, 1) | (5, 0) => Some(0),
            // an atom written aromatic with no lone pair to give: pyridine's n, any c
            _ if ring_aromatic => Some(1),
            _ => None,
        }
    }

    /// The bonds around a ring, from each atom to the next.
    fn ring_edges<'a>(&'a self, ring: &'a [NodeIndex]) -> impl Iterator<Item = EdgeIndex> + 'a {
        ring.iter()
            .zip(ring.iter().cycle().skip(1))
            .filter_map(|(&a, &b)| self.graph.find_edge_undirected(a, b).map(|(edge, _)| edge))
    }

    /// Covalent bonds in cycles: every covalent bond except bridges.
    pub(crate) fn covalent_ring_bonds(&self) -> FixedBitSet {
        let graph = &self.graph;
        let mut ring_bonds = FixedBitSet::with_capacity(graph.edge_bound());
        let covalent =
            EdgeFiltered::from_fn(graph, |edge| edge.weight().order != BondOrder::Dative);
        for edge in graph
            .edge_indices()
            .filter(|&edge| graph[edge].order != BondOrder::Dative)
        {
            ring_bonds.insert(edge.index());
        }
        for edge in bridges(&covalent) {
            ring_bonds.remove(edge.id().index());
        }
        ring_bonds
    }

    /// One shortest ring from each end of each bond, without duplicates. Trying both
    /// ends retains local ties without enumerating combinations around macrocycles.
    // shortcut: representative rings can miss extensive ties; revisit for full corpus order independence.
    fn rings(&self, ring_bonds: &FixedBitSet) -> Vec<Vec<NodeIndex>> {
        let mut seen: HashSet<Vec<NodeIndex>> = HashSet::new();
        let mut rings = Vec::new();
        let mut parent = vec![None; self.graph.node_bound()];
        let mut queue = Vec::new();
        for edge in (&self.graph)
            .edge_references()
            .filter(|edge| ring_bonds.contains(edge.id().index()))
        {
            for (from, to) in [
                (edge.source(), edge.target()),
                (edge.target(), edge.source()),
            ] {
                let Some(ring) = self.shortest_path_without(
                    from,
                    to,
                    edge.id(),
                    ring_bonds,
                    &mut parent,
                    &mut queue,
                ) else {
                    continue;
                };
                let mut key = ring.clone();
                key.sort_unstable();
                if seen.insert(key) {
                    rings.push(ring);
                }
            }
        }
        rings
    }

    /// Breadth-first search from `from` to `to` that may not use `skip`; the path
    /// closes into a ring through `skip`.
    fn shortest_path_without(
        &self,
        from: NodeIndex,
        to: NodeIndex,
        skip: EdgeIndex,
        ring_bonds: &FixedBitSet,
        parent: &mut [Option<NodeIndex>],
        queue: &mut Vec<NodeIndex>,
    ) -> Option<Vec<NodeIndex>> {
        // Retain every queued node, including those left at an early return, so
        // the next search clears only the parent entries this search touched.
        for node in queue.drain(..) {
            parent[node.index()] = None;
        }
        parent[from.index()] = Some(from);
        queue.push(from);
        let mut head = 0;
        while head < queue.len() {
            let node = queue[head];
            head += 1;
            if node == to {
                let mut path = vec![to];
                let mut current = to;
                while current != from {
                    current = parent[current.index()]
                        .unwrap_or_else(|| panic!("visited node has a parent"));
                    path.push(current);
                }
                return Some(path);
            }
            for edge in self
                .graph
                .edges(node)
                .filter(|edge| edge.id() != skip && ring_bonds.contains(edge.id().index()))
            {
                let next = if edge.source() == node {
                    edge.target()
                } else {
                    edge.source()
                };
                if parent[next.index()].is_none() {
                    parent[next.index()] = Some(node);
                    queue.push(next);
                }
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use crate::io::smiles::reader::molecule::parse_molecule;
    use crate::math::similarity::tanimoto::tanimoto_bitset;
    use crate::model::bond::BondOrder;
    use petgraph::stable_graph::NodeIndex;
    use std::fmt::Write;

    fn aromatic_atoms(smiles: &str) -> Vec<bool> {
        let m = parse_molecule(smiles).unwrap().1;
        m.graph.node_weights().map(|atom| atom.aromatic).collect()
    }

    #[test]
    fn cycloparaphenylenes_do_not_enumerate_path_combinations() {
        let label = |n: usize| {
            if n < 10 {
                n.to_string()
            } else {
                format!("%{n}")
            }
        };
        for count in [10, 12, 14] {
            let mut smiles = "c12ccc(".to_owned();
            for i in 1..count - 1 {
                write!(smiles, "-c{}ccc(", label(i + 2)).unwrap();
            }
            write!(smiles, "-c{}ccc-2cc{}", label(count + 1), label(count + 1)).unwrap();
            for i in (1..count - 1).rev() {
                write!(smiles, ")cc{}", label(i + 2)).unwrap();
            }
            smiles.push_str(")cc1");
            let start = std::time::Instant::now();
            let molecule = parse_molecule(&smiles).unwrap().1;
            eprintln!("[{count}]CPP: {:?}", start.elapsed());
            assert_eq!(molecule.graph.node_count(), count * 6);
            assert_eq!(
                molecule
                    .graph
                    .edge_weights()
                    .filter(|b| b.order == BondOrder::Aromatic)
                    .count(),
                count * 6
            );
            assert!(
                molecule.rings(&molecule.covalent_ring_bonds()).len()
                    <= 2 * molecule.graph.edge_count()
            );
        }
    }

    #[test]
    fn explicitly_single_biphenylene_links_stay_single() {
        let molecule = parse_molecule("c1ccc2c(c1)-c1ccccc1-2").unwrap().1;
        assert!(molecule.graph.node_weights().all(|a| a.aromatic));
        assert_eq!(
            molecule
                .graph
                .edge_weights()
                .filter(|b| b.order == BondOrder::Single)
                .count(),
            2
        );
        assert_eq!(
            molecule
                .graph
                .edge_weights()
                .filter(|b| b.order == BondOrder::Aromatic)
                .count(),
            12
        );
    }

    #[test]
    fn test_kekule_rings_become_aromatic() {
        for smiles in [
            "C1=CC=CC=C1",
            "C1=CC=NC=C1",
            "C1=CC=CN1",
            "[N-]1C=CC=C1",
            "[N-]1N=NN=C1",
            "C1=COC=C1",
            "C1=CSC=C1",
            "C1=CN=CN1",
            "C1=CC=C2C=CC=CC2=C1",
            // azulene: neither the 5- nor the 7-ring alone, but both together
            "C1=CC2=CC=CC=CC2=C1",
            // 2-pyridone: the C=O carbon gives no electrons, the NH two
            "O=C1C=CC=CN1",
        ] {
            assert!(
                aromatic_atoms(smiles).iter().filter(|&&a| a).count() >= 5,
                "{smiles}"
            );
        }
    }

    #[test]
    fn test_non_aromatic_rings_stay() {
        for smiles in [
            "C1=CCC=C1",
            "C1=CC=CCC1",
            "C1=CC=CC=CC=C1",
            "C1CCCCC1",
            "C=C1C=CC=C1",
        ] {
            assert!(aromatic_atoms(smiles).iter().all(|&a| !a), "{smiles}");
        }
    }

    #[test]
    fn test_ring_search_with_long_branches_and_disconnected_chains() {
        let chain = "C".repeat(128);
        let smiles = format!("{chain}c1ccccc1{chain}.{chain}");
        let molecule = parse_molecule(&smiles).unwrap().1;
        assert_eq!(
            molecule
                .graph
                .edge_weights()
                .filter(|bond| bond.order == BondOrder::Aromatic)
                .count(),
            6
        );
        assert_eq!(
            molecule.graph.node_weights().filter(|a| a.aromatic).count(),
            6
        );
    }

    #[test]
    fn test_ring_search_reuses_buffers_across_spiro_and_disconnected_rings() {
        let molecule = parse_molecule("C1CC12CC2.C1CC1").unwrap().1;
        let mut rings: Vec<Vec<usize>> = molecule
            .rings(&molecule.covalent_ring_bonds())
            .into_iter()
            .map(|ring| {
                let mut atoms: Vec<usize> = ring.into_iter().map(NodeIndex::index).collect();
                atoms.sort_unstable();
                atoms
            })
            .collect();
        rings.sort_unstable();
        assert_eq!(rings, vec![vec![0, 1, 2], vec![2, 3, 4], vec![5, 6, 7]]);
    }

    #[test]
    fn test_dative_cycle_closures_do_not_create_covalent_ring_bonds() {
        for (smiles, expected_bonds, expected_rings) in [
            ("N->1CC[Cu]1", 0, 0),
            ("N1C=CC=C->[B]<-1", 0, 0),
            ("N->1CC[Cu]1.c1ccccc1", 6, 1),
        ] {
            let molecule = parse_molecule(smiles).unwrap().1;
            let ring_bonds = molecule.covalent_ring_bonds();
            assert_eq!(ring_bonds.count_ones(..), expected_bonds, "{smiles}");
            assert_eq!(
                molecule.rings(&ring_bonds).len(),
                expected_rings,
                "{smiles}"
            );
            for edge in molecule.graph.edge_indices() {
                if molecule.graph[edge].order == BondOrder::Dative {
                    assert!(!ring_bonds.contains(edge.index()), "{smiles}");
                }
            }
        }
    }

    #[test]
    fn test_bonds_outside_aromatic_rings_are_single() {
        let m = parse_molecule("c1ccc2c(c1)Cc1ccccc12").unwrap().1;
        let aromatic = m
            .graph
            .edge_weights()
            .filter(|bond| bond.order == BondOrder::Aromatic)
            .count();
        // two benzene rings; the bond they share with the 5-ring is single
        assert_eq!(aromatic, 12);
        let m = parse_molecule("c1ccccc1:c1ccccc1").unwrap().1;
        let single = m
            .graph
            .edge_weights()
            .filter(|bond| bond.order == BondOrder::Single)
            .count();
        assert_eq!(single, 1);
    }

    #[test]
    fn test_azulene_shared_bond_stays_single() {
        for smiles in ["C1=CC2=CC=CC=CC2=C1", "c1ccc2cccc2cc1", "c1ccc2cccc-2cc1"] {
            let m = parse_molecule(smiles).unwrap().1;
            let count = |order| {
                m.graph
                    .edge_weights()
                    .filter(|bond| bond.order == order)
                    .count()
            };
            // ten aromatic bonds around the outside, the one across the middle single
            assert_eq!(
                (count(BondOrder::Aromatic), count(BondOrder::Single)),
                (10, 1),
                "{smiles}"
            );
            assert!(m.graph.node_weights().all(|atom| atom.aromatic), "{smiles}");
        }
        // naphthalene's rings are aromatic alone, so its shared bond is aromatic
        let m = parse_molecule("C1=CC=C2C=CC=CC2=C1").unwrap().1;
        assert!(m
            .graph
            .edge_weights()
            .all(|bond| bond.order == BondOrder::Aromatic));
    }

    #[test]
    fn test_both_forms_give_the_same_fingerprint() {
        for (kekule, aromatic) in [
            (
                "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
                "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            ),
            ("C1=CC=C2C=CC=CC2=C1", "c1ccc2ccccc2c1"),
            ("C1=CC=NC=C1", "c1ccncc1"),
            ("C1=CC=CN1", "c1cc[nH]c1"),
            ("[N-]1C=CC=C1", "[n-]1cccc1"),
            ("[N-]1N=NN=C1", "[n-]1nnnc1"),
            ("C1=CSC=C1", "c1ccsc1"),
            ("C1=CC=C2C(=C1)C=CN2", "c1ccc2[nH]ccc2c1"),
            ("C1=CC=C(C=C1)C1=CC=CC=C1", "c1ccccc1c1ccccc1"),
            ("C1C2=CC=CC=C2C2=CC=CC=C12", "c1ccc2c(c1)Cc1ccccc1-2"),
            ("c1ccc2c(c1)Cc1ccccc12", "c1ccc2c(c1)Cc1ccccc1-2"),
        ] {
            let a = parse_molecule(kekule).unwrap().1.ecfp(2, 512);
            let b = parse_molecule(aromatic).unwrap().1.ecfp(2, 512);
            assert_eq!(tanimoto_bitset(&a.0, &b.0), 1.0, "{kekule} {aromatic}");
        }
    }
}
