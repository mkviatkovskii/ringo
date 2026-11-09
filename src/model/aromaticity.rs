use crate::model::bond::BondOrder;
use crate::model::molecule::Molecule;
use petgraph::stable_graph::{EdgeIndex, NodeIndex};
use petgraph::visit::{EdgeRef, IntoEdgeReferences};
use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet, VecDeque};

impl Molecule {
    /// Marks every ring that obeys Hückel's rule (4n+2 π electrons) as aromatic: its
    /// atoms aromatic, its bonds `Aromatic`. A ring written with aromatic bonds only
    /// (lowercase atoms, or `:`) that fails the rule is still taken as aromatic. Any other
    /// `Aromatic` bond becomes single. Afterwards a molecule has the same atoms and bonds
    /// whether it was written in Kekulé form ("C1=CC=CC=C1") or aromatic form ("c1ccccc1").
    pub fn perceive_aromaticity(&mut self) {
        let rings = self.rings();
        let ring_bonds: HashSet<EdgeIndex> = rings
            .iter()
            .flat_map(|ring| self.ring_edges(ring))
            .collect();

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
    fn huckel(
        &self,
        atoms: impl Iterator<Item = NodeIndex>,
        ring_bonds: &HashSet<EdgeIndex>,
    ) -> bool {
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
    fn pi_electrons(&self, node: NodeIndex, ring_bonds: &HashSet<EdgeIndex>) -> Option<u32> {
        let atom = &self.graph[node];
        let mut ring_double = false;
        let mut ring_aromatic = false;
        let mut exocyclic_double = None;
        for edge in self.graph.edges(node) {
            let in_ring = ring_bonds.contains(&edge.id());
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

    /// Every shortest ring through each bond, in ring order, without duplicates. Taking all
    /// of them when several tie, as around the central atom of phenalene, keeps the set the
    /// same whatever the atom order. For the small fused rings of drug-like molecules this is
    /// the usual smallest set of rings.
    fn rings(&self) -> Vec<Vec<NodeIndex>> {
        let mut seen: HashSet<Vec<EdgeIndex>> = HashSet::new();
        let mut rings = Vec::new();
        for edge in (&self.graph)
            .edge_references()
            .filter(|edge| edge.weight().order != BondOrder::Dative)
        {
            for ring in self.shortest_paths_without(edge.source(), edge.target(), edge.id()) {
                let mut key: Vec<EdgeIndex> = self.ring_edges(&ring).collect();
                key.sort_unstable();
                if seen.insert(key) {
                    rings.push(ring);
                }
            }
        }
        rings
    }

    /// The bonds a ring search may follow from `node`: all but `skip` and dative bonds.
    fn ring_neighbours(
        &self,
        node: NodeIndex,
        skip: EdgeIndex,
    ) -> impl Iterator<Item = NodeIndex> + '_ {
        self.graph
            .edges(node)
            .filter(move |edge| edge.id() != skip && edge.weight().order != BondOrder::Dative)
            .map(move |edge| {
                if edge.source() == node {
                    edge.target()
                } else {
                    edge.source()
                }
            })
    }

    /// All shortest paths from `to` back to `from` that don't use `skip`; each closes into a
    /// ring through `skip`.
    fn shortest_paths_without(
        &self,
        from: NodeIndex,
        to: NodeIndex,
        skip: EdgeIndex,
    ) -> Vec<Vec<NodeIndex>> {
        let mut distance: HashMap<NodeIndex, usize> = HashMap::from([(from, 0)]);
        let mut queue = VecDeque::from([from]);
        while let Some(node) = queue.pop_front() {
            if node == to {
                break;
            }
            let next_distance = distance[&node] + 1;
            for next in self.ring_neighbours(node, skip) {
                if let Entry::Vacant(entry) = distance.entry(next) {
                    entry.insert(next_distance);
                    queue.push_back(next);
                }
            }
        }
        if !distance.contains_key(&to) {
            return Vec::new();
        }
        // walk back from `to`, one step closer to `from` each time
        let mut paths = Vec::new();
        let mut stack = vec![vec![to]];
        while let Some(path) = stack.pop() {
            let last = path[path.len() - 1];
            if last == from {
                paths.push(path);
                continue;
            }
            for previous in self.ring_neighbours(last, skip) {
                if distance.get(&previous) == Some(&(distance[&last] - 1)) {
                    let mut longer = path.clone();
                    longer.push(previous);
                    stack.push(longer);
                }
            }
        }
        paths
    }
}

#[cfg(test)]
mod tests {
    use crate::io::smiles::reader::molecule::parse_molecule;
    use crate::math::similarity::tanimoto::tanimoto_bitset;
    use crate::model::bond::BondOrder;

    fn aromatic_atoms(smiles: &str) -> Vec<bool> {
        let m = parse_molecule(smiles).unwrap().1;
        m.graph.node_weights().map(|atom| atom.aromatic).collect()
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
    fn test_bonds_outside_aromatic_rings_are_single() {
        let m = parse_molecule("c1ccc2c(c1)Cc1ccccc12").unwrap().1;
        let aromatic = m
            .graph
            .edge_weights()
            .filter(|bond| bond.order == BondOrder::Aromatic)
            .count();
        // two benzene rings; the bond they share with the 5-ring is single
        assert_eq!(aromatic, 12);
        let single = |smiles| {
            let m = parse_molecule(smiles).unwrap().1;
            m.graph
                .edge_weights()
                .filter(|bond| bond.order == BondOrder::Single)
                .count()
        };
        assert_eq!(single("c1ccccc1:c1ccccc1"), 1);
        // biphenylene's 4-ring has aromatic atoms but is not aromatic itself
        assert_eq!(single("c1ccc2c(c1)-c1ccccc1-2"), 2);
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
            // biphenylene: the bonds written single between the benzene rings stay single
            ("C1=CC=C2C(=C1)C1=CC=CC=C12", "c1ccc2c(c1)-c1ccccc1-2"),
        ] {
            let a = parse_molecule(kekule).unwrap().1.ecfp(2, 512);
            let b = parse_molecule(aromatic).unwrap().1.ecfp(2, 512);
            assert_eq!(tanimoto_bitset(&a.0, &b.0), 1.0, "{kekule} {aromatic}");
        }
    }
}
