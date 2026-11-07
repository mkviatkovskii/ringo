use crate::model::bond::{Bond, BondOrder};
use crate::model::fingerprint::Fingerprint;
use crate::model::molecule::Molecule;
use fixedbitset::FixedBitSet;
use petgraph::algo::bridges;
use petgraph::visit::{EdgeIndexable, EdgeRef, NodeIndexable};
use std::collections::HashSet;

impl Molecule {
    /// Morgan/ECFP bit fingerprint; radius 2 is ECFP4. The bits depend neither on atom
    /// order nor on the Rust version. `bits` must be above 0; the index stores both
    /// values in its header.
    pub fn ecfp(&self, radius: usize, bits: usize) -> Fingerprint {
        let graph = &self.graph;
        let nodes: Vec<_> = graph.node_indices().collect();
        let mut fp = FixedBitSet::with_capacity(bits);

        // a bond is in a ring unless it is a bridge
        let mut ring_bonds = FixedBitSet::with_capacity(graph.edge_bound());
        ring_bonds.insert_range(..);
        for edge in bridges(graph) {
            ring_bonds.remove(edge.id().index());
        }

        // per atom, by node index: current identifier, and the bonds its environment covers
        let mut ids = vec![0u64; graph.node_bound()];
        let mut envs = vec![FixedBitSet::with_capacity(graph.edge_bound()); graph.node_bound()];
        for &node in &nodes {
            let atom = &graph[node];
            let in_ring = graph
                .edges(node)
                .any(|edge| ring_bonds.contains(edge.id().index()));
            let id = fnv1a([
                u64::from(atom.element.atomic_number),
                graph.edges(node).count() as u64,
                // every hydrogen, written or implicit, so an acid's =O and -OH differ at radius 0
                u64::from(atom.hs) + u64::from(atom.implicit_hs),
                // sign-extend, so a negative charge keeps the bits it had with `as u64`
                i64::from(atom.charge).cast_unsigned(),
                u64::from(atom.isotope),
                u64::from(in_ring),
            ]);
            ids[node.index()] = id;
            fp.insert(bit(id, bits));
        }

        // bond sets already turned into bits; the empty set is radius 0
        let mut seen = HashSet::from([FixedBitSet::with_capacity(graph.edge_bound())]);
        for layer in 1..=radius {
            let mut next: Vec<(u64, usize, FixedBitSet)> = nodes
                .iter()
                .map(|&node| {
                    let mut env = envs[node.index()].clone();
                    let mut neighbours: Vec<(u64, u64)> = graph
                        .edges(node)
                        .map(|edge| {
                            let other = if edge.source() == node {
                                edge.target()
                            } else {
                                edge.source()
                            };
                            env.insert(edge.id().index());
                            env.union_with(&envs[other.index()]);
                            (bond_code(edge.weight()), ids[other.index()])
                        })
                        .collect();
                    neighbours.sort_unstable();
                    let id = fnv1a(
                        [layer as u64, ids[node.index()]]
                            .into_iter()
                            .chain(neighbours.into_iter().flat_map(|(bond, id)| [bond, id])),
                    );
                    (id, node.index(), env)
                })
                .collect();

            // an environment that covers the same bonds as an earlier one adds no bit;
            // going by identifier keeps that choice independent of atom order
            next.sort_unstable_by_key(|(id, _, _)| *id);
            for (id, index, env) in next {
                if seen.insert(env.clone()) {
                    fp.insert(bit(id, bits));
                }
                ids[index] = id;
                envs[index] = env;
            }
        }

        Fingerprint(fp)
    }
}

/// The splitmix64 finalizer mixes the high bits of `id` into the low ones the modulo
/// keeps; FNV-1a alone only carries changes upward. Modulo in u64, because
/// `id as usize` would drop the high half on 32-bit targets.
#[expect(
    clippy::cast_possible_truncation,
    reason = "the result is below `bits`, a usize"
)]
fn bit(id: u64, bits: usize) -> usize {
    let mut z = id;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^= z >> 31;
    (z % bits as u64) as usize
}

// explicit codes, so reordering BondOrder cannot change stored fingerprints
fn bond_code(bond: &Bond) -> u64 {
    match bond.order {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => 4,
    }
}

/// 64-bit FNV-1a over little-endian words: a fixed, documented hash, unlike
/// `DefaultHasher`, whose algorithm may change between Rust releases.
fn fnv1a(words: impl IntoIterator<Item = u64>) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325;
    for byte in words.into_iter().flat_map(u64::to_le_bytes) {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

#[cfg(test)]
mod test {
    use crate::descriptors::fingerprints::ecfp::fnv1a;
    use crate::io::smiles::reader::molecule::parse_molecule;
    use crate::math::similarity::tanimoto::tanimoto_bitset;
    use crate::model::fingerprint::Fingerprint;

    fn ecfp(smiles: &str, radius: usize, bits: usize) -> Fingerprint {
        parse_molecule(smiles).unwrap().1.ecfp(radius, bits)
    }

    fn ones(fp: &Fingerprint) -> Vec<usize> {
        fp.0.ones().collect()
    }

    #[test]
    fn test_ecfp() {
        let ibuprofen = ecfp("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", 2, 512);
        let naproxen = ecfp("CC(C1=CC2=C(C=C1)C=C(C=C2)OC)C(=O)O", 2, 512);
        let sim = tanimoto_bitset(&ibuprofen.0, &naproxen.0);
        assert!(0.47 < sim && sim < 0.48, "{sim}");
    }

    #[test]
    fn test_ecfp_large_length() {
        // used to panic: bit positions went past the fixed 512-bit set
        let fp = ecfp("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", 2, 2048);
        assert_eq!(fp.0.len(), 2048);
        assert!(fp.0.count_ones(..) > 0);
    }

    #[test]
    fn test_ecfp_ignores_atom_order() {
        for (a, b) in [
            ("CCO", "OCC"),
            ("CC(C)O", "OC(C)C"),
            ("CC(=O)O", "OC(C)=O"),
            ("C1CCOC1", "O1CCCC1"),
            ("c1ccncc1", "n1ccccc1"),
            ("Cc1ccccc1", "c1ccc(C)cc1"),
        ] {
            assert_eq!(ones(&ecfp(a, 2, 512)), ones(&ecfp(b, 2, 512)), "{a} {b}");
        }
    }

    #[test]
    fn test_ecfp_radius() {
        let r0 = ecfp("CCO", 0, 512);
        let r1 = ecfp("CCO", 1, 512);
        // three different atoms: degree 1 C, degree 2 C, O
        assert_eq!(r0.0.count_ones(..), 3);
        assert!(r0.0.is_subset(&r1.0) && r1.0 != r0.0);
        // from radius 2 on, every environment covers the whole molecule seen at radius 1
        assert_eq!(ones(&ecfp("CCO", 2, 512)), ones(&ecfp("CCO", 5, 512)));
    }

    #[test]
    fn test_ecfp_bits_are_stable() {
        // pins the hash and the algorithm: if this changes, existing indexes give wrong scores
        assert_eq!(fnv1a([]), 0xcbf2_9ce4_8422_2325);
        assert_eq!(
            ones(&ecfp("CC(=O)O", 2, 512)),
            [102, 119, 120, 169, 174, 397, 428, 486]
        );
    }

    #[test]
    fn test_ecfp_counts_implicit_hydrogens() {
        // the =O and -OH oxygens of acetic acid used to look the same at radius 0
        assert_eq!(ecfp("CC(=O)O", 0, 512).0.count_ones(..), 4);
        assert_ne!(ones(&ecfp("CC=O", 2, 512)), ones(&ecfp("CCO", 2, 512)));
    }
}
