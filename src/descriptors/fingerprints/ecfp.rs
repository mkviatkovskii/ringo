use crate::model::atom::Atom;
use crate::model::bond::Bond;
use crate::model::fingerprint::Fingerprint;
use crate::model::molecule::Molecule;
use fixedbitset::FixedBitSet;
use petgraph::graph::NodeIndex;
use petgraph::prelude::{EdgeRef, StableGraph};
use petgraph::Undirected;
use std::hash::{DefaultHasher, Hasher};

impl Molecule {
    /// `bits` must be above 0; the index stores both values in its header.
    pub fn ecfp(&self, radius: usize, bits: usize) -> Fingerprint {
        let mut fp = FixedBitSet::with_capacity(bits);

        for node in self.graph.node_indices() {
            ecfp_recursive(
                &self.graph,
                radius,
                1,
                node,
                &mut fp,
                &mut DefaultHasher::new(),
            );
        }

        Fingerprint(fp)
    }
}

fn ecfp_recursive(
    graph: &StableGraph<Atom, Bond, Undirected>,
    radius: usize,
    depth: usize,
    node: NodeIndex,
    fp: &mut FixedBitSet,
    hasher: &mut DefaultHasher,
) {
    if depth > radius {
        return;
    }

    let atom = graph.node_weight(node).unwrap();
    hasher.write_u8(atom.element.atomic_number);
    hasher.write_u8(atom.isotope);
    hasher.write_i8(atom.charge);
    hasher.write_u8(atom.hs);

    let value = hasher.clone().finish();
    fp.insert(value as usize % fp.len());

    for edge in graph.edges(node) {
        let bond = edge.weight();
        hasher.write_u8(bond.order as u8);

        let target = if edge.source() == node {
            edge.target()
        } else {
            edge.source()
        };

        ecfp_recursive(graph, radius, depth + 1, target, fp, hasher);
    }
}

#[cfg(test)]
mod test {
    use crate::io::smiles::reader::molecule::parse_molecule;
    use crate::math::similarity::tanimoto::tanimoto_bitset;
    #[test]
    fn test_ecfp() {
        let ecfp_ibuprofen = parse_molecule("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O")
            .unwrap()
            .1
            .ecfp(2, 512);
        let ecfp_naproxen = parse_molecule("CC(C1=CC2=C(C=C1)C=C(C=C2)OC)C(=O)O")
            .unwrap()
            .1
            .ecfp(2, 512);
        let sim = tanimoto_bitset(&ecfp_ibuprofen.0, &ecfp_naproxen.0);
        assert!(0.51 < sim && sim < 0.52);
    }

    #[test]
    fn test_ecfp_large_length() {
        // used to panic: bit positions went past the fixed 512-bit set
        let fp = parse_molecule("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O")
            .unwrap()
            .1
            .ecfp(2, 2048);
        assert_eq!(fp.0.len(), 2048);
        assert!(fp.0.count_ones(..) > 0);
    }
}
