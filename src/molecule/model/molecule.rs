use std::alloc::System;
use crate::molecule::model::atom::Atom;
use crate::molecule::model::bond::Bond;
use crate::molecule::model::element::atomic_weight;
use petgraph::stable_graph::{EdgeIndex, NodeIndex, StableGraph};
use petgraph::visit::EdgeRef;
use petgraph::Undirected;
use std::borrow::Borrow;
use std::collections::{BTreeSet, HashSet};
use std::hash::Hasher;
use std::iter::FromIterator;
use twox_hash::XxHash64;
use bit_set::BitSet;
use bit_vec::BitVec;

pub struct Molecule {
    graph: StableGraph<Atom, Bond, Undirected>,
}

impl Molecule {
    pub fn new() -> Molecule {
        Molecule {
            graph: Default::default(),
        }
    }

    pub fn add_atom(&mut self, atom: Atom) -> NodeIndex {
        return self.graph.add_node(atom);
    }

    pub fn add_bond(&mut self, atom1: NodeIndex, atom2: NodeIndex, bond: Bond) -> EdgeIndex {
        return self.graph.add_edge(atom1, atom2, bond);
    }

    pub fn get_atom(&self, node: NodeIndex) -> Option<&Atom> {
        return self.graph.node_weight(node);
    }

    pub fn get_bond(&self, edge: EdgeIndex) -> Option<&Bond> {
        return self.graph.edge_weight(edge);
    }

    pub fn count_atoms(&self) -> usize {
        return self.graph.node_count();
    }

    pub fn count_bonds(&self) -> usize {
        return self.graph.edge_count();
    }

    pub fn get_bonds_for_atom(&self, atom: NodeIndex) -> Vec<EdgeIndex> {
        let mut bonds: Vec<EdgeIndex> = Vec::new();
        for edge in self.graph.edges(atom) {
            bonds.push(EdgeIndex::new(edge.id().index()))
        }
        return bonds;
    }

    pub fn get_neighbors_for_atom(&self, atom: NodeIndex) -> BTreeSet<NodeIndex> {
        let mut atoms: BTreeSet<NodeIndex> = BTreeSet::new();
        for edge in self.graph.edges(atom) {
            atoms.insert(NodeIndex::new(edge.source().index()));
            atoms.insert(NodeIndex::new(edge.target().index()));
        }
        atoms.remove(atom.borrow());
        return atoms;
    }

    // TODO: move to Descriptors crate
    pub fn weight(&self) -> f64 {
        let mut weight: f64 = 0.0;
        for atom in self.graph.node_weights() {
            weight += atomic_weight(atom.element.borrow())
        }
        weight
    }

    // TODO: move to Descriptors crate
    pub fn ecfp(&self, radius: usize, fp_length: usize) -> BitVec {
        let mut fp = BitSet::new();

        for node in self.graph.node_indices() {
            ecfp_recursive(&self.graph, radius, 1, node, &mut fp, fp_length);
        }

        BitVec::from_fn(fp_length, |idx| fp.contains(idx))
    }
}

fn ecfp_recursive(
    graph: &StableGraph<Atom, Bond, Undirected>,
    radius: usize,
    depth: usize,
    node: NodeIndex,
    fp: &mut BitSet,
    fp_length: usize,
) {
    //println!("ecfp_recursive({}, {}, {})", radius, depth, node.index());

    if depth > radius {
        return;
    }

    let mut hasher = XxHash64::with_seed(0);

    let atom = graph.node_weight(node).unwrap();
    hasher.write_u8(atom.element.atomic_number);
    hasher.write_u8(atom.isotope);
    hasher.write_i8(atom.charge);
    hasher.write_u8(atom.hs);

    fp.insert(hasher.finish() as usize % fp_length);

    for edge in graph.edges(node) {
        let bond = edge.weight();
        hasher.write_u8(bond.order as u8);

        let target = if edge.source() == node {
            edge.target()
        } else {
            edge.source()
        };

        ecfp_recursive(graph, radius, depth + 1, target, fp, fp_length);
    }
}
