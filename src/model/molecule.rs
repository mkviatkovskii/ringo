use crate::model::atom::Atom;
use crate::model::bond::Bond;
use petgraph::stable_graph::{EdgeIndex, NodeIndex, StableGraph};
use petgraph::visit::EdgeRef;
use petgraph::Undirected;
use std::borrow::Borrow;
use std::collections::BTreeSet;

pub struct Molecule {
    pub graph: StableGraph<Atom, Bond, Undirected>,
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

    pub fn get_bond_by_atoms(&self, ni1: NodeIndex, ni2: NodeIndex) -> Option<&Bond> {
        let e = self.graph.find_edge_undirected(ni1, ni2);
        if e.is_none() {
            return None;
        }
        return self.graph.edge_weight(e.unwrap().0);
    }

    pub fn has_bond(&self, ni1: NodeIndex, ni2: NodeIndex) -> bool {
        return self.graph.find_edge_undirected(ni1, ni2).is_some();
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
}

#[cfg(test)]
mod test {}
