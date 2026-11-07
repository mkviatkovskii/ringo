use crate::model::atom::Atom;
use crate::model::bond::Bond;
use petgraph::stable_graph::{EdgeIndex, NodeIndex, StableGraph};
use petgraph::visit::EdgeRef;
use petgraph::Undirected;
use std::borrow::Borrow;
use std::collections::BTreeSet;

pub struct Molecule {
    // crate-only: the reader computes implicit hydrogens once, so outside edits would leave
    // them stale
    pub(crate) graph: StableGraph<Atom, Bond, Undirected>,
}

impl Default for Molecule {
    fn default() -> Self {
        Self::new()
    }
}

impl Molecule {
    pub fn new() -> Molecule {
        Molecule {
            graph: StableGraph::default(),
        }
    }

    pub(crate) fn add_atom(&mut self, atom: Atom) -> NodeIndex {
        self.graph.add_node(atom)
    }

    pub(crate) fn add_bond(&mut self, atom1: NodeIndex, atom2: NodeIndex, bond: Bond) -> EdgeIndex {
        self.graph.add_edge(atom1, atom2, bond)
    }

    /// The atoms and bonds, read-only.
    pub fn graph(&self) -> &StableGraph<Atom, Bond, Undirected> {
        &self.graph
    }

    pub fn get_atom(&self, node: NodeIndex) -> Option<&Atom> {
        self.graph.node_weight(node)
    }

    pub fn get_bond(&self, edge: EdgeIndex) -> Option<&Bond> {
        self.graph.edge_weight(edge)
    }

    pub fn get_bond_by_atoms(&self, ni1: NodeIndex, ni2: NodeIndex) -> Option<&Bond> {
        let (edge, _) = self.graph.find_edge_undirected(ni1, ni2)?;
        self.graph.edge_weight(edge)
    }

    pub fn has_bond(&self, ni1: NodeIndex, ni2: NodeIndex) -> bool {
        self.graph.find_edge_undirected(ni1, ni2).is_some()
    }

    pub fn count_atoms(&self) -> usize {
        self.graph.node_count()
    }

    pub fn count_bonds(&self) -> usize {
        self.graph.edge_count()
    }

    pub fn get_bonds_for_atom(&self, atom: NodeIndex) -> Vec<EdgeIndex> {
        let mut bonds: Vec<EdgeIndex> = Vec::new();
        for edge in self.graph.edges(atom) {
            bonds.push(EdgeIndex::new(edge.id().index()));
        }
        bonds
    }

    pub fn get_neighbors_for_atom(&self, atom: NodeIndex) -> BTreeSet<NodeIndex> {
        let mut atoms: BTreeSet<NodeIndex> = BTreeSet::new();
        for edge in self.graph.edges(atom) {
            atoms.insert(NodeIndex::new(edge.source().index()));
            atoms.insert(NodeIndex::new(edge.target().index()));
        }
        atoms.remove(atom.borrow());
        atoms
    }
}

#[cfg(test)]
mod test {}
