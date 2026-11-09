use crate::model::atom::Atom;
use crate::model::bond::Bond;
use petgraph::stable_graph::{EdgeIndex, NodeIndex, StableGraph};
use petgraph::visit::EdgeRef;
use petgraph::Undirected;
use std::collections::BTreeSet;

#[derive(Default)]
pub struct Molecule {
    // crate-only: the reader computes implicit hydrogens once, so outside edits would leave
    // them stale
    pub(crate) graph: StableGraph<Atom, Bond, Undirected>,
    pub(crate) stereo: Vec<crate::model::stereo::AtomStereo>,
    pub(crate) bond_stereo: Vec<crate::model::stereo::DoubleBondStereo>,
    pub(crate) allene_stereo: Vec<crate::model::stereo::AlleneStereo>,
    // Aromaticity normalization rewrites atom flags and bond orders; post-load checks read
    // them as written, by node and edge index.
    pub(crate) input_aromatic: Vec<bool>,
    pub(crate) input_orders: Vec<crate::model::bond::BondOrder>,
}

impl Molecule {
    /// Explicit allene/cumulene configurations, with their axes and terminal ligands.
    pub fn allene_stereo(&self) -> &[crate::model::stereo::AlleneStereo] {
        &self.allene_stereo
    }

    pub fn double_bond_stereo(&self) -> &[crate::model::stereo::DoubleBondStereo] {
        &self.bond_stereo
    }
    /// Explicitly supplied configurations. Absence means unspecified.
    pub fn atom_stereo(&self) -> &[crate::model::stereo::AtomStereo] {
        &self.stereo
    }

    pub fn new() -> Molecule {
        Self::default()
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
        self.graph.edges(atom).map(|edge| edge.id()).collect()
    }

    pub fn get_neighbors_for_atom(&self, atom: NodeIndex) -> BTreeSet<NodeIndex> {
        self.graph.neighbors(atom).collect()
    }
}
