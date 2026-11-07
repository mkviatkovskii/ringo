//! Configurations relative to ordered ligand references, never CIP labels.
use petgraph::stable_graph::EdgeIndex;
use petgraph::stable_graph::NodeIndex;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ligand {
    Atom(NodeIndex),
    Hydrogen,
    LonePair,
    /// An unoccupied site in explicitly tagged TB/OH coordination.
    Vacancy,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Configuration {
    Tetrahedral(u8),
    TrigonalBipyramidal(u8),
    Octahedral(u8),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomStereo {
    pub center: NodeIndex,
    /// SMILES encounter order, including hydrogen and ring-closure positions.
    pub ligands: Vec<Ligand>,
    pub configuration: Configuration,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DoubleBondStereo {
    pub bond: EdgeIndex,
    pub ends: [NodeIndex; 2],
    pub references: [NodeIndex; 2],
    pub same_side: bool,
}
