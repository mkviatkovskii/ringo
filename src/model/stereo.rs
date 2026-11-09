//! Configurations relative to ordered ligand references, never CIP labels.
use petgraph::stable_graph::EdgeIndex;
use petgraph::stable_graph::NodeIndex;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Ligand {
    Atom(NodeIndex),
    Hydrogen,
    LonePair,
    /// An unoccupied site in explicitly tagged SP/TB/OH coordination.
    Vacancy,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Configuration {
    Tetrahedral(u8),
    /// SMILES SP1, SP2, SP3: opposite ligand indices are (0, 2), (0, 1), (0, 3).
    SquarePlanar(u8),
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
pub struct AlleneStereo {
    pub center: NodeIndex,
    /// Terminal-to-terminal double-bond chain, ordered by the center's SMILES neighbours.
    /// The center is its middle atom.
    pub axis: Vec<NodeIndex>,
    /// Two SMILES-ordered substituents at axis[0], then two at its last atom.
    pub ligands: [Ligand; 4],
    /// AL1/AL2 (or @/@@), interpreted by collapsing the axis to a tetrahedron.
    pub permutation: u8,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DoubleBondStereo {
    pub bond: EdgeIndex,
    pub ends: [NodeIndex; 2],
    pub references: [NodeIndex; 2],
    pub same_side: bool,
}
