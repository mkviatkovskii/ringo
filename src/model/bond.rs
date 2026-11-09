#[derive(Debug, Eq, PartialEq, Clone, Copy, Hash, PartialOrd, Ord)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
    Dative,
}

pub struct Bond {
    pub order: BondOrder,
    /// Electron-pair donor for a dative bond; None for covalent bonds.
    pub donor: Option<NodeIndex>,
}
use petgraph::stable_graph::NodeIndex;
