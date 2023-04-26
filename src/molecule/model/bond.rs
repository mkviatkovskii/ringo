#[derive(Debug, Eq, PartialEq, Clone, Copy)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
}

pub struct Bond {
    pub order: BondOrder
}
