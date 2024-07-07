#[derive(Debug, Eq, PartialEq, Clone, Copy)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
}

#[derive(Debug, Eq, PartialEq, Clone, Copy)]
pub enum BondDirection {
    Unspecified,
    Up,
    Down,
    Either,
}

pub struct Bond {
    pub order: BondOrder,
    pub direction: BondDirection
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bond_order() {
        let bond = Bond {
            order: BondOrder::Single,
            direction: BondDirection::Unspecified
        };
        assert_eq!(bond.order, BondOrder::Single);
        assert_eq!(bond.direction, BondDirection::Unspecified);
    }
}
