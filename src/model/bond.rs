#[derive(Debug, Eq, PartialEq, Clone, Copy)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
}

pub struct Bond {
    pub order: BondOrder,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bond_order() {
        let bond = Bond {
            order: BondOrder::Single,
        };
        assert_eq!(bond.order, BondOrder::Single);
    }
}
