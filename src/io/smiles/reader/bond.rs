use crate::model::bond::{Bond, BondOrder};
use winnow::combinator::alt;
use winnow::prelude::*;

pub fn parse_bond(input: &mut &str) -> ModalResult<Bond> {
    alt((
        '='.value(BondOrder::Double),
        '#'.value(BondOrder::Triple),
        ':'.value(BondOrder::Aromatic),
        '-'.value(BondOrder::Single),
    ))
    .map(|order| Bond { order })
    .parse_next(input)
}

#[cfg(test)]
mod tests {
    use crate::io::smiles::reader::bond::parse_bond;
    use crate::model::bond::BondOrder;
    use winnow::Parser;

    #[test]
    fn parse_bond_empty() {
        assert!(parse_bond.parse_peek("").is_err());
    }

    #[test]
    fn parse_bond_single() {
        assert_eq!(
            parse_bond.parse_peek("-").unwrap().1.order,
            BondOrder::Single
        );
    }

    #[test]
    fn parse_bond_double() {
        assert_eq!(
            parse_bond.parse_peek("=").unwrap().1.order,
            BondOrder::Double
        );
    }

    #[test]
    fn parse_bond_triple() {
        assert_eq!(
            parse_bond.parse_peek("#").unwrap().1.order,
            BondOrder::Triple
        );
    }

    #[test]
    fn parse_bond_aromatic() {
        assert_eq!(
            parse_bond.parse_peek(":").unwrap().1.order,
            BondOrder::Aromatic
        );
    }
}
