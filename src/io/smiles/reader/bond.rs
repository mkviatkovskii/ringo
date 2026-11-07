use crate::model::bond::{Bond, BondOrder};
use petgraph::stable_graph::NodeIndex;
use winnow::combinator::alt;
use winnow::prelude::*;

#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) struct ParsedBond {
    pub order: BondOrder,
    donor_at_start: Option<bool>,
}

impl ParsedBond {
    pub(crate) fn covalent(order: BondOrder) -> Self {
        Self {
            order,
            donor_at_start: None,
        }
    }

    pub(crate) fn reversed(mut self) -> Self {
        self.donor_at_start = self.donor_at_start.map(|at_start| !at_start);
        self
    }

    pub(crate) fn bind(self, start: NodeIndex, end: NodeIndex) -> Bond {
        Bond {
            order: self.order,
            donor: self
                .donor_at_start
                .map(|at_start| if at_start { start } else { end }),
        }
    }
}

pub(crate) fn parse_bond(input: &mut &str) -> ModalResult<ParsedBond> {
    alt((
        "->".value(ParsedBond {
            order: BondOrder::Dative,
            donor_at_start: Some(true),
        }),
        "<-".value(ParsedBond {
            order: BondOrder::Dative,
            donor_at_start: Some(false),
        }),
        '='.value(ParsedBond::covalent(BondOrder::Double)),
        '#'.value(ParsedBond::covalent(BondOrder::Triple)),
        ':'.value(ParsedBond::covalent(BondOrder::Aromatic)),
        '-'.value(ParsedBond::covalent(BondOrder::Single)),
    ))
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
