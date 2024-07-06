use crate::molecule::model::bond::{Bond, BondOrder};
use nom::branch::alt;
use nom::bytes::complete::tag;
use nom::combinator::map;
use nom::IResult;

pub fn parse_bond(input: &str) -> IResult<&str, Bond> {
    let (input, bond_order) = alt((
        map(tag("="), |_| BondOrder::Double),
        map(tag("#"), |_| BondOrder::Triple),
        map(tag(":"), |_| BondOrder::Aromatic),
        map(tag("-"), |_| BondOrder::Single),
    ))(input)?;
    Ok((input, Bond { order: bond_order }))
}

#[cfg(test)]
mod tests {
    use crate::ringo::molecule::model::bond::BondOrder;
    use crate::ringo::molecule::smiles::reader::bond::parse_bond;

    #[test]
    fn parse_bond_empty() {
        assert!(parse_bond("").is_err())
    }

    #[test]
    fn parse_bond_single() {
        assert_eq!(parse_bond("-").unwrap().1.order, BondOrder::Single);
    }

    #[test]
    fn parse_bond_double() {
        assert_eq!(parse_bond("=").unwrap().1.order, BondOrder::Double);
    }

    #[test]
    fn parse_bond_triple() {
        assert_eq!(parse_bond("#").unwrap().1.order, BondOrder::Triple);
    }

    #[test]
    fn parse_bond_aromatic() {
        assert_eq!(parse_bond(":").unwrap().1.order, BondOrder::Aromatic);
    }
}
