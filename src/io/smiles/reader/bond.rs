use crate::model::bond::{Bond, BondDirection, BondOrder};
use nom::branch::alt;
use nom::bytes::complete::tag;
use nom::combinator::map;
use nom::IResult;

pub fn parse_bond(input: &str) -> IResult<&str, Bond> {
    let (input, (bond_order, bond_direction)) = alt((
        map(tag("="), |_| (BondOrder::Double, BondDirection::Unspecified)),
        map(tag("#"), |_| (BondOrder::Triple, BondDirection::Unspecified)),
        map(tag(":"), |_| (BondOrder::Aromatic, BondDirection::Unspecified)),
        map(tag("-"), |_| (BondOrder::Single, BondDirection::Unspecified)),
        map(tag("/"), |_| (BondOrder::Single, BondDirection::Up)),
        map(tag(r"\"), |_| (BondOrder::Single, BondDirection::Down)),
    ))(input)?;
    Ok((input, Bond { order: bond_order, direction: bond_direction}))
}

#[cfg(test)]
mod tests {
    use crate::io::smiles::reader::bond::parse_bond;
    use crate::model::bond::{BondDirection, BondOrder};

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

    #[test]
    fn parse_bond_cistrans_up() {
        let bond = parse_bond("/").unwrap().1;
        assert_eq!(bond.order, BondOrder::Single);
        assert_eq!(bond.direction, BondDirection::Up);
    }

    #[test]
    fn parse_bond_cistrans_down() {
        let bond = parse_bond(r"\").unwrap().1;
        assert_eq!(bond.order, BondOrder::Single);
        assert_eq!(bond.direction, BondDirection::Down);
    }
}
