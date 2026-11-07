use winnow::ascii::digit1;
use winnow::combinator::{opt, preceded};
use winnow::prelude::*;

/// Parse a hydrogen count. Possible inputs are:
/// * H - 1 hydrogen
/// * Hn - n hydrogens
///
/// Returns the number of hydrogens
pub fn parse_hydrogens(input: &mut &str) -> ModalResult<u8> {
    preceded('H', opt(digit1.try_map(str::parse::<u8>)))
        .map(|count| count.unwrap_or(1))
        .parse_next(input)
}

#[cfg(test)]
mod tests {
    use crate::io::smiles::reader::hydrogens::parse_hydrogens;
    use winnow::Parser;

    #[test]
    fn parse_hydrogens_empty() {
        assert!(parse_hydrogens.parse_peek("").is_err());
    }

    #[test]
    fn parse_hydrogens_h() {
        assert_eq!(parse_hydrogens.parse_peek("H").unwrap().1, 1);
    }

    #[test]
    fn parse_hydrogens_h1() {
        assert_eq!(parse_hydrogens.parse_peek("H1").unwrap().1, 1);
    }

    #[test]
    fn parse_hydrogens_h2() {
        assert_eq!(parse_hydrogens.parse_peek("H2").unwrap().1, 2);
    }

    #[test]
    fn parse_hydrogens_ha() {
        assert_eq!(parse_hydrogens.parse_peek("Ha").unwrap(), ("a", 1));
    }
}
