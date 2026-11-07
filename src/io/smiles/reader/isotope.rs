use winnow::ascii::digit1;
use winnow::prelude::*;

/// Parses isotope value, that should be a number
/// Returns the isotope value
pub fn parse_isotope(input: &mut &str) -> ModalResult<u8> {
    digit1.try_map(str::parse::<u8>).parse_next(input)
}

#[cfg(test)]
mod tests {
    use crate::io::smiles::reader::isotope::parse_isotope;
    use winnow::Parser;

    #[test]
    fn parse_isotope_empty() {
        assert!(parse_isotope.parse_peek("").is_err());
    }

    #[test]
    fn parse_isotope_2() {
        assert_eq!(parse_isotope.parse_peek("2").unwrap().1, 2);
    }

    #[test]
    fn parse_isotope_12() {
        assert_eq!(parse_isotope.parse_peek("12").unwrap().1, 12);
    }
}
