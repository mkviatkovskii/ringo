use nom::{character::complete::digit1, combinator::map_res, IResult, Parser};

/// Parses isotope value, that should be a number
/// Returns the isotope value
pub fn parse_isotope(input: &str) -> IResult<&str, u8> {
    map_res(digit1, str::parse::<u8>).parse(input)
}

#[cfg(test)]
mod tests {
    use crate::io::smiles::reader::isotope::parse_isotope;
    use nom::error::{Error, ErrorKind};

    #[test]
    fn parse_isotope_empty() {
        assert_eq!(
            parse_isotope(""),
            Err(nom::Err::Error(Error::new("", ErrorKind::Digit)))
        );
    }

    #[test]
    fn parse_isotope_2() {
        assert_eq!(parse_isotope("2").unwrap().1, 2);
    }

    #[test]
    fn parse_isotope_12() {
        assert_eq!(parse_isotope("12").unwrap().1, 12);
    }
}
