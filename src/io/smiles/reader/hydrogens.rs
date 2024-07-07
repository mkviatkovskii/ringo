use nom;
use nom::character::complete::{char, digit1};
use nom::combinator::map_res;
use nom::sequence::preceded;
use nom::IResult;

/// Parse a hydrogen count. Possible inputs are:
/// * H - 1 hydrogen
/// * Hn - n hydrogens
/// Returns the number of hydrogens
pub fn parse_hydrogens(input: &str) -> IResult<&str, u8> {
    let single_hydrogen_parser = char('H');
    let mut hydrogen_parser = preceded(char('H'), map_res(digit1, str::parse::<u8>));
    single_hydrogen_parser(input)?;
    let hydrogen_parser_result = hydrogen_parser(input);
    if hydrogen_parser_result.is_ok() {
        return Ok(hydrogen_parser_result?);
    }
    let single_hydrogen_parser_result = single_hydrogen_parser(input);
    if single_hydrogen_parser_result.is_ok() {
        return Ok((single_hydrogen_parser_result?.0, 1));
    }
    Err(nom::Err::Failure(nom::error::Error::new(
        input,
        nom::error::ErrorKind::Alt,
    )))
}

#[cfg(test)]
mod tests {
    use crate::io::smiles::reader::hydrogens::parse_hydrogens;
    use nom::error::{Error, ErrorKind};

    #[test]
    fn parse_hydrogens_empty() {
        assert_eq!(
            parse_hydrogens(""),
            Err(nom::Err::Error(Error::new("", ErrorKind::Char)))
        );
    }

    #[test]
    fn parse_hydrogens_h() {
        assert_eq!(parse_hydrogens("H").unwrap().1, 1);
    }

    #[test]
    fn parse_hydrogens_h1() {
        assert_eq!(parse_hydrogens("H1").unwrap().1, 1);
    }

    #[test]
    fn parse_hydrogens_h2() {
        assert_eq!(parse_hydrogens("H2").unwrap().1, 2);
    }

    #[test]
    fn parse_hydrogens_ha() {
        assert_eq!(parse_hydrogens("Ha").unwrap(), ("a", 1));
    }
}
