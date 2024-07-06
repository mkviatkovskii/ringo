use nom::branch::alt;
use nom::bytes::complete::tag;
use nom::combinator::map_res;
use nom::IResult;

pub fn parse_element(input: &str) -> IResult<&str, u8> {
    map_res(
        alt((
            tag("Cl"),
            tag("Br"),
            tag("C"),
            tag("N"),
            tag("O"),
            tag("F"),
            tag("P"),
            tag("S"),
            tag("I"),
            tag("H"),
        )),
        |symbol| match symbol {
            "Cl" => Ok(17),
            "Br" => Ok(35),
            "C" => Ok(6),
            "N" => Ok(7),
            "O" => Ok(8),
            "F" => Ok(9),
            "P" => Ok(15),
            "S" => Ok(16),
            "I" => Ok(53),
            "H" => Ok(1),
            _ => Err(()),
        },
    )(input)
}

#[cfg(test)]
mod tests {
    use crate::ringo::molecule::smiles::reader::element::parse_element;
    use nom::error::{Error, ErrorKind};

    #[test]
    fn parse_atom_symbol_empty() {
        assert_eq!(
            parse_element(""),
            Err(nom::Err::Error(Error::new("", ErrorKind::Tag)))
        );
    }

    #[test]
    fn parse_atom_symbol_cl() {
        assert_eq!(parse_element("Cl").unwrap().1, 17);
    }

    #[test]
    fn parse_atom_symbol_c() {
        assert_eq!(parse_element("C").unwrap().1, 6);
    }
}
