use nom;
use nom::branch::alt;
use nom::bytes::streaming::tag;
use nom::combinator::{map, value};
use nom::IResult;

use crate::model::chirality::Chirality;

/// Parse a stereocenter chirality configuration. Possible inputs are:
/// @ - anticlockwise
/// @@ - clockwise
/// Returns the configuration
pub fn parse_chirality(input: &str) -> IResult<&str, Chirality> {
    alt((
        map(tag("@@"), |_| Chirality::Clockwise),
        map(tag("@"), |_| Chirality::Anticlockwise),
        value(Chirality::Unspecified, tag("")),
    ))(input)
}

#[cfg(test)]
mod tests {
    use crate::io::smiles::reader::chirality::parse_chirality;
    use crate::model::chirality::Chirality;

    #[test]
    fn test_parse_chirality_unspecified() {
        let (remaining_input, chirality) = parse_chirality("H]").unwrap();
        assert_eq!(chirality, Chirality::Unspecified);
        assert_eq!(remaining_input, "H]");
    }

    #[test]
    fn test_parse_chirality_anticlockwise() {
        let (remaining_input, chirality) = parse_chirality("@H]").unwrap();
        assert_eq!(chirality, Chirality::Anticlockwise);
        assert_eq!(remaining_input, "H]");
    }

    #[test]
    fn test_parse_chirality_clockwise() {
        let (remaining_input, chirality) = parse_chirality("@@H]").unwrap();
        assert_eq!(chirality, Chirality::Clockwise);
        assert_eq!(remaining_input, "H]");
    }
}
