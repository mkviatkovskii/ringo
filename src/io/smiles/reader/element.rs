use crate::model::element::Element;
use winnow::combinator::alt;
use winnow::prelude::*;
use winnow::token::one_of;

/// An element outside brackets: the organic subset, its lowercase aromatic forms and
/// the `*` wildcard. Returns the atomic number and whether the atom is aromatic.
pub fn parse_organic_element(input: &mut &str) -> ModalResult<(u8, bool)> {
    alt((
        // two-letter symbols first, so "Cl" is not read as "C"
        alt((
            "Cl",
            "Br",
            one_of(['B', 'C', 'N', 'O', 'F', 'P', 'S', 'I', 'H', '*']).take(),
        ))
        .verify_map(Element::from_symbol)
        .map(|element| (element.atomic_number, false)),
        one_of(AROMATIC_ORGANIC)
            .take()
            .verify_map(aromatic_number)
            .map(|atomic_number| (atomic_number, true)),
    ))
    .parse_next(input)
}

/// An element inside brackets: any element symbol, an aromatic one such as `c` or
/// `se`, or `*`. Returns the atomic number and whether the atom is aromatic.
pub fn parse_bracket_element(input: &mut &str) -> ModalResult<(u8, bool)> {
    alt((
        alt(("se", "as", "te", one_of(AROMATIC_ORGANIC).take()))
            .verify_map(aromatic_number)
            .map(|atomic_number| (atomic_number, true)),
        // two letters first, so "[Sc]" is scandium; inside brackets nothing else
        // can follow a symbol with a lowercase letter
        alt((
            (one_of('A'..='Z'), one_of('a'..='z')).take(),
            one_of(('A'..='Z', '*')).take(),
        ))
        .verify_map(Element::from_symbol)
        .map(|element| (element.atomic_number, false)),
    ))
    .parse_next(input)
}

const AROMATIC_ORGANIC: [char; 6] = ['b', 'c', 'n', 'o', 'p', 's'];

fn aromatic_number(symbol: &str) -> Option<u8> {
    match symbol {
        "b" => Some(5),
        "c" => Some(6),
        "n" => Some(7),
        "o" => Some(8),
        "p" => Some(15),
        "s" => Some(16),
        "as" => Some(33),
        "se" => Some(34),
        "te" => Some(52),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use crate::io::smiles::reader::element::{parse_bracket_element, parse_organic_element};
    use winnow::Parser;

    #[test]
    fn parse_atom_symbol_empty() {
        assert!(parse_organic_element.parse_peek("").is_err());
        assert!(parse_bracket_element.parse_peek("").is_err());
    }

    #[test]
    fn parse_organic_elements() {
        for (input, expected, rest) in [
            ("Cl", (17, false), ""),
            ("C", (6, false), ""),
            ("Br", (35, false), ""),
            ("B", (5, false), ""),
            ("*", (0, false), ""),
            ("c", (6, true), ""),
            ("n", (7, true), ""),
            // outside brackets "Sc" is sulfur, then an aromatic carbon
            ("Sc", (16, false), "c"),
        ] {
            assert_eq!(
                parse_organic_element.parse_peek(input).unwrap(),
                (rest, expected),
                "{input}"
            );
        }
        // not in the organic subset
        assert!(parse_organic_element
            .parse_peek("Na")
            .is_ok_and(|(rest, _)| rest == "a"));
        assert!(parse_organic_element.parse_peek("Z").is_err());
    }

    #[test]
    fn parse_bracket_elements() {
        for (input, expected) in [
            ("Na", (11, false)),
            ("Sc", (21, false)),
            ("U", (92, false)),
            ("Og", (118, false)),
            ("*", (0, false)),
            ("c", (6, true)),
            ("se", (34, true)),
            ("as", (33, true)),
        ] {
            assert_eq!(
                parse_bracket_element.parse_peek(input).unwrap(),
                ("", expected),
                "{input}"
            );
        }
        // "Xx" is not an element; "X" alone is not either
        assert!(parse_bracket_element.parse_peek("Xx").is_err());
    }
}
