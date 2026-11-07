use winnow::combinator::alt;
use winnow::prelude::*;

pub fn parse_element(input: &mut &str) -> ModalResult<u8> {
    // two-letter symbols first, so "Cl" is not read as "C"
    alt((
        alt(("Cl".value(17), "Br".value(35))),
        alt((
            "C".value(6),
            "N".value(7),
            "O".value(8),
            "F".value(9),
            "P".value(15),
            "S".value(16),
            "I".value(53),
            "H".value(1),
        )),
    ))
    .parse_next(input)
}

#[cfg(test)]
mod tests {
    use crate::io::smiles::reader::element::parse_element;
    use winnow::Parser;

    #[test]
    fn parse_atom_symbol_empty() {
        assert!(parse_element.parse_peek("").is_err());
    }

    #[test]
    fn parse_atom_symbol_cl() {
        assert_eq!(parse_element.parse_peek("Cl").unwrap().1, 17);
    }

    #[test]
    fn parse_atom_symbol_c() {
        assert_eq!(parse_element.parse_peek("C").unwrap().1, 6);
    }
}
