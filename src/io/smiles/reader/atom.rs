use crate::io::smiles::reader::charge::parse_charge;
use crate::io::smiles::reader::element::{parse_bracket_element, parse_organic_element};
use crate::io::smiles::reader::hydrogens::parse_hydrogens;
use crate::io::smiles::reader::isotope::parse_isotope;
use crate::model::atom::Atom;
use crate::model::element::Element;
use winnow::ascii::digit1;
use winnow::combinator::{cut_err, opt, preceded};
use winnow::error::{ContextError, ErrMode};
use winnow::prelude::*;

#[derive(Clone, Copy)]
pub(crate) enum Chirality {
    Default(u8),
    TH(u8),
    AL(u8),
    SP(u8),
    TB(u8),
    OH(u8),
}

fn parse_chirality(input: &mut &str) -> ModalResult<Chirality> {
    '@'.parse_next(input)?;
    if input.starts_with('@') {
        '@'.parse_next(input)?;
        return Ok(Chirality::Default(2));
    }
    for (name, max, make) in [
        ("TH", 2, Chirality::TH as fn(u8) -> Chirality),
        ("AL", 2, Chirality::AL),
        ("SP", 3, Chirality::SP),
        ("TB", 20, Chirality::TB),
        ("OH", 30, Chirality::OH),
    ] {
        if let Some(rest) = input.strip_prefix(name) {
            let n = rest.bytes().take_while(u8::is_ascii_digit).count();
            let value = rest[..n]
                .parse::<u8>()
                .ok()
                .filter(|v| (1..=max).contains(v));
            let Some(value) = value else {
                return Err(ErrMode::Cut(ContextError::new()));
            };
            *input = &rest[n..];
            return Ok(make(value));
        }
    }
    Ok(Chirality::Default(1))
}

pub(crate) fn parse_atom(input: &mut &str) -> ModalResult<(Atom, bool, Option<Chirality>)> {
    let mut isotope: Option<u16> = None;
    let mut charge: Option<i8> = None;
    let mut hs: Option<u8> = None;
    let bracket = opt('[').parse_next(input)?.is_some();
    let mut chirality = None;
    let mut atom_class = 0;
    if bracket {
        isotope = opt(parse_isotope).parse_next(input)?;
    }
    let (atomic_number, aromatic) = if bracket {
        // after "[", anything but an element is a hard error
        cut_err(parse_bracket_element).parse_next(input)?
    } else {
        parse_organic_element(input)?
    };
    if bracket {
        chirality = opt(parse_chirality).parse_next(input)?;
        hs = opt(parse_hydrogens).parse_next(input)?;
        charge = opt(parse_charge).parse_next(input)?;
        atom_class = opt(preceded(':', cut_err(digit1.try_map(str::parse::<u32>))))
            .parse_next(input)?
            .unwrap_or(0);
        // unclosed bracket atom is a hard error, not a reason to stop parsing
        cut_err(']').parse_next(input)?;
    }

    Ok((
        Atom {
            element: Element { atomic_number },
            isotope: isotope.unwrap_or(0),
            atom_class,
            charge: charge.unwrap_or(0),
            hs: hs.unwrap_or(0),
            implicit_hs: 0,
            aromatic,
        },
        bracket,
        chirality,
    ))
}

#[cfg(test)]
mod tests {
    use crate::io::smiles::reader::atom::parse_atom;
    use winnow::Parser;

    fn do_test_parse_atom(input: &str, atomic_number: u8, charge: i8, hs: u8, isotope: u16) {
        let (remaining_input, (atom, _, _)) = parse_atom.parse_peek(input).unwrap();
        assert_eq!(atom.element.atomic_number, atomic_number);
        assert_eq!(atom.isotope, isotope);
        assert_eq!(atom.charge, charge);
        assert_eq!(atom.hs, hs);
        assert_eq!(remaining_input, "");
    }

    #[test]
    fn parse_atom_c() {
        do_test_parse_atom("C", 6, 0, 0, 0);
    }

    #[test]
    fn parse_atom_n_plus() {
        do_test_parse_atom("[N+]", 7, 1, 0, 0);
    }

    #[test]
    fn parse_atom_o_minus_1() {
        do_test_parse_atom("[O-1]", 8, -1, 0, 0);
    }

    #[test]
    fn parse_atom_nh2() {
        do_test_parse_atom("[NH2]", 7, 0, 2, 0);
    }

    #[test]
    fn parse_atom_cl_minus_1() {
        do_test_parse_atom("[Cl-1]", 17, -1, 0, 0);
    }

    #[test]
    fn parse_atom_cl_minus_minus() {
        do_test_parse_atom("[Cl--]", 17, -2, 0, 0);
    }

    #[test]
    fn parse_atom_13_c_h3_minus() {
        do_test_parse_atom("[13CH3-]", 6, -1, 3, 13);
    }

    #[test]
    fn parse_atom_13_c_h3_minus_incorrect_closing() {
        assert!(parse_atom.parse_peek("[13CH3-C").is_err());
    }
}
