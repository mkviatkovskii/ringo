use crate::io::smiles::reader::charge::parse_charge;
use crate::io::smiles::reader::element::parse_element;
use crate::io::smiles::reader::hydrogens::parse_hydrogens;
use crate::io::smiles::reader::isotope::parse_isotope;
use crate::model::atom::Atom;
use crate::model::element::Element;
use nom::combinator::opt;
use nom::IResult;
use crate::io::smiles::reader::chirality::parse_chirality;
use crate::model::chirality::Chirality;

pub fn parse_atom(input: &str) -> IResult<&str, Atom> {
    let mut isotope: Option<u8> = None;
    let mut charge: Option<i8> = None;
    let mut hs: Option<u8> = None;
    let mut chirality: Option<Chirality> = None;
    let atomic_number: u8;

    let (mut input, sqro_found) = opt(nom::character::complete::char('['))(input)?;
    if sqro_found.is_some() {
        (input, isotope) = opt(parse_isotope)(input).unwrap_or((input, None));
    }
    (input, atomic_number) = parse_element(input)?;
    if sqro_found.is_some() {
        (input, chirality) = opt(parse_chirality)(input).unwrap_or((input, None));
        (input, hs) = opt(parse_hydrogens)(input).unwrap_or((input, None));
        (input, charge) = opt(parse_charge)(input).unwrap_or((input, None));
        let sqrc_found: Option<char>;
        (input, sqrc_found) = opt(nom::character::complete::char(']'))(input)?;
        if (sqro_found.is_some() && sqrc_found.is_none())
            || (sqro_found.is_none() && sqrc_found.is_some())
        {
            return Err(nom::Err::Failure(nom::error::Error::new(
                input,
                nom::error::ErrorKind::Verify,
            )));
        }
    }

    Ok((
        input,
        Atom {
            element: Element {
                atomic_number: atomic_number,
            },
            isotope: isotope.unwrap_or(0),
            charge: charge.unwrap_or(0),
            hs: hs.unwrap_or(0),
            chirality: chirality.unwrap_or(Chirality::Unspecified),
        },
    ))
}

#[cfg(test)]
mod tests {
    use crate::io::smiles::reader::atom::parse_atom;
    use crate::model::chirality::Chirality;

    fn do_test_parse_atom(input: &str, atomic_number: u8, charge: i8, hs: u8, isotope: u8, chirality: Chirality) {
        let (remaining_input, atom) = parse_atom(input).unwrap();
        assert_eq!(atom.element.atomic_number, atomic_number);
        assert_eq!(atom.isotope, isotope);
        assert_eq!(atom.charge, charge);
        assert_eq!(atom.hs, hs);
        assert_eq!(remaining_input, "");
    }

    #[test]
    fn parse_atom_c() {
        do_test_parse_atom("C", 6, 0, 0, 0, Chirality::Unspecified);
    }

    #[test]
    fn parse_atom_n_plus() {
        do_test_parse_atom("[N+]", 7, 1, 0, 0, Chirality::Unspecified);
    }

    #[test]
    fn parse_atom_o_minus_1() {
        do_test_parse_atom("[O-1]", 8, -1, 0, 0, Chirality::Unspecified);
    }

    #[test]
    fn parse_atom_nh2() {
        do_test_parse_atom("[NH2]", 7, 0, 2, 0, Chirality::Unspecified);
    }

    #[test]
    fn parse_atom_cl_minus_1() {
        do_test_parse_atom("[Cl-1]", 17, -1, 0, 0, Chirality::Unspecified);
    }

    #[test]
    fn parse_atom_cl_minus_minus() {
        do_test_parse_atom("[Cl--]", 17, -2, 0, 0, Chirality::Unspecified);
    }

    #[test]
    fn parse_atom_13_c_h3_minus() {
        do_test_parse_atom("[13CH3-]", 6, -1, 3, 13, Chirality::Unspecified);
    }

    #[test]
    fn parse_atom_13_c_h3_minus_incorrect_closing() {
        assert!(parse_atom("[13CH3-C").is_err())
    }

    #[test]
    fn parse_atom_chiral() {
        do_test_parse_atom("[C@@H]", 6, 0, 1, 0, Chirality::Clockwise);
    }
}
