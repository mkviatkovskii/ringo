use crate::model::element::Element;

#[derive(Hash, Eq, PartialEq, Debug)]
pub struct Atom {
    pub element: Element,
    pub isotope: u16,
    pub charge: i8,
    pub hs: u8,
    // hydrogens implied by valence for atoms written without brackets
    pub implicit_hs: u8,
    // written in lowercase, such as `c` or `[nH]`
    pub aromatic: bool,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atom() {
        let atom = Atom {
            element: Element { atomic_number: 6 },
            isotope: 12,
            charge: 0,
            hs: 0,
            implicit_hs: 0,
            aromatic: false,
        };
        assert_eq!(atom.element, Element { atomic_number: 6 });
        assert_eq!(atom.isotope, 12);
        assert_eq!(atom.charge, 0);
        assert_eq!(atom.hs, 0);
    }
}
