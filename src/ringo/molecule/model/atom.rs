use crate::ringo::molecule::model::element::Element;

#[derive(Hash, Eq, PartialEq, Debug)]
pub struct Atom {
    pub element: Element,
    pub isotope: u8,
    pub charge: i8,
    pub hs: u8,
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
        };
        assert_eq!(atom.element, Element { atomic_number: 6 });
        assert_eq!(atom.isotope, 12);
        assert_eq!(atom.charge, 0);
        assert_eq!(atom.hs, 0);
    }
}
