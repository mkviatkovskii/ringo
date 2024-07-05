use crate::ringo::molecule::model::element::Element;

#[derive(Hash, Eq, PartialEq)]
pub struct Atom {
    pub element: Element,
    pub isotope: u8,
    pub charge: i8,
    pub hs: u8,
}
