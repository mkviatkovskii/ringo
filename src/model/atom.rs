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
