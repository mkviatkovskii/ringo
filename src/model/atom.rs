use crate::model::element::Element;

#[derive(Hash, Eq, PartialEq, Debug)]
pub struct Atom {
    pub element: Element,
    pub isotope: u16,
    /// SMILES atom class/map label; absent labels default to zero, not chemical identity.
    pub atom_class: u32,
    pub charge: i8,
    pub hs: u8,
    // hydrogens implied by valence for atoms written without brackets
    pub implicit_hs: u8,
    // written in lowercase, such as `c` or `[nH]`
    pub aromatic: bool,
}
