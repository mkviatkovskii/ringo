use crate::model::element::Element;
use crate::model::molecule::Molecule;

impl Molecule {
    pub fn weight(&self) -> f64 {
        let hydrogen = Element { atomic_number: 1 }.atomic_weight();
        let mut weight: f64 = 0.0;
        for atom in self.graph.node_weights() {
            weight += atom.element.atomic_weight();
            weight += (f64::from(atom.hs) + f64::from(atom.implicit_hs)) * hydrogen;
        }
        weight
    }
}
