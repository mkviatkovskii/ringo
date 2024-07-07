use crate::model::molecule::Molecule;

impl Molecule {
    pub fn weight(&self) -> f64 {
        let mut weight: f64 = 0.0;
        for atom in self.graph.node_weights() {
            weight += atom.element.atomic_weight();
        }
        weight
    }
}
