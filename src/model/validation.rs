use crate::model::bond::BondOrder;
use crate::model::molecule::Molecule;
use fixedbitset::FixedBitSet;
use petgraph::algo::maximum_matching;
use petgraph::graph::UnGraph;
use petgraph::stable_graph::{EdgeIndex, NodeIndex};
use petgraph::visit::{EdgeRef, IntoEdgeReferences};
use std::collections::HashMap;

impl Molecule {
    /// Runs optional post-load chemical checks without modifying the molecule: main-group
    /// valence limits, then the written aromatic annotations. Stereogenicity is not checked.
    /// Errors identify the zero-based atom or bond index and what is wrong there.
    pub fn validate(&self) -> Result<(), String> {
        self.validate_valences()?;
        self.validate_aromaticity()
    }

    fn validate_valences(&self) -> Result<(), String> {
        for node in self.graph.node_indices() {
            let atom = &self.graph[node];
            if atom.element.atomic_number != 0
                && i16::from(atom.charge) > i16::from(atom.element.atomic_number)
            {
                return Err(format!(
                    "{}: charge exceeds atomic number",
                    self.describe(node)
                ));
            }
            let valence = self.written_valence(node) + u32::from(self.needs_ring_pi(node));
            // An isolated cation has no bonding valence to limit, even across closed shells.
            if atom.charge > 0 && valence == 0 && self.graph.edges(node).next().is_none() {
                continue;
            }
            let Some(maximum) = maximum_valence(atom.element.atomic_number, atom.charge)
                .map_err(|()| format!("{}: unsupported charge state", self.describe(node)))?
            else {
                continue;
            };
            if valence > u32::from(maximum) {
                return Err(format!(
                    "{}: valence {valence} exceeds maximum {maximum}",
                    self.describe(node)
                ));
            }
        }
        Ok(())
    }

    /// `OpenSMILES` requires that electrons can be assigned to the atoms and bonds written
    /// aromatic, consistent with their hydrogens, charges and external bonds. Ringo checks
    /// that every aromatic atom is in a ring, and that the atoms needing a ring π bond pair
    /// up along aromatic ring bonds: a Kekulé structure exists. An aromatic bond outside
    /// rings, like `c1ccccc1:c1ccccc1`'s link, is single, as loading reads it. Hückel's rule
    /// is not required: `OpenSMILES` also writes anti-aromatic rings in lowercase.
    fn validate_aromaticity(&self) -> Result<(), String> {
        let aromatic: Vec<NodeIndex> = self
            .graph
            .node_indices()
            .filter(|&node| self.aromatic_as_written(node))
            .collect();
        if aromatic.is_empty() {
            return Ok(());
        }
        let ring_bonds = self.covalent_ring_bonds();
        for &node in &aromatic {
            let mut in_ring = false;
            for edge in self.graph.edges(node) {
                let ring_bond = ring_bonds.contains(edge.id().index());
                let other = if edge.source() == node {
                    edge.target()
                } else {
                    edge.source()
                };
                if ring_bond
                    && self.graph[other].element.atomic_number == 0
                    && self.wildcard_in_aromatic_system(edge.id(), &ring_bonds)
                {
                    return Err(format!(
                        "{}: a wildcard atom in an aromatic ring is not supported",
                        self.describe(node)
                    ));
                }
                in_ring |= ring_bond;
            }
            if !in_ring {
                return Err(format!(
                    "{}: aromatic atom outside a ring",
                    self.describe(node)
                ));
            }
        }

        // a perfect matching of the atoms that need a ring π bond is a Kekulé structure
        let mut pairing: UnGraph<NodeIndex, ()> = UnGraph::default();
        let mut index = HashMap::new();
        for &node in aromatic.iter().filter(|&&node| self.needs_ring_pi(node)) {
            index.insert(node, pairing.add_node(node));
        }
        for edge in (&self.graph).edge_references() {
            if !ring_bonds.contains(edge.id().index())
                || self.written_order(edge.id()) != BondOrder::Aromatic
            {
                continue;
            }
            if let (Some(&a), Some(&b)) = (index.get(&edge.source()), index.get(&edge.target())) {
                pairing.add_edge(a, b, ());
            }
        }
        let matching = maximum_matching(&pairing);
        if let Some(unpaired) = pairing
            .node_indices()
            .find(|&node| !matching.contains_node(node))
        {
            return Err(format!(
                "{}: its aromatic system has no Kekulé structure; check its hydrogens, \
                 charges and external bonds",
                self.describe(pairing[unpaired])
            ));
        }
        Ok(())
    }

    /// Whether a ring bond from an aromatic atom to a wildcard puts the wildcard in the
    /// written aromatic system: the bond is written aromatic, or the shortest ring through
    /// it is written aromatic apart from wildcards, as in `c1cc*cc1`. A wildcard in a fused
    /// non-aromatic ring, as in `c1ccc2c(c1)C*2`, is not.
    fn wildcard_in_aromatic_system(&self, edge: EdgeIndex, ring_bonds: &FixedBitSet) -> bool {
        self.written_order(edge) == BondOrder::Aromatic
            || self.ring_through(edge, ring_bonds).is_some_and(|ring| {
                ring.iter().all(|&node| {
                    self.graph[node].element.atomic_number == 0 || self.aromatic_as_written(node)
                })
            })
    }

    fn describe(&self, node: NodeIndex) -> String {
        let atom = &self.graph[node];
        format!(
            "atom {} (atomic number {}, charge {})",
            node.index(),
            atom.element.atomic_number,
            atom.charge,
        )
    }

    /// Bond valence from the current bond orders, as the loader computes hydrogens.
    pub(crate) fn bond_valence(&self, node: NodeIndex) -> u32 {
        self.graph
            .edges(node)
            .map(|edge| order_valence(edge.weight().order, edge.weight().donor, node))
            .sum()
    }

    pub(crate) fn record_input(&mut self) {
        self.input_aromatic = self
            .graph
            .node_indices()
            .map(|node| self.graph[node].aromatic)
            .collect();
        self.input_orders = self
            .graph
            .edge_indices()
            .map(|edge| self.graph[edge].order)
            .collect();
    }

    fn written_order(&self, edge: EdgeIndex) -> BondOrder {
        self.input_orders[edge.index()]
    }

    /// Written bonds and hydrogens, an aromatic bond counting one.
    fn written_valence(&self, node: NodeIndex) -> u32 {
        let atom = &self.graph[node];
        self.graph
            .edges(node)
            .map(|edge| order_valence(self.written_order(edge.id()), edge.weight().donor, node))
            .sum::<u32>()
            + u32::from(atom.hs)
            + u32::from(atom.implicit_hs)
    }

    /// Written in lowercase, or with an aromatic bond.
    fn aromatic_as_written(&self, node: NodeIndex) -> bool {
        self.input_aromatic[node.index()]
            || self
                .graph
                .edges(node)
                .any(|edge| self.written_order(edge.id()) == BondOrder::Aromatic)
    }

    /// Whether an aromatic atom needs one of its aromatic bonds to be double, like
    /// benzene's carbons and pyridine's nitrogen. Lone-pair donors (pyrrole's NH, furan's O),
    /// empty p orbitals (boron, tropylium) and atoms with a written double or triple bond
    /// (2-pyridone's C=O carbon) do not.
    fn needs_ring_pi(&self, node: NodeIndex) -> bool {
        if !self.aromatic_as_written(node)
            || self.graph.edges(node).any(|edge| {
                matches!(
                    self.written_order(edge.id()),
                    BondOrder::Double | BondOrder::Triple
                )
            })
        {
            return false;
        }
        let atom = &self.graph[node];
        let donor_or_empty = match (atom.element.atomic_number, atom.charge) {
            // a lone pair with two neighbours: furan's O, pyrrolide's N-
            (8 | 16 | 34 | 52, 0) | (7, -1) => true,
            // with three neighbours, a lone pair (pyrrole's NH, cyclopentadienide's CH-,
            // the S+ of a thiophenium or of CHEMBL8982's [s+]([O-])) or an empty p orbital
            // (borole's B, tropylium's CH+); with two, borabenzene's B and thiopyrylium's
            // S+ need a ring double bond, like pyridine's N
            (5 | 7 | 15 | 33, 0) | (6, -1 | 1) | (16 | 34 | 52, 1) => {
                self.written_valence(node) == 3
            }
            _ => false,
        };
        !donor_or_empty
    }
}

/// The valence a bond of this order gives `node`: a dative bond counts at its acceptor only.
fn order_valence(order: BondOrder, donor: Option<NodeIndex>, node: NodeIndex) -> u32 {
    match order {
        BondOrder::Dative => u32::from(donor != Some(node)),
        BondOrder::Single | BondOrder::Aromatic => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
    }
}

fn maximum_valence(element: u8, charge: i8) -> Result<Option<u8>, ()> {
    // This bounded policy covers common main-group atoms, not metal coordination.
    if !matches!(element, 1 | 5..=9 | 14..=17 | 33..=35 | 53) {
        return Ok(None);
    }
    if element == 1 {
        return match charge {
            0 => Ok(Some(1)),
            -1 | 1 => Ok(Some(0)),
            _ => Err(()),
        };
    }
    if matches!((element, charge), (15 | 33, -2) | (16 | 34, -1)) {
        return Ok(Some(5));
    }
    // Charged atoms use the isoelectronic element's ceiling in the documented validation policy.
    let isoelectronic = i16::from(element) - i16::from(charge);
    let same_period = match element {
        5..=9 => (3..=10).contains(&isoelectronic),
        14..=17 => (11..=18).contains(&isoelectronic),
        33..=35 => (19..=36).contains(&isoelectronic),
        53 => (37..=54).contains(&isoelectronic),
        _ => false,
    };
    if !same_period {
        return Err(());
    }
    let maximum = match isoelectronic {
        10 | 18 | 36 | 54 => 0,
        4 | 8 => 2,
        5 | 7 | 13 | 31 | 49 => 3,
        6 | 14 | 32 | 50 => 4,
        3 | 9 | 17 | 35 => 1,
        15 | 33 | 51 | 53 => 5,
        16 | 34 | 52 => 6,
        _ => return Err(()),
    };
    Ok(Some(maximum))
}
