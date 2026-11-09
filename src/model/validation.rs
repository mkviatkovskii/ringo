use crate::model::bond::BondOrder;
use crate::model::molecule::Molecule;
use petgraph::stable_graph::NodeIndex;

impl Molecule {
    /// Runs optional post-load chemical checks without modifying the molecule.
    /// Currently checks main-group valence limits, not aromaticity or stereogenicity.
    /// Errors identify the zero-based atom index and the unsupported charge or excessive valence.
    pub fn validate(&self) -> Result<(), String> {
        for node in self.graph.node_indices() {
            let atom = &self.graph[node];
            let context = || {
                format!(
                    "atom {} (atomic number {}, charge {})",
                    node.index(),
                    atom.element.atomic_number,
                    atom.charge,
                )
            };
            let Some(maximum) = maximum_valence(atom.element.atomic_number, atom.charge)
                .map_err(|()| format!("{}: unsupported charge state", context()))?
            else {
                continue;
            };
            let valence = self.input_valences[node.index()];
            if valence > u16::from(maximum) {
                return Err(format!(
                    "{}: valence {valence} exceeds maximum {maximum}",
                    context()
                ));
            }
        }
        Ok(())
    }

    pub(crate) fn bond_valence(&self, node: NodeIndex) -> Option<u8> {
        self.graph.edges(node).try_fold(0u8, |sum, edge| {
            let order = match edge.weight().order {
                BondOrder::Dative => u8::from(edge.weight().donor != Some(node)),
                BondOrder::Single | BondOrder::Aromatic => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
            };
            sum.checked_add(order)
        })
    }

    pub(crate) fn record_valences(&mut self) -> Option<()> {
        self.input_valences = self
            .graph
            .node_indices()
            .map(|node| {
                let atom = &self.graph[node];
                let valence = self
                    .bond_valence(node)?
                    .checked_add(atom.hs)?
                    .checked_add(atom.implicit_hs)?;
                let aromatic = self
                    .graph
                    .edges(node)
                    .any(|e| e.weight().order == BondOrder::Aromatic);
                let has_pi_bond = self
                    .graph
                    .edges(node)
                    .any(|e| matches!(e.weight().order, BondOrder::Double | BondOrder::Triple));
                // Lone-pair donors and empty p orbitals do not need an extra ring bond order.
                let donor_or_empty = match (atom.element.atomic_number, atom.charge) {
                    (5 | 8 | 16 | 34, 0) | (6, -1 | 1) | (7, -1) => true,
                    (7 | 15 | 33, 0) => valence == 3,
                    _ => false,
                };
                let pi = u16::from(aromatic && !has_pi_bond && !donor_or_empty);
                // The written total fits u8; its estimated aromatic contribution can add one.
                Some(u16::from(valence) + pi)
            })
            .collect::<Option<_>>()?;
        Some(())
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
    // Charged atoms use the isoelectronic element's ceiling (RDKit valence policy).
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
        10 | 18 | 36 => 0,
        4 | 8 => 2,
        5 | 7 | 13 | 31 | 49 => 3,
        6 | 14 | 32 | 50 => 4,
        9 | 17 | 35 => 1,
        15 | 33 | 51 | 53 => 5,
        16 | 34 | 52 | 54 => 6,
        _ => return Err(()),
    };
    Ok(Some(maximum))
}
