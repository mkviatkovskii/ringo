use crate::model::bond::BondOrder;
use crate::model::molecule::Molecule;
use petgraph::stable_graph::NodeIndex;
use winnow::error::{ContextError, ErrMode};
use winnow::ModalResult;

pub(super) fn bond_valence(molecule: &Molecule, node: NodeIndex) -> ModalResult<u8> {
    molecule.graph.edges(node).try_fold(0u8, |sum, edge| {
        let order = match edge.weight().order {
            BondOrder::Dative => u8::from(edge.weight().donor != Some(node)),
            BondOrder::Single | BondOrder::Aromatic => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
        };
        sum.checked_add(order).ok_or_else(error)
    })
}

/// Check the written valence before aromaticity can replace double bonds.
/// Arithmetic limits apply even when chemical validation is disabled.
pub(super) fn check_valence(molecule: &Molecule, validate: bool) -> ModalResult<()> {
    for node in molecule.graph.node_indices() {
        let atom = &molecule.graph[node];
        let valence = bond_valence(molecule, node)?
            .checked_add(atom.hs)
            .and_then(|n| n.checked_add(atom.implicit_hs))
            .ok_or_else(error)?;
        if !validate {
            continue;
        }
        let Some(maximum) = maximum_valence(atom.element.atomic_number, atom.charge)? else {
            continue;
        };
        let aromatic = molecule
            .graph
            .edges(node)
            .any(|e| e.weight().order == BondOrder::Aromatic);
        let has_pi_bond = molecule
            .graph
            .edges(node)
            .any(|e| matches!(e.weight().order, BondOrder::Double | BondOrder::Triple));
        // Lone-pair donors and empty p orbitals do not need an extra ring bond order.
        let donor_or_empty = match (atom.element.atomic_number, atom.charge) {
            (5 | 8 | 16 | 34, 0) | (6, -1 | 1) | (7, -1) => true,
            (7 | 15 | 33, 0) => valence == 3,
            _ => false,
        };
        let pi = u8::from(aromatic && !has_pi_bond && !donor_or_empty);
        if valence.checked_add(pi).ok_or_else(error)? > maximum {
            return Err(error());
        }
    }
    Ok(())
}

fn maximum_valence(element: u8, charge: i8) -> ModalResult<Option<u8>> {
    // This bounded policy covers common main-group atoms, not metal coordination.
    if !matches!(element, 1 | 5..=9 | 14..=17 | 33..=35 | 53) {
        return Ok(None);
    }
    if element == 1 {
        return match charge {
            0 => Ok(Some(1)),
            -1 | 1 => Ok(Some(0)),
            _ => Err(error()),
        };
    }
    if !(-2..=2).contains(&charge) {
        return Err(error());
    }
    if matches!((element, charge), (15 | 33, -2) | (16 | 34, -1)) {
        return Ok(Some(5));
    }
    // Charged atoms use the isoelectronic element's ceiling (RDKit valence policy).
    let isoelectronic = i16::from(element) - i16::from(charge);
    let maximum = match isoelectronic {
        2 | 10 | 18 | 36 => 0,
        4 | 8 => 2,
        5 | 7 | 13 | 31 => 3,
        6 | 14 | 32 => 4,
        9 | 17 | 35 => 1,
        15 | 33 | 53 => 5,
        16 | 34 | 52 | 54 => 6,
        _ => return Err(error()),
    };
    Ok(Some(maximum))
}

fn error() -> ErrMode<ContextError> {
    ErrMode::Cut(ContextError::new())
}
