use crate::model::bond::BondOrder;
use crate::model::molecule::Molecule;
use crate::model::stereo::{AlleneStereo, Ligand};
use petgraph::stable_graph::NodeIndex;
use winnow::error::{ContextError, ErrMode};
use winnow::ModalResult;

pub(super) fn resolve_allene(
    molecule: &Molecule,
    center: NodeIndex,
    permutation: u8,
    ligand_order: &[Vec<Option<NodeIndex>>],
    has_previous: &[bool],
) -> ModalResult<AlleneStereo> {
    let error = || ErrMode::Cut(ContextError::new());
    let atom = &molecule.graph[center];
    let neighbors = &ligand_order[center.index()];
    if neighbors.len() != 2 || atom.hs != 0 || atom.implicit_hs != 0 || atom.aromatic {
        return Err(error());
    }
    let mut seen = vec![false; ligand_order.len()];
    seen[center.index()] = true;
    let mut arms = [Vec::new(), Vec::new()];
    let mut ligands = Vec::new();
    for (side, neighbor) in neighbors.iter().enumerate() {
        let mut previous = center;
        let mut current = neighbor.ok_or_else(error)?;
        loop {
            if seen[current.index()]
                || molecule
                    .get_bond_by_atoms(previous, current)
                    .is_none_or(|bond| bond.order != BondOrder::Double)
            {
                return Err(error());
            }
            seen[current.index()] = true;
            arms[side].push(current);
            let atom = &molecule.graph[current];
            let hydrogens = usize::from(atom.hs) + usize::from(atom.implicit_hs);
            let others: Vec<NodeIndex> = ligand_order[current.index()]
                .iter()
                .filter(|node| **node != Some(previous))
                .map(|node| node.ok_or_else(error))
                .collect::<ModalResult<_>>()?;
            let is_double = |node| {
                molecule
                    .get_bond_by_atoms(current, node)
                    .is_some_and(|bond| bond.order == BondOrder::Double)
            };
            if others.iter().copied().any(is_double) {
                if others.len() != 1 || hydrogens != 0 || atom.aromatic {
                    return Err(error());
                }
                previous = current;
                current = others[0];
                continue;
            }
            if others.len() + hydrogens != 2
                || atom.aromatic
                || others.iter().any(|&node| {
                    molecule
                        .get_bond_by_atoms(current, node)
                        .is_none_or(|bond| bond.order != BondOrder::Single)
                })
            {
                return Err(error());
            }
            let mut terminal: Vec<_> = others.into_iter().map(Ligand::Atom).collect();
            // H follows the terminal atom; removing an axis predecessor moves H to slot zero.
            let position = usize::from(
                has_previous[current.index()] && ligand_order[current.index()][0] != Some(previous),
            );
            for _ in 0..hydrogens {
                terminal.insert(position, Ligand::Hydrogen);
            }
            ligands.extend(terminal);
            break;
        }
    }
    // Only an even number of double bonds with the tag at the midpoint is allene-like.
    if arms[0].len() != arms[1].len() {
        return Err(error());
    }
    let mut axis = std::mem::take(&mut arms[0]);
    axis.reverse();
    axis.push(center);
    axis.append(&mut arms[1]);
    Ok(AlleneStereo {
        center,
        axis,
        ligands: ligands.try_into().map_err(|_| error())?,
        permutation,
    })
}
