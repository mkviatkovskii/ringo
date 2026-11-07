use crate::io::smiles::reader::atom::parse_atom;
use crate::io::smiles::reader::bond::parse_bond;
use crate::model::atom::Atom;
use crate::model::bond::{Bond, BondOrder};
use crate::model::molecule::Molecule;
use petgraph::stable_graph::NodeIndex;
use std::collections::HashMap;
use winnow::ascii::digit1;
use winnow::combinator::{alt, repeat};
use winnow::error::{ContextError, ErrMode};
use winnow::prelude::*;

// the atom carries whether it was written in brackets
type Token = (Option<(Atom, bool)>, Option<Bond>, Option<u8>, Option<bool>);

fn parse_cycle_digit(input: &mut &str) -> ModalResult<u8> {
    digit1.try_map(str::parse::<u8>).parse_next(input)
}

/// Hydrogens an atom outside brackets gets: the lowest normal valence that
/// fits its bonds, minus the bonds (OpenSMILES organic subset rules).
fn implicit_hydrogens(molecule: &Molecule, node: NodeIndex, atomic_number: u8) -> u8 {
    let valences: &[u8] = match atomic_number {
        5 => &[3],
        6 => &[4],
        7 => &[3, 5],
        8 => &[2],
        15 => &[3, 5],
        16 => &[2, 4, 6],
        1 | 9 | 17 | 35 | 53 => &[1],
        _ => &[],
    };
    let mut bonds: u8 = 0;
    let mut aromatic = false;
    for edge in molecule.graph.edges(node) {
        bonds += match edge.weight().order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => {
                aromatic = true;
                1
            }
        };
    }
    // an aromatic atom also spends one electron on the ring
    if aromatic {
        bonds += 1;
    }
    valences
        .iter()
        .find(|&&valence| valence >= bonds)
        .map_or(0, |valence| valence - bonds)
}

/// Parses a molecule from the start of `input`; returns the unparsed rest and the molecule.
pub fn parse_molecule(input: &str) -> ModalResult<(&str, Molecule)> {
    parse_molecule_next.parse_peek(input)
}

fn parse_molecule_next(input: &mut &str) -> ModalResult<Molecule> {
    let mut molecule = Molecule::new();
    let mut open_cycles: HashMap<u8, NodeIndex> = HashMap::new();
    let mut stack: Vec<(NodeIndex, BondOrder)> = Vec::new();
    let mut pending_bonds: Vec<(NodeIndex, NodeIndex, Bond)> = Vec::new();
    let mut organic_nodes: Vec<(NodeIndex, u8)> = Vec::new();

    let atoms_and_bonds: Vec<Token> = repeat(
        0..,
        alt((
            parse_atom
                .with_taken()
                .map(|(atom, taken)| (Some((atom, taken.starts_with('['))), None, None, None)),
            parse_bond.map(|bond| (None, Some(bond), None, None)),
            parse_cycle_digit.map(|digit| (None, None, Some(digit), None)),
            '('.map(|_| (None, None, None, Some(true))),
            ')'.map(|_| (None, None, None, Some(false))),
        )),
    )
    .parse_next(input)?;

    let mut prev_node = NodeIndex::end();
    let mut prev_bond = BondOrder::Single;

    for (atom, bond, cycle_digit, open_paren) in atoms_and_bonds {
        if let Some(open) = open_paren {
            if open {
                stack.push((prev_node, prev_bond));
            } else {
                let Some((node, bond)) = stack.pop() else {
                    return Err(ErrMode::Cut(ContextError::new()));
                };
                prev_node = node;
                prev_bond = bond;
            }
        } else if let Some(digit) = cycle_digit {
            // a ring bond needs an atom to start from or end at
            if prev_node == NodeIndex::end() {
                return Err(ErrMode::Cut(ContextError::new()));
            }
            if let Some(open_node) = open_cycles.remove(&digit) {
                pending_bonds.push((prev_node, open_node, Bond { order: prev_bond }));
                prev_bond = BondOrder::Single;
            } else {
                open_cycles.insert(digit, prev_node);
            }
        } else {
            if let Some((atom, bracket)) = atom {
                let atomic_number = atom.element.atomic_number;
                let node = molecule.add_atom(atom);
                if !bracket {
                    organic_nodes.push((node, atomic_number));
                }
                if prev_node != NodeIndex::end() && bond.is_none() {
                    pending_bonds.push((prev_node, node, Bond { order: prev_bond }));
                }
                prev_node = node;
            }

            if let Some(bond) = bond {
                prev_bond = bond.order;
            } else {
                prev_bond = BondOrder::Single;
            }
        }
    }

    if !open_cycles.is_empty() || !stack.is_empty() {
        return Err(ErrMode::Cut(ContextError::new()));
    }

    for (node1, node2, bond) in pending_bonds {
        molecule.add_bond(node1, node2, bond);
    }

    for (node, atomic_number) in organic_nodes {
        let hs = implicit_hydrogens(&molecule, node, atomic_number);
        if let Some(atom) = molecule.graph.node_weight_mut(node) {
            atom.implicit_hs = hs;
        }
    }

    Ok(molecule)
}

#[cfg(test)]
mod tests {
    use crate::io::smiles::reader::molecule::parse_molecule;
    use crate::model::bond::BondOrder;
    use petgraph::stable_graph::{EdgeIndex, NodeIndex};

    #[test]
    fn parse_molecule_empty() {
        let m = parse_molecule("").unwrap().1;
        assert_eq!(m.count_atoms(), 0);
        assert_eq!(m.count_bonds(), 0);
    }

    #[test]
    fn parse_molecule_c() {
        let m = parse_molecule("C").unwrap().1;
        assert_eq!(m.count_atoms(), 1);
        assert_eq!(m.count_bonds(), 0);
        assert_eq!(
            m.get_atom(NodeIndex::new(0)).unwrap().element.atomic_number,
            6
        );
    }

    #[test]
    fn parse_molecule_cc() {
        let m = parse_molecule("CN").unwrap().1;
        assert_eq!(m.count_atoms(), 2);
        assert_eq!(m.count_bonds(), 1);
        assert_eq!(
            m.get_atom(NodeIndex::new(0)).unwrap().element.atomic_number,
            6
        );
        assert_eq!(
            m.get_atom(NodeIndex::new(1)).unwrap().element.atomic_number,
            7
        );
        assert_eq!(
            m.get_bond(EdgeIndex::new(0)).unwrap().order,
            BondOrder::Single
        );
    }

    #[test]
    fn parse_molecule_branch() {
        let m = parse_molecule("C(O)N").unwrap().1;
        assert_eq!(m.count_atoms(), 3);
        assert_eq!(m.count_bonds(), 2);
        assert_eq!(
            m.get_atom(NodeIndex::new(0)).unwrap().element.atomic_number,
            6
        );
        assert_eq!(
            m.get_atom(NodeIndex::new(1)).unwrap().element.atomic_number,
            8
        );
        assert_eq!(
            m.get_atom(NodeIndex::new(2)).unwrap().element.atomic_number,
            7
        );
        assert_eq!(
            m.get_bond(EdgeIndex::new(0)).unwrap().order,
            BondOrder::Single
        );
        assert_eq!(
            m.get_bond(EdgeIndex::new(1)).unwrap().order,
            BondOrder::Single
        );
    }

    #[test]
    fn parse_molecule_c1cc1() {
        let m = parse_molecule("C1P=N#1").unwrap().1;
        assert_eq!(m.count_atoms(), 3);
        assert_eq!(m.count_bonds(), 3);
        assert_eq!(
            m.get_bond(EdgeIndex::new(0)).unwrap().order,
            BondOrder::Single
        );
        assert_eq!(
            m.get_bond(EdgeIndex::new(1)).unwrap().order,
            BondOrder::Double
        );
        assert_eq!(
            m.get_bond(EdgeIndex::new(2)).unwrap().order,
            BondOrder::Triple
        );
        let atom_0_bonds = m.get_bonds_for_atom(NodeIndex::new(0));
        assert_eq!(atom_0_bonds.len(), 2);
        assert_eq!(
            m.get_bond(atom_0_bonds[0]).unwrap().order,
            BondOrder::Single
        );
        assert_eq!(
            m.get_bond(atom_0_bonds[1]).unwrap().order,
            BondOrder::Triple
        );

        let mut atom_0_neighbors = m.get_neighbors_for_atom(NodeIndex::new(0));
        assert_eq!(atom_0_neighbors.len(), 2);
        assert_eq!(
            m.get_atom(atom_0_neighbors.pop_first().unwrap())
                .unwrap()
                .element
                .atomic_number,
            15
        );
        assert_eq!(
            m.get_atom(atom_0_neighbors.pop_first().unwrap())
                .unwrap()
                .element
                .atomic_number,
            7
        );
    }

    #[test]
    fn parse_molecule_branch_double_bond() {
        let m = parse_molecule("C(=O)N").unwrap().1;
        assert_eq!(m.count_atoms(), 3);
        assert_eq!(m.count_bonds(), 2);
        assert_eq!(
            m.get_atom(NodeIndex::new(0)).unwrap().element.atomic_number,
            6
        );
        assert_eq!(
            m.get_atom(NodeIndex::new(1)).unwrap().element.atomic_number,
            8
        );
        assert_eq!(
            m.get_atom(NodeIndex::new(2)).unwrap().element.atomic_number,
            7
        );
        assert_eq!(
            m.get_bond_by_atoms(NodeIndex::new(0), NodeIndex::new(1))
                .unwrap()
                .order,
            BondOrder::Double
        );
        assert_eq!(
            m.get_bond_by_atoms(NodeIndex::new(0), NodeIndex::new(2))
                .unwrap()
                .order,
            BondOrder::Single
        );
    }

    #[test]
    fn parse_molecule_branch_double_bonds() {
        let m = parse_molecule("C(=O)=N").unwrap().1;
        assert_eq!(m.count_atoms(), 3);
        assert_eq!(m.count_bonds(), 2);
        assert_eq!(
            m.get_atom(NodeIndex::new(0)).unwrap().element.atomic_number,
            6
        );
        assert_eq!(
            m.get_atom(NodeIndex::new(1)).unwrap().element.atomic_number,
            8
        );
        assert_eq!(
            m.get_atom(NodeIndex::new(2)).unwrap().element.atomic_number,
            7
        );
        assert_eq!(
            m.get_bond_by_atoms(NodeIndex::new(0), NodeIndex::new(1))
                .unwrap()
                .order,
            BondOrder::Double
        );
        assert_eq!(
            m.get_bond_by_atoms(NodeIndex::new(0), NodeIndex::new(2))
                .unwrap()
                .order,
            BondOrder::Double
        );
    }

    #[test]
    fn parse_molecule_branch_recursive() {
        let m = parse_molecule("C(=S(=O)P)N").unwrap().1;
        assert_eq!(m.count_atoms(), 5);
        assert_eq!(m.count_bonds(), 4);
        assert_eq!(
            m.get_atom(NodeIndex::new(0)).unwrap().element.atomic_number,
            6
        );
        assert_eq!(
            m.get_atom(NodeIndex::new(1)).unwrap().element.atomic_number,
            16
        );
        assert_eq!(
            m.get_atom(NodeIndex::new(2)).unwrap().element.atomic_number,
            8
        );
        assert_eq!(
            m.get_atom(NodeIndex::new(3)).unwrap().element.atomic_number,
            15
        );
        assert_eq!(
            m.get_atom(NodeIndex::new(4)).unwrap().element.atomic_number,
            7
        );
        assert_eq!(
            m.get_bond_by_atoms(NodeIndex::new(0), NodeIndex::new(1))
                .unwrap()
                .order,
            BondOrder::Double
        );
        assert_eq!(
            m.get_bond_by_atoms(NodeIndex::new(1), NodeIndex::new(2))
                .unwrap()
                .order,
            BondOrder::Double
        );
        assert_eq!(
            m.get_bond_by_atoms(NodeIndex::new(1), NodeIndex::new(3))
                .unwrap()
                .order,
            BondOrder::Single
        );
        assert_eq!(
            m.get_bond_by_atoms(NodeIndex::new(0), NodeIndex::new(4))
                .unwrap()
                .order,
            BondOrder::Single
        );
    }

    #[test]
    fn parse_molecule_cycle_double() {
        let m = parse_molecule("N1OC=1S").unwrap().1;
        assert_eq!(m.count_atoms(), 4);
        assert_eq!(m.count_bonds(), 4);
        assert_eq!(
            m.get_atom(NodeIndex::new(0)).unwrap().element.atomic_number,
            7
        );
        assert_eq!(
            m.get_atom(NodeIndex::new(1)).unwrap().element.atomic_number,
            8
        );
        assert_eq!(
            m.get_atom(NodeIndex::new(2)).unwrap().element.atomic_number,
            6
        );
        assert_eq!(
            m.get_atom(NodeIndex::new(3)).unwrap().element.atomic_number,
            16
        );
        assert!(m.has_bond(NodeIndex::new(0), NodeIndex::new(1)));
        assert!(
            m.get_bond_by_atoms(NodeIndex::new(0), NodeIndex::new(1))
                .unwrap()
                .order
                == BondOrder::Single
        );
        assert!(m.has_bond(NodeIndex::new(1), NodeIndex::new(2)));
        assert!(
            m.get_bond_by_atoms(NodeIndex::new(1), NodeIndex::new(2))
                .unwrap()
                .order
                == BondOrder::Single
        );
        assert!(m.has_bond(NodeIndex::new(0), NodeIndex::new(2)));
        assert!(
            m.get_bond_by_atoms(NodeIndex::new(0), NodeIndex::new(2))
                .unwrap()
                .order
                == BondOrder::Double
        );
        assert!(m.has_bond(NodeIndex::new(2), NodeIndex::new(3)));
        assert!(
            m.get_bond_by_atoms(NodeIndex::new(2), NodeIndex::new(3))
                .unwrap()
                .order
                == BondOrder::Single
        );
    }

    #[test]
    fn parse_molecule_cycle_branch() {
        let m = parse_molecule("N1C(=P)S=1O").unwrap().1;
        assert_eq!(m.count_atoms(), 5);
        assert_eq!(m.count_bonds(), 5);
        assert_eq!(
            m.get_atom(NodeIndex::new(0)).unwrap().element.atomic_number,
            7
        );
        assert_eq!(
            m.get_atom(NodeIndex::new(1)).unwrap().element.atomic_number,
            6
        );
        assert_eq!(
            m.get_atom(NodeIndex::new(2)).unwrap().element.atomic_number,
            15
        );
        assert_eq!(
            m.get_atom(NodeIndex::new(3)).unwrap().element.atomic_number,
            16
        );
        assert_eq!(
            m.get_atom(NodeIndex::new(4)).unwrap().element.atomic_number,
            8
        );
        assert!(m.has_bond(NodeIndex::new(0), NodeIndex::new(1)));
        assert!(
            m.get_bond_by_atoms(NodeIndex::new(0), NodeIndex::new(1))
                .unwrap()
                .order
                == BondOrder::Single
        );
        assert!(m.has_bond(NodeIndex::new(1), NodeIndex::new(2)));
        assert!(
            m.get_bond_by_atoms(NodeIndex::new(1), NodeIndex::new(2))
                .unwrap()
                .order
                == BondOrder::Double
        );
        assert!(m.has_bond(NodeIndex::new(1), NodeIndex::new(3)));
        assert!(
            m.get_bond_by_atoms(NodeIndex::new(1), NodeIndex::new(3))
                .unwrap()
                .order
                == BondOrder::Single
        );
        assert!(m.has_bond(NodeIndex::new(3), NodeIndex::new(4)));
        assert!(
            m.get_bond_by_atoms(NodeIndex::new(3), NodeIndex::new(4))
                .unwrap()
                .order
                == BondOrder::Single
        );
        assert!(m.has_bond(NodeIndex::new(3), NodeIndex::new(0)));
        assert!(
            m.get_bond_by_atoms(NodeIndex::new(3), NodeIndex::new(0))
                .unwrap()
                .order
                == BondOrder::Double
        );
    }

    #[test]
    fn parse_molecule_c1cc() {
        assert!(parse_molecule("C1CC").is_err())
    }

    #[test]
    fn parse_molecule_ring_digit_without_atom() {
        // used to panic in add_edge on the placeholder node
        assert!(parse_molecule("1C1").is_err());
        assert!(parse_molecule("(C)1C1").is_err());
    }

    #[test]
    fn molecule_weight() {
        // methane with explicit, bracket and implicit hydrogens
        for smiles in ["C([H])([H])([H])(H)", "[CH4]", "C"] {
            let m = parse_molecule(smiles).unwrap().1;
            assert!((m.weight() - 16.043).abs() < 1e-9, "{}", smiles);
        }
        // ethanol C2H6O
        let m = parse_molecule("CCO").unwrap().1;
        assert!((m.weight() - 46.069).abs() < 1e-9);
        // caffeine C8H10N4O2
        let m = parse_molecule("CN1C=NC2=C1C(=O)N(C(=O)N2C)C").unwrap().1;
        assert!((m.weight() - 194.194).abs() < 1e-9);
        // ammonium NH4 and the methylene radical CH2
        let m = parse_molecule("[NH4+]").unwrap().1;
        assert!((m.weight() - 18.039).abs() < 1e-9);
        let m = parse_molecule("[CH2]").unwrap().1;
        assert!((m.weight() - 14.027).abs() < 1e-9);
    }

    #[test]
    fn implicit_hydrogens() {
        let hs = |smiles: &str| -> Vec<u8> {
            let m = parse_molecule(smiles).unwrap().1;
            m.graph
                .node_weights()
                .map(|atom| atom.implicit_hs)
                .collect()
        };
        assert_eq!(hs("C"), [4]);
        assert_eq!(hs("C=O"), [2, 0]);
        assert_eq!(hs("C#N"), [1, 0]);
        assert_eq!(hs("CCl"), [3, 0]);
        // bracket atoms take only their written hydrogens
        assert_eq!(hs("[CH3]C"), [0, 3]);
        // charges and radicals need brackets, so valence never adds hydrogens to them
        assert_eq!(hs("[NH4+]"), [0]);
        assert_eq!(hs("C[O-]"), [3, 0]);
        assert_eq!(hs("[CH2]"), [0]);
        assert_eq!(hs("[CH2]C"), [0, 3]);
        // next valence when the first one is exceeded
        assert_eq!(hs("CS(=O)C"), [3, 0, 0, 3]);
        assert_eq!(hs("OP(=O)(O)O"), [1, 0, 0, 1, 1]);
        // aromatic bonds spend an extra electron
        assert_eq!(hs("C1:C:C:C:C:C1"), [1, 1, 1, 1, 1, 1]);
    }

    #[test]
    fn ecfp() {
        for smiles in ["C", "CC", "C1COCNC1", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"] {
            println!("{}: ", smiles);
            let m = parse_molecule(smiles).unwrap().1;
            let result = m.ecfp(2, 512);
            println!("{:?}", result);
        }
    }
}
