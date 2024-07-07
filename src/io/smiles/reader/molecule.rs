use crate::io::smiles::reader::atom::parse_atom;
use crate::io::smiles::reader::bond::parse_bond;
use crate::model::bond::{Bond, BondDirection, BondOrder};
use crate::model::molecule::Molecule;
use nom::branch::alt;
use nom::character::complete::{char, digit1};
use nom::combinator::{map, map_res};
use nom::error::ErrorKind;
use nom::multi::many0;
use nom::IResult;
use petgraph::stable_graph::NodeIndex;
use std::collections::HashMap;

fn parse_cycle_digit(input: &str) -> IResult<&str, u8> {
    map_res(digit1, str::parse::<u8>)(input)
}

pub fn parse_molecule(input: &str) -> IResult<&str, Molecule> {
    let mut molecule = Molecule::new();
    let mut open_cycles: HashMap<u8, NodeIndex> = HashMap::new();
    let mut stack: Vec<(NodeIndex, BondOrder)> = Vec::new();
    let mut pending_bonds: Vec<(NodeIndex, NodeIndex, Bond)> = Vec::new();

    let mut parse_atoms_and_bonds = many0(alt((
        map(parse_atom, |atom| (Some(atom), None, None, None)),
        map(parse_bond, |bond| (None, Some(bond), None, None)),
        map(parse_cycle_digit, |digit| (None, None, Some(digit), None)),
        map(char('('), |_| (None, None, None, Some(true))),
        map(char(')'), |_| (None, None, None, Some(false))),
    )));

    let (input, atoms_and_bonds) = parse_atoms_and_bonds(input)?;

    let mut prev_node = NodeIndex::end();
    let mut prev_bond_order = BondOrder::Single;

    for (atom, bond, cycle_digit, open_paren) in atoms_and_bonds {
        if let Some(open) = open_paren {
            if open {
                stack.push((prev_node, prev_bond_order));
            } else {
                if stack.is_empty() {
                    return Err(nom::Err::Failure(nom::error::Error::new(
                        input,
                        ErrorKind::Verify,
                    )));
                }

                let (node, bond) = stack.pop().unwrap();
                prev_node = node;
                prev_bond_order = bond;
            }
        } else if let Some(digit) = cycle_digit {
            if let Some(open_node) = open_cycles.remove(&digit) {
                pending_bonds.push((prev_node, open_node, Bond { order: prev_bond_order, direction: BondDirection::Unspecified }));
                prev_bond_order = BondOrder::Single;
            } else {
                open_cycles.insert(digit, prev_node);
            }
        } else {
            if let Some(atom) = atom {
                let node = molecule.add_atom(atom);
                if prev_node != NodeIndex::end() && bond.is_none() {
                    pending_bonds.push((prev_node, node, Bond { order: prev_bond_order, direction: BondDirection::Unspecified }));
                }
                prev_node = node;
            }

            if let Some(bond) = bond {
                prev_bond_order = bond.order;
            } else {
                prev_bond_order = BondOrder::Single;
            }
        }
    }

    if !open_cycles.is_empty() || !stack.is_empty() {
        return Err(nom::Err::Failure(nom::error::Error::new(
            input,
            ErrorKind::Verify,
        )));
    }

    for (node1, node2, bond) in pending_bonds {
        molecule.add_bond(node1, node2, bond);
    }

    Ok((input, molecule))
}

#[cfg(test)]
mod tests {
    use crate::io::smiles::reader::molecule::parse_molecule;
    use crate::model::bond::BondOrder;
    use petgraph::stable_graph::{EdgeIndex, NodeIndex};
    use crate::model::chirality::Chirality;

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
        assert_eq!(m.get_bond_by_atoms(NodeIndex::new(0), NodeIndex::new(1))
                       .unwrap()
                       .order, BondOrder::Single);
        assert!(m.has_bond(NodeIndex::new(1), NodeIndex::new(2)));
        assert_eq!(m.get_bond_by_atoms(NodeIndex::new(1), NodeIndex::new(2))
                       .unwrap()
                       .order, BondOrder::Double);
        assert!(m.has_bond(NodeIndex::new(1), NodeIndex::new(3)));
        assert_eq!(m.get_bond_by_atoms(NodeIndex::new(1), NodeIndex::new(3))
                       .unwrap()
                       .order, BondOrder::Single);
        assert!(m.has_bond(NodeIndex::new(3), NodeIndex::new(4)));
        assert_eq!(m.get_bond_by_atoms(NodeIndex::new(3), NodeIndex::new(4))
                       .unwrap()
                       .order, BondOrder::Single);
        assert!(m.has_bond(NodeIndex::new(3), NodeIndex::new(0)));
        assert_eq!(m.get_bond_by_atoms(NodeIndex::new(3), NodeIndex::new(0))
                       .unwrap()
                       .order, BondOrder::Double);
    }

    #[test]
    fn test_parse_molecule_c1cc() {
        assert!(parse_molecule("C1CC").is_err())
    }


    #[test]
    fn test_parse_molecule_chiral_anticlockwise() {
        assert!(!parse_molecule("N[C@H](C)C(=O)O").is_err());
    }

    #[test]
    fn test_parse_molecule_chiral_clockwise() {
        assert!(!parse_molecule("N[C@@](F)(C)C(=O)O").is_err());
    }

    #[test]
    fn test_molecule_weight() {
        let m = parse_molecule("C([H])([H])([H])(H)").unwrap().1;
        assert_eq!(m.weight(), 16.0423);
    }

    #[test]
    fn test_molecule_chiral() {
        let m = parse_molecule("N[C@H](C)C(=O)O").unwrap().1;
        assert_eq!(m.get_atom(NodeIndex::new(1)).unwrap().chirality, Chirality::Anticlockwise);
    }
}
