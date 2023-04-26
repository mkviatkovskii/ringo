use crate::molecule::model::bond::{Bond, BondOrder};
use crate::molecule::model::molecule::Molecule;
use crate::molecule::smiles::reader::atom::parse_atom;
use crate::molecule::smiles::reader::bond::parse_bond;
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

fn parse_molecule_old(input: &str) -> IResult<&str, Molecule> {
    let mut molecule = Molecule::new();
    let mut open_cycles: HashMap<u8, NodeIndex> = HashMap::new();

    let mut parse_atoms_and_bonds = many0(alt((
        map(parse_atom, |atom| (Some(atom), None, None)),
        map(parse_bond, |bond| (None, Some(bond), None)),
        map(parse_cycle_digit, |digit| (None, None, Some(digit))),
    )));

    let (input, atoms_and_bonds) = parse_atoms_and_bonds(input)?;

    let mut prev_node = NodeIndex::end();
    let mut prev_bond = BondOrder::Single;

    let mut pending_bonds: Vec<(NodeIndex, NodeIndex, Bond)> = Vec::new();

    for (atom, bond, cycle_digit) in atoms_and_bonds {
        if let Some(digit) = cycle_digit {
            if let Some(open_node) = open_cycles.remove(&digit) {
                pending_bonds.push((prev_node, open_node, Bond { order: prev_bond }));
            } else {
                open_cycles.insert(digit, prev_node);
            }
        } else {
            if let Some(atom) = atom {
                let node = molecule.add_atom(atom);
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

    if !open_cycles.is_empty() {
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

fn parse_molecule(input: &str) -> IResult<&str, Molecule> {
    let mut molecule = Molecule::new();
    let mut open_cycles: HashMap<u8, NodeIndex> = HashMap::new();
    let mut stack: Vec<(NodeIndex, BondOrder)> = Vec::new();

    let mut parse_atoms_and_bonds = many0(alt((
        map(parse_atom, |atom| (Some(atom), None, None, None)),
        map(parse_bond, |bond| (None, Some(bond), None, None)),
        map(parse_cycle_digit, |digit| (None, None, Some(digit), None)),
        map(char('('), |_| (None, None, None, Some(true))),
        map(char(')'), |_| (None, None, None, Some(false))),
    )));

    let (input, atoms_and_bonds) = parse_atoms_and_bonds(input)?;

    let mut prev_node = NodeIndex::end();
    let mut prev_bond = BondOrder::Single;

    let mut pending_bonds: Vec<(NodeIndex, NodeIndex, Bond)> = Vec::new();

    for (atom, bond, cycle_digit, open_paren) in atoms_and_bonds {
        if let Some(open) = open_paren {
            if open {
                stack.push((prev_node, prev_bond));
            } else {
                if stack.is_empty() {
                    return Err(nom::Err::Failure(nom::error::Error::new(
                        input,
                        ErrorKind::Verify,
                    )));
                }

                let (node, bond) = stack.pop().unwrap();
                prev_node = node;
                prev_bond = bond;
            }
        } else if let Some(digit) = cycle_digit {
            if let Some(open_node) = open_cycles.remove(&digit) {
                pending_bonds.push((prev_node, open_node, Bond { order: prev_bond }));
            } else {
                open_cycles.insert(digit, prev_node);
            }
        } else {
            if let Some(atom) = atom {
                let node = molecule.add_atom(atom);
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
    use crate::molecule::model::bond::BondOrder;
    use crate::molecule::smiles::reader::molecule::parse_molecule;
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
    fn parse_molecule_c1cc() {
        assert!(parse_molecule("C1CC").is_err())
    }

    #[test]
    fn molecule_weight() {
        let m = parse_molecule("C([H])([H])([H])(H)").unwrap().1;
        assert_eq!(m.weight(), 16.0423);
    }

    #[test]
    fn ecfp() {
        for smiles in ["C", "CC", "C1COCNC1", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"] {
            println!("{}: ", smiles);
            let m = parse_molecule(smiles).unwrap().1;
            let result = m.ecfp(2, 128);
            for bit in result {
                print!("{}", if bit { 1 } else { 0 });
            }
            println!("");
        }
    }
}
