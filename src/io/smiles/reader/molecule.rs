use crate::io::smiles::reader::atom::parse_atom;
use crate::io::smiles::reader::bond::parse_bond;
use crate::model::atom::Atom;
use crate::model::bond::{Bond, BondOrder};
use crate::model::molecule::Molecule;
use petgraph::stable_graph::NodeIndex;
use std::collections::HashMap;
use winnow::combinator::{alt, preceded, repeat};
use winnow::error::{ContextError, ErrMode};
use winnow::prelude::*;
use winnow::stream::AsChar;
use winnow::token::take_while;

enum Token {
    /// The atom, and whether it was written in brackets.
    Atom(Atom, bool),
    Bond(BondOrder),
    Ring(u8),
    Open,
    Close,
    Dot,
}

// a bond with no symbol is None until both its atoms are known
type PendingBond = (NodeIndex, NodeIndex, Option<BondOrder>);

/// A ring-closure number: one digit, or `%` and exactly two digits, so `12` is rings 1 and 2.
fn parse_cycle_digit(input: &mut &str) -> ModalResult<u8> {
    alt((
        preceded('%', take_while(2, AsChar::is_dec_digit)),
        take_while(1, AsChar::is_dec_digit),
    ))
    .try_map(str::parse::<u8>)
    .parse_next(input)
}

/// Hydrogens an atom outside brackets gets: the lowest normal valence that
/// fits its bonds, minus the bonds (`OpenSMILES` organic subset rules).
fn implicit_hydrogens(
    molecule: &Molecule,
    node: NodeIndex,
    atomic_number: u8,
    aromatic: bool,
) -> u8 {
    let mut valences: &[u8] = match atomic_number {
        5 => &[3],
        6 => &[4],
        7 | 15 => &[3, 5],
        8 => &[2],
        16 => &[2, 4, 6],
        1 | 9 | 17 | 35 | 53 => &[1],
        _ => &[],
    };
    let mut bonds: u8 = 0;
    let mut aromatic = aromatic;
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
    // an aromatic atom also spends one electron on the ring, and keeps its lowest
    // valence: the s in thiophene has no hydrogen
    if aromatic {
        bonds += 1;
        valences = &valences[..valences.len().min(1)];
    }
    valences
        .iter()
        .find(|&&valence| valence >= bonds)
        .map_or(0, |valence| valence - bonds)
}

/// Adds the bonds, giving each one with no symbol its order: aromatic between two
/// aromatic atoms, single otherwise. Aromaticity perception later turns the aromatic
/// ones outside an aromatic ring, like biphenyl's link, into single bonds.
fn add_bonds(molecule: &mut Molecule, pending_bonds: Vec<PendingBond>) {
    for (node1, node2, order) in pending_bonds {
        let aromatic = |node| molecule.get_atom(node).is_some_and(|atom| atom.aromatic);
        let implicit = if aromatic(node1) && aromatic(node2) {
            BondOrder::Aromatic
        } else {
            BondOrder::Single
        };
        let order = order.unwrap_or(implicit);
        molecule.add_bond(node1, node2, Bond { order });
    }
}

/// Parses a molecule from the start of `input`; returns the rest (empty, or whitespace and
/// a title) and the molecule. Anything else left over is SMILES the reader does not
/// support, and an error, so it never yields a partial molecule.
pub fn parse_molecule(input: &str) -> ModalResult<(&str, Molecule)> {
    let (rest, molecule) = parse_molecule_next.parse_peek(input)?;
    if !rest.is_empty() && !rest.starts_with(char::is_whitespace) {
        return Err(ErrMode::Cut(ContextError::new()));
    }
    Ok((rest, molecule))
}

fn parse_molecule_next(input: &mut &str) -> ModalResult<Molecule> {
    let mut molecule = Molecule::new();
    let mut open_cycles: HashMap<u8, (NodeIndex, Option<BondOrder>)> = HashMap::new();
    let mut stack: Vec<NodeIndex> = Vec::new();
    let mut pending_bonds: Vec<PendingBond> = Vec::new();
    let mut organic_nodes: Vec<(NodeIndex, u8, bool)> = Vec::new();

    let tokens: Vec<Token> = repeat(
        0..,
        alt((
            parse_atom
                .with_taken()
                .map(|(atom, taken)| Token::Atom(atom, taken.starts_with('['))),
            parse_bond.map(|bond| Token::Bond(bond.order)),
            parse_cycle_digit.map(Token::Ring),
            '('.map(|_| Token::Open),
            ')'.map(|_| Token::Close),
            '.'.map(|_| Token::Dot),
        )),
    )
    .parse_next(input)?;

    let mut prev_node = NodeIndex::end();
    let mut prev_bond = None;
    let mut after_dot = false;

    for token in tokens {
        // a '.' must stand between two atoms: "CC.", ".CC" and "C..O" are truncated or
        // broken records, not molecules
        if after_dot && !matches!(token, Token::Atom(..)) {
            return Err(ErrMode::Cut(ContextError::new()));
        }
        // a bond symbol must be followed by an atom or a ring digit: "CC=", "C(C=)C" and
        // "C=.C" would otherwise drop it and give a different molecule
        if prev_bond.is_some() && !matches!(token, Token::Atom(..) | Token::Ring(_)) {
            return Err(ErrMode::Cut(ContextError::new()));
        }
        match token {
            Token::Open => stack.push(prev_node),
            Token::Close => {
                let Some(node) = stack.pop() else {
                    return Err(ErrMode::Cut(ContextError::new()));
                };
                prev_node = node;
            }
            Token::Ring(digit) => {
                // a ring bond needs an atom to start from or end at
                if prev_node == NodeIndex::end() {
                    return Err(ErrMode::Cut(ContextError::new()));
                }
                if let Some((open_node, open_bond)) = open_cycles.remove(&digit) {
                    // "C11" would bond an atom to itself, "C1C1" bond two atoms twice
                    let bonded = |(a, b, _): &PendingBond| {
                        (*a, *b) == (prev_node, open_node) || (*a, *b) == (open_node, prev_node)
                    };
                    if open_node == prev_node || pending_bonds.iter().any(bonded) {
                        return Err(ErrMode::Cut(ContextError::new()));
                    }
                    // the bond symbol may stand at either end: "C=1CCCCC1" or "C1CCCCC=1",
                    // but not two different ones: "C=1CCCCC-1"
                    if open_bond.is_some() && prev_bond.is_some() && open_bond != prev_bond {
                        return Err(ErrMode::Cut(ContextError::new()));
                    }
                    pending_bonds.push((prev_node, open_node, prev_bond.or(open_bond)));
                } else {
                    open_cycles.insert(digit, (prev_node, prev_bond));
                }
                prev_bond = None;
            }
            Token::Atom(atom, bracket) => {
                let atomic_number = atom.element.atomic_number;
                let aromatic = atom.aromatic;
                let node = molecule.add_atom(atom);
                if !bracket {
                    organic_nodes.push((node, atomic_number, aromatic));
                }
                if prev_node != NodeIndex::end() {
                    pending_bonds.push((prev_node, node, prev_bond));
                }
                prev_node = node;
                prev_bond = None;
                after_dot = false;
            }
            // "=C", "C.=C" and ".C" have no atom before the bond or dot
            Token::Bond(_) | Token::Dot if prev_node == NodeIndex::end() => {
                return Err(ErrMode::Cut(ContextError::new()));
            }
            Token::Bond(order) => prev_bond = Some(order),
            // the next atom starts a new, unbonded component
            Token::Dot => {
                prev_node = NodeIndex::end();
                after_dot = true;
            }
        }
    }

    if !open_cycles.is_empty() || !stack.is_empty() || prev_bond.is_some() || after_dot {
        return Err(ErrMode::Cut(ContextError::new()));
    }

    add_bonds(&mut molecule, pending_bonds);

    for (node, atomic_number, aromatic) in organic_nodes {
        let hs = implicit_hydrogens(&molecule, node, atomic_number, aromatic);
        if let Some(atom) = molecule.graph.node_weight_mut(node) {
            atom.implicit_hs = hs;
        }
    }
    // after the hydrogens: perception moves electrons around, not hydrogens
    molecule.perceive_aromaticity();

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
        assert_eq!(
            m.get_bond_by_atoms(NodeIndex::new(0), NodeIndex::new(1))
                .unwrap()
                .order,
            BondOrder::Single
        );
        assert!(m.has_bond(NodeIndex::new(1), NodeIndex::new(2)));
        assert_eq!(
            m.get_bond_by_atoms(NodeIndex::new(1), NodeIndex::new(2))
                .unwrap()
                .order,
            BondOrder::Single
        );
        assert!(m.has_bond(NodeIndex::new(0), NodeIndex::new(2)));
        assert_eq!(
            m.get_bond_by_atoms(NodeIndex::new(0), NodeIndex::new(2))
                .unwrap()
                .order,
            BondOrder::Double
        );
        assert!(m.has_bond(NodeIndex::new(2), NodeIndex::new(3)));
        assert_eq!(
            m.get_bond_by_atoms(NodeIndex::new(2), NodeIndex::new(3))
                .unwrap()
                .order,
            BondOrder::Single
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
        assert_eq!(
            m.get_bond_by_atoms(NodeIndex::new(0), NodeIndex::new(1))
                .unwrap()
                .order,
            BondOrder::Single
        );
        assert!(m.has_bond(NodeIndex::new(1), NodeIndex::new(2)));
        assert_eq!(
            m.get_bond_by_atoms(NodeIndex::new(1), NodeIndex::new(2))
                .unwrap()
                .order,
            BondOrder::Double
        );
        assert!(m.has_bond(NodeIndex::new(1), NodeIndex::new(3)));
        assert_eq!(
            m.get_bond_by_atoms(NodeIndex::new(1), NodeIndex::new(3))
                .unwrap()
                .order,
            BondOrder::Single
        );
        assert!(m.has_bond(NodeIndex::new(3), NodeIndex::new(4)));
        assert_eq!(
            m.get_bond_by_atoms(NodeIndex::new(3), NodeIndex::new(4))
                .unwrap()
                .order,
            BondOrder::Single
        );
        assert!(m.has_bond(NodeIndex::new(3), NodeIndex::new(0)));
        assert_eq!(
            m.get_bond_by_atoms(NodeIndex::new(3), NodeIndex::new(0))
                .unwrap()
                .order,
            BondOrder::Double
        );
    }

    #[test]
    fn parse_molecule_c1cc() {
        assert!(parse_molecule("C1CC").is_err());
    }

    #[test]
    fn parse_molecule_ring_digit_without_atom() {
        // used to panic in add_edge on the placeholder node
        assert!(parse_molecule("1C1").is_err());
        assert!(parse_molecule("(C)1C1").is_err());
    }

    #[test]
    fn parse_molecule_rejects_unparsed_rest() {
        // still unsupported, so an error instead of a partial molecule
        for smiles in ["C%1CC%1", "CC$C", "C[C@H](O)N", "Cq", "C[Xx]", "C[c"] {
            assert!(parse_molecule(smiles).is_err(), "{smiles}");
        }
        // a title after whitespace is fine, and comes back as the rest
        let (rest, m) = parse_molecule("CCO ethanol").unwrap();
        assert_eq!((rest, m.count_atoms()), (" ethanol", 3));
        let (rest, m) = parse_molecule("CCO\tethanol").unwrap();
        assert_eq!((rest, m.count_atoms()), ("\tethanol", 3));
    }

    #[test]
    fn parse_aromatic_molecules() {
        let hs = |smiles: &str| -> Vec<u8> {
            let m = parse_molecule(smiles).unwrap().1;
            m.graph
                .node_weights()
                .map(|atom| atom.hs + atom.implicit_hs)
                .collect()
        };
        assert_eq!(hs("c1ccccc1"), [1; 6]);
        // toluene: the substituted carbon has none
        assert_eq!(hs("Cc1ccccc1"), [3, 0, 1, 1, 1, 1, 1]);
        // pyridine, pyrrole, furan and thiophene: heteroatoms with no hydrogen unless written
        assert_eq!(hs("n1ccccc1"), [0, 1, 1, 1, 1, 1]);
        assert_eq!(hs("[nH]1cccc1"), [1, 1, 1, 1, 1]);
        assert_eq!(hs("o1cccc1"), [0, 1, 1, 1, 1]);
        assert_eq!(hs("s1cccc1"), [0, 1, 1, 1, 1]);
        // fused rings: naphthalene's shared carbons have none
        assert_eq!(hs("c1ccc2ccccc2c1"), [1, 1, 1, 0, 1, 1, 1, 1, 0, 1]);
        // aromatic and Kekulé benzene weigh the same
        for smiles in ["c1ccccc1", "C1=CC=CC=C1"] {
            let m = parse_molecule(smiles).unwrap().1;
            assert!((m.weight() - 78.114).abs() < 1e-9, "{smiles}");
        }
    }

    #[test]
    fn parse_aromatic_bond_orders() {
        let orders = |smiles: &str| -> Vec<BondOrder> {
            let m = parse_molecule(smiles).unwrap().1;
            m.graph.edge_weights().map(|bond| bond.order).collect()
        };
        // unwritten bonds between aromatic atoms are aromatic, the ring closure too
        assert_eq!(orders("c1ccccc1"), [BondOrder::Aromatic; 6]);
        // between an aromatic and an aliphatic atom they are single
        assert_eq!(orders("Cc1ccccc1")[0], BondOrder::Single);
        // the bond joining biphenyl's rings is in no ring, so single
        let biphenyl = orders("c1ccccc1c1ccccc1");
        assert_eq!(
            biphenyl.iter().filter(|&&o| o == BondOrder::Single).count(),
            1
        );
        // a written bond wins (edge 5 is the first ring's closure, 6 the "-")
        let written = orders("c1ccccc1-c1ccccc1");
        assert_eq!(
            (written[5], written[6]),
            (BondOrder::Aromatic, BondOrder::Single)
        );
        assert_eq!(orders("c1ccccc1=O")[6], BondOrder::Double);
    }

    #[test]
    fn parse_dots_wildcards_and_any_element() {
        // a salt: two components, no bond between them
        let m = parse_molecule("[Na+].[Cl-]").unwrap().1;
        assert_eq!((m.count_atoms(), m.count_bonds()), (2, 0));
        assert!((m.weight() - 58.44).abs() < 1e-9);
        let m = parse_molecule("CCO.O").unwrap().1;
        assert_eq!((m.count_atoms(), m.count_bonds()), (4, 2));
        // a ring bond may cross a dot
        let m = parse_molecule("C1.C1").unwrap().1;
        assert_eq!(m.count_bonds(), 1);
        // the wildcard gets no hydrogens and no weight
        let m = parse_molecule("*C").unwrap().1;
        assert!((m.weight() - 15.035).abs() < 1e-9);
        assert!(parse_molecule("[*]CC").is_ok());
        // every element in brackets, isotopes past 255, boron outside brackets
        let m = parse_molecule("[238U]").unwrap().1;
        assert_eq!(m.get_atom(NodeIndex::new(0)).unwrap().isotope, 238);
        assert!(parse_molecule("[Og]").is_ok());
        let m = parse_molecule("B(C)(C)C").unwrap().1;
        assert_eq!(m.get_atom(NodeIndex::new(0)).unwrap().implicit_hs, 0);
        // a dot needs an atom on both sides
        for smiles in [
            "CC.", ".CC", "C..O", ".", "[Na+].", "C(C.)C", "C.(C)", ".1C1",
        ] {
            assert!(parse_molecule(smiles).is_err(), "{smiles}");
        }
        assert!(parse_molecule("C(.C)C").is_ok());
    }

    #[test]
    fn parse_rejects_misplaced_bond_symbols() {
        // "CC=" used to parse as ethane: a bond symbol with no atom after it was dropped
        for smiles in [
            "CC=",
            "C(C=)C",
            "C=.C",
            "C=(O)C",
            "=C",
            "C.=C",
            "C==C",
            "C=1CCCCC-1",
        ] {
            assert!(parse_molecule(smiles).is_err(), "{smiles}");
        }
        // bonds where they belong still parse, the same symbol at both ring ends too
        for smiles in [
            "C=C",
            "C(=O)O",
            "C=1CCCCC1",
            "C1CCCCC=1",
            "C=1CCCCC=1",
            "C1.C=1",
            "C=C.O",
        ] {
            assert!(parse_molecule(smiles).is_ok(), "{smiles}");
        }
    }

    #[test]
    fn parse_bond_at_ring_opening() {
        // the "=" belongs to the ring bond, not to the next chain bond
        let m = parse_molecule("C=1CCCCC1").unwrap().1;
        assert_eq!(
            m.get_bond_by_atoms(NodeIndex::new(0), NodeIndex::new(5))
                .unwrap()
                .order,
            BondOrder::Double
        );
        assert_eq!(
            m.get_bond_by_atoms(NodeIndex::new(0), NodeIndex::new(1))
                .unwrap()
                .order,
            BondOrder::Single
        );
    }

    #[test]
    fn parse_molecule_ring_closures() {
        let counts = |smiles: &str| {
            let m = parse_molecule(smiles).unwrap().1;
            (m.count_atoms(), m.count_bonds())
        };
        // each digit is its own ring: bicyclo[2.2.0]hexane opens rings 1 and 2 on one atom
        assert_eq!(counts("C12CCC1CC2"), (6, 7));
        // two-digit ring numbers after %
        assert_eq!(counts("C%10CC%10"), (3, 3));
        assert_eq!(counts("C%12CCC%12"), (4, 4));
        // a ring number can be reused after it closes
        assert_eq!(counts("C1CC1C1CC1"), (6, 7));
        // % with one digit, an atom bonded to itself, and a second bond between two atoms
        for smiles in ["C%1CC1", "C11", "C1C1", "C12CC12"] {
            assert!(parse_molecule(smiles).is_err(), "{smiles}");
        }
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
            println!("{smiles}: ");
            let m = parse_molecule(smiles).unwrap().1;
            let result = m.ecfp(2, 512);
            println!("{result:?}");
        }
    }
}
