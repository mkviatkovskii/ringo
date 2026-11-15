pub mod common;

use ringo::io::smiles::reader::molecule::parse_molecule;
use ringo::model::stereo::{Configuration, Ligand};

#[test]
fn corpus_fixtures_preserve_every_marker() {
    for line in include_str!("data/stereo.smi").lines() {
        let (smiles, _) = line.split_once(' ').unwrap();
        let molecule = parse_molecule(line).unwrap_or_else(|_| panic!("{line}")).1;
        let markers = smiles
            .bytes()
            .enumerate()
            .filter(|&(i, b)| b == b'@' && (i == 0 || smiles.as_bytes()[i - 1] != b'@'))
            .count();
        assert_eq!(
            molecule.atom_stereo().len() + molecule.allene_stereo().len(),
            markers,
            "{line}"
        );
    }
}

#[test]
fn extended_classes_keep_permutation_and_missing_site() {
    assert_eq!(parity("N[S@](Cl)=O"), parity("[S@@](N)(Cl)=O"));
    let molecule = parse_molecule("O[Mn@OH28](Cl)(C)(N)F").unwrap().1;
    let s = &molecule.atom_stereo()[0];
    assert_eq!(s.configuration, Configuration::Octahedral(28));
    assert_eq!(s.ligands.len(), 6);
    assert_eq!(s.ligands[1], Ligand::Vacancy);
    let molecule = parse_molecule("S[As@TB5](F)(N)(Cl)Br").unwrap().1;
    assert_eq!(
        molecule.atom_stereo()[0].configuration,
        Configuration::TrigonalBipyramidal(5)
    );
    let molecule = parse_molecule("O=[S@@](O)C").unwrap().1;
    assert_eq!(molecule.atom_stereo()[0].ligands[1], Ligand::LonePair);
    assert_eq!(parity("[C@H](F)(Cl)Br"), parity("[C@]([H])(F)(Cl)Br"));
    let m = parse_molecule("[2H][C@](F)(Cl)Br").unwrap().1;
    assert_eq!(
        m.atom_stereo()[0].ligands[0],
        Ligand::Atom(petgraph::stable_graph::NodeIndex::new(0))
    );
    assert_eq!(
        m.get_atom(petgraph::stable_graph::NodeIndex::new(0))
            .unwrap()
            .isotope,
        2
    );
}

#[expect(
    clippy::unwrap_used,
    reason = "test helper requires a parsed centre with four valid ligands"
)]
fn parity(smiles: &str) -> u8 {
    let molecule = parse_molecule(smiles).unwrap().1;
    let stereo = &molecule.atom_stereo()[0];
    let Configuration::Tetrahedral(n) = stereo.configuration else {
        panic!()
    };
    let elements: Vec<_> = stereo
        .ligands
        .iter()
        .map(|l| match l {
            Ligand::Atom(a) => molecule.get_atom(*a).unwrap().element.atomic_number,
            Ligand::Hydrogen => 1,
            Ligand::LonePair => 0,
            Ligand::Vacancy => panic!(),
        })
        .collect();
    let inversions = (0..4)
        .flat_map(|i| (i + 1..4).map(move |j| (i, j)))
        .filter(|&(i, j)| elements[i] > elements[j])
        .count();
    (n - 1) ^ u8::try_from(inversions % 2).unwrap()
}

#[test]
fn tetrahedral_permutations_preserve_orientation() {
    let base = parity("N[C@](Br)(O)C");
    assert_eq!(base, parity("N[C@@](O)(Br)C"));
    assert_ne!(base, parity("N[C@@](Br)(O)C"));
    let names = ["F", "Cl", "Br", "I"];
    let mut permutations = 0;
    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                for d in 0..4 {
                    let p = [a, b, c, d];
                    if p.iter()
                        .copied()
                        .collect::<std::collections::BTreeSet<_>>()
                        .len()
                        != 4
                    {
                        continue;
                    }
                    let odd = (0..4)
                        .flat_map(|i| (i + 1..4).map(move |j| (i, j)))
                        .filter(|&(i, j)| p[i] > p[j])
                        .count()
                        % 2
                        == 1;
                    let tag = if odd { "@@" } else { "@" };
                    let smiles = format!(
                        "[C{tag}]({})({})({}){}",
                        names[a], names[b], names[c], names[d]
                    );
                    assert_eq!(parity(&smiles), 0, "{smiles}");
                    permutations += 1;
                }
            }
        }
    }
    assert_eq!(permutations, 24);
}

#[test]
fn hydrogens_and_ring_slots_keep_their_positions() {
    assert_eq!(parity("[C@H](F)(Cl)Br"), parity("[H][C@](F)(Cl)Br"));
    assert_eq!(parity("F[C@H](Cl)Br"), parity("F[C@]([H])(Cl)Br"));
    assert_eq!(parity("N[C@](Br)(O)C"), parity("N[C@TH1](Br)(O)C"));
    let a = parse_molecule("C[C@H]1CCCCO1").unwrap().1;
    let b = parse_molecule("O1CCCC[C@@H]1C").unwrap().1;
    assert_eq!(
        a.atom_stereo()[0].configuration,
        Configuration::Tetrahedral(1)
    );
    assert_eq!(
        b.atom_stereo()[0].configuration,
        Configuration::Tetrahedral(2)
    );
    let indices = |m: &ringo::model::molecule::Molecule| {
        m.atom_stereo()[0]
            .ligands
            .iter()
            .map(|l| match l {
                Ligand::Atom(n) => n.index(),
                Ligand::Hydrogen => usize::MAX,
                _ => panic!(),
            })
            .collect::<Vec<_>>()
    };
    assert_eq!(indices(&a), [0, usize::MAX, 6, 2]);
    assert_eq!(indices(&b), [4, usize::MAX, 0, 6]);
    let m = parse_molecule("F[C@]1(Cl)CCC(Br)C1").unwrap().1;
    let s = &m.atom_stereo()[0];
    let Ligand::Atom(ring) = s.ligands[1] else {
        panic!()
    };
    assert_eq!(ring.index(), 7);
    assert_eq!(s.ligands.len(), 4);
}

#[test]
fn directional_meaning_and_errors() {
    let m = parse_molecule("O=C(/C=C/c1ccccc1)/C=C/c1ccccc1").unwrap().1; // CHEMBL17201
    assert_eq!(m.double_bond_stereo().len(), 2);
    let m = parse_molecule("CC(=O)/N=c1/sc(S(N)(=O)=O)nn1C").unwrap().1; // CHEBI:6822
    assert_eq!(m.double_bond_stereo().len(), 1);
    let a = parse_molecule("F/C=C1CCCCC/1").unwrap().1;
    let b = parse_molecule("F/C=C\\1CCCCC1").unwrap().1;
    assert_eq!(
        a.double_bond_stereo()[0].same_side,
        b.double_bond_stereo()[0].same_side
    );
    assert!(parse_molecule("F/C=C/1CCCCC/1").is_err());
    for (smiles, same_side) in [
        ("F/C=C/F", false),
        ("F\\C=C\\F", false),
        ("C(\\F)=C/F", false),
        ("F/C=C\\F", true),
    ] {
        let m = parse_molecule(smiles).unwrap().1;
        assert_eq!(m.double_bond_stereo().len(), 1);
        assert_eq!(m.double_bond_stereo()[0].same_side, same_side, "{smiles}");
    }
    assert_eq!(
        parse_molecule("F/C=C/C=CC")
            .unwrap()
            .1
            .double_bond_stereo()
            .len(),
        1
    );
    assert_eq!(parse_molecule("FC=CF").unwrap().1.double_bond_stereo(), []);
    for smiles in [
        "C/C",
        "F/C=CF",
        "C/C(\\F)=C/F",
        "F//C=C/F",
        "[C@H2](F)Cl",
        "[C@TH3](F)(Cl)(Br)I",
        "[Co@OH31](F)(Cl)(Br)(I)(N)O",
        "F[C@](Cl)Br",
    ] {
        assert!(parse_molecule(smiles).is_err(), "{smiles}");
    }
}

#[test]
fn macrocyclic_lactone_keeps_its_configured_double_bond() {
    use ringo::model::bond::BondOrder;
    for (smiles, plain, aromatic_atoms) in [
        // CHEMBL2004566, ChEMBL 37 line 2,138,582: the ester O bridges a 10-atom ring
        // with a naive Hückel count of 10, but the three benzene rings alone are aromatic
        (
            "O=C1Nc2ccccc2C(=O)O/C(Cc2ccccc2)=N\\c2ccccc21",
            "O=C1Nc2ccccc2C(=O)OC(Cc2ccccc2)=Nc2ccccc21",
            18,
        ),
        // the same 10-atom ring, reduced to its conjugated atoms
        ("O=C1NC=CC(=O)O/C(C)=N\\C=C1", "O=C1NC=CC(=O)OC(C)=NC=C1", 0),
    ] {
        let m = parse_molecule(smiles)
            .unwrap_or_else(|_| panic!("{smiles}"))
            .1;
        let [stereo] = m.double_bond_stereo() else {
            panic!("{smiles}")
        };
        assert_eq!(m.get_bond(stereo.bond).unwrap().order, BondOrder::Double);
        let elements = stereo
            .ends
            .map(|end| m.get_atom(end).unwrap().element.atomic_number);
        assert_eq!(elements, [6, 7], "{smiles}");
        let aromatic = m.graph().node_weights().filter(|a| a.aromatic).count();
        assert_eq!(aromatic, aromatic_atoms, "{smiles}");
        m.validate().unwrap();
        // the ring is not aromatic because of its constitution, not its stereo marks
        let p = parse_molecule(plain).unwrap().1;
        assert_eq!(p.exact_key(), m.exact_key(), "{plain}");
        assert_eq!(p.ecfp(2, 512).0, m.ecfp(2, 512).0, "{plain}");
    }
}

#[test]
fn indexes_newly_supported_stereo_without_changing_fingerprints() {
    common::assert_search(
        "stereo",
        "N[C@H](C)C(=O)O alanine\nF/C=C/F alkene\n",
        ["N[C@@H](C)C(=O)O", "F/C=C\\F"],
    );
}

#[test]
#[ignore = "requires the downloaded ChEBI/ChEMBL files; set RINGO_SMILES_CORPUS"]
fn corpus_coverage() {
    use std::io::BufRead;
    let root = std::env::var("RINGO_SMILES_CORPUS").unwrap();
    for name in ["chebi_20261105.smi", "chembl_37.smi"] {
        let corpus_path = std::path::Path::new(&root).join(name);
        let before = std::fs::metadata(&corpus_path).unwrap();
        let mut counts = [0usize; 4];
        let mut examples = Vec::new();
        // by kind, so valence and aromaticity changes can be reported apart
        let mut validation_failed: std::collections::BTreeMap<&str, (usize, Vec<String>)> =
            std::collections::BTreeMap::new();
        for line in std::io::BufReader::new(std::fs::File::open(&corpus_path).unwrap()).lines() {
            let line = line.unwrap();
            let smiles = line.split_whitespace().next().unwrap_or("");
            let stereo = smiles.contains(['@', '/', '\\']);
            let success = match parse_molecule(&line) {
                Ok((_, molecule)) => {
                    if let Err(error) = molecule.validate() {
                        let kind = [
                            "no Kekulé structure",
                            "aromatic atom outside a ring",
                            "wildcard atom in an aromatic ring",
                        ]
                        .into_iter()
                        .find(|kind| error.contains(kind))
                        .unwrap_or("valence or charge");
                        let (count, kind_examples) = validation_failed.entry(kind).or_default();
                        *count += 1;
                        if kind_examples.len() < 5 {
                            kind_examples.push(format!("{line}: {error}"));
                        }
                    }
                    true
                }
                Err(_) => false,
            };
            counts[usize::from(stereo) * 2 + usize::from(success)] += 1;
            if stereo && !success && examples.len() < 10 {
                examples.push(line);
            }
        }
        let after = std::fs::metadata(&corpus_path).unwrap();
        assert_eq!(
            (before.len(), before.modified().unwrap()),
            (after.len(), after.modified().unwrap()),
            "corpus changed during scan"
        );
        println!("{name}: [plain rejected, plain accepted, stereo rejected, stereo accepted] = {counts:?}; examples={examples:?}");
        for (kind, (count, kind_examples)) in &validation_failed {
            println!("{name}: loaded then failed validation ({kind}) = {count}; examples={kind_examples:?}");
        }
    }
}
