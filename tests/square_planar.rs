pub mod common;

use ringo::io::smiles::reader::molecule::parse_molecule;
use ringo::model::stereo::{Configuration, Ligand};

#[expect(
    clippy::unwrap_used,
    reason = "test helper requires a parsed square-planar centre and its four ligands"
)]
fn opposite_elements(smiles: &str) -> [[u8; 2]; 2] {
    let molecule = parse_molecule(smiles)
        .unwrap_or_else(|_| panic!("{smiles}"))
        .1;
    let stereo = &molecule.atom_stereo()[0];
    // OpenSMILES U, 4 and Z shapes, in SMILES encounter order.
    let indices = match stereo.configuration {
        Configuration::SquarePlanar(1) => [[0, 2], [1, 3]],
        Configuration::SquarePlanar(2) => [[0, 1], [2, 3]],
        Configuration::SquarePlanar(3) => [[0, 3], [1, 2]],
        _ => panic!("{smiles}"),
    };
    let mut pairs = indices.map(|pair| {
        let mut elements = pair.map(|i| match stereo.ligands[i] {
            Ligand::Atom(node) => molecule.get_atom(node).unwrap().element.atomic_number,
            Ligand::Hydrogen => 1,
            Ligand::Vacancy => 0,
            Ligand::LonePair => panic!("{smiles}"),
        });
        elements.sort_unstable();
        elements
    });
    pairs.sort_unstable();
    pairs
}

#[test]
fn permutations_preserve_opposite_ligands() {
    let expected = [[1, 9], [6, 17]];
    for smiles in [
        "C[Pt@SP1](F)(Cl)[H]",
        "C[Pt@SP2](Cl)(F)[H]",
        "C[Pt@SP3](F)([H])Cl",
        "F[Pt@SP1](Cl)([H])C",
        "[H][Pt@SP1](Cl)(F)C",
    ] {
        assert_eq!(opposite_elements(smiles), expected, "{smiles}");
    }
    for smiles in ["C[Pt@SP2](F)(Cl)[H]", "C[Pt@SP3](F)(Cl)[H]"] {
        assert_ne!(opposite_elements(smiles), expected, "{smiles}");
    }
    for smiles in ["Cl[Pt@SP1]1(F)CCN->1", "Cl[Pt@SP2](F)1CCN->1"] {
        assert_eq!(opposite_elements(smiles), [[6, 7], [9, 17]], "{smiles}");
    }
}

#[test]
fn hydrogens_and_vacancies_keep_their_sites() {
    for (smiles, position, ligand) in [
        ("[Pt@SP1](C)(Cl)F", 0, Ligand::Vacancy),
        ("C[Pt@SP1](Cl)F", 1, Ligand::Vacancy),
        ("[Pt@SP1H](C)(Cl)F", 0, Ligand::Hydrogen),
        ("C[Pt@SP1H](Cl)F", 1, Ligand::Hydrogen),
        ("[P@SP1](F)(Cl)Br", 0, Ligand::Vacancy),
    ] {
        let molecule = parse_molecule(smiles).unwrap().1;
        let stereo = &molecule.atom_stereo()[0];
        assert_eq!(stereo.configuration, Configuration::SquarePlanar(1));
        assert_eq!(stereo.ligands.len(), 4, "{smiles}");
        assert_eq!(stereo.ligands[position], ligand, "{smiles}");
    }
    for (a, b) in [
        ("C[Pt@SP1H](Cl)F", "C[Pt@SP1]([H])(Cl)F"),
        ("[Pt@SP1H](C)(Cl)F", "[H][Pt@SP1](C)(Cl)F"),
        ("C[Pt@SP1H2]Cl", "C[Pt@SP1]([H])([H])Cl"),
        ("C[Pt@SP1](Cl)F", "C[Pt@SP1](*)(Cl)F"),
        ("[Pt@SP1](C)(Cl)F", "*[Pt@SP1](C)(Cl)F"),
    ] {
        assert_eq!(opposite_elements(a), opposite_elements(b), "{a} {b}");
    }
}

#[test]
fn cis_and_trans_platinum_keep_different_arrangements() {
    // RDKit Book examples with repeated ligands and dative bonds.
    let cis = "Cl[Pt@SP1](Cl)(<-[NH3])<-[NH3]";
    let trans = "Cl[Pt@SP2](Cl)(<-[NH3])<-[NH3]";
    assert_eq!(opposite_elements(cis), [[7, 17], [7, 17]]);
    assert_eq!(opposite_elements(trans), [[7, 7], [17, 17]]);
    assert_eq!(
        opposite_elements(cis),
        opposite_elements("[NH3]->[Pt@SP1](Cl)(Cl)<-[NH3]")
    );
    assert_eq!(
        opposite_elements(trans),
        opposite_elements("[NH3]->[Pt@SP3](Cl)(Cl)<-[NH3]")
    );
    common::assert_search("square_planar", &format!("{cis} cisplatin\n"), [trans]);
}

#[test]
fn malformed_tags_and_excess_sites_are_rejected() {
    for smiles in [
        "[Pt@SP](F)(Cl)(Br)I",
        "[Pt@SP0](F)(Cl)(Br)I",
        "[Pt@SP4](F)(Cl)(Br)I",
        "[Pt@SP256](F)(Cl)(Br)I",
        "[Pt@SP-1](F)(Cl)(Br)I",
        "[Pt@SP1foo](F)(Cl)(Br)I",
        "[Pt@SP1](F)(Cl)(Br)(I)N",
        "F[Pt@SP1H2](Cl)Br",
        "[Pt@SP1H5]",
    ] {
        assert!(parse_molecule(smiles).is_err(), "{smiles}");
    }
    for (smiles, permutation) in [("F[Pt@](Cl)(Br)I", 1), ("F[Pt@@](Cl)(Br)I", 2)] {
        let molecule = parse_molecule(smiles).unwrap().1;
        assert_eq!(
            molecule.atom_stereo()[0].configuration,
            Configuration::Tetrahedral(permutation)
        );
    }
}
