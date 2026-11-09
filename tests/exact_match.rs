use ringo::io::smiles::reader::molecule::parse_molecule;
use ringo::model::exact_match::ExactMatchError;
use ringo::model::molecule::Molecule;

fn parse(smiles: &str) -> Molecule {
    parse_molecule(smiles)
        .unwrap_or_else(|_| panic!("{smiles}"))
        .1
}

#[expect(
    clippy::unwrap_used,
    reason = "test helper requires stereo the matcher supports"
)]
fn same(a: &str, b: &str) -> bool {
    let (a_molecule, b_molecule) = (parse(a), parse(b));
    let forward = a_molecule.exact_match(&b_molecule).unwrap();
    assert_eq!(
        forward,
        b_molecule.exact_match(&a_molecule).unwrap(),
        "{a} and {b} match one way only"
    );
    if forward {
        assert_eq!(
            a_molecule.exact_key(),
            b_molecule.exact_key(),
            "{a} and {b} match but have different keys"
        );
    }
    forward
}

#[test]
fn keys_ignore_stereo_and_separate_structures() {
    let key = |smiles| parse(smiles).exact_key();
    for (a, b) in [
        ("C[C@@H](C(=O)O)N", "C[C@H](C(=O)O)N"),
        ("C[C@@H](C(=O)O)N", "CC(N)C(=O)O"),
        ("F/C=C/F", "F/C=C\\F"),
        ("NC(Br)=[C@]=C(O)C", "NC(Br)=[C@@]=C(O)C"),
    ] {
        assert_eq!(key(a), key(b), "{a} {b}");
    }
    for (a, b) in [
        ("CCO", "COC"),
        ("C=C", "CC"),
        ("CC", "C.C"),
        ("N->[Cu]", "N<-[Cu]"),
        ("[2H]C", "C"),
        // alike at every atom's first shell
        ("C1CCCCC1", "C1CCC1.C1CCC1"),
        ("CC(C)CC", "CCCCC"),
    ] {
        assert_ne!(key(a), key(b), "{a} {b}");
    }
    // pins the hash and the refinement: if this changes, existing .key files find nothing
    assert_eq!(key("CC(=O)O"), 0xedda_6aea_80c8_7eac);
}

fn assert_same(pairs: &[(&str, &str)]) {
    for (a, b) in pairs {
        assert!(same(a, b), "{a} should match {b}");
    }
}

fn assert_different(pairs: &[(&str, &str)]) {
    for (a, b) in pairs {
        assert!(!same(a, b), "{a} should not match {b}");
    }
}

#[test]
fn spellings_of_one_molecule_match() {
    assert_same(&[
        ("CCO", "OCC"),
        ("CC(=O)O", "OC(C)=O"),
        ("C1CCCCC1", "C2CCCCC2"),
        ("C1CCCCC1", "C%10CCCCC%10"),
        ("c1ccccc1", "C1=CC=CC=C1"),
        ("Cc1ccccc1", "C1=CC=C(C)C=C1"),
        ("c1ccncc1", "C1=CC=NC=C1"),
        ("c1ccc2ccccc2c1", "C1=CC=C2C=CC=CC2=C1"),
        ("C", "[CH4]"),
        ("C", "[H]C([H])([H])[H]"),
        ("CO", "[H]OC"),
        // atom classes are labels
        ("[CH3:1]O", "CO"),
        ("[Na+].[Cl-]", "[Cl-].[Na+]"),
        ("N->[Cu]", "[Cu]<-N"),
    ]);
}

#[test]
fn different_molecules_dont_match() {
    assert_different(&[
        ("CCO", "COC"),
        ("c1ccccc1", "C1CCCCC1"),
        ("C=C", "CC"),
        ("N", "[NH4+]"),
        ("C", "[13CH4]"),
        // deuterium stays an atom
        ("[2H]C", "C"),
        ("CCO", "CCO.O"),
        ("CC", "C.C"),
        ("N->[Cu]", "N<-[Cu]"),
        // only the direction differs: bracket atoms keep their written hydrogens
        ("[NH3]->[Cu]", "[NH3]<-[Cu]"),
    ]);
}

#[test]
fn tetrahedral_centres_compare_by_parity() {
    // L-alanine three ways, D-alanine
    assert_same(&[
        ("C[C@@H](C(=O)O)N", "N[C@@H](C)C(=O)O"),
        ("C[C@@H](C(=O)O)N", "C[C@H](N)C(=O)O"),
        ("[C@H](F)(Cl)Br", "[C@]([H])(F)(Cl)Br"),
        ("F[C@H](Cl)Br", "F[C@]([H])(Cl)Br"),
        // the ring-closure digit holds its neighbour's place
        ("C[C@H]1CCCCO1", "C[C@@H]1OCCCC1"),
    ]);
    assert_different(&[
        ("C[C@@H](C(=O)O)N", "C[C@H](C(=O)O)N"),
        ("F[C@H](Cl)Br", "F[C@@]([H])(Cl)Br"),
        ("C[C@H]1CCCCO1", "C[C@H]1OCCCC1"),
        // a real centre against an unspecified one
        ("CC[C@H](C)O", "CCC(C)O"),
    ]);
}

#[test]
fn stereo_that_symmetry_undoes_is_ignored() {
    assert_same(&[
        // two methyls, or two hydrogens, on the centre
        ("C[C@H](C)O", "CC(C)O"),
        ("C[C@H](C)O", "C[C@@H](C)O"),
        ("[H][C@]([H])(F)Cl", "FCCl"),
        // two methyls, or two hydrogens, on one end of the double bond
        ("C/C=C(/C)C", "CC=C(C)C"),
        ("[H]/C([H])=C/F", "C=CF"),
        // (2R,4R)-pentane-2,3,4-triol: the middle carbon's branches are alike
        (
            "C[C@@H](O)[C@H](O)[C@@H](C)O",
            "C[C@@H](O)[C@@H](O)[C@@H](C)O",
        ),
        ("C[C@@H](O)[C@H](O)[C@@H](C)O", "C[C@@H](O)C(O)[C@@H](C)O"),
        // the two branches are alike, dative direction included
        (
            "F[C@H]([NH2]->[Cu])[NH2]->[Cu]",
            "FC([NH2]->[Cu])[NH2]->[Cu]",
        ),
    ]);
    // dative direction tells the branches apart
    assert_different(&[
        (
            "F[C@H]([NH2]->[Cu])[NH2]<-[Cu]",
            "F[C@@H]([NH2]->[Cu])[NH2]<-[Cu]",
        ),
        (
            "F[C@H]([NH2]->[Cu])[NH2]<-[Cu]",
            "FC([NH2]->[Cu])[NH2]<-[Cu]",
        ),
    ]);
    // the meso forms: the middle carbon is pseudo-asymmetric, its branches differ only by
    // their own stereo
    assert_different(&[
        (
            "C[C@@H](O)[C@H](O)[C@H](C)O",
            "C[C@@H](O)[C@@H](O)[C@H](C)O",
        ),
        ("C[C@@H](O)[C@H](O)[C@H](C)O", "C[C@@H](O)C(O)[C@H](C)O"),
        ("C[C@@H](O)C(O)[C@@H](C)O", "C[C@@H](O)C(O)[C@H](C)O"),
    ]);
}

#[test]
fn double_bonds_compare_by_side() {
    assert_same(&[
        ("F/C=C/F", "F\\C=C\\F"),
        ("F/C=C/F", "C(\\F)=C/F"),
        // a hydrogen reference stands opposite the end's other substituent
        ("[H]/C(F)=C/F", "F/C=C\\F"),
    ]);
    assert_different(&[("F/C=C/F", "F/C=C\\F"), ("F/C=C/F", "FC=CF")]);
}

#[test]
fn allenes_compare_as_collapsed_tetrahedra() {
    assert_same(&[
        ("NC(Br)=[C@]=C(O)C", "NC(Br)=[C@AL1]=C(O)C"),
        ("NC(Br)=[C@]=C(O)C", "BrC(N)=[C@@]=C(O)C"),
        ("NC(Br)=[C@]=C(O)C", "CC(O)=[C@]=C(Br)N"),
        ("OC=[C@]=CF", "OC([H])=[C@AL1]=C([H])F"),
        ("OC=[C@]=CF", "FC=[C@]=CO"),
        // two methyls, or two hydrogens, on a terminal atom
        ("CC(C)=[C@]=CF", "CC(C)=C=CF"),
        ("C=[C@]=C", "C=C=C"),
    ]);
    assert_different(&[
        ("NC(Br)=[C@]=C(O)C", "NC(Br)=[C@AL2]=C(O)C"),
        ("NC(Br)=[C@]=C(O)C", "NC(Br)=C=C(O)C"),
        // one hydrogen on each terminal is real stereo
        ("OC=[C@]=CF", "OC=[C@@]=CF"),
        ("OC=[C@]=CF", "OC=C=CF"),
    ]);
}

#[test]
fn alike_end_atoms_fold_into_their_neighbour() {
    assert_same(&[
        ("CN(C)C", "N(C)(C)C"),
        ("CS(=O)(=O)C", "O=S(C)(=O)C"),
        ("FC(F)(F)C", "CC(F)(F)F"),
        // two methyls on the centre
        ("C[C@](C)(F)O", "CC(C)(F)O"),
        // L-valine: the isopropyl's methyls don't touch the centre's stereo
        ("CC(C)[C@H](N)C(=O)O", "N[C@@H](C(C)C)C(=O)O"),
    ]);
    assert_different(&[
        ("CN(C)C", "CNCC"),
        ("CS(=O)(=O)C", "CS(=O)C"),
        ("FC(F)(F)C", "FC(F)(Cl)C"),
        ("CC(C)[C@H](N)C(=O)O", "CC(C)[C@@H](N)C(=O)O"),
        ("CC(C)(C)[C@H](O)F", "CC(C)(C)[C@@H](O)F"),
    ]);
}

#[test]
fn components_match_in_any_order() {
    let salt = format!("C[C@H](N)C(=O)O{}", ".Cl".repeat(14));
    let reordered = format!("{}N[C@@H](C)C(=O)O", "Cl.".repeat(14));
    assert_same(&[("O.O.O.CCO", "CCO.O.O.O"), (&salt, &reordered)]);
    assert_different(&[
        ("CCO.O", "CCO.C"),
        ("[Na+].[Cl-].[Na+]", "[Na+].[Cl-].[Cl-]"),
    ]);
}

#[test]
fn other_stereo_is_reported_unsupported() {
    let a = parse("C[Pt@SP1](F)(Cl)[H]");
    assert_eq!(a.exact_match(&a), Err(ExactMatchError::UnsupportedStereo));
}

#[test]
fn corpus_molecules_match_themselves_only() {
    let corpus: Vec<&str> = [
        include_str!("data/molecules.smi"),
        include_str!("data/stereo.smi"),
    ]
    .iter()
    .flat_map(|file| file.lines())
    .map(|line| line.split_whitespace().next().unwrap_or(line))
    .collect();
    for (i, a) in corpus.iter().enumerate() {
        let molecule = parse(a);
        match molecule.exact_match(&molecule) {
            Ok(matched) => assert!(matched, "{a}"),
            Err(error) => assert_eq!(error, ExactMatchError::UnsupportedStereo, "{a}"),
        }
        for b in &corpus[i + 1..] {
            if let Ok(matched) = molecule.exact_match(&parse(b)) {
                assert!(!matched, "{a} matches {b}");
            }
        }
    }
}
