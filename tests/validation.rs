use ringo::io::smiles::reader::molecule::parse_molecule;

#[test]
fn valence_is_optional_and_accounts_for_charge_and_hydrogens() {
    for smiles in [
        "C(C)(C)(C)(C)C",
        "[CH5]",
        "N(C)(C)(C)C",
        "[NH5+]",
        "O(C)(C)C",
        "[OH2-]",
        "[BH5-]",
        "F(F)F",
        "P(F)(F)(F)(F)(F)F",
        "S(F)(F)(F)(F)(F)(F)F",
        "[H]([H])[H]",
        "[H-]C",
        "[C-](C)(C)(C)C",
        "[N+127]",
        "[N-46]",
        "[Cl+3](F)(F)(F)(F)F",
        "CN1=CC=CC=C=1", // aromaticity must not hide the written nitrogen valence of five
        "Br[Br-]Br",
        "I[I-]I",
        "[I-](F)(F)(F)(F)(F)F",
    ] {
        let molecule = parse_molecule(smiles).unwrap().1;
        assert!(molecule.validate().is_err(), "{smiles}");
    }
    for smiles in [
        "CCO",
        "[CH3]",
        "[NH4+]",
        "[OH3+]",
        "[OH-]",
        "[BH4-]",
        "[CH3-]",
        "[CH3+]",
        "[H+]",
        "[H-]",
        "[2H]O[2H]",
        "[SiH4]",
        "P(F)(F)(F)(F)F",
        "S(F)(F)(F)(F)(F)F",
        "[P-2](F)(F)(F)(F)F",
        "[S-](F)(F)(F)(F)F",
        "[As](F)(F)(F)(F)F",
        "[Se](F)(F)(F)(F)(F)F",
        "I(F)(F)(F)(F)F",
        "[Cl-]",
        "[Br-]",
        "[I-]",
        "[O-][Cl+3]([O-])([O-])[O-]", // CHEMBL29216 counterion
        "[Fe+2]",
        "*",
    ] {
        let molecule = parse_molecule(smiles).unwrap().1;
        molecule.validate().unwrap();
    }
}

#[test]
fn ordinary_cations_do_not_require_a_same_period_valence_table_entry() {
    for smiles in [
        "[B+3]",
        "[Si+4]",
        "[As+3]",
        "[As+5]",
        "[Se+4]",
        "[Se+6]",
        "[Si+2]",
        "[C+4]",
        "[O+6]",
        "[P+3]",
        "[P+5]",
        "[S+4]",
        "[S+6]",
        "[C+6]",
        "[12C+4]",
        "[82Se+6]",
        "[Se+4].[Cl-].[Cl-].[Cl-].[Cl-]",
        "[H][B+2]",
        "[B+2][B+2]",
    ] {
        let molecule = parse_molecule(smiles).unwrap().1;
        molecule
            .validate()
            .unwrap_or_else(|error| panic!("{smiles}: {error}"));
    }
    for smiles in [
        "[C+7]",
        "[Si+15]",
        "[N+127]",
        "[Fe+27]",
        "[H+2]",
        "[N-46]",
        "[Si+4]C",
        "[Si+4]->[Cu]",
        "[B+2](F)F",
    ] {
        assert!(
            parse_molecule(smiles).unwrap().1.validate().is_err(),
            "{smiles}"
        );
    }
}

#[test]
fn aromatic_and_dative_valence_keep_supported_spellings() {
    for smiles in [
        "c1ccccc1",
        "C1=CC=CC=C1",
        "c1ccc2ccccc2c1",
        "n1ccccc1",
        "[nH]1cccc1",
        "N1C=CC=C1",
        "[n-]1cccc1",
        "[N-]1C=CC=C1",
        "C[n+]1ccccc1",
        "[nH+]1ccccc1",
        "o1cccc1",
        "s1cccc1",
        "O=c1cccc[nH]1",
        "c1ccccc1:c1ccccc1",
        "[NH3]->[Cu+2]",
        "[Cu+2]<-[NH3]",
        "n1ccccc1->[Cu]",
        "Cl[Pt@SP1](Cl)(<-[NH3])<-[NH3]",
        "FC=[C@]=CCl",
    ] {
        assert!(
            parse_molecule(smiles).unwrap().1.validate().is_ok(),
            "{smiles}"
        );
    }
    for smiles in [
        "Cn1(C)ccccc1",
        "c1(C)(C)ccccc1",
        "[nH2]1cccc1",
        "[NH3]<-[Cu]",
    ] {
        assert!(
            parse_molecule(smiles).unwrap().1.validate().is_err(),
            "{smiles}"
        );
    }
}

#[test]
fn aromatic_annotations_need_a_ring_and_a_kekule_structure() {
    for smiles in [
        "c1ccccc1",
        "c1ccc2ccccc2c1",
        "c1ccc2cccc2cc1",   // azulene
        "c1ccc-2cccc-2cc1", // azulene with its shared bond written single
        "c1ccc2[nH]ccc2c1",
        "c1ccn2cccc2c1", // indolizine: the bridgehead N gives a lone pair
        "Cn1c(=O)c2c(ncn2C)n(C)c1=O", // caffeine: exocyclic C=O carbons need no ring π bond
        "O=c1cccc[nH]1",
        "[cH-]1cccc1",
        "[cH+]1cccccc1",
        "b1ccccc1",  // borabenzene: B=C in its Kekulé structure
        "B1C=CC=C1", // borole
        "[bH]1cccc1",
        "c1cc[o+]cc1",
        "[se]1cccc1",
        "[te]1cccc1",
        "c1cc[n+]([O-])cc1",
        "c1ccc1", // anti-aromatic, but OpenSMILES writes it in lowercase
        "C1:C:C:C:C:C1",
        "c1ccccc1-c1ccccc1",
        "c1ccccc1c1ccccc1", // an aromatic bond outside rings is single
        "*c1ccccc1",
        // the wildcard is in the fused non-aromatic 4-ring, joined by single bonds
        "c1ccc2c(c1)C*2",
        "O=C1Nc2ccccc2C(=O)O/C(Cc2ccccc2)=N\\c2ccccc21", // CHEMBL2004566
    ] {
        let molecule = parse_molecule(smiles).unwrap().1;
        molecule
            .validate()
            .unwrap_or_else(|error| panic!("{smiles}: {error}"));
    }
    for (smiles, error) in [
        ("c1cccc1", "no Kekulé structure"), // cyclopentadiene's CH2 written as c
        ("c1ccnc1", "no Kekulé structure"), // pyrrole without its NH
        ("c1cc[nH]cc1", "no Kekulé structure"), // pyridine with an NH
        ("C=c1ccccc1", "no Kekulé structure"), // an exocyclic C=C takes one ring electron
        ("c1cccc1c1cccc1", "no Kekulé structure"), // the link is single, not a C=C
        ("CccccC", "aromatic atom outside a ring"),
        ("c1ccccc1:C", "aromatic atom outside a ring"),
        ("c", "aromatic atom outside a ring"),
        ("c1cc*cc1", "wildcard atom in an aromatic ring"),
        ("c1cc:*:cc1", "wildcard atom in an aromatic ring"),
        ("c1ccc2c(c1)c*2", "wildcard atom in an aromatic ring"),
        // ambiguous, so unsupported: the 4-ring's only other atoms are aromatic
        ("c1ccc2c(c1)**2", "wildcard atom in an aromatic ring"),
    ] {
        let molecule = parse_molecule(smiles).unwrap().1;
        let message = molecule.validate().unwrap_err();
        assert!(message.contains(error), "{smiles}: {message}");
        assert!(message.starts_with("atom "), "{smiles}: {message}");
    }
}

#[test]
fn chalcogen_cations_give_a_lone_pair_only_with_three_neighbours() {
    for smiles in [
        // three neighbours: a lone pair, like the neutral atom in thiophene
        "C[s+]1cccc1",
        "C[S+]1C=CC=C1",
        "C[se+]1cccc1",
        "C[Se+]1C=CC=C1",
        "C[te+]1cccc1",
        "CN(C)Cc1ccc(CSCCNc2n[s+]([O-])nc2N)o1 CHEMBL8982",
        // two neighbours: a ring double bond, like pyridine's N
        "c1cc[s+]cc1",
        "c1cc[se+]cc1",
    ] {
        let molecule = parse_molecule(smiles).unwrap().1;
        molecule
            .validate()
            .unwrap_or_else(|error| panic!("{smiles}: {error}"));
    }
    for smiles in [
        // a two-neighbour cation still needs a ring double bond, and has no partner here
        "c1c[s+]cc1",
        "c1c[se+]cc1",
        // a three-neighbour cation gives its lone pair, leaving five carbons unpaired
        "C[s+]1ccccc1",
        "C[se+]1ccccc1",
    ] {
        let message = parse_molecule(smiles).unwrap().1.validate().unwrap_err();
        assert!(
            message.contains("no Kekulé structure"),
            "{smiles}: {message}"
        );
    }
}

#[test]
fn loading_always_checks_syntax_and_arithmetic() {
    let overflow = format!("C{}", "(C)".repeat(256));
    let bracket_overflow = format!("[C]{}", "(C)".repeat(256));
    for smiles in ["C(", "C~C", "C1CC", &overflow] {
        assert!(parse_molecule(smiles).is_err(), "{smiles}");
    }
    assert!(parse_molecule("[CH255]").is_ok());
    for smiles in ["[CH255]C", &bracket_overflow] {
        assert!(
            parse_molecule(smiles).unwrap().1.validate().is_err(),
            "{smiles}"
        );
    }
    assert!(parse_molecule(&format!("C{}", "(C)".repeat(255))).is_ok());
}

#[test]
fn cli_validation_runs_after_loading_only_when_requested() {
    use std::process::Command;
    let path =
        std::env::temp_dir().join(format!("ringo-validation-cli-{}.smi", std::process::id()));
    let filename = path.to_str().unwrap();
    let invalid = "C(C)(C)(C)(C)C";
    std::fs::write(
        &path,
        format!("CCO ethanol\n{invalid} overvalent\nC( malformed\n"),
    )
    .unwrap();
    for validate in [false, true] {
        let mut index = Command::new(env!("CARGO_BIN_EXE_ringo-index"));
        if validate {
            index.arg("--validate");
        }
        index.args([filename, "2", "512", "2"]);
        let output = index.output().unwrap();
        assert!(output.status.success(), "{output:?}");
        let expected_skipped = if validate { 2 } else { 1 };
        assert!(
            String::from_utf8_lossy(&output.stderr)
                .contains(&format!("skipped {expected_skipped} line")),
            "{output:?}"
        );
        for query in ["CCO", invalid] {
            let mut search = Command::new(env!("CARGO_BIN_EXE_ringo-search"));
            if validate {
                search.arg("--validate");
            }
            search.args([filename, query, "0", "10"]);
            let output = search.output().unwrap();
            if validate && query == invalid {
                assert!(!output.status.success(), "{output:?}");
                assert!(
                    String::from_utf8_lossy(&output.stderr).contains("query validation failed"),
                    "{output:?}"
                );
            } else {
                assert!(output.status.success(), "{output:?}");
                assert_eq!(
                    String::from_utf8_lossy(&output.stdout).lines().count(),
                    if validate { 1 } else { 2 }
                );
            }
        }
    }
    std::fs::remove_file(filename.to_owned() + ".fp").unwrap();
    std::fs::remove_file(path).unwrap();
}
