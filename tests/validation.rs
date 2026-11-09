use ringo::db::index::{index_file, index_file_with_validation, DEFAULT_BATCH_LINES};
use ringo::db::search::{similarity_search_with_validation, SearchMode};
use ringo::io::smiles::reader::molecule::{parse_molecule, parse_molecule_with_validation};

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
    ] {
        assert!(parse_molecule(smiles).is_ok(), "permissive: {smiles}");
        assert!(
            parse_molecule_with_validation(smiles, false).is_ok(),
            "{smiles}"
        );
        assert!(
            parse_molecule_with_validation(smiles, true).is_err(),
            "validated: {smiles}"
        );
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
        "[O-][Cl+3]([O-])([O-])[O-]", // CHEMBL29216 counterion
        "[Fe+2]",
        "*",
    ] {
        let ordinary = parse_molecule(smiles).unwrap().1;
        let validated = parse_molecule_with_validation(smiles, true).unwrap().1;
        assert_eq!(
            ordinary.ecfp(2, 512).to_bytes(),
            validated.ecfp(2, 512).to_bytes(),
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
        "[NH3]<-[Cu]",
        "Cl[Pt@SP1](Cl)(<-[NH3])<-[NH3]",
        "FC=[C@]=CCl",
    ] {
        let expected = !matches!(smiles, "[NH3]<-[Cu]");
        assert_eq!(
            parse_molecule_with_validation(smiles, true).is_ok(),
            expected,
            "{smiles}"
        );
    }
    for smiles in ["Cn1(C)ccccc1", "c1(C)(C)ccccc1", "[nH2]1cccc1"] {
        assert!(parse_molecule(smiles).is_ok(), "{smiles}");
        assert!(
            parse_molecule_with_validation(smiles, true).is_err(),
            "{smiles}"
        );
    }
}

#[test]
fn syntax_and_arithmetic_errors_are_mandatory_in_both_modes() {
    let overflow = format!("C{}", "(C)".repeat(256));
    let bracket_overflow = format!("[C]{}", "(C)".repeat(256));
    for smiles in [
        "C(",
        "C~C",
        "C1CC",
        "[CH255]C",
        &overflow,
        &bracket_overflow,
    ] {
        for validate in [false, true] {
            assert!(
                parse_molecule_with_validation(smiles, validate).is_err(),
                "{smiles}"
            );
        }
        assert!(parse_molecule(smiles).is_err(), "{smiles}");
    }
    assert!(parse_molecule("[CH255]").is_ok());
    assert!(parse_molecule(&format!("C{}", "(C)".repeat(255))).is_ok());
}

#[test]
fn indexing_and_query_parsing_use_the_selected_policy() {
    let path = std::env::temp_dir().join(format!("ringo-validation-{}.smi", std::process::id()));
    let filename = path.to_str().unwrap();
    let invalid = "C(C)(C)(C)(C)C";
    std::fs::write(
        &path,
        format!("CCO ethanol\n{invalid} overvalent\nC( malformed\n"),
    )
    .unwrap();
    assert_eq!(
        index_file(filename, 2, 512, DEFAULT_BATCH_LINES).unwrap(),
        1
    );
    let search = |query, validate| {
        similarity_search_with_validation(filename, query, 0.0, 1.0, 10, SearchMode::Best, validate)
    };
    assert_eq!(search(invalid, false).unwrap().len(), 2);
    assert_eq!(
        search(invalid, true).err().unwrap().kind(),
        std::io::ErrorKind::InvalidInput
    );
    assert_eq!(
        index_file_with_validation(filename, 2, 512, DEFAULT_BATCH_LINES, true).unwrap(),
        2
    );
    let results = search("CCO", true).unwrap();
    assert_eq!(results.len(), 1);
    assert_eq!(results[0].line, "CCO ethanol\n");
    assert_eq!(search(invalid, false).unwrap().len(), 1);
    std::fs::remove_file(filename.to_owned() + ".fp").unwrap();
    std::fs::remove_file(path).unwrap();
}

#[test]
fn command_line_validation_is_opt_in() {
    use std::process::Command;
    let path =
        std::env::temp_dir().join(format!("ringo-validation-cli-{}.smi", std::process::id()));
    let filename = path.to_str().unwrap();
    let invalid = "C(C)(C)(C)(C)C";
    std::fs::write(&path, format!("CCO ethanol\n{invalid} overvalent\n")).unwrap();
    for validate in [false, true] {
        let mut index = Command::new(env!("CARGO_BIN_EXE_ringo-index"));
        index.arg(filename);
        if validate {
            index.arg("--validate");
        }
        let output = index.output().unwrap();
        assert!(output.status.success(), "{output:?}");
        assert_eq!(
            String::from_utf8_lossy(&output.stderr).contains("skipped 1 line"),
            validate
        );
        let mut search = Command::new(env!("CARGO_BIN_EXE_ringo-search"));
        search.args([filename, invalid, "0", "10"]);
        if validate {
            search.arg("--validate");
        }
        let output = search.output().unwrap();
        assert_eq!(output.status.success(), !validate, "{output:?}");
        if validate {
            assert!(String::from_utf8_lossy(&output.stderr).contains("could not parse query"));
        } else {
            assert_eq!(String::from_utf8_lossy(&output.stdout).lines().count(), 2);
        }
    }
    std::fs::remove_file(filename.to_owned() + ".fp").unwrap();
    std::fs::remove_file(path).unwrap();
}

#[test]
#[ignore = "requires downloaded ChEBI/ChEMBL files; set RINGO_SMILES_CORPUS"]
fn corpus_valence_coverage() {
    use std::io::BufRead;
    let root = std::env::var("RINGO_SMILES_CORPUS").unwrap();
    for name in ["chebi_20261105.smi", "chembl_37.smi"] {
        let path = std::path::Path::new(&root).join(name);
        let before = std::fs::metadata(&path).unwrap();
        let mut counts = [0usize; 3];
        let mut examples = Vec::new();
        for line in std::io::BufReader::new(std::fs::File::open(&path).unwrap()).lines() {
            let line = line.unwrap();
            let category = if parse_molecule_with_validation(&line, true).is_ok() {
                0
            } else if parse_molecule(&line).is_ok() {
                if examples.len() < 10 {
                    examples.push(line.clone());
                }
                1
            } else {
                2
            };
            counts[category] += 1;
        }
        let after = std::fs::metadata(&path).unwrap();
        assert_eq!(
            (before.len(), before.modified().unwrap()),
            (after.len(), after.modified().unwrap())
        );
        println!("{name}: [validated, valence rejected, rejected in both modes] = {counts:?}; examples={examples:?}");
    }
}
