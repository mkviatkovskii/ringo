pub mod common;

use ringo::io::smiles::reader::molecule::parse_molecule;

#[test]
fn atom_classes_preserve_numeric_labels_and_other_properties() {
    let molecule = parse_molecule("[13C@@H:005](F)([1*:4294967295])Cl.[NH4+:12]")
        .unwrap()
        .1;
    let atoms: Vec<_> = molecule.graph().node_weights().collect();
    assert_eq!(
        atoms.iter().map(|a| a.atom_class).collect::<Vec<_>>(),
        [5, 0, u32::MAX, 0, 12]
    );
    assert_eq!((atoms[0].isotope, atoms[0].hs), (13, 1));
    assert_eq!((atoms[2].element.atomic_number, atoms[2].isotope), (0, 1));
    assert_eq!((atoms[4].charge, atoms[4].hs), (1, 4));
    assert_eq!(molecule.atom_stereo().len(), 1);
    for smiles in ["C", "[C]", "[C:0]", "[C:000]"] {
        assert_eq!(
            parse_molecule(smiles)
                .unwrap()
                .1
                .graph()
                .node_weights()
                .next()
                .unwrap()
                .atom_class,
            0
        );
    }
}

#[test]
fn atom_classes_do_not_change_fingerprints_or_stereo() {
    for (mapped, plain) in [
        ("[CH3:12][OH:2]", "[CH3][OH]"),
        ("[*:0]C[1*:4294967295]", "[*]C[1*]"),
        ("[NH4+:005].[Cl-:4]", "[NH4+].[Cl-]"),
        ("[C@H:2]1([F:3])CC[O:4]1", "[C@H]1(F)CC[O]1"),
    ] {
        let mapped = parse_molecule(mapped).unwrap().1;
        let plain = parse_molecule(plain).unwrap().1;
        assert_eq!(mapped.ecfp(2, 2048).0, plain.ecfp(2, 2048).0);
        assert_eq!(mapped.atom_stereo(), plain.atom_stereo());
    }
}

#[test]
fn malformed_atom_classes_are_rejected() {
    for smiles in [
        "[C:]",
        "[C:-1]",
        "[C:+1]",
        "[C:1:2]",
        "[C:1H]",
        "[C:1+]",
        "[C:1@]",
        "[C:4294967296]",
        "C:1",
        "[C: 1]",
        "[C:1.0]",
    ] {
        assert!(parse_molecule(smiles).is_err(), "{smiles}");
    }
}

#[test]
fn mapped_corpus_records_are_indexed_and_searchable() {
    common::assert_self_search("atom_maps", include_str!("data/atom_maps.smi"));
}
