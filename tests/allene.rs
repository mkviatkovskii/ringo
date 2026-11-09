pub mod common;

use ringo::io::smiles::reader::molecule::parse_molecule;
use ringo::model::stereo::Ligand;

#[expect(
    clippy::unwrap_used,
    reason = "test helper requires a resolved allene axis"
)]
fn orientation(smiles: &str) -> u8 {
    let molecule = parse_molecule(smiles)
        .unwrap_or_else(|_| panic!("{smiles}"))
        .1;
    assert!(molecule.atom_stereo().is_empty(), "{smiles}");
    assert!(molecule.double_bond_stereo().is_empty(), "{smiles}");
    assert_eq!(molecule.allene_stereo().len(), 1, "{smiles}");
    let stereo = &molecule.allene_stereo()[0];
    assert_eq!(stereo.axis[stereo.axis.len() / 2], stereo.center);
    assert!(stereo.axis.len() >= 3 && stereo.axis.len() % 2 == 1);
    for pair in stereo.axis.windows(2) {
        assert_eq!(
            molecule.get_bond_by_atoms(pair[0], pair[1]).unwrap().order,
            ringo::model::bond::BondOrder::Double
        );
    }
    let labels = stereo.ligands.map(|ligand| match ligand {
        Ligand::Atom(node) => {
            let atom = molecule.get_atom(node).unwrap();
            (atom.element.atomic_number, atom.isotope)
        }
        Ligand::Hydrogen => (1, 0),
        _ => panic!("{smiles}"),
    });
    // Swapping a terminal pair inverts the sense; reversing the two pairs does not.
    (stereo.permutation - 1) ^ u8::from(labels[0] > labels[1]) ^ u8::from(labels[2] > labels[3])
}

#[test]
fn terminal_permutations_and_axis_reversal_preserve_orientation() {
    // OpenSMILES' collapsed-tetrahedron example; Daylight defines AL1/2 aliases.
    let base = orientation("NC(Br)=[C@]=C(O)C");
    for smiles in [
        "NC(Br)=[C@AL1]=C(O)C",
        "BrC(N)=[C@@]=C(O)C",
        "NC(Br)=[C@@]=C(C)O",
        "CC(O)=[C@]=C(Br)N",
        "[C@](=C(N)Br)=C(O)C",
    ] {
        assert_eq!(orientation(smiles), base, "{smiles}");
    }
    assert_ne!(orientation("NC(Br)=[C@AL2]=C(O)C"), base);
    assert_ne!(orientation("NC(Br)=[C@]=C(C)O"), base);
    let molecule = parse_molecule("NC(Br)=[C@AL2]=C(O)C").unwrap().1;
    assert_eq!(molecule.allene_stereo()[0].permutation, 2);
    assert_eq!(
        molecule.allene_stereo()[0].ligands.map(|ligand| {
            let Ligand::Atom(node) = ligand else { panic!() };
            molecule.get_atom(node).unwrap().element.atomic_number
        }),
        [7, 35, 8, 6]
    );
    assert_eq!(
        parse_molecule("NC(Br)=C=C(O)C").unwrap().1.allene_stereo(),
        []
    );
}

#[test]
fn terminal_hydrogens_keep_smiles_positions_and_isotopes() {
    // Daylight explicitly gives this implicit/explicit-H equivalence.
    let base = orientation("OC=[C@]=CF");
    for smiles in [
        "OC([H])=[C@AL1]=C([H])F",
        "O[CH]=[C@]=[CH]F",
        "C(O)=[C@@]=CF",
        "[CH](O)=[C@@]=CF",
        "FC=[C@]=CO",
    ] {
        assert_eq!(orientation(smiles), base, "{smiles}");
    }
    assert_eq!(
        orientation("[2H]C(F)=[C@]=C(Cl)Br"),
        orientation("FC([2H])=[C@@]=C(Cl)Br")
    );
    assert_ne!(
        orientation("[2H]C([H])=[C@]=C(F)Cl"),
        orientation("[H]C([2H])=[C@]=C(F)Cl")
    );
    let m = parse_molecule("[2H]C(F)=[C@]=C(Cl)Br").unwrap().1;
    let Ligand::Atom(h) = m.allene_stereo()[0].ligands[0] else {
        panic!()
    };
    assert_eq!(m.get_atom(h).unwrap().isotope, 2);
    // Preserve supplied annotations even for chemically degenerate terminal pairs.
    let m = parse_molecule("C=[C@]=C").unwrap().1;
    assert_eq!(m.allene_stereo()[0].ligands, [Ligand::Hydrogen; 4]);
}

#[test]
fn ring_slots_and_disconnected_roots_preserve_order() {
    assert_eq!(
        orientation("NC1=[C@]=C(O)[13CH2]CC1"),
        orientation("OC1=[C@]=C(N)CC[13CH2]1")
    );
    assert_eq!(
        orientation("F1.NC1=[C@]=C(O)Cl"),
        orientation("NC(F)=[C@]=C(O)Cl")
    );
    assert_eq!(
        orientation("F1.C1=[C@]=C(O)Cl"),
        orientation("[CH](F)=[C@]=C(O)Cl")
    );
    assert_eq!(
        orientation("FC(Cl)=[C@]=1.BrC1I"),
        orientation("FC(Cl)=[C@]=C(Br)I")
    );
}

#[test]
fn longer_even_bond_axes_keep_terminal_references() {
    for middle in ["=[C@]=", "=C=[C@]=C=", "=C=C=[C@]=C=C="] {
        let a = format!("FC(Cl){middle}C(Br)I");
        let b = format!("IC(Br){middle}C(Cl)F");
        assert_eq!(orientation(&a), orientation(&b), "{a} {b}");
        assert_ne!(orientation(&a), orientation(&a.replace('@', "@@")));
        let m = parse_molecule(&a).unwrap().1;
        assert_eq!(
            m.allene_stereo()[0].axis.len(),
            middle.matches('=').count() + 1
        );
        let ends = [
            m.allene_stereo()[0].axis[0],
            *m.allene_stereo()[0].axis.last().unwrap(),
        ];
        for (i, ligand) in m.allene_stereo()[0].ligands.iter().enumerate() {
            let Ligand::Atom(node) = ligand else { panic!() };
            assert!(m.has_bond(ends[i / 2], *node));
        }
    }
}

#[test]
fn malformed_tags_and_unsupported_axes_are_rejected() {
    for smiles in [
        "FC(Cl)=[C@AL]=C(Br)I",
        "FC(Cl)=[C@AL0]=C(Br)I",
        "FC(Cl)=[C@AL3]=C(Br)I",
        "FC(Cl)=[C@AL256]=C(Br)I",
        "FC(Cl)=[C@AL-1]=C(Br)I",
        "[C@AL1](F)(Cl)(Br)I",
        "FC(Cl)[C@AL1]=C(Br)I",
        "FC(Cl)=[C@AL1]#C(Br)I",
        "FC(Cl)=[C@H]=C(Br)I",
        "FC(Cl)=[C@](F)=C(Br)I",
        "FC(Cl)=[C@]=C=C(Br)I",
        "FC(Cl)=C=[C@]=C(Br)I",
        "FC(Cl)=[C@]=C(F)=C(Br)I",
        "FC(Cl)=C=[C@]=[CH]=C(Br)I",
        "C1=[C@]=C=1",
        "F[C]=[C@]=C(Br)I",
        "FC(Cl)=[C@]=O",
        "FC(Cl)=[C@]=C(Br)(I)F",
        "FC(Cl)=[C@]=C(Br)->[Cu]",
        "FC(Cl)=[C@]=C#N",
        "F/C(Cl)=[C@]=C(Br)I",
    ] {
        assert!(parse_molecule(smiles).is_err(), "{smiles}");
    }
}

#[test]
fn indexing_keeps_allene_fingerprints_stereo_insensitive() {
    common::assert_search(
        "allene",
        "NC(Br)=[C@AL1]=C(O)C allene\n",
        ["NC(Br)=[C@AL2]=C(O)C", "NC(Br)=C=C(O)C"],
    );
    let m = parse_molecule("NC(Br)=[C@]=C(O)C.N[C@H](C)C(=O)O.F/C=C/F")
        .unwrap()
        .1;
    assert_eq!(m.allene_stereo().len(), 1);
    assert_eq!(m.atom_stereo().len(), 1);
    assert_eq!(m.double_bond_stereo().len(), 1);
}
