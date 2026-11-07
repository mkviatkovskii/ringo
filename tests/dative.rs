pub mod common;

use petgraph::stable_graph::NodeIndex;
use ringo::io::smiles::reader::molecule::parse_molecule;
use ringo::model::bond::BondOrder;

#[test]
fn donor_endpoints_and_hydrogens_are_preserved() {
    for (smiles, donor, hydrogens) in [
        ("N->[Cu]", 0, 3),
        ("N<-[Cu]", 1, 2),
        ("[Cu]<-N", 1, 3),
        ("N->1CC[Cu]1", 0, 2),
    ] {
        let molecule = parse_molecule(smiles).unwrap().1;
        let bond = molecule
            .graph()
            .edge_weights()
            .find(|b| b.order == BondOrder::Dative)
            .unwrap();
        assert_eq!(bond.donor, Some(NodeIndex::new(donor)), "{smiles}");
        let nitrogen = molecule
            .graph()
            .node_weights()
            .find(|a| a.element.atomic_number == 7)
            .unwrap();
        assert_eq!(nitrogen.implicit_hs, hydrogens, "{smiles}");
        assert!(molecule
            .graph()
            .edge_weights()
            .filter(|b| b.order != BondOrder::Dative)
            .all(|b| b.donor.is_none()));
    }
    let molecule = parse_molecule("c1ccn(->[Cu])cc1").unwrap().1;
    let nitrogen = molecule
        .graph()
        .node_weights()
        .find(|a| a.element.atomic_number == 7)
        .unwrap();
    assert_eq!(nitrogen.implicit_hs, 0);
    assert_eq!(
        molecule
            .graph()
            .edge_weights()
            .filter(|b| b.order == BondOrder::Aromatic)
            .count(),
        6
    );
}

#[test]
fn coordination_rings_do_not_become_covalent_aromatic_rings() {
    let molecule = parse_molecule("N1C=CC=C->[B]<-1").unwrap().1;
    assert!(molecule.graph().node_weights().all(|a| !a.aromatic));
    assert_eq!(
        molecule
            .graph()
            .edge_weights()
            .filter(|b| b.order == BondOrder::Dative)
            .count(),
        2
    );
}

#[test]
fn dative_ligands_keep_atom_stereo_sites() {
    for (smiles, sites) in [
        ("[N@](C)(F)(Cl)->[Cu]", 4),
        ("[As@TB5](F)(Cl)(Br)(I)->[Cu]", 5),
        ("[Co@OH1](N)(N)(N)(N)(N)<-N", 6),
    ] {
        let molecule = parse_molecule(smiles).unwrap().1;
        assert_eq!(molecule.atom_stereo()[0].ligands.len(), sites);
        assert_eq!(
            molecule
                .graph()
                .edge_weights()
                .filter(|b| b.order == BondOrder::Dative)
                .count(),
            1
        );
    }
}

#[test]
fn dative_fingerprints_preserve_direction_under_reordering() {
    for (a, b) in [
        ("N->[Cu]", "[Cu]<-N"),
        ("N->1CC[Cu]1", "N1CC[Cu]<-1"),
        ("N->1CC[Cu]<-1", "N1CC[Cu]<-1"),
        (
            "c1cccn->2c1-c1n->3cccc1.[Cu]23(Cl)Cl",
            "[Cl][Cu]1([Cl])<-[n]2ccccc2-c2cccc[n]->12",
        ),
    ] {
        assert_eq!(
            parse_molecule(a).unwrap().1.ecfp(2, 4096).0,
            parse_molecule(b).unwrap().1.ecfp(2, 4096).0,
            "{a} {b}"
        );
    }
    let fingerprint = |s| parse_molecule(s).unwrap().1.ecfp(2, 4096).0;
    assert_ne!(fingerprint("[NH3]->[Cu]"), fingerprint("[NH3]<-[Cu]"));
    assert_ne!(fingerprint("[NH3]->[Cu]"), fingerprint("[NH3]-[Cu]"));
}

#[test]
fn malformed_dative_bonds_are_rejected() {
    for smiles in [
        "N->",
        "N<-",
        "->N",
        "N-><-N",
        "N->(C)",
        "N->1CC[Cu]->1",
        "N->1CC[Cu]-1",
    ] {
        assert!(parse_molecule(smiles).is_err(), "{smiles}");
    }
}

#[test]
fn dative_corpus_records_are_indexed_and_searchable() {
    let records = include_str!("data/dative.smi");
    common::assert_search("dative", records, records.lines());
}
