use super::{write, write_smiles, WriteError};
use crate::io::smiles::reader::molecule::parse_molecule;
use crate::model::molecule::Molecule;
use crate::model::stereo::{Configuration, DoubleBondStereo, Ligand};
use petgraph::stable_graph::NodeIndex;
use petgraph::visit::{EdgeRef, IntoEdgeReferences, NodeIndexable};

fn parse(smiles: &str) -> Molecule {
    parse_molecule(smiles)
        .unwrap_or_else(|e| panic!("{smiles}: {e:?}"))
        .1
}

macro_rules! ensure {
    ($condition:expr, $($message:tt)*) => {
        if !$condition {
            return Err(format!($($message)*));
        }
    };
}

/// Whether reordering `a` into `b` is an odd permutation.
fn odd(a: &[Ligand], b: &[Ligand]) -> Result<bool, String> {
    let positions = a
        .iter()
        .map(|ligand| b.iter().position(|other| other == ligand))
        .collect::<Option<Vec<usize>>>()
        .ok_or_else(|| format!("ligands {a:?} and {b:?}"))?;
    let inversions = (0..positions.len())
        .flat_map(|i| (i + 1..positions.len()).map(move |j| (i, j)))
        .filter(|&(i, j)| positions[i] > positions[j])
        .count();
    Ok(inversions % 2 == 1)
}

/// `stereo`'s configuration for `references` at `ends`, given in either order.
fn same_side(stereo: &DoubleBondStereo, ends: [NodeIndex; 2], references: [NodeIndex; 2]) -> bool {
    let references = if stereo.ends == ends {
        references
    } else {
        [references[1], references[0]]
    };
    let flips = (0..2)
        .filter(|&i| references[i] != stereo.references[i])
        .count();
    stereo.same_side == (flips % 2 == 0)
}

/// Atom `order[i]` of a molecule, by node index, is atom `i` of the one read back.
struct Mapping(Vec<NodeIndex>);

impl Mapping {
    fn new(source: &Molecule, order: &[NodeIndex]) -> Mapping {
        let mut map = vec![NodeIndex::end(); source.graph.node_bound()];
        for (i, node) in order.iter().enumerate() {
            map[node.index()] = NodeIndex::new(i);
        }
        Mapping(map)
    }

    fn node(&self, node: NodeIndex) -> NodeIndex {
        self.0[node.index()]
    }

    fn ligands(&self, ligands: &[Ligand]) -> Vec<Ligand> {
        ligands
            .iter()
            .map(|&ligand| match ligand {
                Ligand::Atom(node) => Ligand::Atom(self.node(node)),
                other => other,
            })
            .collect()
    }
}

/// Whether `b` has `a`'s atoms and bonds under `map`.
fn same_graph(a: &Molecule, b: &Molecule, map: &Mapping) -> Result<(), String> {
    ensure!(b.count_atoms() == a.count_atoms(), "atom count");
    ensure!(b.count_bonds() == a.count_bonds(), "bond count");
    let key = |molecule: &Molecule, node: NodeIndex| {
        let atom = &molecule.graph[node];
        (
            atom.element.atomic_number,
            atom.isotope,
            atom.atom_class,
            atom.charge,
            atom.hs + atom.implicit_hs,
            atom.aromatic,
            molecule.input_valences[node.index()],
        )
    };
    for node in a.graph.node_indices() {
        ensure!(
            key(a, node) == key(b, map.node(node)),
            "atom {}: {:?} became {:?}",
            node.index(),
            key(a, node),
            key(b, map.node(node))
        );
    }
    for edge in a.graph.edge_references() {
        let (from, to) = (edge.source(), edge.target());
        let bond = b.get_bond_by_atoms(map.node(from), map.node(to));
        let expected = (
            edge.weight().order,
            edge.weight().donor.map(|n| map.node(n)),
        );
        ensure!(
            bond.map(|bond| (bond.order, bond.donor)) == Some(expected),
            "bond {}-{}: {expected:?} became {:?}",
            from.index(),
            to.index(),
            bond.map(|bond| (bond.order, bond.donor))
        );
    }
    Ok(())
}

/// Whether `b` has `a`'s configurations under `map`.
fn same_stereo(a: &Molecule, b: &Molecule, map: &Mapping) -> Result<(), String> {
    ensure!(a.stereo.len() == b.stereo.len(), "atom stereo count");
    for s in &a.stereo {
        let center = map.node(s.center);
        let Some(t) = b.stereo.iter().find(|t| t.center == center) else {
            return Err(format!("centre {} lost", s.center.index()));
        };
        let source = map.ligands(&s.ligands);
        let same = match (s.configuration, t.configuration) {
            (Configuration::Tetrahedral(before), Configuration::Tetrahedral(after)) => {
                (before == after) != odd(&source, &t.ligands)?
            }
            (before, after) => before == after && source == t.ligands,
        };
        ensure!(same, "centre {} changed", s.center.index());
    }
    ensure!(
        a.allene_stereo.len() == b.allene_stereo.len(),
        "allene count"
    );
    for s in &a.allene_stereo {
        let center = map.node(s.center);
        let Some(t) = b.allene_stereo.iter().find(|t| t.center == center) else {
            return Err(format!("allene {} lost", s.center.index()));
        };
        // the terminals' hydrogens told apart by putting the terminal in their place
        let tagged = |ligands: Vec<Ligand>, axis: &[NodeIndex]| -> Vec<Ligand> {
            ligands
                .into_iter()
                .enumerate()
                .map(|(i, ligand)| match ligand {
                    Ligand::Hydrogen if i < 2 => Ligand::Atom(axis[0]),
                    Ligand::Hydrogen => Ligand::Atom(axis[axis.len() - 1]),
                    other => other,
                })
                .collect()
        };
        let axis: Vec<NodeIndex> = s.axis.iter().map(|&n| map.node(n)).collect();
        let source = tagged(map.ligands(&s.ligands), &axis);
        let target = tagged(t.ligands.to_vec(), &t.axis);
        ensure!(
            (s.permutation == t.permutation) != odd(&source, &target)?,
            "allene {} changed",
            s.center.index()
        );
    }
    ensure!(
        a.bond_stereo.len() == b.bond_stereo.len(),
        "double-bond stereo count"
    );
    for s in &a.bond_stereo {
        let ends = s.ends.map(|n| map.node(n));
        let Some(t) = b
            .bond_stereo
            .iter()
            .find(|t| t.ends == ends || t.ends == [ends[1], ends[0]])
        else {
            return Err("double bond lost".to_owned());
        };
        ensure!(
            same_side(t, ends, s.references.map(|n| map.node(n))) == s.same_side,
            "double bond changed"
        );
    }
    Ok(())
}

/// Writes `molecule` in the order of `ranks`, and checks that it reads back the same.
fn check(molecule: &Molecule, ranks: &[usize]) -> Result<String, String> {
    let (written, order) = write(molecule, ranks).map_err(|e| format!("{e:?}"))?;
    let context = |e: String| format!("{written}: {e}");
    let read = parse_molecule(&written)
        .map_err(|e| context(format!("unreadable: {e:?}")))?
        .1;
    let map = Mapping::new(molecule, &order);
    same_graph(molecule, &read, &map).map_err(context)?;
    same_stereo(molecule, &read, &map).map_err(context)?;
    Ok(written)
}

/// A small deterministic shuffle, so failures repeat.
fn ranks(count: usize, seed: u64) -> Vec<usize> {
    let mut state = seed.wrapping_mul(0x9e37_79b9_7f4a_7c15) | 1;
    let mut ranks: Vec<usize> = (0..count).collect();
    for i in (1..count).rev() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        ranks.swap(i, usize::try_from(state % (i as u64 + 1)).unwrap());
    }
    ranks
}

fn in_reading_order(molecule: &Molecule) -> Vec<usize> {
    (0..molecule.graph.node_bound()).collect()
}

/// Writes `smiles` in reading order and in shuffled orders, and reads each one back.
fn assert_round_trip(smiles: &str) {
    let molecule = parse(smiles);
    let written = check(&molecule, &in_reading_order(&molecule))
        .unwrap_or_else(|e| panic!("{smiles} -> {e}"));
    // writing what was read gives the same string again
    assert_eq!(write_smiles(&parse(&written)).unwrap(), written, "{smiles}");
    for seed in 0..20 {
        match check(&molecule, &ranks(molecule.graph.node_bound(), seed)) {
            Ok(_) => {}
            // only centres whose permutation number isn't renumbered yet may fail
            Err(e) if e == format!("{:?}", WriteError::UnsupportedStereo) => assert!(
                molecule
                    .stereo
                    .iter()
                    .any(|s| !matches!(s.configuration, Configuration::Tetrahedral(_))),
                "{smiles} (seed {seed})"
            ),
            Err(e) => panic!("{smiles} (seed {seed}) -> {e}"),
        }
    }
}

/// Writes every record of `RINGO_SMILES_CORPUS`'s files in reading order and in one shuffled
/// order, and reads each one back. Missing files are skipped.
#[test]
#[ignore = "needs RINGO_SMILES_CORPUS with the ChEBI and ChEMBL SMILES files"]
fn corpus_round_trip() {
    use rayon::prelude::*;
    use std::collections::BTreeMap;
    use std::io::BufRead;
    let root = std::env::var("RINGO_SMILES_CORPUS").unwrap();
    let mut failed = 0;
    for name in ["chebi_20261105.smi", "chembl_37.smi"] {
        let path = std::path::Path::new(&root).join(name);
        let Ok(file) = std::fs::File::open(&path) else {
            println!("{name}: missing");
            continue;
        };
        let start = std::time::Instant::now();
        let (mut records, mut unreadable) = (0, 0);
        let mut errors: BTreeMap<String, usize> = BTreeMap::new();
        let mut examples = Vec::new();
        let mut lines = std::io::BufReader::new(file)
            .lines()
            .map(Result::unwrap)
            .enumerate();
        loop {
            let batch: Vec<(usize, String)> = lines.by_ref().take(1 << 16).collect();
            if batch.is_empty() {
                break;
            }
            let results: Vec<(String, Option<Vec<String>>)> = batch
                .into_par_iter()
                .map(|(number, line)| {
                    let Ok((_, molecule)) = parse_molecule(&line) else {
                        return (line, None);
                    };
                    let shuffled = ranks(molecule.graph.node_bound(), number as u64);
                    let outcomes = [
                        ("reading order", in_reading_order(&molecule)),
                        ("shuffled", shuffled),
                    ]
                    .into_iter()
                    .filter_map(|(order, ranks)| {
                        check(&molecule, &ranks)
                            .err()
                            .map(|e| format!("{order}: {e}"))
                    })
                    .collect();
                    (line, Some(outcomes))
                })
                .collect();
            for (line, outcomes) in results {
                let Some(outcomes) = outcomes else {
                    unreadable += 1;
                    continue;
                };
                records += 1;
                for e in outcomes {
                    if e.split_once(": ")
                        .is_some_and(|(_, e)| e.chars().all(char::is_alphanumeric))
                    {
                        *errors.entry(e).or_default() += 1;
                    } else {
                        failed += 1;
                        if examples.len() < 200 {
                            examples.push(format!("{line} -> {e}"));
                        }
                    }
                }
            }
        }
        println!(
            "{name}: {records} records ({unreadable} unreadable) in {:.1?}, write errors \
             {errors:?}",
            start.elapsed()
        );
        for example in examples {
            println!("  {example}");
        }
    }
    assert_eq!(failed, 0);
}

#[test]
fn writes_what_it_reads() {
    for (smiles, expected) in [
        ("", ""),
        ("C", "C"),
        ("CCO", "CCO"),
        ("C(C)(C)C", "C(C)(C)C"),
        ("C=C", "C=C"),
        ("C#N", "C#N"),
        ("[Na+].[Cl-]", "[Na+].[Cl-]"),
        ("C1CC1", "C1CC1"),
        ("C1CC1C1CC1", "C1CC1C1CC1"),
        ("C12CCC1CC2", "C12CCC1CC2"),
        ("C1=CC=CC=C1", "c1ccccc1"),
        ("c1ccccc1-c1ccccc1", "c1ccccc1-c1ccccc1"),
        ("c1ccccc1=O", "c1ccccc1=O"),
        ("C1=CC=CN1", "c1ccc[nH]1"),
        ("[nH]1cccc1", "[nH]1cccc1"),
        ("C1=CC=C[Si]=C1", "c1ccc:[Si]:c1"),
        ("[CH2]C", "[CH2]C"),
        ("[13CH3:2][NH3+]", "[13CH3:2][NH3+]"),
        ("[Fe+3]", "[Fe+3]"),
        ("[O-2]", "[O-2]"),
        ("C([H])([H])([H])H", "C([H])([H])([H])[H]"),
        ("[H][H]", "[H][H]"),
        ("N->[Cu]", "N->[Cu]"),
        ("[Cu]<-N", "[Cu]<-N"),
        ("*C", "*C"),
        ("F/C=C/F", "F/C=C/F"),
        ("F/C=C\\F", "F/C=C\\F"),
        ("N[C@@H](C)C(=O)O", "N[C@@H](C)C(=O)O"),
        ("[C@H](F)(Cl)Br", "[C@H](F)(Cl)Br"),
        ("NC(Br)=[C@]=C(O)C", "NC(Br)=[C@]=C(O)C"),
    ] {
        let written = write_smiles(&parse(smiles)).unwrap();
        assert_eq!(written, expected, "{smiles}");
    }
}

#[test]
fn ring_digits_are_reused_and_go_past_nine() {
    // a digit closed on an atom isn't reused on that atom, but is on the next one
    assert_eq!(write_smiles(&parse("C1CC12CC2")).unwrap(), "C1CC12CC2");
    assert_eq!(write_smiles(&parse("C1CC1C2CC2")).unwrap(), "C1CC1C1CC1");
    assert_round_trip("C12C3C4C1C5C2C3C45");
    // two chains of 12 atoms, joined atom for atom
    let rings = "C%10C%11C%12C%13C%14C%15C%16C%17C%18C%19C%20C%21";
    let written = write_smiles(&parse(&rings.repeat(2))).unwrap();
    assert!(written.contains("%10"), "{written}");
    assert_round_trip(&rings.repeat(2));
}

#[test]
fn round_trips() {
    for smiles in [
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "CC(C1=CC2=C(C=C1)C=C(C=C2)OC)C(=O)O",
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "c1ccc2ccccc2c1",
        "c1ccc2c(c1)[nH]c1ccccc12",
        "O=c1cc[nH]cc1",
        "c1ccoc1",
        "c1ccsc1",
        "[se]1cccc1",
        "C1=CC=C[Si]=C1",
        "c1cc2cccc3ccc(c1)c23",
        // the ring pair that shares a single bond must stay the same in any order
        "Cc1cc2cc(=O)cc3oc(C)cc(o1)c23",
        // every shortest ring counts when several tie
        "p12p3p1[P+]23",
        "C1=CC=C2C=CC=C2C=C1",
        // Kekulé forms that perception must treat the same when the output is read back
        "C1=CC2=CC=CC=CC2=C1",
        "C1=CC=C2C(=C1)C1=CC=CC=C12",
        "C1C2=CC=CC=C2C2=CC=CC=C12",
        "C=C1C=CC=C1",
        "C1C=CC=C1",
        "O=C1C=CC=CN1",
        "O=C1C=CC(=O)C=C1",
        "C1=CC=CN1",
        "C1=CC=C2C(=C1)C=CN2",
        "[CH+]1C=CC=CC=C1",
        "[CH-]1C=CC=C1",
        "B1C=CC=C1",
        "C1=CC2=CC3=CC=C(N3)C=C4C=CC(=N4)C=C5C=CC(=N5)C=C1N2",
        "C1=CC=C2C=C3C=CC=CC3=CC2=C1",
        "C1=CC=CC=CC=C1",
        "[Te]1C=CC=C1",
        "O=C1NC=CC2=CC=CC=C12",
        "C12=CC=CC=C1C=CC=C2",
        "c1ccc2c(c1)-c1ccccc1-2",
        "[O-][n+]1ccccc1",
        "Cc1c[nH]c(=O)[nH]c1=O",
        "c1ccc(cc1)-c1ccccc1",
        "C1CC2CCC1C2",
        "OC(=O)CC(O)(CC(O)=O)C(O)=O",
        "[2H]C([2H])([2H])O",
        "[NH4+].[O-]C(=O)C",
        "CCO.O.O",
        "C1.C1",
        "B(C)(C)C",
        "[Na+].[Na+].[O-]S(=O)(=O)[O-]",
        "CN(C)C1[S][Fe]23([S]C(=[S]->2)N(C)C)([S]C(=[S]->3)N(C)C)<-[S]=1",
        "[*:0]c1cc(-c2[o+]c3cc([O-])cc(O)c3cc2O)cc([*:0])c1O",
        "[Cl:1][c:2]1[n:3][c:4]2[c:5]([H:24])[c:6]([H:25])[c:7]([H:26])[c:8]([H:27])[c:9]2[c:10]([H:28])[c:11]1[C:12]([O:13][c:14]1[c:15]([H:31])[c:16]([H:32])[c:17]([H:33])[c:18]2[c:19]([H:34])[c:20]([H:35])[c:21]([H:36])[n:22][c:23]12)([H:29])[H:30]",
        // tetrahedral centres, with H, lone pairs, ring bonds and explicit hydrogens
        "C[C@H](O)CO",
        "C[C@@H](O)CO",
        "OC1C[N@]2CC[C@H]1CC2",
        "CC(C)CN[S@](C)=O",
        "O=[S@@](O)C",
        "[H][C@](F)(Cl)Br",
        "[2H][C@](F)(Cl)Br",
        "F[C@]1(Cl)CCCC1",
        "[C@@H]12CCC[C@H]1CCC2",
        "C[C@]12CC[C@@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC1=CC(=O)CC[C@]31C",
        "ClCC[N@+]12CC[N@@+](CCCl)(CC1)C2",
        "CCO[P@@](=O)(C#N)N(C)C",
        "O=[S@@](O)[S@@](=O)O",
        "C[N+]1(C)C(c2ccccc2)C(=O)O[B@-]1(F)c1ccccc1",
        "O[Si@](CCCN1CCCCC1)(c1ccc(F)cc1)C1CCCCC1",
        // allenes and cumulenes
        "NC(Br)=[C@]=C(O)C",
        "OC=[C@]=CF",
        "OC([H])=[C@AL1]=C([H])F",
        "CC=C=[C@@]=C=CC",
        "C1CC=[C@]=C(C1)C",
        // double bonds, in chains, rings, conjugated and next to aromatic rings
        "F/C=C/F",
        "F/C=C\\F",
        "C/C=C/C=C/C",
        "C/C=C\\C=C/C",
        "F/C=C/C=C",
        "C/C(F)=C(\\Cl)Br",
        "C(/F)=C/F",
        "C1CCCCCC/C=C/1",
        "C/1=C/CCCCCC1",
        "c1ccccc1/C=C/c1ccccc1",
        "O/N=C/c1ccccc1",
        "C/C=C/[C@H](O)C",
        // two configured double bonds with an unconfigured one between them
        "C/C=C(/F)C=CC(/F)=C/C",
        // one mark between two double bonds serves both
        r"O=C/C=C(O)\C=C/C(=O)O",
        r"O=C(O)/C(Cl)=C(Cl)\C(Cl)=C(\Cl)C(=O)O",
        r"CC(/C=C/C=C(\C)C(=O)O)=C\C=C\C=C(C)\C=C\C=C(/C)C(=O)O",
        // a mark that would configure the ring's double bonds is avoided
        r"NC1=C/C(=N\c2ccc(N)cc2)C(N)=C/C1=N\c1ccc(N)cc1",
        // square planar, trigonal bipyramidal and octahedral centres
        "Cl[Pt@SP1](Cl)(<-[NH3])<-[NH3]",
        "Cl[Pt@SP2](Cl)(<-[NH3])<-[NH3]",
        "S[As@TB5](F)(N)(Cl)Br",
        "O[Mn@OH28](Cl)(C)(N)F",
        "CC(C)O[P@](C)(=O)F",
    ] {
        assert_round_trip(smiles);
    }
}

#[test]
fn every_test_fixture_round_trips() {
    for file in [
        include_str!("../../../../tests/data/stereo.smi"),
        include_str!("../../../../tests/data/molecules.smi"),
        include_str!("../../../../tests/data/dative.smi"),
        include_str!("../../../../tests/data/atom_maps.smi"),
    ] {
        for line in file.lines() {
            assert_round_trip(line.split_whitespace().next().unwrap());
        }
    }
}

#[test]
fn reordered_centres_without_a_renumbering_are_errors() {
    let molecule = parse("Cl[Pt@SP1](Cl)(<-[NH3])<-[NH3]");
    let reversed: Vec<usize> = (0..molecule.count_atoms()).rev().collect();
    assert_eq!(
        write(&molecule, &reversed).err(),
        Some(WriteError::UnsupportedStereo)
    );
}

#[test]
fn aromatic_valence_that_notation_cannot_give_back_is_an_error() {
    // the nitrogen has valence four; written `n` it would read back with three
    assert_eq!(
        write_smiles(&parse("[Cu][N]1=CC=CC=C1")).err(),
        Some(WriteError::AromaticValence)
    );
    // as in uranyl phthalocyanine, where the change makes two more rings aromatic
    assert_eq!(
        write_smiles(&parse(
            "[O]=[U]1234(=[O])[N]5=C6N=C7c8ccccc8C(=[N]71)N=c1c7ccccc7c([n]12)=NC1=[N]3\
             C(=Nc2c3ccccc3c([n]24)N=C5c2ccccc26)c2ccccc21"
        ))
        .err(),
        Some(WriteError::AromaticValence)
    );
}
