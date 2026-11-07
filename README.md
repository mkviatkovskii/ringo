# Ringo
File chemical search engine written in Rust.

The SMILES reader preserves tetrahedral (`@`, `@@`, `@TH1`, `@TH2`),
trigonal-bipyramidal (`@TB1`–`@TB20`), octahedral (`@OH1`–`@OH30`),
and ordinary double-bond (`/`, `\`) configurations. Use
`Molecule::atom_stereo()` and `Molecule::double_bond_stereo()` to inspect them.
Atom configurations refer to the recorded ligand order, not CIP R/S labels.
For TB/OH, permutation numbers retain their SMILES meaning; missing coordination
sites occupy the implicit-hydrogen position, following RDKit's SMILES convention.
Three-neighbour N/P/S/As/Se tetrahedral records use a final virtual lone-pair site.
Absence of a record means unspecified stereochemistry. The reader does not
validate chemical plausibility or perceive additional stereocentres.

ECFP fingerprints and similarity search remain **stereo-insensitive**: opposite
configurations can have identical fingerprints. Index format remains version 1.
Rebuild development indexes to include stereo records that older readers skipped.
Allene/cumulene and square-planar stereochemistry remain unsupported. Invalid,
incomplete, unsupported, or contradictory specifications produce parsing errors.
Directional specifications that aromaticity normalization cannot preserve are
also rejected rather than silently discarded.

`tests/data/stereo.smi` contains small examples with source ChEBI/ChEMBL IDs from
the supplied `chebi_20261105.smi` and `chembl_37.smi` downloads. The full downloads
are not required for the normal suite. To run the optional corpus survey, set
`RINGO_SMILES_CORPUS` to their directory and run:

```sh
cargo test --release --locked --test stereo corpus_coverage -- --ignored --nocapture
```
