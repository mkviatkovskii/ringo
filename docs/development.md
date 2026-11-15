# Development

Use stable Rust. The normal checks run without downloading chemical databases:

```sh
cargo +stable fmt --all --check
cargo +stable clippy --all-targets --locked -- -D warnings
cargo +stable test --locked
```

CI runs tests on Linux and Windows. For a faster local build, `cargo test` uses
the debug profile; use `--release` when measuring indexing or search performance.

## Corpus fixtures and coverage

Small SMILES fixtures under `tests/data` retain their ChEBI/ChEMBL identifiers.
The full downloaded databases are not committed and are not required by the normal
test suite.

The optional survey currently expects both `chebi_20261105.smi` and `chembl_37.smi`
in the directory named by `RINGO_SMILES_CORPUS`. Set the variable to your download
directory, then run the ignored survey test.

PowerShell:

```powershell
$env:RINGO_SMILES_CORPUS = 'C:\path\to\smiles'
cargo +stable test --release --locked --test stereo corpus_coverage -- --ignored --nocapture
```

POSIX shell:

```sh
RINGO_SMILES_CORPUS=/path/to/smiles cargo +stable test --release --locked --test stereo corpus_coverage -- --ignored --nocapture
```

The survey loads every record once and reports accepted/rejected counts for plain
and stereo-marked records, plus the count and examples of loaded molecules that
fail optional post-load validation, by kind: valence or charge, and each aromaticity
check. It checks that file size and modification time remain
unchanged during each scan. A passing survey is not a claim of complete parser
coverage or chemical validity. The checksum-pinned baseline and complete rejection
classification are recorded in [issue #8](https://github.com/mkviatkovskii/ringo/issues/8).

A valence failure can reflect an unsupported
representation, such as covalent metal coordination or neutral hypervalent nitrogen.

## Index compatibility during development

Index format remains version 1 for the current development release; version changes
are made once per release. Rebuild development indexes when parser changes admit
previously skipped records or fingerprint semantics change. Atom class/map labels
and the addition of dative-bond support leave existing covalent fingerprints unchanged.
Square-planar and allene/cumulene tags also leave fingerprints unchanged; rebuild
indexes to include previously rejected SP/AL-tagged records. Aromaticity perception
no longer counts a neutral O/S bridge in a ring of nine or more atoms as a π donor.
Rebuild indexes: Kekulé-written macrocycles with such bridges change fingerprints,
and records such as CHEMBL2004566 are now accepted.
