# Ringo

Ringo is a file-based chemical similarity search engine written in Rust. It builds
an ECFP fingerprint index for a SMILES file and searches it using Tanimoto similarity.
It provides two command-line tools and a Rust library.

## Getting started

Use the current stable Rust toolchain. From the repository root, build both tools:

```sh
cargo +stable build --release --locked
```

On x86-64, the repository targets x86-64-v3. For older CPUs, set the `RUSTFLAGS`
environment variable to `-C target-cpu=x86-64` before building or running Cargo.

The repository includes a small SMILES file containing ibuprofen and naproxen.
Index it, then search for compounds similar to ibuprofen:

```sh
cargo +stable run --release --locked --bin ringo-index -- tests/data/molecules.smi
cargo +stable run --release --locked --bin ringo-search -- tests/data/molecules.smi "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O" 0.4 10
```

Indexing creates `tests/data/molecules.smi.fp` and `tests/data/molecules.smi.key`.
Search prints a similarity score followed by the original input record, separated by
a tab. Scores range from 0 to 1; 1 means identical fingerprints. Reindex after editing
the SMILES file.

To find the records that hold one molecule, however they were written:

```sh
cargo +stable run --release --locked --bin ringo-search -- --exact tests/data/molecules.smi "OC(=O)C(C)C1=CC=C(CC(C)C)C=C1"
```

## Input format

Use a UTF-8 text file with one SMILES record per line. An optional name or identifier
may follow the SMILES, separated by whitespace:

```text
CCO ethanol
CC(=O)O acetic_acid
```

The indexer reports how many records could not be parsed or failed enabled validation. Those records are omitted
from search results. An unsupported or malformed query produces an error.

## Command-line options

```text
ringo-index <smiles_file> [radius=2] [bits=512] [batch_lines=65536] [--validate]
ringo-search <smiles_file> <query> <min_similarity> <limit> [max_similarity=1.0] [best|first] [--validate]
ringo-search --exact <smiles_file> <query>
```

The index stores its fingerprint radius and bit count; queries use those same
parameters. Radius 2 corresponds to ECFP4.

Both tools provide `--help` and `--version`; flags may appear before or after
positional arguments.

Add `--validate` to run optional post-load checks before fingerprinting; see the
[validation policy and limits](docs/smiles.md#optional-molecule-validation).
Search checks only the query; build the index with `--validate` to leave out records
that fail.

Search includes scores between the minimum and maximum, inclusive. The default
`best` mode returns the most similar records first, up to the limit; equal scores
keep file order. `first` returns the first matching records in index order and stops
once it reaches the limit.

The index groups records by the number of fingerprint bits set, in file order within
a group, and search reads only the groups that can reach the minimum score. Indexing needs free disk space for about twice the
size of the index while it runs.

Exact search prints the records that hold the same molecule as the query, in file
order. Atom order, ring-closure digits, Kekulé or aromatic form and hydrogens written
as atoms don't matter. Element, isotope, charge, hydrogen count, bonds and
tetrahedral, double-bond and allene stereo must agree. A stereo mark that a symmetry
of the molecule undoes, as in `C[C@H](C)O`, counts as absent. A record that can't be
compared, because it has square-planar, trigonal-bipyramidal or octahedral stereo or
is too symmetric to check, is reported on stderr. The `.key` file holds a structure
key per record, sorted, so a lookup reads only the records with the query's key.

## SMILES support and limits

The reader handles aromatic and Kekulé notation, disconnected components, isotopes,
atom map/class labels, dative bonds, and supported atom and double-bond stereo tags.
See [SMILES support](docs/smiles.md) for the exact syntax, model representation and
unsupported cases.

Fingerprints and similarity search are currently **stereo-insensitive**: opposite
configurations can have identical fingerprints. The reader preserves supported
stereochemistry. A separate optional valence check does not establish chemical
plausibility or assign CIP labels.
Ringo does not yet read every SMILES extension.

## Development

See [development and corpus checks](docs/development.md) for test commands and the
optional ChEBI/ChEMBL coverage survey. Generate the Rust API documentation with
`cargo +stable doc --no-deps`.

## License

Apache-2.0.
