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

Indexing creates `tests/data/molecules.smi.fp`. Search prints a similarity score
followed by the original input record, separated by a tab. Scores range from 0 to 1;
1 means identical fingerprints. Reindex after editing the SMILES file.

## Input format

Use a UTF-8 text file with one SMILES record per line. An optional name or identifier
may follow the SMILES, separated by whitespace:

```text
CCO ethanol
CC(=O)O acetic_acid
```

The indexer reports how many records could not be parsed. Those records are omitted
from search results. An unsupported or malformed query produces an error.

## Command-line options

```text
ringo-index <smiles_file> [radius=2] [bits=512]
ringo-search <smiles_file> <query> <min_similarity> <limit> [max_similarity=1.0] [best|first]
```

The index stores its fingerprint radius and bit count; queries use those same
parameters. Radius 2 corresponds to ECFP4.

Search includes scores between the minimum and maximum, inclusive. The default
`best` mode returns the most similar records first, up to the limit. `first`
returns the first matching records in file order and stops once it reaches the limit.

## SMILES support and limits

The reader handles aromatic and Kekulé notation, disconnected components, isotopes,
atom map/class labels, dative bonds, and supported atom and double-bond stereo tags.
See [SMILES support](docs/smiles.md) for the exact syntax, model representation and
unsupported cases.

Fingerprints and similarity search are currently **stereo-insensitive**: opposite
configurations can have identical fingerprints. The reader preserves supported
stereochemistry but does not validate chemical plausibility or assign CIP labels.
Ringo does not yet read every SMILES extension.

## Development

See [development and corpus checks](docs/development.md) for test commands and the
optional ChEBI/ChEMBL coverage survey. Generate the Rust API documentation with
`cargo +stable doc --no-deps`.

## License

Apache-2.0.
