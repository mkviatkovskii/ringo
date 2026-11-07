# SMILES support

Ringo parses the SMILES at the start of each record; trailing whitespace-separated
text is retained as the record's name or identifier. Unsupported syntax inside the
SMILES is an error rather than a partially parsed molecule.

## Atoms, bonds and components

The reader supports all element symbols in brackets, the organic subset outside
brackets, aromatic atom notation, isotopes, formal charges, bracket hydrogens,
implicit hydrogens, single/double/triple/aromatic bonds, branches, ring closures,
and disconnected components separated by a dot. The wildcard `*` is represented
as atomic number zero; it does not turn a similarity query into a SMARTS pattern.
Ring labels are one digit or `%` followed by two digits.

Supported Kekulé and aromatic spellings are normalized by aromaticity perception.
Hydrogens are counted before normalization. The molecule graph is exposed read-only
so external edits cannot leave computed hydrogen or stereo metadata stale.

## Atom maps and classes

Bracket atoms accept a trailing numeric class/map label, such as `[CH3:12]`,
`[NH4+:005]`, or `[1*:0]`. `Atom::atom_class` stores its numeric value as a
`u32`; an absent label is zero. Labels do not affect fingerprints or stereochemistry.
Malformed labels and values above `u32::MAX` are parsing errors.

## Dative bonds

`->` and `<-` specify electron-pair donation in chains, branches and ring
closures. `BondOrder::Dative` identifies the bond type and `Bond::donor` stores
its donor atom. The molecular graph remains undirected.

For implicit hydrogen calculations, a dative bond contributes one to the acceptor's
valence and zero to the donor's. Covalent aromaticity and fingerprint ring membership
exclude dative bonds. Fingerprint degree and neighbourhoods include them, with
separate codes for donor and acceptor ends. Equivalent reversed spellings, such as
`N->[Cu]` and `[Cu]<-N`, produce the same fingerprint.

## Stereochemistry

| Configuration | Supported notation |
|---|---|
| Tetrahedral | `@`, `@@`, `@TH1`, `@TH2` |
| Trigonal bipyramidal | `@TB1`–`@TB20`, supported coordination-based shorthand |
| Octahedral | `@OH1`–`@OH30`, supported coordination-based shorthand |
| Ordinary double bond | `/`, `\` |

`Molecule::atom_stereo()` and `Molecule::double_bond_stereo()` expose the
preserved configurations. Atom configurations refer to recorded SMILES ligand order,
not CIP R/S labels. TB/OH permutation numbers retain their SMILES meaning; missing
coordination sites occupy the implicit-hydrogen position, following the RDKit SMILES
convention. Supported three-neighbour N/P/S/As/Se tetrahedral configurations use a
virtual lone-pair site in that position, following the Open Babel convention.

Absence of a stereo record means unspecified stereochemistry. Invalid, incomplete,
unsupported or contradictory specifications produce parsing errors. Directional
specifications that aromaticity normalization cannot preserve are also rejected.
The reader does not perceive additional stereocentres or validate chemical
plausibility. Arithmetic and representation limits are checked independently.

ECFP and similarity search are stereo-insensitive. Preserving a configuration in
the molecule model does not yet make opposite configurations score differently.

## Known gaps

- Square-planar and allene/cumulene stereo.
- Enhanced stereo and CXSMILES interpretation. Whitespace-separated extension
  blocks are currently treated as trailing record text, not chemical metadata.
- Tilde (`~`) bonds and extended ring labels such as `%(100)`.
- Some directional/aromatic ring cases; the current parser rejects these rather
  than discarding their stereo annotations.
- CIP assignment and canonical SMILES writing.

[Issue #8](https://github.com/mkviatkovskii/ringo/issues/8) records corpus coverage
and remaining cases. Further core syntax/stereo work is tracked in
[#51](https://github.com/mkviatkovskii/ringo/issues/51), enhanced stereo/CXSMILES in
[#52](https://github.com/mkviatkovskii/ringo/issues/52), and optional chemical
validation in [#44](https://github.com/mkviatkovskii/ringo/issues/44).
