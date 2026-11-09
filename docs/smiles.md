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
| Allene/cumulene | `@AL1`, `@AL2`, `@`, `@@` on the middle atom of an allene-like axis |
| Square planar | `@SP1`–`@SP3`; degree-four shorthand remains tetrahedral |
| Trigonal bipyramidal | `@TB1`–`@TB20`, supported coordination-based shorthand |
| Octahedral | `@OH1`–`@OH30`, supported coordination-based shorthand |
| Ordinary double bond | `/`, `\` |

`Molecule::atom_stereo()`, `Molecule::allene_stereo()` and `Molecule::double_bond_stereo()` expose the
preserved configurations. Atom configurations refer to recorded SMILES ligand order,
not CIP R/S labels. SP/TB/OH permutation numbers retain their SMILES meaning; missing
coordination sites occupy the implicit-hydrogen position, following the RDKit SMILES
convention. Supported three-neighbour N/P/S/As/Se tetrahedral configurations use a
virtual lone-pair site in that position, following the Open Babel convention.

Square-planar tags preserve the four SMILES-ordered ligand sites and the
[OpenSMILES U, 4 and Z permutations](https://opensmiles.org/opensmiles.html#_square_planar_centers).
Using zero-based ligand indices, the opposite pairs are `(0, 2)` and `(1, 3)` for
`SP1`, `(0, 1)` and `(2, 3)` for `SP2`, and `(0, 3)` and `(1, 2)` for `SP3`.
Repeated ligands and dative bonds are supported: `Cl[Pt@SP1](Cl)(<-[NH3])<-[NH3]`
describes cisplatin, while `Cl[Pt@SP2](Cl)(<-[NH3])<-[NH3]` describes transplatin.
Bracket hydrogens and missing SP sites use the same insertion position as TB/OH;
missing sites are `Ligand::Vacancy`, not lone pairs. More than four SP sites is an
error. See the [RDKit conventions](https://www.rdkit.org/docs/RDKit_Book.html#support-for-non-tetrahedral-atomic-stereochemistry)
for hydrogen and missing-site examples.

Allene-like axes contain an even number of consecutive double bonds (three, five,
seven, ... atoms), with the tag on the middle atom. `AlleneStereo` stores that
center, the complete ordered terminal-to-terminal axis, and four terminal ligand
references: two at the first end followed by two at the other end. The first axis
arm is the center's first neighbour in SMILES encounter order. Each pair follows
SMILES encounter order, including ring-closure slots; suppressed terminal hydrogens
immediately follow their terminal atom. Explicit isotopic hydrogens remain atom
references. AL1/AL2 and their `@`/`@@` aliases use the
[OpenSMILES collapsed-tetrahedron convention](https://opensmiles.org/opensmiles.html#_tetrahedral_allene_like_systems).
For example, `NC(Br)=[C@]=C(O)C` records terminal references N, Br, O, C.
The [Daylight hydrogen example](https://daylight.com/dayhtml/doc/theory/theory.smiles.html)
`OC=[C@]=CF` is equivalent to `OC([H])=[C@AL1]=C([H])F`.

Every internal axis atom must have exactly two double-bond neighbours and no H;
each terminal must supply exactly two single-bonded or hydrogen ligand sites.
An off-center tag, cyclic axis, odd number of double bonds, extra axis branches,
missing terminal sites or other terminal bond types is an error. Ligands may
themselves connect through rings. Ordinary double-bond directional records are
separate; `/` and `\` do not encode an allene axis here. These checks preserve the
representation, without requiring chemically distinct terminal substituents.

Absence of a stereo record means unspecified stereochemistry. Invalid, incomplete,
unsupported or contradictory specifications produce parsing errors. Directional
specifications that aromaticity normalization cannot preserve are also rejected.
The reader does not perceive additional stereocentres. Optional valence validation
is described below; arithmetic and representation limits are checked independently.

ECFP and similarity search are stereo-insensitive. Preserving a configuration in
the molecule model does not yet make opposite configurations score differently.

## Optional valence validation

Existing library functions and CLI commands remain permissive by default.
`parse_molecule_with_validation(input, true)` enables valence checks; passing `false`
is equivalent to `parse_molecule(input)`. For example, `C(C)(C)(C)(C)C` parses in
permissive mode and fails with validation. Neither mode changes the structure to
repair its chemistry. Syntax, stereo representation and checked `u8` arithmetic
are mandatory in both modes, including bond-order and hydrogen totals on bracket atoms.

The initial policy checks upper valence limits for H, B, C, N, O, F, Si, P, S,
Cl, As, Se, Br and I. Neutral ceilings are 1 for H/F/Cl/Br, 2 for O, 3 for B/N,
4 for C/Si, 5 for P/As/I and 6 for S/Se. Lower valences are accepted, including
bracket atoms with incomplete valence; this is not radical perception.

For these main-group atoms, charged forms use the ceiling of the isoelectronic
element in the same periodic-table row when the table supplies one. Thus the
common perchlorate spelling `[O-][Cl+3]([O-])([O-])[O-]` is supported. P/As with charge -2 and S/Se
with charge -1 instead have a ceiling of 5. Unsupported charge/table combinations
fail validation. Neutral H allows one bond; H+ and H- allow none. Other hydrogen
charges fail validation. The policy follows the
[RDKit valence ceilings and charge rules](https://www.rdkit.org/docs/RDKit_Book.html#valence-calculation-and-allowed-valences),
without its structure cleanup or full sanitization pipeline.

Bond orders and explicit/implicit hydrogens are counted before aromaticity
normalization. Dative bonds contribute zero at the donor and one at the acceptor.
Aromatic bonds contribute one each, plus one per atom needing a ring pi bond;
lone-pair donors, empty p orbitals and atoms with a written double/triple bond do
not receive that extra contribution. This admits pyridine, pyrrole, pyrrolide,
fused aromatic carbon and supported exocyclic carbonyls without counting every
aromatic bond as a double bond.

Elements outside this list, including transition metals and wildcard atoms, have
no chemical valence limit in this policy. Aromatic kekulizability, isotope
stability, ligand distinctness, stereogenicity and general chemical feasibility
are not validated. Neutral hypervalent nitrogen is rejected rather than rewritten
into a charged representation. Passing validation means only passing these rules.

`index_file_with_validation(..., true)` counts rejected records as skipped and
stores empty fingerprints, just as for syntax errors.
`similarity_search_with_validation(..., true)` rejects an invalid query with
`InvalidInput`. Their existing counterparts select permissive parsing. Append
`--validate` as the last argument to either CLI to select these checks.

The search option checks only the query. Build the index with validation to
exclude invalid dataset records. Validation does not change accepted fingerprints,
and the index does not store the chosen policy; rebuild it when changing policies.
Index format remains version 1.

## Known gaps

- Enhanced stereo and CXSMILES interpretation. Whitespace-separated extension
  blocks are currently treated as trailing record text, not chemical metadata.
- Tilde (`~`) bonds and extended ring labels such as `%(100)`.
- Some directional/aromatic ring cases; the current parser rejects these rather
  than discarding their stereo annotations.
- CIP assignment and canonical SMILES writing.

[Issue #8](https://github.com/mkviatkovskii/ringo/issues/8) records corpus coverage
and remaining cases. Further core syntax/stereo work is tracked in
[#51](https://github.com/mkviatkovskii/ringo/issues/51) (completed), enhanced stereo/CXSMILES in
[#52](https://github.com/mkviatkovskii/ringo/issues/52), and the initial optional
valence policy in [#44](https://github.com/mkviatkovskii/ringo/issues/44).
