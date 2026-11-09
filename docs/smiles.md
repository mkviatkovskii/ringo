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
not CIP R/S labels. SP/TB/OH permutation numbers follow
[OpenSMILES stereochemistry](https://opensmiles.org/opensmiles.html#_stereochemistry).
Ringo additionally supports missing coordination sites at the implicit-hydrogen
position, and virtual lone-pair sites at that position for three-neighbour
N/P/S/As/Se tetrahedral configurations. These are extensions beyond the standard's
explicit stereo rules.

Square-planar tags preserve the four SMILES-ordered ligand sites and the
[OpenSMILES U, 4 and Z permutations](https://opensmiles.org/opensmiles.html#_square_planar_centers).
Using zero-based ligand indices, the opposite pairs are `(0, 2)` and `(1, 3)` for
`SP1`, `(0, 1)` and `(2, 3)` for `SP2`, and `(0, 3)` and `(1, 2)` for `SP3`.
Repeated ligands and dative bonds are supported: `Cl[Pt@SP1](Cl)(<-[NH3])<-[NH3]`
describes cisplatin, while `Cl[Pt@SP2](Cl)(<-[NH3])<-[NH3]` describes transplatin.
Bracket hydrogens and missing SP sites use the same insertion position as TB/OH;
missing sites are `Ligand::Vacancy`, not lone pairs. More than four SP sites is an
error.

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

## Optional molecule validation

Chemical validation is a separate optional step on an already loaded molecule:

```rust
use ringo::io::smiles::reader::molecule::parse_molecule;

let (_, molecule) = parse_molecule("C(C)(C)(C)(C)C").unwrap();
assert!(molecule.validate().is_err());
```

`Molecule::validate()` returns `Result<(), String>`. This first implementation checks
valence; aromaticity consistency and chemical stereogenicity checks are not yet
implemented. A successful result only means the implemented checks passed.
Errors identify the
zero-based atom index and its excessive valence or unsupported charge state.
Validation never modifies the molecule, its stereo records or its fingerprints.
Loading always checks syntax, stereo representation and `u8` arithmetic, including
bond-order and hydrogen totals on bracket atoms, independently of this optional step.

The first implemented check covers upper valence limits for H, B, C, N, O, F, Si, P, S,
Cl, As, Se, Br and I. Neutral ceilings are 1 for H/F/Cl/Br, 2 for O, 3 for B/N,
4 for C/Si, 5 for P/As/I and 6 for S/Se. Lower valences are accepted, including
bracket atoms with incomplete valence; this is not radical perception.

For these main-group atoms, charged forms use the ceiling of the isoelectronic
element in the same periodic-table row when the table supplies one. Thus the
common perchlorate spelling `[O-][Cl+3]([O-])([O-])[O-]` is supported. P/As with charge -2 and S/Se
with charge -1 instead have a ceiling of 5. Unsupported charge/table combinations
fail validation. Neutral H allows one bond; H+ and H- allow none. Other hydrogen
charges fail validation. These chemical validation limits are Ringo's documented
policy; they are distinct from the normal valences used to infer hydrogens in the
[OpenSMILES organic-subset rules](https://opensmiles.org/opensmiles.html#_organic_subset).
For example, Ringo's optional validation rejects neutral nitrogen above valence
three, although the notation's implicit-hydrogen rules include valence five.
Halide anions F-/Cl-/Br-/I- all have a ceiling of zero; polyhalides such as
`Br[Br-]Br` and `I[I-]I` load permissively but fail this bounded covalent-valence policy.

Loading records valence counts before aromaticity normalization replaces written
bond orders, so post-load validation can still detect an over-valent input atom.
Dative bonds contribute zero at the donor and one at the acceptor.
Aromatic bonds contribute one each, plus one per atom needing a ring pi bond;
lone-pair donors, empty p orbitals and atoms with a written double/triple bond do
not receive that extra contribution. This admits pyridine, pyrrole, pyrrolide,
fused aromatic carbon and supported exocyclic carbonyls without counting every
aromatic bond as a double bond.

This is aromatic-aware valence bookkeeping, not aromaticity validation. The
[OpenSMILES aromaticity rules](https://opensmiles.org/opensmiles.html#_aromaticity)
also require a consistent electron assignment for the marked aromatic system,
accounting for hydrogens, external bonds and charges. Those system-wide checks
remain a separate step in #44.

Elements outside this list, including transition metals and wildcard atoms, have
no chemical valence limit in this policy. Aromatic kekulizability, isotope
stability, ligand distinctness, stereogenicity and general chemical feasibility
are not validated. Neutral hypervalent nitrogen is rejected rather than rewritten
into a charged representation. Passing validation means only passing these rules.

Parsing never invokes chemical validation. Library callers can call the method
explicitly after loading. The existing `index_file` and `similarity_search` APIs
remain permissive; `index_file_with_options` and `similarity_search_with_options`
accept `IndexOptions` and `SearchOptions` with a `validate` field to request that
separate post-load step before fingerprinting.

Append `--validate` as the last argument to either CLI to select the same step.
Indexing counts validation failures as skipped records and stores empty fingerprints,
just as for loading failures. Search reports a validation error for an invalid query.
The search switch checks only the loaded query, not records already in the index.
Rebuild an index with `--validate` to exclude failed dataset records. The chosen
validation policy is not stored in the index; accepted fingerprints and index
format version 1 remain unchanged.

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
