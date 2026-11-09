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

## Aromaticity perception

Loading keeps the written aromatic flags and bond orders for optional validation,
then normalizes them. A ring written in lowercase or with `:` bonds stays aromatic
as written, including anti-aromatic rings such as `c1ccc1`. A Kekulé ring becomes
aromatic when its atoms hold 4n+2 π electrons by the
[OpenSMILES electron counts](https://opensmiles.org/opensmiles.html#_aromaticity),
either alone or as a fused pair such as azulene. An aromatic bond outside every
aromatic ring becomes single.

A neutral O or S between two single bonds gives no π electrons when the smallest
ring through it has nine or more atoms: in a macrocycle it is a bridge, not a donor
as in furan. RDKit adopted the same rule in 2026.09
([rdkit/rdkit#9398](https://github.com/rdkit/rdkit/issues/9398)). For example, the
10-atom lactone ring of CHEMBL2004566,
`O=C1Nc2ccccc2C(=O)O/C(Cc2ccccc2)=N\c2ccccc21`, is not aromatic. It keeps its C=N
double bond and that bond's configuration. Its three benzene rings are aromatic.
The rule depends on constitution, not on stereo marks, so the stereo-free spelling
gives the same atoms, bonds and fingerprint. Smaller rings such as furan and 2-pyrone, and nitrogen donors in
macrocycles such as azonine, are unaffected.

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

`Molecule::validate()` returns `Result<(), String>`. It checks valence, then the
written aromatic annotations; chemical stereogenicity checks are not yet
implemented. A successful result only means the implemented checks passed.
Errors identify the zero-based atom index and its excessive valence, unsupported
charge state or aromatic inconsistency.
Loading always checks syntax, stereo representation and arithmetic limits for
implicit-hydrogen calculation, independently of this optional step.

The first implemented check covers upper valence limits for H, B, C, N, O, F, Si, P, S,
Cl, As, Se, Br and I. Neutral ceilings are 1 for H/F/Cl/Br, 2 for O, 3 for B/N,
4 for C/Si, 5 for P/As/I and 6 for S/Se. Lower valences are accepted, including
bracket atoms with incomplete valence; this is not radical perception.

For any known element, a positive charge greater than its atomic number is rejected:
it would require removing more electrons than the neutral atom has. Isolated
cations with no bonds or attached hydrogens pass the valence check when their
charge is within that bound, including `[Si+4]`, `[As+3]` and `[Se+4]`; they do not
need a same-period valence-table entry. This does not check whether an ion is stable.

Bonded main-group atoms use the ceiling of the isoelectronic element in the same
periodic-table row when the table supplies one. The table allows one bond for the
three-electron case, so `[H][B+2]` is supported but `[B+2](F)F` is rejected. The
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

Loading keeps the written bond orders and aromatic flags before aromaticity
normalization replaces them, so post-load validation can still detect an over-valent
input atom. Dative bonds contribute zero at the donor and one at the acceptor.
Aromatic bonds contribute one each, plus one per aromatic atom needing a ring pi
bond (see below). This admits pyridine, pyrrole, pyrrolide, fused aromatic carbon and
supported exocyclic carbonyls without counting every aromatic bond as a double bond.

### Aromatic annotations

The [OpenSMILES aromaticity rules](https://opensmiles.org/opensmiles.html#_aromaticity)
require that electrons can be assigned to the atoms written aromatic, consistent with
their hydrogens, external bonds and charges. Validation checks the annotations as
written, before normalization. An aromatic atom is one written in lowercase or with a
`:` bond. Validation fails when:

- an aromatic atom is in no ring, as in `CccccC`;
- a wildcard `*` is in the written aromatic system, which Ringo does not support:
  joined to an aromatic atom by a `:` bond, or on a ring whose other atoms are all
  aromatic, as in `c1cc*cc1`. A wildcard in a fused non-aromatic ring, as in
  `c1ccc2c(c1)C*2`, is allowed;
- the aromatic atoms needing a ring double bond cannot pair up along aromatic ring
  bonds, so there is no Kekulé structure. Examples are `c1cccc1`, pyrrole written
  `c1ccnc1` without its NH, and `C=c1ccccc1`.

An aromatic atom needs a ring double bond unless it has a written double or triple
bond (2-pyridone's C=O carbon), or gives a lone pair or an empty p orbital. Neutral
O/S/Se/Te and N- give a lone pair. Neutral B/N/P/As and charged carbon with three
neighbours or hydrogens give a lone pair or an empty p orbital: pyrrole's NH,
indolizine's bridgehead N, borole's B, `[cH-]` and `[cH+]`. With two neighbours they
need a double bond, like pyridine's N or borabenzene's B. Charged N/O/S, such as
pyridinium and pyrylium, need one.

An aromatic bond outside rings, such as the link in `c1ccccc1c1ccccc1`, is single,
as loading reads it. It cannot supply a ring double bond, so `c1cccc1c1cccc1` fails.
Hückel's rule is not required of written aromatic rings: OpenSMILES also writes
anti-aromatic rings such as `c1ccc1` in lowercase. Kekulé input has no aromatic
annotations to check, even where perception later makes its rings aromatic.
Radicals are not perceived: an atom with too few written hydrogens can still pass if
a Kekulé structure exists.

Elements outside the valence list, including transition metals and wildcard atoms,
have no chemical valence limit in this policy. Isotope stability, ligand
distinctness, stereogenicity and general chemical feasibility are not validated.
Neutral hypervalent nitrogen is rejected rather than rewritten into a charged
representation.

## Known gaps

- Enhanced stereo and CXSMILES interpretation. Whitespace-separated extension
  blocks are currently treated as trailing record text, not chemical metadata.
- Tilde (`~`) bonds and extended ring labels such as `%(100)`.
- Some directional/aromatic ring cases, such as the expanded porphyrin CHEBI:47035;
  the current parser rejects these rather than discarding their stereo annotations.
- CIP assignment and canonical SMILES writing.

[Issue #8](https://github.com/mkviatkovskii/ringo/issues/8) records corpus coverage
and remaining cases. Further core syntax/stereo work is tracked in
[#51](https://github.com/mkviatkovskii/ringo/issues/51) (completed), enhanced stereo/CXSMILES in
[#52](https://github.com/mkviatkovskii/ringo/issues/52), and optional
validation in [#44](https://github.com/mkviatkovskii/ringo/issues/44).
