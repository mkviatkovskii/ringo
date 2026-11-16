//! Configurations relative to ordered ligand references, never CIP labels.
use petgraph::stable_graph::EdgeIndex;
use petgraph::stable_graph::NodeIndex;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Ligand {
    Atom(NodeIndex),
    Hydrogen,
    LonePair,
    /// An unoccupied site in explicitly tagged SP/TB/OH coordination.
    Vacancy,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Configuration {
    Tetrahedral(u8),
    /// SMILES SP1, SP2, SP3: opposite ligand indices are (0, 2), (0, 1), (0, 3).
    SquarePlanar(u8),
    TrigonalBipyramidal(u8),
    Octahedral(u8),
}

/// Normalized coordination arrangement invariant under the polyhedron's symmetry group.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum CoordinationGeometry<T> {
    SquarePlanar([[T; 2]; 2]),
    TrigonalBipyramidal((T, [T; 3], T)),
    Octahedral((T, [T; 4], T)),
}

#[derive(Clone, Copy)]
enum OhShape {
    U,
    Z,
    Four,
}

// OpenSMILES 3.8.6: (from, towards, anti-clockwise `@`).
const TB_SPEC: [(usize, usize, bool); 20] = [
    (0, 4, true),  // TB1
    (0, 4, false), // TB2
    (0, 3, true),  // TB3
    (0, 3, false), // TB4
    (0, 2, true),  // TB5
    (0, 2, false), // TB6
    (0, 1, true),  // TB7
    (0, 1, false), // TB8
    (1, 4, true),  // TB9
    (1, 3, true),  // TB10
    (1, 4, false), // TB11
    (1, 3, false), // TB12
    (1, 2, true),  // TB13
    (1, 2, false), // TB14
    (2, 4, true),  // TB15
    (2, 3, true),  // TB16
    (3, 4, true),  // TB17
    (3, 4, false), // TB18
    (2, 3, false), // TB19
    (2, 4, false), // TB20
];

// OpenSMILES 3.8.7: (equatorial shape, from, towards, anti-clockwise `@`).
const OH_SPEC: [(OhShape, usize, usize, bool); 30] = [
    (OhShape::U, 0, 5, true),     // OH1
    (OhShape::U, 0, 5, false),    // OH2
    (OhShape::U, 0, 4, true),     // OH3
    (OhShape::Z, 0, 5, true),     // OH4
    (OhShape::Z, 0, 4, true),     // OH5
    (OhShape::U, 0, 3, true),     // OH6
    (OhShape::Z, 0, 3, true),     // OH7
    (OhShape::Four, 0, 5, false), // OH8
    (OhShape::Four, 0, 4, false), // OH9
    (OhShape::Four, 0, 5, true),  // OH10
    (OhShape::Four, 0, 4, true),  // OH11
    (OhShape::Four, 0, 3, false), // OH12
    (OhShape::Four, 0, 3, true),  // OH13
    (OhShape::Z, 0, 5, false),    // OH14
    (OhShape::Z, 0, 4, false),    // OH15
    (OhShape::U, 0, 4, false),    // OH16
    (OhShape::Z, 0, 3, false),    // OH17
    (OhShape::U, 0, 3, false),    // OH18
    (OhShape::U, 0, 2, true),     // OH19
    (OhShape::Z, 0, 2, true),     // OH20
    (OhShape::Four, 0, 2, false), // OH21
    (OhShape::Four, 0, 2, true),  // OH22
    (OhShape::Z, 0, 2, false),    // OH23
    (OhShape::U, 0, 2, false),    // OH24
    (OhShape::U, 0, 1, true),     // OH25
    (OhShape::Z, 0, 1, true),     // OH26
    (OhShape::Four, 0, 1, false), // OH27
    (OhShape::Four, 0, 1, true),  // OH28
    (OhShape::Z, 0, 1, false),    // OH29
    (OhShape::U, 0, 1, false),    // OH30
];

// All 24 proper octahedral rotations of (top, e0, e1, e2, e3, bottom).
const OH_ROTATIONS: [[usize; 6]; 24] = [
    [0, 1, 2, 3, 4, 5],
    [0, 2, 3, 4, 1, 5],
    [0, 3, 4, 1, 2, 5],
    [0, 4, 1, 2, 3, 5],
    [1, 0, 4, 5, 2, 3],
    [1, 2, 0, 4, 5, 3],
    [1, 4, 5, 2, 0, 3],
    [1, 5, 2, 0, 4, 3],
    [2, 0, 1, 5, 3, 4],
    [2, 1, 5, 3, 0, 4],
    [2, 3, 0, 1, 5, 4],
    [2, 5, 3, 0, 1, 4],
    [3, 0, 2, 5, 4, 1],
    [3, 2, 5, 4, 0, 1],
    [3, 4, 0, 2, 5, 1],
    [3, 5, 4, 0, 2, 1],
    [4, 0, 3, 5, 1, 2],
    [4, 1, 0, 3, 5, 2],
    [4, 3, 5, 1, 0, 2],
    [4, 5, 1, 0, 3, 2],
    [5, 1, 4, 3, 2, 0],
    [5, 2, 1, 4, 3, 0],
    [5, 3, 2, 1, 4, 0],
    [5, 4, 3, 2, 1, 0],
];

fn ordered_pair<T: Ord + Copy>(a: T, b: T) -> [T; 2] {
    if a <= b {
        [a, b]
    } else {
        [b, a]
    }
}

impl Configuration {
    pub(crate) fn coordination_variant_count(self) -> u8 {
        match self {
            Configuration::Tetrahedral(_) => 0,
            Configuration::SquarePlanar(_) => 3,
            Configuration::TrigonalBipyramidal(_) => 20,
            Configuration::Octahedral(_) => 30,
        }
    }

    pub(crate) fn with_coordination_number(self, number: u8) -> Option<Configuration> {
        match self {
            Configuration::SquarePlanar(_) if (1..=3).contains(&number) => {
                Some(Configuration::SquarePlanar(number))
            }
            Configuration::TrigonalBipyramidal(_) if (1..=20).contains(&number) => {
                Some(Configuration::TrigonalBipyramidal(number))
            }
            Configuration::Octahedral(_) if (1..=30).contains(&number) => {
                Some(Configuration::Octahedral(number))
            }
            _ => None,
        }
    }

    pub(crate) fn coordination_geometry<T: Ord + Copy>(
        self,
        sites: &[T],
    ) -> Option<CoordinationGeometry<T>> {
        match self {
            Configuration::Tetrahedral(_) => None,
            Configuration::SquarePlanar(n) => {
                let s: [T; 4] = sites.try_into().ok()?;
                let (p0, p1) = match n {
                    1 => (ordered_pair(s[0], s[2]), ordered_pair(s[1], s[3])),
                    2 => (ordered_pair(s[0], s[1]), ordered_pair(s[2], s[3])),
                    3 => (ordered_pair(s[0], s[3]), ordered_pair(s[1], s[2])),
                    _ => return None,
                };
                Some(CoordinationGeometry::SquarePlanar(if p0 <= p1 {
                    [p0, p1]
                } else {
                    [p1, p0]
                }))
            }
            Configuration::TrigonalBipyramidal(n) => {
                let s: [T; 5] = sites.try_into().ok()?;
                let &(from, towards, ccw) = TB_SPEC.get(usize::from(n).checked_sub(1)?)?;
                let mut eq = [s[0]; 3];
                let mut k = 0;
                for (i, &site) in s.iter().enumerate() {
                    if i != from && i != towards {
                        eq[k] = site;
                        k += 1;
                    }
                }
                let [e0, e1, e2] = if ccw { eq } else { [eq[0], eq[2], eq[1]] };
                let (top, bot) = (s[from], s[towards]);
                let rotations = [
                    (top, [e0, e1, e2], bot),
                    (top, [e1, e2, e0], bot),
                    (top, [e2, e0, e1], bot),
                    (bot, [e0, e2, e1], top),
                    (bot, [e2, e1, e0], top),
                    (bot, [e1, e0, e2], top),
                ];
                let mut best = rotations[0];
                for candidate in &rotations[1..] {
                    if *candidate < best {
                        best = *candidate;
                    }
                }
                Some(CoordinationGeometry::TrigonalBipyramidal(best))
            }
            Configuration::Octahedral(n) => {
                let s: [T; 6] = sites.try_into().ok()?;
                let &(shape, from, towards, ccw) = OH_SPEC.get(usize::from(n).checked_sub(1)?)?;
                let mut eq = [s[0]; 4];
                let mut k = 0;
                for (i, &site) in s.iter().enumerate() {
                    if i != from && i != towards {
                        eq[k] = site;
                        k += 1;
                    }
                }
                let cyc = match shape {
                    OhShape::U => [eq[0], eq[1], eq[2], eq[3]],
                    OhShape::Z => [eq[0], eq[1], eq[3], eq[2]],
                    OhShape::Four => [eq[0], eq[3], eq[1], eq[2]],
                };
                let cyc = if ccw {
                    cyc
                } else {
                    [cyc[0], cyc[3], cyc[2], cyc[1]]
                };
                let ordered = [s[from], cyc[0], cyc[1], cyc[2], cyc[3], s[towards]];
                let tuple = |r: [usize; 6]| {
                    (
                        ordered[r[0]],
                        [ordered[r[1]], ordered[r[2]], ordered[r[3]], ordered[r[4]]],
                        ordered[r[5]],
                    )
                };
                let mut best = tuple(OH_ROTATIONS[0]);
                for &rotation in &OH_ROTATIONS[1..] {
                    let candidate = tuple(rotation);
                    if candidate < best {
                        best = candidate;
                    }
                }
                Some(CoordinationGeometry::Octahedral(best))
            }
        }
    }

    pub(crate) fn coordination_tag<T: Ord + Copy>(
        self,
        source: &[T],
        target: &[T],
    ) -> Option<String> {
        let goal = self.coordination_geometry(source)?;
        let number = (1..=self.coordination_variant_count()).find(|&number| {
            self.with_coordination_number(number)
                .and_then(|cfg| cfg.coordination_geometry(target))
                .as_ref()
                == Some(&goal)
        })?;
        Some(match self {
            Configuration::Tetrahedral(_) => return None,
            Configuration::SquarePlanar(_) => format!("@SP{number}"),
            Configuration::TrigonalBipyramidal(_) => format!("@TB{number}"),
            Configuration::Octahedral(_) => format!("@OH{number}"),
        })
    }

    pub(crate) fn all_coordination_trivial<T: Ord + Copy>(self, sites: &[T]) -> bool {
        let Some(first) = self.coordination_geometry(sites) else {
            return false;
        };
        (1..=self.coordination_variant_count()).all(|number| {
            self.with_coordination_number(number)
                .and_then(|cfg| cfg.coordination_geometry(sites))
                .as_ref()
                == Some(&first)
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomStereo {
    pub center: NodeIndex,
    /// SMILES encounter order, including hydrogen and ring-closure positions.
    pub ligands: Vec<Ligand>,
    pub configuration: Configuration,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlleneStereo {
    pub center: NodeIndex,
    /// Terminal-to-terminal double-bond chain, ordered by the center's SMILES neighbours.
    /// The center is its middle atom.
    pub axis: Vec<NodeIndex>,
    /// Two SMILES-ordered substituents at axis[0], then two at its last atom.
    pub ligands: [Ligand; 4],
    /// AL1/AL2 (or @/@@), interpreted by collapsing the axis to a tetrahedron.
    pub permutation: u8,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DoubleBondStereo {
    pub bond: EdgeIndex,
    pub ends: [NodeIndex; 2],
    pub references: [NodeIndex; 2],
    pub same_side: bool,
}
