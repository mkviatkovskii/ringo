//! Whether two molecules are the same, however they were written (#58).
use crate::model::bond::{Bond, BondOrder};
use crate::model::molecule::Molecule;
use crate::model::stereo::{Configuration, Ligand};
use petgraph::algo::subgraph_isomorphisms_iter;
use petgraph::graph::{Graph, NodeIndex as GraphIndex};
use petgraph::stable_graph::{EdgeReference, NodeIndex};
use petgraph::visit::{EdgeRef, IntoEdgeReferences, NodeIndexable};
use petgraph::Undirected;
use std::cell::OnceCell;
use std::collections::HashMap;

// shortcut: large symmetric molecules give up with TooSymmetric (49 of ChEMBL's 1.09M
// stereo molecules, all dendrimers and multivalent conjugates, against their stereo-free
// form; up to 1.2 s each). It caps the mappings tried, not VF2's own backtracking. A
// stereo-aware canonical form (#63) would remove the limit.
const MAX_MAPPINGS: usize = 10_000;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExactMatchError {
    /// Square-planar, trigonal-bipyramidal or octahedral stereo, which isn't compared yet.
    UnsupportedStereo,
    /// More than `MAX_MAPPINGS` atom mappings to try.
    TooSymmetric,
}

#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
struct AtomKey {
    atomic_number: u8,
    isotope: u16,
    charge: i8,
    hydrogens: usize,
    aromatic: bool,
    /// Folded twin end atoms, with their bonds, sorted.
    twins: Vec<(AtomKey, BondOrder)>,
}

/// What became of a source atom.
#[derive(Clone, Copy)]
enum Fate {
    Kept(usize),
    /// Folded into its neighbour's hydrogen count.
    Hydrogen,
    /// One of two or more alike end atoms on one neighbour, folded into its key.
    Twin,
}

/// A stereo ligand, by its atom in the folded graph.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Site {
    Atom(usize),
    Hydrogen,
    LonePair,
    Twin,
}

impl Site {
    fn mapped(self, mapping: &[usize]) -> Site {
        match self {
            Site::Atom(atom) => Site::Atom(mapping[atom]),
            other => other,
        }
    }
}

enum Stereo {
    /// A tetrahedral centre, or an allene's middle atom with its four terminal ligands.
    Parity {
        center: usize,
        ligands: [Site; 4],
        clockwise: bool,
    },
    DoubleBond {
        ends: [usize; 2],
        references: [Site; 2],
        same_side: bool,
    },
}

impl Stereo {
    fn mapped(&self, mapping: &[usize]) -> Stereo {
        match *self {
            Stereo::Parity {
                center,
                ligands,
                clockwise,
            } => Stereo::Parity {
                center: mapping[center],
                ligands: ligands.map(|site| site.mapped(mapping)),
                clockwise,
            },
            Stereo::DoubleBond {
                ends,
                references,
                same_side,
            } => Stereo::DoubleBond {
                ends: ends.map(|end| mapping[end]),
                references: references.map(|site| site.mapped(mapping)),
                same_side,
            },
        }
    }

    /// The atoms it belongs to, sorted, so a mapped element can be looked up.
    fn anchor(&self, mapping: impl Fn(usize) -> usize) -> [usize; 2] {
        match self {
            Stereo::Parity { center, .. } => [mapping(*center); 2],
            Stereo::DoubleBond { ends, .. } => {
                let [a, b] = ends.map(mapping);
                [a.min(b), a.max(b)]
            }
        }
    }

    /// Whether this element, carried over by `mapping`, has `target`'s configuration.
    fn agrees(&self, mapping: &[usize], target: &Stereo) -> bool {
        match (self, target) {
            (
                Stereo::Parity {
                    ligands, clockwise, ..
                },
                Stereo::Parity {
                    ligands: target_ligands,
                    clockwise: target_clockwise,
                    ..
                },
            ) => {
                let mut positions = [0; 4];
                for (position, site) in positions.iter_mut().zip(ligands) {
                    let site = site.mapped(mapping);
                    let Some(found) = target_ligands.iter().position(|&other| other == site) else {
                        return false;
                    };
                    *position = found;
                }
                let inversions = (0..4)
                    .flat_map(|i| (i + 1..4).map(move |j| (i, j)))
                    .filter(|&(i, j)| positions[i] > positions[j])
                    .count();
                (inversions % 2 == 0) == (clockwise == target_clockwise)
            }
            (
                Stereo::DoubleBond {
                    ends,
                    references,
                    same_side,
                },
                Stereo::DoubleBond {
                    ends: target_ends,
                    references: target_references,
                    same_side: target_same_side,
                },
            ) => {
                let target_references = if *target_ends == ends.map(|end| mapping[end]) {
                    *target_references
                } else {
                    [target_references[1], target_references[0]]
                };
                // each end's other substituent is on the opposite side
                let flips = (0..2)
                    .filter(|&i| references[i].mapped(mapping) != target_references[i])
                    .count();
                (*same_side == (flips % 2 == 0)) == *target_same_side
            }
            _ => false,
        }
    }
}

/// A molecule as compared. Hydrogens written as atoms are folded into their neighbour's
/// count, so `[H]C([H])([H])[H]` is the same as `C`. Alike end atoms on one neighbour,
/// such as the methyls of `N(C)C` or the oxygens of `S(=O)(=O)`, are folded into its key:
/// swapping them changes nothing, but VF2 would try both orders of every such pair.
struct Folded {
    graph: FoldedGraph,
    /// (donor, acceptor) of each dative bond.
    datives: Vec<(usize, usize)>,
    stereo: Vec<Stereo>,
    /// Index into `stereo` by anchor.
    anchors: HashMap<[usize; 2], usize>,
    fake: OnceCell<Result<Vec<bool>, ExactMatchError>>,
}

type FoldedGraph = Graph<AtomKey, BondOrder, Undirected>;

impl Folded {
    fn new(graph: FoldedGraph, datives: Vec<(usize, usize)>, stereo: Vec<Stereo>) -> Folded {
        let anchors = stereo
            .iter()
            .enumerate()
            .map(|(i, element)| (element.anchor(|atom| atom), i))
            .collect();
        Folded {
            graph,
            datives,
            stereo,
            anchors,
            fake: OnceCell::new(),
        }
    }

    /// One `Folded` per connected component, so k alike salt or water parts cost k
    /// matches instead of k! mappings. No stereo element spans two components.
    fn components(molecule: &Molecule) -> Result<Vec<Folded>, ExactMatchError> {
        let (graph, fates, datives) = Folded::fold(molecule);
        let stereo = Folded::read_stereo(molecule, &graph, &fates)?;
        let mut component = vec![usize::MAX; graph.node_count()];
        // each atom's index within its component
        let mut local = vec![0; graph.node_count()];
        let mut parts: Vec<FoldedGraph> = Vec::new();
        for start in graph.node_indices() {
            if component[start.index()] != usize::MAX {
                continue;
            }
            let mut part = Graph::new_undirected();
            let mut stack = vec![start];
            component[start.index()] = parts.len();
            while let Some(node) = stack.pop() {
                local[node.index()] = part.add_node(graph[node].clone()).index();
                for next in graph.neighbors(node) {
                    if component[next.index()] == usize::MAX {
                        component[next.index()] = parts.len();
                        stack.push(next);
                    }
                }
            }
            parts.push(part);
        }
        for edge in graph.edge_references() {
            let (a, b) = (edge.source().index(), edge.target().index());
            parts[component[a]].add_edge(
                GraphIndex::new(local[a]),
                GraphIndex::new(local[b]),
                *edge.weight(),
            );
        }
        let mut part_datives = vec![Vec::new(); parts.len()];
        for (donor, acceptor) in datives {
            part_datives[component[donor]].push((local[donor], local[acceptor]));
        }
        let mut part_stereo: Vec<Vec<Stereo>> = parts.iter().map(|_| Vec::new()).collect();
        for element in stereo {
            let [atom, _] = element.anchor(|atom| atom);
            part_stereo[component[atom]].push(element.mapped(&local));
        }
        Ok(parts
            .into_iter()
            .zip(part_datives)
            .zip(part_stereo)
            .map(|((graph, datives), stereo)| Folded::new(graph, datives, stereo))
            .collect())
    }

    /// Whether `target` is the same component.
    fn matches(&self, target: &Folded) -> Result<bool, ExactMatchError> {
        any_mapping(self, target, |mapping| {
            let datives = self.datives.iter().all(|&(donor, acceptor)| {
                target
                    .datives
                    .contains(&(mapping[donor], mapping[acceptor]))
            });
            Ok(datives && self.stereo_agrees(target, mapping)?)
        })
    }

    /// The folded graph, what became of each source atom, and the dative bonds.
    fn fold(molecule: &Molecule) -> (FoldedGraph, Vec<Fate>, Vec<(usize, usize)>) {
        let source = molecule.graph();
        let other = |edge: EdgeReference<'_, Bond>, node: NodeIndex| {
            if edge.source() == node {
                edge.target()
            } else {
                edge.source()
            }
        };
        // a plain hydrogen on one heavy atom; [2H], [H+] and H2 stay atoms
        let hydrogen = |node: NodeIndex| {
            let atom = &source[node];
            let mut edges = source.edges(node);
            atom.element.atomic_number == 1
                && atom.isotope == 0
                && atom.charge == 0
                && atom.hs + atom.implicit_hs == 0
                && matches!((edges.next(), edges.next()), (Some(edge), None)
                    if edge.weight().order == BondOrder::Single
                        && source[other(edge, node)].element.atomic_number != 1)
        };
        let mut keys: Vec<Option<AtomKey>> = vec![None; source.node_bound()];
        for node in source.node_indices().filter(|&node| !hydrogen(node)) {
            let atom = &source[node];
            let folded = source
                .neighbors(node)
                .filter(|&other| hydrogen(other))
                .count();
            keys[node.index()] = Some(AtomKey {
                atomic_number: atom.element.atomic_number,
                isotope: atom.isotope,
                charge: atom.charge,
                hydrogens: usize::from(atom.hs) + usize::from(atom.implicit_hs) + folded,
                aromatic: atom.aromatic,
                twins: Vec::new(),
            });
        }

        // end atoms: one covalent bond to a kept atom, and no other
        let mut ends: HashMap<(NodeIndex, BondOrder, &AtomKey), Vec<NodeIndex>> = HashMap::new();
        for node in source.node_indices() {
            let Some(key) = &keys[node.index()] else {
                continue;
            };
            let mut bonds = source
                .edges(node)
                .filter(|&edge| keys[other(edge, node).index()].is_some());
            if let (Some(edge), None) = (bonds.next(), bonds.next()) {
                if edge.weight().order != BondOrder::Dative {
                    let group = (other(edge, node), edge.weight().order, key);
                    ends.entry(group).or_default().push(node);
                }
            }
        }
        let twins: Vec<(NodeIndex, BondOrder, Vec<NodeIndex>)> = ends
            .into_iter()
            .filter(|(_, group)| group.len() > 1)
            .map(|((parent, order, _), group)| (parent, order, group))
            .collect();
        let mut fates = vec![Fate::Hydrogen; source.node_bound()];
        for (parent, order, group) in twins {
            for twin in group {
                fates[twin.index()] = Fate::Twin;
                if let (Some(key), Some(parent)) =
                    (keys[twin.index()].take(), &mut keys[parent.index()])
                {
                    parent.twins.push((key, order));
                }
            }
        }

        let mut graph = Graph::new_undirected();
        for node in source.node_indices() {
            if let Some(mut key) = keys[node.index()].take() {
                key.twins.sort();
                fates[node.index()] = Fate::Kept(graph.add_node(key).index());
            }
        }
        let mut datives = Vec::new();
        for edge in source.edge_references() {
            let (Fate::Kept(a), Fate::Kept(b)) =
                (fates[edge.source().index()], fates[edge.target().index()])
            else {
                continue;
            };
            graph.add_edge(GraphIndex::new(a), GraphIndex::new(b), edge.weight().order);
            if let Some(donor) = edge.weight().donor {
                datives.push(if donor == edge.source() {
                    (a, b)
                } else {
                    (b, a)
                });
            }
        }
        (graph, fates, datives)
    }

    /// The stereo elements in folded-graph atoms, without those that two alike ligands
    /// rule out.
    fn read_stereo(
        molecule: &Molecule,
        graph: &FoldedGraph,
        fates: &[Fate],
    ) -> Result<Vec<Stereo>, ExactMatchError> {
        let site_of = |node: NodeIndex| match fates[node.index()] {
            Fate::Kept(atom) => Site::Atom(atom),
            Fate::Hydrogen => Site::Hydrogen,
            Fate::Twin => Site::Twin,
        };
        // a hydrogen or twin is an end atom, so it is never a centre or a double-bond end
        let atom = |node: NodeIndex| match fates[node.index()] {
            Fate::Kept(atom) => Ok(atom),
            _ => Err(ExactMatchError::UnsupportedStereo),
        };
        let site = |ligand: Ligand| match ligand {
            Ligand::Atom(node) => Ok(site_of(node)),
            Ligand::Hydrogen => Ok(Site::Hydrogen),
            Ligand::LonePair => Ok(Site::LonePair),
            Ligand::Vacancy => Err(ExactMatchError::UnsupportedStereo),
        };
        // tag 1 is @, 2 is @@
        let parity = |center: NodeIndex, ligands: &[Ligand], tag: u8| {
            let ligands: Vec<Site> = ligands.iter().map(|&l| site(l)).collect::<Result<_, _>>()?;
            Ok(Stereo::Parity {
                center: atom(center)?,
                ligands: ligands
                    .try_into()
                    .map_err(|_| ExactMatchError::UnsupportedStereo)?,
                clockwise: match tag {
                    1 => false,
                    2 => true,
                    _ => return Err(ExactMatchError::UnsupportedStereo),
                },
            })
        };
        let mut stereo = Vec::new();
        for center in molecule.atom_stereo() {
            let Configuration::Tetrahedral(tag) = center.configuration else {
                return Err(ExactMatchError::UnsupportedStereo);
            };
            stereo.push(parity(center.center, &center.ligands, tag)?);
        }
        for allene in molecule.allene_stereo() {
            stereo.push(parity(allene.center, &allene.ligands, allene.permutation)?);
        }
        for bond in molecule.double_bond_stereo() {
            let ends = [atom(bond.ends[0])?, atom(bond.ends[1])?];
            let mut same_side = bond.same_side;
            let mut references = [Site::Hydrogen; 2];
            for side in 0..2 {
                references[side] = site_of(bond.references[side]);
                if references[side] != Site::Hydrogen {
                    continue;
                }
                // a folded hydrogen: refer to the end's other substituent, opposite it
                let other = graph
                    .neighbors(GraphIndex::new(ends[side]))
                    .map(petgraph::graph::NodeIndex::index)
                    .find(|&node| node != ends[1 - side]);
                if let Some(other) = other {
                    references[side] = Site::Atom(other);
                    same_side = !same_side;
                }
            }
            stereo.push(Stereo::DoubleBond {
                ends,
                references,
                same_side,
            });
        }
        // two hydrogens or twins on a centre or a double-bond end make it no stereo at all
        stereo.retain(|element| match element {
            Stereo::Parity { ligands, .. } => {
                !ligands.contains(&Site::Twin)
                    && ligands
                        .iter()
                        .filter(|&&site| site == Site::Hydrogen)
                        .count()
                        < 2
            }
            Stereo::DoubleBond { ends, .. } => ends.iter().all(|&end| {
                let key = &graph[GraphIndex::new(end)];
                key.hydrogens < 2 && key.twins.is_empty()
            }),
        });
        Ok(stereo)
    }

    /// Whether every stereo element agrees with `target`'s under `mapping`. An element
    /// with no counterpart, or one that disagrees, passes only if it is fake.
    fn stereo_agrees(&self, target: &Folded, mapping: &[usize]) -> Result<bool, ExactMatchError> {
        let mut paired = vec![false; target.stereo.len()];
        for (i, element) in self.stereo.iter().enumerate() {
            let counterpart = target.anchors.get(&element.anchor(|atom| mapping[atom]));
            if let Some(&j) = counterpart {
                paired[j] = true;
            }
            let agrees = counterpart.is_some_and(|&j| element.agrees(mapping, &target.stereo[j]));
            if !agrees && !self.fake()?[i] {
                return Ok(false);
            }
        }
        for (j, _) in paired.iter().enumerate().filter(|(_, &paired)| !paired) {
            if !target.fake()?[j] {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn fake(&self) -> Result<&[bool], ExactMatchError> {
        self.fake
            .get_or_init(|| self.find_fake())
            .as_deref()
            .map_err(|error| *error)
    }

    /// Marks each stereo element that a symmetry of the molecule inverts while keeping all
    /// the other stereo, such as `C[C@H](C)O`, whose methyls swap. It describes no
    /// stereoisomer, so it matches a molecule written without it.
    fn find_fake(&self) -> Result<Vec<bool>, ExactMatchError> {
        let mut fake = vec![false; self.stereo.len()];
        if self.stereo.is_empty() {
            return Ok(fake);
        }
        any_mapping(self, self, |mapping| {
            for (i, element) in self.stereo.iter().enumerate() {
                if fake[i]
                    || element.anchor(|atom| mapping[atom]) != element.anchor(|atom| atom)
                    || element.agrees(mapping, element)
                {
                    continue;
                }
                fake[i] = self.stereo.iter().enumerate().all(|(j, other)| {
                    j == i
                        || self
                            .anchors
                            .get(&other.anchor(|atom| mapping[atom]))
                            .is_some_and(|&k| other.agrees(mapping, &self.stereo[k]))
                });
            }
            // every symmetry may mark another element
            Ok(false)
        })?;
        Ok(fake)
    }
}

/// Calls `found` with each mapping of `a`'s atoms onto `b`'s, until it returns true.
fn any_mapping(
    a: &Folded,
    b: &Folded,
    mut found: impl FnMut(&[usize]) -> Result<bool, ExactMatchError>,
) -> Result<bool, ExactMatchError> {
    if a.graph.node_count() != b.graph.node_count() || a.graph.edge_count() != b.graph.edge_count()
    {
        return Ok(false);
    }
    let mut atoms = |x: &AtomKey, y: &AtomKey| x == y;
    let mut bonds = |x: &BondOrder, y: &BondOrder| x == y;
    let (a, b) = (&a.graph, &b.graph);
    // equal atom and bond counts make every subgraph mapping a full isomorphism
    let Some(mappings) = subgraph_isomorphisms_iter(&a, &b, &mut atoms, &mut bonds) else {
        return Ok(false);
    };
    for (tried, mapping) in mappings.enumerate() {
        if tried == MAX_MAPPINGS {
            return Err(ExactMatchError::TooSymmetric);
        }
        if found(&mapping)? {
            return Ok(true);
        }
    }
    Ok(false)
}

impl Molecule {
    /// Whether `other` is the same molecule, however either was written: atom order,
    /// ring-closure digits, Kekulé or aromatic form, hydrogens as atoms or counts. Atoms
    /// must agree on element, isotope, charge, hydrogens and aromaticity, bonds on order
    /// and dative direction, and tetrahedral, allene and double-bond stereo must agree.
    /// Atom classes are labels, so they are ignored.
    pub fn exact_match(&self, other: &Molecule) -> Result<bool, ExactMatchError> {
        let query = Folded::components(self)?;
        let mut target: Vec<Option<Folded>> =
            Folded::components(other)?.into_iter().map(Some).collect();
        if query.len() != target.len() {
            return Ok(false);
        }
        // the first alike target component will do: being the same is an equivalence
        'parts: for part in &query {
            for slot in &mut target {
                if slot.as_ref().map(|other| part.matches(other)).transpose()? == Some(true) {
                    *slot = None;
                    continue 'parts;
                }
            }
            return Ok(false);
        }
        Ok(true)
    }
}
