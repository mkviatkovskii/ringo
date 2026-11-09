//! Writes a molecule as SMILES that `parse_molecule` reads back as the same molecule (#63).
use crate::io::smiles::reader::molecule::{double_bond_stereo, implicit_hydrogens};
use crate::model::bond::{Bond, BondOrder};
use crate::model::molecule::Molecule;
use crate::model::stereo::{Configuration, DoubleBondStereo, Ligand};
use petgraph::stable_graph::{EdgeIndex, EdgeReference, NodeIndex};
use petgraph::visit::{EdgeIndexable, EdgeRef, IntoEdgeReferences, NodeIndexable};
use std::fmt::Write;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WriteError {
    /// A square-planar, trigonal-bipyramidal or octahedral centre whose ligands come out in
    /// another order: renumbering its permutation isn't supported yet.
    UnsupportedStereo,
    /// Double-bond configurations that `/` and `\` can't express in this atom order, such as
    /// one with no single bond to mark at an end, or two that need opposite marks on one bond.
    DoubleBondStereo,
    /// An aromatic atom whose valence aromatic notation can't give back, such as a neutral
    /// nitrogen with a double bond and two single bonds in the ring: written `n`, it reads
    /// back as a lone-pair donor, and the rings around it change.
    AromaticValence,
    /// More than 99 ring bonds open at once.
    TooManyRings,
}

/// Writes `molecule` in reading order: each component starts at its first atom, and each
/// atom's neighbours follow in reading order.
pub fn write_smiles(molecule: &Molecule) -> Result<String, WriteError> {
    let ranks: Vec<usize> = (0..molecule.graph.node_bound()).collect();
    write_smiles_ranked(molecule, &ranks)
}

/// Writes `molecule` in the order of `ranks`, indexed by node index: each component starts
/// at its lowest-ranked atom, and each atom's neighbours follow in rank order. Equal ranks
/// fall back to reading order. Branches come first and the highest-ranked neighbour
/// continues the chain.
pub fn write_smiles_ranked(molecule: &Molecule, ranks: &[usize]) -> Result<String, WriteError> {
    write(molecule, ranks).map(|(smiles, _)| smiles)
}

/// The SMILES, and the atoms in the order written.
pub(crate) fn write(
    molecule: &Molecule,
    ranks: &[usize],
) -> Result<(String, Vec<NodeIndex>), WriteError> {
    assert!(
        ranks.len() >= molecule.graph.node_bound(),
        "a rank for every atom"
    );
    let changed = |node: NodeIndex| {
        molecule
            .input_valences
            .get(node.index())
            .is_some_and(|&valence| valence != molecule.written_valence(node))
    };
    if molecule
        .graph
        .node_indices()
        .any(|node| molecule.graph[node].aromatic && changed(node))
    {
        return Err(WriteError::AromaticValence);
    }
    let writer = Writer::new(molecule, ranks)?;
    let mut out = String::new();
    writer.write(&mut out)?;
    Ok((out, writer.order))
}

/// Whether an aromatic atom of this element is written in lowercase: `b`, `c`, `n`, `o`,
/// `p`, `s`, and in brackets `as`, `se` and `te`.
fn has_lowercase(atomic_number: u8) -> bool {
    matches!(atomic_number, 5 | 6 | 7 | 8 | 15 | 16 | 33 | 34 | 52)
}

/// Whether the element may stand outside brackets (the organic subset), written in
/// lowercase or not.
fn organic(atomic_number: u8, lowercase: bool) -> bool {
    if lowercase {
        matches!(atomic_number, 5 | 6 | 7 | 8 | 15 | 16)
    } else {
        matches!(
            atomic_number,
            0 | 5 | 6 | 7 | 8 | 9 | 15 | 16 | 17 | 35 | 53
        )
    }
}

/// A stereo ligand, with hydrogens told apart by the atom they are on.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Site {
    Ligand(Ligand),
    Hydrogen(NodeIndex),
}

/// Whether reordering `from` into `to` takes an odd number of swaps, or None if they don't
/// hold the same sites.
fn odd_permutation(from: &[Site], to: &[Site]) -> Option<bool> {
    let positions = from
        .iter()
        .map(|site| to.iter().position(|other| other == site))
        .collect::<Option<Vec<usize>>>()?;
    let mut seen = vec![false; to.len()];
    for &position in &positions {
        if std::mem::replace(&mut seen[position], true) {
            return None;
        }
    }
    let inversions = (0..positions.len())
        .flat_map(|i| (i + 1..positions.len()).map(move |j| (i, j)))
        .filter(|&(i, j)| positions[i] > positions[j])
        .count();
    Some(inversions % 2 == 1)
}

/// Whether `stereo` puts `a` at its first end and `b` at its second on the same side.
fn same_side(stereo: &DoubleBondStereo, a: NodeIndex, b: NodeIndex) -> bool {
    let flips = usize::from(a != stereo.references[0]) + usize::from(b != stereo.references[1]);
    stereo.same_side == (flips % 2 == 0)
}

/// A depth-first step: the atom's neighbours with their bonds, the next one to try, and the atom.
type Frame = (Vec<(NodeIndex, EdgeIndex)>, usize, NodeIndex);

struct Writer<'a> {
    molecule: &'a Molecule,
    /// Atoms in the order written.
    order: Vec<NodeIndex>,
    /// Each atom's place in `order`, by node index.
    position: Vec<usize>,
    parent: Vec<Option<NodeIndex>>,
    /// Neighbours reached from each atom, by node index: branches, then the chain.
    children: Vec<Vec<NodeIndex>>,
    /// Ring-bond partners of each atom, by node index, in the order their digits follow it.
    rings: Vec<Vec<NodeIndex>>,
    /// Whether each atom, by node index, is written in lowercase.
    lowercase: Vec<bool>,
    /// The `/` (true) or `\` (false) on each bond, by edge index, written from the atom
    /// written first.
    directions: Vec<Option<bool>>,
}

impl<'a> Writer<'a> {
    fn new(molecule: &'a Molecule, ranks: &[usize]) -> Result<Writer<'a>, WriteError> {
        let graph = &molecule.graph;
        let bound = graph.node_bound();
        let key = |node: &NodeIndex| (ranks[node.index()], node.index());
        let mut writer = Writer {
            molecule,
            order: Vec::with_capacity(graph.node_count()),
            position: vec![usize::MAX; bound],
            parent: vec![None; bound],
            children: vec![Vec::new(); bound],
            rings: vec![Vec::new(); bound],
            lowercase: vec![false; bound],
            directions: vec![None; graph.edge_bound()],
        };
        for node in graph.node_indices() {
            let atom = &graph[node];
            writer.lowercase[node.index()] =
                atom.aromatic && has_lowercase(atom.element.atomic_number);
        }

        // depth first, so every bond outside the tree joins an atom to one written before it
        let mut starts: Vec<NodeIndex> = graph.node_indices().collect();
        starts.sort_by_key(key);
        let mut used = vec![false; graph.edge_bound()];
        let mut closings: Vec<Vec<NodeIndex>> = vec![Vec::new(); bound];
        let mut openings: Vec<Vec<NodeIndex>> = vec![Vec::new(); bound];
        for start in starts {
            if writer.position[start.index()] != usize::MAX {
                continue;
            }
            let mut stack: Vec<Frame> = Vec::new();
            let enter = |writer: &mut Writer<'_>, node: NodeIndex| {
                writer.position[node.index()] = writer.order.len();
                writer.order.push(node);
                let mut next: Vec<(NodeIndex, EdgeIndex)> = graph
                    .edges(node)
                    .map(|edge| {
                        let other = if edge.source() == node {
                            edge.target()
                        } else {
                            edge.source()
                        };
                        (other, edge.id())
                    })
                    .collect();
                next.sort_by_key(|(other, _)| key(other));
                (next, 0, node)
            };
            stack.push(enter(&mut writer, start));
            while let Some((next, i, node)) = stack.last_mut() {
                let node = *node;
                let Some(&(other, edge)) = next.get(*i) else {
                    stack.pop();
                    continue;
                };
                *i += 1;
                if std::mem::replace(&mut used[edge.index()], true) {
                    continue;
                }
                if writer.position[other.index()] == usize::MAX {
                    writer.parent[other.index()] = Some(node);
                    writer.children[node.index()].push(other);
                    stack.push(enter(&mut writer, other));
                } else {
                    closings[node.index()].push(other);
                    openings[other.index()].push(node);
                }
            }
        }
        // ring digits on an atom: the rings it closes, then the ones it opens
        for node in graph.node_indices() {
            let index = node.index();
            let by_position = |nodes: &mut Vec<NodeIndex>| {
                nodes.sort_by_key(|node| writer.position[node.index()]);
                std::mem::take(nodes)
            };
            let mut rings = by_position(&mut closings[index]);
            rings.extend(by_position(&mut openings[index]));
            writer.rings[index] = rings;
        }
        writer.assign_directions()?;
        Ok(writer)
    }

    /// Neighbours in the order the reader lists them: the previous atom, the ring bonds, then
    /// the branches and the chain.
    fn neighbours(&self, node: NodeIndex) -> impl Iterator<Item = NodeIndex> + '_ {
        let index = node.index();
        self.parent[index]
            .into_iter()
            .chain(self.rings[index].iter().copied())
            .chain(self.children[index].iter().copied())
    }

    fn hydrogens(&self, node: NodeIndex) -> u8 {
        let atom = &self.molecule.graph[node];
        atom.hs + atom.implicit_hs
    }

    /// Whether `node` is written after the other atom of a bond.
    fn written_second(&self, node: NodeIndex, other: NodeIndex) -> bool {
        self.position[node.index()] > self.position[other.index()]
    }

    /// Picks `/` and `\` bonds for every double-bond configuration, then checks that the
    /// reader gets the same configurations, and no others, back from them.
    #[allow(
        clippy::too_many_lines,
        reason = "direction selection and parity constraints form one checked assignment"
    )]
    fn assign_directions(&mut self) -> Result<(), WriteError> {
        let molecule = self.molecule;
        if molecule.bond_stereo.is_empty() {
            return Ok(());
        }
        let mut stereo: Vec<&DoubleBondStereo> = molecule.bond_stereo.iter().collect();
        stereo.sort_by_key(|s| {
            s.ends
                .map(|end| self.position[end.index()])
                .into_iter()
                .min()
        });
        let edges = molecule.graph.edge_bound();
        let other_end = |edge: EdgeReference<'_, Bond>, end: NodeIndex| {
            if edge.source() == end {
                edge.target()
            } else {
                edge.source()
            }
        };
        let configured_end =
            |node: NodeIndex| molecule.bond_stereo.iter().any(|s| s.ends.contains(&node));
        let unconfigured = |node: NodeIndex| {
            molecule.graph.edges(node).any(|edge| {
                edge.weight().order == BondOrder::Double
                    && !molecule.bond_stereo.iter().any(|s| s.bond == edge.id())
            })
        };
        // a mark at `node` and one already at the far end of an unconfigured double bond from
        // it would configure that bond
        let completes = |marked: &[bool], node: NodeIndex| {
            molecule.graph.edges(node).any(|double| {
                double.weight().order == BondOrder::Double
                    && !molecule.bond_stereo.iter().any(|s| s.bond == double.id())
                    && molecule
                        .graph
                        .edges(other_end(double, node))
                        .any(|edge| edge.id() != double.id() && marked[edge.id().index()])
            })
        };
        // one mark at each end: a bond already marked, else one that configures no other
        // double bond with the marks so far, else one shared with another configured double
        // bond, else one leading to no unconfigured double bond, which a later mark could
        // configure
        let mut marked = vec![false; edges];
        let mut chosen = Vec::new();
        for s in &stereo {
            for side in 0..2 {
                let (end, partner) = (s.ends[side], s.ends[1 - side]);
                let best = molecule
                    .graph
                    .edges(end)
                    .filter(|&edge| {
                        let other = other_end(edge, end);
                        let markable = match edge.weight().order {
                            BondOrder::Single => true,
                            // `/` between two lowercase atoms reads back as aromatic
                            BondOrder::Aromatic => {
                                self.lowercase[end.index()] && self.lowercase[other.index()]
                            }
                            _ => false,
                        };
                        other != partner && markable
                    })
                    .min_by_key(|&edge| {
                        let other = other_end(edge, end);
                        (
                            !marked[edge.id().index()],
                            completes(&marked, other),
                            !configured_end(other),
                            unconfigured(other),
                            self.position[other.index()],
                        )
                    })
                    .ok_or(WriteError::DoubleBondStereo)?;
                if !marked[best.id().index()] {
                    marked[best.id().index()] = true;
                    chosen.push(best.id().index());
                }
            }
        }
        // Every two marks around a configured double bond fix how their directions relate,
        // `direction[a] != direction[b] == parity`; a mark seen from its end is up when its
        // direction differs from whether that end is written second.
        let mut constraints: Vec<Vec<(usize, bool)>> = vec![Vec::new(); edges];
        for s in &stereo {
            let marks = [0, 1].map(|side| {
                let (end, partner) = (s.ends[side], s.ends[1 - side]);
                molecule
                    .graph
                    .edges(end)
                    .filter(|&edge| marked[edge.id().index()] && other_end(edge, end) != partner)
                    .map(|edge| {
                        let other = other_end(edge, end);
                        (edge.id().index(), other, self.written_second(end, other))
                    })
                    .collect::<Vec<_>>()
            });
            let mut link = |(a, a_second): (usize, bool), (b, b_second): (usize, bool), apart| {
                let parity = apart != (a_second != b_second);
                constraints[a].push((b, parity));
                constraints[b].push((a, parity));
            };
            // two substituents at one end are on opposite sides
            for end in &marks {
                for pair in end.windows(2) {
                    link((pair[0].0, pair[0].2), (pair[1].0, pair[1].2), true);
                }
            }
            let ((a, x, a_second), (b, y, b_second)) = (marks[0][0], marks[1][0]);
            link((a, a_second), (b, b_second), !same_side(s, x, y));
        }
        for &start in &chosen {
            if self.directions[start].is_some() {
                continue;
            }
            self.directions[start] = Some(true);
            let mut stack = vec![start];
            while let Some(a) = stack.pop() {
                let direction = self.directions[a] == Some(true);
                for &(b, parity) in &constraints[a] {
                    match self.directions[b] {
                        None => {
                            self.directions[b] = Some(direction != parity);
                            stack.push(b);
                        }
                        Some(other) if other != (direction != parity) => {
                            return Err(WriteError::DoubleBondStereo);
                        }
                        Some(_) => {}
                    }
                }
            }
        }
        // the reader decides from these, so it has the last word
        let mut written = Vec::new();
        for edge in molecule.graph.edge_references() {
            if let Some(direction) = self.directions[edge.id().index()] {
                let (a, b) = (edge.source(), edge.target());
                written.push(if self.written_second(a, b) {
                    (b, a, direction)
                } else {
                    (a, b, direction)
                });
            }
        }
        let read = double_bond_stereo(molecule, &written).ok_or(WriteError::DoubleBondStereo)?;
        let agrees = |s: &DoubleBondStereo| {
            read.iter().any(|r| {
                if r.ends == s.ends {
                    same_side(r, s.references[0], s.references[1]) == s.same_side
                } else if r.ends == [s.ends[1], s.ends[0]] {
                    same_side(r, s.references[1], s.references[0]) == s.same_side
                } else {
                    false
                }
            })
        };
        if read.len() != molecule.bond_stereo.len() || !molecule.bond_stereo.iter().all(agrees) {
            return Err(WriteError::DoubleBondStereo);
        }
        Ok(())
    }

    fn write(&self, out: &mut String) -> Result<(), WriteError> {
        enum Step {
            Atom(NodeIndex),
            Open,
            Close,
        }
        let mut digits: Vec<Option<u8>> = vec![None; self.molecule.graph.edge_bound()];
        let mut in_use = [false; 100];
        for (i, &start) in self
            .order
            .iter()
            .filter(|node| self.parent[node.index()].is_none())
            .enumerate()
        {
            if i > 0 {
                out.push('.');
            }
            let mut stack = vec![Step::Atom(start)];
            while let Some(step) = stack.pop() {
                let node = match step {
                    Step::Open => {
                        out.push('(');
                        continue;
                    }
                    Step::Close => {
                        out.push(')');
                        continue;
                    }
                    Step::Atom(node) => node,
                };
                if let Some(parent) = self.parent[node.index()] {
                    self.write_bond(out, parent, node);
                }
                self.write_atom(out, node)?;
                let mut closed = Vec::new();
                for &partner in &self.rings[node.index()] {
                    let Some((edge, _)) = self.molecule.graph.find_edge_undirected(node, partner)
                    else {
                        continue;
                    };
                    let digit = if self.written_second(node, partner) {
                        let digit = digits[edge.index()].unwrap_or_default();
                        closed.push(digit);
                        digit
                    } else {
                        // digits closed on this atom stay taken until after it
                        let digit = (1..100u8)
                            .find(|&digit| !in_use[usize::from(digit)])
                            .ok_or(WriteError::TooManyRings)?;
                        in_use[usize::from(digit)] = true;
                        digits[edge.index()] = Some(digit);
                        self.write_bond(out, node, partner);
                        digit
                    };
                    if digit < 10 {
                        out.push(char::from(b'0' + digit));
                    } else {
                        let _ = write!(out, "%{digit}");
                    }
                }
                for digit in closed {
                    in_use[usize::from(digit)] = false;
                }
                if let Some((last, branches)) = self.children[node.index()].split_last() {
                    stack.push(Step::Atom(*last));
                    for &branch in branches.iter().rev() {
                        stack.extend([Step::Close, Step::Atom(branch), Step::Open]);
                    }
                }
            }
        }
        Ok(())
    }

    /// The symbol for the bond from `from` to `to`, `from` written first.
    fn write_bond(&self, out: &mut String, from: NodeIndex, to: NodeIndex) {
        let Some((edge, _)) = self.molecule.graph.find_edge_undirected(from, to) else {
            return;
        };
        let bond = &self.molecule.graph[edge];
        let lowercase = self.lowercase[from.index()] && self.lowercase[to.index()];
        let symbol = match (self.directions[edge.index()], bond.order) {
            (Some(true), _) => "/",
            (Some(false), _) => "\\",
            (None, BondOrder::Single) if lowercase => "-",
            (None, BondOrder::Aromatic) if !lowercase => ":",
            (None, BondOrder::Single | BondOrder::Aromatic) => "",
            (None, BondOrder::Double) => "=",
            (None, BondOrder::Triple) => "#",
            (None, BondOrder::Dative) if bond.donor == Some(from) => "->",
            (None, BondOrder::Dative) => "<-",
        };
        out.push_str(symbol);
    }

    fn write_atom(&self, out: &mut String, node: NodeIndex) -> Result<(), WriteError> {
        let atom = &self.molecule.graph[node];
        let atomic_number = atom.element.atomic_number;
        let lowercase = self.lowercase[node.index()];
        let chirality = self.chirality(node)?;
        let hydrogens = self.hydrogens(node);
        let bare = chirality.is_none()
            && atom.isotope == 0
            && atom.charge == 0
            && atom.atom_class == 0
            && organic(atomic_number, lowercase)
            && implicit_hydrogens(self.molecule, node, atomic_number, lowercase)
                .is_ok_and(|implicit| implicit == hydrogens);
        let symbol = atom.element.symbol();
        if bare {
            if lowercase {
                out.push_str(&symbol.to_ascii_lowercase());
            } else {
                out.push_str(symbol);
            }
            return Ok(());
        }
        out.push('[');
        if atom.isotope > 0 {
            let _ = write!(out, "{}", atom.isotope);
        }
        if lowercase {
            out.push_str(&symbol.to_ascii_lowercase());
        } else {
            out.push_str(symbol);
        }
        out.push_str(&chirality.unwrap_or_default());
        match hydrogens {
            0 => {}
            1 => out.push('H'),
            n => {
                let _ = write!(out, "H{n}");
            }
        }
        match atom.charge {
            0 => {}
            1 => out.push('+'),
            -1 => out.push('-'),
            charge => {
                let _ = write!(out, "{charge:+}");
            }
        }
        if atom.atom_class > 0 {
            let _ = write!(out, ":{}", atom.atom_class);
        }
        out.push(']');
        Ok(())
    }

    /// The `@...` tag for `node`'s configuration in the ligand order of the output.
    fn chirality(&self, node: NodeIndex) -> Result<Option<String>, WriteError> {
        let molecule = self.molecule;
        if let Some(stereo) = molecule.stereo.iter().find(|s| s.center == node) {
            let ligands = self.atom_ligands(node, stereo.configuration);
            let source: Vec<Site> = stereo.ligands.iter().copied().map(Site::Ligand).collect();
            let tag = match stereo.configuration {
                Configuration::Tetrahedral(n) => {
                    let odd = odd_permutation(&source, &ligands)
                        .filter(|_| ligands.len() == 4)
                        .ok_or(WriteError::UnsupportedStereo)?;
                    if (n == 1) == odd { "@@" } else { "@" }.to_owned()
                }
                // shortcut: kept only when the ligands come out in the same order
                _ if source != ligands => return Err(WriteError::UnsupportedStereo),
                Configuration::SquarePlanar(n) => format!("@SP{n}"),
                Configuration::TrigonalBipyramidal(n) => format!("@TB{n}"),
                Configuration::Octahedral(n) => format!("@OH{n}"),
            };
            return Ok(Some(tag));
        }
        if let Some(stereo) = molecule.allene_stereo.iter().find(|s| s.center == node) {
            let (first, last) = (stereo.axis[0], stereo.axis[stereo.axis.len() - 1]);
            let tagged = |ligand: Ligand, terminal: NodeIndex| match ligand {
                Ligand::Hydrogen => Site::Hydrogen(terminal),
                other => Site::Ligand(other),
            };
            let source: Vec<Site> = stereo
                .ligands
                .iter()
                .enumerate()
                .map(|(i, &ligand)| tagged(ligand, if i < 2 { first } else { last }))
                .collect();
            // the arm of the centre's first neighbour comes first
            let middle = stereo.axis.len() / 2;
            let mut ends = [(first, stereo.axis[1]), (last, stereo.axis[middle * 2 - 1])];
            if self.neighbours(node).next() != Some(stereo.axis[middle - 1]) {
                ends.swap(0, 1);
            }
            let mut ligands = Vec::new();
            for (terminal, previous) in ends {
                let mut terminal_ligands: Vec<Site> = self
                    .neighbours(terminal)
                    .filter(|&other| other != previous)
                    .map(|other| Site::Ligand(Ligand::Atom(other)))
                    .collect();
                let position = usize::from(
                    self.parent[terminal.index()].is_some()
                        && self.neighbours(terminal).next() != Some(previous),
                );
                for _ in 0..self.hydrogens(terminal) {
                    terminal_ligands.insert(position, Site::Hydrogen(terminal));
                }
                ligands.extend(terminal_ligands);
            }
            let odd = odd_permutation(&source, &ligands).ok_or(WriteError::UnsupportedStereo)?;
            let tag = if (stereo.permutation == 1) == odd {
                "@@"
            } else {
                "@"
            };
            return Ok(Some(tag.to_owned()));
        }
        Ok(None)
    }

    /// The ligands the reader lists for a centre written at `node` with this kind of tag.
    fn atom_ligands(&self, node: NodeIndex, configuration: Configuration) -> Vec<Site> {
        let mut ligands: Vec<Site> = self
            .neighbours(node)
            .map(|other| Site::Ligand(Ligand::Atom(other)))
            .collect();
        let position = usize::from(self.parent[node.index()].is_some());
        for _ in 0..self.hydrogens(node) {
            ligands.insert(position, Site::Ligand(Ligand::Hydrogen));
        }
        let coordination = match configuration {
            Configuration::Tetrahedral(_) => ligands.len(),
            Configuration::SquarePlanar(_) => 4,
            Configuration::TrigonalBipyramidal(_) => 5,
            Configuration::Octahedral(_) => 6,
        };
        while ligands.len() < coordination {
            ligands.insert(position, Site::Ligand(Ligand::Vacancy));
        }
        let atomic_number = self.molecule.graph[node].element.atomic_number;
        if ligands.len() == 3 && matches!(atomic_number, 7 | 15 | 16 | 33 | 34) {
            ligands.insert(position, Site::Ligand(Ligand::LonePair));
        }
        ligands
    }
}

#[cfg(test)]
mod tests;
