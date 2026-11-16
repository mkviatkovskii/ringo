use crate::model::bond::BondOrder;
use crate::model::matcher::{mix64, mix_word};
use crate::model::molecule::Molecule;
use fixedbitset::FixedBitSet;
use petgraph::stable_graph::{EdgeIndex, NodeIndex};
use petgraph::visit::{EdgeIndexable, EdgeRef, IntoEdgeReferences, NodeIndexable};

const MAX_SHORTEST_RING_PATHS: u8 = 16;

#[derive(Clone, Copy, Default)]
struct RingVisit {
    parent_node: usize,
    parent_edge: usize,
    distance: usize,
    paths: u8,
    stamp: u32,
}

struct Ring {
    atoms: Vec<NodeIndex>,
    edges: Vec<EdgeIndex>,
}

struct CovalentCsr {
    n: usize,
    offsets: Vec<usize>,
    adj: Vec<(usize, usize)>,
    any_written_aromatic: bool,
}

impl CovalentCsr {
    fn neighbors(&self, u: usize) -> &[(usize, usize)] {
        &self.adj[self.offsets[u]..self.offsets[u + 1]]
    }
}

/// What Hückel counting needs to know about the rings an atom is in.
struct RingContext {
    /// Covalent bonds in any ring, by edge index.
    bonds: FixedBitSet,
    /// Size of the smallest ring through each ring atom (0 if not in any ring).
    smallest: Vec<usize>,
}

impl RingContext {
    fn new(node_bound: usize, rings: &[Ring], bonds: FixedBitSet) -> RingContext {
        let mut smallest = vec![0usize; node_bound];
        for ring in rings {
            let len = ring.atoms.len();
            for &node in &ring.atoms {
                let slot = &mut smallest[node.index()];
                if *slot == 0 || len < *slot {
                    *slot = len;
                }
            }
        }
        RingContext { bonds, smallest }
    }
}

/// Hash-indexed deduplication table for sorted ring-edge slices; avoids per-duplicate `Vec` and
/// `RandomState` `HashSet` allocations in `rings_with_edges`.
struct SeenRings {
    entries: Vec<(u64, usize, usize)>,
    pool: Vec<usize>,
}

impl SeenRings {
    fn new() -> Self {
        Self {
            entries: Vec::new(),
            pool: Vec::new(),
        }
    }

    fn insert(&mut self, sorted_edges: &[usize]) -> bool {
        let mut h = sorted_edges.len() as u64;
        for &e in sorted_edges {
            h = mix_word(h, e as u64);
        }
        let hash = mix64(h);
        for &(entry_hash, start, len) in &self.entries {
            if entry_hash == hash && &self.pool[start..start + len] == sorted_edges {
                return false;
            }
        }
        let start = self.pool.len();
        self.pool.extend_from_slice(sorted_edges);
        self.entries.push((hash, start, sorted_edges.len()));
        true
    }
}

impl Molecule {
    /// Marks every ring that obeys Hückel's rule (4n+2 π electrons) as aromatic: its
    /// atoms aromatic, its bonds `Aromatic`. A ring written entirely with aromatic
    /// bonds that fails the rule is still taken as aromatic, as written. Any other
    /// `Aromatic` bond becomes single. Afterwards a molecule has the same atoms and bonds
    /// whether it was written in Kekulé form ("C1=CC=CC=C1") or aromatic form ("c1ccccc1").
    pub fn perceive_aromaticity(&mut self) {
        self.exact_cache.take();
        let (ring_bonds, csr) = self.covalent_ring_data();
        if ring_bonds.is_clear() {
            if csr.any_written_aromatic {
                for edge in self.graph.edge_indices().collect::<Vec<_>>() {
                    if let Some(bond) = self.graph.edge_weight_mut(edge) {
                        if bond.order == BondOrder::Aromatic {
                            bond.order = BondOrder::Single;
                        }
                    }
                }
            }
            return;
        }
        let rings = self.rings_with_edges(&ring_bonds, &csr);
        let context = RingContext::new(self.graph.node_bound(), &rings, ring_bonds);

        let written: Vec<bool> = rings
            .iter()
            .map(|ring| self.written_aromatic(ring))
            .collect();
        let protected: Vec<bool> = rings
            .iter()
            .map(|ring| self.protect_configured_double_bond(ring))
            .collect();
        // Aromatic by itself, so all its bonds are aromatic.
        let mut alone: Vec<bool> = rings
            .iter()
            .zip(&protected)
            .map(|(ring, &prot)| !prot && self.huckel(ring.atoms.iter().copied(), &context))
            .collect();
        let mut aromatic = alone.clone();
        let mut shared_bonds = FixedBitSet::with_capacity(self.graph.edge_bound());
        self.mark_fused_aromatic_pairs(
            &rings,
            &written,
            &protected,
            &alone,
            &context,
            &mut aromatic,
            &mut shared_bonds,
        );
        for (i, &is_written) in written.iter().enumerate() {
            if !aromatic[i] && is_written {
                aromatic[i] = true;
                alone[i] = true;
            }
        }

        let mut aromatic_bonds = FixedBitSet::with_capacity(self.graph.edge_bound());
        for (i, ring) in rings.iter().enumerate().filter(|&(i, _)| aromatic[i]) {
            for &edge in &ring.edges {
                if alone[i] || !shared_bonds.contains(edge.index()) {
                    aromatic_bonds.insert(edge.index());
                }
            }
        }
        self.recover_hidden_aromatic_bonds(&context.bonds, &shared_bonds, &mut aromatic_bonds);
        self.apply_aromaticity(&aromatic_bonds);
    }

    #[expect(
        clippy::too_many_arguments,
        reason = "internal helper splitting fused-pair logic cleanly out of perceive_aromaticity"
    )]
    fn mark_fused_aromatic_pairs(
        &self,
        rings: &[Ring],
        written: &[bool],
        protected: &[bool],
        alone: &[bool],
        context: &RingContext,
        aromatic: &mut [bool],
        shared_bonds: &mut FixedBitSet,
    ) {
        if rings.len() < 2 {
            return;
        }
        let mut atom_mark = vec![0u32; self.graph.node_bound()];
        let mut epoch = 0u32;
        let mut pairs = Vec::new();
        for i in 0..rings.len() {
            if alone[i] || protected[i] {
                continue;
            }
            epoch = next_epoch(&mut atom_mark, epoch);
            for &a in &rings[i].atoms {
                atom_mark[a.index()] = epoch;
            }
            for j in i + 1..rings.len() {
                if alone[j] || protected[j] || written[i] != written[j] {
                    continue;
                }
                let mut shared_count = 0usize;
                let mut shared_a = NodeIndex::end();
                let mut shared_b = NodeIndex::end();
                for &node in &rings[j].atoms {
                    if atom_mark[node.index()] == epoch {
                        match shared_count {
                            0 => shared_a = node,
                            1 => shared_b = node,
                            _ => {}
                        }
                        shared_count += 1;
                        if shared_count > 2 {
                            break;
                        }
                    }
                }
                if shared_count == 2 {
                    if let Some((edge, _)) = self.graph.find_edge_undirected(shared_a, shared_b) {
                        pairs.push((i, j, edge));
                    }
                }
            }
        }
        if pairs.is_empty() {
            return;
        }
        // when one ring could pair with either of two others, as in phenalene, a shared bond
        // written `-` picks the pair, so a written molecule reads back the same in any order
        pairs.sort_by_key(|&(_, _, edge)| self.graph[edge].order == BondOrder::Aromatic);
        let mut union = Vec::new();
        for (i, j, edge) in pairs {
            if aromatic[i] || aromatic[j] {
                continue;
            }
            epoch = next_epoch(&mut atom_mark, epoch);
            union.clear();
            for &node in rings[i].atoms.iter().chain(&rings[j].atoms) {
                let slot = &mut atom_mark[node.index()];
                if *slot != epoch {
                    *slot = epoch;
                    union.push(node);
                }
            }
            if self.huckel(union.iter().copied(), context) {
                aromatic[i] = true;
                aromatic[j] = true;
                shared_bonds.insert(edge.index());
            }
        }
    }

    fn recover_hidden_aromatic_bonds(
        &self,
        ring_bonds: &FixedBitSet,
        shared_bonds: &FixedBitSet,
        aromatic_bonds: &mut FixedBitSet,
    ) {
        let missed: Vec<EdgeIndex> = self
            .graph
            .edge_indices()
            .filter(|&edge| {
                let idx = edge.index();
                ring_bonds.contains(idx)
                    && self.graph[edge].order == BondOrder::Aromatic
                    && !aromatic_bonds.contains(idx)
                    && !shared_bonds.contains(idx)
            })
            .collect();
        if missed.is_empty() {
            return;
        }
        let written_ring_bonds = self.bridges_filtered(|bond| bond == BondOrder::Aromatic);
        for edge in missed {
            if written_ring_bonds.contains(edge.index()) {
                aromatic_bonds.insert(edge.index());
            }
        }
    }

    fn apply_aromaticity(&mut self, aromatic_bonds: &FixedBitSet) {
        for idx in aromatic_bonds.ones() {
            if let Some((a, b)) = self.graph.edge_endpoints(EdgeIndex::new(idx)) {
                self.graph[a].aromatic = true;
                self.graph[b].aromatic = true;
            }
        }
        for edge in self.graph.edge_indices().collect::<Vec<_>>() {
            if let Some(bond) = self.graph.edge_weight_mut(edge) {
                if aromatic_bonds.contains(edge.index()) {
                    bond.order = BondOrder::Aromatic;
                } else if bond.order == BondOrder::Aromatic {
                    // like the bond joining biphenyl's rings, or fluorene's in its 5-ring
                    bond.order = BondOrder::Single;
                }
            }
        }
    }

    /// Whether every bond of the ring was written aromatic: with no symbol between lowercase
    /// atoms, or with `:`. A `-` or `=` stands as written, so the 4-ring of biphenylene,
    /// `c1ccc2c(c1)-c1ccccc1-2`, isn't aromatic although its atoms are.
    fn written_aromatic(&self, ring: &Ring) -> bool {
        ring.edges
            .iter()
            .all(|&edge| self.graph[edge].order == BondOrder::Aromatic)
    }

    // Mixed aromatic and Kekulé notation can explicitly retain a configured double bond.
    fn protect_configured_double_bond(&self, ring: &Ring) -> bool {
        !self.bond_stereo.is_empty()
            && ring
                .edges
                .iter()
                .any(|&edge| self.bond_stereo.iter().any(|s| s.bond == edge))
            && (ring.atoms.iter().any(|&node| self.graph[node].aromatic)
                || ring
                    .edges
                    .iter()
                    .any(|&edge| self.graph[edge].order == BondOrder::Aromatic))
    }

    /// Whether these atoms hold 4n+2 π electrons, each of them giving some.
    fn huckel(&self, atoms: impl Iterator<Item = NodeIndex>, context: &RingContext) -> bool {
        let mut electrons = 0;
        for node in atoms {
            match self.pi_electrons(node, context) {
                Some(count) => electrons += count,
                None => return false,
            }
        }
        electrons % 4 == 2
    }

    /// π electrons an atom gives to an aromatic ring, or None when it cannot be in one
    /// (an sp3 carbon such as the CH2 of cyclopentadiene).
    fn pi_electrons(&self, node: NodeIndex, context: &RingContext) -> Option<u32> {
        let atom = &self.graph[node];
        let mut ring_doubles = 0u8;
        let mut ring_aromatic = false;
        let mut exocyclic_double = None;
        let mut sigma_valence = u32::from(atom.hs) + u32::from(atom.implicit_hs);
        let mut neighbours = usize::from(atom.hs) + usize::from(atom.implicit_hs);
        for edge in self.graph.edges(node) {
            let in_ring = context.bonds.contains(edge.id().index());
            let order = edge.weight().order;
            match order {
                BondOrder::Triple => return None,
                BondOrder::Double if in_ring => ring_doubles += 1,
                BondOrder::Aromatic if in_ring => ring_aromatic = true,
                // an aromatic bond outside rings, like biphenyl's link, ends up single
                BondOrder::Single | BondOrder::Aromatic | BondOrder::Dative => {}
                BondOrder::Double => {
                    let other = if edge.source() == node {
                        edge.target()
                    } else {
                        edge.source()
                    };
                    exocyclic_double = Some(self.graph[other].element.atomic_number);
                }
            }
            if order != BondOrder::Dative {
                neighbours += 1;
                sigma_valence += 1;
            } else if edge.weight().donor != Some(node) {
                sigma_valence += 1;
            }
        }
        // An atom with cumulated or mixed ring/exocyclic double bonds (`=C=`, `=S(=O)`, `=[N+]=`)
        // cannot donate a single π electron without losing a π bond and changing valence.
        if ring_doubles > 1 || (ring_doubles > 0 && exocyclic_double.is_some()) {
            return None;
        }
        match exocyclic_double {
            // the C of a C=O, as in 2-pyridone, gives none
            Some(7 | 8 | 16) => return Some(0),
            // an exocyclic C=C, as in fulvene, takes the electrons out of the ring
            Some(_) => return None,
            None => {}
        }
        if ring_doubles == 1 {
            // A ring double-bond atom whose σ-bond valence already reads back as an aromatic
            // donor or empty-p-orbital atom (`=[S]-`, `=[PH]-`, `=[As]-`, 4-valent `=[N]-`)
            // would lose a valence unit under aromatic notation, so keep the ring Kekulé.
            return (!Self::aromatic_donor_or_empty(
                atom.element.atomic_number,
                atom.charge,
                sigma_valence,
            ))
            .then_some(1);
        }
        match (atom.element.atomic_number, atom.charge) {
            // a lone pair: pyrrole's N, furan's O, thiophene's S, cyclopentadienyl anion
            (7 | 15, 0) | (6, -1) if neighbours == 3 => Some(2),
            // but a neutral O or S between two single bonds only bridges a macrocycle
            // (9 or more atoms, as in RDKit since 2026.09, rdkit/rdkit#9398): the
            // lactone ring of CHEMBL2004566 keeps its C=N double bond
            (8 | 16, 0) if neighbours == 2 && context.smallest[node.index()] >= 9 => None,
            (7, -1) | (8 | 16 | 34 | 52, 0) if neighbours == 2 => Some(2),
            // an empty p orbital: tropylium, boron
            (6, 1) | (5, 0) if neighbours == 3 => Some(0),
            // an atom written aromatic with no lone pair to give: pyridine's n, any c
            _ if ring_aromatic => Some(1),
            _ => None,
        }
    }

    fn build_csr_filtered(&self, keep: impl Fn(BondOrder) -> bool) -> (FixedBitSet, CovalentCsr) {
        let graph = &self.graph;
        let n = graph.node_bound();
        let mut deg = vec![0usize; n];
        let mut ring_bonds = FixedBitSet::with_capacity(graph.edge_bound());
        let mut any_written_aromatic = false;
        let mut edge_count = 0usize;
        for edge in graph.edge_references() {
            let order = edge.weight().order;
            if order == BondOrder::Aromatic {
                any_written_aromatic = true;
            }
            if !keep(order) {
                continue;
            }
            let (u, v) = (edge.source().index(), edge.target().index());
            deg[u] += 1;
            deg[v] += 1;
            ring_bonds.insert(edge.id().index());
            edge_count += 1;
        }
        if edge_count == 0 {
            return (
                ring_bonds,
                CovalentCsr {
                    n,
                    offsets: vec![0; n + 1],
                    adj: Vec::new(),
                    any_written_aromatic,
                },
            );
        }
        let mut offsets = Vec::with_capacity(n + 1);
        offsets.push(0);
        for &d in &deg {
            offsets.push(offsets.last().copied().unwrap_or(0) + d);
        }
        let mut cursor = offsets[..n].to_vec();
        let mut adj = vec![(0usize, 0usize); offsets[n]];
        for edge in graph.edge_references() {
            if !keep(edge.weight().order) {
                continue;
            }
            let e_idx = edge.id().index();
            let (u, v) = (edge.source().index(), edge.target().index());
            adj[cursor[u]] = (e_idx, v);
            cursor[u] += 1;
            adj[cursor[v]] = (e_idx, u);
            cursor[v] += 1;
        }
        (
            ring_bonds,
            CovalentCsr {
                n,
                offsets,
                adj,
                any_written_aromatic,
            },
        )
    }

    fn bridges_filtered(&self, keep: impl Fn(BondOrder) -> bool) -> FixedBitSet {
        let (ring_bonds, csr) = self.covalent_ring_data_with(keep);
        let _ = csr;
        ring_bonds
    }

    fn covalent_ring_data_with(
        &self,
        keep: impl Fn(BondOrder) -> bool,
    ) -> (FixedBitSet, CovalentCsr) {
        let (mut ring_bonds, csr) = self.build_csr_filtered(keep);
        let n = csr.n;
        if csr.adj.is_empty() {
            return (ring_bonds, csr);
        }
        let mut tin = vec![0u32; n];
        let mut low = vec![0u32; n];
        let mut timer = 0u32;
        let mut stack: Vec<(usize, usize, usize, usize)> = Vec::new();
        for start in 0..n {
            if csr.offsets[start] == csr.offsets[start + 1] || tin[start] != 0 {
                continue;
            }
            timer += 1;
            tin[start] = timer;
            low[start] = timer;
            stack.push((
                start,
                csr.offsets[start],
                csr.offsets[start + 1],
                usize::MAX,
            ));
            while let Some(top) = stack.last_mut() {
                let u = top.0;
                if top.1 < top.2 {
                    let (e_idx, v) = csr.adj[top.1];
                    top.1 += 1;
                    if e_idx == top.3 {
                        continue;
                    }
                    if tin[v] != 0 {
                        low[u] = low[u].min(tin[v]);
                    } else {
                        timer += 1;
                        tin[v] = timer;
                        low[v] = timer;
                        stack.push((v, csr.offsets[v], csr.offsets[v + 1], e_idx));
                    }
                } else {
                    let child_low = low[u];
                    let p_edge = top.3;
                    stack.pop();
                    if let Some(&(p, _, _, _)) = stack.last() {
                        low[p] = low[p].min(child_low);
                        if child_low > tin[p] {
                            ring_bonds.set(p_edge, false);
                        }
                    }
                }
            }
        }
        (ring_bonds, csr)
    }

    fn covalent_ring_data(&self) -> (FixedBitSet, CovalentCsr) {
        self.covalent_ring_data_with(|order| order != BondOrder::Dative)
    }

    /// Covalent bonds in cycles: every covalent bond except bridges.
    pub(crate) fn covalent_ring_bonds(&self) -> FixedBitSet {
        self.covalent_ring_data().0
    }

    /// Shortest rings through every covalent ring bond, without duplicates.
    #[cfg(test)]
    fn rings(&self, ring_bonds: &FixedBitSet) -> Vec<Vec<NodeIndex>> {
        let (_, csr) = self.build_csr_filtered(|order| order != BondOrder::Dative);
        self.rings_with_edges(ring_bonds, &csr)
            .into_iter()
            .map(|r| r.atoms)
            .collect()
    }

    fn rings_with_edges(&self, ring_bonds: &FixedBitSet, csr: &CovalentCsr) -> Vec<Ring> {
        let n = csr.n;
        let mut ring_deg = vec![0u16; n];
        for idx in ring_bonds.ones() {
            if let Some((a, b)) = self.graph.edge_endpoints(EdgeIndex::new(idx)) {
                ring_deg[a.index()] += 1;
                ring_deg[b.index()] += 1;
            }
        }

        let mut seen = SeenRings::new();
        let mut rings = Vec::new();
        let mut key_buf = Vec::with_capacity(16);
        let mut record = |atoms: Vec<NodeIndex>, edges: Vec<EdgeIndex>| {
            key_buf.clear();
            key_buf.extend(edges.iter().map(|e| e.index()));
            key_buf.sort_unstable();
            if seen.insert(&key_buf) {
                rings.push(Ring { atoms, edges });
            }
        };

        let mut handled = FixedBitSet::with_capacity(self.graph.edge_bound());
        let mut visited = vec![RingVisit::default(); n];
        let mut stamp = 0u32;
        let mut queue = Vec::new();
        let mut path_atoms = Vec::new();
        let mut path_edges = Vec::new();

        for e_idx in ring_bonds.ones() {
            if handled.contains(e_idx) {
                continue;
            }
            let Some((from_node, to_node)) = self.graph.edge_endpoints(EdgeIndex::new(e_idx))
            else {
                continue;
            };
            let (from, to) = (from_node.index(), to_node.index());
            handled.insert(e_idx);
            mark_deg2_chain(csr, ring_bonds, &ring_deg, &mut handled, from, e_idx);
            mark_deg2_chain(csr, ring_bonds, &ring_deg, &mut handled, to, e_idx);

            let Some(paths) = shortest_path_csr(
                csr,
                from,
                to,
                e_idx,
                ring_bonds,
                &mut visited,
                &mut stamp,
                &mut queue,
                &mut path_atoms,
                &mut path_edges,
            ) else {
                continue;
            };
            match paths {
                1 => {
                    record(path_atoms.clone(), path_edges.clone());
                }
                2..=MAX_SHORTEST_RING_PATHS => {
                    visit_shortest_rings_csr(
                        csr,
                        from,
                        to,
                        e_idx,
                        ring_bonds,
                        &visited,
                        stamp,
                        &mut record,
                    );
                }
                // shortcut: keep two representatives for many tied paths, revisit if these
                // macrocycles need order-independent perception.
                _ => {
                    record(path_atoms.clone(), path_edges.clone());
                    if shortest_path_csr(
                        csr,
                        to,
                        from,
                        e_idx,
                        ring_bonds,
                        &mut visited,
                        &mut stamp,
                        &mut queue,
                        &mut path_atoms,
                        &mut path_edges,
                    )
                    .is_some()
                    {
                        record(path_atoms.clone(), path_edges.clone());
                    }
                }
            }
        }
        rings
    }

    /// The shortest ring of covalent bonds through a bond, as its atoms; `ring_bonds`
    /// comes from `covalent_ring_bonds`.
    pub(crate) fn ring_through(
        &self,
        edge: EdgeIndex,
        ring_bonds: &FixedBitSet,
    ) -> Option<Vec<NodeIndex>> {
        let (a, b) = self.graph.edge_endpoints(edge)?;
        if !ring_bonds.contains(edge.index()) {
            return None;
        }
        let (_, csr) = self.build_csr_filtered(|order| order != BondOrder::Dative);
        let mut visited = vec![RingVisit::default(); csr.n];
        let mut stamp = 0u32;
        let mut queue = Vec::new();
        let mut atoms = Vec::new();
        let mut edges = Vec::new();
        shortest_path_csr(
            &csr,
            a.index(),
            b.index(),
            edge.index(),
            ring_bonds,
            &mut visited,
            &mut stamp,
            &mut queue,
            &mut atoms,
            &mut edges,
        )
        .map(|_| atoms)
    }
}

fn next_epoch(marks: &mut [u32], epoch: u32) -> u32 {
    let next = epoch.wrapping_add(1);
    if next == 0 {
        marks.fill(0);
        return 1;
    }
    next
}

/// Bonds along a maximal degree-2 chain in `ring_bonds` lie in the exact same set of simple
/// cycles, so only one representative bond in the chain needs a shortest-path search.
fn mark_deg2_chain(
    csr: &CovalentCsr,
    ring_bonds: &FixedBitSet,
    ring_deg: &[u16],
    handled: &mut FixedBitSet,
    start: usize,
    start_edge: usize,
) {
    let mut cur = start;
    let mut prev_edge = start_edge;
    while ring_deg[cur] == 2 {
        let Some(&(next_edge, next_node)) = csr
            .neighbors(cur)
            .iter()
            .find(|&&(e_idx, _)| e_idx != prev_edge && ring_bonds.contains(e_idx))
        else {
            break;
        };
        if handled.contains(next_edge) {
            break;
        }
        handled.insert(next_edge);
        prev_edge = next_edge;
        cur = next_node;
    }
}

/// Breadth-first search from `from` to `to` in `csr` that may not use `skip`; writes the
/// shortest cycle's atoms and edges into `out_atoms` / `out_edges` and returns the path
/// multiplicity at `to`.
#[expect(
    clippy::too_many_arguments,
    reason = "reuses caller scratch buffers across all ring bonds without heap allocation"
)]
fn shortest_path_csr(
    csr: &CovalentCsr,
    from: usize,
    to: usize,
    skip: usize,
    ring_bonds: &FixedBitSet,
    visited: &mut [RingVisit],
    stamp: &mut u32,
    queue: &mut Vec<usize>,
    out_atoms: &mut Vec<NodeIndex>,
    out_edges: &mut Vec<EdgeIndex>,
) -> Option<u8> {
    *stamp = stamp.wrapping_add(1);
    if *stamp == 0 {
        visited.fill(RingVisit::default());
        *stamp = 1;
    }
    let cur_stamp = *stamp;
    queue.clear();
    visited[from] = RingVisit {
        parent_node: from,
        parent_edge: skip,
        distance: 0,
        paths: 1,
        stamp: cur_stamp,
    };
    queue.push(from);
    let mut head = 0usize;
    while head < queue.len() {
        let u = queue[head];
        head += 1;
        if u == to {
            out_atoms.clear();
            out_edges.clear();
            out_atoms.push(NodeIndex::new(to));
            out_edges.push(EdgeIndex::new(skip));
            let mut cur = to;
            while cur != from {
                let v = visited[cur];
                out_edges.push(EdgeIndex::new(v.parent_edge));
                cur = v.parent_node;
                out_atoms.push(NodeIndex::new(cur));
            }
            return Some(visited[to].paths);
        }
        let visit = visited[u];
        let distance = visit.distance + 1;
        for &(e_idx, next) in csr.neighbors(u) {
            if e_idx == skip || !ring_bonds.contains(e_idx) {
                continue;
            }
            let slot = &mut visited[next];
            if slot.stamp != cur_stamp {
                *slot = RingVisit {
                    parent_node: u,
                    parent_edge: e_idx,
                    distance,
                    paths: visit.paths,
                    stamp: cur_stamp,
                };
                queue.push(next);
            } else if slot.distance == distance {
                slot.paths = (slot.paths + visit.paths).min(MAX_SHORTEST_RING_PATHS + 1);
            }
        }
    }
    None
}

/// Walk the shortest-path DAG backward, visiting only nodes on a shortest path.
#[expect(
    clippy::too_many_arguments,
    reason = "passes caller scratch state without per-call wrapper struct"
)]
fn visit_shortest_rings_csr(
    csr: &CovalentCsr,
    from: usize,
    to: usize,
    skip: usize,
    ring_bonds: &FixedBitSet,
    visited: &[RingVisit],
    stamp: u32,
    record: &mut impl FnMut(Vec<NodeIndex>, Vec<EdgeIndex>),
) {
    let mut path: Vec<(usize, usize, usize)> = vec![(to, skip, 0)];
    while let Some(&mut (node, _, ref mut cursor)) = path.last_mut() {
        if node == from {
            let atoms = path.iter().map(|&(u, _, _)| NodeIndex::new(u)).collect();
            let edges = path
                .iter()
                .map(|&(_, e_idx, _)| EdgeIndex::new(e_idx))
                .collect();
            record(atoms, edges);
            path.pop();
            continue;
        }
        let nbrs = csr.neighbors(node);
        let distance = visited[node].distance;
        let mut pushed = false;
        while *cursor < nbrs.len() {
            let (e_idx, next) = nbrs[*cursor];
            *cursor += 1;
            if e_idx == skip || !ring_bonds.contains(e_idx) {
                continue;
            }
            let nv = visited[next];
            if nv.stamp == stamp && nv.distance + 1 == distance {
                path.push((next, e_idx, 0));
                pushed = true;
                break;
            }
        }
        if !pushed {
            path.pop();
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::io::smiles::reader::molecule::parse_molecule;
    use crate::math::similarity::tanimoto::tanimoto_bitset;
    use crate::model::bond::BondOrder;
    use petgraph::stable_graph::NodeIndex;
    use std::fmt::Write;

    fn aromatic_atoms(smiles: &str) -> Vec<bool> {
        let m = parse_molecule(smiles).unwrap().1;
        m.graph.node_weights().map(|atom| atom.aromatic).collect()
    }

    #[test]
    fn cycloparaphenylenes_do_not_enumerate_path_combinations() {
        let label = |n: usize| {
            if n < 10 {
                n.to_string()
            } else {
                format!("%{n}")
            }
        };
        for count in [10, 12, 14] {
            let mut smiles = "c12ccc(".to_owned();
            for i in 1..count - 1 {
                write!(smiles, "-c{}ccc(", label(i + 2)).unwrap();
            }
            write!(smiles, "-c{}ccc-2cc{}", label(count + 1), label(count + 1)).unwrap();
            for i in (1..count - 1).rev() {
                write!(smiles, ")cc{}", label(i + 2)).unwrap();
            }
            smiles.push_str(")cc1");
            let start = std::time::Instant::now();
            let molecule = parse_molecule(&smiles).unwrap().1;
            eprintln!("[{count}]CPP: {:?}", start.elapsed());
            assert_eq!(molecule.graph.node_count(), count * 6);
            assert_eq!(
                molecule
                    .graph
                    .edge_weights()
                    .filter(|b| b.order == BondOrder::Aromatic)
                    .count(),
                count * 6
            );
            assert!(
                molecule.rings(&molecule.covalent_ring_bonds()).len()
                    <= 2 * molecule.graph.edge_count()
            );
        }
    }

    #[test]
    fn explicitly_single_biphenylene_links_stay_single() {
        let molecule = parse_molecule("c1ccc2c(c1)-c1ccccc1-2").unwrap().1;
        assert!(molecule.graph.node_weights().all(|a| a.aromatic));
        assert_eq!(
            molecule
                .graph
                .edge_weights()
                .filter(|b| b.order == BondOrder::Single)
                .count(),
            2
        );
        assert_eq!(
            molecule
                .graph
                .edge_weights()
                .filter(|b| b.order == BondOrder::Aromatic)
                .count(),
            12
        );
    }

    #[test]
    fn test_kekule_rings_become_aromatic() {
        for smiles in [
            "C1=CC=CC=C1",
            "C1=CC=NC=C1",
            "C1=CC=CN1",
            "[N-]1C=CC=C1",
            "[N-]1N=NN=C1",
            "C1=COC=C1",
            "C1=CSC=C1",
            "C1=CN=CN1",
            "C1=CC=C2C=CC=CC2=C1",
            // azulene: neither the 5- nor the 7-ring alone, but both together
            "C1=CC2=CC=CC=CC2=C1",
            // 2-pyridone: the C=O carbon gives no electrons, the NH two
            "O=C1C=CC=CN1",
        ] {
            assert!(
                aromatic_atoms(smiles).iter().filter(|&&a| a).count() >= 5,
                "{smiles}"
            );
        }
    }

    #[test]
    fn test_non_aromatic_rings_stay() {
        for smiles in [
            "C1=CCC=C1",
            "C1=CC=CCC1",
            "C1=CC=CC=CC=C1",
            "C1CCCCC1",
            "C=C1C=CC=C1",
            "C1=CC=C=CC=1",
            "FC1=C=C=C(F)C(F)=C1F",
            "C1=NC=[N+]=[N+]=C1",
            "N1=S=NSN=S=NS1",
        ] {
            assert!(aromatic_atoms(smiles).iter().all(|&a| !a), "{smiles}");
        }
        // Fused aromatic rings stay aromatic while adjacent cumulated rings keep their double bonds.
        // Mixed ring/exocyclic double bonds and non-3-coordinate charged C or neutral B also stay non-aromatic.
        for (smiles, expected_aromatic) in [
            ("Cc1ccc2cc1N=C=N2", 6),
            (
                "C=Cc1ccc(N(c2ccc(C=C)cc2)C2C=CC(C)C3=c4c(CC)cccc4=C=C=C32)cc1",
                18,
            ),
            ("CS1(=O)=Nc2nc(C(=O)O)cn2-c2ccccc21", 11),
            ("[B]1Oc2ccccc2O1", 6),
            ("F[B-](F)(F)F.[C+]1Sc2ccccc2S1", 6),
            ("CN1C=CN(C)[CH2+]1Cl.F[B-](F)(F)F", 0),
            ("C1=CC=[S]C=C1", 0),
            ("C1=[PH]c2ccccc2N1", 6),
            ("N[As]1=CC=CC=C1", 0),
            ("[Cu][N]1=CC=CC=C1", 0),
            (
                "C[C@@H]1CN(c2cc(F)c(C3=[S]C(C(=O)N4CCN(C)CC4)=CN3)cc2NC(=O)c2c[nH]c(=O)cc2C(F)(F)F)C[C@H](C)N1C",
                12,
            ),
            (
                "C1=CC=[N]([Fe+2]([N]2=CC=CC=C2)([N]2=CC=CC=C2)([N]2=CC=CC=C2)([N]2=CC=CC=C2)[N]2=CC=CC=C2)C=C1",
                0,
            ),
        ] {
            let count = aromatic_atoms(smiles).iter().filter(|&&a| a).count();
            assert_eq!(count, expected_aromatic, "{smiles}");
        }
    }

    #[test]
    fn test_macrocycle_oxygen_and_sulfur_bridges_give_no_electrons() {
        // oxonin and thionin: 10 by a naive count, but the O or S only bridges a 9-ring
        for smiles in [
            "C1=CC=CC=COC=C1",
            "C1=CC=CC=CSC=C1",
            "O=C1NC=CC(=O)OC(C)=NC=C1",
        ] {
            assert!(aromatic_atoms(smiles).iter().all(|&a| !a), "{smiles}");
        }
        // furan, 2-pyrone and azonine keep their donors; written aromatic stays aromatic
        for (smiles, ring) in [
            ("C1=COC=C1", 5),
            ("O=C1OC=CC=C1", 6),
            ("C1=CC=CC=C[NH]C=C1", 9),
            ("c1ccccoccc1", 9),
        ] {
            let count = aromatic_atoms(smiles).iter().filter(|&&a| a).count();
            assert_eq!(count, ring, "{smiles}");
        }
    }

    #[test]
    fn test_ring_search_with_long_branches_and_disconnected_chains() {
        let chain = "C".repeat(128);
        let smiles = format!("{chain}c1ccccc1{chain}.{chain}");
        let molecule = parse_molecule(&smiles).unwrap().1;
        assert_eq!(
            molecule
                .graph
                .edge_weights()
                .filter(|bond| bond.order == BondOrder::Aromatic)
                .count(),
            6
        );
        assert_eq!(
            molecule.graph.node_weights().filter(|a| a.aromatic).count(),
            6
        );
    }

    #[test]
    fn test_ring_search_reuses_buffers_across_spiro_and_disconnected_rings() {
        let molecule = parse_molecule("C1CC12CC2.C1CC1").unwrap().1;
        let mut rings: Vec<Vec<usize>> = molecule
            .rings(&molecule.covalent_ring_bonds())
            .into_iter()
            .map(|ring| {
                let mut atoms: Vec<usize> = ring.into_iter().map(NodeIndex::index).collect();
                atoms.sort_unstable();
                atoms
            })
            .collect();
        rings.sort_unstable();
        assert_eq!(rings, vec![vec![0, 1, 2], vec![2, 3, 4], vec![5, 6, 7]]);
    }

    #[test]
    fn test_dative_cycle_closures_do_not_create_covalent_ring_bonds() {
        for (smiles, expected_bonds, expected_rings) in [
            ("N->1CC[Cu]1", 0, 0),
            ("N1C=CC=C->[B]<-1", 0, 0),
            ("N->1CC[Cu]1.c1ccccc1", 6, 1),
        ] {
            let molecule = parse_molecule(smiles).unwrap().1;
            let ring_bonds = molecule.covalent_ring_bonds();
            assert_eq!(ring_bonds.count_ones(..), expected_bonds, "{smiles}");
            assert_eq!(
                molecule.rings(&ring_bonds).len(),
                expected_rings,
                "{smiles}"
            );
            for edge in molecule.graph.edge_indices() {
                if molecule.graph[edge].order == BondOrder::Dative {
                    assert!(!ring_bonds.contains(edge.index()), "{smiles}");
                }
            }
        }
    }

    #[test]
    fn test_bonds_outside_aromatic_rings_are_single() {
        let m = parse_molecule("c1ccc2c(c1)Cc1ccccc12").unwrap().1;
        let aromatic = m
            .graph
            .edge_weights()
            .filter(|bond| bond.order == BondOrder::Aromatic)
            .count();
        // two benzene rings; the bond they share with the 5-ring is single
        assert_eq!(aromatic, 12);
        let m = parse_molecule("c1ccccc1:c1ccccc1").unwrap().1;
        let single = m
            .graph
            .edge_weights()
            .filter(|bond| bond.order == BondOrder::Single)
            .count();
        assert_eq!(single, 1);
    }

    #[test]
    fn test_azulene_shared_bond_stays_single() {
        for smiles in ["C1=CC2=CC=CC=CC2=C1", "c1ccc2cccc2cc1", "c1ccc2cccc-2cc1"] {
            let m = parse_molecule(smiles).unwrap().1;
            let count = |order| {
                m.graph
                    .edge_weights()
                    .filter(|bond| bond.order == order)
                    .count()
            };
            // ten aromatic bonds around the outside, the one across the middle single
            assert_eq!(
                (count(BondOrder::Aromatic), count(BondOrder::Single)),
                (10, 1),
                "{smiles}"
            );
            assert!(m.graph.node_weights().all(|atom| atom.aromatic), "{smiles}");
        }
        // naphthalene's rings are aromatic alone, so its shared bond is aromatic
        let m = parse_molecule("C1=CC=C2C=CC=CC2=C1").unwrap().1;
        assert!(m
            .graph
            .edge_weights()
            .all(|bond| bond.order == BondOrder::Aromatic));
    }

    #[test]
    fn test_both_forms_give_the_same_fingerprint() {
        for (kekule, aromatic) in [
            (
                "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
                "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            ),
            ("C1=CC=C2C=CC=CC2=C1", "c1ccc2ccccc2c1"),
            ("C1=CC=NC=C1", "c1ccncc1"),
            ("C1=CC=CN1", "c1cc[nH]c1"),
            ("[N-]1C=CC=C1", "[n-]1cccc1"),
            ("[N-]1N=NN=C1", "[n-]1nnnc1"),
            ("C1=CSC=C1", "c1ccsc1"),
            ("C1=CC=C2C(=C1)C=CN2", "c1ccc2[nH]ccc2c1"),
            ("C1=CC=C(C=C1)C1=CC=CC=C1", "c1ccccc1c1ccccc1"),
            ("C1C2=CC=CC=C2C2=CC=CC=C12", "c1ccc2c(c1)Cc1ccccc1-2"),
            ("c1ccc2c(c1)Cc1ccccc12", "c1ccc2c(c1)Cc1ccccc1-2"),
        ] {
            let a = parse_molecule(kekule).unwrap().1.ecfp(2, 512);
            let b = parse_molecule(aromatic).unwrap().1.ecfp(2, 512);
            assert_eq!(tanimoto_bitset(&a.0, &b.0), 1.0, "{kekule} {aromatic}");
        }
    }
}
