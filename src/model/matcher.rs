//! Symmetry-pruned graph and subgraph matching primitives shared by exact search and future
//! substructure search (#57).

use crate::model::bond::BondOrder;
use std::collections::BTreeMap;

/// Fast, platform-independent 64-bit integer mixer (`SplitMix64` finalizer).
#[inline]
pub(crate) const fn mix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9e37_79b9_7f4a_7c15);
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

/// Folds a 64-bit word into a running 64-bit hash state without byte-by-byte iteration.
#[inline]
pub(crate) const fn mix_word(state: u64, word: u64) -> u64 {
    state.rotate_left(13).wrapping_mul(0x517c_c1b7_2722_0a95) ^ word
}

/// Packs bond order and donor orientation into a single comparable byte.
#[must_use]
pub(crate) fn pack_bond(order: BondOrder, donor_is_self: bool) -> u8 {
    let base = match order {
        BondOrder::Single => 0,
        BondOrder::Double => 1,
        BondOrder::Triple => 2,
        BondOrder::Aromatic => 3,
        BondOrder::Dative => 4,
    };
    (base << 1) | u8::from(donor_is_self)
}

/// Dense CSR graph with precomputed per-vertex base hashes for fast refinement and matching.
pub(crate) struct CsrGraph {
    pub n: usize,
    pub offsets: Vec<usize>,
    pub edges: Vec<(u8, usize)>,
    pub base_hashes: Vec<u64>,
}

impl CsrGraph {
    #[allow(
        dead_code,
        reason = "convenience builder for #57 substructure query/target construction"
    )]
    pub(crate) fn from_adj(n: usize, adj: &[Vec<(u8, usize)>], base_hashes: Vec<u64>) -> Self {
        let mut offsets = Vec::with_capacity(n + 1);
        let total_edges: usize = adj.iter().map(Vec::len).sum();
        let mut edges = Vec::with_capacity(total_edges);
        offsets.push(0);
        for nbrs in adj {
            edges.extend_from_slice(nbrs);
            offsets.push(edges.len());
        }
        Self {
            n,
            offsets,
            edges,
            base_hashes,
        }
    }

    pub(crate) fn neighbors(&self, u: usize) -> &[(u8, usize)] {
        &self.edges[self.offsets[u]..self.offsets[u + 1]]
    }

    pub(crate) fn edge_between(&self, u: usize, v: usize) -> Option<u8> {
        self.neighbors(u)
            .iter()
            .find_map(|&(code, nbr)| (nbr == v).then_some(code))
    }
}

/// An automorphism permutation together with the subset of vertices it moves.
#[derive(Clone)]
pub(crate) struct Symmetry {
    pub map: Vec<usize>,
    pub support: Vec<usize>,
}

impl Symmetry {
    pub(crate) fn from_map(map: Vec<usize>) -> Self {
        let support = map
            .iter()
            .enumerate()
            .filter_map(|(i, &target)| (target != i).then_some(i))
            .collect();
        Self { map, support }
    }

    pub(crate) fn preserves_fixed(&self, fixed: &[bool]) -> bool {
        !self.support.iter().any(|&idx| fixed[idx])
    }

    pub(crate) fn preserves_colors(&self, colors: &[usize]) -> bool {
        self.support
            .iter()
            .all(|&idx| colors[idx] == colors[self.map[idx]])
    }
}

/// Automorphism generators of a graph, providing global symmetry classes and dynamic
/// stabilizer orbits during backtracking.
pub(crate) struct SymmetryGroup {
    pub n: usize,
    pub generators: Vec<Symmetry>,
}

impl SymmetryGroup {
    /// Computes symmetry generators for `graph` using partition refinement and nauty-style
    /// automorphism pruning. `stereo_sig` may inject optional per-vertex stereo keys (ordered
    /// slice per vertex) into refinement and leaf comparison.
    pub(crate) fn compute<S, L>(
        graph: &CsrGraph,
        stereo_free_refs: &[bool],
        mut stereo_keys: S,
        mut leaf_stereo: L,
    ) -> Self
    where
        S: FnMut(&[usize], &mut Vec<Vec<u64>>),
        L: FnMut(&[usize], &mut Vec<u64>),
    {
        let n = graph.n;
        if n <= 1 {
            return Self {
                n,
                generators: Vec::new(),
            };
        }
        let mut engine = SymEngine::new(graph);
        let (mut colors, mut distinct) = engine.initial_colors();
        distinct = engine.refine(&mut colors, distinct, &mut stereo_keys);
        engine.seed_twin_symmetries(&colors, stereo_free_refs);
        engine.visit(colors, distinct, &mut stereo_keys, &mut leaf_stereo);
        Self {
            n,
            generators: engine.symmetries,
        }
    }

    /// Global vertex symmetry classes (`0..num_orbits`).
    #[allow(
        dead_code,
        reason = "public building block for #57 substructure search"
    )]
    pub(crate) fn vertex_orbits(&self) -> Vec<usize> {
        let fixed = vec![false; self.n];
        let mut orbits = vec![0usize; self.n];
        let mut parent: Vec<usize> = (0..self.n).collect();
        self.stabilizer_orbits(&fixed, &mut orbits, &mut parent);
        orbits
    }

    /// Computes connected components (orbits) of vertices under generators that fix every
    /// vertex marked `true` in `fixed`. Reuses caller-provided `orbits` and `parent` buffers.
    pub(crate) fn stabilizer_orbits(
        &self,
        fixed: &[bool],
        orbits: &mut [usize],
        parent: &mut [usize],
    ) {
        for (i, slot) in parent.iter_mut().enumerate().take(self.n) {
            *slot = i;
        }
        for sym in &self.generators {
            if !sym.preserves_fixed(fixed) {
                continue;
            }
            for &u in &sym.support {
                let v = sym.map[u];
                union_roots(parent, u, v);
            }
        }
        for (i, orbit) in orbits.iter_mut().enumerate().take(self.n) {
            *orbit = find_root(parent, i);
        }
    }
}

fn find_root(parent: &mut [usize], mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    x
}

fn union_roots(parent: &mut [usize], a: usize, b: usize) {
    let ra = find_root(parent, a);
    let rb = find_root(parent, b);
    if ra != rb {
        if ra < rb {
            parent[rb] = ra;
        } else {
            parent[ra] = rb;
        }
    }
}

/// Refines a `CsrGraph` to stable color hashes `(colors, color_hashes, color_hist)`, returning
/// dense vertex color IDs `0..distinct` sorted canonically so isomorphic graphs get identical
/// `color_hist`.
pub(crate) fn refine_hash_signature(graph: &CsrGraph) -> (Vec<usize>, Vec<u64>, Vec<(u64, usize)>) {
    let n = graph.n;
    if n == 0 {
        return (Vec::new(), Vec::new(), Vec::new());
    }
    let mut hashes = graph.base_hashes.clone();
    let mut scratch = Vec::with_capacity(n);
    let mut distinct = count_distinct_u64(&mut scratch, &hashes);
    if distinct < n {
        let mut nbrs: Vec<u64> = Vec::with_capacity(8);
        let mut next_hashes = vec![0u64; n];
        for _ in 0..n {
            for u in 0..n {
                nbrs.clear();
                for &(code, v) in graph.neighbors(u) {
                    let edge_mix =
                        mix64(hashes[v] ^ (u64::from(code).wrapping_mul(0x9e37_79b9_7f4a_7c15)));
                    nbrs.push(edge_mix);
                }
                if nbrs.len() > 1 {
                    nbrs.sort_unstable();
                }
                let mut h = mix_word(hashes[u], nbrs.len() as u64);
                for &edge_mix in &nbrs {
                    h = mix_word(h, edge_mix);
                }
                next_hashes[u] = mix64(h);
            }
            std::mem::swap(&mut hashes, &mut next_hashes);
            let next_distinct = count_distinct_u64(&mut scratch, &hashes);
            if next_distinct <= distinct {
                break;
            }
            distinct = next_distinct;
            if distinct >= n {
                break;
            }
        }
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_unstable_by_key(|&u| hashes[u]);
    let mut colors = vec![0usize; n];
    let mut color_hashes = Vec::with_capacity(distinct);
    let mut color_hist: Vec<(u64, usize)> = Vec::with_capacity(distinct);
    for (pos, &u) in order.iter().enumerate() {
        let h = hashes[u];
        if pos == 0 || h != hashes[order[pos - 1]] {
            color_hashes.push(h);
            color_hist.push((h, 1));
        } else if let Some(last) = color_hist.last_mut() {
            last.1 += 1;
        }
        colors[u] = color_hashes.len() - 1;
    }
    (colors, color_hashes, color_hist)
}

fn count_distinct_u64(scratch: &mut Vec<u64>, values: &[u64]) -> usize {
    scratch.clear();
    scratch.extend_from_slice(values);
    scratch.sort_unstable();
    scratch.dedup();
    scratch.len()
}

struct SymEngine<'a> {
    graph: &'a CsrGraph,
    symmetries: Vec<Symmetry>,
    leaves: BTreeMap<Vec<u64>, Vec<usize>>,
    color_counts: Vec<usize>,
    cell_starts: Vec<usize>,
    order: Vec<usize>,
    next_colors: Vec<usize>,
    nbr_buf: Vec<u64>,
    nbr_ranges: Vec<(usize, usize)>,
    stereo_buf: Vec<Vec<u64>>,
    covered_mark: Vec<u32>,
    covered_epoch: u32,
    covered_stack: Vec<usize>,
}

impl<'a> SymEngine<'a> {
    fn new(graph: &'a CsrGraph) -> Self {
        let n = graph.n;
        Self {
            graph,
            symmetries: Vec::new(),
            leaves: BTreeMap::new(),
            color_counts: vec![0; n + 1],
            cell_starts: vec![0; n + 1],
            order: (0..n).collect(),
            next_colors: vec![0; n],
            nbr_buf: Vec::new(),
            nbr_ranges: vec![(0, 0); n],
            stereo_buf: (0..n).map(|_| Vec::new()).collect(),
            covered_mark: vec![0; n],
            covered_epoch: 0,
            covered_stack: Vec::new(),
        }
    }

    fn initial_colors(&self) -> (Vec<usize>, usize) {
        let n = self.graph.n;
        let mut codes: Vec<Vec<u8>> = Vec::with_capacity(n);
        for u in 0..n {
            let mut row: Vec<u8> = self.graph.neighbors(u).iter().map(|&(c, _)| c).collect();
            row.sort_unstable();
            codes.push(row);
        }
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_unstable_by(|&a, &b| {
            (self.graph.base_hashes[a], &codes[a]).cmp(&(self.graph.base_hashes[b], &codes[b]))
        });
        let mut colors = vec![0usize; n];
        let mut c = 0usize;
        for i in 0..n {
            if i > 0
                && (self.graph.base_hashes[order[i]], &codes[order[i]])
                    != (self.graph.base_hashes[order[i - 1]], &codes[order[i - 1]])
            {
                c += 1;
            }
            colors[order[i]] = c;
        }
        (colors, c + 1)
    }

    fn seed_twin_symmetries(&mut self, colors: &[usize], stereo_free: &[bool]) {
        let n = self.graph.n;
        let mut candidates: Vec<((usize, usize, u8), usize)> = Vec::new();
        for u in 0..n {
            if !stereo_free[u] {
                continue;
            }
            let nbrs = self.graph.neighbors(u);
            if nbrs.len() == 1 {
                let (code, parent) = nbrs[0];
                candidates.push(((colors[u], parent, code), u));
            }
        }
        candidates.sort_unstable();
        for pair in candidates.windows(2) {
            if pair[0].0 == pair[1].0 {
                let (a, b) = (pair[0].1, pair[1].1);
                let mut map: Vec<usize> = (0..n).collect();
                map.swap(a, b);
                self.symmetries.push(Symmetry::from_map(map));
            }
        }
    }

    fn refine<S>(&mut self, colors: &mut [usize], mut distinct: usize, stereo_keys: &mut S) -> usize
    where
        S: FnMut(&[usize], &mut Vec<Vec<u64>>),
    {
        let n = self.graph.n;
        loop {
            if distinct >= n {
                return distinct;
            }
            self.color_counts[..distinct].fill(0);
            for &c in colors.iter() {
                self.color_counts[c] += 1;
            }
            for row in &mut self.stereo_buf {
                row.clear();
            }
            stereo_keys(colors, &mut self.stereo_buf);

            self.nbr_buf.clear();
            for u in 0..n {
                if self.color_counts[colors[u]] <= 1 {
                    self.nbr_ranges[u] = (0, 0);
                    continue;
                }
                let start = self.nbr_buf.len();
                for &(code, v) in self.graph.neighbors(u) {
                    self.nbr_buf
                        .push((u64::from(code) << 56) | (colors[v] as u64));
                }
                self.nbr_buf[start..].sort_unstable();
                self.nbr_ranges[u] = (start, self.nbr_buf.len());
            }

            let mut offset = 0usize;
            for c in 0..distinct {
                self.cell_starts[c] = offset;
                offset += self.color_counts[c];
            }
            for (u, &c) in colors.iter().enumerate().take(n) {
                let slot = self.cell_starts[c];
                self.cell_starts[c] = slot + 1;
                self.order[slot] = u;
            }
            let mut start = 0usize;
            for c in 0..distinct {
                let count = self.color_counts[c];
                let end = start + count;
                if count > 1 {
                    let nbr_buf = &self.nbr_buf;
                    let nbr_ranges = &self.nbr_ranges;
                    let stereo_buf = &self.stereo_buf;
                    self.order[start..end].sort_unstable_by(|&a, &b| {
                        let (sa, ea) = nbr_ranges[a];
                        let (sb, eb) = nbr_ranges[b];
                        nbr_buf[sa..ea]
                            .cmp(&nbr_buf[sb..eb])
                            .then_with(|| stereo_buf[a].cmp(&stereo_buf[b]))
                    });
                }
                start = end;
            }

            let mut next_color = 0usize;
            self.next_colors[self.order[0]] = 0;
            for idx in 1..n {
                let prev = self.order[idx - 1];
                let cur = self.order[idx];
                let differs = colors[prev] != colors[cur]
                    || self.nbr_buf[self.nbr_ranges[prev].0..self.nbr_ranges[prev].1]
                        != self.nbr_buf[self.nbr_ranges[cur].0..self.nbr_ranges[cur].1]
                    || self.stereo_buf[prev] != self.stereo_buf[cur];
                if differs {
                    next_color += 1;
                }
                self.next_colors[cur] = next_color;
            }
            let next_distinct = next_color + 1;
            colors.copy_from_slice(&self.next_colors);
            if next_distinct == distinct {
                return next_distinct;
            }
            distinct = next_distinct;
        }
    }

    fn covered(&mut self, node: usize, explored: &[usize], preserving: &[usize]) -> bool {
        self.covered_epoch = self.covered_epoch.wrapping_add(1);
        if self.covered_epoch == 0 {
            self.covered_mark.fill(0);
            self.covered_epoch = 1;
        }
        let epoch = self.covered_epoch;
        self.covered_stack.clear();
        self.covered_stack.push(node);
        self.covered_mark[node] = epoch;
        while let Some(u) = self.covered_stack.pop() {
            if explored.contains(&u) {
                return true;
            }
            for &gen_idx in preserving {
                let next = self.symmetries[gen_idx].map[u];
                if self.covered_mark[next] != epoch {
                    self.covered_mark[next] = epoch;
                    self.covered_stack.push(next);
                }
            }
        }
        false
    }

    fn leaf_signature<L>(&self, rank_of: &[usize], leaf_stereo: &mut L) -> (Vec<u64>, Vec<usize>)
    where
        L: FnMut(&[usize], &mut Vec<u64>),
    {
        let n = self.graph.n;
        let mut order = vec![0usize; n];
        for (u, &r) in rank_of.iter().enumerate() {
            order[r] = u;
        }
        let mut sig = Vec::with_capacity(n + self.graph.edges.len());
        for &u in &order {
            sig.push(self.graph.base_hashes[u]);
        }
        let mut edges_sig = Vec::new();
        for &u in &order {
            let ru = rank_of[u];
            let start = edges_sig.len();
            for &(code, v) in self.graph.neighbors(u) {
                let rv = rank_of[v];
                edges_sig.push(((ru as u64) << 40) | ((rv as u64) << 16) | u64::from(code));
            }
            edges_sig[start..].sort_unstable();
        }
        sig.extend_from_slice(&edges_sig);
        leaf_stereo(rank_of, &mut sig);
        (sig, order)
    }

    fn visit<S, L>(
        &mut self,
        mut colors: Vec<usize>,
        distinct: usize,
        stereo_keys: &mut S,
        leaf_stereo: &mut L,
    ) where
        S: FnMut(&[usize], &mut Vec<Vec<u64>>),
        L: FnMut(&[usize], &mut Vec<u64>),
    {
        let distinct = self.refine(&mut colors, distinct, stereo_keys);
        let n = self.graph.n;
        if distinct >= n {
            let (sig, order) = self.leaf_signature(&colors, leaf_stereo);
            if let Some(previous) = self.leaves.get(&sig) {
                let mut map: Vec<usize> = (0..n).collect();
                for (&a, &b) in previous.iter().zip(&order) {
                    map[a] = b;
                }
                self.symmetries.push(Symmetry::from_map(map));
            } else {
                self.leaves.insert(sig, order);
            }
            return;
        }

        self.color_counts[..distinct].fill(0);
        for &c in &colors {
            self.color_counts[c] += 1;
        }
        let Some((target_color, &cell_len)) = self.color_counts[..distinct]
            .iter()
            .enumerate()
            .filter(|&(_, &count)| count > 1)
            .min_by_key(|&(_, &count)| count)
        else {
            return;
        };
        let mut cell = Vec::with_capacity(cell_len);
        for (u, &c) in colors.iter().enumerate() {
            if c == target_color {
                cell.push(u);
            }
        }

        let mut explored = Vec::new();
        let mut checked = 0usize;
        let mut preserving = Vec::new();
        for node in cell {
            for (idx, sym) in self.symmetries.iter().enumerate().skip(checked) {
                if sym.preserves_colors(&colors) {
                    preserving.push(idx);
                }
            }
            checked = self.symmetries.len();
            if self.covered(node, &explored, &preserving) {
                continue;
            }
            explored.push(node);
            let mut child = colors.clone();
            for c in &mut child {
                if *c > target_color {
                    *c += 1;
                }
            }
            child[node] += 1;
            self.visit(child, distinct + 1, stereo_keys, leaf_stereo);
        }
    }
}

/// Strategy hooks allowing [`find_embedding`] to serve both exact isomorphism and future
/// substructure matching (#57).
pub(crate) trait EmbeddingSpec {
    fn atom_compatible(&self, q: usize, t: usize) -> bool;
    fn bond_compatible(&self, q_code: u8, t_code: u8) -> bool;
    /// Called immediately after `mapping[q]` has been tentatively set (`target_used[t] = true`).
    fn step_ok(&self, q: usize, mapping: &[usize]) -> bool;
    /// Optional final check once all query vertices are mapped.
    fn leaf_ok(&self, mapping: &[usize]) -> bool;
}

/// Computes a connected, MRV-ordered traversal sequence `(q_vertex, mapped_parent_idx)` over
/// `q_graph` starting from vertices with the fewest candidates (`domain_sizes`).
/// Computes a connected, MRV-ordered traversal sequence `(q_vertex, mapped_parent_idx)` over
/// `q_graph` starting from vertices with the fewest candidates (`domain_sizes`), breaking ties
/// by highest degree and most-recently-reached neighbor (depth-first within a domain tier) so
/// cyclic/chiral substructures close before sibling symmetric branches fan out.
pub(crate) fn plan_query_order(
    q_graph: &CsrGraph,
    domain_sizes: &[usize],
) -> Vec<(usize, Option<usize>)> {
    let n = q_graph.n;
    let mut visited = vec![false; n];
    let mut plan = Vec::with_capacity(n);
    let mut frontier: Vec<(usize, usize, usize)> = Vec::new();
    let mut seq = 0usize;
    while plan.len() < n {
        let start = (0..n)
            .filter(|&u| !visited[u])
            .min_by_key(|&u| {
                (
                    domain_sizes[u],
                    std::cmp::Reverse(q_graph.neighbors(u).len()),
                )
            })
            .unwrap_or(0);
        visited[start] = true;
        plan.push((start, None));
        frontier.clear();
        for &(_, v) in q_graph.neighbors(start) {
            if !visited[v] {
                visited[v] = true;
                frontier.push((v, start, seq));
                seq += 1;
            }
        }
        while !frontier.is_empty() {
            let best_idx = (0..frontier.len())
                .min_by_key(|&i| {
                    let (v, _, s) = frontier[i];
                    (
                        domain_sizes[v],
                        std::cmp::Reverse(q_graph.neighbors(v).len()),
                        std::cmp::Reverse(s),
                    )
                })
                .unwrap_or(0);
            let (u, parent, _) = frontier.swap_remove(best_idx);
            plan.push((u, Some(parent)));
            for &(_, v) in q_graph.neighbors(u) {
                if !visited[v] {
                    visited[v] = true;
                    frontier.push((v, u, seq));
                    seq += 1;
                }
            }
        }
    }
    plan
}

pub(crate) enum EmbedOutcome {
    Matched(Vec<usize>),
    NoMatch,
    BudgetExceeded,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum SearchStatus {
    Found,
    Exhausted,
    BudgetExceeded,
}

/// Searches for an injective embedding of `q_graph` into `t_graph` honoring `spec`, pruning
/// redundant target branches via `t_sym` stabilizer orbits.
pub(crate) fn find_embedding<S: EmbeddingSpec>(
    q_graph: &CsrGraph,
    t_graph: &CsrGraph,
    t_sym: &SymmetryGroup,
    plan: &[(usize, Option<usize>)],
    spec: &S,
) -> Option<Vec<usize>> {
    match run_embedding(q_graph, t_graph, Some(t_sym), plan, spec, None) {
        EmbedOutcome::Matched(mapping) => Some(mapping),
        EmbedOutcome::NoMatch | EmbedOutcome::BudgetExceeded => None,
    }
}

/// Fast unpruned search with a failed-branch budget, avoiding `SymmetryGroup` construction on
/// typical hits and simple non-hits.
pub(crate) fn find_embedding_unpruned<S: EmbeddingSpec>(
    q_graph: &CsrGraph,
    t_graph: &CsrGraph,
    plan: &[(usize, Option<usize>)],
    spec: &S,
    backtrack_budget: usize,
) -> EmbedOutcome {
    run_embedding(q_graph, t_graph, None, plan, spec, Some(backtrack_budget))
}

fn run_embedding<S: EmbeddingSpec>(
    q_graph: &CsrGraph,
    t_graph: &CsrGraph,
    t_sym: Option<&SymmetryGroup>,
    plan: &[(usize, Option<usize>)],
    spec: &S,
    budget: Option<usize>,
) -> EmbedOutcome {
    let has_sym = t_sym.is_some_and(|sym| !sym.generators.is_empty());
    let mut ctx = EmbedState {
        q_graph,
        t_graph,
        t_sym,
        has_sym,
        plan,
        spec,
        mapping: vec![usize::MAX; q_graph.n],
        target_used: vec![false; t_graph.n],
        orbits: if has_sym {
            vec![0; t_graph.n]
        } else {
            Vec::new()
        },
        uf_parent: if has_sym {
            vec![0; t_graph.n]
        } else {
            Vec::new()
        },
        cand_buf: Vec::with_capacity(q_graph.n),
        tried_buf: Vec::new(),
        budget,
    };
    match ctx.search(0) {
        SearchStatus::Found => EmbedOutcome::Matched(ctx.mapping),
        SearchStatus::Exhausted => EmbedOutcome::NoMatch,
        SearchStatus::BudgetExceeded => EmbedOutcome::BudgetExceeded,
    }
}

struct EmbedState<'a, S: EmbeddingSpec> {
    q_graph: &'a CsrGraph,
    t_graph: &'a CsrGraph,
    t_sym: Option<&'a SymmetryGroup>,
    has_sym: bool,
    plan: &'a [(usize, Option<usize>)],
    spec: &'a S,
    mapping: Vec<usize>,
    target_used: Vec<bool>,
    orbits: Vec<usize>,
    uf_parent: Vec<usize>,
    cand_buf: Vec<usize>,
    tried_buf: Vec<usize>,
    budget: Option<usize>,
}

impl<S: EmbeddingSpec> EmbedState<'_, S> {
    fn push_candidates(&mut self, q: usize, parent: Option<usize>) {
        if let Some(pq) = parent {
            let pt = self.mapping[pq];
            let Some(q_edge_from_parent) = self.q_graph.edge_between(pq, q) else {
                return;
            };
            for &(t_code, t) in self.t_graph.neighbors(pt) {
                if !self.target_used[t]
                    && self.spec.bond_compatible(q_edge_from_parent, t_code)
                    && self.edges_compatible(q, t)
                {
                    self.cand_buf.push(t);
                }
            }
        } else {
            for t in 0..self.t_graph.n {
                if !self.target_used[t] && self.edges_compatible(q, t) {
                    self.cand_buf.push(t);
                }
            }
        }
    }

    fn edges_compatible(&self, q: usize, t: usize) -> bool {
        if !self.spec.atom_compatible(q, t) {
            return false;
        }
        if self.q_graph.neighbors(q).len() > self.t_graph.neighbors(t).len() {
            return false;
        }
        for &(q_code, q_nbr) in self.q_graph.neighbors(q) {
            let t_nbr = self.mapping[q_nbr];
            if t_nbr != usize::MAX {
                let Some(t_code) = self.t_graph.edge_between(t, t_nbr) else {
                    return false;
                };
                if !self.spec.bond_compatible(q_code, t_code) {
                    return false;
                }
            }
        }
        true
    }

    fn search(&mut self, depth: usize) -> SearchStatus {
        if depth == self.plan.len() {
            return if self.spec.leaf_ok(&self.mapping) {
                SearchStatus::Found
            } else {
                SearchStatus::Exhausted
            };
        }
        let (q, parent) = self.plan[depth];
        let cand_base = self.cand_buf.len();
        self.push_candidates(q, parent);
        let cand_end = self.cand_buf.len();
        if cand_base == cand_end {
            return SearchStatus::Exhausted;
        }
        let use_orbits = self.has_sym && (cand_end - cand_base) > 1;
        if use_orbits {
            if let Some(sym) = self.t_sym {
                sym.stabilizer_orbits(&self.target_used, &mut self.orbits, &mut self.uf_parent);
            }
        }
        let tried_base = self.tried_buf.len();
        for idx in cand_base..cand_end {
            let t = self.cand_buf[idx];
            let orb = if use_orbits { self.orbits[t] } else { t };
            if use_orbits && self.tried_buf[tried_base..].contains(&orb) {
                continue;
            }
            self.mapping[q] = t;
            self.target_used[t] = true;
            if self.spec.step_ok(q, &self.mapping) {
                match self.search(depth + 1) {
                    SearchStatus::Found => {
                        self.cand_buf.truncate(cand_base);
                        self.tried_buf.truncate(tried_base);
                        return SearchStatus::Found;
                    }
                    SearchStatus::BudgetExceeded => {
                        self.target_used[t] = false;
                        self.mapping[q] = usize::MAX;
                        self.cand_buf.truncate(cand_base);
                        self.tried_buf.truncate(tried_base);
                        return SearchStatus::BudgetExceeded;
                    }
                    SearchStatus::Exhausted => {}
                }
            }
            self.target_used[t] = false;
            self.mapping[q] = usize::MAX;
            if let Some(rem) = &mut self.budget {
                if *rem == 0 {
                    self.cand_buf.truncate(cand_base);
                    self.tried_buf.truncate(tried_base);
                    return SearchStatus::BudgetExceeded;
                }
                *rem -= 1;
            }
            if use_orbits {
                self.tried_buf.push(orb);
            }
        }
        self.cand_buf.truncate(cand_base);
        self.tried_buf.truncate(tried_base);
        SearchStatus::Exhausted
    }
}
