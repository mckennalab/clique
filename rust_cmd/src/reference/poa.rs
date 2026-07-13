#![allow(dead_code)]
//! Partial-order alignment (POA) graph over a reference panel — Milestone 1:
//! construction only.
//!
//! A POA graph is a DAG of single-base nodes. Each node records the set of
//! references passing through it; shared bases become nodes supported by all
//! references (the backbone), and divergences branch. It is built progressively:
//! the first reference is a linear chain, and each subsequent reference is
//! aligned to the current graph (a Needleman–Wunsch generalised over the DAG)
//! and fused in — matches merge into existing nodes, substitutions add a sibling
//! node in the same column (an `aligned_to` ring), insertions add a new node, and
//! deletions leave a "skip" edge over the deleted node.
//!
//! This module implements construction and recovery only. Aligning a *read* to
//! the graph and classifying it by which branch it traverses is Milestone 2.

use std::collections::VecDeque;

use crate::reference::fasta_reference::ReferenceManager;

/// Virtual start / end nodes have this sentinel base (never a real DNA base).
const SENTINEL: u8 = 0;

// Simple integer scores. Match beats a substitution, and a substitution beats an
// indel pair, so near-identical equal-length references align ungapped.
const MATCH: i32 = 2;
const MISMATCH: i32 = -2;
const GAP: i32 = -3;

/// Half-width of the banded classification DP (in read positions). Comfortably
/// covers small indels plus error jitter for near-identical panels.
const CLASSIFY_BAND: usize = 25;

/// One node of the POA graph: a single base plus the references through it.
struct PoaNode {
    base: u8,
    supporting_refs: Vec<u16>,       // reference ids passing through this node (unique)
    out_edges: Vec<(usize, Vec<u16>)>, // (successor node id, refs traversing this edge)
    in_edges: Vec<usize>,            // predecessor node ids
    aligned_to: Vec<usize>,          // sibling nodes occupying the same column
}

impl PoaNode {
    fn new(base: u8, ref_id: Option<u16>) -> PoaNode {
        PoaNode {
            base,
            supporting_refs: ref_id.into_iter().collect(),
            out_edges: Vec::new(),
            in_edges: Vec::new(),
            aligned_to: Vec::new(),
        }
    }
}

/// A partial-order alignment graph built from a reference panel.
pub struct PoaGraph {
    nodes: Vec<PoaNode>,
    start: usize,
    end: usize,
    ref_names: Vec<Vec<u8>>,
    min_margin: usize,
}

/// The outcome of classifying a read against the graph.
#[derive(Debug, Clone, PartialEq)]
pub struct PoaClassification {
    /// Best-matching reference name.
    pub best: Vec<u8>,
    /// Branch-votes for the best reference.
    pub best_score: usize,
    /// Runner-up reference name (if the panel has >= 2 references).
    pub second: Option<Vec<u8>>,
    /// Branch-votes for the runner-up.
    pub second_score: usize,
    /// `best_score - second_score` — the confidence margin.
    pub margin: usize,
    /// Discriminating columns the read matched (carried a vote).
    pub informative_columns: usize,
    /// True when the margin is below the configured minimum (a near-tie call).
    pub ambiguous: bool,
}

/// How the sequence ends are treated in a graph alignment.
#[derive(Clone, Copy, PartialEq)]
enum AlignMode {
    /// Both the sequence and the graph are aligned end-to-end (leading/trailing
    /// gaps penalised). Used to fuse a full reference into the graph.
    Global,
    /// Free gaps at the graph ends: the read may cover only a substring of the
    /// reference (a "fitting" alignment — the read is fully consumed, but the
    /// uncovered graph prefix/suffix costs nothing). Used to classify reads.
    SemiGlobal,
}

/// One step of a sequence-to-graph alignment (from the DP traceback).
#[derive(Debug, Clone, PartialEq)]
enum AlignOp {
    /// The read base aligned to graph node `node` (a match if bases are equal,
    /// otherwise a substitution), arriving from graph node `from` — i.e. the read
    /// traversed the graph edge `from -> node` (used for edge-based voting).
    MatchOrSub { from: usize, node: usize, base: u8 },
    /// The read base is an insertion relative to the graph.
    Insert { base: u8 },
    /// The graph node is deleted (skipped) by the read.
    Delete { node: usize },
}

/// A traceback pointer for one DP cell.
#[derive(Clone, Copy)]
enum Back {
    None,
    Match(usize),  // came from predecessor topo-index (consumed a read base)
    Delete(usize), // came from predecessor topo-index (node deleted)
    Insert,        // came from the same node at j-1 (read base inserted)
    FromPred(usize), // zero-cost transition into a virtual node (END)
}

impl PoaGraph {
    /// Start a graph from the first reference as a linear chain.
    pub fn new(first_name: &[u8], first_seq: &[u8]) -> PoaGraph {
        let mut graph = PoaGraph {
            nodes: vec![PoaNode::new(SENTINEL, None), PoaNode::new(SENTINEL, None)],
            start: 0,
            end: 1,
            ref_names: vec![first_name.to_vec()],
            min_margin: 1,
        };
        let mut prev = graph.start;
        for &base in first_seq {
            let node = graph.push_node(base, 0);
            graph.add_edge(prev, node, 0);
            prev = node;
        }
        let end = graph.end;
        graph.add_edge(prev, end, 0);
        graph
    }

    /// Build a graph from a panel of `(name, sequence)` references.
    pub fn from_references(references: &[(Vec<u8>, Vec<u8>)]) -> Result<PoaGraph, String> {
        let (first_name, first_seq) = references
            .first()
            .ok_or_else(|| "cannot build a POA graph from an empty panel".to_string())?;
        let mut graph = PoaGraph::new(first_name, first_seq);
        for (name, seq) in &references[1..] {
            graph.add_reference(name, seq)?;
        }
        Ok(graph)
    }

    /// Build a graph from a loaded [`ReferenceManager`]. References are added in
    /// name order for determinism. `min_margin` is the smallest branch-vote gap
    /// for a confident (non-ambiguous) call.
    pub fn from_reference_manager(rm: &ReferenceManager, min_margin: usize) -> Result<PoaGraph, String> {
        let mut refs: Vec<(Vec<u8>, Vec<u8>)> = rm
            .references
            .values()
            .map(|r| (r.name.clone(), r.sequence.clone()))
            .collect();
        refs.sort_by(|a, b| a.0.cmp(&b.0));
        let mut graph = PoaGraph::from_references(&refs)?;
        graph.min_margin = min_margin.max(1);
        Ok(graph)
    }

    /// Whether the panel has any branch (discriminating) nodes to classify on.
    pub fn is_discriminable(&self) -> bool {
        self.branch_node_count() > 0
    }

    /// Classify a read by aligning it to the graph and voting, at each
    /// discriminating column the read matches, for the references whose branch it
    /// took. Best reference = most votes, with a top-2 margin confidence.
    pub fn classify_read(&self, read: &[u8]) -> PoaClassification {
        // Band the DP for near-full-length reads (their optimal path stays near
        // the depth diagonal). Short/partial reads can enter the graph anywhere,
        // so they use the exact full DP.
        let graph_len = self.node_count();
        // Band only genuinely near-full-length reads (>= 85% of the graph), whose
        // optimal path hugs the depth diagonal. A shorter fragment can enter the
        // graph anywhere, so its diagonal is unknown -> exact full DP.
        let band = if read.len() * 100 >= graph_len * 85 && read.len() <= graph_len * 3 / 2 {
            Some(CLASSIFY_BAND)
        } else {
            None
        };
        let ops = self.align_sequence(read, AlignMode::SemiGlobal, band).unwrap_or_default();
        let n_refs = self.n_refs();

        // Nodes the read truly matches (its base equals the node's base).
        let mut matched = vec![false; self.nodes.len()];
        for op in &ops {
            if let AlignOp::MatchOrSub { node, base, .. } = op {
                if self.nodes[*node].base == *base {
                    matched[*node] = true;
                }
            }
        }

        let mut scores = vec![0usize; n_refs];
        let mut informative = 0usize;
        for op in &ops {
            if let AlignOp::MatchOrSub { from, node, base } = op {
                // Edge-based vote: the read truly matches `node` and arrived via
                // edge `from -> node`. We vote on the EDGE's ref set (more specific
                // than the node's — it captures which branch the read took, so it
                // discriminates refs that merely share `node` and votes correctly
                // on a deletion's skip edge). Require the source to be reliable too
                // (the graph start, or another node the read truly matched) so a
                // transition out of an error/N-masked position casts no spurious
                // vote based on how the DP happened to align that position.
                let source_reliable = *from == self.start || matched[*from];
                if self.nodes[*node].base == *base && source_reliable {
                    if let Some((_, refs)) =
                        self.nodes[*from].out_edges.iter().find(|(t, _)| t == node)
                    {
                        if refs.len() < n_refs {
                            informative += 1;
                            for &r in refs {
                                scores[r as usize] += 1;
                            }
                        }
                    }
                }
            }
        }

        // Rank by votes desc, then reference id asc for a stable tie order.
        let mut ranked: Vec<(usize, usize)> = scores.iter().copied().enumerate().collect();
        ranked.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));

        let (best_id, best_score) = ranked[0];
        let (second, second_score) = match ranked.get(1) {
            Some(&(id, s)) => (Some(self.ref_names[id].clone()), s),
            None => (None, 0),
        };
        let margin = best_score - second_score;
        let ambiguous = n_refs >= 2 && margin < self.min_margin;

        PoaClassification {
            best: self.ref_names[best_id].clone(),
            best_score,
            second,
            second_score,
            margin,
            informative_columns: informative,
            ambiguous,
        }
    }

    /// Align a reference to the current graph and fuse it in.
    pub fn add_reference(&mut self, name: &[u8], seq: &[u8]) -> Result<(), String> {
        let ref_id = self.ref_names.len() as u16;
        // Construction aligns full references end-to-end and is one-time, so it
        // uses the exact (unbanded) global DP.
        let ops = self.align_sequence(seq, AlignMode::Global, None)?;
        self.incorporate(&ops, ref_id, seq);
        self.ref_names.push(name.to_vec());
        Ok(())
    }

    // ---- accessors ----------------------------------------------------------

    /// Number of references fused into the graph.
    pub fn n_refs(&self) -> usize {
        self.ref_names.len()
    }

    /// Number of real (non-virtual) base nodes.
    pub fn node_count(&self) -> usize {
        self.nodes.len() - 2
    }

    /// Nodes supported by every reference (the shared backbone).
    pub fn backbone_node_count(&self) -> usize {
        let n = self.n_refs();
        self.real_nodes().filter(|&i| self.nodes[i].supporting_refs.len() == n).count()
    }

    /// Nodes supported by only some references (the divergences).
    pub fn branch_node_count(&self) -> usize {
        let n = self.n_refs();
        self.real_nodes().filter(|&i| self.nodes[i].supporting_refs.len() < n).count()
    }

    /// Recover a reference's sequence by walking its supporting path from the
    /// start. Returns `None` if the path is broken (should not happen).
    pub fn reference_path(&self, ref_id: u16) -> Option<Vec<u8>> {
        let mut seq = Vec::new();
        let mut cur = self.start;
        loop {
            // Follow the out-edge this reference traverses (edge-level support
            // uniquely determines the path even when refs share nodes).
            let next = self.nodes[cur]
                .out_edges
                .iter()
                .find(|(to, refs)| *to != self.end && refs.contains(&ref_id))
                .map(|(to, _)| *to);
            match next {
                Some(s) => {
                    seq.push(self.nodes[s].base);
                    cur = s;
                }
                None => break,
            }
        }
        Some(seq)
    }

    // ---- construction internals --------------------------------------------

    fn real_nodes(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.nodes.len()).filter(move |&i| i != self.start && i != self.end)
    }

    fn push_node(&mut self, base: u8, ref_id: u16) -> usize {
        self.nodes.push(PoaNode::new(base, Some(ref_id)));
        self.nodes.len() - 1
    }

    fn add_edge(&mut self, from: usize, to: usize, ref_id: u16) {
        match self.nodes[from].out_edges.iter_mut().find(|(t, _)| *t == to) {
            Some((_, refs)) => {
                if !refs.contains(&ref_id) {
                    refs.push(ref_id);
                }
            }
            None => self.nodes[from].out_edges.push((to, vec![ref_id])),
        }
        if !self.nodes[to].in_edges.contains(&from) {
            self.nodes[to].in_edges.push(from);
        }
    }

    fn add_ref_support(&mut self, node: usize, ref_id: u16) {
        if !self.nodes[node].supporting_refs.contains(&ref_id) {
            self.nodes[node].supporting_refs.push(ref_id);
        }
    }

    /// Kahn topological sort; errors if the graph is not a DAG.
    fn topological_order(&self) -> Result<Vec<usize>, String> {
        let mut in_degree: Vec<usize> = self.nodes.iter().map(|n| n.in_edges.len()).collect();
        let mut queue: VecDeque<usize> =
            (0..self.nodes.len()).filter(|&i| in_degree[i] == 0).collect();
        let mut order = Vec::with_capacity(self.nodes.len());
        while let Some(u) = queue.pop_front() {
            order.push(u);
            for &(v, _) in &self.nodes[u].out_edges {
                in_degree[v] -= 1;
                if in_degree[v] == 0 {
                    queue.push_back(v);
                }
            }
        }
        if order.len() != self.nodes.len() {
            return Err("POA graph contains a cycle (not a DAG)".to_string());
        }
        Ok(order)
    }

    /// Needleman–Wunsch of `seq` against the graph (over the DAG), returning the
    /// traceback as a list of alignment operations in read order.
    fn align_sequence(
        &self,
        seq: &[u8],
        mode: AlignMode,
        band: Option<usize>,
    ) -> Result<Vec<AlignOp>, String> {
        let topo = self.topological_order()?;
        let n = topo.len();
        let m = seq.len();

        let mut pos_of = vec![usize::MAX; self.nodes.len()];
        for (ti, &nid) in topo.iter().enumerate() {
            pos_of[nid] = ti;
        }

        // Banding: restrict each node's read-position range to a window around
        // the topological-depth diagonal. `depth[v]` is the longest path from the
        // start; a node at depth d maps to read position ~ d*m/D. A read that
        // aligns near-diagonally (a near-full-length read) keeps its optimal path
        // inside the band; short/partial reads must use the full DP (band = None).
        // `range[ti]` is the inclusive `[lo, hi]` of read positions to fill.
        let range: Vec<(usize, usize)> = match band {
            None => vec![(0, m); n],
            Some(w) => {
                let mut depth = vec![0usize; self.nodes.len()];
                for &u in &topo {
                    for &(v, _) in &self.nodes[u].out_edges {
                        depth[v] = depth[v].max(depth[u] + 1);
                    }
                }
                let d_end = depth[self.end].max(1);
                topo.iter()
                    .map(|&nid| {
                        if nid == self.start {
                            (0, m) // the start row is a universal predecessor
                        } else {
                            let center = depth[nid] * m / d_end;
                            (center.saturating_sub(w), (center + w).min(m))
                        }
                    })
                    .collect()
            }
        };

        let neg = i32::MIN / 4;
        // Banded j-ranges prune the columns filled per row; unwritten (out-of-band)
        // cells stay `neg`, so predecessor reads outside the band read as neg. A
        // compact flat layout was tried but was measurably SLOWER on short
        // amplicons (the per-access bound check outweighs the allocation savings);
        // it would only pay off for long references, so we keep the simple layout.
        let mut score = vec![vec![neg; m + 1]; n];
        let mut back = vec![vec![Back::None; m + 1]; n];

        for ti in 0..n {
            let nid = topo[ti];

            if nid == self.start {
                score[ti][0] = 0;
                for j in 1..=m {
                    score[ti][j] = score[ti][j - 1] + GAP;
                    back[ti][j] = Back::Insert;
                }
                continue;
            }

            let preds: Vec<usize> = self.nodes[nid].in_edges.iter().map(|&u| pos_of[u]).collect();

            if nid == self.end {
                let (lo, hi) = range[ti];
                for j in lo..=hi {
                    let mut best = neg;
                    let mut b = Back::None;
                    for &pt in &preds {
                        if score[pt][j] > best {
                            best = score[pt][j];
                            b = Back::FromPred(pt);
                        }
                    }
                    if j >= 1 && score[ti][j - 1] + GAP > best {
                        best = score[ti][j - 1] + GAP;
                        b = Back::Insert;
                    }
                    score[ti][j] = best;
                    back[ti][j] = b;
                }
                continue;
            }

            // real node
            let base = self.nodes[nid].base;
            let (lo, hi) = range[ti];
            for j in lo..=hi {
                // Semi-global: a read may enter the graph at any node for free
                // (the uncovered graph prefix costs nothing).
                if j == 0 && mode == AlignMode::SemiGlobal {
                    score[ti][0] = 0;
                    back[ti][0] = Back::None;
                    continue;
                }

                let mut best = neg;
                let mut b = Back::None;

                // delete this node (gap in the read)
                for &pt in &preds {
                    let v = score[pt][j] + GAP;
                    if v > best {
                        best = v;
                        b = Back::Delete(pt);
                    }
                }

                if j >= 1 {
                    // match / substitution
                    let s = if base == seq[j - 1] { MATCH } else { MISMATCH };
                    for &pt in &preds {
                        let v = score[pt][j - 1] + s;
                        if v > best {
                            best = v;
                            b = Back::Match(pt);
                        }
                    }
                    // insertion (gap in the graph)
                    let v = score[ti][j - 1] + GAP;
                    if v > best {
                        best = v;
                        b = Back::Insert;
                    }
                }

                score[ti][j] = best;
                back[ti][j] = b;
            }
        }

        // Choose the cell to trace back from.
        // - Global: (END, m) — both fully consumed.
        // - Semi-global: the best node at column m (read fully consumed; the
        //   uncovered graph suffix is free), tracing back until the read is
        //   exhausted (j == 0), leaving the graph prefix free.
        let start_pos = pos_of[self.start];
        let (mut ti, mut j) = match mode {
            AlignMode::Global => (pos_of[self.end], m),
            AlignMode::SemiGlobal => {
                let mut best_ti = pos_of[self.end];
                let mut best = score[best_ti][m];
                for cand in 0..n {
                    if topo[cand] == self.start {
                        continue;
                    }
                    if score[cand][m] > best {
                        best = score[cand][m];
                        best_ti = cand;
                    }
                }
                (best_ti, m)
            }
        };
        let mut ops = Vec::new();
        loop {
            let done = match mode {
                AlignMode::Global => ti == start_pos && j == 0,
                AlignMode::SemiGlobal => j == 0,
            };
            if done {
                break;
            }
            let nid = topo[ti];
            match back[ti][j] {
                Back::FromPred(pt) => ti = pt,
                Back::Insert => {
                    ops.push(AlignOp::Insert { base: seq[j - 1] });
                    j -= 1;
                }
                Back::Delete(pt) => {
                    ops.push(AlignOp::Delete { node: nid });
                    ti = pt;
                }
                Back::Match(pt) => {
                    ops.push(AlignOp::MatchOrSub { from: topo[pt], node: nid, base: seq[j - 1] });
                    ti = pt;
                    j -= 1;
                }
                Back::None => {
                    return Err("POA traceback reached a dead end".to_string());
                }
            }
        }
        ops.reverse();
        Ok(ops)
    }

    /// Fuse an aligned reference into the graph.
    fn incorporate(&mut self, ops: &[AlignOp], ref_id: u16, _seq: &[u8]) {
        let mut prev = self.start;
        for op in ops {
            match op {
                AlignOp::MatchOrSub { node, base, .. } => {
                    let v = *node;
                    if self.nodes[v].base == *base {
                        // true match: merge into the existing node
                        self.add_ref_support(v, ref_id);
                        self.add_edge(prev, v, ref_id);
                        prev = v;
                    } else {
                        // substitution: reuse a same-column sibling with this base,
                        // else create one and link it into the column ring.
                        let sibling = self.nodes[v]
                            .aligned_to
                            .iter()
                            .copied()
                            .find(|&w| self.nodes[w].base == *base);
                        let w = match sibling {
                            Some(w) => {
                                self.add_ref_support(w, ref_id);
                                w
                            }
                            None => {
                                let w = self.push_node(*base, ref_id);
                                let mut ring = self.nodes[v].aligned_to.clone();
                                ring.push(v);
                                for x in ring {
                                    self.nodes[w].aligned_to.push(x);
                                    self.nodes[x].aligned_to.push(w);
                                }
                                w
                            }
                        };
                        self.add_edge(prev, w, ref_id);
                        prev = w;
                    }
                }
                AlignOp::Insert { base } => {
                    let w = self.push_node(*base, ref_id);
                    self.add_edge(prev, w, ref_id);
                    prev = w;
                }
                AlignOp::Delete { .. } => {
                    // node skipped by this reference: leave prev where it is so the
                    // next edge jumps over the deleted node (a deletion branch).
                }
            }
        }
        let end = self.end;
        self.add_edge(prev, end, ref_id);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const UNEDITED: &str = "GCCTCCACGGCCACTAGTATTATGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACGTATTAGTCATCGCTATTACCATGTACTCATCCTGTCATCTTAGCTAAGATGACAGGTAATTCGAATTTAAATCGGATCCGCGGCC";
    const LEFT: &str = "GCCTCCACGGCCACTAGTATTATGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACGTATTAGTCATCGCTATTACCATGCGATGCTCCTGTCGTCTTAGCTAAGATGACAGGTAATTCGAATTTAAATCGGATCCGCGGCC";
    const RIGHT: &str = "GCCTCCACGGCCACTAGTATTATGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACGTATTAGTCATCGCTATTACCATGCATACTTCCTGTCATCTTAGCTAAGACGACAGGTAATTCGAATTTAAATCGGATCCGCGGCC";

    fn rnf2_panel() -> Vec<(Vec<u8>, Vec<u8>)> {
        vec![
            (b"unedited".to_vec(), UNEDITED.as_bytes().to_vec()),
            (b"left".to_vec(), LEFT.as_bytes().to_vec()),
            (b"right".to_vec(), RIGHT.as_bytes().to_vec()),
        ]
    }

    #[test]
    fn test_single_reference_chain() {
        let g = PoaGraph::new(b"a", b"ACGTACGT");
        assert_eq!(g.node_count(), 8);
        assert_eq!(g.n_refs(), 1);
        assert_eq!(g.reference_path(0).unwrap(), b"ACGTACGT".to_vec());
        assert!(g.topological_order().is_ok());
    }

    #[test]
    fn test_identical_references_fully_merge() {
        let refs = vec![(b"a".to_vec(), b"ACGTACGT".to_vec()), (b"b".to_vec(), b"ACGTACGT".to_vec())];
        let g = PoaGraph::from_references(&refs).unwrap();
        // Identical sequences share every node: no new nodes, all backbone.
        assert_eq!(g.node_count(), 8);
        assert_eq!(g.backbone_node_count(), 8);
        assert_eq!(g.branch_node_count(), 0);
        assert_eq!(g.reference_path(0).unwrap(), b"ACGTACGT".to_vec());
        assert_eq!(g.reference_path(1).unwrap(), b"ACGTACGT".to_vec());
    }

    #[test]
    fn test_single_substitution_makes_a_column_ring() {
        // Differ at position 2 only (A vs G).
        let refs = vec![(b"a".to_vec(), b"ACAT".to_vec()), (b"b".to_vec(), b"ACGT".to_vec())];
        let g = PoaGraph::from_references(&refs).unwrap();
        assert!(g.topological_order().is_ok());
        // Shared A,C,_,T positions are backbone (3 nodes); the middle column has
        // two nodes (A for a, G for b) -> 2 branch nodes.
        assert_eq!(g.backbone_node_count(), 3);
        assert_eq!(g.branch_node_count(), 2);
        assert_eq!(g.node_count(), 5);
        assert_eq!(g.reference_path(0).unwrap(), b"ACAT".to_vec());
        assert_eq!(g.reference_path(1).unwrap(), b"ACGT".to_vec());
    }

    #[test]
    fn test_insertion_adds_a_branch_node() {
        // b has an extra base inserted in the middle.
        let refs = vec![(b"a".to_vec(), b"ACGT".to_vec()), (b"b".to_vec(), b"ACAGT".to_vec())];
        let g = PoaGraph::from_references(&refs).unwrap();
        assert!(g.topological_order().is_ok());
        assert_eq!(g.reference_path(0).unwrap(), b"ACGT".to_vec());
        assert_eq!(g.reference_path(1).unwrap(), b"ACAGT".to_vec());
        // The inserted base is supported only by b.
        assert!(g.branch_node_count() >= 1);
    }

    #[test]
    fn test_deletion_leaves_a_skip_edge() {
        // b is missing a base relative to a.
        let refs = vec![(b"a".to_vec(), b"ACGT".to_vec()), (b"b".to_vec(), b"AGT".to_vec())];
        let g = PoaGraph::from_references(&refs).unwrap();
        assert!(g.topological_order().is_ok());
        assert_eq!(g.reference_path(0).unwrap(), b"ACGT".to_vec());
        assert_eq!(g.reference_path(1).unwrap(), b"AGT".to_vec());
    }

    #[test]
    fn test_rnf2_panel_construction() {
        let g = PoaGraph::from_references(&rnf2_panel()).unwrap();
        assert_eq!(g.n_refs(), 3);
        assert!(g.topological_order().is_ok(), "graph must be a DAG");
        eprintln!(
            "RNF2 POA: {} nodes ({} backbone, {} branch) over {} refs",
            g.node_count(), g.backbone_node_count(), g.branch_node_count(), g.n_refs()
        );

        // Every reference is recovered exactly from its supporting path.
        assert_eq!(g.reference_path(0).unwrap(), UNEDITED.as_bytes().to_vec());
        assert_eq!(g.reference_path(1).unwrap(), LEFT.as_bytes().to_vec());
        assert_eq!(g.reference_path(2).unwrap(), RIGHT.as_bytes().to_vec());

        // The panel is ~156 bp and differs by only a few bases, so the backbone
        // dominates and there are only a handful of branch nodes.
        assert!(g.backbone_node_count() > 140, "backbone too small: {}", g.backbone_node_count());
        assert!(g.branch_node_count() < 20, "too many branch nodes: {}", g.branch_node_count());
    }

    #[test]
    fn test_deletion_only_distinction_via_edge_voting() {
        // `del` is `has` minus one base. The two share every node, so node-level
        // voting cannot tell them apart; edge-level voting uses the skip edge that
        // `del` takes across the deleted base.
        let refs = vec![
            (b"has".to_vec(), b"ACGTACGT".to_vec()),
            (b"del".to_vec(), b"ACGACGT".to_vec()), // missing the T at index 3
        ];
        let g = PoaGraph::from_references(&refs).unwrap();
        let has = g.classify_read(b"ACGTACGT");
        assert_eq!(has.best, b"has".to_vec());
        assert!(!has.ambiguous, "the base-carrying read should be confident");
        let del = g.classify_read(b"ACGACGT");
        assert_eq!(del.best, b"del".to_vec(), "deletion read must classify to the deleting ref");
        assert!(!del.ambiguous, "the skip edge should make the deletion read confident");
    }

    #[test]
    fn test_classify_each_exact_reference_to_itself() {
        let g = PoaGraph::from_references(&rnf2_panel()).unwrap();
        for (id, seq) in [(0u16, UNEDITED), (1, LEFT), (2, RIGHT)] {
            let c = g.classify_read(seq.as_bytes());
            assert_eq!(c.best, g.ref_names[id as usize], "ref {} misclassified", id);
            assert!(!c.ambiguous, "exact reference {} should be unambiguous (margin {})", id, c.margin);
            assert!(c.margin >= 1);
        }
    }

    #[test]
    fn test_backbone_only_read_is_ambiguous() {
        // Mask the columns where references differ; the read then carries no
        // branch signal, so every reference ties.
        let g = PoaGraph::from_references(&rnf2_panel()).unwrap();
        // The RNF2 differences sit around cols 95-100 and 108/121; N-mask them.
        let mut read = UNEDITED.as_bytes().to_vec();
        for &p in &[95usize, 96, 97, 98, 99, 100, 108, 121] {
            read[p] = b'N';
        }
        let c = g.classify_read(&read);
        assert!(c.ambiguous, "backbone-only read should be ambiguous (margin {})", c.margin);
    }

    #[test]
    fn test_partial_read_classifies_via_semi_global() {
        // A truncated read covering only the middle of the reference (including
        // the discriminating region ~95-121) must still classify correctly —
        // this is the semi-global (fitting) win: the uncovered reference ends
        // are free rather than penalised as end gaps.
        let g = PoaGraph::from_references(&rnf2_panel()).unwrap();
        for (id, seq) in [(0u16, UNEDITED), (1, LEFT), (2, RIGHT)] {
            let fragment = &seq.as_bytes()[60..145]; // covers the divergences, not the ends
            let c = g.classify_read(fragment);
            assert_eq!(c.best, g.ref_names[id as usize],
                "partial read of ref {} misclassified as {:?}", id, String::from_utf8_lossy(&c.best));
            assert!(!c.ambiguous, "partial ref {} should still be a confident call (margin {})", id, c.margin);
        }
    }

    #[test]
    fn test_offset_partial_read_is_not_penalised() {
        // A short fragment from near one end (positions 90-130) still lands on
        // its true substring under semi-global and calls the right reference.
        let g = PoaGraph::from_references(&rnf2_panel()).unwrap();
        let fragment = &LEFT.as_bytes()[90..130];
        let c = g.classify_read(fragment);
        assert_eq!(c.best, b"left".to_vec());
    }

    #[test]
    fn test_backbone_error_does_not_flip_the_call() {
        let g = PoaGraph::from_references(&rnf2_panel()).unwrap();
        let mut read = LEFT.as_bytes().to_vec();
        // corrupt a backbone base near the start (position 5, well away from the
        // discriminating region) -> must still call "left".
        read[5] = if read[5] == b'A' { b'T' } else { b'A' };
        let c = g.classify_read(&read);
        assert_eq!(c.best, b"left".to_vec(), "backbone error flipped the call");
        assert!(!c.ambiguous);
    }

    /// Validation against real reads: classify each read and compare to clique's
    /// alignment-based reference call (the BAM RNAME) from a TSV
    /// (`read_seq \t clique_ref`) in `CLIQUE_POA_TSV`. Ignored by default; run:
    ///   CLIQUE_POA_TSV=reads_ref.tsv cargo test --bin clique \
    ///     reference::poa::tests::validate_on_real_reads -- --ignored --nocapture
    #[test]
    #[ignore]
    fn validate_on_real_reads() {
        use std::io::BufRead;
        let path = std::env::var("CLIQUE_POA_TSV").expect("set CLIQUE_POA_TSV");
        let g = PoaGraph::from_references(&rnf2_panel()).unwrap();
        eprintln!(
            "POA graph: {} nodes ({} backbone, {} branch)",
            g.node_count(), g.backbone_node_count(), g.branch_node_count()
        );

        let file = std::fs::File::open(&path).unwrap();
        let (mut total, mut agree, mut ambiguous, mut disagree) = (0, 0, 0, 0);
        for line in std::io::BufReader::new(file).lines() {
            let line = line.unwrap();
            let mut it = line.split('\t');
            let (seq, clique_ref) = (it.next().unwrap(), it.next().unwrap().as_bytes());
            let c = g.classify_read(seq.as_bytes());
            total += 1;
            if c.ambiguous {
                ambiguous += 1;
            } else if c.best == clique_ref {
                agree += 1;
            } else {
                disagree += 1;
            }
        }
        let confident = agree + disagree;
        let conc = if confident == 0 { 1.0 } else { agree as f64 / confident as f64 };
        eprintln!(
            "reads={} confident-agree={} ({:.1}%) ambiguous={} disagree={} | confident-concordance={:.1}%",
            total, agree, 100.0 * agree as f64 / total as f64,
            ambiguous, disagree, 100.0 * conc,
        );
        assert!(conc > 0.95, "POA confident-concordance with the aligner too low: {:.3}", conc);
    }

    #[test]
    fn test_recovery_holds_under_any_build_order() {
        // Progressive POA is order-DEPENDENT (node counts may differ by build
        // order), but the invariant that matters holds regardless: every
        // reference is recovered exactly from its edge-supported path.
        let mut panel = rnf2_panel();
        panel.swap(0, 2); // build in a different order
        let g = PoaGraph::from_references(&panel).unwrap();
        assert!(g.topological_order().is_ok());
        // panel[0] is now "right", panel[1] "left", panel[2] "unedited".
        assert_eq!(g.reference_path(0).unwrap(), RIGHT.as_bytes().to_vec());
        assert_eq!(g.reference_path(1).unwrap(), LEFT.as_bytes().to_vec());
        assert_eq!(g.reference_path(2).unwrap(), UNEDITED.as_bytes().to_vec());
        // The backbone still dominates.
        assert!(g.backbone_node_count() > 140);
    }
}
