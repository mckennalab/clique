//! IDF-weighted k-mer reference router.
//!
//! The default router keeps only k-mers unique to a single reference and votes;
//! on near-identical panels that index empties out, and raw vote counts carry a
//! size bias. This module instead weights every panel k-mer by its
//! **inverse document frequency**:
//!
//! ```text
//! weight(kmer) = ln(N / df(kmer))
//! ```
//!
//! where `N` is the number of references and `df` is how many contain the k-mer.
//! A k-mer in one reference gets `ln N` (maximally discriminating); one shared by
//! all references gets `ln 1 = 0` and contributes nothing — so the shared backbone
//! is automatically ignored and the score is driven by the k-mers that actually
//! distinguish references. A read is scored by summing the weights of its distinct
//! k-mers against each reference; the best reference is returned with a top-2
//! **margin ratio** confidence and an `ambiguous` flag. This runs at the routing
//! stage (no alignment), so it stays fast on large panels.

use std::collections::{HashMap, HashSet};

use crate::reference::fasta_reference::ReferenceManager;

/// Per-k-mer index entry: its IDF weight and the references that contain it.
struct KmerEntry {
    weight: f64,
    refs: Vec<Vec<u8>>, // reference names
}

/// An IDF-weighted k-mer index over a reference panel.
pub struct IdfIndex {
    kmers: HashMap<Vec<u8>, KmerEntry>,
    kmer_size: usize,
    kmer_skip: usize,
    n_refs: usize,
    min_margin_ratio: f64,
}

/// The outcome of scoring one read.
#[derive(Debug, Clone, PartialEq)]
pub struct IdfScore {
    /// Best reference name, or `None` when no read k-mer matched a discriminating
    /// panel k-mer (the caller should fall back to a full search).
    pub best: Option<Vec<u8>>,
    /// Summed IDF weight for the best reference.
    pub best_score: f64,
    /// Runner-up reference name (if any).
    pub second: Option<Vec<u8>>,
    /// Summed IDF weight for the runner-up.
    pub second_score: f64,
    /// `(best - second) / best`, in `[0, 1]` (0 when `best_score == 0`).
    pub margin_ratio: f64,
    /// Distinct read k-mers that hit a discriminating (weight > 0) panel k-mer.
    pub informative_kmers: usize,
    /// True when there is no discriminating signal or the margin is below the floor.
    pub ambiguous: bool,
}

impl IdfIndex {
    /// Build an IDF index from a loaded [`ReferenceManager`]. `min_margin_ratio`
    /// is the smallest top-2 relative gap for a confident (non-ambiguous) call.
    pub fn from_reference_manager(rm: &ReferenceManager, min_margin_ratio: f64) -> IdfIndex {
        let n_refs = rm.references.len();

        // Collect, per k-mer, the set of references that contain it.
        let mut kmer_to_refs: HashMap<Vec<u8>, Vec<Vec<u8>>> = HashMap::new();
        for reference in rm.references.values() {
            let distinct: HashSet<Vec<u8>> =
                ReferenceManager::sequence_to_kmers(&reference.sequence, &rm.kmer_size, &rm.kmer_skip)
                    .into_iter()
                    .map(|(kmer, _count)| kmer)
                    .collect();
            for kmer in distinct {
                kmer_to_refs.entry(kmer).or_default().push(reference.name.clone());
            }
        }

        // Convert document frequencies into IDF weights.
        let n = n_refs as f64;
        let kmers = kmer_to_refs
            .into_iter()
            .map(|(kmer, refs)| {
                let df = refs.len() as f64;
                let weight = (n / df).ln(); // 0 when df == N (shared by all)
                (kmer, KmerEntry { weight, refs })
            })
            .collect();

        IdfIndex {
            kmers,
            kmer_size: rm.kmer_size,
            kmer_skip: rm.kmer_skip,
            n_refs,
            min_margin_ratio,
        }
    }

    /// Number of references the index was built over.
    pub fn n_refs(&self) -> usize {
        self.n_refs
    }

    /// Whether any panel k-mer is discriminating (weight > 0). If not, IDF cannot
    /// route (all references share every k-mer) and the caller should fall back.
    pub fn has_discriminating_kmers(&self) -> bool {
        self.kmers.values().any(|e| e.weight > 0.0)
    }

    /// Score a read: sum IDF weights over its distinct k-mers per reference and
    /// return the best reference with a top-2 margin-ratio confidence.
    pub fn score(&self, read: &[u8]) -> IdfScore {
        let distinct: HashSet<Vec<u8>> =
            ReferenceManager::sequence_to_kmers(&read.to_vec(), &self.kmer_size, &self.kmer_skip)
                .into_iter()
                .map(|(kmer, _count)| kmer)
                .collect();

        let mut scores: HashMap<Vec<u8>, f64> = HashMap::new();
        let mut informative = 0usize;
        for kmer in &distinct {
            if let Some(entry) = self.kmers.get(kmer) {
                if entry.weight > 0.0 {
                    informative += 1;
                    for name in &entry.refs {
                        *scores.entry(name.clone()).or_insert(0.0) += entry.weight;
                    }
                }
            }
        }

        // Rank by score desc, then name asc for a stable tie order.
        let mut ranked: Vec<(Vec<u8>, f64)> = scores.into_iter().collect();
        ranked.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap().then(a.0.cmp(&b.0)));

        let (best, best_score) = match ranked.first() {
            Some((n, s)) => (Some(n.clone()), *s),
            None => (None, 0.0),
        };
        let (second, second_score) = match ranked.get(1) {
            Some((n, s)) => (Some(n.clone()), *s),
            None => (None, 0.0),
        };

        let margin_ratio = if best_score > 0.0 {
            (best_score - second_score) / best_score
        } else {
            0.0
        };
        let ambiguous = best_score <= 0.0 || margin_ratio < self.min_margin_ratio;

        IdfScore {
            best,
            best_score,
            second,
            second_score,
            margin_ratio,
            informative_kmers: informative,
            ambiguous,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::fasta_reference::{Reference, ReferenceManager};

    fn manager(refs: Vec<(&str, &str)>, k: usize, skip: usize) -> ReferenceManager<'static, 'static, 'static> {
        let structs: Vec<Reference> = refs
            .iter()
            .map(|(name, seq)| Reference {
                sequence: seq.as_bytes().to_vec(),
                name: name.as_bytes().to_vec(),
                suffix_table: ReferenceManager::find_seeds(&seq.as_bytes().to_vec(), k),
            })
            .collect();
        ReferenceManager::from_fasta_vec(structs, k, skip)
    }

    #[test]
    fn test_backbone_kmers_get_zero_weight() {
        // Two references sharing a long backbone, differing only in the middle.
        let a = "AAAAAAAAAACGACAAAAAAAAAA";
        let b = "AAAAAAAAAATGTCAAAAAAAAAA";
        let idx = IdfIndex::from_reference_manager(&manager(vec![("a", a), ("b", b)], 6, 1), 0.05);
        assert!(idx.has_discriminating_kmers());
        // A backbone k-mer present in both refs must have weight 0.
        let backbone = b"AAAAAA".to_vec();
        assert!(idx.kmers.get(&backbone).map_or(true, |e| e.weight == 0.0));
    }

    #[test]
    fn test_scores_read_to_its_own_reference() {
        let a = "AAAAAAAAAACGACAAAAAAAAAA";
        let b = "AAAAAAAAAATGTCAAAAAAAAAA";
        let idx = IdfIndex::from_reference_manager(&manager(vec![("a", a), ("b", b)], 6, 1), 0.05);

        let sa = idx.score(a.as_bytes());
        assert_eq!(sa.best.as_deref(), Some(b"a".as_ref()));
        assert!(sa.best_score > sa.second_score);
        assert!(!sa.ambiguous, "an exact reference should be a confident call");

        let sb = idx.score(b.as_bytes());
        assert_eq!(sb.best.as_deref(), Some(b"b".as_ref()));
    }

    #[test]
    fn test_backbone_only_read_is_ambiguous_or_unrouted() {
        // A read that is pure backbone carries no discriminating k-mers.
        let a = "AAAAAAAAAACGACAAAAAAAAAA";
        let b = "AAAAAAAAAATGTCAAAAAAAAAA";
        let idx = IdfIndex::from_reference_manager(&manager(vec![("a", a), ("b", b)], 6, 1), 0.05);
        let s = idx.score(b"AAAAAAAAAAAAAA");
        assert!(s.ambiguous);
        assert!(s.best.is_none() || s.best_score == 0.0);
    }

    #[test]
    fn test_three_reference_idf_weights() {
        // k-mer in 1 of 3 refs -> ln(3); in 2 of 3 -> ln(1.5); in 3 -> 0.
        let refs = vec![
            ("r1", "GGGGGGGGAAAAAAAA"),
            ("r2", "GGGGGGGGAAAAAAAA"), // identical to r1 -> shares everything
            ("r3", "CCCCCCCCTTTTTTTT"), // fully distinct
        ];
        let idx = IdfIndex::from_reference_manager(&manager(refs, 8, 1), 0.05);
        // A k-mer only in r3 should weight ln(3); one shared by r1+r2 weights ln(1.5).
        let r3_kmer = b"CCCCCCCC".to_vec();
        assert!((idx.kmers[&r3_kmer].weight - (3.0f64).ln()).abs() < 1e-9);
        let shared_kmer = b"GGGGGGGG".to_vec();
        assert!((idx.kmers[&shared_kmer].weight - (3.0f64 / 2.0).ln()).abs() < 1e-9);
    }

    /// Validation on the diverse 180-guide library (the regime IDF targets):
    /// build the index from the real reference FASTA, route the real reads, and
    /// report the confident-routing rate. Ignored by default (streams a 9 MB
    /// FASTQ); run explicitly:
    ///   cargo test --bin clique reference::idf::tests::validate_on_large_library \
    ///     -- --ignored --nocapture
    #[test]
    #[ignore]
    fn validate_on_large_library() {
        use crate::read_strategies::read_set::ReadIterator;
        use std::path::PathBuf;

        let rm = ReferenceManager::from_fa_file(
            &"test_data/18guide1_pcr_sequence.fasta".to_string(),
            8,
            4,
        );
        let idx = IdfIndex::from_reference_manager(&rm, 0.05);
        eprintln!("panel: {} references, discriminating k-mers: {}", idx.n_refs(), idx.has_discriminating_kmers());
        assert!(idx.has_discriminating_kmers());

        let reads = ReadIterator::new(
            PathBuf::from("test_data/PAM_TWIST_1_018_S20_merged_001.fastq.gz"),
            None,
            None,
            None,
        );

        let (mut total, mut routed, mut confident, mut unrouted) = (0, 0, 0, 0);
        for read in reads.take(5000) {
            total += 1;
            let sc = idx.score(&read.read_one.seq().to_vec());
            match sc.best {
                None => unrouted += 1,
                Some(_) => {
                    routed += 1;
                    if !sc.ambiguous {
                        confident += 1;
                    }
                }
            }
        }
        eprintln!(
            "reads={} routed={} ({:.1}%) confident={} ({:.1}%) unrouted={}",
            total,
            routed, 100.0 * routed as f64 / total as f64,
            confident, 100.0 * confident as f64 / total as f64,
            unrouted,
        );
        assert!(routed as f64 / total as f64 > 0.8, "IDF routed too few reads on a diverse panel");
    }

    #[test]
    fn test_reads_sharing_discriminating_kmers_give_middling_margin() {
        // r1 shares its AAAACCCC half with r2; a read = r1 therefore gives r2 a
        // partial (non-zero) score, so the margin ratio is middling (~0.7 here),
        // and the threshold flips the call between confident and ambiguous.
        let refs = vec![("r1", "AAAACCCCGGGG"), ("r2", "AAAACCCCTTTT"), ("r3", "TTTTGGGGCCCC")];
        let lax = IdfIndex::from_reference_manager(&manager(refs.clone(), 4, 1), 0.5);
        let s = lax.score(refs[0].1.as_bytes());
        assert_eq!(s.best.as_deref(), Some(b"r1".as_ref()));
        assert!(s.second_score > 0.0, "runner-up must score via shared discriminating k-mers");
        assert!(s.margin_ratio > 0.5 && s.margin_ratio < 1.0, "margin {} not middling", s.margin_ratio);
        assert!(!s.ambiguous);

        let strict = IdfIndex::from_reference_manager(&manager(refs.clone(), 4, 1), 0.9);
        assert!(strict.score(refs[0].1.as_bytes()).ambiguous, "high threshold should flag the middling margin");
    }
}
