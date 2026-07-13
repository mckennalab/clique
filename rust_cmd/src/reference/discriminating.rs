//! Discriminating-position classifier for near-identical reference panels.
//!
//! When a panel of references differs by only a handful of bases (e.g. the RNF2
//! `unedited` / `left` / `right` amplicons, ~2 bp apart), the unique-k-mer router
//! empties out and a full-length alignment score is dominated by the shared
//! backbone, so the distinguishing bases are effectively ignored. This module
//! instead:
//!
//! 1. Projects every reference onto a common **anchor** coordinate frame (the
//!    longest reference), by aligning each reference to the anchor.
//! 2. Finds the **discriminating positions** — anchor columns where the panel
//!    references disagree (a differing base, or a gap vs. a base).
//! 3. Classifies a read by projecting it onto the anchor and comparing its bases
//!    at *only* those positions against each reference's signature, returning the
//!    best reference plus a top-2 **margin** confidence and an `ambiguous` flag.
//!
//! Scope/limits (prototype): discrimination is done in anchor coordinates, so
//! differences that are pure *insertions* relative to the anchor are not scored
//! (substitutions and deletions are). Reference bases are treated as concrete
//! ACGT; IUPAC degenerate codes are compared byte-exact.

use std::collections::BTreeMap;

use crate::alignment::scoring_functions::AffineScoring;
use crate::alignment_manager::align_two_strings;

const GAP: u8 = b'-';

/// A classifier built from a reference panel. Cheap to query once built.
pub struct DiscriminatingClassifier {
    /// Name of the anchor reference the coordinate frame is built on.
    pub anchor_name: String,
    /// The anchor reference sequence (reads are projected onto its coordinates).
    anchor_seq: Vec<u8>,
    /// Discriminating anchor columns (0-based indices into the anchor sequence).
    pub positions: Vec<usize>,
    /// Per-reference expected base at each discriminating position (`-` = gap).
    /// Each vector is aligned with `positions`.
    pub signatures: BTreeMap<String, Vec<u8>>,
    /// Panel reference names, in a stable order.
    reference_names: Vec<String>,
    /// Minimum top-2 score gap required to call a read unambiguously.
    min_margin: usize,
}

/// The outcome of classifying one read.
#[derive(Debug, Clone, PartialEq)]
pub struct Classification {
    /// Best-matching reference name.
    pub best: String,
    /// Matches at discriminating positions for the best reference.
    pub best_score: usize,
    /// Runner-up reference name (if the panel has >= 2 references).
    pub second: Option<String>,
    /// Matches for the runner-up.
    pub second_score: usize,
    /// `best_score - second_score` — the confidence margin.
    pub margin: usize,
    /// Discriminating positions the read actually covered with a real (ACGT) base.
    pub informative_positions: usize,
    /// Total number of discriminating positions in the panel.
    pub total_positions: usize,
    /// True when the margin is below `min_margin` (a forced/near-tie call).
    pub ambiguous: bool,
}

impl DiscriminatingClassifier {
    /// Build a classifier from `(name, sequence)` reference pairs. `min_margin`
    /// is the smallest top-2 score gap that counts as an unambiguous call
    /// (1 = require a strict winner).
    pub fn new(references: &[(String, Vec<u8>)], min_margin: usize) -> Result<Self, String> {
        DiscriminatingClassifier::new_excluding(references, min_margin, &[])
    }

    /// Like [`new`](Self::new) but drops the given anchor positions from the
    /// discriminating set. Use this to exclude CRISPR edit / target-window
    /// positions in a base-editing panel, where the columns that distinguish the
    /// reference *alleles* are the same columns the editor changes — otherwise an
    /// edited read is pulled toward whichever allele carries the edited base.
    pub fn new_excluding(
        references: &[(String, Vec<u8>)],
        min_margin: usize,
        exclude_positions: &[usize],
    ) -> Result<Self, String> {
        if references.is_empty() {
            return Err("cannot build a classifier from an empty reference panel".to_string());
        }
        let mut names: Vec<String> = references.iter().map(|(n, _)| n.clone()).collect();
        names.sort();
        names.dedup();
        if names.len() != references.len() {
            return Err("reference panel has duplicate names".to_string());
        }

        // Anchor = longest reference; ties broken by name for determinism.
        let anchor = references
            .iter()
            .max_by(|a, b| a.1.len().cmp(&b.1.len()).then(b.0.cmp(&a.0)))
            .unwrap();
        let anchor_name = anchor.0.clone();
        let anchor_seq = anchor.1.clone();
        let anchor_len = anchor_seq.len();

        let scoring = AffineScoring::default_dna();

        // Project each reference onto anchor coordinates: proj[p] is that
        // reference's base at anchor position p (or `-` for a deletion / no cover).
        let mut projections: BTreeMap<String, Vec<u8>> = BTreeMap::new();
        for (name, seq) in references {
            let proj = if name == &anchor_name {
                anchor_seq.clone()
            } else {
                project_onto_anchor(&anchor_seq, &anchor_name, seq, name, &scoring)
            };
            projections.insert(name.clone(), proj);
        }

        // A position is discriminating if the references show >1 distinct base
        // there (and it is not in the excluded set, e.g. an edit/target window).
        let excluded: std::collections::BTreeSet<usize> = exclude_positions.iter().copied().collect();
        let reference_names: Vec<String> = references.iter().map(|(n, _)| n.clone()).collect();
        let mut positions = Vec::new();
        for p in 0..anchor_len {
            if excluded.contains(&p) {
                continue;
            }
            let mut seen: Option<u8> = None;
            let mut differs = false;
            for name in &reference_names {
                let base = projections[name][p];
                match seen {
                    None => seen = Some(base),
                    Some(s) if s != base => {
                        differs = true;
                        break;
                    }
                    _ => {}
                }
            }
            if differs {
                positions.push(p);
            }
        }

        let signatures = reference_names
            .iter()
            .map(|name| {
                let proj = &projections[name];
                (name.clone(), positions.iter().map(|&p| proj[p]).collect::<Vec<u8>>())
            })
            .collect();

        Ok(DiscriminatingClassifier {
            anchor_name,
            anchor_seq,
            positions,
            signatures,
            reference_names,
            min_margin: min_margin.max(1),
        })
    }

    /// Build from a loaded [`crate::reference::fasta_reference::ReferenceManager`].
    pub fn from_reference_manager(
        rm: &crate::reference::fasta_reference::ReferenceManager,
        min_margin: usize,
    ) -> Result<Self, String> {
        let refs: Vec<(String, Vec<u8>)> = rm
            .references
            .values()
            .map(|r| (String::from_utf8_lossy(&r.name).into_owned(), r.sequence.clone()))
            .collect();
        DiscriminatingClassifier::new(&refs, min_margin)
    }

    /// Number of discriminating positions the panel has (0 = references are
    /// identical in anchor coordinates, i.e. not discriminable here).
    pub fn n_positions(&self) -> usize {
        self.positions.len()
    }

    /// Classify a read sequence: align it to the anchor, then score its bases at
    /// the discriminating positions against every reference's signature.
    pub fn classify_read(&self, read_seq: &[u8]) -> Classification {
        let scoring = AffineScoring::default_dna();
        let projection = project_onto_anchor(
            &self.anchor_seq,
            &self.anchor_name,
            read_seq,
            &"read".to_string(),
            &scoring,
        );
        self.classify_projection(&projection)
    }

    /// Classify a read already projected onto anchor coordinates (`read_proj[p]`
    /// is the read's base at anchor position `p`, `-` for gap). Exposed mainly
    /// for testing and for callers that already hold a read↔anchor alignment.
    pub fn classify_projection(&self, read_proj: &[u8]) -> Classification {
        // Score each reference by matches at the discriminating positions.
        let mut scores: Vec<(String, usize)> = self
            .reference_names
            .iter()
            .map(|name| {
                let sig = &self.signatures[name];
                let matches = self
                    .positions
                    .iter()
                    .enumerate()
                    .filter(|(i, &p)| read_proj.get(p).copied().unwrap_or(GAP) == sig[*i])
                    .count();
                (name.clone(), matches)
            })
            .collect();

        // Rank by score desc, name asc for a stable tie order.
        scores.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));

        let informative_positions = self
            .positions
            .iter()
            .filter(|&&p| is_real_base(read_proj.get(p).copied().unwrap_or(GAP)))
            .count();

        let (best, best_score) = scores[0].clone();
        let (second, second_score) = match scores.get(1) {
            Some((n, s)) => (Some(n.clone()), *s),
            None => (None, 0),
        };
        let margin = best_score - second_score;
        let ambiguous = self.reference_names.len() >= 2 && margin < self.min_margin;

        Classification {
            best,
            best_score,
            second,
            second_score,
            margin,
            informative_positions,
            total_positions: self.positions.len(),
            ambiguous,
        }
    }

}

/// Whether a projected base is a concrete ACGT base (carries classification signal).
fn is_real_base(b: u8) -> bool {
    matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

/// Align `other` to `anchor` and return `other`'s base at each anchor position
/// (`-` where `other` has a deletion). Insertions in `other` relative to the
/// anchor fall on anchor-gap columns and are dropped. Result length == anchor.len().
fn project_onto_anchor(
    anchor: &[u8],
    anchor_name: &String,
    other: &[u8],
    other_name: &String,
    scoring: &AffineScoring,
) -> Vec<u8> {
    let aln = align_two_strings(
        &anchor.to_vec(),
        &other.to_vec(),
        None,
        scoring,
        false, // global
        anchor_name,
        other_name,
        None,
    );

    let mut proj = vec![GAP; anchor.len()];
    let mut anchor_pos = 0usize;
    for (a, o) in aln.reference_aligned.iter().zip(aln.read_aligned.iter()) {
        if *a != GAP {
            if anchor_pos < proj.len() {
                proj[anchor_pos] = *o;
            }
            anchor_pos += 1;
        }
        // *a == GAP: an insertion in `other` relative to the anchor -> skip.
    }
    proj
}

#[cfg(test)]
mod tests {
    use super::*;

    // The real RNF2 near-identical panel: three 156 bp amplicons that differ by
    // ~2 bp each (a base-editing palindrome experiment). The unique-k-mer router
    // empties out on these, which is exactly the case this classifier targets.
    const UNEDITED: &str = "GCCTCCACGGCCACTAGTATTATGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACGTATTAGTCATCGCTATTACCATGTACTCATCCTGTCATCTTAGCTAAGATGACAGGTAATTCGAATTTAAATCGGATCCGCGGCC";
    const LEFT: &str = "GCCTCCACGGCCACTAGTATTATGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACGTATTAGTCATCGCTATTACCATGCGATGCTCCTGTCGTCTTAGCTAAGATGACAGGTAATTCGAATTTAAATCGGATCCGCGGCC";
    const RIGHT: &str = "GCCTCCACGGCCACTAGTATTATGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACGTATTAGTCATCGCTATTACCATGCATACTTCCTGTCATCTTAGCTAAGACGACAGGTAATTCGAATTTAAATCGGATCCGCGGCC";

    fn panel() -> Vec<(String, Vec<u8>)> {
        vec![
            ("unedited".to_string(), UNEDITED.as_bytes().to_vec()),
            ("left".to_string(), LEFT.as_bytes().to_vec()),
            ("right".to_string(), RIGHT.as_bytes().to_vec()),
        ]
    }

    #[test]
    fn test_finds_discriminating_positions_on_rnf2_panel() {
        let clf = DiscriminatingClassifier::new(&panel(), 1).unwrap();
        // The panel differs at only a few columns; must be > 0 and small.
        assert!(clf.n_positions() > 0, "expected discriminating positions");
        assert!(clf.n_positions() < 10, "panel should differ at only a few positions, got {}", clf.n_positions());
        // Each reference has a signature the length of `positions`.
        for name in ["unedited", "left", "right"] {
            assert_eq!(clf.signatures[name].len(), clf.positions.len());
        }
    }

    #[test]
    fn test_classifies_each_exact_reference_to_itself() {
        let clf = DiscriminatingClassifier::new(&panel(), 1).unwrap();
        for (name, seq) in panel() {
            let c = clf.classify_read(&seq);
            assert_eq!(c.best, name, "reference {} misclassified as {}", name, c.best);
            assert!(!c.ambiguous, "exact reference {} should be unambiguous (margin {})", name, c.margin);
            assert!(c.margin >= 1);
        }
    }

    #[test]
    fn test_backbone_only_read_is_ambiguous() {
        // A read identical to the shared backbone but N-masked at the
        // discriminating positions carries no signal -> ambiguous, all tie.
        let clf = DiscriminatingClassifier::new(&panel(), 1).unwrap();
        let mut read = UNEDITED.as_bytes().to_vec();
        for &p in &clf.positions {
            read[p] = b'N';
        }
        let c = clf.classify_read(&read);
        assert!(c.ambiguous, "backbone-only read should be ambiguous");
        assert_eq!(c.margin, 0);
        assert_eq!(c.informative_positions, 0);
    }

    #[test]
    fn test_single_reference_is_never_ambiguous() {
        let single = vec![("only".to_string(), UNEDITED.as_bytes().to_vec())];
        let clf = DiscriminatingClassifier::new(&single, 1).unwrap();
        let c = clf.classify_read(UNEDITED.as_bytes());
        assert_eq!(c.best, "only");
        assert!(!c.ambiguous);
        assert!(c.second.is_none());
    }

    #[test]
    fn test_classify_projection_direct() {
        // A minimal 3-reference panel differing at positions 2 and 5.
        //           0123456
        let a = b"ACAGTAA".to_vec(); // ref_a
        let b = b"ACCGTGA".to_vec(); // ref_b (pos2 A->C, pos5 A->G)
        let c = b"ACAGTGA".to_vec(); // ref_c (pos5 A->G only)
        let refs = vec![
            ("a".to_string(), a.clone()),
            ("b".to_string(), b.clone()),
            ("c".to_string(), c.clone()),
        ];
        let clf = DiscriminatingClassifier::new(&refs, 1).unwrap();
        assert_eq!(clf.positions, vec![2, 5]);

        // A read matching ref_b exactly at both discriminating positions.
        let cls = clf.classify_projection(&b);
        assert_eq!(cls.best, "b");
        assert_eq!(cls.best_score, 2);
        assert!(!cls.ambiguous);

        // A read that is C at pos2 (only b) but A at pos5 (a and c): b scores 1,
        // c scores 1 (pos5 A? c has G at 5) ... check the tie logic explicitly.
        // read: pos2=C (b), pos5=A -> a matches pos5(A) + pos2? a has A at 2 =>
        // a: pos2 A!=C, pos5 A==A -> 1 ; b: pos2 C==C, pos5 G!=A -> 1 ; c: pos2 A!=C, pos5 G!=A -> 0
        let read = b"ACCGTAA".to_vec();
        let cls2 = clf.classify_projection(&read);
        assert_eq!(cls2.best_score, 1);
        assert_eq!(cls2.second_score, 1);
        assert!(cls2.ambiguous, "1-vs-1 tie must be ambiguous");
    }

    #[test]
    fn test_margin_threshold_respected() {
        // Two references differing at exactly one position -> max margin is 1.
        let refs = vec![
            ("x".to_string(), b"AAAA".to_vec()),
            ("y".to_string(), b"AACA".to_vec()),
        ];
        // With min_margin = 2, even a perfect match can't clear the bar.
        let strict = DiscriminatingClassifier::new(&refs, 2).unwrap();
        let c = strict.classify_read(b"AACA");
        assert_eq!(c.best, "y");
        assert_eq!(c.margin, 1);
        assert!(c.ambiguous, "margin 1 < min_margin 2 must be ambiguous");

        // With min_margin = 1 the same read is a confident call.
        let lax = DiscriminatingClassifier::new(&refs, 1).unwrap();
        assert!(!lax.classify_read(b"AACA").ambiguous);
    }

    #[test]
    fn test_empty_and_duplicate_panels_error() {
        assert!(DiscriminatingClassifier::new(&[], 1).is_err());
        let dup = vec![
            ("a".to_string(), b"ACGT".to_vec()),
            ("a".to_string(), b"ACGA".to_vec()),
        ];
        assert!(DiscriminatingClassifier::new(&dup, 1).is_err());
    }

    /// Validation against real reads. Reads a TSV of `read_seq \t clique_ref`
    /// (clique's alignment-based reference call) from the path in the
    /// `CLIQUE_DISCRIM_TSV` env var, classifies each read with this module, and
    /// prints concordance + ambiguity stats. Ignored by default (needs data):
    ///   CLIQUE_DISCRIM_TSV=reads_ref.tsv cargo test --bin clique \
    ///     discriminating::tests::validate_on_real_reads -- --ignored --nocapture
    #[test]
    #[ignore]
    fn validate_on_real_reads() {
        use std::io::BufRead;
        let path = std::env::var("CLIQUE_DISCRIM_TSV").expect("set CLIQUE_DISCRIM_TSV");
        let reads: Vec<(String, String)> = std::io::BufReader::new(std::fs::File::open(&path).unwrap())
            .lines()
            .map(|l| {
                let l = l.unwrap();
                let mut it = l.split('\t');
                (it.next().unwrap().to_string(), it.next().unwrap().to_string())
            })
            .collect();

        // The ABE editing window around the palindrome target (observed edit
        // sites from the `ce` tags). Excluding these leaves only the structural
        // allele differences (~95-100). NOTE: in this particular palindrome panel
        // the edit sites (108/121) are also reference-defining, so excluding them
        // removes real signal -- the exclusion mechanism is meant for panels where
        // edit windows do NOT coincide with the allele differences.
        let edit_window: Vec<usize> = (101..=125).collect();

        let run = |clf: &DiscriminatingClassifier, label: &str| {
            let (mut agree, mut ambiguous, mut disagree) = (0, 0, 0);
            for (seq, clique_ref) in &reads {
                let c = clf.classify_read(seq.as_bytes());
                if c.ambiguous { ambiguous += 1; }
                else if &c.best == clique_ref { agree += 1; }
                else { disagree += 1; }
            }
            let confident = agree + disagree;
            let conc = if confident == 0 { 1.0 } else { agree as f64 / confident as f64 };
            eprintln!(
                "[{}] positions={:?}\n    reads={} agree={} ({:.1}%) ambiguous={} disagree={} | confident-concordance={:.1}%",
                label, clf.positions, reads.len(),
                agree, 100.0 * agree as f64 / reads.len() as f64,
                ambiguous, disagree, 100.0 * conc,
            );
            conc
        };

        let all = DiscriminatingClassifier::new(&panel(), 1).unwrap();
        let no_edits = DiscriminatingClassifier::new_excluding(&panel(), 1, &edit_window).unwrap();
        eprintln!("anchor={}", all.anchor_name);
        let conc_all = run(&all, "all discriminating positions");
        let _conc_no_edits = run(&no_edits, "structural positions only (edit window excluded)");

        // The aligner is the established (not ground-truth) method; on confident
        // calls the classifier should agree with it strongly.
        assert!(conc_all > 0.95, "confident-call concordance with the aligner too low: {:.3}", conc_all);
    }

    #[test]
    fn test_read_with_error_outside_discriminating_positions_still_correct() {
        // An error in the backbone (not at a discriminating position) must not
        // change the call: this is the whole point vs. a full-length score.
        let clf = DiscriminatingClassifier::new(&panel(), 1).unwrap();
        let mut read = LEFT.as_bytes().to_vec();
        // corrupt a backbone base that is NOT a discriminating position
        let backbone_pos = (0..read.len()).find(|p| !clf.positions.contains(p)).unwrap();
        read[backbone_pos] = if read[backbone_pos] == b'A' { b'T' } else { b'A' };
        let c = clf.classify_read(&read);
        assert_eq!(c.best, "left", "backbone error should not flip the call");
        assert!(!c.ambiguous);
    }
}
