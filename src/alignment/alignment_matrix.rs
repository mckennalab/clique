use ndarray::prelude::*;
use ndarray::Array;
use fasta_comparisons::*;

pub enum AlignmentType {
    SIMPLE,
    AFFINE,
    CONVEX,
    ANCHORED_CONVEX,
    INVERSION,
    ANCHORED_INVERSION,
}

pub enum AlignmentDirection {
    UP,
    LEFT,
    DIAG,
}

/// The basic container for our alignments
struct AlignmentMatrix<F: ScoringFunction> {
    pub alignment_matrix: Array::<i64, _>,
    pub traceback_matrix: Array::<i64, _>,
    pub alignment_type: AlignmentType,
    pub scoring_function: F,
}

/// Trait required to instantiate a Scoring instance
pub trait ScoringFunction {
    fn score(&self, a: u8, b: u8, gapSize: u32, alignmentDirection: AlignmentDirection) -> i64;
}


/// Create new aligner instance. The size hints help to
/// avoid unnecessary memory allocations.
///
/// # Arguments
///
/// * `m` - the expected size of x
/// * `n` - the expected size of y
/// * `gap_open` - the score for opening a gap (should be negative)
/// * `gap_extend` - the score for extending a gap (should be negative)
/// * `match_fn` - function that returns the score for substitutions
///    (see also [`bio::alignment::pairwise::Scoring`](struct.Scoring.html))
pub fn with_capacity<F: ScoringFunction> (m: usize, n: usize, alignment_type: AlignmentType, scoringFunction: F) -> Self {
    assert!(gap_open <= 0, "gap_open can't be positive");
    assert!(gap_extend <= 0, "gap_extend can't be positive");

    Aligner {
        I: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
        D: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
        S: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
        Lx: Vec::with_capacity(n + 1),
        Ly: Vec::with_capacity(m + 1),
        Sn: Vec::with_capacity(m + 1),
        traceback: Traceback::with_capacity(m, n),
        scoring: Scoring::new(gap_open, gap_extend, match_fn),
    }
}


impl<F: MatchFunc> AlignmentMatrix<F> {
    /// Create new aligner instance with given gap open and gap extend penalties
    /// and the score function.
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should be negative)
    /// * `gap_extend` - the score for extending a gap (should be negative)
    /// * `match_fn` - function that returns the score for substitutions
    ///    (see also [`bio::alignment::pairwise::Scoring`](struct.Scoring.html))
    pub fn new(gap_open: i32, gap_extend: i32, match_fn: F) -> Self {
        Aligner::with_capacity(
            DEFAULT_ALIGNER_CAPACITY,
            DEFAULT_ALIGNER_CAPACITY,
            gap_open,
            gap_extend,
            match_fn,
        )
    }
}