use std::cmp::{max, min};
use std::collections::HashMap;
use std::fmt;
use std::ops::Add;

use ndarray::Array;
use ndarray::prelude::*;
use num_traits::identities::Zero;
use log::{trace};
use serde::{Deserialize, Serialize};
use crate::alignment::fasta_bit_encoding::{FASTA_UNSET, FastaBase, reverse_complement};

use crate::alignment::scoring_functions::*;

use noodles_sam;


use noodles_sam::alignment::RecordBuf;

use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::{Cigar, Data, Name, QualityScores};
use crate::consensus::consensus_builders::get_reference_alignment_rate;

use noodles_sam::alignment::record::cigar::op::Kind;

use crate::alignment_manager::simplify_cigar_string;
use noodles_sam::alignment::record::cigar::Op as Op;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record_buf::Cigar as CigarBuf;

pub const MAX_NEG_SCORE: f64 = -100000.0;


#[derive(Copy, Clone, Debug)]
pub struct MatchedPosition {
    pub search_start: usize,
    pub ref_start: usize,
    pub length: usize,
}

/// Where our alignment starts, how the sequences align, and the sequences themselves.
#[derive(Clone, Debug)]
pub struct SharedSegments {
    pub start_position: usize,
    pub alignment_segments: Vec<MatchedPosition>,
}

pub struct AlignmentCigar {
    pub alignment_start: usize,
    pub alignment_tags: Vec<AlignmentTag>,
}

/// Our alignment tags -- we don't support clipping for our global alignment approach(es)
#[derive(Eq, PartialEq, Debug, Copy, Clone, Hash, Serialize, Deserialize)]
pub enum AlignmentTag {
    MatchMismatch(usize),
    Del(usize),
    Ins(usize),
    SoftClip(usize),
    HardClip(usize),
    InversionOpen,
    InversionClose,
}

impl From<u8> for AlignmentTag {
    fn from(value: u8) -> Self {
        match value {
            b'X' | b'M' => AlignmentTag::MatchMismatch(1),
            b'D' => AlignmentTag::Del(1),
            b'I' => AlignmentTag::Ins(1),
            _ => { panic!("Cannot convert {} to AlignmentTag", value) }
        }
    }
}
impl From<Op> for AlignmentTag {
    fn from(value: Op) -> Self {
        match value.kind() {
            Kind::Match => {AlignmentTag::MatchMismatch(value.len())}
            Kind::Deletion => {AlignmentTag::Del(value.len())}
            Kind::Insertion => {AlignmentTag::Ins(value.len())}
            Kind::SoftClip => {AlignmentTag::SoftClip(value.len())}
            Kind::HardClip => {AlignmentTag::HardClip(value.len())}
            _ => {panic!("Unmatched op tag")}
        }
    }
}


impl AlignmentTag {
    pub fn to_op(&self) -> Op {
        match self {
            AlignmentTag::MatchMismatch(x) => {Op::new(Kind::Match, *x)}
            AlignmentTag::Del(x) => {Op::new(Kind::Deletion, *x)}
            AlignmentTag::Ins(x) => {Op::new(Kind::Insertion, *x)}
            AlignmentTag::SoftClip(x) => {Op::new(Kind::SoftClip, *x)}
            AlignmentTag::InversionOpen => {panic!("Not sure how to translate the inversion open tag")}
            AlignmentTag::InversionClose => {panic!("Not sure how to translate the inversion close tag")}
            AlignmentTag::HardClip(x) => {Op::new(Kind::HardClip, *x)}
        }
    }
}

impl fmt::Display for AlignmentTag {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            AlignmentTag::MatchMismatch(size) => write!(f, "{}M", size),
            AlignmentTag::Del(size) => write!(f, "{}D", size),
            AlignmentTag::Ins(size) => write!(f, "{}I", size),
            AlignmentTag::InversionOpen => write!(f, "<"),
            AlignmentTag::InversionClose => write!(f, ">"),
            AlignmentTag::SoftClip(size) => write!(f, "{}S", size),
            AlignmentTag::HardClip(size) => write!(f, "{}H", size),
        }
    }
}

#[allow(dead_code)]
#[derive(Eq, PartialEq, Debug, Clone, Hash, Copy)]
pub enum AlignmentType {
    Simple,
    Affine,
    Convex,
    AnchoredConvex,
    InversionAwareSimple,
    AnchoredInversionAwareSimple,
    AnchoredAffineInversion,
}

#[derive(Eq, PartialEq, Debug, Copy, Clone, Hash)]
pub enum InvMove {
    Up(usize),
    Left(usize),
    Diag(usize),
}
#[allow(dead_code)]
#[derive(Eq, PartialEq, Debug, Copy, Clone, Hash)]
pub enum AlignmentDirection {
    Up(usize),
    Left(usize),
    Diag(usize),
    Inv(AlignmentLocation, AlignmentLocation, InvMove),
}

impl fmt::Display for AlignmentDirection {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            AlignmentDirection::Up(_size) => write!(f, "| "),
            AlignmentDirection::Left(_size) => write!(f, "<-"),
            AlignmentDirection::Diag(_size) => write!(f, "\\ "),
            AlignmentDirection::Inv(pos1, pos2, jump_to) => write!(f, "I {:?},{:?},{:?}", pos1, pos2, jump_to),
        }
    }
}

impl Add for AlignmentDirection {
    type Output = AlignmentDirection;

    fn add(self, rhs: Self) -> Self::Output {
        match self {
            AlignmentDirection::Up(size) => {
                match rhs {
                    AlignmentDirection::Up(size2) => AlignmentDirection::Up(size + size2),
                    AlignmentDirection::Left(_size2) => panic!("adding discordant types: UP AND LEFT"),
                    AlignmentDirection::Diag(_size2) => panic!("adding discordant types: UP AND DIAG"),
                    AlignmentDirection::Inv(_x1, _y1, _jump_to) => panic!("adding discordant types: UP AND INV"),
                }
            }
            AlignmentDirection::Diag(size) => {
                match rhs {
                    AlignmentDirection::Up(_size2) => panic!("adding discordant types: DIAG AND UP"),
                    AlignmentDirection::Left(_size2) => panic!("adding discordant types: DIAG AND LEFT"),
                    AlignmentDirection::Diag(size2) => AlignmentDirection::Diag(size + size2),
                    AlignmentDirection::Inv(_x1, _y1, _jump_to) => panic!("adding discordant types: DIAG AND INV"),
                }
            }
            AlignmentDirection::Left(size) => {
                match rhs {
                    AlignmentDirection::Up(_size2) => panic!("adding discordant types: LEFT AND UP"),
                    AlignmentDirection::Left(size2) => AlignmentDirection::Left(size + size2),
                    AlignmentDirection::Diag(_size2) => panic!("adding discordant types: LEFT AND DIAG"),
                    AlignmentDirection::Inv(_x1, _y1, _jump_to) => panic!("adding discordant types: LEFT AND INV"),
                }
            }
            AlignmentDirection::Inv(_x1, _y1, _jump_to) => {
                match rhs {
                    AlignmentDirection::Up(_size2) => panic!("adding discordant types: INV AND UP"),
                    AlignmentDirection::Left(_size2) => panic!("adding discordant types: INV AND LEFT"),
                    AlignmentDirection::Diag(_size2) => panic!("adding discordant types: INV AND DIAG"),
                    AlignmentDirection::Inv(_x1, _y1, _jump_to) => panic!("CANT ADD TWO INVERSIONS"),
                }
            }
        }
    }
}

impl Zero for AlignmentDirection {
    fn zero() -> Self {
        AlignmentDirection::Up(0)
    }

    fn is_zero(&self) -> bool {
        match self {
            AlignmentDirection::Up(size) => *size == 0,
            AlignmentDirection::Left(size) => *size == 0,
            AlignmentDirection::Diag(size) => *size == 0,
            AlignmentDirection::Inv(_x1, _y1, _jump_to) => false, // any inversion has a size != 0
        }
    }
}

#[derive(Clone)]
pub struct Alignment<K> {
    pub scores: Array<f64, K>,
    pub traceback: Array<AlignmentDirection, K>,
    pub alignment_type: AlignmentType,
    pub is_local: bool,
}

pub fn create_scoring_record_3d(hint_seq_a_len: usize, hint_seq_b_len: usize, alignment_type: AlignmentType, local_alignment: bool) -> Alignment<Ix3> {
    Alignment {
        scores: Array::<f64, Ix3>::zeros((hint_seq_a_len, hint_seq_b_len, 3).f()),
        traceback: Array::<AlignmentDirection, Ix3>::zeros((hint_seq_a_len, hint_seq_b_len, 3).f()),
        alignment_type,
        is_local: local_alignment,
    }
}
/// For inversions we need a way to update the alignment matrix, removing 'used' alignment to find the next best alignment
///
/// This function iterates over the alignment matrix starting from the specified `row` and `column`, updating
/// scores using a specified scoring function until no more updates can be made or the edge of the matrix is reached.
/// The function can iterate either by rows or columns based on the `by_row` flag.
///
/// # Parameters
/// - `alignment`: A mutable reference to an `Alignment<Ix3>` object representing the 3D alignment matrix.
/// - `sequence1`: A slice of `FastaBase` representing the first sequence involved in the alignment.
/// - `sequence2`: A slice of `FastaBase` representing the second sequence involved in the alignment.
/// - `scoring_function`: A reference to an `AffineScoring` object used to calculate score updates during the alignment.
/// - `row`: The starting row index for the update process.
/// - `column`: The starting column index for the update process.
/// - `by_row`: A boolean flag indicating whether the update should proceed by rows (`true`) or by columns (`false`).
///
/// # Returns
/// Returns the number of updates made to the alignment matrix.
///
/// # Examples
/// ```
/// // Example usage of `update_sub_vector3d`
/// let mut alignment = Alignment::new(...); // Assume a valid 3D alignment matrix is initialized here
/// let sequence1 = ...; // Assume sequence1 is initialized here
/// let sequence2 = ...; // Assume sequence2 is initialized here
/// let scoring_function = AffineScoring::new(...); // Assume a scoring function is initialized here
/// let row = 0;
/// let column = 0;
/// let by_row = true;
///
/// let updates = update_sub_vector3d(&mut alignment, &sequence1, &sequence2, &scoring_function, row, column, by_row);
/// println!("Number of updates made: {}", updates);
/// ```
#[allow(dead_code)]
#[inline(always)]
fn update_sub_vector3d(alignment: &mut Alignment<Ix3>,
                       sequence1: &[FastaBase],
                       sequence2: &[FastaBase],
                       scoring_function: &AffineScoring,
                       row: usize,
                       column: usize,
                       by_row: bool) -> usize {
    let mut row_pos = if by_row { row + 1 } else { row };
    let mut col_pos = if by_row { column } else { column + 1 };
    let mut still_updating = true;
    let mut update_count = 0;

    while row_pos < alignment.scores.shape()[0] && col_pos < alignment.scores.shape()[1] && still_updating {
        let updates = update_3d_score_local(alignment, sequence1, sequence2, scoring_function, row_pos, col_pos);

        let any_update = updates.0 || updates.1 || updates.2;
        if any_update {
            if by_row {
                row_pos += 1;
            } else {
                col_pos += 1;
            }

            update_count += 1;
            still_updating = true;
        } else {
            still_updating = false;
        }
    }
    update_count
}

/// Cleans up the alignment based on the previous result and attempts to find the next best match
/// in a 3D alignment matrix.
///
/// This function iterates through the path of a previous alignment result, updating the 3D alignment matrix
/// for each position in the path by attempting to find better matches both by rows and by columns. After
/// completing updates for the path, it continues to search for better matches beyond the path's end, updating
/// the alignment matrix until no further improvements can be made or until it reaches the edge of the matrix.
///
/// # Parameters
/// - `alignment`: A mutable reference to an `Alignment<Ix3>` object representing the 3D alignment matrix.
/// - `sequence1`: A slice of `FastaBase` representing the first sequence involved in the alignment.
/// - `sequence2`: A slice of `FastaBase` representing the second sequence involved in the alignment.
/// - `scoring_function`: A reference to an `AffineScoring` object used to calculate score updates during the alignment.
/// - `previous_result`: A reference to an `AlignmentResult` object containing the path of the previously found alignment.
///
/// # Side Effects
/// - Updates the alignment matrix within the `alignment` parameter based on the path provided in `previous_result`
///   and continues to search for and apply updates beyond the path's end.
///
/// # Examples
/// ```
/// // Example usage of `clean_and_find_next_best_match_3d`
/// let mut alignment = Alignment::new(...); // Assume a valid 3D alignment matrix is initialized here
/// let sequence1 = ...; // Assume sequence1 is initialized here
/// let sequence2 = ...; // Assume sequence2 is initialized here
/// let scoring_function = AffineScoring::new(...); // Assume a scoring function is initialized here
/// let previous_result = ...; // Assume a previous alignment result is available
///
/// clean_and_find_next_best_match_3d(&mut alignment, &sequence1, &sequence2, &scoring_function, &previous_result);
/// // The alignment matrix within `alignment` is now updated based on `previous_result` and further refinement.
/// ```
#[inline(always)]
#[allow(dead_code)]
fn clean_and_find_next_best_match_3d(alignment: &mut Alignment<Ix3>,
                                     sequence1: &[FastaBase],
                                     sequence2: &[FastaBase],
                                     scoring_function: &AffineScoring,
                                     previous_result: &AlignmentResult) {
    let mut current_row = 0;
    let mut current_col = 0;
    for path_entry in &previous_result.path {
        current_row = path_entry.x;
        current_col = path_entry.y;
        for _z in 0..3 {
            update_sub_vector3d(alignment, sequence1, sequence2, scoring_function, current_row, current_col, true);
            update_sub_vector3d(alignment, sequence1, sequence2, scoring_function, current_row, current_col, false);
        }
    }

    let mut still_updating_rows = true;
    let mut still_updating_cols = true;
    while (still_updating_rows || still_updating_cols) &&
        current_row < alignment.scores.shape()[0] &&
        current_col < alignment.scores.shape()[1] {
        let row_update_count = update_sub_vector3d(alignment, sequence1, sequence2, scoring_function, current_row, current_col, true);
        let col_update_count = update_sub_vector3d(alignment, sequence1, sequence2, scoring_function, current_row, current_col, false);
        still_updating_rows = row_update_count > 0;
        still_updating_cols = col_update_count > 0;
        current_row += 1;
        current_col += 1;
    }
}


/// Affine matrix dimensions are row,column,dimension, where dim 1 is match, dim 2 is deletion (relative to read, sequence2) and dim 3 is insertion
pub fn perform_affine_alignment(alignment: &mut Alignment<Ix3>,
                                sequence1: &[FastaBase],
                                sequence2: &[FastaBase],
                                scoring_function: &AffineScoring) {
    let new_bandwidth = max(sequence1.len(), sequence2.len());
    perform_affine_alignment_bandwidth(alignment, sequence1, sequence2, scoring_function, &new_bandwidth)
}


/// Affine matrix dimensions are row,column,dimension, where dim 1 is match, dim 2 is deletion (relative to read, sequence2) and dim 3 is insertion
pub fn perform_affine_alignment_bandwidth(alignment: &mut Alignment<Ix3>,
                                          sequence1: &[FastaBase],
                                          sequence2: &[FastaBase],
                                          scoring_function: &AffineScoring,
                                          bandwidth: &usize) {
    assert_eq!(alignment.scores.shape()[2], 3);
    assert!(alignment.scores.shape()[0] > sequence1.len(), "Asked to align sequence 1 with length {} in a matrix sized {} in that dimension, sequence {}", sequence1.len() + 1, alignment.scores.shape()[0], FastaBase::string_from_slice(sequence1));
    assert!(alignment.scores.shape()[1] > sequence2.len(), "Asked to align sequence 2 with length {} in a matrix sized {} in that dimension, sequence {}", sequence2.len() + 1, alignment.scores.shape()[1], FastaBase::string_from_slice(sequence2));

    alignment.scores[[0, 0, 0]] = 0.0;
    alignment.scores[[0, 0, 1]] = MAX_NEG_SCORE;
    alignment.scores[[0, 0, 2]] = MAX_NEG_SCORE;

    // first column (going down)
    for x in 1..(sequence1.len() + 1) {
        alignment.scores[[x, 0, 0]] = MAX_NEG_SCORE;
        alignment.traceback[[x, 0, 0]] = AlignmentDirection::Up(1);
        alignment.scores[[x, 0, 1]] = (scoring_function.gap_open() + (x as f64 * scoring_function.gap_extend())) * scoring_function.final_gap_multiplier();
        alignment.traceback[[x, 0, 1]] = AlignmentDirection::Up(1);
        alignment.scores[[x, 0, 2]] = (scoring_function.gap_open() + (x as f64 * scoring_function.gap_extend())) * scoring_function.final_gap_multiplier();
        alignment.traceback[[x, 0, 2]] = AlignmentDirection::Up(1);
    }
    // top row
    for y in 1..(sequence2.len() + 1) {
        alignment.scores[[0, y, 0]] = MAX_NEG_SCORE;
        alignment.traceback[[0, y, 0]] = AlignmentDirection::Left(1);
        alignment.scores[[0, y, 1]] = (scoring_function.gap_open() + (y as f64 * scoring_function.gap_extend())) * scoring_function.final_gap_multiplier();
        alignment.traceback[[0, y, 1]] = AlignmentDirection::Left(1);
        alignment.scores[[0, y, 2]] = (scoring_function.gap_open() + (y as f64 * scoring_function.gap_extend())) * scoring_function.final_gap_multiplier();
        alignment.traceback[[0, y, 2]] = AlignmentDirection::Left(1);
    }

    let update_function = match alignment.is_local {
        true => { update_3d_score_local }
        false => { update_3d_score }
    };

    for x in 1..(sequence1.len() + 1) {
        let y_bounds = ((x as f64 / (sequence1.len() + 1) as f64) * (sequence2.len() + 1) as f64) as i64;
        //println!("y bounds {} {} {}",y_bounds,(x as f64 / (sequence1.len() + 1) as f64),(sequence2.len() + 1) as f64);

        let y_bounds = (max(1, y_bounds - (*bandwidth as i64)), min(sequence2.len() as i64 + 1, y_bounds + (*bandwidth as i64)));
        //println!("y bounds {} {} from x = {} bandwidth = {}",y_bounds.0,y_bounds.1,x,bandwidth);
        assert!(y_bounds.0 >= 0);
        assert!(y_bounds.1 >= 0);
        for y in (y_bounds.0)..(y_bounds.1) {
            update_function(alignment, sequence1, sequence2, scoring_function, x, y as usize);
        }
    }
}

/// Affine matrix dimensions are row,column,dimension, where dim 1 is match, dim 2 is deletion (relative to read, sequence2) and dim 3 is insertion
#[allow(dead_code)]
fn perform_inversion_aware_alignment(alignment: &mut Alignment<Ix3>,
                                     alignment_inversion: &HashMap<AlignmentLocation, BoundedAlignment>,
                                     sequence1: &Vec<FastaBase>,
                                     sequence2: &Vec<FastaBase>,
                                     scoring_function: &InversionScoring) {
    assert_eq!(alignment.scores.shape()[2], 3);
    assert_eq!(alignment.scores.shape()[0], sequence1.len() + 1);
    assert_eq!(alignment.scores.shape()[1], sequence2.len() + 1);

    alignment.scores[[0, 0, 0]] = 0.0;
    alignment.scores[[0, 0, 1]] = MAX_NEG_SCORE;
    alignment.scores[[0, 0, 2]] = MAX_NEG_SCORE;

    // first column (going down)
    for x in 1..alignment.scores.shape()[0] {
        alignment.scores[[x, 0, 0]] = MAX_NEG_SCORE;
        alignment.traceback[[x, 0, 0]] = AlignmentDirection::Up(1);
        alignment.scores[[x, 0, 1]] = scoring_function.gap_open() + (x as f64 * scoring_function.gap_extend());
        alignment.traceback[[x, 0, 1]] = AlignmentDirection::Up(1);
        alignment.scores[[x, 0, 2]] = scoring_function.gap_open() + (x as f64 * scoring_function.gap_extend());
        alignment.traceback[[x, 0, 2]] = AlignmentDirection::Up(1);
    }
    // top row
    for y in 1..alignment.scores.shape()[1] {
        alignment.scores[[0, y, 0]] = MAX_NEG_SCORE;
        alignment.traceback[[0, y, 0]] = AlignmentDirection::Left(1);
        alignment.scores[[0, y, 1]] = scoring_function.gap_open() + (y as f64 * scoring_function.gap_extend());
        alignment.traceback[[0, y, 1]] = AlignmentDirection::Left(1);
        alignment.scores[[0, y, 2]] = scoring_function.gap_open() + (y as f64 * scoring_function.gap_extend());
        alignment.traceback[[0, y, 2]] = AlignmentDirection::Left(1);
    }

    for x in 1..alignment.scores.shape()[0] {
        for y in 1..alignment.scores.shape()[1] {
            update_inversion_alignment(alignment, alignment_inversion, sequence1, sequence2, scoring_function, x, y);
        }
    }
}

#[allow(dead_code)]
fn update_inversion_alignment(alignment: &mut Alignment<Ix3>,
                              alignment_inversion: &HashMap<AlignmentLocation, BoundedAlignment>,
                              sequence1: &[FastaBase],
                              sequence2: &[FastaBase],
                              scoring_function: &InversionScoring,
                              x: usize,
                              y: usize) -> (bool, bool, bool) {
    let mut _update_x = false;
    let mut _update_y = false;
    let mut _update_z = false;

    {
        let match_score = scoring_function.match_mismatch(&sequence1[x - 1], &sequence2[y - 1]);

        let pos = AlignmentLocation { x, y };

        let max_match_mismatch = [
            if alignment.is_local { 0.0 } else { MAX_NEG_SCORE },
            alignment.scores[[x - 1, y - 1, 0]] + match_score,
            if alignment.is_local { match_score } else { MAX_NEG_SCORE },
        ];

        let max_match_mismatch = max_match_mismatch.iter().max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();

        let best_match = [
            if alignment_inversion.contains_key(&pos.clone()) {
                let inv = alignment_inversion.get(&pos.clone());
                if let Some(inversion) = inv {
                    let first_path = &inversion.bounding_box.0;
                    let last_path = inversion.bounding_box.1;
                    assert_eq!(&last_path, &pos);

                    let inv_best_match = [
                        (alignment.scores[[first_path.x - 1, first_path.y - 1, 1]], InvMove::Up(1)),
                        (alignment.scores[[first_path.x - 1, first_path.y - 1, 2]], InvMove::Left(1)),
                        (alignment.scores[[first_path.x - 1, first_path.y - 1, 0]], InvMove::Diag(1)), ];
                    let inv_best_match = inv_best_match.iter().max_by(|x, y| x.0.partial_cmp(&y.0).unwrap()).unwrap();

                    (inversion.alignment_result.score + inv_best_match.0 + scoring_function.inversion_cost(), AlignmentDirection::Inv(*first_path, last_path, inv_best_match.1))
                } else {
                    (MAX_NEG_SCORE, AlignmentDirection::Up(1))
                }
            } else { (MAX_NEG_SCORE, AlignmentDirection::Up(1)) },
            (*max_match_mismatch, AlignmentDirection::Diag(1)),
            (alignment.scores[[x - 1, y - 1, 1]] + match_score, AlignmentDirection::Up(1)),
            (alignment.scores[[x - 1, y - 1, 2]] + match_score, AlignmentDirection::Left(1)),
        ];

        let best_match = best_match.iter().max_by(|x, y| x.0.partial_cmp(&y.0).unwrap()).unwrap();

        let best_match = match &best_match.1 {
            AlignmentDirection::Inv(_pos1, pos2, jump) => {
                let aln = &alignment_inversion.get(&pos.clone());
                if let Some(align) = aln {
                    let start_pos = &align.bounding_box.0.clone();
                    let inv = AlignmentDirection::Inv(*start_pos, *pos2, *jump);
                    (best_match.0, inv)
                } else {
                    panic!("Unable to unwrap alignment for inversion!");
                }
            }
            _ => {
                *best_match
            }
        };

        _update_x = alignment.scores[[x, y, 0]] != best_match.0;
        alignment.scores[[x, y, 0]] = best_match.0;
        alignment.traceback[[x, y, 0]] = best_match.1;
    }
    {
        let best_gap_x = [
            (alignment.scores[[x - 1, y, 1]] + scoring_function.gap_extend(), AlignmentDirection::Up(1)),
            (alignment.scores[[x - 1, y, 2]] + scoring_function.gap_open() + scoring_function.gap_extend(), AlignmentDirection::Left(1)),
            (alignment.scores[[x - 1, y, 0]] + scoring_function.gap_open() + scoring_function.gap_extend(), AlignmentDirection::Diag(1))];
        let best_gap_x = best_gap_x.iter().max_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
        _update_y = alignment.scores[[x, y, 1]] != best_gap_x.unwrap().0;
        alignment.scores[[x, y, 1]] = best_gap_x.unwrap().0;
        alignment.traceback[[x, y, 1]] = best_gap_x.unwrap().1;
    }
    {
        let best_gap_y = [
            (alignment.scores[[x, y - 1, 1]] + scoring_function.gap_open() + scoring_function.gap_extend(), AlignmentDirection::Up(1)),
            (alignment.scores[[x, y - 1, 2]] + scoring_function.gap_extend(), AlignmentDirection::Left(1)),
            (alignment.scores[[x, y - 1, 0]] + scoring_function.gap_open() + scoring_function.gap_extend(), AlignmentDirection::Diag(1))];
        let best_gap_y = best_gap_y.iter().max_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
        _update_z = alignment.scores[[x, y, 2]] != best_gap_y.unwrap().0;
        alignment.scores[[x, y, 2]] = best_gap_y.unwrap().0;
        alignment.traceback[[x, y, 2]] = best_gap_y.unwrap().1;
    }
    (_update_x, _update_y, _update_z)
}

#[inline(always)]
fn update_3d_score_local(alignment: &mut Alignment<Ix3>, sequence1: &[FastaBase], sequence2: &[FastaBase], scoring_function: &AffineScoring, x: usize, y: usize) -> (bool, bool, bool) {
    let mut _update_x = false;
    let mut _update_y = false;
    let mut _update_z = false;

    let gap_multiplier = if x == sequence1.len() || y == sequence2.len() { scoring_function.final_gap_multiplier() } else { 1.0 };
    let x1 = scoring_function.gap_open() + (scoring_function.gap_extend() * gap_multiplier);

    {
        // match-mismatch matrix update
        let match_score = scoring_function.match_mismatch(&sequence1[x - 1], &sequence2[y - 1]);

        let max_match_mismatch = three_way_max_and_direction(
            &(if alignment.is_local { 0.0 } else { MAX_NEG_SCORE }),
            &(alignment.scores[[x - 1, y - 1, 0]] + match_score),
            &(if alignment.is_local { match_score } else { MAX_NEG_SCORE })).0;

        let best_match = three_way_max_and_direction(
            &(alignment.scores[[x - 1, y - 1, 1]] + match_score),
            &(alignment.scores[[x - 1, y - 1, 2]] + match_score),
            &max_match_mismatch);

        _update_x = alignment.scores[[x, y, 0]] != best_match.0;
        alignment.scores[[x, y, 0]] = best_match.0;
        alignment.traceback[[x, y, 0]] = best_match.1;
    }
    {
        let best_gap_x = three_way_max_and_direction(
            &(alignment.scores[[x - 1, y, 1]] + scoring_function.gap_extend()),
            &(alignment.scores[[x - 1, y, 2]] + x1),
            &(alignment.scores[[x - 1, y, 0]] + x1));

        _update_y = alignment.scores[[x, y, 1]] != best_gap_x.0;
        alignment.scores[[x, y, 1]] = best_gap_x.0;
        alignment.traceback[[x, y, 1]] = best_gap_x.1;
    }
    {
        let best_gap_y = three_way_max_and_direction(
            &(alignment.scores[[x, y - 1, 1]] + x1),
            &(alignment.scores[[x, y - 1, 2]] + (scoring_function.gap_extend())),
            &(alignment.scores[[x, y - 1, 0]] + x1));


        _update_z = alignment.scores[[x, y, 2]] != best_gap_y.0;
        alignment.scores[[x, y, 2]] = best_gap_y.0;
        alignment.traceback[[x, y, 2]] = best_gap_y.1;
    }
    (_update_x, _update_y, _update_z)
}

#[inline(always)]
fn update_3d_score(alignment: &mut Alignment<Ix3>, sequence1: &[FastaBase], sequence2: &[FastaBase], scoring_function: &AffineScoring, x: usize, y: usize) -> (bool, bool, bool) {
    let mut _update_x = false;
    let mut _update_y = false;
    let mut _update_z = false;
    assert!(!alignment.is_local);
    let gap_multiplier = if x == sequence1.len() || y == sequence2.len() { scoring_function.final_gap_multiplier() } else { 1.0 };
    let x1 = scoring_function.gap_open() + (scoring_function.gap_extend() * gap_multiplier);
    let local_gap_ext = scoring_function.gap_extend() * gap_multiplier;

    {
        // match-mismatch matrix update
        let match_score = scoring_function.match_mismatch(&sequence1[x - 1], &sequence2[y - 1]);

        let best_match = three_way_max_and_direction(
            &(alignment.scores[[x - 1, y - 1, 1]] + match_score),
            &(alignment.scores[[x - 1, y - 1, 2]] + match_score),
            &(alignment.scores[[x - 1, y - 1, 0]] + match_score));

        _update_x = alignment.scores[[x, y, 0]] != best_match.0;
        alignment.scores[[x, y, 0]] = best_match.0;
        alignment.traceback[[x, y, 0]] = best_match.1;
    }
    {
        let best_gap_x = three_way_max_and_direction(
            &(alignment.scores[[x - 1, y, 1]] + local_gap_ext),
            &(alignment.scores[[x - 1, y, 2]] + x1),
            &(alignment.scores[[x - 1, y, 0]] + x1));

        _update_y = alignment.scores[[x, y, 1]] != best_gap_x.0;
        alignment.scores[[x, y, 1]] = best_gap_x.0;
        alignment.traceback[[x, y, 1]] = best_gap_x.1;
    }
    {
        let best_gap_y = three_way_max_and_direction(
            &(alignment.scores[[x, y - 1, 1]] + x1),
            &(alignment.scores[[x, y - 1, 2]] + local_gap_ext),
            &(alignment.scores[[x, y - 1, 0]] + x1));


        _update_z = alignment.scores[[x, y, 2]] != best_gap_y.0;
        alignment.scores[[x, y, 2]] = best_gap_y.0;
        alignment.traceback[[x, y, 2]] = best_gap_y.1;
    }
    (_update_x, _update_y, _update_z)
}

/// This was a hot point in the function above; we did a 3-way comparison using vectors which was really costly,
/// mostly due to mallocs and frees. This does a three-way comparison and returns both the value and the
/// traversal direction for alignment. It's about a 2X speedup over the previous version
#[inline(always)]
pub fn three_way_max_and_direction(up_value: &f64, left_value: &f64, diag_value: &f64) -> (f64, AlignmentDirection) {
    if up_value > left_value {
        if up_value > diag_value {
            (*up_value, AlignmentDirection::Up(1))
        } else {
            (*diag_value, AlignmentDirection::Diag(1))
        }
    } else if left_value > diag_value {
        (*left_value, AlignmentDirection::Left(1))
    } else {
        (*diag_value, AlignmentDirection::Diag(1))
    }
}


#[derive(Serialize, Deserialize, Hash, PartialEq, Eq, Debug, Copy, Clone)]
pub struct AlignmentLocation {
    pub x: usize,
    pub y: usize,
}


#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct AlignmentResult {
    pub reference_name: String,
    pub read_name: String,
    pub reference_aligned: Vec<FastaBase>,
    pub read_aligned: Vec<FastaBase>,
    pub read_quals: Option<Vec<u8>>,
    pub cigar_string: Vec<AlignmentTag>,
    pub path: Vec<AlignmentLocation>,
    pub score: f64,
    pub reference_start: usize,
    pub read_start: usize,
    pub bounding_box: Option<(AlignmentLocation, AlignmentLocation)>,
}

#[allow(dead_code)]
impl AlignmentResult {
    pub fn from_match_segment(str1: &[FastaBase],
                              str2: &[FastaBase],
                              reference_name: &String,
                              read_name: &String,
                              start_x: usize,
                              start_y: usize,
                              af_score: &AffineScoring) -> AlignmentResult {
        let cigar_string: Vec<AlignmentTag> = vec![AlignmentTag::MatchMismatch(str1.len())];
        let path = (start_x..(start_x + str1.len())).zip(start_y..(start_y + str1.len())).map(|(x, y)| AlignmentLocation { x, y }).collect();
        let score = str1.iter().zip(str2.iter()).map(|(xb, yb)| af_score.match_mismatch(xb, yb)).sum();

        AlignmentResult {
            reference_name: reference_name.clone(),
            read_name: read_name.clone(),
            reference_aligned: str1.to_vec(),
            read_aligned: str2.to_vec(),
            read_quals: None,
            cigar_string,
            path,
            score,
            reference_start: start_x,
            read_start: start_y,
            bounding_box: None,
        }
    }

    pub fn to_cigar_string(&self) -> Box<dyn noodles_sam::alignment::record::Cigar> {
        let cigar: CigarBuf = self.cigar_string.iter().map(|t|t.to_op()).collect();
        Box::new(cigar)
    }

    pub fn to_sam_record(&self, reference_id: &i32, extra_tags: &HashMap<[u8; 2],String>, read_names: Option<Vec<String>>) -> noodles_sam::alignment::RecordBuf {

        // set the aux. data with alignments stats and the extracted tags
        let mut data = Data::default();
        extra_tags.iter().for_each(|(k,v)| {data.insert(Tag::from(k.clone()), Value::from(v.clone()));});

        data.insert(Tag::from([b'r', b'm']), Value::from(get_reference_alignment_rate(&self.reference_aligned, &self.read_aligned).to_string()));
        data.insert(Tag::from([b'r', b's']), Value::from(self.score.to_string()));

        match read_names {
            None => {}
            Some(x) => {data.insert(Tag::from([b'a', b'r']), Value::from(x.join(",")));}
        }

        data.insert(Tag::from([b'a', b's']), Value::from(self.score.to_string()));

        let seq_len = self.read_aligned.iter().filter(|b| **b != FASTA_UNSET).collect::<Vec<&FastaBase>>().len();
        
        // set the read name
        RecordBuf::builder()
            .set_name(Name::from(self.read_name.as_bytes()))
            .set_sequence(FastaBase::vec_u8(&self.read_aligned.clone().iter().cloned().filter(|b| *b != FASTA_UNSET).collect::<Vec<FastaBase>>()).into())
            .set_cigar(Cigar::from_iter(self.cigar_string.iter().map(|m| m.to_op()).into_iter()).clone())
            .set_alignment_start(noodles_core::Position::new(self.reference_start+1).unwrap())
            .set_quality_scores(match &self.read_quals {
                Some(x) => {QualityScores::from(x.clone())},
                None => {QualityScores::from(vec![b'H'; seq_len])}})
            .set_reference_sequence_id(*reference_id as usize)
            .set_flags(Flags::PROPERLY_ALIGNED)
            .set_data(data).build()
    }

    #[allow(dead_code)]
    fn slice_out_inversions(&self) -> Vec<AlignmentResult> {
        let mut alignment_string1: Vec<FastaBase> = Vec::new();
        let mut alignment_string2: Vec<FastaBase> = Vec::new();
        let mut cigar_string: Vec<AlignmentTag> = Vec::new();
        let mut path = Vec::new();
        let mut score: f64 = 0.0;
        let mut return_vec = Vec::new();

        let mut position = 0;

        for cigar in &self.cigar_string {
            match cigar {
                AlignmentTag::MatchMismatch(size) | AlignmentTag::Ins(size) | AlignmentTag::Del(size) => {
                    alignment_string1.extend(self.reference_aligned[position..(position + size)].iter());
                    alignment_string2.extend(self.read_aligned[position..(position + size)].iter());
                    cigar_string.push(*cigar);
                    path.extend_from_slice(&self.path[position..(position + size)]);
                    score = -1.0;
                    position += size;
                }
                AlignmentTag::InversionOpen | AlignmentTag::InversionClose => {
                    if !alignment_string1.is_empty() || !alignment_string2.is_empty() {
                        return_vec.push(AlignmentResult {
                            reference_name: self.reference_name.clone(),
                            read_name: self.read_name.clone(),
                            reference_aligned: alignment_string1,
                            read_aligned: alignment_string2,
                            read_quals: None,
                            cigar_string,
                            path,
                            score,
                            reference_start: self.reference_start,
                            read_start: self.read_start,
                            bounding_box: None,
                        });
                        alignment_string1 = Vec::new();
                        alignment_string2 = Vec::new();
                        cigar_string = Vec::new();
                        path = Vec::new();
                        score = 0.0;
                    }
                }
                _ => { panic!("Unknown tag type") }
            }
        }
        if !alignment_string1.is_empty() || !alignment_string2.is_empty() {
            return_vec.push(AlignmentResult {
                reference_name: self.reference_name.clone(),
                read_name: self.read_name.clone(),
                reference_aligned: alignment_string1,
                read_aligned: alignment_string2,
                read_quals: None,
                cigar_string,
                path,
                score,
                reference_start: self.reference_start,
                read_start: self.read_start,
                bounding_box: None,
            });
        }
        return_vec
    }

    fn convert_inverted_path(&self, total_string_length: usize) -> AlignmentResult {
        let mut new_path = Vec::new();

        let half_point = total_string_length as f64 / 2.0;

        for path_obj in &self.path {
            let new_y = (1.0 + half_point + (half_point - (path_obj.y as f64))).round() as usize;
            trace!("old {},{} new {}",path_obj.x,path_obj.y,new_y);
            new_path.push(AlignmentLocation { x: path_obj.x, y: new_y });
        }
        new_path.reverse();
        let bounds = (
            AlignmentLocation { x: new_path[new_path.len() - 1].x, y: new_path[0].y },
            AlignmentLocation { x: new_path[0].x, y: new_path[new_path.len() - 1].y });
        AlignmentResult {
            reference_name: self.reference_name.clone(),
            read_name: self.read_name.clone(),
            reference_aligned: self.reference_aligned.clone(),
            read_aligned: self.read_aligned.clone(),
            read_quals: None,
            cigar_string: self.cigar_string.clone(),
            path: new_path,
            score: self.score,
            reference_start: self.reference_start,
            read_start: self.read_start,
            bounding_box: Some(bounds),
        }
    }
}

pub fn find_max_value_3d_array(matrix: &Array<f64, Ix3>) -> Option<(AlignmentLocation, f64)> {
    let mut _g_max_row = 0;
    let mut _g_max_col = 0;
    let mut _g_max_z = 0;
    let mut _g_max_val = MAX_NEG_SCORE;

    for x in 0..matrix.shape()[0] {
        for y in 0..matrix.shape()[1] {
            for z in 0..matrix.shape()[2] {
                // normal update combined with update rule #5 from Waterman and Eggart and
                // rule #6 from Waterman and Eggart
                if matrix[[x, y, z]] > _g_max_val ||
                    matrix[[x, y, z]] == _g_max_val && (x + y) < (_g_max_row + _g_max_col) ||
                    (matrix[[x, y, z]] == _g_max_val &&
                        (_g_max_row + _g_max_col) == (x + y) &&
                        x < _g_max_row) {
                    _g_max_row = x;
                    _g_max_col = y;
                    _g_max_z = z;
                    _g_max_val = matrix[[x, y, z]];
                }
            }
        }
    }

    if _g_max_val > MAX_NEG_SCORE {
        //info!("3d max is {},{},{}={}", g_max_row, g_max_col, g_max_z, g_max_val);
        Some((AlignmentLocation { x: _g_max_row, y: _g_max_col }, _g_max_val))
    } else {
        None
    }
}

pub struct BoundedAlignment {
    alignment_result: AlignmentResult,
    bounding_box: (AlignmentLocation, AlignmentLocation),
}

#[allow(dead_code)]
pub(crate) fn inversion_alignment(reference: &Vec<FastaBase>, read: &Vec<FastaBase>, reference_name: &String, read_name: &String, inversion_score: &InversionScoring, my_aff_score: &AffineScoring, local: bool) -> AlignmentResult {
    let mut alignment_mat = create_scoring_record_3d(reference.len() + 1, read.len() + 1, AlignmentType::Affine, local);
    let mut inversion_mat = create_scoring_record_3d(reference.len() + 1, read.len() + 1, AlignmentType::Affine, true);

    let mut long_enough_hits: HashMap<AlignmentLocation, BoundedAlignment> = HashMap::new();
    let rev_comp_read = &reverse_complement(read);
    perform_affine_alignment(&mut inversion_mat, reference, rev_comp_read, my_aff_score);

    let mut aligned_inv: Option<AlignmentResult> = Some(perform_3d_global_traceback(&mut inversion_mat, None, reference, rev_comp_read, reference_name, read_name, None));

    while aligned_inv.is_some() {
        let aligned_inv_local = aligned_inv.unwrap();
        let length = aligned_inv_local.path.len();
        if length > 1 {
            let converted_path = aligned_inv_local.convert_inverted_path(read.len());
            let bounding = converted_path.bounding_box.unwrap();
            let b_alignment = BoundedAlignment { alignment_result: converted_path, bounding_box: bounding };
            let true_position = AlignmentLocation { x: bounding.1.x, y: bounding.1.y };
            aligned_inv = if length >= inversion_score.min_inversion_length {
                clean_and_find_next_best_match_3d(&mut inversion_mat, reference, rev_comp_read, my_aff_score, &aligned_inv_local);
                long_enough_hits.insert(true_position, b_alignment);
                Some(perform_3d_global_traceback(&mut inversion_mat, None, reference, rev_comp_read, reference_name, read_name, None))
            } else {
                None
            }
        } else {
            aligned_inv = None
        }
    }
    perform_inversion_aware_alignment(&mut alignment_mat, &long_enough_hits, reference, read, inversion_score);
    perform_3d_global_traceback(&mut alignment_mat, Some(&long_enough_hits), reference, read, reference_name, read_name,None)
}

#[allow(dead_code)]
pub fn  perform_3d_global_traceback(alignment: &mut Alignment<Ix3>,
                                   inversion_mapping: Option<&HashMap<AlignmentLocation, BoundedAlignment>>,
                                   sequence1: &[FastaBase],
                                   sequence2: &[FastaBase],
                                   sequence1_name: &String,
                                   sequence2_name: &String,
                                   starting_position: Option<(usize, usize)>,
) -> AlignmentResult {
    let mut alignment1: Vec<FastaBase> = Vec::with_capacity(sequence1.len() * 2);
    let mut alignment2: Vec<FastaBase> = Vec::with_capacity(sequence1.len() * 2);
    let mut cigars: Vec<AlignmentTag> = Vec::with_capacity(100);

    let mut _starting_x = sequence1.len();
    let mut _starting_y = sequence2.len();

    if starting_position.is_some() {
        _starting_x = starting_position.unwrap().0;
        _starting_y = starting_position.unwrap().1;
    } else if alignment.is_local {
        //println!("trying to find max for seq {}",String::from_utf8(sequence2.clone()).unwrap());
        let max_value_tuple = find_max_value_3d_array(&alignment.scores).unwrap();
        //println!("done! {}",String::from_utf8(sequence2.clone()).unwrap());
        _starting_x = max_value_tuple.0.x;
        _starting_y = max_value_tuple.0.y;
    }
    let starting_z = [(alignment.scores[[_starting_x, _starting_y, 0]], 0),
        (alignment.scores[[_starting_x, _starting_y, 1]], 1),
        (alignment.scores[[_starting_x, _starting_y, 2]], 2)];

    let mut _starting_z = starting_z.iter().max_by(|x, y| x.0.partial_cmp(&y.0).unwrap()).unwrap().1;
    let score = alignment.scores[[_starting_x, _starting_y, _starting_z]];
    //info!("Score is {} from {},{},{}", score, starting_x, starting_y, starting_z);

    let mut path = Vec::new();

    while _starting_x > 0 && _starting_y > 0 && ((alignment.is_local && alignment.scores[[_starting_x, _starting_y, _starting_z]] != 0.0) || !alignment.is_local) {
        alignment.scores[[_starting_x, _starting_y, 0]] = 0.0;
        alignment.scores[[_starting_x, _starting_y, 1]] = 0.0;
        alignment.scores[[_starting_x, _starting_y, 2]] = 0.0;

        //println!("Position: {},{},{}",_starting_x,_starting_y,_starting_z);
        path.push(AlignmentLocation { x: _starting_x, y: _starting_y });

        let movement_delta = match alignment.traceback[[_starting_x, _starting_y, _starting_z]] {
            AlignmentDirection::Diag(size) => (0, size),
            AlignmentDirection::Up(size) => (1, size),
            AlignmentDirection::Left(size) => (2, size),
            AlignmentDirection::Inv(pos1, pos2, inv_move) => {
                let inversion_alignment = inversion_mapping.unwrap().get(&AlignmentLocation { x: pos2.x, y: pos2.y }).unwrap();
                for p in &inversion_alignment.alignment_result.path {
                    path.push(*p);
                }

                //println!("INVERSOIN adding {} and {}", String::from_utf8(inversion_alignment.alignment_result.alignment_string1.clone()).unwrap(),String::from_utf8(inversion_alignment.alignment_result.alignment_string2.clone()).unwrap());
                let mut a1_rev = inversion_alignment.alignment_result.reference_aligned.clone();
                a1_rev.reverse();
                let mut a2_rev = inversion_alignment.alignment_result.read_aligned.clone();
                a2_rev.reverse();

                //println!("pushing -----> {} and {} ", String::from_utf8(a1_rev.clone()).unwrap(),String::from_utf8(a2_rev.clone()).unwrap());
                alignment1.extend(&a1_rev);
                alignment2.extend(&a2_rev);
                _starting_x = pos1.x - 1;
                _starting_y = pos1.y - 1;
                let matrix_move = match inv_move {
                    InvMove::Diag(_) => { 0 }
                    InvMove::Up(_) => { 1 }
                    InvMove::Left(_) => { 2 }
                };
                cigars.push(AlignmentTag::InversionClose); // it'll get reversed at the end
                cigars.extend(&inversion_alignment.alignment_result.cigar_string);
                cigars.push(AlignmentTag::InversionOpen);
                (matrix_move, 0)
            }
        };

        match _starting_z {
            0 => {
                if movement_delta.1 > 0 { cigars.push(AlignmentTag::MatchMismatch(1)) };
                trace!("PUSH MM");
                for _i in 0..(movement_delta.1) {
                    alignment1.push(*sequence1.get(_starting_x - 1).unwrap());
                    alignment2.push(*sequence2.get(_starting_y - 1).unwrap());
                    _starting_x -= 1;
                    _starting_y -= 1;
                }
            }
            1 => {
                if movement_delta.1 > 0 { cigars.push(AlignmentTag::Del(1)) };
                trace!("PUSH INS");
                for _i in 0..(movement_delta.1) {
                    alignment1.push(*sequence1.get(_starting_x - 1).unwrap());
                    alignment2.push(FASTA_UNSET);
                    _starting_x -= 1;
                }
            }
            2 => {
                if movement_delta.1 > 0 { cigars.push(AlignmentTag::Ins(1)) };
                trace!("PUSH DEL");
                for _i in 0..(movement_delta.1) {
                    alignment1.push(FASTA_UNSET);
                    alignment2.push(*sequence2.get(_starting_y - 1).unwrap());
                    _starting_y -= 1;
                }
            }
            _ => panic!("Unknown matrix")
        }
        trace!("looping");
        _starting_z = movement_delta.0;
    }
    //info!("OUT: {},{},{}", starting_x, starting_y, starting_z);
    while _starting_x > 0 && !alignment.is_local {
        alignment1.push(*sequence1.get(_starting_x - 1).unwrap());
        alignment2.push(FASTA_UNSET);
        _starting_x -= 1;
        cigars.push(AlignmentTag::Del(1));
    }
    while _starting_y > 0 && !alignment.is_local {
        alignment1.push(FASTA_UNSET);
        alignment2.push(*sequence2.get(_starting_y - 1).unwrap());
        _starting_y -= 1;
        cigars.push(AlignmentTag::Ins(1));
    }
    alignment1.reverse();
    alignment2.reverse();
    path.reverse();
    cigars.reverse();

    AlignmentResult {
        reference_name: sequence1_name.clone(),
        read_name: sequence2_name.clone(),
        reference_aligned: alignment1,
        read_aligned: alignment2,
        read_quals: None,
        cigar_string: simplify_cigar_string(&cigars),
        path,
        score,
        reference_start: 0, // we're global
        read_start: 0, // we're global
        bounding_box: None,
    }
}

#[allow(dead_code)]
pub fn pretty_print_3d_matrix(alignment: &Alignment<Ix3>, sequence1: &Vec<u8>, sequence2: &Vec<u8>) {
    println!("DIM: {:?}", alignment.scores.shape());
    print!("     ");
    for i in 0..sequence2.len() + 1 {
        print!("{: >9}", if i > 0 { sequence2[i - 1] as char } else { b'x' as char });
    }
    println!();

    for z in 0..3 {
        for x in 0..sequence1.len() + 1 {
            print!("{:>7?},", if x > 0 { sequence1[x - 1] as char } else { b' ' as char });
            print!("[");
            for y in 0..sequence2.len() + 1 {
                print!("{:>8},", alignment.scores[[x, y, z]]);
            }
            println!("]");
        }
        println!();
    }

    print!("      ");
    for i in 0..sequence2.len() + 1 {
        print!("{: >9}", if i > 0 { sequence2[i - 1] as char } else { b' ' as char });
    }
    println!();

    for z in 0..3 {
        for x in 0..sequence1.len() + 1 {
            print!("{:>7?},", if x > 0 { sequence1[x - 1] as char } else { b' ' as char });
            print!("[");
            for y in 0..sequence2.len() + 1 {
                print!("{:>7},", alignment.traceback[[x, y, z]]);
            }
            println!("]");
        }
        println!();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::fasta_bit_encoding::{reverse_complement};

    fn str_to_fasta_vec(input: &str) -> Vec<FastaBase> {
        FastaBase::from_vec_u8_default_ns(&input.as_bytes().to_vec())
    }

    fn fasta_vec_to_string(input: &Vec<FastaBase>) -> String {
        FastaBase::string(input)
    }

    #[test]
    fn waterman_eggart_affine_test_case_2nds() {
        let reference = str_to_fasta_vec("CCAATCTACTACTGCTTGCAGTAC");
        let test_read = str_to_fasta_vec("AGTCCGAGGGCTACTCTACTGAAC");

        let my_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -9.0,
            special_character_score: 8.0,
            gap_open: -10.0,
            gap_extend: -10.0,
            final_gap_multiplier: 1.0,
        };

        let mut alignment_mat = create_scoring_record_3d(reference.len() + 1, test_read.len() + 1, AlignmentType::Affine, true);
        perform_affine_alignment(&mut alignment_mat, &reference, &test_read, &my_score);
        //pretty_print_3d_matrix(&alignment_mat, &reference, &test_read);

        let results = perform_3d_global_traceback(&mut alignment_mat, None, &reference, &test_read, &"reference_name".to_ascii_uppercase(), &"read_name".to_ascii_uppercase(), None);

        assert_eq!(results.reference_aligned, str_to_fasta_vec("CCAATCTACT"));
        assert_eq!(results.read_aligned, str_to_fasta_vec("CTACTCTACT"));

        clean_and_find_next_best_match_3d(&mut alignment_mat, &reference, &test_read, &my_score, &results);
        //pretty_print_3d_matrix(&alignment_mat, &reference, &test_read);

        let results = perform_3d_global_traceback(&mut alignment_mat, None, &reference, &test_read, &"reference_name".to_ascii_uppercase(), &"read_name".to_ascii_uppercase(), None);
        assert_eq!(results.reference_aligned, FastaBase::from_vec_u8(&"CTACTACTGCT".as_bytes().to_vec()));
        assert_eq!(results.read_aligned, str_to_fasta_vec("CTACT-CTACT"));
    }

    #[test]
    fn waterman_eggart_affine_test_case() {
        let reference = str_to_fasta_vec("CCAATCTACTACTGCTTGCAGTAC");
        let test_read = str_to_fasta_vec("AGTCCGAGGGCTACTCTACTGAAC");

        let my_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -9.0,
            special_character_score: 8.0,
            gap_open: -20.0,
            gap_extend: -10.0,
            final_gap_multiplier: 1.0,
        };

        let mut alignment_mat = create_scoring_record_3d(reference.len() + 1, test_read.len() + 1, AlignmentType::Affine, true);
        perform_affine_alignment(&mut alignment_mat, &reference, &test_read, &my_score);

        let results = perform_3d_global_traceback(&mut alignment_mat, None, &reference, &test_read, &"reference_name".to_ascii_uppercase(), &"read_name".to_ascii_uppercase(), None);
        assert_eq!(results.reference_aligned, str_to_fasta_vec("CCAATCTACT"));
        assert_eq!(results.read_aligned, str_to_fasta_vec("CTACTCTACT"));
    }

    #[test]
    fn affine_special_scoring_test() {
        let reference = str_to_fasta_vec("AAAANAAAA");
        let test_read = str_to_fasta_vec("AAAAAAAA");

        let my_score = AffineScoring {
            match_score: 6.0,
            mismatch_score: -6.0,
            special_character_score: 5.0,
            gap_open: -10.0,
            gap_extend: -10.0,
            final_gap_multiplier: 1.0,
        };

        let mut alignment_mat = create_scoring_record_3d(reference.len() + 1, test_read.len() + 1, AlignmentType::Affine, false);
        perform_affine_alignment(&mut alignment_mat, &reference, &test_read, &my_score);

        let results = perform_3d_global_traceback(&mut alignment_mat, None, &reference, &test_read, &"reference_name".to_ascii_uppercase(), &"read_name".to_ascii_uppercase(), None);
        assert_eq!(results.reference_aligned, str_to_fasta_vec("AAAANAAAA"));
        assert_eq!(results.read_aligned, str_to_fasta_vec("AAAA-AAAA"));
    }

    #[test]
    fn affine_loose_ends() {
        let reference = str_to_fasta_vec("ACGTACGTACGT");
        let test_read = str_to_fasta_vec("ACGTACGTT");

        let my_score = AffineScoring {
            match_score: 6.0,
            mismatch_score: -6.0,
            special_character_score: 5.0,
            gap_open: -10.0,
            gap_extend: -10.0,
            final_gap_multiplier: 1.0,
        };

        let mut alignment_mat = create_scoring_record_3d(reference.len() + 1, test_read.len() + 1, AlignmentType::Affine, false);

        let iterations = 50000;
        use std::time::Instant;
        let now = Instant::now();
        for _i in 0..iterations {
            perform_affine_alignment(&mut alignment_mat, &reference, &test_read, &my_score);
            let _results = perform_3d_global_traceback(&mut alignment_mat, None, &reference, &test_read, &"reference_name".to_ascii_uppercase(), &"read_name".to_ascii_uppercase(), None);
        }
        let elapsed = now.elapsed();
        println!("Elapsed affine: {:.2?}", elapsed);
/*
        let now2 = Instant::now();
        for _i in 0..iterations {
            perform_rust_bio_alignment(&reference, &test_read);
        }
        let elapsed2 = now2.elapsed();

        println!("Elapsed biorust: {:.2?}", elapsed2);*/
    }


    #[test]
    fn affine_special_practical_test() {
        let reference = str_to_fasta_vec("AAAAAAAA############################AGATCGGAAGAGCGTCGTGTAGGGAAAGA");
        let test_read = str_to_fasta_vec("AAAAAAAAAAAAAAAAAAAAAAAAATATCTCGTTTAATTGACTCTGAAATCAAGATCGGAAGAGCGTCGTGTAGGGAAAGA");

        let my_score = AffineScoring {
            match_score: 6.0,
            mismatch_score: -6.0,
            special_character_score: 5.0,
            gap_open: -10.0,
            gap_extend: -10.0,
            final_gap_multiplier: 1.0,
        };

        let mut alignment_mat = create_scoring_record_3d(reference.len() + 1, test_read.len() + 1, AlignmentType::Affine, false);
        perform_affine_alignment(&mut alignment_mat, &reference, &test_read, &my_score);

        let results = perform_3d_global_traceback(&mut alignment_mat, None, &reference, &test_read, &"reference_name".to_ascii_uppercase(), &"read_name".to_ascii_uppercase(), None);
        assert_eq!(results.reference_aligned, str_to_fasta_vec("----------------AAAAAAAA############################AGATCGGAAGAGCGTCGTGTAGGGAAAGA"));
        assert_eq!(results.read_aligned, str_to_fasta_vec("AAAAAAAAAAAAAAAAAAAAAAAAATATCTCGTTTAATTGACTCTGAAATCAAGATCGGAAGAGCGTCGTGTAGGGAAAGA"));
    }


    #[test]
    fn affine_alignment_test() {
        let reference = str_to_fasta_vec("AAAA");
        let test_read = str_to_fasta_vec("AATAA");

        let my_score = AffineScoring {
            match_score: 6.0,
            mismatch_score: -6.0,
            special_character_score: 8.0,
            gap_open: -10.0,
            gap_extend: -10.0,
            final_gap_multiplier: 1.0,
        };

        let mut alignment_mat = create_scoring_record_3d(reference.len() + 1, test_read.len() + 1, AlignmentType::Affine, false);
        perform_affine_alignment(&mut alignment_mat, &reference, &test_read, &my_score);

        let results = perform_3d_global_traceback(&mut alignment_mat, None, &reference, &test_read, &"reference_name".to_ascii_uppercase(), &"read_name".to_ascii_uppercase(), None);
        assert_eq!(results.reference_aligned, str_to_fasta_vec("AA-AA"));
        assert_eq!(results.read_aligned, str_to_fasta_vec("AATAA"));
    }

    #[test]
    fn affine_alignment_test_favor_non_special_characters() {
        let reference = str_to_fasta_vec("TTAAGCAGTGGTATCAACGCAGAGTACGCCTTAGGTTAACTTGCTATTTCTAGCTCTAACCCCACCCACGATTGCCGCCGACCCCCATATAAGAAANNNNNNNNNNNNNNNNNNNNNNNNNNAGAT");
        let test_read = str_to_fasta_vec("TTAAGCAGTGGTATCAACGCAGAGTACGCCTTAGGTTAACTTGCTAGTTCTAGCTCTAACCCCACCAACAAGTTTTTCAACACCTAGCGTGT");

        let my_score = AffineScoring::default_dna();

        let mut alignment_mat = create_scoring_record_3d(reference.len() + 1, test_read.len() + 1, AlignmentType::Affine, false);
        perform_affine_alignment(&mut alignment_mat, &reference, &test_read, &my_score);

        //pretty_print_3d_matrix(&alignment_mat, &FastaBase::vec_u8(&reference), &FastaBase::vec_u8(&test_read));
        let results = perform_3d_global_traceback(&mut alignment_mat, None, &reference, &test_read, &"reference_name".to_ascii_uppercase(), &"read_name".to_ascii_uppercase(), None);
        println!("{}\n{}\n{}", FastaBase::string(results.reference_aligned.as_slice()), FastaBase::string(results.read_aligned.as_slice()),results.score);
        assert_eq!(results.reference_aligned, str_to_fasta_vec("TTAAGCAGTGGTATCAACGCAGAGTACGCCTTAGGTTAACTTGCTATTTCTAGCTCTAACCCCACCCACGATTGCCGCCGACCCCCATATAAGAAANNNNNNNNNNNNNNNNNNNNNNNNNNAGAT"));
        //                                                            TTAAGCAGTGGTATCAACGCAGAGTACGCCTTAGGTTAACTTGCTAGTTCTAGCTCTAACCCCACC----------------------------AACAAGTTTTTCAACACCTAGCGTG------T
        assert_eq!(results.read_aligned, str_to_fasta_vec("TTAAGCAGTGGTATCAACGCAGAGTACGCCTTAGGTTAACTTGCTAGTTCTAGCTCTAACCCCACC----------------------------AACAAGTTTTTCAACACCTAGCGTGT------"));
    }

    #[test]
    fn affine_alignment_cigar_test() {
        let reference = str_to_fasta_vec("AAAA");
        let test_read = str_to_fasta_vec("AATAA");

        let my_score = AffineScoring {
            match_score: 6.0,
            mismatch_score: -6.0,
            special_character_score: 8.0,
            gap_open: -10.0,
            gap_extend: -10.0,
            final_gap_multiplier: 1.0,
        };

        let mut alignment_mat = create_scoring_record_3d(reference.len() + 1, test_read.len() + 1, AlignmentType::Affine, false);
        perform_affine_alignment(&mut alignment_mat, &reference, &test_read, &my_score);

        let results = perform_3d_global_traceback(&mut alignment_mat, None, &reference, &test_read, &"reference_name".to_ascii_uppercase(), &"read_name".to_ascii_uppercase(), None);
        assert_eq!(results.reference_aligned, str_to_fasta_vec("AA-AA"));
        assert_eq!(results.read_aligned, str_to_fasta_vec("AATAA"));

        println!("CIGAR : {:?}", results.cigar_string);
    }

    #[test]
    fn affine_alignment_test2() {
        let reference = str_to_fasta_vec("CCAATCTACTACTGCTTGCA");
        let test_read = reverse_complement(&str_to_fasta_vec("GCCACTCTCGCTGTACTGTG"));
        println!("Aligned {} and {} pre {} --",
                 fasta_vec_to_string(&reference),
                 fasta_vec_to_string(&test_read),
                 fasta_vec_to_string(&str_to_fasta_vec("GCCACTCTCGCTGTACTGTG")));

        let my_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -11.0,
            special_character_score: 8.0,
            gap_open: -15.0,
            gap_extend: -5.0,
            final_gap_multiplier: 1.0,
        };

        let mut alignment_mat = create_scoring_record_3d(reference.len() + 1, test_read.len() + 1, AlignmentType::Affine, true);
        perform_affine_alignment(&mut alignment_mat, &reference, &test_read, &my_score);
        //pretty_print_3d_matrix(&alignment_mat, &reference, &test_read);

        let results = perform_3d_global_traceback(&mut alignment_mat, None, &reference, &test_read, &"reference_name".to_ascii_uppercase(), &"read_name".to_ascii_uppercase(), None);
        println!("Aligned {} and {}; from {} and {}",
                 fasta_vec_to_string(&results.reference_aligned),
                 fasta_vec_to_string(&results.read_aligned),
                 fasta_vec_to_string(&reference),
                 fasta_vec_to_string(&test_read));
        assert_eq!(fasta_vec_to_string(&results.reference_aligned), "TACTGC");
        assert_eq!(fasta_vec_to_string(&results.read_aligned), "TACAGC");
    }

    #[test]
    fn inversion_alignment_setup_test() {
        let reference = str_to_fasta_vec("CCAATCTACTACTGCTTGCA");
        let test_read = reverse_complement(&str_to_fasta_vec("GCCACTCTCGCTGTACTGTG"));
        //let reference = String::from("CCAAT").as_bytes().to_owned();
        //let test_read = reverse_complement(&String::from("CTGTG").as_bytes().to_owned());

        let my_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -11.0,
            special_character_score: 8.0,
            gap_open: -15.0,
            gap_extend: -5.0,
            final_gap_multiplier: 1.0,
        };

        let mut alignment_mat = create_scoring_record_3d(reference.len() + 1, test_read.len() + 1, AlignmentType::Affine, true);
        perform_affine_alignment(&mut alignment_mat, &reference, &test_read, &my_score);

        let results = perform_3d_global_traceback(&mut alignment_mat, None, &reference, &test_read, &"reference_name".to_ascii_uppercase(), &"read_name".to_ascii_uppercase(), None);

        assert_eq!(fasta_vec_to_string(&results.reference_aligned), "TACTGC");
        assert_eq!(fasta_vec_to_string(&results.read_aligned), "TACAGC");
    }

    #[test]
    fn inversion_alignment_test() {
        let reference = str_to_fasta_vec("CCAATCTACTACTGCTTGCA");
        let test_read = str_to_fasta_vec("GCCACTCTCGCTGTACTGTG");

        let my_score = InversionScoring {
            match_score: 10.0,
            mismatch_score: -11.0,
            //special_character_score: 9.0,
            gap_open: -15.0,
            gap_extend: -5.0,
            inversion_penalty: -2.0,
            min_inversion_length: 4,
        };


        let my_aff_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -11.0,
            special_character_score: 8.0,
            gap_open: -15.0,
            gap_extend: -5.0,
            final_gap_multiplier: 1.0,
        };

        let results = inversion_alignment(&reference, &test_read, &"reference_name".to_ascii_uppercase(), &"read_name".to_ascii_uppercase(), &my_score, &my_aff_score, true);

        println!("Aligned {} and {} from {} and {}",
                 fasta_vec_to_string(&results.reference_aligned),
                 fasta_vec_to_string(&results.read_aligned),
                 fasta_vec_to_string(&reference),
                 fasta_vec_to_string(&test_read));
        //                                                                  CCAATCTACtactgcTTG
        //                                                                  CCACTCT-CTCTCGCCTG
        assert_eq!(fasta_vec_to_string(&results.reference_aligned), "CCAATCTACTACTGCTTG");
        assert_eq!(fasta_vec_to_string(&results.read_aligned), "CCACTCT-CTACAGCCTG");
    }

    #[test]
    fn inversion_alignment_global_test() {
        let reference = str_to_fasta_vec("CCAATCTACTACTGCTTGCA");
        let test_read = str_to_fasta_vec("CCGTAGATTTACTGCTTGCA");

        let my_score = InversionScoring {
            match_score: 10.0,
            mismatch_score: -11.0,
            //special_character_score: 10.0,
            gap_open: -15.0,
            gap_extend: -5.0,
            inversion_penalty: -2.0,
            min_inversion_length: 2,
        };


        let my_aff_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -11.0,
            special_character_score: 8.0,
            gap_open: -15.0,
            gap_extend: -5.0,
            final_gap_multiplier: 1.0,
        };

        let results = inversion_alignment(&reference, &test_read, &"reference_name".to_ascii_uppercase(), &"read_name".to_ascii_uppercase(), &my_score, &my_aff_score, false);


        println!("Aligned {} and {} from {} and {}",
                 fasta_vec_to_string(&results.reference_aligned),
                 fasta_vec_to_string(&results.read_aligned),
                 fasta_vec_to_string(&reference),
                 fasta_vec_to_string(&test_read));
        //                                                                  CCAATCTACtactgcTTG
        //                                                                  CCACTCT-CTCTCGCCTG
        assert_eq!(fasta_vec_to_string(&results.reference_aligned), "CCAATCTACTACTGCTTGCA");
        assert_eq!(fasta_vec_to_string(&results.read_aligned), "CCAATCTACTACTGCTTGCA");
    }

    #[test]
    fn inversion_alignment_cigar_test() {
        let reference = str_to_fasta_vec("CCAATCTACTACTGCTTGCA");
        let test_read = str_to_fasta_vec("CCGTAGATTTACTGCTTGCA");

        let my_score = InversionScoring {
            match_score: 10.0,
            mismatch_score: -11.0,
            //special_character_score: 9.0,
            gap_open: -15.0,
            gap_extend: -5.0,
            inversion_penalty: -2.0,
            min_inversion_length: 4,
        };


        let my_aff_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -11.0,
            special_character_score: 8.0,
            gap_open: -15.0,
            gap_extend: -5.0,
            final_gap_multiplier: 1.0,
        };

        let results = inversion_alignment(&reference, &test_read, &"reference_name".to_ascii_uppercase(), &"read_name".to_ascii_uppercase(), &my_score, &my_aff_score, false);

        println!("Aligned {} and {} from {} and {}",
                 fasta_vec_to_string(&results.reference_aligned),
                 fasta_vec_to_string(&results.read_aligned),
                 fasta_vec_to_string(&reference),
                 fasta_vec_to_string(&test_read));

        println!("CIGAR: {:?}", results.cigar_string);
    }

}

