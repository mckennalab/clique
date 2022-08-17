use std::fmt;
use std::mem;
use std::ops::Add;
use std::str;

use bio::alignment::pairwise::MatchFunc;
use ndarray::Array;
use ndarray::prelude::*;
use num_traits::identities::Zero;
use suffix::SuffixTable;

use float_cmp::approx_eq;
use priority_queue::PriorityQueue;

use crate::alignment::scoring_functions::*;

// used by some alignments to quickly match and extend kmer anchors
pub struct SuffixTableLookup<'s, 't> {
    pub suffix_table: SuffixTable<'s, 't>,
    pub seed_size: usize,
}

pub const MAX_NEG_SCORE: f64 = -100000.0;

/// Find the suffix array 'seeds' given a reference sequence
///
/// # Arguments
///
/// * `name` - a u8 Vec representing the reference sequence
/// * `seed_size` - not used in the suffix array creation, but tracked for each analysis
pub fn find_seeds(reference: &Vec<u8>, seed_size: usize) -> SuffixTableLookup {
    return SuffixTableLookup { suffix_table: SuffixTable::new(String::from_utf8(reference.clone()).unwrap()), seed_size };
}


#[derive(Copy, Clone, Debug)]
pub struct MatchedPosition {
    pub search_start: usize,
    pub ref_start: usize,
    pub length: usize,
}

/// Where our alignment starts, how the sequences align, and the sequences themselves.
pub struct SharedSegments {
    pub start_position: usize,
    pub alignment_cigar: Vec<MatchedPosition>,
}

pub struct AlignmentCigar {
    pub alignment_start: usize,
    pub alignment_tags: Vec<AlignmentTags>,
}

/// Our alignment tags -- we don't support clipping for our global alignment approach(es)
#[derive(Eq, PartialEq, Debug, Copy, Clone, Hash)]
pub enum AlignmentTags {
    MatchMismatch(usize),
    Del(usize),
    Ins(usize),
    Inversion(usize),
}

pub fn simplify_cigar_string(cigar_tokens: &Vec<AlignmentTags>) -> Vec<AlignmentTags> {
    let mut new_cigar = Vec::new();

    let mut last_token: Option<AlignmentTags> = None; // zero length, so combining won't affect the final cigar string

    for token in cigar_tokens.into_iter() {
        last_token = match token {
            AlignmentTags::MatchMismatch(size) => {
                match last_token {
                    Some(AlignmentTags::MatchMismatch(size_old)) => Some(AlignmentTags::MatchMismatch(size + size_old)),
                    _ => {
                        if last_token != None { new_cigar.push(last_token.unwrap()); }

                        Some(AlignmentTags::MatchMismatch(*size))
                    }
                }
            }
            AlignmentTags::Del(size) => {
                match last_token {
                    Some(AlignmentTags::Del(size_old)) => Some(AlignmentTags::Del(size + size_old)),
                    _ => {
                        if last_token != None { new_cigar.push(last_token.unwrap()); }
                        Some(AlignmentTags::Del(*size))
                    }
                }
            }
            AlignmentTags::Ins(size) => {
                match last_token {
                    Some(AlignmentTags::Ins(size_old)) => Some(AlignmentTags::Ins(size + size_old)),
                    _ => {
                        if last_token != None { new_cigar.push(last_token.unwrap()); }
                        Some(AlignmentTags::Ins(*size))
                    }
                }
            }
            AlignmentTags::Inversion(size) => {
                match last_token {
                    _ => {
                        // we don't combine inversions
                        if last_token != None { new_cigar.push(last_token.unwrap()); }
                        Some(AlignmentTags::Ins(*size))
                    }
                }
            }
        }
    }
    new_cigar
}

impl fmt::Display for AlignmentTags {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            AlignmentTags::MatchMismatch(size) => write!(f, "{}M", size),
            AlignmentTags::Del(size) => write!(f, "{}D", size),
            AlignmentTags::Ins(size) => write!(f, "{}I", size),
            AlignmentTags::Inversion(size) => write!(f, "{}F", size),
        }
    }
}

#[derive(Eq, PartialEq, Debug, Copy, Clone, Hash)]
pub enum AlignmentType {
    SIMPLE,
    AFFINE,
    CONVEX,
    ANCHORED_CONVEX,
    INVERSION,
    ANCHORED_INVERSION,
    ANCHORED_AFFINE_INVERSION,
}

#[derive(Eq, PartialEq, Debug, Copy, Clone, Hash)]
pub enum AlignmentDirection {
    UP(usize),
    LEFT(usize),
    DIAG(usize),
}

impl fmt::Display for AlignmentDirection {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            AlignmentDirection::UP(size) => write!(f, "| "),
            AlignmentDirection::LEFT(size) => write!(f, "<-"),
            AlignmentDirection::DIAG(size) => write!(f, "\\ "),
        }
    }
}

impl Add for AlignmentDirection {
    type Output = AlignmentDirection;

    fn add(self, rhs: Self) -> Self::Output {
        match self {
            AlignmentDirection::UP(size) => {
                match rhs {
                    AlignmentDirection::UP(size2) => AlignmentDirection::UP(size + size2),
                    AlignmentDirection::LEFT(size2) => panic!("adding discordant types: UP AND LEFT"),
                    AlignmentDirection::DIAG(size2) => panic!("adding discordant types: UP AND DIAG"),
                }
            }
            AlignmentDirection::DIAG(size) => {
                match rhs {
                    AlignmentDirection::UP(size2) => panic!("adding discordant types: DIAG AND UP"),
                    AlignmentDirection::LEFT(size2) => panic!("adding discordant types: DIAG AND LEFT"),
                    AlignmentDirection::DIAG(size2) => AlignmentDirection::DIAG(size + size2),
                }
            }
            AlignmentDirection::LEFT(size) => {
                match rhs {
                    AlignmentDirection::UP(size2) => panic!("adding discordant types: LEFT AND UP"),
                    AlignmentDirection::LEFT(size2) => AlignmentDirection::LEFT(size + size2),
                    AlignmentDirection::DIAG(size2) => panic!("adding discordant types: LEFT AND DIAG"),
                }
            }
        }
    }
}

impl Zero for AlignmentDirection {
    fn zero() -> Self {
        AlignmentDirection::UP(0)
    }

    fn is_zero(&self) -> bool {
        match self {
            AlignmentDirection::UP(size) => *size == 0 as usize,
            AlignmentDirection::LEFT(size) => *size == 0 as usize,
            AlignmentDirection::DIAG(size) => *size == 0 as usize,
        }
    }
}


pub struct Alignment<K> {
    pub scores: Array::<f64, K>,
    pub traceback: Array::<AlignmentDirection, K>,
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

pub fn create_scoring_record_2d(hint_seq_a_len: usize, hint_seq_b_len: usize, alignment_type: AlignmentType, local_alignment: bool) -> Alignment<Ix2> {
    Alignment {
        scores: Array::<f64, Ix2>::zeros((hint_seq_a_len, hint_seq_b_len).f()),
        traceback: Array::<AlignmentDirection, Ix2>::zeros((hint_seq_a_len, hint_seq_b_len).f()),
        alignment_type,
        is_local: local_alignment,
    }
}

/// Performs alignment inline on the passed alignment structure
/// our assumed structure is:
///    Seq2 (columns) ->
/// S  [    matrix    ]
/// e  [    matrix    ]
/// q  [    matrix    ]
/// 1  [    matrix    ]
/// |  [    matrix    ]
/// v  [    matrix    ]
///
/// so addressing matrix[1,2] would be the alignment of first base in sequence 1 to the second base in sequence 2
///
///
fn perform_simple_alignment(alignment: &mut Alignment<Ix2>, sequence1: &Vec<u8>, sequence2: &Vec<u8>, scoring_function: &dyn ScoringFunction) {
    assert_eq!(alignment.alignment_type, AlignmentType::SIMPLE);
    assert_eq!(alignment.scores.shape()[0], sequence1.len() + 1);
    assert_eq!(alignment.scores.shape()[1], sequence2.len() + 1);

    // first column (going down)
    for x in 0..alignment.scores.shape()[0] {
        alignment.scores[[x, 0]] = if alignment.is_local { 0.0 } else { scoring_function.gap(x) };
        alignment.traceback[[x, 0]] = AlignmentDirection::UP(1);
    }
    // top row
    for y in 0..alignment.scores.shape()[1] {
        alignment.scores[[0, y]] = if alignment.is_local { 0.0 } else { scoring_function.gap(y) };
        alignment.traceback[[0, y]] = AlignmentDirection::LEFT(1);
    }

    for x in 1..alignment.scores.shape()[0] {
        for y in 1..alignment.scores.shape()[1] {
            let match_score = scoring_function.match_mismatch(&sequence1[x - 1], &sequence2[y - 1]);
            let comp_values = vec![(if alignment.is_local { 0.0 } else { MAX_NEG_SCORE }, AlignmentDirection::DIAG(1)),
                                   (alignment.scores[[x - 1, y - 1]] + match_score, AlignmentDirection::DIAG(1)),
                                   (alignment.scores[[x, y - 1]] + scoring_function.gap(1), AlignmentDirection::LEFT(1)),
                                   (alignment.scores[[x - 1, y]] + scoring_function.gap(1), AlignmentDirection::UP(1))];
            let best_score = comp_values.iter().max_by(|x, y| x.0.partial_cmp(&y.0).unwrap());

            alignment.scores[[x, y]] = best_score.unwrap().0;
            alignment.traceback[[x, y]] = best_score.unwrap().1;
        }
    }
}

fn update_sub_vector(alignment: &mut Alignment<Ix2>,
                     sequence1: &Vec<u8>,
                     sequence2: &Vec<u8>,
                     scoring_function: &dyn ScoringFunction,
                     row: usize,
                     column: usize,
                     by_row: bool) -> usize {
    let mut row_pos = if by_row { row + 1 } else { row };
    let mut col_pos = if by_row { column } else { column + 1 };
    let mut still_updating = true;
    let mut update_count = 0;

    while row_pos < alignment.scores.shape()[0] && col_pos < alignment.scores.shape()[1] && still_updating {
        let match_score = scoring_function.match_mismatch(&sequence1[row_pos - 1], &sequence2[col_pos - 1]);
        let comp_values = vec![(alignment.scores[[row_pos - 1, col_pos - 1]] + match_score, AlignmentDirection::DIAG(1)),
                               (alignment.scores[[row_pos, col_pos - 1]] + scoring_function.gap(1), AlignmentDirection::LEFT(1)),
                               (alignment.scores[[row_pos - 1, col_pos]] + scoring_function.gap(1), AlignmentDirection::UP(1))];

        let best_score = comp_values.iter().max_by(|x, y| x.0.partial_cmp(&y.0).unwrap()).unwrap();

        if best_score.0 != alignment.scores[[row_pos, col_pos]] {
            alignment.scores[[row_pos, col_pos]] = best_score.0;
            alignment.traceback[[row_pos, col_pos]] = best_score.1;
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

fn clean_and_find_next_best_match_2d(alignment: &mut Alignment<Ix2>,
                                     sequence1: &Vec<u8>,
                                     sequence2: &Vec<u8>,
                                     scoring_function: &dyn ScoringFunction,
                                     previous_result: &AlignmentResult) {
    let mut update_count = 0;

    let mut current_row = 0;
    let mut current_col = 0;
    for path_entry in &previous_result.path {
        current_row = path_entry.positions[0];
        current_col = path_entry.positions[1];
        let row_update_count = update_sub_vector(alignment, sequence1, sequence2, scoring_function, current_row, current_col, true);
        let col_update_count = update_sub_vector(alignment, sequence1, sequence2, scoring_function, current_row, current_col, false);
        update_count += row_update_count + col_update_count;
    }

    let mut still_updating_rows = true;
    let mut still_updating_cols = true;
    while (still_updating_rows || still_updating_rows) &&
        current_row < alignment.scores.shape()[0] &&
        current_col < alignment.scores.shape()[1] {
        let row_update_count = update_sub_vector(alignment, sequence1, sequence2, scoring_function, current_row, current_col, true);
        let col_update_count = update_sub_vector(alignment, sequence1, sequence2, scoring_function, current_row, current_col, false);
        update_count += row_update_count + col_update_count;
        still_updating_rows = row_update_count > 0;
        still_updating_cols = col_update_count > 0;
        current_row += 1;
        current_col += 1;
    }

    println!("UPDATED {}", update_count);
}


/// Affine matrix dimensions are row,column,dimension, where dim 1 is match, dim 2 is deletion (relative to read, sequence2) and dim 3 is insertion
fn perform_affine_alignment(alignment: &mut Alignment<Ix3>,
                            sequence1: &Vec<u8>,
                            sequence2: &Vec<u8>,
                            scoring_function: &dyn AffineScoringFunction) {
    assert_eq!(alignment.scores.shape()[2], 3);
    assert_eq!(alignment.scores.shape()[0], sequence1.len() + 1);
    assert_eq!(alignment.scores.shape()[1], sequence2.len() + 1);

    alignment.scores[[0, 0, 0]] = 0.0;
    alignment.scores[[0, 0, 1]] = MAX_NEG_SCORE;
    alignment.scores[[0, 0, 2]] = MAX_NEG_SCORE;

    // first column (going down)
    for x in 1..alignment.scores.shape()[0] {
        alignment.scores[[x, 0, 0]] = MAX_NEG_SCORE;
        alignment.traceback[[x, 0, 0]] = AlignmentDirection::UP(1);
        alignment.scores[[x, 0, 1]] = if alignment.is_local { 0.0 } else { scoring_function.gap_open() + (x as f64 * scoring_function.gap_extend())};
        alignment.traceback[[x, 0, 1]] = AlignmentDirection::UP(1);
        alignment.scores[[x, 0, 2]] = if alignment.is_local { 0.0 } else { scoring_function.gap_open() + (x as f64 * scoring_function.gap_extend())};
        alignment.traceback[[x, 0, 2]] = AlignmentDirection::UP(1);
    }
    // top row
    for y in 1..alignment.scores.shape()[1] {
        alignment.scores[[0, y, 0]] = MAX_NEG_SCORE;
        alignment.traceback[[0, y, 0]] = AlignmentDirection::LEFT(1);
        alignment.scores[[0, y, 1]] = if alignment.is_local { 0.0 } else { scoring_function.gap_open() + (y as f64 * scoring_function.gap_extend())};
        alignment.traceback[[0, y, 1]] = AlignmentDirection::LEFT(1);
        alignment.scores[[0, y, 2]] = if alignment.is_local { 0.0 } else { scoring_function.gap_open() + (y as f64 * scoring_function.gap_extend())};
        alignment.traceback[[0, y, 2]] = AlignmentDirection::LEFT(1);
    }

    for x in 1..alignment.scores.shape()[0] {
        for y in 1..alignment.scores.shape()[1] {
            {
                let match_score = scoring_function.match_mismatch(&sequence1[x - 1], &sequence2[y - 1]);
                let best_match = vec![
                    (if alignment.is_local { 0.0 } else { MAX_NEG_SCORE }, AlignmentDirection::DIAG(1)),
                    (alignment.scores[[x - 1, y - 1, 1]] + match_score, AlignmentDirection::UP(1)),
                    (alignment.scores[[x - 1, y - 1, 2]] + match_score, AlignmentDirection::LEFT(1)),
                    (alignment.scores[[x - 1, y - 1, 0]] + match_score, AlignmentDirection::DIAG(1)), ];
                let best_match = best_match.iter().max_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
                alignment.scores[[x, y, 0]] = best_match.unwrap().0;
                alignment.traceback[[x, y, 0]] = best_match.unwrap().1;
            }
            {
                let best_gap_x = vec![
                    (alignment.scores[[x - 1, y, 1]] + scoring_function.gap_extend(), AlignmentDirection::UP(1)),
                    (alignment.scores[[x - 1, y, 2]] + scoring_function.gap_open(), AlignmentDirection::LEFT(1)),
                    (alignment.scores[[x - 1, y, 0]] + scoring_function.gap_open(), AlignmentDirection::DIAG(1))];
                let best_gap_x = best_gap_x.iter().max_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
                alignment.scores[[x, y, 1]] = best_gap_x.unwrap().0;
                alignment.traceback[[x, y, 1]] = best_gap_x.unwrap().1;
            }
            {
                let best_gap_y = vec![
                    (alignment.scores[[x, y - 1, 1]] + scoring_function.gap_open(), AlignmentDirection::UP(1)),
                    (alignment.scores[[x, y - 1, 2]] + scoring_function.gap_extend(), AlignmentDirection::LEFT(1)),
                    (alignment.scores[[x, y - 1, 0]] + scoring_function.gap_open(), AlignmentDirection::DIAG(1))];
                let best_gap_y = best_gap_y.iter().max_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
                alignment.scores[[x, y, 2]] = best_gap_y.unwrap().0;
                alignment.traceback[[x, y, 2]] = best_gap_y.unwrap().1;
            }
        }
    }
}

pub struct AlignmentLocation {
    dims: usize,
    positions: Vec<usize>,
}

pub struct AlignmentResult {
    alignment_string1: Vec<u8>,
    alignment_string2: Vec<u8>,
    cigar_string: Vec<AlignmentTags>,
    path: Vec<AlignmentLocation>,
    score: f64,
}

pub fn find_max_value_2d_array(matrix: &Array::<f64, Ix2>) -> Option<(AlignmentLocation, f64)> {
    let mut g_max_row = 0;
    let mut g_max_col = 0;
    let mut g_max_val = MAX_NEG_SCORE;

    for x in 0..matrix.shape()[0] {
        for y in 0..matrix.shape()[1] {
            if matrix[[x, y]] > g_max_val {
                g_max_row = x;
                g_max_col = y;
                g_max_val = matrix[[x, y]];
            }
            // update rule #5 from Waterman and Eggart
            else if matrix[[x, y]] == g_max_val && (x + y) < (g_max_row + g_max_col) {
                g_max_row = x;
                g_max_col = y;
                g_max_val = matrix[[x, y]];
            }
            // update rule #6 from Waterman and Eggart
            else if matrix[[x, y]] == g_max_val &&
                (g_max_row + g_max_col) == (x + y) &&
                x < g_max_row {
                g_max_row = x;
                g_max_col = y;
                g_max_val = matrix[[x, y]];
            }
        }
    }
    if g_max_val > MAX_NEG_SCORE {
        Some((AlignmentLocation { dims: 2, positions: vec![g_max_row, g_max_col] }, g_max_val))
    } else {
        None
    }
}

pub fn find_max_value_3d_array(matrix: &Array::<f64, Ix3>) -> Option<(AlignmentLocation, f64)> {
    let mut g_max_row = 0;
    let mut g_max_col = 0;
    let mut g_max_z = 0;
    let mut g_max_val = MAX_NEG_SCORE;

    for x in 0..matrix.shape()[0] {
        for y in 0..matrix.shape()[1] {
            for z in 0..matrix.shape()[2] {
                if matrix[[x, y, z]] > g_max_val {
                    g_max_row = x;
                    g_max_col = y;
                    g_max_z = z;
                    g_max_val = matrix[[x, y, z]];
                }
                // update rule #5 from Waterman and Eggart
                else if matrix[[x, y, z]] == g_max_val && (x + y) < (g_max_row + g_max_col) {
                    g_max_row = x;
                    g_max_col = y;
                    g_max_z = z;
                    g_max_val = matrix[[x, y, z]];
                }
                // update rule #6 from Waterman and Eggart
                else if matrix[[x, y, z]] == g_max_val &&
                    (g_max_row + g_max_col) == (x + y) &&
                    x < g_max_row {
                    g_max_row = x;
                    g_max_col = y;
                    g_max_z = z;
                    g_max_val = matrix[[x, y, z]];
                }
            }
        }
    }
    if g_max_val > MAX_NEG_SCORE {
        Some((AlignmentLocation { dims: 2, positions: vec![g_max_row, g_max_col, g_max_z] }, g_max_val))
    } else {
        None
    }
}


fn perform_2d_global_traceback(alignment: &mut Alignment<Ix2>, sequence1: &Vec<u8>, sequence2: &Vec<u8>) -> AlignmentResult {
    let mut alignment1 = Vec::new();
    let mut alignment2 = Vec::new();
    let mut cigars: Vec<AlignmentTags> = Vec::new();

    let mut starting_x = sequence1.len();
    let mut starting_y = sequence2.len();

    if alignment.is_local {
        let max_value_tuple = find_max_value_2d_array(&alignment.scores).unwrap();
        starting_x = max_value_tuple.0.positions[0];
        starting_y = max_value_tuple.0.positions[1];
    }
    let score = alignment.scores[[starting_x, starting_y]];

    let mut path = Vec::new();

    while (starting_x > 0 && starting_y > 0) &&
        ((alignment.is_local && !(alignment.scores[[starting_x, starting_y]] == 0.0)) || !alignment.is_local)
    {
        // for Waterman-Eggart, make sure we zero out values and record the path
        alignment.scores[[starting_x, starting_y]] = 0.0;
        path.push(AlignmentLocation { dims: 2, positions: vec![starting_x, starting_y] });

        match alignment.traceback[[starting_x, starting_y]] {
            AlignmentDirection::DIAG(size) => {
                for index in 0..size {
                    alignment1.push(*sequence1.get(starting_x - 1).unwrap());
                    alignment2.push(*sequence2.get(starting_y - 1).unwrap());
                    starting_x -= 1;
                    starting_y -= 1;
                    cigars.push(AlignmentTags::MatchMismatch(1));
                }
            }
            AlignmentDirection::LEFT(size) => {
                for index in 0..size {
                    alignment1.push(b'-');
                    alignment2.push(*sequence2.get(starting_y - 1).unwrap());
                    starting_y -= 1;
                    cigars.push(AlignmentTags::Ins(1));
                }
            }
            AlignmentDirection::UP(size) => {
                for index in 0..size {
                    alignment1.push(*sequence1.get(starting_x - 1).unwrap());
                    alignment2.push(b'-');
                    starting_x -= 1;
                    cigars.push(AlignmentTags::Del(1));
                }
            }
        }
    }

    while starting_x > 0 && !alignment.is_local {
        alignment1.push(*sequence1.get(starting_x).unwrap());
        alignment2.push(b'-');
        starting_x -= 1;
        cigars.push(AlignmentTags::Del(1));
    }
    while starting_y > 0 && !alignment.is_local {
        alignment1.push(b'-');
        alignment2.push(*sequence2.get(starting_y).unwrap());
        starting_y -= 1;
        cigars.push(AlignmentTags::Ins(1));
    }

    alignment1.reverse();
    alignment2.reverse();
    path.reverse();

    AlignmentResult {
        alignment_string1: alignment1,
        alignment_string2: alignment2,
        cigar_string: simplify_cigar_string(&cigars),
        path,
        score,
    }
}


fn perform_3d_global_traceback(alignment: &mut Alignment<Ix3>, sequence1: &Vec<u8>, sequence2: &Vec<u8>, stop_at_zero: bool) -> AlignmentResult {
    let mut alignment1 = Vec::new();
    let mut alignment2 = Vec::new();
    let mut cigars: Vec<AlignmentTags> = Vec::new();

    let mut starting_x = sequence1.len();
    let mut starting_y = sequence2.len();
    let mut starting_z = vec![(alignment.scores[[starting_x, starting_y, 0]], 0),
                              (alignment.scores[[starting_x, starting_y, 1]], 1),
                              (alignment.scores[[starting_x, starting_y, 2]], 2)];

    if alignment.is_local {
        let max_value_tuple = find_max_value_3d_array(&alignment.scores).unwrap();
        starting_x = max_value_tuple.0.positions[0];
        starting_y = max_value_tuple.0.positions[1];
    }

    let mut starting_z = starting_z.iter().max_by(|x, y| x.0.partial_cmp(&y.0).unwrap()).unwrap().1;
    let score = alignment.scores[[starting_x, starting_y, starting_z]];
    let mut path = Vec::new();


    while (starting_x > 0 && starting_y > 0) &&
        ((alignment.is_local && !(alignment.scores[[starting_x, starting_y, starting_z]] == 0.0)) || !alignment.is_local) {
        alignment.scores[[starting_x, starting_y, starting_z]] = 0.0;
        path.push(AlignmentLocation { dims: 3, positions: vec![starting_x, starting_y, starting_z] });

        let movement_delta = match alignment.traceback[[starting_x, starting_y, starting_z]] {
            AlignmentDirection::DIAG(size) => (0, size),
            AlignmentDirection::UP(size) => (1, size),
            AlignmentDirection::LEFT(size) => (2, size),
        };

        match starting_z {
            0 => {
                cigars.push(AlignmentTags::MatchMismatch(1));
                for i in 0..(movement_delta.1) {
                    alignment1.push(*sequence1.get(starting_x - 1).unwrap());
                    alignment2.push(*sequence2.get(starting_y - 1).unwrap());
                    starting_x -= 1;
                    starting_y -= 1;
                }
            }
            1 => {
                cigars.push(AlignmentTags::Ins(1));
                for i in 0..(movement_delta.1) {
                    alignment1.push(*sequence1.get(starting_x - 1).unwrap());
                    alignment2.push(b'-');
                    starting_x -= 1;
                }
                starting_z = movement_delta.0;
            }
            2 => {
                cigars.push(AlignmentTags::Del(1));
                for i in 0..(movement_delta.1) {
                    alignment1.push(b'-');
                    alignment2.push(*sequence2.get(starting_y - 1).unwrap());
                    starting_y -= 1;
                }
            }
            _ => panic!("Unknown matrix")
        }
        starting_z = movement_delta.0;
    }
    println!("OUT: {},{},{}", starting_x, starting_y, starting_z);
    while starting_x > 0 && !alignment.is_local {
        alignment1.push(*sequence1.get(starting_x - 1).unwrap());
        alignment2.push(b'-');
        starting_x -= 1;
        cigars.push(AlignmentTags::Del(1));
    }
    while starting_y > 0 && !alignment.is_local {
        alignment1.push(b'-');
        alignment2.push(*sequence2.get(starting_y - 1).unwrap());
        starting_y -= 1;
        cigars.push(AlignmentTags::Ins(1));
    }
    alignment1.reverse();
    alignment2.reverse();
    path.reverse();

    AlignmentResult {
        alignment_string1: alignment1,
        alignment_string2: alignment2,
        cigar_string: simplify_cigar_string(&cigars),
        path,
        score,
    }
}

pub fn pretty_print_3d_matrix(alignment: &Alignment<Ix3>, sequence1: &Vec<u8>, sequence2: &Vec<u8>) {
    let mut starting_x = sequence1.len();
    let mut starting_y = sequence2.len();
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
                print!("{:>4?},", alignment.traceback[[x, y, z]]);
            }
            println!("]");
        }
        println!();
    }
}


pub fn pretty_print_2d_matrix(alignment: &Alignment<Ix2>, sequence1: &Vec<u8>, sequence2: &Vec<u8>) {
    let mut starting_x = sequence1.len();
    let mut starting_y = sequence2.len();
    println!("DIM: {:?}", alignment.scores.shape());
    print!("     ");
    for i in 0..sequence2.len() + 1 {
        print!("{: >9}", if i > 0 { sequence2[i - 1] as char } else { b'x' as char });
    }
    println!();

    for x in 0..sequence1.len() + 1 {
        print!("{:>7?},", if x > 0 { sequence1[x - 1] as char } else { b' ' as char });
        print!("[");
        for y in 0..sequence2.len() + 1 {
            print!("{:>8},", alignment.scores[[x, y]]);
        }
        println!("]");
    }
    println!();

    print!("      ");
    for i in 0..sequence2.len() + 1 {
        print!("{: >9}", if i > 0 { sequence2[i - 1] as char } else { b' ' as char });
    }
    println!();

    for x in 0..sequence1.len() + 1 {
        print!("{:>7?},", if x > 0 { sequence1[x - 1] as char } else { b' ' as char });
        print!("[");
        for y in 0..sequence2.len() + 1 {
            print!("{:>4?},", alignment.traceback[[x, y]]);
        }
        println!("]");
    }
    println!();
}

#[cfg(test)]
mod tests {
    use super::*;

// find_max_value_2d_array(matrix: Array::<f64, Ix2>) -> Option<(AlignmentLocation,f64)>

    #[test]
    fn find_max_value() {
        let test_array: Array::<f64, Ix2> = array![[1.,2.,3.], [4.,5.,6.]];
        let best_loc = find_max_value_2d_array(&test_array);
        assert!(!best_loc.is_none());
        let unwrapped = best_loc.unwrap();
        assert_eq!(unwrapped.1, 6.0);
        assert_eq!(unwrapped.0.positions[0], 1);
        assert_eq!(unwrapped.0.positions[1], 2);
    }

    #[test]
    fn simple_alignment_test() {
        let reference = String::from("AAAA").as_bytes().to_owned();
        let test_read = String::from("AATAA").as_bytes().to_owned();

        let my_score = SimpleScoring {
            match_score: 6.0,
            mismatch_score: -4.0,
            gap_score: -10.0,
        };

        let mut alignment_mat = create_scoring_record_2d(reference.len() + 1, test_read.len() + 1, AlignmentType::SIMPLE, false);
        let aligned = perform_simple_alignment(&mut alignment_mat, &reference, &test_read, &my_score);

        let results = perform_2d_global_traceback(&mut alignment_mat, &reference, &test_read);
        assert_eq!(str::from_utf8(&results.alignment_string1).unwrap(), "AA-AA");
        assert_eq!(str::from_utf8(&results.alignment_string2).unwrap(), "AATAA");
    }

    #[test]
    fn waterman_eggart_test_case() {
        let reference = String::from("CCAATCTACTACTGCTTGCAGTAC").as_bytes().to_owned();
        let test_read = String::from("AGTCCGAGGGCTACTCTACTGAAC").as_bytes().to_owned();

        let my_score = SimpleScoring {
            match_score: 10.0,
            mismatch_score: -9.0,
            gap_score: -20.0,
        };

        let mut alignment_mat = create_scoring_record_2d(reference.len() + 1, test_read.len() + 1, AlignmentType::SIMPLE, true);
        let aligned = perform_simple_alignment(&mut alignment_mat, &reference, &test_read, &my_score);
        pretty_print_2d_matrix(&alignment_mat, &reference, &test_read);

        let results = perform_2d_global_traceback(&mut alignment_mat, &reference, &test_read);
        assert_eq!(str::from_utf8(&results.alignment_string1).unwrap(), "CCAATCTACT");
        assert_eq!(str::from_utf8(&results.alignment_string2).unwrap(), "CTACTCTACT");
    }

    #[test]
    fn waterman_eggart_affine_test_case() {
        let reference = String::from("CCAATCTACTACTGCTTGCAGTAC").as_bytes().to_owned();
        let test_read = String::from("AGTCCGAGGGCTACTCTACTGAAC").as_bytes().to_owned();

        let my_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -9.0,
            gap_open: -20.0,
            gap_extend: -10.0,
        };

        let mut alignment_mat = create_scoring_record_3d(reference.len() + 1, test_read.len() + 1, AlignmentType::AFFINE, true);
        let aligned = perform_affine_alignment(&mut alignment_mat, &reference, &test_read, &my_score);
        pretty_print_3d_matrix(&alignment_mat, &reference, &test_read);

        let results = perform_3d_global_traceback(&mut alignment_mat, &reference, &test_read, true);
        println!("SSSSSSSSSSSSS = {},{}",str::from_utf8(&results.alignment_string1).unwrap(),str::from_utf8(&results.alignment_string2).unwrap());
        assert_eq!(str::from_utf8(&results.alignment_string1).unwrap(), "CCAATCTACT");
        assert_eq!(str::from_utf8(&results.alignment_string2).unwrap(), "CTACTCTACT");
    }

    #[test]
    fn waterman_eggart_test_case_round2() {
        let reference = String::from("CCAATCTACTACTGCTTGCAGTAC").as_bytes().to_owned();
        let test_read = String::from("AGTCCGAGGGCTACTCTACTGAAC").as_bytes().to_owned();

        let my_score = SimpleScoring {
            match_score: 10.0,
            mismatch_score: -9.0,
            gap_score: -20.0,
        };

        let mut alignment_mat = create_scoring_record_2d(reference.len() + 1, test_read.len() + 1, AlignmentType::SIMPLE, true);
        let aligned = perform_simple_alignment(&mut alignment_mat, &reference, &test_read, &my_score);
        pretty_print_2d_matrix(&alignment_mat, &reference, &test_read);

        let results = perform_2d_global_traceback(&mut alignment_mat, &reference, &test_read);
        pretty_print_2d_matrix(&alignment_mat, &reference, &test_read);
        clean_and_find_next_best_match_2d(&mut alignment_mat, &reference, &test_read, &my_score, &results);
        pretty_print_2d_matrix(&alignment_mat, &reference, &test_read);
        assert_eq!(str::from_utf8(&results.alignment_string1).unwrap(), "CCAATCTACT");
        assert_eq!(str::from_utf8(&results.alignment_string2).unwrap(), "CTACTCTACT");

        let results = perform_2d_global_traceback(&mut alignment_mat, &reference, &test_read);
        pretty_print_2d_matrix(&alignment_mat, &reference, &test_read);
        assert_eq!(str::from_utf8(&results.alignment_string1).unwrap(), "CTACTACTGCT");
        assert_eq!(str::from_utf8(&results.alignment_string2).unwrap(), "CTACT-CTACT");
    }


    #[test]
    fn affine_alignment_test() {
        let reference = String::from("AAAA").as_bytes().to_owned();
        let test_read = String::from("AATAA").as_bytes().to_owned();

        let my_score = AffineScoring {
            match_score: 6.0,
            mismatch_score: -6.0,
            gap_open: -10.0,
            gap_extend: -10.0,
        };

        let mut alignment_mat = create_scoring_record_3d(reference.len() + 1, test_read.len() + 1, AlignmentType::AFFINE, false);
        let aligned = perform_affine_alignment(&mut alignment_mat, &reference, &test_read, &my_score);
        //pretty_print_3d_matrix(&alignment_mat, &reference, &test_read);

        let results = perform_3d_global_traceback(&mut alignment_mat, &reference, &test_read, true);
        assert_eq!(str::from_utf8(&results.alignment_string1).unwrap(), "AA-AA");
        assert_eq!(str::from_utf8(&results.alignment_string2).unwrap(), "AATAA");
    }
    /*
        #[test]
        fn test_basic_align_with_anchors() {
            let reference = String::from("NNNNNNNNCATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCATATGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCGTAGAAGAAAGTTAGCTCTTGTGCGAGCTACAGGAACGATGTTTGATTAGAGTAAGCAGAGGACAAGGGCTCGCGTGCAGCCGAAGTTTGGCCGGTACTCTCCAACCGTTAACAACAACACCTTTCATCGAAATCCGCTTGGTAACAACACTAGGATGTTTCGACGCACTCGATAACCGGGAAACCAAGAGAAGTTTCCGAGCGCCACAGCCCAACTTAACTTCGCCATGTTTGAGACACGCGATCGGGACCACAAGACGTTACTCTTTGGGACCGCCGTAAGCGAGTATAAGCAGATTGTGTTTCGGCGTCCAAGTTGCCGTCAAAAGCTTACTGAGTTTCGCTGCCGGCGGAATGCTATTAATCGCGCCTACTTTCATGGAAGACGTTTCCGCTAAAATGGACTGGTGTTTCATGTCGGGAGCCGCTTTGTAAGATGGAGCTACTTTCCAGTCTGAGTTCATGCGAGAACACCACGAGAGTTTCATGAGTGGCCTCTCGAATCAACAGTCTACAAGTTTGGAGTTATCCGACACATCAAAACCAGCCATTCGTTTCATGAGATGGATCGCATACTAACCTTAACGGAGTTTGTAGTCACGGACGAACGATAAACGAGCATAATCTTTCGAGTCGGAGGCGAGTACTTAACGGATATAACGTTTCGTGCCAATGTTAACCGGTCAACTACCACTCAGTTTCTTGTTCATCATAACACTGAAACTGAGATCGTCTTTGGTGCAATTCCAATACGGCTAACTTACGCATACTTTGATGACGCCGTGATTATATCAAGAACCTACCGCTTTCATGGCGGTAACGGTATCCAAAGAATTGGTGTGTTTCGTGCATGCAGTGTCGGACTAAGACTAGGAATGTTTGCAGTGGCCGCAGTTATTCCAAGAGGATGCTTCTTTCCAGCTAACGGTCGCCTAATAAGATGTAACTGGTTTCTTGAGGAGGCATGTACCCGAAGCTTAGAGTAGTCTCCTCTATGAATACTGAAGGACTTGCGTAGTTATGTACAAGCTCACCAACGGACGGGTGCTTCCACATATAACGTTAGCATCTCGTGTGCTATTCGTAAGAGTTTCTAGTCACGGACGAACGATAAAGTACCAACGCCTTTCATGAGTGGCCTCTCGAATCAAGTGATCGGACCTTTGGACGCACTCGATAACCGGGAAGTTATCCAGACTTTCGTGCCAATGTTAACCGGTCAATAAGAGCTACCTTTGATGACGCCGTGATTATATCAATACGCTTCTGGTTTGGGCGTCCAAGTTGCCGTCAAATAGTAGTGACCTTTGCAGTCTGAGTTCATGCGAGAATCACCGCGAAGTTTGTTGTTCATCATAACACTGAAATCCGCAATTAGTTTCCAGCTAACGGTCGCCTAATAATCGGTAGCACGTTTCCGAGCGCCACAGCCCAACTTATCTTACACACGTTTCATCGAAATCCGCTTGGTAACATGCAAGTGTAGTTTGCAGTGGCCGCAGTTATTCCAATGTGTGTGAGCTTTCGAGTTATCCGACACATCAAACAACCGATTAACTTTCGTGCATGCAGTGTCGGACTACAAGAATAGTGCTTTGATGGCGGTAACGGTATCCAACACACTATTACCTTTCATGTCGGGAGCCGCTTTGTACAGTGTGCTACCTTTGGTGCAATTCCAATACGGCTACATAGCTAACGGTTTCTTGAGGAGGCATGTACCCGACCGCCTTATTCGTTTGATGGAAGACGTTTCCGCTAACCGGTTCTAAGGTTTCGGACCGCCGTAAGCGAGTATCCTAGTACATGGTTTCGCTGCCGGCGGAATGCTATTCCTGCTAACGAGTTTCATGAGATGGATCGCATACTACGAAGCGTCATGTTTGGCCGGTACTCTCCAACCGTTCGAGCTTCTTCCTTTCGAGTCGGAGGCGAGTACTTACGATGATAGAGCTTTCAGACACGCGATCGGGACCACCGCAGCTAACAGTAATAGGACCGACCGACCGTTCGATTCAGGGAGATTGCCCTACACTATGCGGCAGCTGGCATAGACTCCTAAGGAGATGCGTACTTGTTAAATAGGACTCTTTCATCGAAATCCGCTTGGTAACCGCTAGGTTACGTTTGTTGTTCATCATAACACTGAACGTAACTATGTCTTTCGAGTTATCCGACACATCAAACTAAGTATGAGCTTTCCGAGCGCCACAGCCCAACTTCTAGCTAATCTCTTTGGCCGGTACTCTCCAACCGTTCTATTATGCCTGTTTGGCTGCCGGCGGAATGCTATTCTCCTGCTACACTTTGTAGTCACGGACGAACGATAACTGCTTAGAACCTTTGCAGTGGCCGCAGTTATTCCACTTAACGCGGAGTTTGGACGCACTCGATAACCGGGAGAACATTAGCTCTTTCATGAGATGGATCGCATACTAGAAGACATTAGGTTTCATGACGCCGTGATTATATCAGAAGTGTTACGGTTTCGGCGTCCAAGTTGCCGTCAAGAATTGATGAGCTTTCGTGCCAATGTTAACCGGTCAGACATTATAGCCTTTGGTGCAATTCCAATACGGCTAGACTCACACGACTTTGGAGTCGGAGGCGAGTACTTAGAGTTGTTGACGTTTCCAGCTAACGGTCGCCTAATAGATGATAGAACGTTTGCAGTCTGAGTTCATGCGAGAGCAATAAGCTACTTTCGGACCGCCGTAAGCGAGTATGCGATTAAGTAGTTTGATGTCGGGAGCCGCTTTGTAGGATACTCGACGTTTCAGACACGCGATCGGGACCACGGTAATATGACGTTTCTTGAGGAGGCATGTACCCGAGGTGTTGCATGGTTTCATGGAAGACGTTTCCGCTAAGTAACGGTATTCTTTGATGAGTGGCCTCTCGAATCAGTACCAGACTTGTTTCATGGCGGTAACGGTATCCAAGTAGCATAATCGTTTGGTGCATGCAGTGTCGGACTAGTATGCTATCGCAGGAGGATGGGGCAGGACAGGACGCGGCCACCCCAGGCCTTTCCAGAGCAAACCTGGAGAAGATTCACAATAGACAGGCCAAGAAACCCGGTGCTTCCTCCAGAGCCGTTTAAAGCTGATATGAGGAAATAAAGAGTGAACTGGAAGGATCCGATATCGCCACCGTGGCTGAATGAGACTGGTGTCGACCTGTGCCT").as_bytes().to_owned();
            let test_read = String::from("GTATTGCTCATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCAAATGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCATAGAAGAAGCTTAGCTCTTGTGCGAGCTACAGGAACGATGTTTGATTAGAGTAAGCAGAGGACAAGGGCTCGCGTGCAGCCGAAGTTTGGCCGGTACTCTCCATACAGTGTGCTACCTTTGGTGCAATTCCAATACGGCTACATAGCTAACGGTTTCTTGAGGAGGCATGTACCCGACCGCCTTATTCGTTTGATGGAAGACGTTTCCGCTAACCGGTTCTAAGGTTTCGGGACCGCCGTAAGCGATTGATGAGCTTTCGTGCCAATGTTAACCGGTCAGACATTATAGCCTTTGGTGCAATTCCAATACGGTAATATGACGTTTCTTGAGGAGGCATGTACCCGAGGTGTTGCATGGTTTCATGGAAGACGTTTCCGCTAAGTAACGGTATTCTTTGATGAGTGGCCTCTCGAATCAGTACCAGACTTGTTTCATGGCGGTAACGGTATCCAAGTAGCATAATCGTTTGGTGCATGCAGTGTCGGACTAGTATGCTATCGCAGGAGGATGGGGCAGGACAGGACGCGGCCACCCCAGGCCTTTCCAGAGCAAACCTGGAGAAGATTCACAATAGACAGGCCAAGAAACCCGGTGCTTCCTCCAGAGCCGTTTAAAGCTGATATGAGGAAATAAAGAGTGAACTGGAAGGATCCCATATCGACAATACGTAACTGAACGAAGTACACCAGTATT").as_bytes().to_owned();

            let my_score = AffineScoring {
                match_score: 6.0,
                mismatch_score: -6.0,
                gap_open: -10.0,
                gap_extend: -1.0,
            };
            let mut alignment_mat = create_scoring_record_3d(reference.len() + 1, test_read.len() + 1, AlignmentType::AFFINE);
            let aligned = perform_affine_alignment(&mut alignment_mat, &reference, &test_read, &my_score);
            let fwd_score_mp = perform_3d_global_traceback(&mut alignment_mat, &reference, &test_read, true);


        }*/
}

