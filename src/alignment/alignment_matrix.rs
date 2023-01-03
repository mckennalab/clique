use std::collections::HashMap;
use std::fmt;
use std::ops::Add;

use ndarray::Array;
use ndarray::prelude::*;
use num_traits::identities::Zero;
use suffix::SuffixTable;
use log::{info, trace, warn};

use crate::alignment::scoring_functions::*;
use bio::alignment::distance::simd::bounded_levenshtein;

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
    pub alignment_segments: Vec<MatchedPosition>,
}

pub struct AlignmentCigar {
    pub alignment_start: usize,
    pub alignment_tags: Vec<AlignmentTag>,
}

/// Our alignment tags -- we don't support clipping for our global alignment approach(es)
#[derive(Eq, PartialEq, Debug, Copy, Clone, Hash)]
pub enum AlignmentTag {
    MatchMismatch(usize),
    Del(usize),
    Ins(usize),
    InversionOpen,
    InversionClose,
}


pub fn simplify_cigar_string(cigar_tokens: &Vec<AlignmentTag>) -> Vec<AlignmentTag> {
    let mut new_cigar = Vec::new();

    let mut last_token: Option<AlignmentTag> = None; // zero length, so combining won't affect the final cigar string

    for token in cigar_tokens.into_iter() {
        last_token = match token {
            AlignmentTag::MatchMismatch(size) => {
                match last_token {
                    Some(AlignmentTag::MatchMismatch(size_old)) => Some(AlignmentTag::MatchMismatch(size + size_old)),
                    _ => {
                        if last_token != None { new_cigar.push(last_token.unwrap()); }

                        Some(AlignmentTag::MatchMismatch(*size))
                    }
                }
            }
            AlignmentTag::Del(size) => {
                match last_token {
                    Some(AlignmentTag::Del(size_old)) => Some(AlignmentTag::Del(size + size_old)),
                    _ => {
                        if last_token != None { new_cigar.push(last_token.unwrap()); }
                        Some(AlignmentTag::Del(*size))
                    }
                }
            }
            AlignmentTag::Ins(size) => {
                match last_token {
                    Some(AlignmentTag::Ins(size_old)) => Some(AlignmentTag::Ins(size + size_old)),
                    _ => {
                        if last_token != None { new_cigar.push(last_token.unwrap()); }
                        Some(AlignmentTag::Ins(*size))
                    }
                }
            }
            AlignmentTag::InversionOpen => {
                match last_token {
                    _ => {
                        // we don't combine inversions
                        if last_token != None { new_cigar.push(last_token.unwrap()); }
                        Some(AlignmentTag::InversionOpen)
                    }
                }
            }
            AlignmentTag::InversionClose => {
                match last_token {
                    _ => {
                        // we don't combine inversions
                        if last_token != None { new_cigar.push(last_token.unwrap()); }
                        Some(AlignmentTag::InversionClose)
                    }
                }
            }
        }
    }
    if let Some(x) = last_token {
        new_cigar.push(x);
    }
    new_cigar.reverse();
    new_cigar
}

impl fmt::Display for AlignmentTag {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            AlignmentTag::MatchMismatch(size) => write!(f, "{}M", size),
            AlignmentTag::Del(size) => write!(f, "{}D", size),
            AlignmentTag::Ins(size) => write!(f, "{}I", size),
            AlignmentTag::InversionOpen => write!(f, "<"),
            AlignmentTag::InversionClose => write!(f, ">"),
        }
    }
}
#[allow(dead_code)]
#[derive(Eq, PartialEq, Debug, Clone, Hash)]
pub enum AlignmentType {
    SIMPLE,
    AFFINE,
    CONVEX,
    ANCHOREDCONVEX,
    INVERSION,
    ANCHOREDINVERSION,
    ANCHOREDAFFINEINVERSION,
}

#[derive(Eq, PartialEq, Debug, Copy, Clone, Hash)]
pub enum InvMove {
    UP(usize),
    LEFT(usize),
    DIAG(usize),
}

#[derive(Eq, PartialEq, Debug, Copy, Clone, Hash)]
pub enum AlignmentDirection {
    UP(usize),
    LEFT(usize),
    DIAG(usize),
    INV(AlignmentLocation, AlignmentLocation, InvMove),
}

impl fmt::Display for AlignmentDirection {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            AlignmentDirection::UP(_size) => write!(f, "| "),
            AlignmentDirection::LEFT(_size) => write!(f, "<-"),
            AlignmentDirection::DIAG(_size) => write!(f, "\\ "),
            AlignmentDirection::INV(pos1, pos2, jumpTo) => write!(f, "I {:?},{:?},{:?}", pos1, pos2, jumpTo),
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
                    AlignmentDirection::LEFT(_size2) => panic!("adding discordant types: UP AND LEFT"),
                    AlignmentDirection::DIAG(_size2) => panic!("adding discordant types: UP AND DIAG"),
                    AlignmentDirection::INV(_x1, _y1, jumpTo) => panic!("adding discordant types: UP AND INV"),
                }
            }
            AlignmentDirection::DIAG(size) => {
                match rhs {
                    AlignmentDirection::UP(_size2) => panic!("adding discordant types: DIAG AND UP"),
                    AlignmentDirection::LEFT(_size2) => panic!("adding discordant types: DIAG AND LEFT"),
                    AlignmentDirection::DIAG(size2) => AlignmentDirection::DIAG(size + size2),
                    AlignmentDirection::INV(_x1, _y1, jumpTo) => panic!("adding discordant types: DIAG AND INV"),
                }
            }
            AlignmentDirection::LEFT(size) => {
                match rhs {
                    AlignmentDirection::UP(_size2) => panic!("adding discordant types: LEFT AND UP"),
                    AlignmentDirection::LEFT(size2) => AlignmentDirection::LEFT(size + size2),
                    AlignmentDirection::DIAG(_size2) => panic!("adding discordant types: LEFT AND DIAG"),
                    AlignmentDirection::INV(_x1, _y1, jumpTo) => panic!("adding discordant types: LEFT AND INV"),
                }
            }
            AlignmentDirection::INV(_x1, _y1, jumpTo) => {
                match rhs {
                    AlignmentDirection::UP(_size2) => panic!("adding discordant types: INV AND UP"),
                    AlignmentDirection::LEFT(_size2) => panic!("adding discordant types: INV AND LEFT"),
                    AlignmentDirection::DIAG(_size2) => panic!("adding discordant types: INV AND DIAG"),
                    AlignmentDirection::INV(_x1, _y1, jumpTo) => panic!("CANT ADD TWO INVERSIONS"),
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
            AlignmentDirection::INV(_x1, _y1, jumpTo) => false, // any inversion has a size != 0
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


pub fn alignment_distance(str1: &Vec<u8>, str2: &Vec<u8>) -> usize {
    let mut mat = create_scoring_record_2d(str1.len(), str2.len(), AlignmentType::SIMPLE, false);
    perform_unit_alignment(&mut mat,str1, str2);

    let alignment = perform_2d_global_traceback(&mut mat, str1, str2);

    (f64::round(f64::abs(f64::max(str1.len() as f64,str2.len() as f64) - alignment.score))) as usize
}


pub fn create_scoring_record_2d(hint_seq_a_len: usize, hint_seq_b_len: usize, alignment_type: AlignmentType, local_alignment: bool) -> Alignment<Ix2> {
    Alignment {
        scores: Array::<f64, Ix2>::zeros((hint_seq_a_len, hint_seq_b_len).f()),
        traceback: Array::<AlignmentDirection, Ix2>::zeros((hint_seq_a_len, hint_seq_b_len).f()),
        alignment_type,
        is_local: local_alignment,
    }
}



fn perform_unit_alignment(alignment: &mut Alignment<Ix2>,sequence1: &Vec<u8>, sequence2: &Vec<u8>) {
    let gap_cost = -1.0;
    let mismatch_cost = -1.0;
    let match_value = 1.0;

    // first column (going down)
    for x in 0..alignment.scores.shape()[0] {
        alignment.scores[[x, 0]] = if alignment.is_local { 0.0 } else { gap_cost * (x as f64)};
        alignment.traceback[[x, 0]] = AlignmentDirection::UP(1);
    }
    // top row
    for y in 0..alignment.scores.shape()[1] {
        alignment.scores[[0, y]] = if alignment.is_local { 0.0 } else { gap_cost * (y as f64) };
        alignment.traceback[[0, y]] = AlignmentDirection::LEFT(1);
    }

    for x in 1..alignment.scores.shape()[0] {
        for y in 1..alignment.scores.shape()[1] {
            let match_score = if &sequence1[x - 1] == &sequence2[y - 1] {match_value} else {mismatch_cost};
            let comp_values = vec![(if alignment.is_local { 0.0 } else { MAX_NEG_SCORE }, AlignmentDirection::DIAG(1)),
                                   (alignment.scores[[x - 1, y - 1]] + match_score, AlignmentDirection::DIAG(1)),
                                   (alignment.scores[[x, y - 1]] + gap_cost, AlignmentDirection::LEFT(1)),
                                   (alignment.scores[[x - 1, y]] + gap_cost, AlignmentDirection::UP(1))];
            let best_score = comp_values.iter().max_by(|x, y| x.0.partial_cmp(&y.0).unwrap());

            alignment.scores[[x, y]] = best_score.unwrap().0;
            alignment.traceback[[x, y]] = best_score.unwrap().1;
        }
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
#[allow(dead_code)]
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

#[allow(dead_code)]
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

fn update_sub_vector3d(alignment: &mut Alignment<Ix3>,
                       sequence1: &Vec<u8>,
                       sequence2: &Vec<u8>,
                       scoring_function: &dyn AffineScoringFunction,
                       row: usize,
                       column: usize,
                       by_row: bool) -> usize {
    let mut row_pos = if by_row { row + 1 } else { row };
    let mut col_pos = if by_row { column } else { column + 1 };
    let mut still_updating = true;
    let mut update_count = 0;

    while row_pos < alignment.scores.shape()[0] && col_pos < alignment.scores.shape()[1] && still_updating {
        let updates = update_3d_score(alignment, sequence1, sequence2, scoring_function, row_pos, col_pos);

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


#[allow(dead_code)]
fn clean_and_find_next_best_match_2d(alignment: &mut Alignment<Ix2>,
                                     sequence1: &Vec<u8>,
                                     sequence2: &Vec<u8>,
                                     scoring_function: &dyn ScoringFunction,
                                     previous_result: &AlignmentResult) {
    let mut current_row = 0;
    let mut current_col = 0;
    for path_entry in &previous_result.path {
        current_row = path_entry.x;
        current_col = path_entry.y;
        update_sub_vector(alignment, sequence1, sequence2, scoring_function, current_row, current_col, true);
        update_sub_vector(alignment, sequence1, sequence2, scoring_function, current_row, current_col, false);
    }

    let mut _still_updating_rows = true;
    let mut _still_updating_cols = true;
    while (_still_updating_rows || _still_updating_rows) &&
        current_row < alignment.scores.shape()[0] &&
        current_col < alignment.scores.shape()[1] {
        let row_update_count = update_sub_vector(alignment, sequence1, sequence2, scoring_function, current_row, current_col, true);
        let col_update_count = update_sub_vector(alignment, sequence1, sequence2, scoring_function, current_row, current_col, false);
        _still_updating_rows = row_update_count > 0;
        _still_updating_cols = col_update_count > 0;
        current_row += 1;
        current_col += 1;
    }
}


fn clean_and_find_next_best_match_3d(alignment: &mut Alignment<Ix3>,
                                     sequence1: &Vec<u8>,
                                     sequence2: &Vec<u8>,
                                     scoring_function: &dyn AffineScoringFunction,
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

    let mut _still_updating_rows = true;
    let mut _still_updating_cols = true;
    while (_still_updating_rows || _still_updating_rows) &&
        current_row < alignment.scores.shape()[0] &&
        current_col < alignment.scores.shape()[1] {
        let row_update_count = update_sub_vector3d(alignment, sequence1, sequence2, scoring_function, current_row, current_col, true);
        let col_update_count = update_sub_vector3d(alignment, sequence1, sequence2, scoring_function, current_row, current_col, false);
        _still_updating_rows = row_update_count > 0;
        _still_updating_cols = col_update_count > 0;
        current_row += 1;
        current_col += 1;
    }
}


/// Affine matrix dimensions are row,column,dimension, where dim 1 is match, dim 2 is deletion (relative to read, sequence2) and dim 3 is insertion
pub fn perform_affine_alignment(alignment: &mut Alignment<Ix3>,
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
        alignment.scores[[x, 0, 1]] = scoring_function.gap_open() + (x as f64 * scoring_function.gap_extend());
        alignment.traceback[[x, 0, 1]] = AlignmentDirection::UP(1);
        alignment.scores[[x, 0, 2]] = scoring_function.gap_open() + (x as f64 * scoring_function.gap_extend());
        alignment.traceback[[x, 0, 2]] = AlignmentDirection::UP(1);
    }
    // top row
    for y in 1..alignment.scores.shape()[1] {
        alignment.scores[[0, y, 0]] = MAX_NEG_SCORE;
        alignment.traceback[[0, y, 0]] = AlignmentDirection::LEFT(1);
        alignment.scores[[0, y, 1]] = scoring_function.gap_open() + (y as f64 * scoring_function.gap_extend());
        alignment.traceback[[0, y, 1]] = AlignmentDirection::LEFT(1);
        alignment.scores[[0, y, 2]] = scoring_function.gap_open() + (y as f64 * scoring_function.gap_extend());
        alignment.traceback[[0, y, 2]] = AlignmentDirection::LEFT(1);
    }

    for x in 1..alignment.scores.shape()[0] {
        for y in 1..alignment.scores.shape()[1] {
            update_3d_score(alignment, &sequence1, &sequence2, scoring_function, x, y);
        }
    }
}


/// Affine matrix dimensions are row,column,dimension, where dim 1 is match, dim 2 is deletion (relative to read, sequence2) and dim 3 is insertion
fn perform_inversion_aware_alignment(alignment: &mut Alignment<Ix3>,
                                     alignment_inversion: &HashMap<AlignmentLocation, BoundedAlignment>,
                                     sequence1: &Vec<u8>,
                                     sequence2: &Vec<u8>,
                                     scoring_function: &dyn InversionScoringFunction) {

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
        alignment.scores[[x, 0, 1]] = scoring_function.gap_open() + (x as f64 * scoring_function.gap_extend());
        alignment.traceback[[x, 0, 1]] = AlignmentDirection::UP(1);
        alignment.scores[[x, 0, 2]] = scoring_function.gap_open() + (x as f64 * scoring_function.gap_extend());
        alignment.traceback[[x, 0, 2]] = AlignmentDirection::UP(1);
    }
    // top row
    for y in 1..alignment.scores.shape()[1] {
        alignment.scores[[0, y, 0]] = MAX_NEG_SCORE;
        alignment.traceback[[0, y, 0]] = AlignmentDirection::LEFT(1);
        alignment.scores[[0, y, 1]] = scoring_function.gap_open() + (y as f64 * scoring_function.gap_extend());
        alignment.traceback[[0, y, 1]] = AlignmentDirection::LEFT(1);
        alignment.scores[[0, y, 2]] = scoring_function.gap_open() + (y as f64 * scoring_function.gap_extend());
        alignment.traceback[[0, y, 2]] = AlignmentDirection::LEFT(1);
    }

    for x in 1..alignment.scores.shape()[0] {
        for y in 1..alignment.scores.shape()[1] {
            update_inversion_alignment(alignment, alignment_inversion, &sequence1, &sequence2, scoring_function, x, y);
        }
    }
}

fn update_inversion_alignment(alignment: &mut Alignment<Ix3>,
                              alignment_inversion: &HashMap<AlignmentLocation, BoundedAlignment>,
                              sequence1: &Vec<u8>,
                              sequence2: &Vec<u8>,
                              scoring_function: &dyn InversionScoringFunction,
                              x: usize,
                              y: usize) -> (bool, bool, bool) {
    let mut _update_x = false;
    let mut _update_y = false;
    let mut _update_z = false;

    if x == 0 && y == 0 {
        trace!("length 1 {} length 2 {}", sequence1.len(), sequence2.len());
        alignment_inversion.iter().for_each(|x| info!("------------> {:?}",x.0.clone()));
    }

    {
        let match_score = scoring_function.match_mismatch(&sequence1[x - 1], &sequence2[y - 1]);

        let pos = AlignmentLocation { x, y };

        let max_match_mismatch = vec![
            if alignment.is_local {0.0} else {MAX_NEG_SCORE},
            alignment.scores[[x - 1, y - 1, 0]] + match_score,
            if alignment.is_local {match_score} else {MAX_NEG_SCORE},
        ];
        let max_match_mismatch = max_match_mismatch.iter().max_by(|x, y| x.partial_cmp(&y).unwrap()).unwrap();

        let best_match = vec![
            if alignment_inversion.contains_key(&pos.clone()) {
                trace!("pos {} {} {:?}", pos.x, pos.y, alignment_inversion.get(&pos.clone()).unwrap().alignment_result.score);
                let inv = alignment_inversion.get(&pos.clone());
                if let Some(inversion) = inv {
                    let bounding = &inversion.bounding_box.0;
                    let first_path = &inversion.bounding_box.0.clone();
                    let last_path = inversion.bounding_box.1.clone();
                    assert_eq!(&last_path,&pos);
                    trace!("Score {} {} {} {} {} {}",first_path.x, first_path.y, x,y, &inversion.alignment_result.score, alignment.scores[[first_path.x - 1, first_path.y - 1, 0]]);

                    let inv_best_match = vec![
                        (alignment.scores[[first_path.x - 1, first_path.y - 1, 1]] , InvMove::UP(1)),
                        (alignment.scores[[first_path.x - 1, first_path.y - 1, 2]], InvMove::LEFT(1)),
                        (alignment.scores[[first_path.x - 1, first_path.y - 1, 0]], InvMove::DIAG(1)),];
                    trace!("+++++++++ Positions {:?}",&inv_best_match);
                    let inv_best_match = inv_best_match.iter().max_by(|x, y| x.0.partial_cmp(&y.0).unwrap()).unwrap();
                    trace!("+++++++++ Positions {:?} {}",&inv_best_match,&inversion.alignment_result.score);

                    (inversion.alignment_result.score + inv_best_match.0 + scoring_function.inversion_cost(), AlignmentDirection::INV(first_path.clone(), last_path.clone(), inv_best_match.1))
                } else {
                    (MAX_NEG_SCORE, AlignmentDirection::UP(1))
                }
            } else { (MAX_NEG_SCORE, AlignmentDirection::UP(1)) },
            (*max_match_mismatch, AlignmentDirection::DIAG(1)),
            (alignment.scores[[x - 1, y - 1, 1]] + match_score, AlignmentDirection::UP(1)),
            (alignment.scores[[x - 1, y - 1, 2]] + match_score, AlignmentDirection::LEFT(1)),
        ];

        let mut best_match = best_match.iter().max_by(|x, y| x.0.partial_cmp(&y.0).unwrap()).unwrap().clone();

        best_match = match best_match.1 {
            AlignmentDirection::INV(pos1, pos2, jump) => {

                let aln = &alignment_inversion.get(&pos.clone());
                if let Some(align) = aln {
                    let start_pos = &align.bounding_box.0.clone();
                    let end_pos = &align.bounding_box.1.clone();
                    let inv = AlignmentDirection::INV(start_pos.clone(), pos2.clone(), jump);
                    //println!("ADDING INVERSION {},{},{},{} score {}  ---- {},{}", start_pos.x, start_pos.y, x, y, align.alignment_result.score, end_pos.x, end_pos.y);
                    (best_match.0, inv)
                } else {
                    info!("SCORES: {},{},{},{}",alignment.scores[[x - 1, y - 1, 1]] + match_score,alignment.scores[[x - 1, y - 1, 2]] + match_score,max_match_mismatch,best_match.0);
                    panic!("Unable to unwrap alignment for inversion!");
                }
            }
            _ => {
                best_match
            }
        };

        _update_x = alignment.scores[[x, y, 0]] != best_match.0;
        alignment.scores[[x, y, 0]] = best_match.0;
        alignment.traceback[[x, y, 0]] = best_match.1;
    }
    {
        let best_gap_x = vec![
            (alignment.scores[[x - 1, y, 1]] + scoring_function.gap_extend(), AlignmentDirection::UP(1)),
            (alignment.scores[[x - 1, y, 2]] + scoring_function.gap_open() + scoring_function.gap_extend(), AlignmentDirection::LEFT(1)),
            (alignment.scores[[x - 1, y, 0]] + scoring_function.gap_open() + scoring_function.gap_extend(), AlignmentDirection::DIAG(1))];
        let best_gap_x = best_gap_x.iter().max_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
        _update_y = alignment.scores[[x, y, 1]] != best_gap_x.unwrap().0;
        alignment.scores[[x, y, 1]] = best_gap_x.unwrap().0;
        alignment.traceback[[x, y, 1]] = best_gap_x.unwrap().1;
    }
    {
        let best_gap_y = vec![
            (alignment.scores[[x, y - 1, 1]] + scoring_function.gap_open() + scoring_function.gap_extend(), AlignmentDirection::UP(1)),
            (alignment.scores[[x, y - 1, 2]] + scoring_function.gap_extend(), AlignmentDirection::LEFT(1)),
            (alignment.scores[[x, y - 1, 0]] + scoring_function.gap_open() + scoring_function.gap_extend(), AlignmentDirection::DIAG(1))];
        let best_gap_y = best_gap_y.iter().max_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
        _update_z = alignment.scores[[x, y, 2]] != best_gap_y.unwrap().0;
        alignment.scores[[x, y, 2]] = best_gap_y.unwrap().0;
        alignment.traceback[[x, y, 2]] = best_gap_y.unwrap().1;
    }
    (_update_x, _update_y, _update_z)
}

fn update_3d_score(alignment: &mut Alignment<Ix3>, sequence1: &Vec<u8>, sequence2: &Vec<u8>, scoring_function: &dyn AffineScoringFunction, x: usize, y: usize) -> (bool, bool, bool) {
    let mut _update_x = false;
    let mut _update_y = false;
    let mut _update_z = false;
    {
        let match_score = scoring_function.match_mismatch(&sequence1[x - 1], &sequence2[y - 1]);
        let max_match_mismatch = vec![
            if alignment.is_local {0.0} else {MAX_NEG_SCORE},
            alignment.scores[[x - 1, y - 1, 0]] + match_score,
            if alignment.is_local {match_score} else {MAX_NEG_SCORE},
        ];
        let max_match_mismatch = max_match_mismatch.iter().max_by(|x, y| x.partial_cmp(&y).unwrap()).unwrap();
        let best_match = vec![
            (alignment.scores[[x - 1, y - 1, 1]] + match_score, AlignmentDirection::UP(1)),
            (alignment.scores[[x - 1, y - 1, 2]] + match_score, AlignmentDirection::LEFT(1)),
            (*max_match_mismatch, AlignmentDirection::DIAG(1)), ];
        let best_match = best_match.iter().max_by(|x, y| x.0.partial_cmp(&y.0).unwrap());


        _update_x = alignment.scores[[x, y, 0]] != best_match.unwrap().0;
        alignment.scores[[x, y, 0]] = best_match.unwrap().0;
        alignment.traceback[[x, y, 0]] = best_match.unwrap().1;
    }
    {
        let best_gap_x = vec![
            (alignment.scores[[x - 1, y, 1]] + scoring_function.gap_extend(), AlignmentDirection::UP(1)),
            (alignment.scores[[x - 1, y, 2]] + scoring_function.gap_open() + scoring_function.gap_extend(), AlignmentDirection::LEFT(1)),
            (alignment.scores[[x - 1, y, 0]] + scoring_function.gap_open() + scoring_function.gap_extend(), AlignmentDirection::DIAG(1))];
        let best_gap_x = best_gap_x.iter().max_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
        _update_y = alignment.scores[[x, y, 1]] != best_gap_x.unwrap().0;
        alignment.scores[[x, y, 1]] = best_gap_x.unwrap().0;
        alignment.traceback[[x, y, 1]] = best_gap_x.unwrap().1;
    }
    {
        let best_gap_y = vec![
            (alignment.scores[[x, y - 1, 1]] + scoring_function.gap_open() + scoring_function.gap_extend(), AlignmentDirection::UP(1)),
            (alignment.scores[[x, y - 1, 2]] + scoring_function.gap_extend(), AlignmentDirection::LEFT(1)),
            (alignment.scores[[x, y - 1, 0]] + scoring_function.gap_open() + scoring_function.gap_extend(), AlignmentDirection::DIAG(1))];
        let best_gap_y = best_gap_y.iter().max_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
        _update_z = alignment.scores[[x, y, 2]] != best_gap_y.unwrap().0;
        alignment.scores[[x, y, 2]] = best_gap_y.unwrap().0;
        alignment.traceback[[x, y, 2]] = best_gap_y.unwrap().1;
    }
    (_update_x, _update_y, _update_z)
}

#[derive(Hash, std::cmp::Eq, Debug, Copy, Clone)]
pub struct AlignmentLocation {
    pub x: usize,
    pub y: usize,
}


impl PartialEq for AlignmentLocation {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x &&
            self.y == other.y
    }
}


#[derive(Clone)]
pub struct AlignmentResult {
    pub alignment_string1: Vec<u8>,
    pub alignment_string2: Vec<u8>,
    pub cigar_string: Vec<AlignmentTag>,
    pub path: Vec<AlignmentLocation>,
    pub score: f64,
    pub bounding_box : Option<(AlignmentLocation,AlignmentLocation)>,
}

impl AlignmentResult {

    #[allow(dead_code)]
    fn slice_out_inversions(&self) -> Vec<AlignmentResult> {
        let mut alignment_string1: Vec<u8> = Vec::new();
        let mut alignment_string2: Vec<u8> = Vec::new();
        let mut cigar_string: Vec<AlignmentTag> = Vec::new();
        let mut path = Vec::new();
        let mut score: f64 = 0.0;
        let mut return_vec = Vec::new();

        let mut position = 0;

        for cigar in &self.cigar_string {
            match cigar {
                AlignmentTag::MatchMismatch(size) | AlignmentTag::Ins(size) | AlignmentTag::Del(size) => {
                    alignment_string1.extend(self.alignment_string1[position..(position + size)].iter());
                    alignment_string2.extend(self.alignment_string2[position..(position + size)].iter());
                    cigar_string.push(*cigar);
                    path.extend_from_slice(&self.path[position..(position + size)]);
                    score = -1.0;
                    position = position + size;
                }
                AlignmentTag::InversionOpen => {
                    if alignment_string1.len() > 0 || alignment_string2.len() > 0 {
                        return_vec.push(AlignmentResult{
                            alignment_string1,
                            alignment_string2,
                            cigar_string,
                            path,
                            score,
                            bounding_box: None
                        });
                        alignment_string1 = Vec::new();
                        alignment_string2 = Vec::new();
                        cigar_string = Vec::new();
                        path = Vec::new();
                        score = 0.0;
                    }
                }
                AlignmentTag::InversionClose => {
                    if alignment_string1.len() > 0 || alignment_string2.len() > 0 {
                        return_vec.push(AlignmentResult{
                            alignment_string1,
                            alignment_string2,
                            cigar_string,
                            path,
                            score,
                            bounding_box: None
                        });
                        alignment_string1 = Vec::new();
                        alignment_string2 = Vec::new();
                        cigar_string = Vec::new();
                        path = Vec::new();
                        score = 0.0;
                    }
                }
            }
        }
        if alignment_string1.len() > 0 || alignment_string2.len() > 0 {
            return_vec.push(AlignmentResult{
                alignment_string1,
                alignment_string2,
                cigar_string,
                path,
                score,
                bounding_box: None
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
            new_path.push(AlignmentLocation { x: path_obj.x, y: new_y});
        }
        new_path.reverse();
        let bounds = (
            AlignmentLocation{x: new_path[new_path.len()-1].x, y: new_path[0].y},
            AlignmentLocation{x: new_path[0].x, y: new_path[new_path.len()-1].y});
        AlignmentResult {
            alignment_string1: self.alignment_string1.clone(),
            alignment_string2: self.alignment_string2.clone(),
            cigar_string: self.cigar_string.clone(),
            path: new_path,
            score: self.score,
            bounding_box: Some(bounds)
        }
    }
}

#[allow(dead_code)]
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
        Some((AlignmentLocation { x: g_max_row, y: g_max_col }, g_max_val))
    } else {
        None
    }
}

pub fn find_max_value_3d_array(matrix: &Array::<f64, Ix3>) -> Option<(AlignmentLocation, f64)> {
    let mut _g_max_row = 0;
    let mut _g_max_col = 0;
    let mut _g_max_z = 0;
    let mut _g_max_val = MAX_NEG_SCORE;

    for x in 0..matrix.shape()[0] {
        for y in 0..matrix.shape()[1] {
            for z in 0..matrix.shape()[2] {
                if matrix[[x, y, z]] > _g_max_val {
                    _g_max_row = x;
                    _g_max_col = y;
                    _g_max_z = z;
                    _g_max_val = matrix[[x, y, z]];
                }
                // update rule #5 from Waterman and Eggart
                else if matrix[[x, y, z]] == _g_max_val && (x + y) < (_g_max_row + _g_max_col) {
                    _g_max_row = x;
                    _g_max_col = y;
                    _g_max_z = z;
                    _g_max_val = matrix[[x, y, z]];
                }
                // update rule #6 from Waterman and Eggart
                else if matrix[[x, y, z]] == _g_max_val &&
                    (_g_max_row + _g_max_col) == (x + y) &&
                    x < _g_max_row {
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
    bounding_box: (AlignmentLocation,AlignmentLocation)
}

pub(crate) fn inversion_alignment(reference: &Vec<u8>, read: &Vec<u8>, my_score: &InversionScoring, my_aff_score: &AffineScoring, local: bool) -> AlignmentResult {
    let mut alignment_mat = create_scoring_record_3d(reference.len() + 1, read.len() + 1, AlignmentType::AFFINE, local);
    let mut inversion_mat = create_scoring_record_3d(reference.len() + 1, read.len() + 1, AlignmentType::AFFINE, true);

    info!("GGGGGGGGGG --------> ref: {} read: {}", String::from_utf8(reference.clone()).unwrap(), String::from_utf8(read.clone()).unwrap());
    let mut long_enough_hits: HashMap<AlignmentLocation, BoundedAlignment> = HashMap::new();
    let rev_comp_read = &reverse_complement(&read);
    perform_affine_alignment(&mut inversion_mat, &reference, rev_comp_read, my_aff_score);
    //pretty_print_3d_matrix(&inversion_mat, &reference, &reverse_complement(&test_read));
    //info!("Aligning {} and {} ({})",str::from_utf8(reference).unwrap(),str::from_utf8(test_read).unwrap(),str::from_utf8(&reverse_complement(&test_read)).unwrap());

    let mut aligned_inv: Option<AlignmentResult> = Some(perform_3d_global_traceback(&mut inversion_mat, None, &reference, &rev_comp_read, None));

    //println!("-------------------------------Adding INVERSIONS");
    while aligned_inv.is_some() {
        let aligned_inv_local = aligned_inv.unwrap();
        let length = aligned_inv_local.path.len();
        if length > 1 {
            let first_pos = &aligned_inv_local.path[0].clone();
            let position = &aligned_inv_local.path[length - 1].clone();
            //println!("============={},{}============={},{}================",first_pos.x,first_pos.y,position.x,position.y);
            let converted_path = aligned_inv_local.convert_inverted_path(read.len());
            let bounding = converted_path.bounding_box.unwrap().clone();
            //println!("============={},{}============={},{}================",bounding.0.x,bounding.0.y,bounding.1.x,bounding.1.y);
            let bAlignment = BoundedAlignment{ alignment_result: converted_path, bounding_box: bounding };
            let true_position = AlignmentLocation { x: bounding.1.x, y: bounding.1.y };
            aligned_inv = if length >= my_score.min_inversion_length {
                //println!("INVERSION put in buffer: {} with score {} position {} {} truepos = {} {} ",String::from_utf8(aligned_inv_local.alignment_string1.clone()).unwrap(),aligned_inv_local.score, position.x, position.y,true_position.x, true_position.y);
                clean_and_find_next_best_match_3d(&mut inversion_mat, &reference, &rev_comp_read, my_aff_score, &aligned_inv_local);
                long_enough_hits.insert(true_position, bAlignment);
                Some(perform_3d_global_traceback(&mut inversion_mat, None, &reference, &rev_comp_read,  None))
            } else {
                None
            }
        } else {
            aligned_inv = None
        }
    }
    perform_inversion_aware_alignment(&mut alignment_mat, &mut long_enough_hits, &reference, &read, my_score);
    //pretty_print_3d_matrix(&alignment_mat, &reference, &read);
    //pretty_print_3d_matrix(&alignment_mat, &reference, &read);
    let res = perform_3d_global_traceback(&mut alignment_mat, Some(&long_enough_hits), &reference, &read, None);
    res
}

#[allow(dead_code)]
fn perform_2d_global_traceback(alignment: &mut Alignment<Ix2>, sequence1: &Vec<u8>, sequence2: &Vec<u8>) -> AlignmentResult {
    let mut alignment1 = Vec::new();
    let mut alignment2 = Vec::new();
    let mut cigars: Vec<AlignmentTag> = Vec::new();

    let mut starting_x = sequence1.len();
    let mut starting_y = sequence2.len();

    if alignment.is_local {
        let max_value_tuple = find_max_value_2d_array(&alignment.scores).unwrap();
        starting_x = max_value_tuple.0.x;
        starting_y = max_value_tuple.0.y;
    }
    let score = alignment.scores[[starting_x, starting_y]];

    let mut path = Vec::new();

    while (starting_x > 0 && starting_y > 0) &&
        ((alignment.is_local && !(alignment.scores[[starting_x, starting_y]] == 0.0)) || !alignment.is_local)
    {
        // for Waterman-Eggart, make sure we zero out values and record the path
        alignment.scores[[starting_x, starting_y]] = 0.0;
        path.push(AlignmentLocation { x: starting_x, y: starting_y });

        match alignment.traceback[[starting_x, starting_y]] {
            AlignmentDirection::DIAG(size) => {
                for _index in 0..size {
                    alignment1.push(*sequence1.get(starting_x - 1).unwrap());
                    alignment2.push(*sequence2.get(starting_y - 1).unwrap());
                    starting_x -= 1;
                    starting_y -= 1;
                    cigars.push(AlignmentTag::MatchMismatch(1));
                }
            }
            AlignmentDirection::LEFT(size) => {
                for _index in 0..size {
                    alignment1.push(b'-');
                    alignment2.push(*sequence2.get(starting_y - 1).unwrap());
                    starting_y -= 1;
                    cigars.push(AlignmentTag::Ins(1));
                }
            }
            AlignmentDirection::UP(size) => {
                for _index in 0..size {
                    alignment1.push(*sequence1.get(starting_x - 1).unwrap());
                    alignment2.push(b'-');
                    starting_x -= 1;
                    cigars.push(AlignmentTag::Del(1));
                }
            }
            AlignmentDirection::INV(pos1, pos2, invMove) => {}
        }
    }

    while starting_x > 0 && !alignment.is_local {
        alignment1.push(*sequence1.get(starting_x).unwrap());
        alignment2.push(b'-');
        starting_x -= 1;
        cigars.push(AlignmentTag::Del(1));
    }
    while starting_y > 0 && !alignment.is_local {
        alignment1.push(b'-');
        alignment2.push(*sequence2.get(starting_y).unwrap());
        starting_y -= 1;
        cigars.push(AlignmentTag::Ins(1));
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
        bounding_box: None
    }
}

#[allow(dead_code)]
pub fn perform_3d_global_traceback(alignment: &mut Alignment<Ix3>,
                               inversion_mapping: Option<&HashMap<AlignmentLocation, BoundedAlignment>>,
                               sequence1: &Vec<u8>,
                               sequence2: &Vec<u8>,
                               starting_position: Option<(usize, usize)>,
) -> AlignmentResult {
    info!("STARTING {:?}",&starting_position);

    let mut alignment1 = Vec::new();
    let mut alignment2 = Vec::new();
    let mut cigars: Vec<AlignmentTag> = Vec::new();

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
    let starting_z = vec![(alignment.scores[[_starting_x, _starting_y, 0]], 0),
                              (alignment.scores[[_starting_x, _starting_y, 1]], 1),
                              (alignment.scores[[_starting_x, _starting_y, 2]], 2)];

    let mut _starting_z = starting_z.iter().max_by(|x, y| x.0.partial_cmp(&y.0).unwrap()).unwrap().1;
    let score = alignment.scores[[_starting_x, _starting_y, _starting_z]];
    //info!("Score is {} from {},{},{}", score, starting_x, starting_y, starting_z);

    let mut path = Vec::new();

    while (_starting_x > 0 && _starting_y > 0) &&
        ((alignment.is_local && !(alignment.scores[[_starting_x, _starting_y, _starting_z]] == 0.0)) || !alignment.is_local) {
        alignment.scores[[_starting_x, _starting_y, 0]] = 0.0;
        alignment.scores[[_starting_x, _starting_y, 1]] = 0.0;
        alignment.scores[[_starting_x, _starting_y, 2]] = 0.0;

        //println!("Position: {},{},{}",_starting_x,_starting_y,_starting_z);
        path.push(AlignmentLocation { x: _starting_x, y: _starting_y });

        let movement_delta = match alignment.traceback[[_starting_x, _starting_y, _starting_z]] {
            AlignmentDirection::DIAG(size) => (0, size),
            AlignmentDirection::UP(size) => (1, size),
            AlignmentDirection::LEFT(size) => (2, size),
            AlignmentDirection::INV(pos1, pos2, invMove) => {

                //println!("INVERSION ADD: positions {}.{} {:?}",pos2.x,pos2.y, invMove);
                let inversion_alignment= inversion_mapping.unwrap().get(&AlignmentLocation{x: pos2.x, y: pos2.y}).unwrap();
                for p in &inversion_alignment.alignment_result.path {
                    path.push(p.clone());
                }

                //println!("INVERSOIN adding {} and {}", String::from_utf8(inversion_alignment.alignment_result.alignment_string1.clone()).unwrap(),String::from_utf8(inversion_alignment.alignment_result.alignment_string2.clone()).unwrap());
                let mut a1_rev = inversion_alignment.alignment_result.alignment_string1.clone();
                a1_rev.reverse();
                let mut a2_rev = inversion_alignment.alignment_result.alignment_string2.clone();
                a2_rev.reverse();

                //println!("pushing -----> {} and {} ", String::from_utf8(a1_rev.clone()).unwrap(),String::from_utf8(a2_rev.clone()).unwrap());
                alignment1.extend(&a1_rev);
                alignment2.extend(&a2_rev);
                _starting_x = pos1.x - 1;
                _starting_y = pos1.y - 1;
                let matrix_move = match invMove {
                    InvMove::DIAG(_) => {0}
                    InvMove::UP(_) => {1}
                    InvMove::LEFT(_) => {2}
                };
                cigars.push(AlignmentTag::InversionClose); // it'll get reversed at the end
                cigars.extend(&inversion_alignment.alignment_result.cigar_string);
                cigars.push(AlignmentTag::InversionOpen);
                (matrix_move, 0)
            }
        };

        match _starting_z {
            0 => {
                if movement_delta.1 > 0 {cigars.push(AlignmentTag::MatchMismatch(1))};
                trace!("PUSH MM");
                for _i in 0..(movement_delta.1) {
                    alignment1.push(*sequence1.get(_starting_x - 1).unwrap());
                    alignment2.push(*sequence2.get(_starting_y - 1).unwrap());
                    _starting_x -= 1;
                    _starting_y -= 1;
                }
            }
            1 => {
                if movement_delta.1 > 0 {cigars.push(AlignmentTag::Del(1))};
                trace!("PUSH INS");
                for _i in 0..(movement_delta.1) {
                    alignment1.push(*sequence1.get(_starting_x - 1).unwrap());
                    alignment2.push(b'-');
                    _starting_x -= 1;
                }
            }
            2 => {
                if movement_delta.1 > 0 {cigars.push(AlignmentTag::Ins(1))};
                trace!("PUSH DEL");
                for _i in 0..(movement_delta.1) {
                    alignment1.push(b'-');
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
        alignment2.push(b'-');
        _starting_x -= 1;
        cigars.push(AlignmentTag::Del(1));
    }
    while _starting_y > 0 && !alignment.is_local {
        alignment1.push(b'-');
        alignment2.push(*sequence2.get(_starting_y - 1).unwrap());
        _starting_y -= 1;
        cigars.push(AlignmentTag::Ins(1));
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
        bounding_box: None
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


#[allow(dead_code)]
pub fn pretty_print_2d_matrix(alignment: &Alignment<Ix2>, sequence1: &Vec<u8>, sequence2: &Vec<u8>) {
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

pub(crate) fn reverse_complement(bases: &Vec<u8>) -> Vec<u8> {
    let mut new_bases = bases.clone().iter().map(|&b| reverse_base(b)).collect::<Vec<u8>>();
    new_bases.reverse();
    new_bases
}

fn reverse_base(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'a' => b't',
        b'C' => b'G',
        b'c' => b'g',
        b'G' => b'C',
        b'g' => b'c',
        b'T' => b'A',
        b't' => b'a',
        _ => b'N',
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str;
// find_max_value_2d_array(matrix: Array::<f64, Ix2>) -> Option<(AlignmentLocation,f64)>

    #[test]
    fn find_max_value() {
        let test_array: Array::<f64, Ix2> = array![[1.,2.,3.], [4.,5.,6.]];
        let best_loc = find_max_value_2d_array(&test_array);
        assert!(!best_loc.is_none());
        let unwrapped = best_loc.unwrap();
        assert_eq!(unwrapped.1, 6.0);
        assert_eq!(unwrapped.0.x, 1);
        assert_eq!(unwrapped.0.y, 2);
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
        perform_simple_alignment(&mut alignment_mat, &reference, &test_read, &my_score);

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
        perform_simple_alignment(&mut alignment_mat, &reference, &test_read, &my_score);
        pretty_print_2d_matrix(&alignment_mat, &reference, &test_read);

        let results = perform_2d_global_traceback(&mut alignment_mat, &reference, &test_read);
        assert_eq!(str::from_utf8(&results.alignment_string1).unwrap(), "CCAATCTACT");
        assert_eq!(str::from_utf8(&results.alignment_string2).unwrap(), "CTACTCTACT");
    }

    #[test]
    fn waterman_eggart_affine_test_case_2nds() {
        let reference = String::from("CCAATCTACTACTGCTTGCAGTAC").as_bytes().to_owned();
        let test_read = String::from("AGTCCGAGGGCTACTCTACTGAAC").as_bytes().to_owned();

        let my_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -9.0,
            special_character_score: 8.0,
            gap_open: -10.0,
            gap_extend: -10.0,
        };

        let mut alignment_mat = create_scoring_record_3d(reference.len() + 1, test_read.len() + 1, AlignmentType::AFFINE, true);
        perform_affine_alignment(&mut alignment_mat, &reference, &test_read, &my_score);
        //pretty_print_3d_matrix(&alignment_mat, &reference, &test_read);

        let results = perform_3d_global_traceback(&mut alignment_mat, None, &reference, &test_read, None);
        assert_eq!(str::from_utf8(&results.alignment_string1).unwrap(), "CCAATCTACT");
        assert_eq!(str::from_utf8(&results.alignment_string2).unwrap(), "CTACTCTACT");

        clean_and_find_next_best_match_3d(&mut alignment_mat, &reference, &test_read, &my_score, &results);
        //pretty_print_3d_matrix(&alignment_mat, &reference, &test_read);

        let results = perform_3d_global_traceback(&mut alignment_mat, None, &reference, &test_read, None);
        assert_eq!(str::from_utf8(&results.alignment_string1).unwrap(), "CTACTACTGCT");
        assert_eq!(str::from_utf8(&results.alignment_string2).unwrap(), "CTACT-CTACT");
    }

    #[test]
    fn waterman_eggart_affine_test_case() {
        let reference = String::from("CCAATCTACTACTGCTTGCAGTAC").as_bytes().to_owned();
        let test_read = String::from("AGTCCGAGGGCTACTCTACTGAAC").as_bytes().to_owned();

        let my_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -9.0,
            special_character_score: 8.0,
            gap_open: -20.0,
            gap_extend: -10.0,
        };

        let mut alignment_mat = create_scoring_record_3d(reference.len() + 1, test_read.len() + 1, AlignmentType::AFFINE, true);
        perform_affine_alignment(&mut alignment_mat, &reference, &test_read, &my_score);

        let results = perform_3d_global_traceback(&mut alignment_mat, None, &reference, &test_read, None);
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
        perform_simple_alignment(&mut alignment_mat, &reference, &test_read, &my_score);
        //pretty_print_2d_matrix(&alignment_mat, &reference, &test_read);

        let results = perform_2d_global_traceback(&mut alignment_mat, &reference, &test_read);
        assert_eq!(str::from_utf8(&results.alignment_string1).unwrap(), "CCAATCTACT");
        assert_eq!(str::from_utf8(&results.alignment_string2).unwrap(), "CTACTCTACT");

        clean_and_find_next_best_match_2d(&mut alignment_mat, &reference, &test_read, &my_score, &results);


        let results = perform_2d_global_traceback(&mut alignment_mat, &reference, &test_read);
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
            special_character_score: 8.0,
            gap_open: -10.0,
            gap_extend: -10.0,
        };

        let mut alignment_mat = create_scoring_record_3d(reference.len() + 1, test_read.len() + 1, AlignmentType::AFFINE, false);
        perform_affine_alignment(&mut alignment_mat, &reference, &test_read, &my_score);

        let results = perform_3d_global_traceback(&mut alignment_mat, None, &reference, &test_read, None);
        assert_eq!(str::from_utf8(&results.alignment_string1).unwrap(), "AA-AA");
        assert_eq!(str::from_utf8(&results.alignment_string2).unwrap(), "AATAA");
    }

    #[test]
    fn affine_alignment_cigar_test() {
        let reference = String::from("AAAA").as_bytes().to_owned();
        let test_read = String::from("AATAA").as_bytes().to_owned();

        let my_score = AffineScoring {
            match_score: 6.0,
            mismatch_score: -6.0,
            special_character_score: 8.0,
            gap_open: -10.0,
            gap_extend: -10.0,
        };

        let mut alignment_mat = create_scoring_record_3d(reference.len() + 1, test_read.len() + 1, AlignmentType::AFFINE, false);
        perform_affine_alignment(&mut alignment_mat, &reference, &test_read, &my_score);

        let results = perform_3d_global_traceback(&mut alignment_mat, None, &reference, &test_read, None);
        assert_eq!(str::from_utf8(&results.alignment_string1).unwrap(), "AA-AA");
        assert_eq!(str::from_utf8(&results.alignment_string2).unwrap(), "AATAA");

        println!("CIGAR : {:?}", results.cigar_string);
    }

    #[test]
    fn affine_alignment_test2() {
        let reference = String::from("CCAATCTACTACTGCTTGCA").as_bytes().to_owned();
        let test_read = reverse_complement(&String::from("GCCACTCTCGCTGTACTGTG").as_bytes().to_owned());

        let my_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -11.0,
            special_character_score: 8.0,
            gap_open: -15.0,
            gap_extend: -5.0,
        };

        let mut alignment_mat = create_scoring_record_3d(reference.len() + 1, test_read.len() + 1, AlignmentType::AFFINE, true);
        perform_affine_alignment(&mut alignment_mat, &reference, &test_read, &my_score);
        //pretty_print_3d_matrix(&alignment_mat, &reference, &test_read);

        let results = perform_3d_global_traceback(&mut alignment_mat, None, &reference, &test_read, None);
        /*println!("Aligned {} and {} from {} and {}",
                 str::from_utf8(&results.alignment_string1).unwrap(),
                 str::from_utf8(&results.alignment_string2).unwrap(),
                 str::from_utf8(&reference).unwrap(),
                 str::from_utf8(&test_read).unwrap(), );*/
        assert_eq!(str::from_utf8(&results.alignment_string1).unwrap(), "TACTGC");
        assert_eq!(str::from_utf8(&results.alignment_string2).unwrap(), "TACAGC");
    }

    #[test]
    fn inversion_alignment_setup_test() {
        let reference = String::from("CCAATCTACTACTGCTTGCA").as_bytes().to_owned();
        let test_read = reverse_complement(&String::from("GCCACTCTCGCTGTACTGTG").as_bytes().to_owned());
        //let reference = String::from("CCAAT").as_bytes().to_owned();
        //let test_read = reverse_complement(&String::from("CTGTG").as_bytes().to_owned());

        let my_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -11.0,
            special_character_score: 8.0,
            gap_open: -15.0,
            gap_extend: -5.0,
        };

        let mut alignment_mat = create_scoring_record_3d(reference.len() + 1, test_read.len() + 1, AlignmentType::AFFINE, true);
        perform_affine_alignment(&mut alignment_mat, &reference, &test_read, &my_score);
        //pretty_print_3d_matrix(&alignment_mat, &reference, &test_read);

        let results = perform_3d_global_traceback(&mut alignment_mat, None, &reference, &test_read, None);

        assert_eq!(str::from_utf8(&results.alignment_string1).unwrap(), "TACTGC");
        assert_eq!(str::from_utf8(&results.alignment_string2).unwrap(), "TACAGC");
    }

    #[test]
    fn inversion_alignment_test() {
        let reference = String::from("CCAATCTACTACTGCTTGCA").as_bytes().to_owned();
        let test_read = String::from("GCCACTCTCGCTGTACTGTG").as_bytes().to_owned();

        let my_score = InversionScoring {
            match_score: 10.0,
            mismatch_score: -11.0,
            special_character_score: 9.0,
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
        };

        let results = inversion_alignment(&reference,&test_read,&my_score,&my_aff_score,true);

        println!("Aligned {} and {} from {} and {}",
                 str::from_utf8(&results.alignment_string1).unwrap(),
                 str::from_utf8(&results.alignment_string2).unwrap(),
                 str::from_utf8(&reference).unwrap(),
                 str::from_utf8(&test_read).unwrap());
        //                                                                  CCAATCTACtactgcTTG
        //                                                                  CCACTCT-CTCTCGCCTG
        assert_eq!(str::from_utf8(&results.alignment_string1).unwrap(), "CCAATCTACTACTGCTTG");
        assert_eq!(str::from_utf8(&results.alignment_string2).unwrap(), "CCACTCT-CTACAGCCTG");
    }

    #[test]
    fn inversion_alignment_global_test() {
        let reference = String::from("CCAATCTACTACTGCTTGCA").as_bytes().to_owned();
        let test_read = String::from("CCGTAGATTTACTGCTTGCA").as_bytes().to_owned();

        let my_score = InversionScoring {
            match_score: 10.0,
            mismatch_score: -11.0,
            special_character_score: 10.0,
            gap_open: -15.0,
            gap_extend: -5.0,
            inversion_penalty: -2.0,
            min_inversion_length: 2
        };


        let my_aff_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -11.0,
            special_character_score: 8.0,
            gap_open: -15.0,
            gap_extend: -5.0,
        };

        let results = inversion_alignment(&reference,&test_read,&my_score,&my_aff_score,false);


        println!("Aligned {} and {} from {} and {}",
                 str::from_utf8(&results.alignment_string1).unwrap(),
                 str::from_utf8(&results.alignment_string2).unwrap(),
                 str::from_utf8(&reference).unwrap(),
                 str::from_utf8(&test_read).unwrap());
        //                                                                  CCAATCTACtactgcTTG
        //                                                                  CCACTCT-CTCTCGCCTG
        assert_eq!(str::from_utf8(&results.alignment_string1).unwrap(), "CCAATCTACTACTGCTTGCA");
        assert_eq!(str::from_utf8(&results.alignment_string2).unwrap(), "CCAATCTACTACTGCTTGCA");
    }
    #[test]
    fn inversion_alignment_cigar_test() {
        let reference = String::from("CCAATCTACTACTGCTTGCA").as_bytes().to_owned();
        let test_read = String::from("CCGTAGATTTACTGCTTGCA").as_bytes().to_owned();

        let my_score = InversionScoring {
            match_score: 10.0,
            mismatch_score: -11.0,
            special_character_score: 9.0,
            gap_open: -15.0,
            gap_extend: -5.0,
            inversion_penalty: -2.0,
            min_inversion_length: 4
        };


        let my_aff_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -11.0,
            special_character_score: 8.0,
            gap_open: -15.0,
            gap_extend: -5.0,
        };

        let results = inversion_alignment(&reference,&test_read,&my_score,&my_aff_score,false);

        println!("Aligned {} and {} from {} and {}",
                 str::from_utf8(&results.alignment_string1).unwrap(),
                 str::from_utf8(&results.alignment_string2).unwrap(),
                 str::from_utf8(&reference).unwrap(),
                 str::from_utf8(&test_read).unwrap());

        println!("CIGAR: {:?}",results.cigar_string);
    }

}

