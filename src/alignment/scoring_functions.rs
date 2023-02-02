use crate::alignment::fasta_bit_encoding::{FastaBase, FASTA_N, u8_to_encoding_defaulted_to_N};

/// Trait required to instantiate a Scoring instance
pub trait ScoringFunction {
    fn match_mismatch(&self, a: &u8, b: &u8) -> f64;
    fn gap(&self, length: usize) -> f64;
}

#[allow(dead_code)]
pub struct SimpleScoring {
    pub(crate) match_score: f64,
    pub(crate) mismatch_score: f64,
    pub(crate) gap_score: f64,
}


impl ScoringFunction for SimpleScoring {
    fn match_mismatch(&self, a: &u8, b: &u8) -> f64 {
        if a == b { self.match_score } else { self.mismatch_score }
    }

    fn gap(&self, length: usize) -> f64 {
        self.gap_score * length as f64
    }
}
/// Trait required to instantiate a Scoring instance
pub trait ConvexScoringFunction {
    fn match_mismatch(&self, a: &u8, b: &u8) -> f64;
    fn gap(&self, length: usize) -> f64;
}

#[allow(dead_code)]
pub struct ConvexScoring {
    pub(crate) match_score: f64,
    pub(crate) mismatch_score: f64,
    pub(crate) gap_score: f64,
    pub(crate) gap_open: f64,
    pub(crate) gap_extend: f64,
}


impl ConvexScoringFunction for ConvexScoring {
    fn match_mismatch(&self, a: &u8, b: &u8) -> f64 {
        if a == b { self.match_score } else { self.mismatch_score }
    }

    fn gap(&self, length: usize) -> f64 {
        self.gap_open + f64::log10(length as f64)
    }
}


/// Trait required to instantiate a Scoring instance
pub trait AffineScoringFunction {
    fn match_mismatch(&self, a: &u8, b: &u8) -> f64;
    fn gap_open(&self) -> f64;
    fn gap_extend(&self) -> f64;
    fn final_gap_multiplier(&self) -> f64;
}

pub struct AffineScoring {
    pub(crate) match_score: f64,
    pub(crate) mismatch_score: f64,
    pub(crate) special_character_score: f64,
    pub(crate) gap_open: f64,
    pub(crate) gap_extend: f64,
    pub(crate) final_gap_multiplier: f64,

}

impl AffineScoringFunction for AffineScoring {
    fn match_mismatch(&self, a: &u8, b: &u8) -> f64 {
        let bit_a = u8_to_encoding_defaulted_to_N(a);
        let bit_b = u8_to_encoding_defaulted_to_N(b);
        if bit_a == bit_b && (bit_a .identity(&FASTA_N) || bit_b.identity(&FASTA_N)) { self.special_character_score }
        else if bit_a == bit_b {self.match_score}
        else { self.mismatch_score }
    }

    fn gap_open(&self) -> f64 {
        self.gap_open as f64
    }

    fn gap_extend(&self) -> f64 {
        self.gap_extend as f64
    }

    fn final_gap_multiplier(&self) -> f64 {self.final_gap_multiplier}
}


pub trait InversionScoringFunction {
    fn match_mismatch(&self, a: &u8, b: &u8) -> f64;
    fn gap_open(&self) -> f64;
    fn gap_extend(&self) -> f64;
    fn inversion_cost(&self) -> f64;
}
pub struct InversionScoring {
    pub(crate) match_score: f64,
    pub(crate) mismatch_score: f64,
    pub(crate) special_character_score: f64,
    pub(crate) gap_open: f64,
    pub(crate) gap_extend: f64,
    pub(crate) inversion_penalty: f64,
    pub(crate) min_inversion_length: usize,
}


impl InversionScoringFunction for InversionScoring {

    fn match_mismatch(&self, a: &u8, b: &u8) -> f64 {
        let bit_a = u8_to_encoding_defaulted_to_N(a);
        let bit_b = u8_to_encoding_defaulted_to_N(b);
        if bit_a == bit_b { self.match_score } else { self.mismatch_score }
    }
    fn gap_open(&self) -> f64 {
        self.gap_open as f64
    }

    fn gap_extend(&self) -> f64 {
        self.gap_extend as f64
    }
    fn inversion_cost(&self) -> f64 {
        self.inversion_penalty as f64
    }
}