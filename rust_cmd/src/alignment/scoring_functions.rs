use FASTA_N;

/// Trait required to instantiate a Scoring instance
#[allow(dead_code)]
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
#[allow(dead_code)]
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
/// TODO: this was removed for performance reasons
//pub trait AffineScoringFunction {
//    fn match_mismatch(&self, a: &FastaBase, b: &FastaBase) -> f64;
//    fn gap_open(&self) -> f64;
//    fn gap_extend(&self) -> f64;
//    fn final_gap_multiplier(&self) -> f64;
//}

pub struct AffineScoring {
    pub(crate) match_score: f64,
    pub(crate) mismatch_score: f64,
    pub(crate) special_character_score: f64,
    pub(crate) gap_open: f64,
    pub(crate) gap_extend: f64,
    pub(crate) final_gap_multiplier: f64,

}

impl AffineScoring {
    /// this is the default, matching DNAFull from Emboss' WATER
    pub fn default_dna() -> AffineScoring {
        AffineScoring {
            match_score: 5.0,
            mismatch_score: -4.0,
            special_character_score: 4.0,
            gap_open: -10.0,
            gap_extend: -0.5,
            final_gap_multiplier: 0.5,
        }
    }

    /// inverted for a distance metric
    pub fn distance_dna() -> AffineScoring {
        AffineScoring {
            match_score: 0.0,
            mismatch_score: -1.0,
            special_character_score: -1.0,
            gap_open: 0.0, // gap open costs gap_open + gap_extend, we automatically get -1.0 for the start of a gap
            gap_extend: -1.0,
            final_gap_multiplier: 1.0,
        }
    }

    pub fn match_mismatch(&self, bit_a: &u8, bit_b: &u8) -> f64 {
        if *bit_a == FASTA_N || *bit_b == FASTA_N || *bit_a < 58 || *bit_b < 58 { self.special_character_score } else if bit_a == bit_b { self.match_score } else { self.mismatch_score }
    }

    pub fn gap_open(&self) -> f64 {
        self.gap_open
    }

    pub fn gap_extend(&self) -> f64 {
        self.gap_extend
    }

    pub fn final_gap_multiplier(&self) -> f64 { self.final_gap_multiplier }
}


pub struct InversionScoring {
    pub match_score: f64,
    pub mismatch_score: f64,
    pub gap_open: f64,
    pub gap_extend: f64,
    pub inversion_penalty: f64,
    pub min_inversion_length: usize,
}


impl InversionScoring {

    #[allow(dead_code)]
    pub fn default() -> InversionScoring {
        InversionScoring {
            match_score: 9.0,
            mismatch_score: -21.0,
            gap_open: -25.0,
            gap_extend: -1.0,
            inversion_penalty: -40.0,
            min_inversion_length: 20,
        }
    }

    pub fn match_mismatch(&self, a: &u8, b: &u8) -> f64 {
        if a == b { self.match_score } else { self.mismatch_score }
    }
    pub fn gap_open(&self) -> f64 {
        self.gap_open
    }

    pub fn gap_extend(&self) -> f64 {
        self.gap_extend
    }
    pub fn inversion_cost(&self) -> f64 {
        self.inversion_penalty
    }
}

// TODO: BUG - ConvexScoring::gap() ignores the `gap_score` and `gap_extend` fields entirely.
// It only uses `gap_open + log10(length)`. The `gap_score` and `gap_extend` fields are defined
// but never used, which is likely an oversight. Also, `gap(0)` will return `-Infinity`
// because `log10(0.0) = -Infinity`.

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_scoring_match() {
        let scoring = SimpleScoring {
            match_score: 5.0,
            mismatch_score: -4.0,
            gap_score: -2.0,
        };
        assert_eq!(scoring.match_mismatch(&b'A', &b'A'), 5.0);
        assert_eq!(scoring.match_mismatch(&b'A', &b'T'), -4.0);
    }

    #[test]
    fn test_simple_scoring_gap() {
        let scoring = SimpleScoring {
            match_score: 5.0,
            mismatch_score: -4.0,
            gap_score: -2.0,
        };
        assert_eq!(scoring.gap(1), -2.0);
        assert_eq!(scoring.gap(3), -6.0);
        assert_eq!(scoring.gap(0), 0.0);
    }

    #[test]
    fn test_convex_scoring_match() {
        let scoring = ConvexScoring {
            match_score: 5.0,
            mismatch_score: -4.0,
            gap_score: -2.0,
            gap_open: -10.0,
            gap_extend: -1.0,
        };
        assert_eq!(scoring.match_mismatch(&b'A', &b'A'), 5.0);
        assert_eq!(scoring.match_mismatch(&b'A', &b'T'), -4.0);
    }

    #[test]
    fn test_convex_scoring_gap() {
        let scoring = ConvexScoring {
            match_score: 5.0,
            mismatch_score: -4.0,
            gap_score: -2.0,
            gap_open: -10.0,
            gap_extend: -1.0,
        };
        // gap(1) = -10.0 + log10(1) = -10.0 + 0.0 = -10.0
        assert_eq!(scoring.gap(1), -10.0);
        // gap(10) = -10.0 + log10(10) = -10.0 + 1.0 = -9.0
        assert_eq!(scoring.gap(10), -9.0);
    }

    #[test]
    fn test_affine_scoring_default_dna() {
        let scoring = AffineScoring::default_dna();
        assert_eq!(scoring.match_score, 5.0);
        assert_eq!(scoring.mismatch_score, -4.0);
        assert_eq!(scoring.gap_open, -10.0);
        assert_eq!(scoring.gap_extend, -0.5);
        assert_eq!(scoring.final_gap_multiplier, 0.5);
    }

    #[test]
    fn test_affine_scoring_distance_dna() {
        let scoring = AffineScoring::distance_dna();
        assert_eq!(scoring.match_score, 0.0);
        assert_eq!(scoring.mismatch_score, -1.0);
    }

    #[test]
    fn test_affine_scoring_match_mismatch_regular() {
        let scoring = AffineScoring::default_dna();
        assert_eq!(scoring.match_mismatch(&b'A', &b'A'), 5.0);
        assert_eq!(scoring.match_mismatch(&b'A', &b'T'), -4.0);
        assert_eq!(scoring.match_mismatch(&b'G', &b'G'), 5.0);
        assert_eq!(scoring.match_mismatch(&b'C', &b'T'), -4.0);
    }

    #[test]
    fn test_affine_scoring_match_mismatch_n_bases() {
        let scoring = AffineScoring::default_dna();
        // N (FASTA_N = b'N' = 78) should trigger special character score
        assert_eq!(scoring.match_mismatch(&b'N', &b'A'), 4.0);
        assert_eq!(scoring.match_mismatch(&b'A', &b'N'), 4.0);
        assert_eq!(scoring.match_mismatch(&b'N', &b'N'), 4.0);
    }

    #[test]
    fn test_affine_scoring_special_characters() {
        let scoring = AffineScoring::default_dna();
        // Characters with ASCII < 58 (digits, special chars) get special score
        // ASCII '0' = 48, '9' = 57, '#' = 35
        assert_eq!(scoring.match_mismatch(&b'0', &b'A'), 4.0);
        assert_eq!(scoring.match_mismatch(&b'A', &b'#'), 4.0);
        assert_eq!(scoring.match_mismatch(&b'1', &b'2'), 4.0);
    }

    #[test]
    fn test_affine_scoring_gap_accessors() {
        let scoring = AffineScoring::default_dna();
        assert_eq!(scoring.gap_open(), -10.0);
        assert_eq!(scoring.gap_extend(), -0.5);
        assert_eq!(scoring.final_gap_multiplier(), 0.5);
    }

    #[test]
    fn test_inversion_scoring_default() {
        let scoring = InversionScoring::default();
        assert_eq!(scoring.match_score, 9.0);
        assert_eq!(scoring.mismatch_score, -21.0);
        assert_eq!(scoring.gap_open, -25.0);
        assert_eq!(scoring.gap_extend, -1.0);
        assert_eq!(scoring.inversion_penalty, -40.0);
        assert_eq!(scoring.min_inversion_length, 20);
    }

    #[test]
    fn test_inversion_scoring_match_mismatch() {
        let scoring = InversionScoring::default();
        assert_eq!(scoring.match_mismatch(&b'A', &b'A'), 9.0);
        assert_eq!(scoring.match_mismatch(&b'A', &b'T'), -21.0);
    }

    #[test]
    fn test_inversion_scoring_gap_accessors() {
        let scoring = InversionScoring::default();
        assert_eq!(scoring.gap_open(), -25.0);
        assert_eq!(scoring.gap_extend(), -1.0);
        assert_eq!(scoring.inversion_cost(), -40.0);
    }
}