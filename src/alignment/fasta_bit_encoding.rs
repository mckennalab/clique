extern crate derive_more;

use derive_more::{From, Display, Add};



/// our core Fasta base - representing a single fasta character in a u8 data store. We actually pack
/// everything into a u4 (if that was a type) but FastaString below handles the more dense packing and
/// unpacking into u64 structures
#[derive(Shrinkwrap)]
#[derive(Clone, Copy, From, Display, Add, Debug)]
#[display(fmt = "{}", _0)]
pub struct FastaBase(u8);

//pub struct FastaMask(u8);

/// our comparisons are done using logical AND operations for speed (the layout of bits is really important).
/// Each degenerate base should also be 'equal' to the correct 'ACGT' combination and be equal to other degenerate
/// bases that share any overlap in their base patterns. E.g. R == K, but R != C (N == everything)
impl PartialEq for FastaBase {
    fn eq(&self, other: &Self) -> bool {
        self.0 & other.0 > (0 as u8)
    }
}
impl Eq for FastaBase {}

use std::ops::BitAnd;
impl BitAnd for FastaBase {
    type Output = FastaBase;
    fn bitand(self, rhs: FastaBase) -> FastaBase {
        FastaBase(self.0.bitand(rhs.0))
    }
}
impl BitAnd for &FastaBase {
    type Output = FastaBase;
    fn bitand(self, rhs: &FastaBase) -> FastaBase {
        FastaBase(self.0.bitand(rhs.0))
    }
}

/// TODO: so I couldn't figure out how to create the const version as a standard trait, so I did it here (self note: just look at u8 impl you dummy)
const fn add_two_string_encodings(a: FastaBase, b: FastaBase) -> FastaBase {
    FastaBase((a.0 + b.0) as u8)
}

const FASTA_N: FastaBase = FastaBase(0xF as u8);
const FASTA_A: FastaBase = FastaBase(0x1 as u8);
const FASTA_C: FastaBase = FastaBase(0x2 as u8);
const FASTA_G: FastaBase = FastaBase(0x4 as u8);
const FASTA_T: FastaBase = FastaBase(0x8 as u8);
const FASTA_R: FastaBase = add_two_string_encodings(FASTA_A, FASTA_G);
const FASTA_Y: FastaBase = add_two_string_encodings(FASTA_C, FASTA_T);
const FASTA_K: FastaBase = add_two_string_encodings(FASTA_G, FASTA_T);
const FASTA_M: FastaBase = add_two_string_encodings(FASTA_A, FASTA_C);
const FASTA_S: FastaBase = add_two_string_encodings(FASTA_C, FASTA_G);
const FASTA_W: FastaBase = add_two_string_encodings(FASTA_A, FASTA_T);
const FASTA_B: FastaBase = add_two_string_encodings(FASTA_C, add_two_string_encodings(FASTA_G, FASTA_T));
const FASTA_D: FastaBase = add_two_string_encodings(FASTA_A, add_two_string_encodings(FASTA_G, FASTA_T));
const FASTA_H: FastaBase = add_two_string_encodings(FASTA_A, add_two_string_encodings(FASTA_C, FASTA_T));
const FASTA_V: FastaBase = add_two_string_encodings(FASTA_A, add_two_string_encodings(FASTA_C, FASTA_G));

#[allow(dead_code)]
pub fn char_to_encoding(base: &char) -> Option<FastaBase> {
    match base {
        'A' | 'a' => Some(FASTA_A),
        'C' | 'c' => Some(FASTA_C),
        'G' | 'g' => Some(FASTA_G),
        'T' | 't' => Some(FASTA_T),

        'R' | 'r' => Some(FASTA_R),
        'Y' | 'y' => Some(FASTA_C),
        'K' | 'k' => Some(FASTA_K),
        'M' | 'm' => Some(FASTA_M),

        'S' | 's' => Some(FASTA_S),
        'W' | 'w' => Some(FASTA_W),
        'B' | 'b' => Some(FASTA_B),
        'D' | 'd' => Some(FASTA_D),

        'H' | 'h' => Some(FASTA_H),
        'V' | 'v' => Some(FASTA_V),

        'N' | 'n' => Some(FASTA_N),
        _ => None,
    }
}

#[allow(dead_code)]
pub fn u8_to_encoding(base: &u8) -> Option<FastaBase> {
    match base {
        b'A' | b'a' => Some(FASTA_A),
        b'C' | b'c' => Some(FASTA_C),
        b'G' | b'g' => Some(FASTA_G),
        b'T' | b't' => Some(FASTA_T),

        b'R' | b'r' => Some(FASTA_R),
        b'Y' | b'y' => Some(FASTA_C),
        b'K' | b'k' => Some(FASTA_K),
        b'M' | b'm' => Some(FASTA_M),

        b'S' | b's' => Some(FASTA_S),
        b'W' | b'w' => Some(FASTA_W),
        b'B' | b'b' => Some(FASTA_B),
        b'D' | b'd' => Some(FASTA_D),

        b'H' | b'h' => Some(FASTA_H),
        b'V' | b'v' => Some(FASTA_V),

        b'N' | b'n' => Some(FASTA_N),
        _ => None,
    }
}

// see this RFC about why we have to do the match guard approach: https://rust-lang.github.io/rfcs/1445-restrict-constants-in-patterns.html
pub fn complement(base: &FastaBase) -> FastaBase {
    match base {
        &x if x == FASTA_N => FASTA_N,
        &x if x == FASTA_A => FASTA_T,
        &x if x == FASTA_T => FASTA_A,
        &x if x == FASTA_G => FASTA_C,
        &x if x == FASTA_R => FASTA_Y,
        &x if x == FASTA_Y => FASTA_R,
        &x if x == FASTA_K => FASTA_M,
        &x if x == FASTA_M => FASTA_K,
        &x if x == FASTA_S => FASTA_W,
        &x if x == FASTA_W => FASTA_S,
        &x if x == FASTA_B => FASTA_V,
        &x if x == FASTA_V => FASTA_B,
        &x if x == FASTA_D => FASTA_H,
        &x if x == FASTA_H => FASTA_D,
        _ => panic!("Unknown base {}",base),
    }
}



/// We have a private type of FastaBase packed into a u64, right aligned
#[derive(Shrinkwrap)]
pub struct Fasta64(u64);

#[derive(Clone, From, Debug, PartialEq, Eq)]
pub struct FastaString {
    packed_bases: Vec<FastaBase>,
}

impl<Idx> std::ops::Index<Idx> for FastaString
    where
        Idx: std::slice::SliceIndex<[FastaBase]>,
{
    type Output = Idx::Output;

    fn index(&self, index: Idx) -> &Self::Output {
        &self.packed_bases[index]
    }
}

impl FastaString {
    pub fn from(string: &str) -> FastaString {
        FastaString{packed_bases: string.chars().map(|s| {
            match char_to_encoding(&s) {
                None => {panic!("unable to convert char {} to a FastaBase from string {}",s,string)}
                Some(t) => {t}
            }
        }).collect()}
    }

    #[allow(dead_code)]
    pub fn reverse_complement(&self) -> FastaString {
        FastaString{packed_bases: self.packed_bases.iter().map(|b|complement(b)).rev().collect()}
    }

    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.packed_bases.len()
    }
}


/*
pub fn hamming_bit_strings(a: &BitEncodedFasta, b: &BitEncodedFasta, width: usize) -> u64 {
    assert_eq!(a.string.len(), b.string.len());
    assert_eq!(a.length, b.length);
    assert!(width <= CHAR_PER_ENCODING);

    a.string.iter().zip(b.string.iter()).map(|(base_a, base_b)|
        (0..width).map(|i| {
            let mask = SINGLE_BIT_MASK << (((CHAR_PER_ENCODING - i) -1) * BIT_WIDTH);
            let and = base_a & base_b;
            let result_and = and & mask;
            if result_and > 0 {0} else {1}
        }).sum::<u64>()
    ).sum()
}


pub fn string_to_bit(input: &Vec<u8>) -> BitEncodedFasta {
    let mut string = Vec::new();

    for chunk in input.chunks(CHAR_PER_ENCODING) {
        let mut current_encoding: u64 = 0x0;
        for index in 0..chunk.len() {

            let converted_base = base_to_encoding(&chunk[index]);
            current_encoding = current_encoding + (converted_base << ((CHAR_PER_ENCODING - index) -1) * BIT_WIDTH);
        }
        // pad the remainder
        for index in chunk.len()..CHAR_PER_ENCODING {
            current_encoding = current_encoding + (SINGLE_BIT_MASK << ((CHAR_PER_ENCODING - index) -1) * BIT_WIDTH);
        }
        string.push(current_encoding);
    }
    BitEncodedFasta { string, length: input.len() }
}*/


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reverse_comp() {
        let fasta_string = FastaString::from("ACGT");
        let rev_comp = fasta_string.reverse_complement();
        assert_eq!(fasta_string, rev_comp);

    }

    #[test]
    fn bit_compare_simple() {
        let bit_one = FASTA_A; // StringEncodingPair { bases: FASTA_A, mask: SINGLE_BIT_MASK };
        let bit_two = FASTA_A; // StringEncodingPair { bases: FASTA_A, mask: SINGLE_BIT_MASK };
        assert_eq!(bit_one, bit_two);

        let bit_one = FASTA_A; // StringEncodingPair { bases: FASTA_A, mask: SINGLE_BIT_MASK };
        let bit_two = FASTA_C; // StringEncodingPair { bases: FASTA_C, mask: SINGLE_BIT_MASK };
        assert_ne!(bit_one, bit_two);

        let bit_one = FASTA_A; // StringEncodingPair { bases: FASTA_A, mask: SINGLE_BIT_MASK };
        let bit_two = FASTA_N; // StringEncodingPair { bases: FASTA_N, mask: SINGLE_BIT_MASK };
        assert_eq!(bit_one, bit_two);

        let bit_one = FASTA_A; // StringEncodingPair { bases: FASTA_A, mask: SINGLE_BIT_MASK };
        let bit_two = FASTA_R; // StringEncodingPair { bases: FASTA_R, mask: SINGLE_BIT_MASK };
        assert_eq!(bit_one, bit_two);
    }

    #[test]
    fn bit_compare_complex() {
        let bit_one = FASTA_R; // StringEncodingPair { bases: FASTA_A, mask: SINGLE_BIT_MASK };
        let bit_two = FASTA_A; // StringEncodingPair { bases: FASTA_A, mask: SINGLE_BIT_MASK };
        assert_eq!(bit_one, bit_two);

        let bit_one = FASTA_R; // StringEncodingPair { bases: FASTA_A, mask: SINGLE_BIT_MASK };
        let bit_two = FASTA_M; // StringEncodingPair { bases: FASTA_C, mask: SINGLE_BIT_MASK };
        assert_eq!(bit_one, bit_two);

        let bit_one = FASTA_M; // StringEncodingPair { bases: FASTA_A, mask: SINGLE_BIT_MASK };
        let bit_two = FASTA_N; // StringEncodingPair { bases: FASTA_N, mask: SINGLE_BIT_MASK };
        assert_eq!(bit_one, bit_two);

        let bit_one = FASTA_V; // StringEncodingPair { bases: FASTA_A, mask: SINGLE_BIT_MASK };
        let bit_two = FASTA_T; // StringEncodingPair { bases: FASTA_R, mask: SINGLE_BIT_MASK };
        assert_ne!(bit_one, bit_two);
    }


/*
    #[test]
    fn bit_compare_common_types() {
        let test_read = String::from("AAAAA").as_bytes().to_owned();

        let aligned_string1 = string_to_bit(&test_read);
        let aligned_string2 = string_to_bit(&test_read);

        assert_eq!(hamming_bit_strings(&aligned_string1, &aligned_string2, CHAR_PER_ENCODING), 0);

        let test_read2 = String::from("CAAAA").as_bytes().to_owned();
        let aligned_string2 = string_to_bit(&test_read2);
        assert_eq!(hamming_bit_strings(&aligned_string1, &aligned_string2, CHAR_PER_ENCODING), 1);

        let test_read2 = String::from("CCCCC").as_bytes().to_owned();
        let aligned_string2 = string_to_bit(&test_read2);
        assert_eq!(hamming_bit_strings(&aligned_string1, &aligned_string2, CHAR_PER_ENCODING), 5);

        // do 100K comparisons
        for i in 0..100000 {
            assert_eq!(hamming_bit_strings(&aligned_string1, &aligned_string2, CHAR_PER_ENCODING), 5);
        }
    }*/
}
