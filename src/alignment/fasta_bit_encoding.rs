extern crate derive_more;

use derive_more::{From, Display, Add};
use std::ops::{BitAnd, BitOr, BitXor, Shl, Shr};


/// our core Fasta base - representing a single fasta character in a u8 data store. We actually pack
/// everything into a u4 (if that was a type) but FastaString below handles the more dense packing and
/// unpacking into u64 structures
#[derive(Shrinkwrap)]
#[derive(Clone, Copy, From, Display, Add, Debug)]
#[display(fmt = "{}", _0)]
pub struct FastaBase(u8);

impl FastaBase {
    pub fn identity(&self, other: &FastaBase) -> bool {
        (*self ^ *other).0 == 0
    }
}
/// our comparisons are done using logical AND operations for speed (the layout of bits is really important).
/// Each degenerate base should also be 'equal' to the correct 'ACGT' combination and be equal to other degenerate
/// bases that share any overlap in their base patterns. E.g. R == K, but R != C (N == everything)
impl PartialEq for FastaBase {
    fn eq(&self, other: &Self) -> bool {
        self.0 & other.0 > (0 as u8)
    }
}

impl Eq for FastaBase {}

impl BitOr for FastaBase {
    type Output = FastaBase;
    fn bitor(self, rhs: FastaBase) -> FastaBase {
        FastaBase(self.0.bitor(rhs.0))
    }
}

impl BitXor for FastaBase {
    type Output = FastaBase;
    fn bitxor(self, rhs: FastaBase) -> FastaBase {
        FastaBase(self.0.bitxor(rhs.0))
    }
}

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

impl Shl<usize> for FastaBase {
    type Output = FastaBase;

    fn shl(self, rhs: usize) -> FastaBase {
        FastaBase(self.0 << rhs)
    }
}

impl Shr<usize> for FastaBase {
    type Output = FastaBase;

    fn shr(self, rhs: usize) -> FastaBase {
        FastaBase(self.0 >> rhs)
    }
}


/// TODO: so I couldn't figure out how to create the const version as a standard trait, so I did it here (self note: just look at u8 impl you dummy)
const fn add_two_string_encodings(a: FastaBase, b: FastaBase) -> FastaBase {
    FastaBase((a.0 + b.0) as u8)
}

pub const FASTA_N: FastaBase = FastaBase(0xF as u8);
pub const FASTA_A: FastaBase = FastaBase(0x1 as u8);
pub const FASTA_C: FastaBase = FastaBase(0x2 as u8);
pub const FASTA_G: FastaBase = FastaBase(0x4 as u8);
pub const FASTA_T: FastaBase = FastaBase(0x8 as u8);
pub const FASTA_R: FastaBase = add_two_string_encodings(FASTA_A, FASTA_G);
pub const FASTA_Y: FastaBase = add_two_string_encodings(FASTA_C, FASTA_T);
pub const FASTA_K: FastaBase = add_two_string_encodings(FASTA_G, FASTA_T);
pub const FASTA_M: FastaBase = add_two_string_encodings(FASTA_A, FASTA_C);
pub const FASTA_S: FastaBase = add_two_string_encodings(FASTA_C, FASTA_G);
pub const FASTA_W: FastaBase = add_two_string_encodings(FASTA_A, FASTA_T);
pub const FASTA_B: FastaBase = add_two_string_encodings(FASTA_C, add_two_string_encodings(FASTA_G, FASTA_T));
pub const FASTA_D: FastaBase = add_two_string_encodings(FASTA_A, add_two_string_encodings(FASTA_G, FASTA_T));
pub const FASTA_H: FastaBase = add_two_string_encodings(FASTA_A, add_two_string_encodings(FASTA_C, FASTA_T));
pub const FASTA_V: FastaBase = add_two_string_encodings(FASTA_A, add_two_string_encodings(FASTA_C, FASTA_G));

#[allow(dead_code)]
pub fn char_to_encoding(base: &char) -> Option<FastaBase> {
    match base {
        'A' | 'a' => Some(FASTA_A),
        'C' | 'c' => Some(FASTA_C),
        'G' | 'g' => Some(FASTA_G),
        'T' | 't' => Some(FASTA_T),

        'R' | 'r' => Some(FASTA_R),
        'Y' | 'y' => Some(FASTA_Y),
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
pub fn u8_to_encoding_defaulted_to_N(base: &u8) -> FastaBase {
    match base {
        b'A' | b'a' => FASTA_A,
        b'C' | b'c' => FASTA_C,
        b'G' | b'g' => FASTA_G,
        b'T' | b't' => FASTA_T,

        b'R' | b'r' => FASTA_R,
        b'Y' | b'y' => FASTA_Y,
        b'K' | b'k' => FASTA_K,
        b'M' | b'm' => FASTA_M,

        b'S' | b's' => FASTA_S,
        b'W' | b'w' => FASTA_W,
        b'B' | b'b' => FASTA_B,
        b'D' | b'd' => FASTA_D,

        b'H' | b'h' => FASTA_H,
        b'V' | b'v' => FASTA_V,

        _ => FASTA_N,
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
        b'Y' | b'y' => Some(FASTA_Y),
        b'K' | b'k' => Some(FASTA_K),
        b'M' | b'm' => Some(FASTA_M),

        b'S' | b's' => Some(FASTA_S),
        b'W' | b'w' => Some(FASTA_W),
        b'B' | b'b' => Some(FASTA_B),
        b'D' | b'd' => Some(FASTA_D),

        b'H' | b'h' => Some(FASTA_H),
        b'V' | b'v' => Some(FASTA_V),

        b'N' | b'n' => Some(FASTA_N),
        _ => Some(FASTA_N),
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
        _ => panic!("Unknown base {}", base),
    }
}

#[derive(Clone, From, Debug, PartialEq, Eq)]
pub struct FastaString {
    packed_bases: Vec<u64>,
    character_length: usize,
}
/*
impl<Idx> std::ops::Index<Idx> for FastaString
    where
        Idx: std::slice::SliceIndex<[FastaBase]>,
{
    type Output = Idx::Output;

    fn index(&self, index: Idx) -> &Self::Output {

        //let bin = index / FastaString::fasta_base_per_u64;
        //let offset = index % FastaString::fasta_base_per_u64;
        &self.packed_bases[index]
    }
}
*/
impl FastaString {
    const fasta_base_per_u64: usize = 64 / 4;
    const FULL_SHIFT: u64 = 0xFFFFFFFFFFFFFFFF;
    const OFFSET_SHIFT: [u64; 16] = [
        ((0xF as u64) << 60) as u64, ((0xF as u64) << 56) as u64, ((0xF as u64) << 52) as u64, ((0xF as u64) << 48) as u64,
        ((0xF as u64) << 44) as u64, ((0xF as u64) << 40) as u64, ((0xF as u64) << 36) as u64, ((0xF as u64) << 32) as u64,
        (0xF << 28) as u64, (0xF << 24) as u64, (0xF << 20) as u64, (0xF << 16) as u64,
        (0xF << 12) as u64, (0xF << 8) as u64, (0xF << 4) as u64, (0xF) as u64];


    /*
    /// a helper method to set a specific
    fn set_offset(u64_index: usize, offset: usize, value: FastaBase, full_array: &mut Vec<u64>) {

    }
    pub fn from(string: &str) -> FastaString {
        let final_bases: Vec<u64> = Vec::with_capacity((string.len() as f64 / FastaString::fasta_base_per_u64 as f64).ceil() as usize);

        string.chars().enumerate().for_each(|(index, s)| {
            let final_base_index = (index / FastaString::fasta_base_per_u64) as usize;

            match char_to_encoding(&s) {
                None => { panic!("unable to convert char {} to a FastaBase within string {}", s, string) }
                Some(t) => {
                    let offset = (FastaString::fasta_base_per_u64 - (index % FastaString::fasta_base_per_u64)) - 1;

                    final_bases[final_base_index] =
                        (u64::from(t) << (offset * 4)) as u64 | final_bases[final_base_index];
                }
            }
        });

        FastaString { packed_bases: final_bases, character_length: string.len() }
    }

    #[allow(dead_code)]
    pub fn reverse_complement(&self) -> FastaString {
        /*let final_bases: Vec<u64> = Vec::with_capacity((string.len() as f64 / FastaString::fasta_base_per_u64 as f64).ceil() as usize);
        // bases are laid out in order, we want to swap the order and switch to complement bases
        let mut forward_position = 0;
        let mut reverse_position = self.character_length - 1;
        while forward_index <= reverse_index {
            let forward_offset  = (FastaString::fasta_base_per_u64 - (forward_position % FastaString::fasta_base_per_u64)) - 1;
            let forward_index =  (forward_position / FastaString::fasta_base_per_u64) as usize;
            let reverse_offset  = (FastaString::fasta_base_per_u64 - (reverse_position % FastaString::fasta_base_per_u64)) - 1;
            let reverse_index =  (reverse_position / FastaString::fasta_base_per_u64) as usize;


        }*/

    }
*/
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
    use std::collections::HashMap;
    use super::*;

    #[test]
    fn test_identity() {
        assert_eq!(FASTA_N.identity(&FASTA_N),true);
        assert_eq!(FASTA_N.identity(&FASTA_A),false);
    }

    #[test]
    fn test_u8_to_encoding_defaulted_to_N() {
        assert_eq!(FASTA_A, u8_to_encoding_defaulted_to_N(&b'A'));
        assert_eq!(FASTA_C, u8_to_encoding_defaulted_to_N(&b'C'));
        assert_eq!(FASTA_G, u8_to_encoding_defaulted_to_N(&b'G'));
        assert_eq!(FASTA_T, u8_to_encoding_defaulted_to_N(&b'T'));

        assert_eq!(FASTA_A, u8_to_encoding_defaulted_to_N(&b'a'));
        assert_eq!(FASTA_C, u8_to_encoding_defaulted_to_N(&b'c'));
        assert_eq!(FASTA_G, u8_to_encoding_defaulted_to_N(&b'g'));
        assert_eq!(FASTA_T, u8_to_encoding_defaulted_to_N(&b't'));

        assert_eq!(FASTA_N, u8_to_encoding_defaulted_to_N(&b'N'));
        assert_eq!(FASTA_N, u8_to_encoding_defaulted_to_N(&b'n'));

        assert_eq!(FASTA_B, u8_to_encoding_defaulted_to_N(&b'B'));
        assert_eq!(FASTA_B, u8_to_encoding_defaulted_to_N(&b'b'));

        assert_eq!(FASTA_D, u8_to_encoding_defaulted_to_N(&b'D'));
        assert_eq!(FASTA_D, u8_to_encoding_defaulted_to_N(&b'd'));

        assert_eq!(FASTA_R, u8_to_encoding_defaulted_to_N(&b'R'));
        assert_eq!(FASTA_R, u8_to_encoding_defaulted_to_N(&b'r'));

        assert_eq!(FASTA_Y, u8_to_encoding_defaulted_to_N(&b'Y'));
        assert_eq!(FASTA_Y, u8_to_encoding_defaulted_to_N(&b'y'));

        assert_eq!(FASTA_K, u8_to_encoding_defaulted_to_N(&b'K'));
        assert_eq!(FASTA_K, u8_to_encoding_defaulted_to_N(&b'k'));

        assert_eq!(FASTA_M, u8_to_encoding_defaulted_to_N(&b'M'));
        assert_eq!(FASTA_M, u8_to_encoding_defaulted_to_N(&b'm'));

        assert_eq!(FASTA_S, u8_to_encoding_defaulted_to_N(&b'S'));
        assert_eq!(FASTA_S, u8_to_encoding_defaulted_to_N(&b's'));

        assert_eq!(FASTA_W, u8_to_encoding_defaulted_to_N(&b'W'));
        assert_eq!(FASTA_W, u8_to_encoding_defaulted_to_N(&b'w'));

        assert_eq!(FASTA_H, u8_to_encoding_defaulted_to_N(&b'H'));
        assert_eq!(FASTA_H, u8_to_encoding_defaulted_to_N(&b'h'));

        assert_eq!(FASTA_V, u8_to_encoding_defaulted_to_N(&b'V'));
        assert_eq!(FASTA_V, u8_to_encoding_defaulted_to_N(&b'v'));
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
    fn bit_compare_all() {
        let mut known_mapping = HashMap::new();
        let mut known_mismapping = HashMap::new();

        known_mapping.insert(b'A', vec![b'A']);
        known_mapping.insert(b'C', vec![b'C']);
        known_mapping.insert(b'G', vec![b'G']);
        known_mapping.insert(b'T', vec![b'T']);
        known_mismapping.insert(b'A', vec![b'C', b'G', b'T']);
        known_mismapping.insert(b'C', vec![b'A', b'G', b'T']);
        known_mismapping.insert(b'G', vec![b'A', b'C', b'T']);
        known_mismapping.insert(b'T', vec![b'A', b'C', b'G']);

        known_mapping.insert(b'R', vec![b'A', b'G']);
        known_mapping.insert(b'Y', vec![b'C', b'T']);
        known_mapping.insert(b'K', vec![b'G', b'T']);
        known_mapping.insert(b'M', vec![b'A', b'C']);
        known_mapping.insert(b'S', vec![b'C', b'G']);
        known_mapping.insert(b'W', vec![b'A', b'T']);

        known_mismapping.insert(b'R', vec![b'T', b'C', b'Y']);
        known_mismapping.insert(b'Y', vec![b'G', b'A', b'R']);
        known_mismapping.insert(b'K', vec![b'A', b'C', b'M']);
        known_mismapping.insert(b'M', vec![b'G', b'T', b'K']);
        known_mismapping.insert(b'S', vec![b'A', b'T', b'W']);
        known_mismapping.insert(b'W', vec![b'C', b'G', b'S']);

        known_mapping.insert(b'B', vec![b'C', b'G', b'T']);
        known_mapping.insert(b'D', vec![b'A', b'G', b'T']);
        known_mapping.insert(b'H', vec![b'A', b'C', b'T']);
        known_mapping.insert(b'V', vec![b'A', b'C', b'G']);

        known_mismapping.insert(b'B', vec![b'A']);
        known_mismapping.insert(b'D', vec![b'C']);
        known_mismapping.insert(b'H', vec![b'G']);
        known_mismapping.insert(b'V', vec![b'T']);

        known_mapping.insert(b'N', vec![b'A', b'C', b'G', b'T']);

        known_mapping.iter().for_each(|(x, y)| {
            y.iter().for_each(|z| {
                // u8_to_encoding
                assert_eq!(u8_to_encoding(x).unwrap(),
                           u8_to_encoding(z).unwrap(),
                           "Testing {} and {}",
                           String::from_utf8(vec![*x]).unwrap(),
                           String::from_utf8(vec![*z]).unwrap());

                assert_eq!(u8_to_encoding_defaulted_to_N(x),
                           u8_to_encoding_defaulted_to_N(z),
                           "Testing {} and {}",
                           String::from_utf8(vec![*x]).unwrap(),
                           String::from_utf8(vec![*z]).unwrap());
            })
        });

        known_mismapping.iter().for_each(|(x, y)| {
            y.iter().for_each(|z| {
                assert_ne!(u8_to_encoding(x).unwrap(),
                           u8_to_encoding(z).unwrap(),
                           "Testing {} and {}",
                           String::from_utf8(vec![*x]).unwrap(),
                           String::from_utf8(vec![*z]).unwrap());

                assert_ne!(u8_to_encoding_defaulted_to_N(x),
                           u8_to_encoding_defaulted_to_N(z),
                           "Testing {} and {}",
                           String::from_utf8(vec![*x]).unwrap(),
                           String::from_utf8(vec![*z]).unwrap());
            })
        });
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
