extern crate derive_more;

use std::cmp::Ordering;
use std::fmt;
use derive_more::{From, Add};
use std::ops::{BitAnd, BitOr, BitXor, Shl, Shr};
use serde::{Serialize,Deserialize};


/// our core Fasta base - representing a single fasta character in a u8 data store. We actually pack
/// everything into a u4 (if that was a type) but FastaString below handles the more dense packing and
/// unpacking into u64 structures
#[derive(Clone, Copy, From, Add,Serialize, Deserialize)]
pub struct FastaBase(u8);

impl FastaBase {
    pub fn identity(&self, other: &FastaBase) -> bool {
        (*self ^ *other).0 == 0
    }
}


impl fmt::Debug for FastaBase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", encoding_to_u8(self) as char)
    }
}


/// our comparisons are done using logical AND operations for speed (the layout of bits is really important).
/// Each degenerate base should also be 'equal' to the correct 'ACGT' matches and be equal to other degenerate
/// bases that share any overlap in their base patterns. E.g. R == K, but R != C (and N == everything)
impl PartialEq for FastaBase {
    fn eq(&self, other: &Self) -> bool {
        self.0 & other.0 > (0 as u8)
    }
}

impl Eq for FastaBase {}

impl PartialOrd for FastaBase {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Right now we simply offer a stable sort for bases -- there's no natural ordering to nucleotides;
/// you could argue alphabetical, but we simply sort on their underlying bit encoding
impl Ord for FastaBase {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.cmp(&other.0)
    }
}

impl fmt::Display for FastaBase {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", encoding_to_u8(&self))
    }
}


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

pub const FASTA_UNSET: FastaBase = FastaBase(0x10 as u8);
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
        '-' => Some(FASTA_UNSET),

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
pub fn encoding_to_u8(base: &FastaBase) -> u8 {
    match base {
        &x if x.identity(&FASTA_UNSET) => {b'-'},
        &x if x.identity(&FASTA_A) => {b'A'},
        &x if x.identity(&FASTA_C) => {b'C'},
        &x if x.identity(&FASTA_G) => {b'G'},
        &x if x.identity(&FASTA_T) => {b'T'},

        &x if x.identity(&FASTA_R) => {b'R'},
        &x if x.identity(&FASTA_Y) => {b'Y'},
        &x if x.identity(&FASTA_K) => {b'K'},
        &x if x.identity(&FASTA_M) => {b'M'},

        &x if x.identity(&FASTA_S) => {b'S'},
        &x if x.identity(&FASTA_W) => {b'W'},
        &x if x.identity(&FASTA_B) => {b'B'},
        &x if x.identity(&FASTA_D) => {b'D'},

        &x if x.identity(&FASTA_H) => {b'H'},
        &x if x.identity(&FASTA_V) => {b'V'},
        &x if x.identity(&FASTA_N) => {b'N'},
        _ => {panic!("Unable to convert {:?}",base)},
    }
}

//pub fn compare_fasta_vec(vec1: &[FastaBase], vec2: &[FastaBase]) ->

#[allow(dead_code)]
pub fn u8_to_encoding_defaulted_to_n(base: &u8) -> FastaBase {
    match base {
        b'-' => FASTA_UNSET,

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
        b'-' => Some(FASTA_UNSET),

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
        _ => None,
    }
}

// see this RFC about why we have to do the match guard approach: https://rust-lang.github.io/rfcs/1445-restrict-constants-in-patterns.html
pub fn complement(base: &FastaBase) -> FastaBase {
    match base {
        &x if x.identity(&FASTA_UNSET) => FASTA_UNSET,
        &x if x.identity(&FASTA_N) => FASTA_N,
        &x if x.identity(&FASTA_A) => FASTA_T,
        &x if x.identity(&FASTA_C) => FASTA_G,
        &x if x.identity(&FASTA_T) => FASTA_A,
        &x if x.identity(&FASTA_G) => FASTA_C,
        &x if x.identity(&FASTA_R) => FASTA_Y,
        &x if x.identity(&FASTA_Y) => FASTA_R,
        &x if x.identity(&FASTA_K) => FASTA_M,
        &x if x.identity(&FASTA_M) => FASTA_K,
        &x if x.identity(&FASTA_S) => FASTA_W,
        &x if x.identity(&FASTA_W) => FASTA_S,
        &x if x.identity(&FASTA_B) => FASTA_V,
        &x if x.identity(&FASTA_V) => FASTA_B,
        &x if x.identity(&FASTA_D) => FASTA_H,
        &x if x.identity(&FASTA_H) => FASTA_D,
        _ => panic!("Unknown base {}", base),
    }
}

pub fn str_to_fasta_vec(st: &str) -> Vec<FastaBase> {
    let mut ret_vec = Vec::with_capacity(st.len());
    for b in st.as_bytes() {
        ret_vec.push(u8_to_encoding_defaulted_to_n(b));
    }
    ret_vec
}

#[allow(dead_code)]
pub fn fasta_vec_to_string(vc: &Vec<FastaBase>) -> String {
    String::from_utf8(fasta_vec_to_vec_u8(vc)).unwrap()
}

pub fn fasta_vec_to_vec_u8(vc: &Vec<FastaBase>) -> Vec<u8> {
    let mut ret_vec = Vec::with_capacity(vc.len());
    for b in vc {
        ret_vec.push(encoding_to_u8(b));
    }
    ret_vec
}


pub(crate) fn reverse_complement(bases: &Vec<FastaBase>) -> Vec<FastaBase> {
    let mut new_bases = bases.clone().iter().map(|b| {
        complement(b)
    }).collect::<Vec<FastaBase>>();
    new_bases.reverse();
    new_bases
}

#[derive(Clone, From, Debug, PartialEq, Eq)]
pub struct FastaString {
    packed_bases: Vec<u64>,
    character_length: usize,
}
/*
impl FastaString {
    const fasta_base_per_u64: usize = 64 / 4;
    const FULL_SHIFT: u64 = 0xFFFFFFFFFFFFFFFF;
    const OFFSET_POSITIONS: [u64; 16] = [60,56,52,48,44,40,36,32,28,24,20,16,12,8,4,0];
    const OFFSET_SHIFT: [u64; 16] = [
        ((0xF as u64) << 60) as u64, ((0xF as u64) << 56) as u64, ((0xF as u64) << 52) as u64, ((0xF as u64) << 48) as u64,
        ((0xF as u64) << 44) as u64, ((0xF as u64) << 40) as u64, ((0xF as u64) << 36) as u64, ((0xF as u64) << 32) as u64,
        (0xF << 28) as u64, (0xF << 24) as u64, (0xF << 20) as u64, (0xF << 16) as u64,
        (0xF << 12) as u64, (0xF << 8) as u64, (0xF << 4) as u64, (0xF) as u64];

    pub fn index(&self, index: usize) -> FastaBase {
        let u64_index = (index / FastaString::fasta_base_per_u64) as usize;
        let u64_internal_offset = (index % FastaString::fasta_base_per_u64) as usize;

        u8_to_encoding_defaulted_to_N(&u8::try_from((((*(&self.packed_bases.get(u64_index).unwrap()) & (0xF << &u64_internal_offset) as u64) >> &u64_internal_offset))).unwrap())
    }

    pub fn from(string: &str) -> FastaString {
        let mut final_bases: Vec<u64> = vec![0; (string.len() as f64 / FastaString::fasta_base_per_u64 as f64).ceil() as usize];

        string.chars().enumerate().for_each(|(index, s)| {
            let final_base_index = (index / FastaString::fasta_base_per_u64) as usize;

            match char_to_encoding(&s) {
                None => { panic!("unable to convert char {} to a FastaBase within string {}", s, string) }
                Some(t) => {
                    let offset = (FastaString::fasta_base_per_u64 - (index % FastaString::fasta_base_per_u64)) - 1;

                    final_bases[final_base_index] =
                        (u64::from(t.0) << (offset * 4)) as u64 | final_bases[final_base_index];
                }
            }
        });

        FastaString { packed_bases: final_bases, character_length: string.len() }
    }

    #[allow(dead_code)]
    pub fn reverse_complement(&self) { //-> FastaString {
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

    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.packed_bases.len()
    }
}
*/


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
        assert_eq!(FASTA_A, u8_to_encoding_defaulted_to_n(&b'A'));
        assert_eq!(FASTA_C, u8_to_encoding_defaulted_to_n(&b'C'));
        assert_eq!(FASTA_G, u8_to_encoding_defaulted_to_n(&b'G'));
        assert_eq!(FASTA_T, u8_to_encoding_defaulted_to_n(&b'T'));

        assert_eq!(FASTA_A, u8_to_encoding_defaulted_to_n(&b'a'));
        assert_eq!(FASTA_C, u8_to_encoding_defaulted_to_n(&b'c'));
        assert_eq!(FASTA_G, u8_to_encoding_defaulted_to_n(&b'g'));
        assert_eq!(FASTA_T, u8_to_encoding_defaulted_to_n(&b't'));

        assert_eq!(FASTA_N, u8_to_encoding_defaulted_to_n(&b'N'));
        assert_eq!(FASTA_N, u8_to_encoding_defaulted_to_n(&b'n'));

        assert_eq!(FASTA_B, u8_to_encoding_defaulted_to_n(&b'B'));
        assert_eq!(FASTA_B, u8_to_encoding_defaulted_to_n(&b'b'));

        assert_eq!(FASTA_D, u8_to_encoding_defaulted_to_n(&b'D'));
        assert_eq!(FASTA_D, u8_to_encoding_defaulted_to_n(&b'd'));

        assert_eq!(FASTA_R, u8_to_encoding_defaulted_to_n(&b'R'));
        assert_eq!(FASTA_R, u8_to_encoding_defaulted_to_n(&b'r'));

        assert_eq!(FASTA_Y, u8_to_encoding_defaulted_to_n(&b'Y'));
        assert_eq!(FASTA_Y, u8_to_encoding_defaulted_to_n(&b'y'));

        assert_eq!(FASTA_K, u8_to_encoding_defaulted_to_n(&b'K'));
        assert_eq!(FASTA_K, u8_to_encoding_defaulted_to_n(&b'k'));

        assert_eq!(FASTA_M, u8_to_encoding_defaulted_to_n(&b'M'));
        assert_eq!(FASTA_M, u8_to_encoding_defaulted_to_n(&b'm'));

        assert_eq!(FASTA_S, u8_to_encoding_defaulted_to_n(&b'S'));
        assert_eq!(FASTA_S, u8_to_encoding_defaulted_to_n(&b's'));

        assert_eq!(FASTA_W, u8_to_encoding_defaulted_to_n(&b'W'));
        assert_eq!(FASTA_W, u8_to_encoding_defaulted_to_n(&b'w'));

        assert_eq!(FASTA_H, u8_to_encoding_defaulted_to_n(&b'H'));
        assert_eq!(FASTA_H, u8_to_encoding_defaulted_to_n(&b'h'));

        assert_eq!(FASTA_V, u8_to_encoding_defaulted_to_n(&b'V'));
        assert_eq!(FASTA_V, u8_to_encoding_defaulted_to_n(&b'v'));
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
    fn bit_complement() {
        let reference = str_to_fasta_vec("CCAATCTACTACTGCTTGCA");
        let ref_rev = reverse_complement(&reference);
        let ref_rev2 = reverse_complement(&ref_rev);
        assert_eq!(reference, ref_rev2);
    }

    #[test]
    fn bit_ordering() {
        assert!(FASTA_N > FASTA_A);
        assert!(FASTA_A < FASTA_N);
    }

    fn bit_ordering_vec() {
        let str1 = vec![FASTA_N,FASTA_N,FASTA_N,FASTA_N,FASTA_N];
        let str2 = vec![FASTA_A,FASTA_C,FASTA_G,FASTA_N,FASTA_N];
        assert!(str1 > str2);
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

                assert_eq!(u8_to_encoding_defaulted_to_n(x),
                           u8_to_encoding_defaulted_to_n(z),
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

                assert_ne!(u8_to_encoding_defaulted_to_n(x),
                           u8_to_encoding_defaulted_to_n(z),
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
    fn basic_bitstring_conversion() {
        let base = "A";
        let converted = FastaString::from(&base);
        let str_version = format!("{:#64b}",converted.packed_bases[0]);

        assert_eq!(str_version," 0b1000000000000000000000000000000000000000000000000000000000000");

        let base = "AN";
        let converted = FastaString::from(&base);
        let str_version = format!("{:#64b}",converted.packed_bases[0]);

        assert_eq!(str_version," 0b1111100000000000000000000000000000000000000000000000000000000");

        let base = "NNNNNNNNNNNNNNNN";
        let converted = FastaString::from(&base);
        let str_version = format!("{:#64b}",converted.packed_bases[0]);

        assert_eq!(str_version,"0b1111111111111111111111111111111111111111111111111111111111111111");

        let base = "NNNNNNNNNNNNNNNNN";
        let converted = FastaString::from(&base);
        let str_version = format!("{:#64b}",converted.packed_bases[0]);

        assert_eq!(str_version,"0b1111111111111111111111111111111111111111111111111111111111111111");
        let str_version = format!("{:#64b}",converted.packed_bases[1]);

        assert_eq!(str_version,"0b1111000000000000000000000000000000000000000000000000000000000000");

        let base = "NNNNNNNNNNNNNNNNA";
        let converted = FastaString::from(&base);
        let str_version = format!("{:#64b}",converted.packed_bases[0]);

        assert_eq!(str_version,"0b1111111111111111111111111111111111111111111111111111111111111111");
        let str_version = format!("{:#64b}",converted.packed_bases[1]);

        assert_eq!(str_version," 0b1000000000000000000000000000000000000000000000000000000000000");
    }


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
