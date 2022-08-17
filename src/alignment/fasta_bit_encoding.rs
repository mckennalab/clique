use std::ops::BitXor;

/// Encodings for all FASTA bases
const FASTA_N: u64 = 0xF;
const FASTA_A: u64 = 0x1;
const FASTA_C: u64 = 0x2;
const FASTA_G: u64 = 0x4;
const FASTA_T: u64 = 0x8;
const FASTA_R: u64 = FASTA_A + FASTA_G;
const FASTA_Y: u64 = FASTA_C + FASTA_T;
const FASTA_K: u64 = FASTA_G + FASTA_T;
const FASTA_M: u64 = FASTA_A + FASTA_C;
const FASTA_S: u64 = FASTA_C + FASTA_G;
const FASTA_W: u64 = FASTA_A + FASTA_T;
const FASTA_B: u64 = FASTA_C + FASTA_G + FASTA_T;
const FASTA_D: u64 = FASTA_A + FASTA_G + FASTA_T;
const FASTA_H: u64 = FASTA_A + FASTA_C + FASTA_T;
const FASTA_V: u64 = FASTA_A + FASTA_C + FASTA_G; // r == ACG


type StringEncoding = u64;
type StringMask = u64;

const CHAR_PER_ENCODING: usize = 16;
const SINGLE_BIT_MASK: u64 = 0xF;
const FULLSET: u64 = 0xFFFFFFFFFFFFFFFF;

const BIT_WIDTH: usize = 0x4;

/*
#[derive(Debug)]
pub struct StringEncodingPair {
    bases: StringEncoding,
    mask: StringMask,
}*/

pub fn eq_bases(a:&StringEncoding,b:&StringEncoding) -> bool {
        let comp = a & b;
        comp > 0
    }


pub struct BitEncodedFasta {
    string: Vec<StringEncoding>,
    length: usize, // with bit packing sometimes
}

pub fn base_to_encoding(base: &u8) -> StringEncoding {
    match base {
        b'A' | b'a' => FASTA_A,
        b'C' | b'c' => FASTA_C,
        b'G' | b'g' => FASTA_G,
        b'T' | b't' => FASTA_T,

        b'R' | b'r' => FASTA_R,
        b'Y' | b'y' => FASTA_C,
        b'K' | b'k' => FASTA_K,
        b'M' | b'm' => FASTA_M,

        b'S' | b's' => FASTA_S,
        b'W' | b'w' => FASTA_W,
        b'B' | b'b' => FASTA_B,
        b'D' | b'd' => FASTA_D,

        b'H' | b'h' => FASTA_H,
        b'V' | b'v' => FASTA_V,

        _ => FASTA_N, // I'm not sure if this should be the default, or we should guard against invalid characters...
    }
}

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
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bit_compare_simple() {
        let bit_one = FASTA_A; // StringEncodingPair { bases: FASTA_A, mask: SINGLE_BIT_MASK };
        let bit_two = FASTA_A; // StringEncodingPair { bases: FASTA_A, mask: SINGLE_BIT_MASK };
        assert_eq!(eq_bases(&bit_one, &bit_two), true);

        let bit_one = FASTA_A; // StringEncodingPair { bases: FASTA_A, mask: SINGLE_BIT_MASK };
        let bit_two = FASTA_C; // StringEncodingPair { bases: FASTA_C, mask: SINGLE_BIT_MASK };
        assert_eq!(eq_bases(&bit_one, &bit_two),false);

        let bit_one = FASTA_A; // StringEncodingPair { bases: FASTA_A, mask: SINGLE_BIT_MASK };
        let bit_two = FASTA_N; // StringEncodingPair { bases: FASTA_N, mask: SINGLE_BIT_MASK };
        assert_eq!(eq_bases(&bit_one, &bit_two),true);

        let bit_one = FASTA_A; // StringEncodingPair { bases: FASTA_A, mask: SINGLE_BIT_MASK };
        let bit_two = FASTA_R; // StringEncodingPair { bases: FASTA_R, mask: SINGLE_BIT_MASK };
        assert_eq!(eq_bases(&bit_one, &bit_two),true);
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
    }
}

/*

 */