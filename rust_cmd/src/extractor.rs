use std::collections::{BTreeMap, HashMap};
use std::{hash::BuildHasherDefault};
use nohash_hasher::NoHashHasher;

use crate::fasta_comparisons::DEGENERATEBASES;
use crate::fasta_comparisons::KNOWNBASES;
use crate::alignment::fasta_bit_encoding::{u8_to_encoding};


pub const REFERENCE_CHAR: u8 = b'R';
pub const READ_CHAR: u8 = b'E';

lazy_static! {

    pub static ref SPECIAL_CHARACTERS: HashMap::<u8, bool, BuildHasherDefault<NoHashHasher<u8>>> = {
        let mut hashedvalues : HashMap::<u8, bool, BuildHasherDefault<NoHashHasher<u8>>> = HashMap::with_capacity_and_hasher(10, BuildHasherDefault::default());
        hashedvalues.insert(b'!',true);
        hashedvalues.insert(b'@',true);
        hashedvalues.insert(b'#',true);
        hashedvalues.insert(b'$',true);
        hashedvalues.insert(b'%',true);
        hashedvalues.insert(b'^',true);
        hashedvalues.insert(b'&',true);
        hashedvalues.insert(b'*',true);
        hashedvalues.insert(b'(',true);
        hashedvalues.insert(b')',true);
        hashedvalues
        };
}

pub fn error_out_on_unknown_base(base: &u8) {
    panic!("Unknown base (not a FASTA base or degenerate) seen in reference {}", base);
}

pub fn stretch_sequence_to_alignment(aligned_version: &Vec<u8>, native_version: &Vec<u8>) -> Vec<u8> {
    assert!(aligned_version.len() >= native_version.len());
    let mut native_result = Vec::new();
    let mut native_index = 0;
    let mut aligned_index = 0;
    while aligned_index < aligned_version.len() && native_index < native_version.len() {
        if aligned_version.get(aligned_index).unwrap() == &b'-' {
            aligned_index += 1;
            native_result.push(b'-');
        }
        else {
            aligned_index += 1;
            native_index += 1;
            native_result.push(aligned_version.get(aligned_index).unwrap().clone());
        }
    }
    assert!(aligned_index >= aligned_version.len() -1 && native_index >= native_version.len() -1);
    native_result
}

pub fn extract_tagged_sequences(aligned_read: &Vec<u8>, aligned_ref: &Vec<u8>) -> BTreeMap<u8, String> {
    let mut special_values: BTreeMap<u8, Vec<u8>> = BTreeMap::new();
    let mut in_extractor = false;

    for (reference_base, read_base) in std::iter::zip(aligned_ref, aligned_read) {
        match (u8_to_encoding(reference_base).is_some(), reference_base.is_ascii_uppercase() || (*reference_base == b'-' && in_extractor), in_extractor) {
            (_x, true, _z) => {
                in_extractor = true;
                special_values.entry(REFERENCE_CHAR).or_insert_with(Vec::new).push(reference_base.clone());
                special_values.entry(READ_CHAR).or_insert_with(Vec::new).push(read_base.clone());
            }
            (false, _y, false) if SPECIAL_CHARACTERS.contains_key(reference_base) => {
                special_values.entry(*reference_base).or_insert_with(Vec::new).push(read_base.clone());
            }
            (false, _y, true) if SPECIAL_CHARACTERS.contains_key(reference_base) => {
                special_values.entry(REFERENCE_CHAR).or_insert_with(Vec::new).push(reference_base.clone());
                special_values.entry(READ_CHAR).or_insert_with(Vec::new).push(read_base.clone());
                special_values.entry(*reference_base).or_insert_with(Vec::new).push(read_base.clone());
            }
            (_x, false, _z) => {
                in_extractor = false;
            }
        }
    }

    special_values.iter().map(|(key, value)| {
        (key.clone(), String::from_utf8(value.clone()).unwrap())
    }).collect()
}

// a custom scoring function for matching nucleotide bases to each other or
// zero-cost matches to other characters. This allows us to encode
pub fn custom_umi_score(a: u8, b: u8) -> i32 {
    match (a, b) {
        (a, b) if KNOWNBASES.contains_key(&a) && KNOWNBASES.contains_key(&b) && KNOWNBASES[&a] == KNOWNBASES[&b] => { 10 }
        (a, b) if KNOWNBASES.contains_key(&a) && KNOWNBASES.contains_key(&b) && DEGENERATEBASES.contains_key(&a) && DEGENERATEBASES[&a].contains_key(&b) => { 10 }
        (a, b) if KNOWNBASES.contains_key(&a) && KNOWNBASES.contains_key(&b) && DEGENERATEBASES.contains_key(&b) && DEGENERATEBASES[&b].contains_key(&a) => { 10 }
        (a, b) if KNOWNBASES.contains_key(&a) && KNOWNBASES.contains_key(&b) => { -8 }
        _ => { 7 } // special characters here
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tagged_sequence_test() {
        let reference = String::from("AATGATACGGCGACCACCGAGATCTACAC##########ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNNN##########CTGTAGGTAGTTTGTC").as_bytes().to_owned();
        let test_read = String::from("AATGATACGGCGACCAGATCTACACACCCCTTTGCACACTCTTTCCCTACACGACGCTCTTCCGATCTAAAAAAAATTTTTTTTTTCTGTAGGTAGTTTGTC").as_bytes().to_owned();

        let keyvalues = extract_tagged_sequences(&test_read, &reference);
        for (key, value) in keyvalues {
            println!("{} / {}", key, value);
        }
    }

    #[test]
    fn tagged_sequence_test_space() {
        let reference = String::from("AAATACTTGTACTTCGTTCAGTTACGTATTGCTAAGCAGTGGTAT*********GAGTACC------TTA--CAGTTCGATCTA").as_bytes().to_owned();
        let test_read = String::from("-------------------------------CT-AGCAG----ATCACCGTAAGGACTACCAGACGTTTAGCC-----------").as_bytes().to_owned();

        let keyvalues = extract_tagged_sequences(&test_read, &reference);

        assert_eq!(keyvalues.get(&b'*').unwrap(), "CACCGTAAG");
    }

    #[test]
    fn test_real_example() {
        let reference = String::from("tcgtcggcagcgtcagatgtgtataagagacagctagcagATCACCGTAAGGACTACCAGACGTTTAGCTGCCGGCGGAATGCTATTACTGCATTTAATGGAAGACGTTTCCGCTAAGCTCTATTTAATGTCGGGAGCCGCTTTGTAACCTGATTTACAGTCTGAGTTCATGCGAGAGAACTCTTTAATGAGTGGCCTCTCGAATCACTGAGATTTAGAGTTATCCGACACATCAAAAGGATCTTTAATGAGATGGATCGCATACTAGACAGTTGCCANNNNNNNNNNNNgcttgcactgtactctacgcgactc############agatcg").as_bytes().to_owned();
        let test_read = String::from("-----------------------------------AGCAGATCACCGTAAGGACTACCAGACGTTTAGCTGCCGGCGGAATGCTATTACTGCATTTAATGGAAGACGTTTCCGCTAAGCTCTATTTAATGTCGGGAGCCGCTTTGTAACCTGATTTACAGTCTGAGTTCATGCGAGAGAACTCTTTAATGAGTGGCCTCTCGAATCACTGAGATTTAGAGTTATCCGACA-------AGGATCTTTAATGAGATG--------------------CCACCTAGTCTCCAGGCTTGCACTGTACTCTACGCGACTCTCACCAACCGAAA----").as_bytes().to_owned();

        let keyvalues = extract_tagged_sequences(&test_read, &reference);

        println!("{:?}", &keyvalues);
        assert_eq!(keyvalues.get(&b'#').unwrap(), "TCACCAACCGAA");
    }

    #[test]
    fn lower_and_uppercase_test() {
        let reference = String::from("aaatacttgtacttcgttcaGTTACGTATTGCTAAGCAGTGGTAT*********GAGTACC------TTA--caaaaaaaaaaa").as_bytes().to_owned();
        let test_read = String::from("AAATACTTGTACTTCGTTCA-----------CT-AGCAG----ATCACCGTAAGGACTACCAGACGTTTAGCC-----------").as_bytes().to_owned();

        let keyvalues = extract_tagged_sequences(&test_read, &reference);

        assert_eq!(keyvalues.get(&READ_CHAR).unwrap(), "-----------CT-AGCAG----ATCACCGTAAGGACTACCAGACGTTTAGC");
    }
}