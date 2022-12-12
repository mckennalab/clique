use std::collections::{HashMap, BTreeMap};
use crate::fasta_comparisons::DEGENERATEBASES;
use crate::fasta_comparisons::KNOWNBASES;
use crate::fasta_comparisons::KNOWNBASESPLUSINSERT;

use bio::alignment::{Alignment, AlignmentMode, AlignmentOperation, TextSlice};
use bio::alignment::pairwise::{Aligner, MIN_SCORE, Scoring};


pub fn extract_tagged_sequences(aligned_read: &Vec<u8>, aligned_ref: &Vec<u8>) -> BTreeMap<String, String> {
    let mut special_values: BTreeMap<String, Vec<u8>> = BTreeMap::new();
    let empty = &Vec::new();

    let mut current_extractor = -1;
    let mut in_extractor = false;

    for (reference_base, read_base) in std::iter::zip(aligned_ref, aligned_read) {
        if !KNOWNBASESPLUSINSERT.contains_key(&reference_base) && !DEGENERATEBASES.contains_key(&reference_base) {
            let mut current_code = special_values.get(&format!("{}",reference_base)).unwrap_or(empty).clone();
            current_code.push(read_base.clone());

            special_values.insert(format!("{}",reference_base), current_code.clone());
        } else if reference_base.is_ascii_uppercase() {
            if !in_extractor {
                current_extractor += 1;
                in_extractor = true;
            }

            let mut current_code = match special_values.get(&(format!("r{}",current_extractor))) {
                None => {
                    let new_vec: Vec<u8> = Vec::new();
                    new_vec
                }
                Some(x) => { x.clone() }
            };

            current_code.push(read_base.clone());

            special_values.insert(format!("r{}",current_extractor), current_code.clone());

            let mut current_code = match special_values.get(&(format!("e{}",current_extractor))) {
                None => {
                    let new_vec: Vec<u8> = Vec::new();
                    new_vec
                }
                Some(x) => { x.clone() }
            };

            current_code.push(reference_base.clone());

            special_values.insert(format!("e{}",current_extractor), current_code.clone());
        } else {
            in_extractor = false;
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
        for _n in 1..100 {
            let reference = String::from("AAATACTTGTACTTCGTTCAGTTACGTATTGCTAAGCAGTGGTAT*********GAGTACC------TTA--CAGTTCGATCTA").as_bytes().to_owned();
            let test_read = String::from("                               CT-AGCAG----ATCACCGTAAGGACTACCAGACGTTTAGCC           ").as_bytes().to_owned();

            let keyvalues = extract_tagged_sequences(&test_read, &reference);
            assert_eq!(keyvalues.get("*").unwrap(), "CACCGTAAG");
        }
    }
}