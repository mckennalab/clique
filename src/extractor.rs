use std::collections::{HashMap, BTreeMap};
use crate::fasta_comparisons::DEGENERATEBASES;
use crate::fasta_comparisons::KNOWNBASES;
use crate::fasta_comparisons::KNOWNBASESPLUSINSERT;

use bio::alignment::{Alignment, AlignmentMode, AlignmentOperation, TextSlice};
use bio::alignment::pairwise::{Aligner, MIN_SCORE, Scoring};


pub fn extract_tagged_sequences(aligned_read: &Vec<u8>, aligned_ref: &Vec<u8>) -> BTreeMap<String, String> {
    let mut special_values: BTreeMap<u8, Vec<u8>> = BTreeMap::new();
    let empty = &Vec::new(); // ugh this is dumb
    for (reference_base, read_base) in std::iter::zip(aligned_ref, aligned_read) {
        if !KNOWNBASESPLUSINSERT.contains_key(&reference_base) && !DEGENERATEBASES.contains_key(&reference_base) {
            let mut current_code = special_values.get(&reference_base).unwrap_or(empty).clone();
            current_code.push(read_base.clone());

            special_values.insert(*reference_base, current_code.clone());
        } else if reference_base.is_ascii_uppercase() {
            let mut current_code = special_values.get(&('r' as u8)).unwrap_or(empty).clone();
            current_code.push(read_base.clone());

            special_values.insert('r' as u8, current_code.clone());

            let mut current_code = special_values.get(&('e' as u8)).unwrap_or(empty).clone();
            current_code.push(reference_base.clone());

            special_values.insert('e' as u8, current_code.clone());
        }
    }
    special_values.iter().map(|(key, value)| {
        (String::from_utf8(vec![*key]).unwrap(), String::from_utf8(value.clone()).unwrap())
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
            assert_eq!(keyvalues.get("*").unwrap(),"CACCGTAAG");

        }
    }

}