use std::cmp::{min};
use nohash_hasher::NoHashHasher;
use noodles_sam::alignment::record::cigar::op::*;
use std::collections::{BTreeMap, HashMap, VecDeque};
use std::hash::BuildHasherDefault;
use itertools::Itertools;
use FASTA_UNSET;
use utils::base_utils::is_valid_fasta_base;
use utils::read_utils::u8s;
use crate::alignment::scoring_functions::AffineScoring;
use crate::alignment_functions::align_two_strings;

use crate::fasta_comparisons::DEGENERATEBASES;
use crate::fasta_comparisons::KNOWNBASES;
use crate::read_strategies::sequence_layout::{ReferenceRecord, SequenceLayout, UMIConfiguration};

pub const REFERENCE_CHAR: u8 = b'R';
pub const READ_CHAR: u8 = b'E';

lazy_static! {
    pub static ref SPECIAL_CHARACTERS: HashMap::<u8, bool, BuildHasherDefault<NoHashHasher<u8>>> = {
        let mut hashedvalues: HashMap<u8, bool, BuildHasherDefault<NoHashHasher<u8>>> =
            HashMap::with_capacity_and_hasher(10, BuildHasherDefault::default());
        hashedvalues.insert(b'0', true);
        hashedvalues.insert(b'1', true);
        hashedvalues.insert(b'2', true);
        hashedvalues.insert(b'3', true);
        hashedvalues.insert(b'4', true);
        hashedvalues.insert(b'5', true);
        hashedvalues.insert(b'6', true);
        hashedvalues.insert(b'7', true);
        hashedvalues.insert(b'8', true);
        hashedvalues.insert(b'9', true);
        hashedvalues
    };
}

pub fn error_out_on_unknown_base(base: &u8) {
    panic!(
        "Unknown base (not a FASTA base or degenerate) seen in reference {}",
        base
    );
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RecoveredAlignedSeq {
    pub aligned_read: Vec<u8>,
    pub aligned_ref: Vec<u8>,
}

pub enum SoftClipResolution {
    Clip,
    MatchMismatch,
    Realign,
}

pub fn recover_soft_clipped_align_sequences(
    unaligned_read: &[u8],
    one_based_start_pos: usize,
    cigar: &Vec<Op>,
    soft_clip_as_match: &SoftClipResolution,
    reference: &[u8],
) -> RecoveredAlignedSeq {
    let mut aligned_read = Vec::new();
    let mut aligned_ref = Vec::new();
    let mut read_pos = 0;
    let mut ref_pos = one_based_start_pos - 1;
    //println!("START read pos {} ref pos: {} ", read_pos, ref_pos);

    if ref_pos > 0 && cigar.len() > 0 && cigar[0].kind() != Kind::SoftClip {
        aligned_read.extend(b"-".repeat(ref_pos));
        aligned_ref.extend(reference[0..ref_pos].to_vec());
        //println!("ADDING");
        //ref_pos += ref_pos;
    };
    //println!("START read pos {} ref pos: {} ", read_pos, ref_pos);


    //println!("\n\nread length {} cigar {:?} ref length {:?}", unaligned_read.len(), cigar, reference.len());
    //println!("read {} ref: {}", String::from_utf8(unaligned_read.clone().to_vec()).unwrap(), String::from_utf8(reference.clone().to_vec()).unwrap(), );
    for (cigar_index, cigar_sub) in cigar.iter().enumerate() {
        let cigar_v = cigar_sub;
        let len = cigar_v.len();
        //println!("START read pos {} ref pos: {} kind {:?} length {}", read_pos, ref_pos, cigar_v.kind(), cigar_v.len());

        match cigar_v.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                // println!("Ref length {} ({} {} {:?})",(ref_pos + len) - ref_pos,ref_pos,ref_pos+len,cigar_sub);
                aligned_read.extend(unaligned_read[read_pos..read_pos + len].to_vec());
                aligned_ref.extend(reference[ref_pos..ref_pos + len].to_vec());
                read_pos += len;
                ref_pos += len;
            }
            Kind::Insertion => {
                aligned_read.extend(unaligned_read[read_pos..read_pos + len].to_vec());
                aligned_ref.extend(b"-".repeat(len));
                read_pos += len;
            }
            Kind::Deletion | Kind::Skip => {
                aligned_read.extend(b"-".repeat(len));
                aligned_ref.extend(reference[ref_pos..ref_pos + len].to_vec());
                ref_pos += len;
            }
            Kind::SoftClip => {
                match soft_clip_as_match {
                    SoftClipResolution::Clip => {
                        // just add dashes to the
                        aligned_ref.extend(unaligned_read[ref_pos..ref_pos + len].to_vec());
                        aligned_read.extend(b"-".repeat(len));
                        read_pos += len;
                        ref_pos += len;
                    }
                    SoftClipResolution::MatchMismatch => {
                        if cigar_index == 0 {
                            if ref_pos >= len {
                                let dashes = ref_pos - len;
                                aligned_ref.extend(reference[0..ref_pos].to_vec());
                                aligned_read.extend(b"-".repeat(dashes));
                                aligned_read.extend(unaligned_read[0..len].to_vec());
                            } else {
                                let ref_dashes = len - ref_pos;
                                aligned_ref.extend(b"-".repeat(ref_dashes));
                                aligned_ref.extend(reference[0..ref_pos].to_vec());
                                aligned_read.extend(unaligned_read[0..len].to_vec());
                            }

                            read_pos += len;
                        } else if ref_pos + len >= reference.len() {
                            let dashes = ref_pos + len - reference.len();
                            aligned_ref.extend(reference[ref_pos..].to_vec());
                            aligned_ref.extend(b"-".repeat(dashes));
                            aligned_read.extend(unaligned_read[read_pos..read_pos + len].to_vec());
                            read_pos += len;
                        } else {
                            aligned_read.extend(unaligned_read[read_pos..read_pos + len].to_vec());
                            aligned_ref.extend(reference[ref_pos..ref_pos + len].to_vec());
                            read_pos += len;
                        }
                    }
                    SoftClipResolution::Realign => {
                        if cigar_index == 0 {
                            // grab the reference up to this point plus the trimmed read sequence and realign
                            let clipped_read = &unaligned_read[0..len];
                            let clipped_ref = &reference[0..ref_pos];
                            println!("Ref length {}",ref_pos);
                            let aligned_sequence = align_two_strings(&"read".to_string(), clipped_ref, clipped_read, &AffineScoring::default_dna(), false, &"ref".to_string(), None);
                            aligned_ref.extend(aligned_sequence.reference_aligned.clone());
                            aligned_read.extend(aligned_sequence.read_aligned.clone());
                            read_pos += len;
                            println!("Go1 {}",u8s(&aligned_sequence.read_aligned));
                            println!("Go1 {}",u8s(&aligned_sequence.reference_aligned));

                        } else if cigar_index == cigar.len() - 1 {

                            // take the rest of the reference and the trimmed read segment and realign
                            let max_right = min(read_pos + len,unaligned_read.len());
                            let clipped_read = &unaligned_read[read_pos..max_right];
                            let clipped_ref = &reference[ref_pos..];
                            let aligned_sequence = align_two_strings(&"read".to_string(), clipped_ref, clipped_read, &AffineScoring::default_dna(), false, &"ref".to_string(), None);
                            aligned_ref.extend(aligned_sequence.reference_aligned.clone());
                            aligned_read.extend(aligned_sequence.read_aligned.clone());
                            read_pos += len;
                            ref_pos += reference[ref_pos..].len();
                            println!("Go2 {}",u8s(&aligned_sequence.read_aligned));
                            println!("Go2 {}",u8s(&aligned_sequence.reference_aligned));

                        }
                    }
                }
            }
            Kind::HardClip | Kind::Pad => {
                // do nothing
            }
        }
        //println!("END read pos {} ref pos: {} kind {:?} length {}", read_pos, ref_pos, cigar_v.kind(), cigar_v.len());
        //println!("read {} {} \nref  {} {}\n", String::from_utf8(aligned_read.clone()).unwrap(),aligned_read.len(), String::from_utf8(aligned_ref.clone()).unwrap(), aligned_ref.len());
    }
    if ref_pos < reference.len() {
        aligned_ref.extend(reference[ref_pos..].to_vec());
        aligned_read.extend(b"-".repeat(reference.len() - ref_pos));
    }

    RecoveredAlignedSeq {
        aligned_read,
        aligned_ref,
    }
}

pub fn stretch_sequence_to_alignment(aligned_version: &[u8], native_version: &[u8]) -> Vec<u8> {
    assert!(
        aligned_version.len() >= native_version.len(),
        "The aligned version {} is shorter than the native (unaligned) version {}",
        String::from_utf8(aligned_version.to_vec()).unwrap(),
        String::from_utf8(native_version.to_vec()).unwrap()
    );

    let mut native_result = Vec::new();
    let mut native_index = 0;
    let mut aligned_index = 0;
    while aligned_index < aligned_version.len() && native_index < native_version.len() {
        if aligned_version.get(aligned_index).unwrap() == &b'-' {
            aligned_index += 1;
            native_result.push(b'-');
        } else {
            native_result.push(native_version.get(native_index).unwrap().clone());
            aligned_index += 1;
            native_index += 1;
        }
    }
    //assert!(aligned_index >= aligned_version.len() -1 && native_index >= native_version.len() -1);
    native_result
}

pub fn gap_proportion_per_tag(tags: &BTreeMap<u8, String>) -> Vec<f64> {
    let mut gap_proportions = Vec::new();
    for (key, value) in tags {
        if key != &REFERENCE_CHAR && key != &READ_CHAR && key <= &b'9' && key >= &b'0' {
            let mut gap_count = 0;
            let mut total_count = 0;
            for base in value.as_bytes() {
                if base == &b'-' {
                    gap_count += 1;
                }
                total_count += 1;
            }
            gap_proportions.push(gap_count as f64 / total_count as f64);
        }
    }
    gap_proportions
}

pub fn extract_tagged_sequences(aligned_read: &[u8], aligned_ref: &[u8]) -> BTreeMap<u8, String> {
    let mut special_values: BTreeMap<u8, Vec<u8>> = BTreeMap::new();
    let mut in_extractor = false;
    let mut next_extractor_read = b'a';
    let mut next_extractor_ref = b'A';

    for (reference_base, read_base) in std::iter::zip(aligned_ref, aligned_read) {
        match (
            is_valid_fasta_base(reference_base),
            reference_base.is_ascii_uppercase() || (*reference_base == b'-' && in_extractor),
            in_extractor,
        ) {
            (_x, true, _z) => {
                in_extractor = true;
                special_values
                    .entry(next_extractor_ref)
                    .or_insert_with(Vec::new)
                    .push(reference_base.clone());
                special_values
                    .entry(next_extractor_read)
                    .or_insert_with(Vec::new)
                    .push(read_base.clone());
            }
            (false, _y, false) if SPECIAL_CHARACTERS.contains_key(reference_base) => {
                special_values
                    .entry(*reference_base)
                    .or_insert_with(Vec::new)
                    .push(read_base.clone());
            }
            (false, _y, true) if SPECIAL_CHARACTERS.contains_key(reference_base) => {
                special_values
                    .entry(next_extractor_ref)
                    .or_insert_with(Vec::new)
                    .push(reference_base.clone());
                special_values
                    .entry(next_extractor_read)
                    .or_insert_with(Vec::new)
                    .push(read_base.clone());
                special_values
                    .entry(*reference_base)
                    .or_insert_with(Vec::new)
                    .push(read_base.clone());
            }
            (_x, false, _z) => {
                if in_extractor {
                    next_extractor_read += 1;
                    next_extractor_ref += 1;
                }
                in_extractor = false;
            }
        }
    }


    special_values
        .iter()
        .map(|(key, value)| {
            //println!("key {} value {}",key.clone(), String::from_utf8(value.clone()).unwrap());
            (key.clone(), String::from_utf8(value.clone()).unwrap())
        })
        .collect()
}

pub fn get_sorting_order(read_structure: &SequenceLayout, reference_name: &String) -> Vec<char> {
    match read_structure.references.get(reference_name) {
        None => {
            panic!("Unable to find reference {} ", reference_name);
        }
        Some(reference) => {
            let mut sorting_order = reference
                .umi_configurations
                .iter()
                .map(|x| x.1.clone())
                .collect::<Vec<UMIConfiguration>>();
            sorting_order.sort_by(|a, b| a.order.cmp(&b.order));
            let sorted_tags = sorting_order
                .iter()
                .map(|x| x.symbol)
                .collect::<Vec<char>>();
            sorted_tags
        }
    }
}

pub fn extract_tag_sequences(
    reference_tags: &ReferenceRecord,
    ets: BTreeMap<u8, String>,
) -> (bool, VecDeque<(char, Vec<u8>)>) {
    let mut invalid_read = false;

    let mut vec_tags =    reference_tags
            .umi_configurations
            .iter()
            .map(|(_umi_name, umi_obj)| {
                let ets_hit = ets.get(&umi_obj.symbol.to_string().as_bytes()[0]);

                match ets_hit {
                    Some(e) => {
                        if e.len() != umi_obj.length {
                            invalid_read = true;
                        };
                        let mut gaps = 0;
                        let str = e.as_bytes()
                            .into_iter()
                            .map(|f| {
                                if *f == FASTA_UNSET {
                                    gaps += 1;
                                }
                                *f
                            })
                            .collect::<Vec<u8>>();

                        if gaps > umi_obj.max_gaps.unwrap_or(gaps) {
                            println!("tossing reads {} {}", gaps, e);
                            invalid_read = true;
                        } //else {
                        //println!("keeping reads {} {} level {}",gaps, e, umi_obj.max_gaps.unwrap_or(gaps));
                        //}

                        Some((umi_obj.order, (
                            umi_obj.symbol,
                            str
                        )))
                    }

                    None => {
                        println!("None read");

                        invalid_read = true;
                        None
                    }
                }
            })
            .filter(|x| x.is_some())
            .map(|x| x.unwrap())
            .collect::<Vec<(usize,(char, Vec<u8>))>>();
        vec_tags.sort_by(|x,y| x.0.cmp(&y.0));
        let queue = VecDeque::from(vec_tags.into_iter().map(|(x,y)| y).collect::<Vec<(char, Vec<u8>)>>());

    (invalid_read, queue)
}

// a custom scoring function for matching nucleotide bases to each other or
// zero-cost matches to other characters. This allows us to encode
pub fn custom_umi_score(a: u8, b: u8) -> i32 {
    match (a, b) {
        (a, b)
        if KNOWNBASES.contains_key(&a)
            && KNOWNBASES.contains_key(&b)
            && KNOWNBASES[&a] == KNOWNBASES[&b] =>
            {
                10
            }
        (a, b)
        if KNOWNBASES.contains_key(&a)
            && KNOWNBASES.contains_key(&b)
            && DEGENERATEBASES.contains_key(&a)
            && DEGENERATEBASES[&a].contains_key(&b) =>
            {
                10
            }
        (a, b)
        if KNOWNBASES.contains_key(&a)
            && KNOWNBASES.contains_key(&b)
            && DEGENERATEBASES.contains_key(&b)
            && DEGENERATEBASES[&b].contains_key(&a) =>
            {
                10
            }
        (a, b) if KNOWNBASES.contains_key(&a) && KNOWNBASES.contains_key(&b) => -8,
        _ => 7, // special characters here
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gap_proportion_per_tag_test() {
        let mut tags = BTreeMap::new();
        tags.insert(b'0', String::from("ACGT"));
        tags.insert(b'1', String::from("ACGT"));
        assert_eq!(
            gap_proportion_per_tag(&tags)
                .iter()
                .max_by(|a, b| a.total_cmp(b))
                .unwrap(),
            &0.0
        );

        tags.insert(b'1', String::from("AC--"));
        assert_eq!(
            gap_proportion_per_tag(&tags)
                .iter()
                .max_by(|a, b| a.total_cmp(b))
                .unwrap(),
            &0.5
        );

        tags.insert(b'1', String::from("----"));
        assert_eq!(
            gap_proportion_per_tag(&tags)
                .iter()
                .max_by(|a, b| a.total_cmp(b))
                .unwrap(),
            &1.0
        );
    }

    #[test]
    fn tagged_sequence_test() {
        let reference = String::from("AATGATACGGCGACCACCGAGATCTACAC0000000000ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNNN1111111111CTGTAGGTAGTTTGTC").as_bytes().to_owned();
        let test_read = String::from("AATGATACGGCGACCAGATCTACACACCCCTTTGCACACTCTTTCCCTACACGACGCTCTTCCGATCTAAAAAAAATTTTTTTTTTCTGTAGGTAGTTTGTC").as_bytes().to_owned();

        let keyvalues = extract_tagged_sequences(&test_read, &reference);
        for (key, value) in keyvalues {
            println!("{} / {}", key, value);
        }
    }

    #[test]
    fn tagged_sequence_test_space() {
        let reference = String::from(
            "AAATACTTGTACTTCGTTCAGTTACGTATTGCTAAGCAGTGGTAT111111111GAGTACC------TTA--CAGTTCGATCTA",
        )
            .as_bytes()
            .to_owned();
        let test_read = String::from(
            "-------------------------------CT-AGCAG----ATCACCGTAAGGACTACCAGACGTTTAGCC-----------",
        )
            .as_bytes()
            .to_owned();

        let keyvalues = extract_tagged_sequences(&test_read, &reference);

        assert_eq!(keyvalues.get(&b'1').unwrap(), "CACCGTAAG");
    }

    #[test]
    fn test_real_example() {
        let reference = String::from("tcgtcggcagcgtcagatgtgtataagagacagctagcagATCACCGTAAGGACTACCAGACGTTTAGCTGCCGGCGGAATGCTATTACTGCATTTAATGGAAGACGTTTCCGCTAAGCTCTATTTAATGTCGGGAGCCGCTTTGTAACCTGATTTACAGTCTGAGTTCATGCGAGAGAACTCTTTAATGAGTGGCCTCTCGAATCACTGAGATTTAGAGTTATCCGACACATCAAAAGGATCTTTAATGAGATGGATCGCATACTAGACAGTTGCCANNNNNNNNNNNNgcttgcactgtactctacgcgactc111111111111agatcg").as_bytes().to_owned();
        let test_read = String::from("-----------------------------------AGCAGATCACCGTAAGGACTACCAGACGTTTAGCTGCCGGCGGAATGCTATTACTGCATTTAATGGAAGACGTTTCCGCTAAGCTCTATTTAATGTCGGGAGCCGCTTTGTAACCTGATTTACAGTCTGAGTTCATGCGAGAGAACTCTTTAATGAGTGGCCTCTCGAATCACTGAGATTTAGAGTTATCCGACA-------AGGATCTTTAATGAGATG--------------------CCACCTAGTCTCCAGGCTTGCACTGTACTCTACGCGACTCTCACCAACCGAAA----").as_bytes().to_owned();

        let keyvalues = extract_tagged_sequences(&test_read, &reference);

        println!("{:?}", &keyvalues);
        assert_eq!(keyvalues.get(&b'1').unwrap(), "TCACCAACCGAA");
    }

    #[test]
    fn lower_and_uppercase_test() {
        let reference = String::from(
            "aaatacttgtacttcgttcaGTTACGTATTGCTAAGCAGTGGTAT111111111GAGTACC------TTA--caaaaaaaaaaa",
        )
            .as_bytes()
            .to_owned();
        let test_read = String::from(
            "AAATACTTGTACTTCGTTCA-----------CT-AGCAG----ATCACCGTAAGGACTACCAGACGTTTAGCC-----------",
        )
            .as_bytes()
            .to_owned();

        let keyvalues = extract_tagged_sequences(&test_read, &reference);

        keyvalues
            .iter()
            .for_each(|(key, value)| println!("key-value {} / {}", *key as char, value));

        assert_eq!(
            keyvalues.get(&b'A').unwrap(),
            "GTTACGTATTGCTAAGCAGTGGTAT111111111GAGTACC------TTA--"
        );
        assert_eq!(
            keyvalues.get(&b'a').unwrap(),
            "-----------CT-AGCAG----ATCACCGTAAGGACTACCAGACGTTTAGC"
        );
    }

    pub fn strip_gaps(str: &Vec<u8>) -> Vec<u8> {
        str.iter().filter(|x| **x != b'-').map(|x|*x).collect()
    }

    #[test]
    fn test_recover_align_sequences() {
        let read =                   "TTCCGATCTGTCATAACACCACACTAGAATCACGCGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT".as_bytes();
        let reference = "CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTAGGAAAGGACAGTGGGAGTGGCACC".as_bytes();
        let start_pos: usize = 23;
        let mut cigar = Vec::new();
        cigar.push(Op::new(Kind::SoftClip, 9));
        cigar.push(Op::new(Kind::Match, 58));


        let recovered = recover_soft_clipped_align_sequences(
            read,
            start_pos,
            &cigar,
            &SoftClipResolution::Realign,
            reference,
        );
        println!("Ref\n{}\nread\n{} ",String::from_utf8(recovered.aligned_ref.clone()).unwrap(),String::from_utf8(recovered.aligned_read.clone()).unwrap());
        assert_eq!(u8s(&"-------------TTCCGATCTGTCATAACACCACACTAGAATCACGCGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT-----------------------------".as_bytes().to_vec()),u8s(&recovered.aligned_read));
        assert_eq!(read,strip_gaps(&recovered.aligned_read));
        assert_eq!(     u8s(&"CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTAGGAAAGGACAGTGGGAGTGGCACC".as_bytes().to_vec()),u8s(&recovered.aligned_ref));


        let read =                        "TTCCGATCTGTCATAACACCACACTAGAATCACGCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTAGGAAAGGACAGTGGGAGTGGCA".as_bytes();
        let reference =      "CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTAGGAAAGGACAGTGGGAGTGGCACC".as_bytes();
        let start_pos: usize = 14;
        let mut cigar = Vec::new();
        cigar.push(Op::new(Kind::Match, 38));
        cigar.push(Op::new(Kind::Insertion, 4));
        cigar.push(Op::new(Kind::Match, 54));
        cigar.push(Op::new(Kind::SoftClip, 2));


        let recovered = recover_soft_clipped_align_sequences(
            read,
            start_pos,
            &cigar,
            &SoftClipResolution::Realign,
            reference,
        );

        println!("Ref\n{}\nread\n{} ",String::from_utf8(recovered.aligned_ref.clone()).unwrap(),String::from_utf8(recovered.aligned_read.clone()).unwrap());
        assert_eq!(u8s(&"-------------TTCCGATCTGTCATAACACCACACTAGAATCACGCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTAGGAAAGGACAGTGGGAGTGGCA--".as_bytes().to_vec()),u8s(&recovered.aligned_read));
        assert_eq!(u8s(&"CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNT----TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTAGGAAAGGACAGTGGGAGTGGCACC".as_bytes().to_vec()),u8s(&recovered.aligned_ref));
        assert_eq!(read,strip_gaps(&recovered.aligned_read));
        assert_eq!(reference,strip_gaps(&recovered.aligned_ref));

        /*
        //                      TTCCGATCTGTCATAACACCACACTAGAATCACGCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAGCGATGCAATTTCCTCATTTTATTAGGAAAGGACAGTGGGAGTGGCACCTTCCAGGGTCAAGGAAGGCACGGGGGAGGGGCAAACAACAGATGGCTGGCAACTAGAAGGCACAGTGAGCTTGTACATAACTACGCAAGTCCTTGCTAGGACCGGCCTTAAAGCCACGTGGCGGCCGCCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCGTGGTACAATATGCGTCTCCGAAATTAACCCGGTGCGTTTAAACGAAAAGGACCGACTACTACCTCGCGAAAGCTCTAAGCGTCGTGTCAGCGAAACTTCGCGGAGGTTCGACATCGAAAGACACGCGGGTGTATGTGGCGAAAGCAGCAACCTGATCTGGGGTGAAAAGCCATGGACGCCGGGACGAGAAAGGTCTAGGACTGTTTTGCGAGAAAAGGATTAGAGTTAGAATCGCGAAACGCTCGCGTTCCACCGCTCCGAAAGATCCCGAGGTCGTTTTACCGAAAGCGACGACTTCTGTCATAGTGAAACGATTGGACGTCTCTGGTGCGAAATCGCGGGTTGTACAACATACGAAACCGAGGCTACAACCCCGGACGAAAAGGTATAGGTAGCTAACACGCGAAACCCTAGGGATCGTGCTAGCCGAAAGCCCTATCACGCAGGGGACTGAAAAACATGGGCACGCCCCCGATGAAACGCTGCTTGTCTGGCCTCGCGAAAGAATGAGCAGAGCGTGAGGCGAAAAGCTTAAGCTGTGCACTCTCGAAAGTCGGTGTCCATCAGTGGATGAAACAGCGGGTTCCTGCTCCCGCGAAACGCCACCTGTACGTTACTTCGAAAATGAAGGGACAGCGGCGGACGAAAGTCATATTCCGTTGTGGTACGAAATTGGTCCTGATGCACGCACAGAAAAGATTGACCTCCGTTCGTACGAAAGCTCGGCCTCTGGGAGTCGTGAAAGACTCGGATCCGCACCAGATGAAAGGCACACCCACGCCCGTCACGAAAACCCAAACCTTGTATGTATGGAAATCTTCTGCGTCCGGGCCGCGGAAAAGCGTATACCTATCTCGCATGAAAGTCTCTCACCTCGTCTACGCGAAACGCTCGTACGCGTACGGGCTGAAAGCGATACACCGCTCGCCCCTGAAACCCTCTAGTTACGCGCCAGTGAAAGAGTCGCGTAGAGTACAGTGCAAGGTCGACAATCAACCTCTGGATTACATCCGATTGCCTCACTGTGCGAAAGTACTCGATGGCGTGGCTTAGAAAGCGTACAGTCTCCGTGCCGGGAAAATAAGAGCGCCTGCGGTTATGAAATCGTGGGCTACTCCTGGGTGGAAAGCTATCCTGCACATTAGTACGAAAGGTGCCAGGTTGCTTCGATCGAAAGCCCGAGAGATCACTCGTAGGAAACTACGCCGGTCACGACGGGCGAAACGACATGAACTCATCCGGACGAAAGGTAGTCCTTACGGTGATCTGCTAGGGTCTCTCCTAGCAACGGTTACTCCATCTGGTACACCCCCTGCTCGGGGCAAGTACCTGATGCGGCACAATGTCTAGCAGGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTTAACCTTAGCAATACGTCTGTCGAAGCAGCTACAA
        let read =      "TTTGTCCTGTACTTCGTTCAGCGTATTGCTAAGGTTAACACAAAGACACCGACAACTTTCTTCAGCACCTCTACACGACGCTCTTCCGATCTAGCACCACACATACCCCCGTGCGAGCGCTTTTTTTTTTTTTTTTTTTTTTTTTTGGCGATGCAATTTCCTCATTTTATTAGGAAAGGACAGTGGGAGTGGCACCTTCCAGGGTCAAGGAAGGCACGGGGGAGGGGCAAACAACAGATGGCTGGCAACTAGAAGGCACAGTGAGCTTGTACATAACTACGCAAGTCCTTGCTAGGACCGGCCTTAAAGCCACGTGGCGGCCGCCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCGTGGCACAACACGCGTCTCCGAAATTAACCCGGCGCGTTTAAACGAAAAGGACCGACCACCACCTCGCGAAAGCTCTAAGCGTCGTGTCAGCGAAACTTCGCGGAGGCCCGACATCGAAAGACACGCGGGCGCATATGGCGAAAGCAGCAACCTGACCTGGGGTGAAAAGCCATGGACGCCGGGACGAGAAAGGTCTAGGACCGTTTTGCGAGAAAAGGATTAGAGCTAGAATCGCGAAACGCTCGCGTCCCACCGCTCCGAAAGATCCCGAGGCCGTTTTACCGAAAGCGACGACTTCCGTCATAGTGAAACGATTGGACGCCTCTGGTGCGAAATCGCGGGTCGCACAACATACGAAAACCGAGGCTACAACCCCGGACGAAAGGTATAGGCAGCCAACACGCGAAACCCTAGGGATCGCGCTAGCCGAAAGCCCTATCACGGAGGGGGCTGAAAAACATGGGCACGCCCCCGATGAAACGCTGCTTGCCCGGCCTCGCGAAAGAATGAGCTGAGCGTGAGGCGAAAAGCTTAAGCTGCGCACTCTCGAAAGTCGGTGTCCATCAGTGGATGAAACAGCGGGTTCCTGCTCCCGCGAAACGCCACCTGCACGTTACTTCGAAAATGAAGGGACAGCGGCGGACGAAAGTCACATCCCGCCGTGGTACGAAATTGGTCCTGATGCACGCACAGAAAAGATTGACCTCTGTTCGTACGAAAGCTCGGCCTCTGGGAGTCGTGAAAGACTCGGATCCGTACCAGATGAAAGGCACACCCACGTCCGTCACGAAAACCCAAACCTTGTATGTATGGAAATCTTCTGCGTTCGGGCCGCGGAAAAGCGTATACCTATCTCGCATGAAAGTCTCTCATCCCGTCTACGCGAAACGCTCGTATGCGTACGGGCTGAAAGTTGATATACTGTTCGCCCCTGAAACCCTCTAGTCACGCGCCAGTGAAAGAGTCGCGTAGAGTACAGTGCAAGGTCGACAATCAACCTCTGGATTACATCCGATTGCCTCACTGTGCGAAAGTACTCGATGGCGTGGCTTAGAAAGCGTACAGTCTCTGTGCCGGGAAAATAAGAGCGTCTGCGGTTATGAAATCGTGGGCTACCCCTGGGTGGAAAGCTATCCCGCACATTAGTACGAAAGGTGCCAGGTTGCTTCGATCGAAAGCCCGAGAGATCACTCGTAGGAAACTACGCCGGTTACGACGGGCGAAACGACATGAACTCATCCGGACGAAAGGTAGTCCTTACGGTGATCTGCTAGGGTCTCTCCTAGCAACGGTTACTCGATTTGGTACGCAAACGCCCTGCCTTCGGTACCTGATGCGGCACAATGTCTAGCAGGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTTAACCTTAGCAATACG".as_bytes();
        let reference = "CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTAGGAAAGGACAGTGGGAGTGGCACCTTCCAGGGTCAAGGAAGGCACGGGGGAGGGGCAAACAACAGATGGCTGGCAACTAGAAGGCACAGTGAGCTTGTACATAACTACGCAAGTCCTTGCTAGGACCGGCCTTAAAGCCACGTGGCGGCCGCCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCGTGGTACAATATGCGTCTCCGAAATTAACCCGGTGTGTTTAAACGAAAAGGACCGACTACTACCTCGCGAAAGCTCTAAGTGTTGTGTCAGCGAAACTTCGCGGAGGTTCGACATCGAAAGACACGCGGGTGTATATGGCGAAAGCAGCAACCTGATCTGGGGTGAAAAGCCATGGATGTCGGGACGAGAAAGGTCTAGGACTGTTTTGCGAGAAAAGGATTAGAGTTAGAATCGCGAAACGCTCGCGTTCTACCGCTCCGAAAGATCCCGAGGTTGTTTTACCGAAAGCGACGACTTCTGTCATAGTGAAACGATTGGACGTCTCTGGTGCGAAATCGCGGGTTGTACAACATACGAAACCGAGGCTATAATCCCGGACGAAAGGTATAGGTAGCTAACACGCGAAACCCTAGGGATCGTGCTAGCCGAAAGCCCTATTATGTAGGGGACTGAAAAACATGGGTACGTCCCCGATGAAACGCTGCTTGTCTGGCCTCGCGAAAGAATGAGCTGAGTGTGAGGCGAAAAGCTTAAGCTGTGCACTCTCGAAAGTCGGTGTCTATTAGTGGATGAAACAGCGGGTTCCTGCTCCCGCGAAACGCCACCTGTATGTTACTTCGAAAATGAAGGGATAGTGGCGGACGAAAGTCATATTCCGTTGTGGTACGAAATTGGTCCTGATGTACGCACAGAAAAGATTGACCTCTGTTCGTACGAAAGCTCGGCCTCTGGGAGTCGTGAAAGACTCGGATCCGTACCAGATGAAAGGCACACCCATGTCCGTCACGAAAACCCAAACCTTGTATGTATGGAAATCTTCTGCGTTCGGGCCGCGGAAAAGCGTATACCTATCTCGCATGAAAGTCTCTTATCTTGTCTACGCGAAACGCTCGTATGCGTACGGGCTGAAAGCGATATACTGTTCGCCCCTGAAACCCTCTAGTTATGCGCCAGTGAAAGAGTCGCGTAGAGTACAGTGCAAGGTCGACAATCAACCTCTGGATTACATCCGATTGCCTTACTGTGCGAAAGTACTCGATGGTGTGGCTTAGAAAGCGTACAGTCTCTGTGCCGGGAAAATAAGAGCGTCTGCGGTTATGAAATCGTGGGCTACTCCTGGGTGGAAAGCTATCCTGTATATTAGTACGAAAGGTGCCAGGTTGCTTCGATCGAAAGCCCGAGAGATTACTCGTAGGAAACTACGCCGGTTACGACGGGCGAAACGACATGAACTTATCCGGACGAAAGGTAGTCCTTACGGTGATCTGCTAGGGTCTCTCCTAGCAACGGTTACTCGATTTGGTACNNNNNNNNNNNNNNNNNNGTACCTGATGCGGCACAATGTCTAGC".as_bytes();
        let start_pos: usize = 1;
        let mut cigar = Vec::new();
        // 42S 24M 10I 3M 5I 591M 1I 970M 68S
        // 42S 24M 10I 3M 5I 591M 1I 970M 68S
        // 42S 24M 10I 3M 5I 591M 1I 970M 68S
        // 70S 74M 12I 3M 5I 567M 1I 556M 1I 438M 49S
        cigar.push(Op::new(Kind::SoftClip, 70));
        cigar.push(Op::new(Kind::Match, 74));
        cigar.push(Op::new(Kind::Insertion, 12));
        cigar.push(Op::new(Kind::Match, 3));
        cigar.push(Op::new(Kind::Insertion, 5));
        cigar.push(Op::new(Kind::Match, 567));
        cigar.push(Op::new(Kind::Insertion, 1));
        cigar.push(Op::new(Kind::Match, 556));
        cigar.push(Op::new(Kind::Insertion, 1));
        cigar.push(Op::new(Kind::Match, 438));
        cigar.push(Op::new(Kind::SoftClip, 49));

        let recovered = recover_align_sequences(
            read,
            start_pos,
            &cigar,
            &SoftClipResolution::Realign,
            reference,
        );
        println!("*\n*\n*\n*\n*\n*\nRef\n{}\nread\n{} ",String::from_utf8(recovered.aligned_ref.clone()).unwrap(),String::from_utf8(recovered.aligned_read.clone()).unwrap());
        assert_eq!("-------------TTCCGATCTGTCATAACACCACACTAGAATCACGCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTAGGAAAGGACAGTGGGAGTGGCA--".as_bytes().to_vec(),recovered.aligned_read);
        assert_eq!("CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNTT---TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTAGGAAAGGACAGTGGGAGTGGCACC".as_bytes().to_vec(),recovered.aligned_ref);

        println!("Ref\n{}\nread\n{} ",String::from_utf8(recovered.aligned_ref.clone()).unwrap(),String::from_utf8(recovered.aligned_read.clone()).unwrap());
        assert_eq!("-------------TTCCGATCTGTCATAACACCACACTAGAATCACGCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTAGGAAAGGACAGTGGGAGTGGCA--".as_bytes().to_vec(),recovered.aligned_read);
        assert_eq!("CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNT----TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTAGGAAAGGACAGTGGGAGTGGCACC".as_bytes().to_vec(),recovered.aligned_ref);
        assert_eq!(read,strip_gaps(&recovered.aligned_read));
        assert_eq!(reference,strip_gaps(&recovered.aligned_ref));

        //                      TTCCGATCTGTCATAACACCACACTAGAATCACGCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAGCGATGCAATTTCCTCATTTTATTAGGAAAGGACAGTGGGAGTGGCACCTTCCAGGGTCAAGGAAGGCACGGGGGAGGGGCAAACAACAGATGGCTGGCAACTAGAAGGCACAGTGAGCTTGTACATAACTACGCAAGTCCTTGCTAGGACCGGCCTTAAAGCCACGTGGCGGCCGCCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCGTGGTACAATATGCGTCTCCGAAATTAACCCGGTGCGTTTAAACGAAAAGGACCGACTACTACCTCGCGAAAGCTCTAAGCGTCGTGTCAGCGAAACTTCGCGGAGGTTCGACATCGAAAGACACGCGGGTGTATGTGGCGAAAGCAGCAACCTGATCTGGGGTGAAAAGCCATGGACGCCGGGACGAGAAAGGTCTAGGACTGTTTTGCGAGAAAAGGATTAGAGTTAGAATCGCGAAACGCTCGCGTTCCACCGCTCCGAAAGATCCCGAGGTCGTTTTACCGAAAGCGACGACTTCTGTCATAGTGAAACGATTGGACGTCTCTGGTGCGAAATCGCGGGTTGTACAACATACGAAACCGAGGCTACAACCCCGGACGAAAAGGTATAGGTAGCTAACACGCGAAACCCTAGGGATCGTGCTAGCCGAAAGCCCTATCACGCAGGGGACTGAAAAACATGGGCACGCCCCCGATGAAACGCTGCTTGTCTGGCCTCGCGAAAGAATGAGCAGAGCGTGAGGCGAAAAGCTTAAGCTGTGCACTCTCGAAAGTCGGTGTCCATCAGTGGATGAAACAGCGGGTTCCTGCTCCCGCGAAACGCCACCTGTACGTTACTTCGAAAATGAAGGGACAGCGGCGGACGAAAGTCATATTCCGTTGTGGTACGAAATTGGTCCTGATGCACGCACAGAAAAGATTGACCTCCGTTCGTACGAAAGCTCGGCCTCTGGGAGTCGTGAAAGACTCGGATCCGCACCAGATGAAAGGCACACCCACGCCCGTCACGAAAACCCAAACCTTGTATGTATGGAAATCTTCTGCGTCCGGGCCGCGGAAAAGCGTATACCTATCTCGCATGAAAGTCTCTCACCTCGTCTACGCGAAACGCTCGTACGCGTACGGGCTGAAAGCGATACACCGCTCGCCCCTGAAACCCTCTAGTTACGCGCCAGTGAAAGAGTCGCGTAGAGTACAGTGCAAGGTCGACAATCAACCTCTGGATTACATCCGATTGCCTCACTGTGCGAAAGTACTCGATGGCGTGGCTTAGAAAGCGTACAGTCTCCGTGCCGGGAAAATAAGAGCGCCTGCGGTTATGAAATCGTGGGCTACTCCTGGGTGGAAAGCTATCCTGCACATTAGTACGAAAGGTGCCAGGTTGCTTCGATCGAAAGCCCGAGAGATCACTCGTAGGAAACTACGCCGGTCACGACGGGCGAAACGACATGAACTCATCCGGACGAAAGGTAGTCCTTACGGTGATCTGCTAGGGTCTCTCCTAGCAACGGTTACTCCATCTGGTACACCCCCTGCTCGGGGCAAGTACCTGATGCGGCACAATGTCTAGCAGGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTTAACCTTAGCAATACGTCTGTCGAAGCAGCTACAA
        let read =      "TTCCGATCTGTCATAACACCACACTAGAATCACGCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAGCGATGCAATTTCCTCATTTTATTAGGAAAGGACAGTGGGAGTGGCACCTTCCAGGGTCAAGGAAGGCACGGGGGAGGGGCAAACAACAGATGGCTGGCAACTAGAAGGCACAGTGAGCTTGTACATAACTACGCAAGTCCTTGCTAGGACCGGCCTTAAAGCCACGTGGCGGCCGCCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCGTGGTACAATATGCGTCTCCGAAATTAACCCGGTGCGTTTAAACGAAAAGGACCGACTACTACCTCGCGAAAGCTCTAAGCGTCGTGTCAGCGAAACTTCGCGGAGGTTCGACATCGAAAGACACGCGGGTGTATGTGGCGAAAGCAGCAACCTGATCTGGGGTGAAAAGCCATGGACGCCGGGACGAGAAAGGTCTAGGACTGTTTTGCGAGAAAAGGATTAGAGTTAGAATCGCGAAACGCTCGCGTTCCACCGCTCCGAAAGATCCCGAGGTCGTTTTACCGAAAGCGACGACTTCTGTCATAGTGAAACGATTGGACGTCTCTGGTGCGAAATCGCGGGTTGTACAACATACGAAACCGAGGCTACAACCCCGGACGAAAAGGTATAGGTAGCTAACACGCGAAACCCTAGGGATCGTGCTAGCCGAAAGCCCTATCACGCAGGGGACTGAAAAACATGGGCACGCCCCCGATGAAACGCTGCTTGTCTGGCCTCGCGAAAGAATGAGCAGAGCGTGAGGCGAAAAGCTTAAGCTGTGCACTCTCGAAAGTCGGTGTCCATCAGTGGATGAAACAGCGGGTTCCTGCTCCCGCGAAACGCCACCTGTACGTTACTTCGAAAATGAAGGGACAGCGGCGGACGAAAGTCATATTCCGTTGTGGTACGAAATTGGTCCTGATGCACGCACAGAAAAGATTGACCTCCGTTCGTACGAAAGCTCGGCCTCTGGGAGTCGTGAAAGACTCGGATCCGCACCAGATGAAAGGCACACCCACGCCCGTCACGAAAACCCAAACCTTGTATGTATGGAAATCTTCTGCGTCCGGGCCGCGGAAAAGCGTATACCTATCTCGCATGAAAGTCTCTCACCTCGTCTACGCGAAACGCTCGTACGCGTACGGGCTGAAAGCGATACACCGCTCGCCCCTGAAACCCTCTAGTTACGCGCCAGTGAAAGAGTCGCGTAGAGTACAGTGCAAGGTCGACAATCAACCTCTGGATTACATCCGATTGCCTCACTGTGCGAAAGTACTCGATGGCGTGGCTTAGAAAGCGTACAGTCTCCGTGCCGGGAAAATAAGAGCGCCTGCGGTTATGAAATCGTGGGCTACTCCTGGGTGGAAAGCTATCCTGCACATTAGTACGAAAGGTGCCAGGTTGCTTCGATCGAAAGCCCGAGAGATCACTCGTAGGAAACTACGCCGGTCACGACGGGCGAAACGACATGAACTCATCCGGACGAAAGGTAGTCCTTACGGTGATCTGCTAGGGTCTCTCCTAGCAACGGTTACTCCATCTGGTACACCCCCTGCTCGGGGCAAGTACCTGATGCGGCACAATGTCTAGCAGGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTTAACCTTAGCAATACGTCTGTCGAAGCAGCTACAA".as_bytes();
        let reference = "CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTAGGAAAGGACAGTGGGAGTGGCACCTTCCAGGGTCAAGGAAGGCACGGGGGAGGGGCAAACAACAGATGGCTGGCAACTAGAAGGCACAGTGAGCTTGTACATAACTACGCAAGTCCTTGCTAGGACCGGCCTTAAAGCCACGTGGCGGCCGCCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCGTGGTACAATATGCGTCTCCGAAATTAACCCGGTGTGTTTAAACGAAAAGGACCGACTACTACCTCGCGAAAGCTCTAAGTGTTGTGTCAGCGAAACTTCGCGGAGGTTCGACATCGAAAGACACGCGGGTGTATATGGCGAAAGCAGCAACCTGATCTGGGGTGAAAAGCCATGGATGTCGGGACGAGAAAGGTCTAGGACTGTTTTGCGAGAAAAGGATTAGAGTTAGAATCGCGAAACGCTCGCGTTCTACCGCTCCGAAAGATCCCGAGGTTGTTTTACCGAAAGCGACGACTTCTGTCATAGTGAAACGATTGGACGTCTCTGGTGCGAAATCGCGGGTTGTACAACATACGAAACCGAGGCTATAATCCCGGACGAAAGGTATAGGTAGCTAACACGCGAAACCCTAGGGATCGTGCTAGCCGAAAGCCCTATTATGTAGGGGACTGAAAAACATGGGTACGTCCCCGATGAAACGCTGCTTGTCTGGCCTCGCGAAAGAATGAGCTGAGTGTGAGGCGAAAAGCTTAAGCTGTGCACTCTCGAAAGTCGGTGTCTATTAGTGGATGAAACAGCGGGTTCCTGCTCCCGCGAAACGCCACCTGTATGTTACTTCGAAAATGAAGGGATAGTGGCGGACGAAAGTCATATTCCGTTGTGGTACGAAATTGGTCCTGATGTACGCACAGAAAAGATTGACCTCTGTTCGTACGAAAGCTCGGCCTCTGGGAGTCGTGAAAGACTCGGATCCGTACCAGATGAAAGGCACACCCATGTCCGTCACGAAAACCCAAACCTTGTATGTATGGAAATCTTCTGCGTTCGGGCCGCGGAAAAGCGTATACCTATCTCGCATGAAAGTCTCTTATCTTGTCTACGCGAAACGCTCGTATGCGTACGGGCTGAAAGCGATATACTGTTCGCCCCTGAAACCCTCTAGTTATGCGCCAGTGAAAGAGTCGCGTAGAGTACAGTGCAAGGTCGACAATCAACCTCTGGATTACATCCGATTGCCTTACTGTGCGAAAGTACTCGATGGTGTGGCTTAGAAAGCGTACAGTCTCTGTGCCGGGAAAATAAGAGCGTCTGCGGTTATGAAATCGTGGGCTACTCCTGGGTGGAAAGCTATCCTGTATATTAGTACGAAAGGTGCCAGGTTGCTTCGATCGAAAGCCCGAGAGATTACTCGTAGGAAACTACGCCGGTTACGACGGGCGAAACGACATGAACTTATCCGGACGAAAGGTAGTCCTTACGGTGATCTGCTAGGGTCTCTCCTAGCAACGGTTACTCGATTTGGTACNNNNNNNNNNNNNNNNNNGTACCTGATGCGGCACAATGTCTAGC".as_bytes();
        let start_pos: usize = 51;
        let mut cigar = Vec::new();
        // 42S 24M 10I 3M 5I 591M 1I 970M 68S
        // 42S 24M 10I 3M 5I 591M 1I 970M 68S
        // 42S 24M 10I 3M 5I 591M 1I 970M 68S
        cigar.push(Op::new(Kind::SoftClip, 42));
        cigar.push(Op::new(Kind::Match, 24));
        cigar.push(Op::new(Kind::Insertion, 10));
        cigar.push(Op::new(Kind::Match, 3));
        cigar.push(Op::new(Kind::Insertion, 6));
        cigar.push(Op::new(Kind::Match, 591));
        cigar.push(Op::new(Kind::Insertion, 1));
        cigar.push(Op::new(Kind::Match, 970));
        cigar.push(Op::new(Kind::SoftClip, 68));

        let recovered = recover_align_sequences(
            read,
            start_pos,
            &cigar,
            &SoftClipResolution::Realign,
            reference,
        );
        println!("Ref\n{}\nread\n{} ",String::from_utf8(recovered.aligned_ref.clone()).unwrap(),String::from_utf8(recovered.aligned_read.clone()).unwrap());
        assert_eq!("-------------TTCCGATCTGTCATAACACCACACTAGAATCACGCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTAGGAAAGGACAGTGGGAGTGGCA--".as_bytes().to_vec(),recovered.aligned_read);
        assert_eq!("CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNTT---TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTAGGAAAGGACAGTGGGAGTGGCACC".as_bytes().to_vec(),recovered.aligned_ref);
*/

    }
}
