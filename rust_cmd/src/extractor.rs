use crate::alignment::fasta_bit_encoding::FastaBase;
use nohash_hasher::NoHashHasher;
use noodles_bam::record::Cigar;
use noodles_sam::alignment::record::cigar::op::*;
use std::collections::{BTreeMap, HashMap, VecDeque};
use std::hash::BuildHasherDefault;

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

pub fn recover_align_sequences(
    unaligned_read: &[u8],
    start_pos: usize,
    cigar: &Cigar,
    soft_clip_as_match: &bool,
    reference: &[u8],
) -> RecoveredAlignedSeq {
    
    let mut aligned_read = Vec::new();
    let mut aligned_ref = Vec::new();
    let mut read_pos = 0;
    let mut ref_pos = start_pos - 1;

    if start_pos > 1 {
        // usually we fill in with dashes -- unless it's a clipped read, we deal with that down below
        let first_cigar = cigar.iter().next().unwrap().unwrap();
        if first_cigar.kind() != Kind::SoftClip  {
            aligned_ref.extend(reference[0..start_pos - 1 ].to_vec());
            aligned_read.extend(b"-".repeat(start_pos - 1));
        };
    };
        
    //println!("read length {} cigar: {:?}", unaligned_read.len(), cigar, );
    for (cigar_index, cigar_sub) in cigar.iter().enumerate() {
        let cigar_v = cigar_sub.unwrap();
        let len = cigar_v.len();
        //println!("read pos {} ref pos: {}", read_pos, ref_pos);


        match cigar_v.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
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
                if *soft_clip_as_match {
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
                } else {
                    aligned_read.extend(unaligned_read[read_pos..read_pos + len].to_vec());
                    aligned_ref.extend(b"-".repeat(len));
                    read_pos += len;
                }
            }
            Kind::HardClip | Kind::Pad => {
                // do nothing
            }
        }
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
            FastaBase::valid(reference_base),
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
        .map(|(key, value)| (key.clone(), String::from_utf8(value.clone()).unwrap()))
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
) -> (bool, VecDeque<(char, Vec<FastaBase>)>) {
    let mut invalid_read = false;
    let queue = VecDeque::from(
        reference_tags
            .umi_configurations
            .iter()
            .map(|(_umi_name, umi_obj)| {
                let ets_hit = ets.get(&umi_obj.symbol.to_string().as_bytes()[0]);
                match ets_hit {
                    Some(e) => {
                        if e.len() != umi_obj.length {
                            invalid_read = true;
                        };

                        Some((
                            umi_obj.symbol,
                            e.as_bytes()
                                .iter()
                                .map(|f| FastaBase::from(f.clone()))
                                .collect::<Vec<FastaBase>>(),
                        ))
                    }
                    None => {
                        invalid_read = true;
                        None
                    }
                }
            })
            .filter(|x| x.is_some())
            .map(|x| x.unwrap())
            .collect::<Vec<(char, Vec<FastaBase>)>>(),
    );
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
}
