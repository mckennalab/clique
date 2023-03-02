use bio::io::fastq::Record;
use needletail::Sequence;
use crate::alignment::fasta_bit_encoding::{FASTA_UNSET, FastaBase, str_to_fasta_vec};
use crate::alignment::scoring_functions::{AffineScoring};
use crate::alignment_functions::align_two_strings;
use crate::utils::read_utils::combine_phred_scores;

pub fn read_merger(read1: &Record, read2: &Record) -> MergedSequence {
    let read1_seq = str_to_fasta_vec(String::from_utf8(read1.seq().to_vec()).unwrap().as_str());
    let rev_comp_read2 = str_to_fasta_vec(String::from_utf8(read2.seq().reverse_complement().to_vec()).unwrap().as_str());
    let mut rev_comp_read2_qual = read2.qual().to_vec();
    rev_comp_read2_qual.reverse();

    let merge_initial_scoring = AffineScoring{
        match_score: 10.0,
        mismatch_score: -15.0,
        special_character_score: 8.0,
        gap_open: -15.0,
        gap_extend: -15.0,
        final_gap_multiplier: 0.1,
    };

    let results = align_two_strings(&read1_seq, &rev_comp_read2, &merge_initial_scoring, false);

    alignment_rate_and_consensus(
        &results.alignment_string1,
        read1.qual(),
        &results.alignment_string2,
        rev_comp_read2_qual.as_slice())
}

pub struct MergedSequence {
    read_bases: Vec<FastaBase>,
    read_quals: Vec<u8>,
    mismatch_rate: f64,
}
pub fn alignment_rate_and_consensus(alignment_1: &Vec<FastaBase>, qual_scores1: &[u8], alignment_2: &Vec<FastaBase>, qual_scores2: &[u8]) -> MergedSequence {
    let mut resulting_alignment = Vec::new();
    let mut resulting_quality_scores = Vec::new();
    let mut alignment_1_qual_position = 0;
    let mut alignment_2_qual_position = 0;
    let mut agreed_bases = 0;

    assert_eq!(alignment_1.len(), alignment_2.len());

    for i in 0..alignment_1.len() {
        match (alignment_1[i],alignment_2[i]) {
            (a, b) if a == b => {
                resulting_alignment.push(a.clone());
                resulting_quality_scores.push(combine_phred_scores(&qual_scores1[alignment_1_qual_position],&qual_scores2[alignment_2_qual_position],true));
                agreed_bases += 1;
                alignment_1_qual_position += 1;
                alignment_2_qual_position += 1;

            },
            (a, b) if a == FASTA_UNSET => {
                resulting_alignment.push(b.clone());
                resulting_quality_scores.push(qual_scores2[alignment_2_qual_position].clone());
                alignment_2_qual_position += 1;
            },
            (a, b) if b == FASTA_UNSET => {
                resulting_alignment.push(a.clone());
                resulting_quality_scores.push(qual_scores1[alignment_1_qual_position].clone());
                alignment_1_qual_position += 1;
            },
            (a,b) => {
                // bases disagree -- we take the higher quality base
                let final_base = if qual_scores1[alignment_1_qual_position] > qual_scores2[alignment_2_qual_position] {
                    a.clone() } else {b.clone()};

                resulting_alignment.push(final_base);
                resulting_quality_scores.push(combine_phred_scores(&qual_scores1[alignment_1_qual_position],&qual_scores2[alignment_2_qual_position],false));

                alignment_1_qual_position += 1;
                alignment_2_qual_position += 1;
            },
        }
    }

    MergedSequence {
        read_bases: resulting_alignment,
        read_quals: resulting_quality_scores,
        mismatch_rate: agreed_bases as f64 / alignment_1.len() as f64,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
/*
    #[test]
    fn read_merger_simple() {
        let read1_fwd = "AAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGG".as_bytes();
        let read1_qls = "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH".as_bytes();
        let read2_fwd = "AAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGCCCCC".as_bytes();
        let read2_qls = "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH".as_bytes();

        let record1 = bio::io::fastq::Record::with_attrs("fakeRead", None, read1_fwd, read1_qls);
        let record2 = bio::io::fastq::Record::with_attrs("fakeRead", None, read2_fwd, read2_qls);

        let merged = read_merger(&record1, &record2);
        assert_eq!(String::from_utf8(merged.read_bases).unwrap(),String::from("AAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGCCCCCCCCCTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    }

    #[test]
    fn read_merger_simple_no_merge() {
        let read1_fwd = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".as_bytes();
        let read1_qls = "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH".as_bytes();
        let read2_fwd = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG".as_bytes();
        let read2_qls = "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH".as_bytes();

        let record1 = bio::io::fastq::Record::with_attrs("fakeRead", None, read1_fwd, read1_qls);
        let record2 = bio::io::fastq::Record::with_attrs("fakeRead", None, read2_fwd, read2_qls);


        let merged = read_merger(&record1, &record2);
        assert_eq!(String::from_utf8(merged.read_bases).unwrap(),String::from("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"));
    }

    #[test]
    fn read_merger_broken_tail() {
        let read1_fwd = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGCGGAA".as_bytes();
        let read1_qls = "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH".as_bytes();
        let read2_fwd = "AAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGCCCCC".as_bytes();
        let read2_qls = "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH".as_bytes();

        let record1 = bio::io::fastq::Record::with_attrs("fakeRead", None, read1_fwd, read1_qls);
        let record2 = bio::io::fastq::Record::with_attrs("fakeRead", None, read2_fwd, read2_qls);


        let merged = read_merger(&record1, &record2);
        assert_eq!(String::from_utf8(merged.read_bases).unwrap(),String::from("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGCCCCCCCCCTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    }


    #[test]
    fn read_merger_real_reads() {
        let read1_fwd = "GTGGAAAGGACGAAACACCGACGTCTACGTAGACGTACGTTGGAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAACAAGTGGCACCGAGTCGGTGCTTTTTTCGCATTCTACCGTCTTTAGCAAGNTGA".as_bytes();
        let read1_qls = "FFFFFFFFFFAFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF/FFAFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF/F/FFFFFFFFFFAFFFFFFFFFFFFFFFA/F=FFFFFFFFFFFFFFFAFFF#FFF".as_bytes();
        let read2_fwd = "TTGGCCGCGGATCCGATTTAAATTCGAATTCAAACATCGACCTGTACGTCTACGTAGACGTACAGGTCGATACTGTTGCGAATGATCACCTTGCTAAAGTCACGGTAGAATGCGAAAAAAGCACCGACTCGGTGCCACTTTTTCAAGTT".as_bytes();
        let read2_qls = "FFFFFFFFAFFAFFFFF/FFFFFFFFAFFFFFFFFFFFF/FFFFAFFFFFFFFFFFAFFFF/FFFFFFFFAAFFFFFFAFF/FF=FFFFFFFAFFFFFFFFFFFFFFFFFFFFF=FAFFFFFFFFFFFFFFFFFFFFFF=F=FFF=FFF".as_bytes();

        let record1 = bio::io::fastq::Record::with_attrs("fakeRead", None, read1_fwd, read1_qls);
        let record2 = bio::io::fastq::Record::with_attrs("fakeRead", None, read2_fwd, read2_qls);


        let merged = read_merger(&record1, &record2);
        assert_eq!(String::from_utf8(merged.read_bases).unwrap(),String::from("GTGGAAAGGACGAAACACCGACGTCTACGTAGACGTACGTTGGAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTCGCATTCTACCGTGACTTTAGCAAGGTGATCATTCGCAACAGTATCGACCTGTACGTCTACGTAGACGTACAGGTCGATGTTTGAATTCGAATTTAAATCGGATCCGCGGCCAA"));
    }

    use std::time::{Duration, Instant};

    #[test]
    fn read_merger_many_real_reads() {
        let read1_fwd = "GTGGAAAGGACGAAACACCGACGTCTACGTAGACGTACGTTGGAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAACAAGTGGCACCGAGTCGGTGCTTTTTTCGCATTCTACCGTCTTTAGCAAGNTGA".as_bytes();
        let read1_qls = "FFFFFFFFFFAFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF/FFAFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF/F/FFFFFFFFFFAFFFFFFFFFFFFFFFA/F=FFFFFFFFFFFFFFFAFFF#FFF".as_bytes();
        let read2_fwd = "TTGGCCGCGGATCCGATTTAAATTCGAATTCAAACATCGACCTGTACGTCTACGTAGACGTACAGGTCGATACTGTTGCGAATGATCACCTTGCTAAAGTCACGGTAGAATGCGAAAAAAGCACCGACTCGGTGCCACTTTTTCAAGTT".as_bytes();
        let read2_qls = "FFFFFFFFAFFAFFFFF/FFFFFFFFAFFFFFFFFFFFF/FFFFAFFFFFFFFFFFAFFFF/FFFFFFFFAAFFFFFFAFF/FF=FFFFFFFAFFFFFFFFFFFFFFFFFFFFF=FAFFFFFFFFFFFFFFFFFFFFFF=F=FFF=FFF".as_bytes();

        let record1 = bio::io::fastq::Record::with_attrs("fakeRead", None, read1_fwd, read1_qls);
        let record2 = bio::io::fastq::Record::with_attrs("fakeRead", None, read2_fwd, read2_qls);

        let start = Instant::now();

        for i in 0..1000 {
            read_merger(&record1, &record2);
        }
        let duration = start.elapsed();
        println!("Time elapsed in expensive_function() is: {:?}", duration);

    }*/
}