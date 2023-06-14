use bio::io::fastq::Record;
use needletail::Sequence;
use serde::{Deserialize, Serialize};
use crate::alignment::fasta_bit_encoding::{FASTA_UNSET, FastaBase, str_to_fasta_vec};
use crate::alignment::scoring_functions::{AffineScoring, AffineScoringFunction};
use crate::alignment_functions::align_two_strings;
use crate::read_strategies::read_set::{ReadIterator, ReadSetContainer};
use crate::read_strategies::sequence_layout::{MergeStrategy, ReadPosition, SequenceLayoutDesign};
use crate::utils::read_utils::combine_phred_scores;

/// Merges two DNA sequences (read1 and read2) into a single sequence, filling the gap with undefined bases, based on the given reference length.
///
/// This function takes two DNA sequences in the form of `Record` objects (read1 and read2) and a reference length.
/// It first converts the sequences into FASTA vectors and creates the reverse complement of read2.
/// Then, it calculates the gap size based on the reference length and the lengths of the input sequences.
/// Finally, it merges the sequences by extending read1, filling the gap with undefined bases, and extending the reverse complement of read2.
///
/// # Arguments
///
/// * `read1` - A reference to a `Record` object representing the first DNA sequence.
/// * `read2` - A reference to a `Record` object representing the second DNA sequence.
/// * `reference` - A reference sequence
///
/// # Returns
///
/// * `MergedSequence` - A struct containing the merged sequence, quality scores, and mismatch rate.
///
/// # Panics
///
/// * If the gap size is negative, the function will panic with an error message. The prereq for using this function is to not call it with overlapping reads,
/// as that takes additional processing
///
/// # Example
///
/// ```
/// let merged_sequence = merge_reads_by_concatenation(&read1, &read2, &reference, &scoring);
/// ```
pub fn merge_reads_by_concatenation(read1: &Record, read2: &Record, reference_size: usize) -> MergedSequence {
    let read1_seq = FastaBase::from_vec_u8(&read1.seq().to_vec());
    let rev_comp_read2 = FastaBase::from_vec_u8(&read2.seq().reverse_complement().to_vec());

    // simply 'align' them to the ends
    match &reference_size {
        x if x > &0 => {
            let gap_size: i64 = reference_size as i64 - read1_seq.len() as i64 - rev_comp_read2.len() as i64;
            if gap_size < 0 {
                panic!("Our gap size cannot be less than 0 ({} here, ref length {}), for sequence {:?} and {:?}", gap_size, reference_size, read1_seq, rev_comp_read2);
            }

            let mut ret_vec = Vec::with_capacity(*x);
            ret_vec.extend(read1_seq);
            ret_vec.extend(vec![FASTA_UNSET; gap_size.abs() as usize]);
            ret_vec.extend(rev_comp_read2);
            let quals = vec![b'H'; *x];

            MergedSequence {
                read_bases: ret_vec,
                read_quals: quals,
                mismatch_rate: 0.0,
            }
        }
        _ => {
            let mut ret_vec = Vec::with_capacity(read1_seq.len() + rev_comp_read2.len());
            ret_vec.extend(read1_seq);
            ret_vec.extend(rev_comp_read2);
            let quals = vec![b'H'; ret_vec.len()];

            MergedSequence {
                read_bases: ret_vec,
                read_quals: quals,
                mismatch_rate: 0.0,
            }
        }
    }
}
lazy_static! {
    pub static ref DEFAULT_ALIGNMENT_AFFINE_SCORING: Box<dyn AffineScoringFunction + Send + Sync> =
        Box::new(AffineScoring {
            match_score: 10.0,
            mismatch_score: -15.0,
            special_character_score: 8.0,
            gap_open: -15.0,
            gap_extend: -15.0,
            final_gap_multiplier: 0.1,
        });
}

pub struct MergedReadSequence {
    underlying_iterator: ReadIterator,
    read_structure: SequenceLayoutDesign,
    read_content: (bool, bool, bool, bool), // a shortcut to save figuring this out for every iteration
}

impl MergedReadSequence {
    pub fn new(iterator: ReadIterator, read_structure: &SequenceLayoutDesign) -> MergedReadSequence {
        MergedReadSequence {
            underlying_iterator: iterator,
            read_structure: read_structure.clone(),
            read_content: (read_structure.reads.contains(&ReadPosition::READ1),
                           read_structure.reads.contains(&ReadPosition::READ2),
                           read_structure.reads.contains(&ReadPosition::INDEX1),
                           read_structure.reads.contains(&ReadPosition::INDEX2)),
        }
    }

    fn decision_tree(&self, it: ReadSetContainer) -> NamedRead {
        match (&self.read_content, &self.read_structure.merge) {
            ((true, true, false, false), Some(MergeStrategy::Align)) => {
                NamedRead{
                    name: it.read_one.id().as_bytes().to_vec(),
                    seq: merge_reads_by_alignment(&it.read_one, &it.read_two.unwrap(),DEFAULT_ALIGNMENT_AFFINE_SCORING.as_ref()).read_bases,
                }
            }
            ((true, true, false, false), Some(MergeStrategy::Concatenate)) => {
                NamedRead{
                    name: it.read_one.id().as_bytes().to_vec(),
                    seq: merge_reads_by_concatenation(&it.read_one, &it.read_two.unwrap(),0).read_bases
                }
            }
            ((true, false, false, false), _) => {
                NamedRead {
                    name: it.read_one.id().as_bytes().to_vec(),
                    seq: it.read_one.seq().iter().map(|x| FastaBase::from(*x)).collect(),
                }
            }
            _ => {
                panic!("We don't support this read structure yet: {:?}", self.read_structure);
            }
        }
    }
}

#[derive(Debug, PartialEq, Serialize, Deserialize, Clone)]
pub struct NamedRead {
    pub name: Vec<u8>,
    pub seq: Vec<FastaBase>,
}

impl Iterator for MergedReadSequence {
    type Item = NamedRead;

    fn next(&mut self) -> Option<Self::Item> {
        let reads = self.underlying_iterator.next();
        match reads {
            Some(reads) => {
                Some(self.decision_tree(reads))
            }
            None => None,
        }
    }
}


pub fn merged_iterator(iterator: ReadIterator, read_structure: &SequenceLayoutDesign) -> MergedReadSequence {
    MergedReadSequence::new(iterator, read_structure)
}


/// Merges two DNA sequences (read1 and read2) using an alignment algorithm to find the optimal merging point.
///
/// This function takes two DNA sequences in the form of `Record` objects (read1 and read2).
/// It first converts the sequences into FASTA vectors and creates the reverse complement of read2.
/// Then, it calculates the quality scores for read2 and reverses them.
///
/// The function uses the affine alignment algorithm with custom scoring parameters to find the optimal alignment between the two sequences.
/// Finally, it returns a `MergedSequence` object containing the consensus sequence, quality scores, and mismatch rate for the aligned sequences.
///
/// # Arguments
///
/// * `read1` - A reference to a `Record` object representing the first DNA sequence.
/// * `read2` - A reference to a `Record` object representing the second DNA sequence.
/// * `scoring` - The affine scoring scheme to use
///
/// # Returns
///
/// * `MergedSequence` - A struct containing the merged sequence, quality scores, and mismatch rate.
///
/// # Example
///
/// ```
/// let merged_sequence = read_merger(&read1, &read2);
/// ```
pub fn merge_reads_by_alignment(read1: &Record, read2: &Record, merge_initial_scoring: &dyn AffineScoringFunction) -> MergedSequence {
    let read1_seq = str_to_fasta_vec(String::from_utf8(read1.seq().to_vec()).unwrap().as_str());
    let rev_comp_read2 = str_to_fasta_vec(String::from_utf8(read2.seq().reverse_complement().to_vec()).unwrap().as_str());
    let mut rev_comp_read2_qual = read2.qual().to_vec();
    rev_comp_read2_qual.reverse();

    let results = align_two_strings(&read1_seq, &rev_comp_read2, merge_initial_scoring, false);

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

/// Computes the consensus sequence, quality scores, and mismatch rate for two aligned DNA sequences.
///
/// This function takes two aligned DNA sequences (alignment_1 and alignment_2) and their corresponding quality scores (qual_scores1 and qual_scores2).
/// It iterates through the aligned sequences, building the consensus sequence by comparing the bases at each position.
///
/// # Arguments
///
/// * `alignment_1` - A reference to a `Vec<FastaBase>` representing the first aligned DNA sequence.
/// * `qual_scores1` - A reference to a slice of u8 values representing the quality scores for the first sequence.
/// * `alignment_2` - A reference to a `Vec<FastaBase>` representing the second aligned DNA sequence.
/// * `qual_scores2` - A reference to a slice of u8 values representing the quality scores for the second sequence.
///
/// # Returns
///
/// * `MergedSequence` - A struct containing the consensus sequence, quality scores, and mismatch rate.
///
/// # Example
///
/// ```
/// let merged_sequence = alignment_rate_and_consensus(&alignment_1, &qual_scores1, &alignment_2, &qual_scores2);
/// ```
///
/// # Note
///
/// The function assumes that the input sequences have the same length, and it will panic if they do not.
pub fn alignment_rate_and_consensus(alignment_1: &Vec<FastaBase>, qual_scores1: &[u8], alignment_2: &Vec<FastaBase>, qual_scores2: &[u8]) -> MergedSequence {
    let mut resulting_alignment = Vec::new();
    let mut resulting_quality_scores = Vec::new();
    let mut alignment_1_qual_position = 0;
    let mut alignment_2_qual_position = 0;
    let mut agreed_bases = 0;

    assert_eq!(alignment_1.len(), alignment_2.len());

    for i in 0..alignment_1.len() {
        match (alignment_1[i], alignment_2[i]) {
            (a, b) if a == b => {
                resulting_alignment.push(a.clone());
                resulting_quality_scores.push(combine_phred_scores(&qual_scores1[alignment_1_qual_position], &qual_scores2[alignment_2_qual_position], true));
                agreed_bases += 1;
                alignment_1_qual_position += 1;
                alignment_2_qual_position += 1;
            }
            (a, b) if a == FASTA_UNSET => {
                resulting_alignment.push(b.clone());
                resulting_quality_scores.push(qual_scores2[alignment_2_qual_position].clone());
                alignment_2_qual_position += 1;
            }
            (a, b) if b == FASTA_UNSET => {
                resulting_alignment.push(a.clone());
                resulting_quality_scores.push(qual_scores1[alignment_1_qual_position].clone());
                alignment_1_qual_position += 1;
            }
            (a, b) => {
                // bases disagree -- we take the higher quality base
                let final_base = if qual_scores1[alignment_1_qual_position] > qual_scores2[alignment_2_qual_position] {
                    a.clone()
                } else { b.clone() };

                resulting_alignment.push(final_base);
                resulting_quality_scores.push(combine_phred_scores(&qual_scores1[alignment_1_qual_position], &qual_scores2[alignment_2_qual_position], false));

                alignment_1_qual_position += 1;
                alignment_2_qual_position += 1;
            }
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

    fn get_scoring_scheme() -> AffineScoring {
        AffineScoring {
            match_score: 10.0,
            mismatch_score: -15.0,
            special_character_score: 8.0,
            gap_open: -15.0,
            gap_extend: -15.0,
            final_gap_multiplier: 0.1,
        }
    }

    #[test]
    fn read_merger_simple() {
        let read1_fwd = "AAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGG".as_bytes();
        let read1_qls = "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH".as_bytes();
        let read2_fwd = "AAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGCCCCC".as_bytes();
        let read2_qls = "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH".as_bytes();

        let record1 = bio::io::fastq::Record::with_attrs("fakeRead", None, read1_fwd, read1_qls);
        let record2 = bio::io::fastq::Record::with_attrs("fakeRead", None, read2_fwd, read2_qls);

        let merged = merge_reads_by_alignment(&record1, &record2, &get_scoring_scheme());
        assert_eq!(merged.read_bases, str_to_fasta_vec("AAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGCCCCCCCCCTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    }

    #[test]
    fn read_merger_simple_no_merge() {
        let read1_fwd = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".as_bytes();
        let read1_qls = "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH".as_bytes();
        let read2_fwd = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG".as_bytes();
        let read2_qls = "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH".as_bytes();

        let record1 = bio::io::fastq::Record::with_attrs("fakeRead", None, read1_fwd, read1_qls);
        let record2 = bio::io::fastq::Record::with_attrs("fakeRead", None, read2_fwd, read2_qls);


        let merged = merge_reads_by_alignment(&record1, &record2, &get_scoring_scheme());
        assert_eq!(merged.read_bases, str_to_fasta_vec("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"));
    }

    #[test]
    fn read_merger_broken_tail() {
        let read1_fwd = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGCGGAA".as_bytes();
        let read1_qls = "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH".as_bytes();
        let read2_fwd = "AAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGCCCCC".as_bytes();
        let read2_qls = "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH".as_bytes();

        let record1 = bio::io::fastq::Record::with_attrs("fakeRead", None, read1_fwd, read1_qls);
        let record2 = bio::io::fastq::Record::with_attrs("fakeRead", None, read2_fwd, read2_qls);


        let merged = merge_reads_by_alignment(&record1, &record2, &get_scoring_scheme());
        assert_eq!(merged.read_bases, str_to_fasta_vec("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGCCCCCCCCCTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    }


    #[test]
    fn read_merger_real_reads() {
        let read1_qls = "FFFFFFFFFFAFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF/FFAFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF/F/FFFFFFFFFFAFFFFFFFFFFFFFFFA/F=FFFFFFFFFFFFFFFAFFF#FFF".as_bytes();
        let read1_fwd = "GTGGAAAGGACGAAACACCGACGTCTACGTAGACGTACGTTGGAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAACAAGTGGCACCGAGTCGGTGCTTTTTTCGCATTCTACCGTCTTTAGCAAGNTGA".as_bytes();
        let read2_fwd = "TTGGCCGCGGATCCGATTTAAATTCGAATTCAAACATCGACCTGTACGTCTACGTAGACGTACAGGTCGATACTGTTGCGAATGATCACCTTGCTAAAGTCACGGTAGAATGCGAAAAAAGCACCGACTCGGTGCCACTTTTTCAAGTT".as_bytes();
        let read2_qls = "FFFFFFFFAFFAFFFFF/FFFFFFFFAFFFFFFFFFFFF/FFFFAFFFFFFFFFFFAFFFF/FFFFFFFFAAFFFFFFAFF/FF=FFFFFFFAFFFFFFFFFFFFFFFFFFFFF=FAFFFFFFFFFFFFFFFFFFFFFF=F=FFF=FFF".as_bytes();

        let record1 = bio::io::fastq::Record::with_attrs("fakeRead", None, read1_fwd, read1_qls);
        let record2 = bio::io::fastq::Record::with_attrs("fakeRead", None, read2_fwd, read2_qls);


        let merged = merge_reads_by_alignment(&record1, &record2, &get_scoring_scheme());
        assert_eq!(merged.read_bases, str_to_fasta_vec("GTGGAAAGGACGAAACACCGACGTCTACGTAGACGTACGTTGGAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTCGCATTCTACCGTGACTTTAGCAAGGTGATCATTCGCAACAGTATCGACCTGTACGTCTACGTAGACGTACAGGTCGATGTTTGAATTCGAATTTAAATCGGATCCGCGGCCAA"));
    }

    use std::time::{Instant};

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
            merge_reads_by_alignment(&record1, &record2, &get_scoring_scheme());
        }
        let duration = start.elapsed();
        println!("Time elapsed in expensive_function() is: {:?}", duration);
    }

    #[test]
    fn space_to_reference_tester() {
        let reference = str_to_fasta_vec("TTACAAAATTTGTGAAAGATTGACTGGTATTCTTAACTATGTTGCTCCTTTTACGCTATGTGGATACGCTGCTTTAATGCCTTTGTATCATGCGTT");
        let reads_reference_aligned = str_to_fasta_vec("TTACAAAATTTGTGAAAGATTGACTGGTATTCTTAACTATGTTGCTCCTTTT---------------------TTAATGCCTTTGTATCATGCGTT");
        let read1_fwd = "TTACAAAATTTGTGAAAGATTGACTGGTATTCTTAACTATGTTGCTCCTTTT".as_bytes();
        let read1_qls = "FFFFFFFFFFAFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF/FFAFFFFFFF".as_bytes();
        let read2_fwd = "AACGCATGATACAAAGGCATTAA".as_bytes();
        let read2_qls = "FFFFFFFFAFFAFFFFF/FFFFF".as_bytes();

        let record1 = bio::io::fastq::Record::with_attrs("fakeRead", None, read1_fwd, read1_qls);
        let record2 = bio::io::fastq::Record::with_attrs("fakeRead", None, read2_fwd, read2_qls);

        let fake_aligned = merge_reads_by_concatenation(&record1, &record2, &reference);
        assert_eq!(fake_aligned.read_bases, reads_reference_aligned);
    }
}