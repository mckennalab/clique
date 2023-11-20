use bio::io::fastq::Record;
use needletail::Sequence;
use serde::{Deserialize, Serialize};
use crate::alignment::fasta_bit_encoding::{encoding_to_u8, FASTA_UNSET, FastaBase, reverse_complement};
use crate::alignment::scoring_functions::{AffineScoring, AffineScoringFunction};
use crate::alignment_functions::align_two_strings;
use crate::read_strategies::read_set::{ReadIterator, ReadSetContainer};
use crate::read_strategies::sequence_layout::{AlignedReadOrientation, MergeStrategy, ReadPosition, SequenceLayoutDesign};
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
/// * `reads` - A reference to a `ReadSetContainer` object representing the sequence reads
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
pub fn merge_reads_by_concatenation(reads: &ReadSetContainer, sequence_layout: &SequenceLayoutDesign) -> MergedSequence {

    // aim to have too much capacity, it's faster this way, assuming the same for Vec instead of String (https://github.com/hoodie/concatenation_benchmarks-rs)
    let mut final_sequence: Vec<FastaBase> = Vec::with_capacity(reads.read_one.seq().len() * 4);
    let mut final_sequence_quals: Vec<u8> = Vec::with_capacity(reads.read_one.seq().len() * 4);

    let mut chain_aligned_seq: Option<Vec<FastaBase>> = None;
    let mut chain_aligned_quals: Option<Vec<u8>> = None;


    for read_layout in sequence_layout.reads.clone() {
        match read_layout {
            ReadPosition::READ1 { chain_align, orientation } => {
                match chain_align {
                    Some(x) if x => {
                        let ret = chain_align_update(&reads.read_one, &mut chain_aligned_seq, &mut chain_aligned_quals, &orientation);
                        chain_aligned_seq = Some(ret.0);
                        chain_aligned_quals = Some(ret.1);
                    }
                    _ => {
                        if chain_aligned_seq.is_some() {
                            final_sequence.extend(chain_aligned_seq.unwrap());
                            final_sequence_quals.extend(chain_aligned_quals.unwrap());
                            chain_aligned_seq = None;
                            chain_aligned_quals = None;
                        }
                        final_sequence.extend(sequence_to_fasta_vec(reads.read_one.seq(), &orientation));
                        final_sequence_quals.extend(reads.read_one.qual());
                    }
                }
            }
            ReadPosition::READ2 { chain_align, orientation } => {
                assert!(reads.read_two.is_some());
                match chain_align {
                    Some(x) if x => {
                        let ret = chain_align_update(&reads.read_two.as_ref().unwrap(), &mut chain_aligned_seq, &mut chain_aligned_quals, &orientation);
                        chain_aligned_seq = Some(ret.0);
                        chain_aligned_quals = Some(ret.1);
                    }
                    _ => {
                        if chain_aligned_seq.is_some() {
                            final_sequence.extend(chain_aligned_seq.unwrap());
                            final_sequence_quals.extend(chain_aligned_quals.unwrap());
                            chain_aligned_seq = None;
                            chain_aligned_quals = None;
                        }
                        final_sequence.extend(sequence_to_fasta_vec(reads.read_two.as_ref().unwrap().seq(), &orientation));
                        final_sequence_quals.extend(reads.read_two.as_ref().unwrap().qual());
                    }
                }
            }
            ReadPosition::INDEX1 { chain_align, orientation } => {
                assert!(reads.index_one.is_some());
                match chain_align {
                    Some(x) if x => {
                        let ret = chain_align_update(&reads.index_one.as_ref().unwrap(), &mut chain_aligned_seq, &mut chain_aligned_quals, &orientation);
                        chain_aligned_seq = Some(ret.0);
                        chain_aligned_quals = Some(ret.1);
                    }
                    _ => {
                        if chain_aligned_seq.is_some() {
                            final_sequence.extend(chain_aligned_seq.unwrap());
                            final_sequence_quals.extend(chain_aligned_quals.unwrap());
                            chain_aligned_seq = None;
                            chain_aligned_quals = None;
                        }
                        final_sequence.extend(sequence_to_fasta_vec(reads.index_one.as_ref().unwrap().seq(), &orientation));
                        final_sequence_quals.extend(reads.index_one.as_ref().unwrap().qual());
                    }
                }
            }
            ReadPosition::INDEX2 { chain_align, orientation } => {
                assert!(reads.index_two.is_some());
                match chain_align {
                    Some(x) if x => {
                        let ret = chain_align_update(&reads.index_two.as_ref().unwrap(), &mut chain_aligned_seq, &mut chain_aligned_quals, &orientation);
                        chain_aligned_seq = Some(ret.0);
                        chain_aligned_quals = Some(ret.1);
                    }
                    _ => {
                        if chain_aligned_seq.is_some() {
                            final_sequence.extend(chain_aligned_seq.unwrap());
                            final_sequence_quals.extend(chain_aligned_quals.unwrap());
                            chain_aligned_seq = None;
                            chain_aligned_quals = None;
                        }
                        final_sequence.extend(sequence_to_fasta_vec(reads.index_two.as_ref().unwrap().seq(), &orientation));
                        final_sequence_quals.extend(reads.index_two.as_ref().unwrap().qual());
                    }
                }
            }
            ReadPosition::SPACER { spacer_sequence } => {
                final_sequence.extend(sequence_to_fasta_vec(spacer_sequence.as_bytes(), &AlignedReadOrientation::Forward));
            }
        }
    }
    if chain_aligned_seq.is_some() {
        final_sequence.extend(chain_aligned_seq.unwrap());
        final_sequence_quals.extend(chain_aligned_quals.unwrap());
    }
    //println!("Final sequence: {}", FastaBase::to_string(&final_sequence));
    MergedSequence { read_bases: final_sequence, read_quals: final_sequence_quals }
}

fn chain_align_update(record: &Record, chain_aligned_seq: &Option<Vec<FastaBase>>, chain_aligned_quals: &Option<Vec<u8>>, orientation: &AlignedReadOrientation) -> (Vec<FastaBase>, Vec<u8>) {
    let seq: Vec<u8> = record.seq().to_vec();
    let qual: Vec<u8> = record.qual().to_vec();

    match chain_aligned_seq {
        None => {
            (sequence_to_fasta_vec(&seq, orientation), qual.clone())
        }
        Some(x) => {
            let merged = merge_fasta_bases_by_alignment(&x,
                                                        &chain_aligned_quals.as_ref().unwrap(),
                                                        &sequence_to_fasta_vec(seq.as_slice(), orientation),
                                                        &qual,
                                                        DEFAULT_ALIGNMENT_AFFINE_SCORING.as_ref());
            (merged.read_bases, merged.read_quals)
        }
    }
}

pub fn sequence_to_fasta_vec(sequence: &[u8], orientation: &AlignedReadOrientation) -> Vec<FastaBase> {
    match orientation {
        AlignedReadOrientation::Forward => { FastaBase::from_u8_slice(sequence) }
        AlignedReadOrientation::Reverse => {
            let mut bases = FastaBase::from_u8_slice(sequence);
            bases.reverse();
            bases
        }
        AlignedReadOrientation::ReverseComplement => {
            reverse_complement(&FastaBase::from_u8_slice(sequence))
        }
    }
}

lazy_static! {

    pub static ref DEFAULT_ALIGNMENT_AFFINE_SCORING: Box<dyn AffineScoringFunction + Send + Sync> =
        Box::new(AffineScoring {
            match_score: 10.0,
            mismatch_score: -5.0,
            special_character_score: 8.0,
            gap_open: -15.0,
            gap_extend: -1.0,
            final_gap_multiplier: 0.25,
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
            read_content: (MergedReadSequence::contains_read1(read_structure),
                           MergedReadSequence::contains_read2(read_structure),
                           MergedReadSequence::contains_index1(read_structure),
                           MergedReadSequence::contains_index2(read_structure)),
            read_structure: read_structure.clone(),
        }
    }

    // these functions are very dump, but save us from using a macro crate to do this for us
    fn contains_read1(read_structure: &SequenceLayoutDesign) -> bool {
        read_structure.reads.iter().any(|s| match s {
            ReadPosition::READ1 { chain_align: _, orientation: _ } => true,
            _ => false
        })
    }
    fn contains_read2(read_structure: &SequenceLayoutDesign) -> bool {
        read_structure.reads.iter().any(|s| match s {
            ReadPosition::READ2 { chain_align: _, orientation: _ } => true,
            _ => false
        })
    }
    fn contains_index1(read_structure: &SequenceLayoutDesign) -> bool {
        read_structure.reads.iter().any(|s| match s {
            ReadPosition::INDEX1 { chain_align: _, orientation: _ } => true,
            _ => false
        })
    }
    fn contains_index2(read_structure: &SequenceLayoutDesign) -> bool {
        read_structure.reads.iter().any(|s| match s {
            ReadPosition::INDEX2 { chain_align: _, orientation: _ } => true,
            _ => false
        })
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct UnifiedRead {
    name: Option<Vec<u8>>,
    seq: Option<Vec<FastaBase>>,
    read_structure: SequenceLayoutDesign,
    read_pattern: (bool, bool, bool, bool),
    underlying_reads: ReadSetContainer,

}

impl UnifiedRead {
    pub fn new(read_structure: SequenceLayoutDesign, underlying_reads: ReadSetContainer) -> UnifiedRead {
        UnifiedRead {
            name: None,
            seq: None,
            underlying_reads,
            read_pattern: (MergedReadSequence::contains_read1(&read_structure),
                            MergedReadSequence::contains_read2(&read_structure),
                            MergedReadSequence::contains_index1(&read_structure),
                            MergedReadSequence::contains_index2(&read_structure)),
            read_structure,

        }
    }
    pub fn name(&mut self) -> &Vec<u8> {
        if !self.name.is_some() {
            let _ = &mut self.decision_tree();
        }
        self.name.as_ref().unwrap()
    }
    pub fn seq(&mut self) -> &Vec<FastaBase> {
        if !self.seq.is_some() {
            let _ = &mut self.decision_tree();
        }
        self.seq.as_ref().unwrap()
    }

    fn decision_tree(&mut self) {
        match (self.read_pattern, &self.read_structure.merge) {
            ((true, true, false, false), Some(MergeStrategy::Align)) => {
                self.name = Some(self.underlying_reads.read_one.id().as_bytes().to_vec());
                self.seq = Some(merge_reads_by_alignment(&self.underlying_reads.read_one,
                                                         &self.underlying_reads.read_two.as_ref().unwrap(),
                                                         DEFAULT_ALIGNMENT_AFFINE_SCORING.as_ref()).read_bases);
            }
            ((true, true, false, false), Some(strat)) if strat == &MergeStrategy::Concatenate || strat == &MergeStrategy::ConcatenateBothForward => {
                self.name = Some(self.underlying_reads.read_one.id().as_bytes().to_vec());
                self.seq = Some(merge_reads_by_concatenation(&self.underlying_reads, &self.read_structure).read_bases);
            }
            ((true, true, true, false), Some(strat)) if strat == &MergeStrategy::Concatenate || strat == &MergeStrategy::ConcatenateBothForward => {
                self.name = Some(self.underlying_reads.read_one.id().as_bytes().to_vec());
                self.seq = Some(merge_reads_by_concatenation(&self.underlying_reads, &self.read_structure).read_bases);
            }
            ((true, false, false, false), _) => {
                self.name = Some(self.underlying_reads.read_one.id().as_bytes().to_vec());
                self.seq = Some(self.underlying_reads.read_one.seq().iter().map(|x| FastaBase::from(*x)).collect());
            }
            _ => {
                panic!("We don't support this read structure yet: {:?}", self.read_structure);
            }
        }
    }
}


impl Iterator for MergedReadSequence {
    type Item = UnifiedRead;

    fn next(&mut self) -> Option<Self::Item> {
        let reads = self.underlying_iterator.next();
        match reads {
            Some(reads) => {
                Some(UnifiedRead::new(self.read_structure.clone(), reads))
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
    let read1_seq = FastaBase::from_vec_u8(&read1.seq().to_vec());
    let rev_comp_read2 = FastaBase::from_vec_u8(&read2.seq().reverse_complement().to_vec());
    let mut rev_comp_qual_read2 = read2.qual().to_vec();
    rev_comp_qual_read2.reverse();
    merge_fasta_bases_by_alignment(&read1_seq, &read1.qual().to_vec(), &rev_comp_read2, &rev_comp_qual_read2, merge_initial_scoring)
}

pub fn merge_fasta_bases_by_alignment(read1_seq: &Vec<FastaBase>,
                                      read1_quals: &Vec<u8>,
                                      read2_seq: &Vec<FastaBase>,
                                      read2_quals: &Vec<u8>,
                                      merge_initial_scoring: &dyn AffineScoringFunction) -> MergedSequence {

    let results = align_two_strings(&read1_seq, &read2_seq, merge_initial_scoring, false);

    alignment_rate_and_consensus(
        &results.reference_aligned,
        read1_quals,
        &results.read_aligned,
        read2_quals)
}

pub struct MergedSequence {
    read_bases: Vec<FastaBase>,
    read_quals: Vec<u8>,
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

    debug!("{} {}",FastaBase::to_string(alignment_1),FastaBase::to_string(alignment_2));
    assert_eq!(alignment_1.len(), alignment_2.len());

    for i in 0..alignment_1.len() {
        match (alignment_1[i], alignment_2[i]) {
            (a, b) if a == b => {
                resulting_alignment.push(a.clone());
                resulting_quality_scores.push(combine_phred_scores(&qual_scores1[alignment_1_qual_position], &qual_scores2[alignment_2_qual_position], true));
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
                let final_base = if qual_scores1[alignment_1_qual_position] >= qual_scores2[alignment_2_qual_position] {
                    debug!("A Base {} {} {} {} {} {}", encoding_to_u8(&a) as char, encoding_to_u8(&b) as char, qual_scores1[alignment_1_qual_position] as char, qual_scores2[alignment_2_qual_position] as char, alignment_1_qual_position, alignment_2_qual_position);
                    a.clone()
                } else {
                    debug!("B Base {} {} {} {} {} {}", encoding_to_u8(&a) as char, encoding_to_u8(&b) as char, qual_scores1[alignment_1_qual_position] as char, qual_scores2[alignment_2_qual_position] as char, alignment_1_qual_position, alignment_2_qual_position);
                    b.clone()
                };

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
    }
}

#[cfg(test)]
mod tests {
    use std::cmp::Ordering;
    use std::collections::BTreeMap;
    use super::*;

    fn str_to_fasta_vec(input: &str) -> Vec<FastaBase> {
        FastaBase::from_vec_u8(&input.as_bytes().to_vec())
    }

    fn get_scoring_scheme() -> AffineScoring {
        AffineScoring {
            match_score: 10.0,
            mismatch_score: -5.0,
            special_character_score: 8.0,
            gap_open: -15.0,
            gap_extend: -1.0,
            final_gap_multiplier: 0.25,
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
        zip_and_convert(&FastaBase::to_string(&merged.read_bases),&String::from("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGCCCCCCCCCTTTTTTTTTTTTTTTTTTTTTTTTTT"));

        assert_eq!(merged.read_bases, str_to_fasta_vec("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGCGGAACCCCCCCTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    }

    fn zip_and_convert(str1: &String, str2: &String) {
        assert_eq!(str1.len(),str2.len());
        let mut res1 = Vec::with_capacity(str1.len());
        let mut res2 = Vec::with_capacity(str2.len());
        str1.as_bytes().iter().zip(str2.as_bytes().iter()).for_each(|(x,y)| {
            match (x,y) {
                (x,y) if x == y => {
                    let xUpper = x.to_ascii_uppercase();
                    let yUpper = y.to_ascii_uppercase();
                    res1.push(xUpper);
                    res2.push(yUpper);
                }
                (x,y) if x != y => {
                    let xUpper = x.to_ascii_lowercase();
                    let yUpper = y.to_ascii_lowercase();
                    res1.push(xUpper);
                    res2.push(yUpper);
                }
                _ => {
                    panic!("what?");
                }
            }
        });
        println!("{}\n{}",String::from_utf8(res1).unwrap(),String::from_utf8(res2).unwrap());
    }

    #[test]
    fn read_merger_real_reads_from_meisam() {
        let read1_fwd = "CGAATGTCAAAGTCAATGCGTTAGGGTTTCTTATATGGTGGTTTCTAACATTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCGTACTCTGCGTTGATACCACTGCTTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATG".as_bytes();
        let read1_qls = "AAFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF".as_bytes();
        let read2_fwd = "AATCAGTGGTATAAAAGCAGAGTACTCCTTAGGTTAACTTTCTATTTCTAGCTCTAACCCCAATGTTAGAAACCCCCATATAAGAAACCCTAACGCATTGACTTTGACATTCGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGAT".as_bytes();
        let read2_qls = "=FAF6FFFFFFF//FFFFFFFFFF//FAAAFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFAFFFFF/FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFAFFAFFFFFF//FF/FA/F/F=F//=/".as_bytes();

        let record1 = bio::io::fastq::Record::with_attrs("fakeRead", None, read1_fwd, read1_qls);
        let record2 = bio::io::fastq::Record::with_attrs("fakeRead", None, read2_fwd, read2_qls);


        let merged = merge_reads_by_alignment(&record1, &record2, &get_scoring_scheme());
        let real = "ATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTCGAATGTCAAAGTCAATGCGTTAGGGTTTCTTATATGGTGGTTTCTAACATTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCGTACTCTGCGTTGATACCACTGCTTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATG";
        println!("{}\n{}",FastaBase::to_string(&merged.read_bases),real);
        zip_and_convert(&FastaBase::to_string(&merged.read_bases),&String::from(real.clone()));
        assert_eq!(merged.read_bases, str_to_fasta_vec(real));
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
    use crate::read_strategies::sequence_layout::ReadPosition::{READ1, READ2, SPACER};

    #[test]
    fn read_merger_many_real_reads() {
        let read1_fwd = "GTGGAAAGGACGAAACACCGACGTCTACGTAGACGTACGTTGGAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAACAAGTGGCACCGAGTCGGTGCTTTTTTCGCATTCTACCGTCTTTAGCAAGNTGA".as_bytes();
        let read1_qls = "FFFFFFFFFFAFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF/FFAFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF/F/FFFFFFFFFFAFFFFFFFFFFFFFFFA/F=FFFFFFFFFFFFFFFAFFF#FFF".as_bytes();
        let read2_fwd = "TTGGCCGCGGATCCGATTTAAATTCGAATTCAAACATCGACCTGTACGTCTACGTAGACGTACAGGTCGATACTGTTGCGAATGATCACCTTGCTAAAGTCACGGTAGAATGCGAAAAAAGCACCGACTCGGTGCCACTTTTTCAAGTT".as_bytes();
        let read2_qls = "FFFFFFFFAFFAFFFFF/FFFFFFFFAFFFFFFFFFFFF/FFFFAFFFFFFFFFFFAFFFF/FFFFFFFFAAFFFFFFAFF/FF=FFFFFFFAFFFFFFFFFFFFFFFFFFFFF=FAFFFFFFFFFFFFFFFFFFFFFF=F=FFF=FFF".as_bytes();

        let record1 = bio::io::fastq::Record::with_attrs("fakeRead", None, read1_fwd, read1_qls);
        let record2 = bio::io::fastq::Record::with_attrs("fakeRead", None, read2_fwd, read2_qls);

        let start = Instant::now();

        for _i in 0..1000 {
            merge_reads_by_alignment(&record1, &record2, &get_scoring_scheme());
        }
        let duration = start.elapsed();
        println!("Time elapsed in expensive_function() is: {:?}", duration);
    }

    #[test]
    fn check_orientation() {
        let read1_fwd = "AAAAAAAAAA".as_bytes();
        let read1_qls = "FFFFFFFFFF".as_bytes();
        let read2_fwd = "TTTTTTTTTT".as_bytes();
        let read2_qls = "FFFFFFFFFF".as_bytes();

        let both_fwd = FastaBase::from_vec_u8(&"AAAAAAAAAATTTTTTTTTT".as_bytes().to_vec());
        let both_rc = FastaBase::from_vec_u8(&"AAAAAAAAAAAAAAAAAAAA".as_bytes().to_vec());

        let record1 = bio::io::fastq::Record::with_attrs("fakeRead", None, read1_fwd, read1_qls);
        let record2 = bio::io::fastq::Record::with_attrs("fakeRead", None, read2_fwd, read2_qls);

        // ReadSetContainer, sequence_layout: &SequenceLayoutDesig
        let read_set = ReadSetContainer {
            read_one: record1,
            read_two: Some(record2),
            index_one: None,
            index_two: None,
        };

        let sequence_layout = SequenceLayoutDesign {
            aligner: None,
            merge: None,
            reads: vec![READ1 { chain_align: None, orientation: AlignedReadOrientation::Forward }, READ2 { chain_align: None, orientation: AlignedReadOrientation::ReverseComplement }],
            known_strand: true,
            umi_configurations: BTreeMap::new(),
        };

        let fake_aligned = merge_reads_by_concatenation(&read_set, &sequence_layout);
        assert_eq!(fake_aligned.read_bases.cmp(&both_rc) == Ordering::Equal, true);

        let sequence_layout = SequenceLayoutDesign {
            aligner: None,
            merge: None,
            reads: vec![READ1 { chain_align: None, orientation: AlignedReadOrientation::Forward }, READ2 { chain_align: None, orientation: AlignedReadOrientation::Reverse }],
            known_strand: true,
            umi_configurations: BTreeMap::new(),
        };

        let fake_aligned = merge_reads_by_concatenation(&read_set, &sequence_layout);
        assert_ne!(fake_aligned.read_bases.cmp(&both_rc) == Ordering::Equal, true);


        let sequence_layout = SequenceLayoutDesign {
            aligner: None,
            merge: None,
            reads: vec![READ1 { chain_align: None, orientation: AlignedReadOrientation::Forward }, READ2 { chain_align: None, orientation: AlignedReadOrientation::Forward }],
            known_strand: true,
            umi_configurations: BTreeMap::new(),
        };
        let fake_aligned = merge_reads_by_concatenation(&read_set, &sequence_layout);
        assert_eq!(fake_aligned.read_bases.cmp(&both_fwd) == Ordering::Equal, true);
    }


    #[test]
    fn check_spacers() {
        let read1_fwd = "AAAAAAAAAA".as_bytes();
        let read1_qls = "FFFFFFFFFF".as_bytes();
        let read2_fwd = "TTTTTTTTTT".as_bytes();
        let read2_qls = "FFFFFFFFFF".as_bytes();

        let both_fwd = FastaBase::from_vec_u8(&"AAAAAAAAAAACGTACGTACGTTTTTTTTTTTGGGG".as_bytes().to_vec());

        let record1 = bio::io::fastq::Record::with_attrs("fakeRead", None, read1_fwd, read1_qls);
        let record2 = bio::io::fastq::Record::with_attrs("fakeRead", None, read2_fwd, read2_qls);

        // ReadSetContainer, sequence_layout: &SequenceLayoutDesig
        let read_set = ReadSetContainer {
            read_one: record1,
            read_two: Some(record2),
            index_one: None,
            index_two: None,
        };

        let sequence_layout = SequenceLayoutDesign {
            aligner: None,
            merge: None,
            reads: vec![READ1 { chain_align: None, orientation: AlignedReadOrientation::Forward },
                        SPACER { spacer_sequence: "ACGTACGTACGT".to_string() },
                        READ2 { chain_align: None, orientation: AlignedReadOrientation::Forward },
                        SPACER { spacer_sequence: "GGGG".to_string() }, ],
            known_strand: true,
            umi_configurations: BTreeMap::new(),
        };

        let fake_aligned = merge_reads_by_concatenation(&read_set, &sequence_layout);
        assert_eq!(fake_aligned.read_bases.cmp(&both_fwd) == Ordering::Equal, true);
    }
}