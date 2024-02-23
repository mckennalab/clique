use std::{str, cmp};

use crate::fasta_comparisons::DEGENERATEBASES;

use std::convert::TryFrom;
use itertools::Itertools;
use crate::alignment::alignment_matrix::{AlignmentCigar, AlignmentTag, MatchedPosition, SharedSegments};
use crate::alignment::fasta_bit_encoding::{FASTA_UNSET, FastaBase};
use crate::alignment::scoring_functions::AffineScoring;
use crate::reference::fasta_reference::SuffixTableLookup;

/// find a read's orientation using the greatest total number of matching bases
///
/// # Arguments
///
/// * `search_string` - a u8 Vec representing the search string
/// * `reference` - a u8 Vec representing the reference string
/// * `seeds` - a suffix array lookup object
pub fn orient_by_longest_segment(search_string: &Vec<FastaBase>, reference: &Vec<u8>, seeds: &SuffixTableLookup) -> (bool, SharedSegments, SharedSegments) {
    let search_string = FastaBase::vec_u8(search_string);
    let fwd_score_mp = find_greedy_non_overlapping_segments(&search_string, reference, seeds);
    let fwd_score: usize = fwd_score_mp.alignment_segments.clone().into_iter().map(|p| p.length).sum();

    let rev_score_mp =  find_greedy_non_overlapping_segments(&bio::alphabets::dna::revcomp(search_string), reference, seeds);
    let rev_score: usize = rev_score_mp.alignment_segments.clone().into_iter().map(|p| p.length).sum();

    (fwd_score > rev_score,fwd_score_mp,rev_score_mp)
}


/// create an alignment string, with matching lengths, from a read, reference, and their CIGAR string
///
/// # Arguments
///
/// * `search_string` - a u8 Vec representing the search string
/// * `reference` - a u8 Vec representing the reference string
/// * `alignment` - contains the starting position and the alignment CIGAR strings
#[allow(dead_code)]
pub fn cigar_alignment_to_full_string(read: &Vec<u8>, reference: &Vec<u8>, alignment: &AlignmentCigar) -> (String, String) {
    let mut read_align = String::new();
    let mut ref_align = String::new();
    let mut current_read_pos = 0;
    let mut current_ref_pos = 0;

    // pad the beginning of the aligment up until the start point
    read_align.push_str(String::from_utf8(vec![b'-'; alignment.alignment_start]).unwrap().as_str());
    ref_align.push_str(str::from_utf8(&reference[0..alignment.alignment_start]).unwrap());
    current_ref_pos += alignment.alignment_start;


    // now process the CIGAR string's individual tokens
    for token in alignment.alignment_tags.clone() {
        match token {
            AlignmentTag::MatchMismatch(size) => {
                read_align.push_str(str::from_utf8(&read[current_read_pos..(current_read_pos + size)]).unwrap());
                ref_align.push_str(str::from_utf8(&reference[current_ref_pos..(current_ref_pos + size)]).unwrap());
                current_read_pos += size;
                current_ref_pos += size;
            }
            AlignmentTag::Del(size) => {
                read_align.push_str(String::from_utf8(vec![b'-'; size]).unwrap().as_str());
                ref_align.push_str(str::from_utf8(&reference[current_ref_pos..(current_ref_pos + size)]).unwrap());
                current_ref_pos += size;
            }
            AlignmentTag::Ins(size) => {
                ref_align.push_str(String::from_utf8(vec![b'-'; size]).unwrap().as_str());
                read_align.push_str(str::from_utf8(&read[current_read_pos..(current_read_pos + size)]).unwrap());
                current_read_pos += size;
            }
            AlignmentTag::InversionOpen => {
                panic!("unclear how to handle InversionOpen");
            }
            AlignmentTag::InversionClose => {
                panic!("unclear how to handle InversionClose");
            }
            AlignmentTag::SoftClip(_) => {
                panic!("unclear how to handle SoftClip");
            }
        }
    }
    (read_align,ref_align)

}


/// find a series of exact matches between the search string and the reference
///
/// # Arguments
///
/// * `search_string` - a u8 Vec representing the search string
/// * `reference` - a u8 Vec representing the reference string
/// * `seeds` - a suffix array lookup object
pub fn find_greedy_non_overlapping_segments(search_string: &[u8], reference: &[u8], seeds: &SuffixTableLookup) -> SharedSegments {
    let mut return_hits: Vec<MatchedPosition> = Vec::new();
    let mut position = 0;
    let mut least_ref_pos = reference.len() as usize;
    let mut greatest_ref_pos = 0;

    while (position as i64) < (search_string.len() as i64 - seeds.seed_size as i64) {
        let ref_positions = seeds.suffix_table.positions_internal(&search_string[position..(position + seeds.seed_size)]);
        let mut longest_hit = 0;
        for ref_position in ref_positions {
            //println!("testing position {}", ref_position);
            if ref_position >= &greatest_ref_pos {
                let extended_hit_size = extend_hit(search_string, position, reference, *ref_position as usize);
                if extended_hit_size > longest_hit {
                    return_hits.push(MatchedPosition { search_start: position, ref_start: *ref_position as usize, length: extended_hit_size });
                    position += extended_hit_size;
                    least_ref_pos = cmp::min(usize::try_from(*ref_position).unwrap(),least_ref_pos);
                    greatest_ref_pos = cmp::max(ref_position + &(extended_hit_size as u32),greatest_ref_pos);
                    longest_hit = extended_hit_size;
                    //println!("taking position {} with greatest_ref_pos {} least {} longest {}", ref_position, greatest_ref_pos, least_ref_pos, longest_hit);
                }
            }
        }
        position += 1;
    }
    SharedSegments { start_position: least_ref_pos as usize, alignment_segments: return_hits}
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlignmentResults {
    pub aligned_read: Vec<FastaBase>,
    pub aligned_ref: Vec<FastaBase>,
    pub cigar_tags: Vec<AlignmentTag>,
}


/*
/// find a series of exact matches between the search string and the reference, and then align the
/// sequences between those exact matches using an inversion aware aligner
///
/// # Arguments
///
/// * `search_string` - a u8 Vec representing the search string
/// * `reference` - a u8 Vec representing the reference string
/// * `seeds` - a suffix array lookup object
#[allow(dead_code)]
pub fn align_string_with_anchors(search_string: &Vec<FastaBase>,
                                 reference: &Vec<FastaBase>,
                                 overlaps: &SharedSegments,
                                 my_inv_score: Option<&InversionScoring>,
                                 my_aff_score: &AffineScoring,
                                 alignment_mat: &mut Alignment<Ix3>) -> AlignmentResult {

    debug!("ref {} read {}",FastaBase::string(reference),FastaBase::string(search_string));

    let mut alignment_ref: Vec<FastaBase> = Vec::new();
    let mut alignment_read: Vec<FastaBase> = Vec::new();
    let mut alignment_cigar = Vec::new();
    let mut read_alignment_last_position: usize = 0;
    let mut ref_alignment_last_position: usize = 0;

    for overlap in &overlaps.alignment_segments {
        assert!(read_alignment_last_position <= overlap.search_start,"READ START FAILURE: {} and {}",read_alignment_last_position,overlap.search_start);
        assert!(ref_alignment_last_position <= overlap.ref_start,"REF START FAILURE: {} and {} from {}",ref_alignment_last_position,overlap.ref_start,overlap.length);

        // look back to see what segment we haven't aligned in the read
        let read_slice = slice_for_alignment(search_string, read_alignment_last_position, overlap.search_start);
        let ref_slice  = slice_for_alignment(reference, ref_alignment_last_position, overlap.ref_start);

        let alignment =
            match (read_slice.len(),ref_slice.len(), my_inv_score) {
                (x ,y, None) if x < 5 && y < 5 && x == y => {
                    AlignmentResult::from_match_segment(&ref_slice, &read_slice, ref_alignment_last_position,read_alignment_last_position,my_aff_score)
                },
                (_x, _y ,Some(inv_score)) => inversion_alignment(&ref_slice, &read_slice, inv_score, my_aff_score, false),
                (_x, _y ,None) => {
                    perform_affine_alignment(alignment_mat, &ref_slice, &read_slice, my_aff_score);

                    perform_3d_global_traceback(alignment_mat, None, &ref_slice, &read_slice, None)
                },
            };
        read_alignment_last_position += read_slice.len();
        ref_alignment_last_position += ref_slice.len();

        debug!("Pushing {:?}",alignment.cigar_string.clone());
        alignment_ref.extend(alignment.reference_aligned);
        alignment_read.extend(alignment.read_aligned);
        alignment_cigar.extend(alignment.cigar_string.into_iter().rev().collect::<Vec<AlignmentTag>>());

        alignment_ref.extend_from_slice(&reference[overlap.ref_start..overlap.ref_start+overlap.length]);
        alignment_read.extend_from_slice(&search_string[overlap.search_start..overlap.search_start+overlap.length]);
        // now add the matching segment

        debug!("Pushing {:?}",AlignmentTag::MatchMismatch(overlap.length));
        read_alignment_last_position += overlap.length;
        ref_alignment_last_position += overlap.length;
        alignment_cigar.push(AlignmentTag::MatchMismatch(overlap.length));
    }

    // process the last segment
    if overlaps.alignment_segments.len() > 0 {
        let read_stop = overlaps.alignment_segments[overlaps.alignment_segments.len() - 1].search_start + overlaps.alignment_segments[overlaps.alignment_segments.len() - 1].length;
        if read_stop < search_string.len() {

            // look back to see what segment we haven't aligned in the read
            let read_slice = slice_for_alignment(&search_string, read_alignment_last_position, search_string.len());
            let ref_slice = slice_for_alignment(&reference, ref_alignment_last_position, reference.len());
            let alignment=
            match my_inv_score {
                Some(x) => inversion_alignment(&ref_slice, &read_slice, x, my_aff_score,false),
                None => {
                    let mut alignment_mat = create_scoring_record_3d(ref_slice.len() + 1, read_slice.len() + 1, AlignmentType::Affine, false);
                    perform_affine_alignment(&mut alignment_mat, &ref_slice, &read_slice, my_aff_score);

                    perform_3d_global_traceback(&mut alignment_mat, None, &ref_slice, &read_slice, None)
                }
            };
            debug!("Pushing {:?}",alignment.cigar_string.clone());

            alignment_ref.extend(alignment.reference_aligned);
            alignment_read.extend(alignment.read_aligned);
            alignment_cigar.extend(alignment.cigar_string.into_iter().rev().collect::<Vec<AlignmentTag>>());

        } else if ref_alignment_last_position < reference.len() {
            let gap_len = reference.len() - ref_alignment_last_position;
            alignment_ref.extend(reference[ref_alignment_last_position..reference.len()].to_vec());
            alignment_read.extend(vec![FASTA_UNSET; gap_len]);
            debug!("Pushing {:?}",Del(gap_len));

            alignment_cigar.push(Del(gap_len));
        }
    } else {
        let alignment=
            match my_inv_score {
                Some(x) => {
                    inversion_alignment(&reference, &search_string, x, my_aff_score,true)
                }
                None => {
                    let mut alignment_mat = create_scoring_record_3d(reference.len() + 1, search_string.len() + 1, AlignmentType::Affine, false);
                    perform_affine_alignment(&mut alignment_mat, &reference, &search_string, my_aff_score);

                    perform_3d_global_traceback(&mut alignment_mat, None, &reference, &search_string, None)
                }
            };
        alignment_ref.extend(alignment.reference_aligned);
        alignment_read.extend(alignment.read_aligned);
        debug!("Pushing {:?}",alignment.cigar_string.clone());
        alignment_cigar.extend(alignment.cigar_string.into_iter().rev().collect::<Vec<AlignmentTag>>());
    }

    let score = calculate_score_from_strings(&alignment_ref, &alignment_read, my_aff_score);
    validate_cigar_string(&alignment_ref, &alignment_read, &alignment_cigar);

    AlignmentResult {
        reference_aligned: alignment_ref,
        read_aligned: alignment_read,
        cigar_string: simplify_cigar_string(&alignment_cigar),
        path: vec![],
        score,
        reference_start: 0,
        read_start: 0,
        bounding_box: None,
    }
}
*/
pub fn validate_cigar_string(reference: &Vec<FastaBase>, read: &Vec<FastaBase>, cigars: &Vec<AlignmentTag>) {
    assert_eq!(reference.len(),read.len());
    debug!("CIGARS: {:?}",cigars);

    let mut cigar_pos = 0;
    cigars.iter().for_each(|c| {
        debug!("CIGAR: {}",c.clone());
       match c {
           AlignmentTag::MatchMismatch(length) => {
               assert_eq!(reference[cigar_pos..cigar_pos + length].iter().counts().get(&FASTA_UNSET).unwrap_or(&0),&0,"CIGAR failure on reference (M({})): {}",length, FastaBase::string(&reference[cigar_pos..cigar_pos + length].to_vec()));
               assert_eq!(read[cigar_pos..cigar_pos + length].iter().counts().get(&FASTA_UNSET).unwrap_or(&0),&0,"CIGAR failure on read (M({})): {}",length, FastaBase::string(&read[cigar_pos..cigar_pos + length].to_vec()));
               cigar_pos += length;
           }
           AlignmentTag::Del(length) => {
               debug!("bit {}",FastaBase::string(&read[cigar_pos..cigar_pos + length].to_vec()));
               assert_eq!(reference[cigar_pos..cigar_pos + length].iter().counts().get(&FASTA_UNSET).unwrap_or(&0),&0,"CIGAR failure on reference (D({})): {}",length, FastaBase::string(&reference[cigar_pos..cigar_pos + length].to_vec()));
               assert_eq!(read[cigar_pos..cigar_pos + length].iter().counts().get(&FASTA_UNSET).unwrap_or(&0),length,"CIGAR failure on read (D({})): {}",length, FastaBase::string(&read[cigar_pos..cigar_pos + length].to_vec()));
               cigar_pos += length;
           }
           AlignmentTag::Ins(length) => {
               assert_eq!(reference[cigar_pos..cigar_pos + length].iter().counts().get(&FASTA_UNSET).unwrap_or(&0),length,"CIGAR failure on reference (I({})): {}",length, FastaBase::string(&reference[cigar_pos..cigar_pos + length].to_vec()));
               assert_eq!(read[cigar_pos..cigar_pos + length].iter().counts().get(&FASTA_UNSET).unwrap_or(&0),&0,"CIGAR failure on reference (I({})): {}",length, FastaBase::string(&read[cigar_pos..cigar_pos + length].to_vec()));
               cigar_pos += length;
           }
           AlignmentTag::SoftClip(length) => {
               assert_eq!(reference[cigar_pos..cigar_pos + length].iter().counts().get(&FASTA_UNSET).unwrap_or(&0),&0,"CIGAR failure on reference (S({})): {}",length, FastaBase::string(&reference[cigar_pos..cigar_pos + length].to_vec()));
               assert_eq!(read[cigar_pos..cigar_pos + length].iter().counts().get(&FASTA_UNSET).unwrap_or(&0),length,"CIGAR failure on reference (S({})): {}",length, FastaBase::string(&read[cigar_pos..cigar_pos + length].to_vec()));
               cigar_pos += length;
           }
           AlignmentTag::InversionOpen => {}
           AlignmentTag::InversionClose => {}
       }
    });
    assert_eq!(cigar_pos,reference.len());
}

pub fn calculate_score_from_strings(reference: &Vec<FastaBase>, read: &Vec<FastaBase>, my_aff_score: &AffineScoring) -> f64 {
    assert_eq!(reference.len(),read.len());
    let mut in_indel = false;
    reference.iter().zip(read.iter()).map(|(a,b)| {
        match (a,b, in_indel) {
            (a, b, _) if !FASTA_UNSET.identity(a) && FASTA_UNSET.identity(b) => {
                in_indel = false;
                my_aff_score.match_mismatch(a, b)
            },
            (_, _, true) => {
                my_aff_score.gap_extend()
            },
            (_, _, false) => {
                in_indel = true;
                my_aff_score.gap_open()
            },
        }
    }).sum()
}

pub fn slice_for_alignment(read: &Vec<FastaBase>, read_start: usize, read_stop: usize) -> Vec<FastaBase> {
    assert!(read_stop <= read.len(),"Read position requested {} when our length is only {} for read {} ",read_stop,read.len(), FastaBase::string(read));
    read[read_start..read_stop].to_vec()
}

/// Extend a seed hit within the reference to its maximum length, using degenerate base matching
pub fn extend_hit(search_string: &[u8], search_location: usize, reference: &[u8], reference_location: usize) -> usize {
    let mut current_length = 0;
    while current_length + search_location < search_string.len() && current_length + reference_location < reference.len() {
        let search_loc = current_length + search_location;
        let ref_loc = current_length + reference_location;

        match (DEGENERATEBASES.get(&search_string[search_loc]),
            DEGENERATEBASES.get(&reference[ref_loc])) {
            (None, None) => {return current_length},
            (None, _) => {return current_length},
            (_, None) => {return current_length},
            (x,y) => {
                match (x.unwrap().contains_key(&reference[ref_loc]),
                       y.unwrap().contains_key(&search_string[search_loc])) {
                    (true,true) => {current_length += 1;}
                    (_, _) => { return current_length },
                }
            }
        }
    }
    current_length
}

#[cfg(test)]
mod tests {
    use crate::reference::fasta_reference::ReferenceManager;
    use super::*;


    #[test]
    fn orient_by_longest_segment_test() {
        let reference = String::from("AAAAATATATATATATAT").as_bytes().to_owned();
        let test_read = String::from("AAAAAGGGGGGGGGGGGG").as_bytes().to_owned();
        let reference_lookup = ReferenceManager::find_seeds(&reference, 5);

        let aligned_string = orient_by_longest_segment(&FastaBase::from_vec_u8(&test_read), &reference, &reference_lookup);
        assert_eq!(aligned_string.1.alignment_segments.len(),1);
        assert_eq!(aligned_string.1.alignment_segments.get(0).unwrap().search_start,0);

        let reference = String::from("AAAAATATATATATATATCCACC").as_bytes().to_owned();
        let test_read = String::from("AAAAAGGGGGGGGGGGGGCCACC").as_bytes().to_owned();
        let reference_lookup = ReferenceManager::find_seeds(&reference, 5);

        let aligned_string = orient_by_longest_segment(&FastaBase::from_vec_u8(&test_read), &reference, &reference_lookup);

        println!("{:?}",reference_lookup.suffix_table);

        assert_eq!(aligned_string.1.alignment_segments.len(),2);
        assert_eq!(aligned_string.1.alignment_segments.get(0).unwrap().search_start,0);
        assert_eq!(aligned_string.1.alignment_segments.get(1).unwrap().search_start,18);
    }

    #[test]
    fn simple_extend_test() {
        let reference = String::from("AATGATACGG").as_bytes().to_owned();
        let test_read = String::from("AATGATACGG").as_bytes().to_owned();

        let aligned_string = extend_hit(&test_read, 0, &reference, 0);
        print!("BLAH {}", aligned_string);
        assert_eq!(aligned_string, 10);
    }

    #[test]
    fn simple_extend_stop_before_end_test() {
        let reference = String::from("AATGATACGGAAA").as_bytes().to_owned();
        let test_read = String::from("AATGATACGG").as_bytes().to_owned();

        let aligned_string = extend_hit(&test_read, 0, &reference, 0);
        print!("BLAH {}", aligned_string);
        assert_eq!(aligned_string, 10);
    }

    #[test]
    fn simple_extend_internal_test() {
        let reference = String::from("GGAATGATACGGAAA").as_bytes().to_owned();
        let test_read = String::from("AATGATACGG").as_bytes().to_owned();

        let aligned_string = extend_hit(&test_read, 0, &reference, 2);
        print!("BLAH {}", aligned_string);
        assert_eq!(aligned_string, 10);
    }
    fn str_to_fasta_vec(input: &str) -> Vec<FastaBase> {
        FastaBase::from_vec_u8(&input.as_bytes().to_vec())
    }

    #[test]
    fn simple_extend_short_test() {
        let reference = String::from("AAA").as_bytes().to_owned();
        let test_read = String::from("AATGATACGG").as_bytes().to_owned();

        let aligned_string = extend_hit(&test_read, 0, &reference, 0);
        print!("BLAH {}", aligned_string);
        assert_eq!(aligned_string, 2);
    }

    #[test]
    fn suffix_array_test() {
        let refseq = String::from("AATGATACGG").as_bytes().to_owned();
        let reference = ReferenceManager::find_seeds(&refseq, 20);
        assert!(reference.suffix_table.contains("AAT"));
        assert!(!reference.suffix_table.contains("TAAT"));
    }

    #[test]
    fn find_greedy_non_overlapping_segments_test() {
        let refseq = String::from("NNNNNNNNCATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCATATGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCGTAGAAGAAAGTTAGCTCTTGTGCGAGCTACAGGAACGATGTTTGATTAGAGTAAGCAGAGGACAAGGGCTCGCGTGCAGCCGAAGTTTGGCCGGTACTCTCCAACCGTTAACAACAACACCTTTCATCGAAATCCGCTTGGTAACAACACTAGGATGTTTCGACGCACTCGATAACCGGGAAACCAAGAGAAGTTTCCGAGCGCCACAGCCCAACTTAACTTCGCCATGTTTGAGACACGCGATCGGGACCACAAGACGTTACTCTTTGGGACCGCCGTAAGCGAGTATAAGCAGATTGTGTTTCGGCGTCCAAGTTGCCGTCAAAAGCTTACTGAGTTTCGCTGCCGGCGGAATGCTATTAATCGCGCCTACTTTCATGGAAGACGTTTCCGCTAAAATGGACTGGTGTTTCATGTCGGGAGCCGCTTTGTAAGATGGAGCTACTTTCCAGTCTGAGTTCATGCGAGAACACCACGAGAGTTTCATGAGTGGCCTCTCGAATCAACAGTCTACAAGTTTGGAGTTATCCGACACATCAAAACCAGCCATTCGTTTCATGAGATGGATCGCATACTAACCTTAACGGAGTTTGTAGTCACGGACGAACGATAAACGAGCATAATCTTTCGAGTCGGAGGCGAGTACTTAACGGATATAACGTTTCGTGCCAATGTTAACCGGTCAACTACCACTCAGTTTCTTGTTCATCATAACACTGAAACTGAGATCGTCTTTGGTGCAATTCCAATACGGCTAACTTACGCATACTTTGATGACGCCGTGATTATATCAAGAACCTACCGCTTTCATGGCGGTAACGGTATCCAAAGAATTGGTGTGTTTCGTGCATGCAGTGTCGGACTAAGACTAGGAATGTTTGCAGTGGCCGCAGTTATTCCAAGAGGATGCTTCTTTCCAGCTAACGGTCGCCTAATAAGATGTAACTGGTTTCTTGAGGAGGCATGTACCCGAAGCTTAGAGTAGTCTCCTCTATGAATACTGAAGGACTTGCGTAGTTATGTACAAGCTCACCAACGGACGGGTGCTTCCACATATAACGTTAGCATCTCGTGTGCTATTCGTAAGAGTTTCTAGTCACGGACGAACGATAAAGTACCAACGCCTTTCATGAGTGGCCTCTCGAATCAAGTGATCGGACCTTTGGACGCACTCGATAACCGGGAAGTTATCCAGACTTTCGTGCCAATGTTAACCGGTCAATAAGAGCTACCTTTGATGACGCCGTGATTATATCAATACGCTTCTGGTTTGGGCGTCCAAGTTGCCGTCAAATAGTAGTGACCTTTGCAGTCTGAGTTCATGCGAGAATCACCGCGAAGTTTGTTGTTCATCATAACACTGAAATCCGCAATTAGTTTCCAGCTAACGGTCGCCTAATAATCGGTAGCACGTTTCCGAGCGCCACAGCCCAACTTATCTTACACACGTTTCATCGAAATCCGCTTGGTAACATGCAAGTGTAGTTTGCAGTGGCCGCAGTTATTCCAATGTGTGTGAGCTTTCGAGTTATCCGACACATCAAACAACCGATTAACTTTCGTGCATGCAGTGTCGGACTACAAGAATAGTGCTTTGATGGCGGTAACGGTATCCAACACACTATTACCTTTCATGTCGGGAGCCGCTTTGTACAGTGTGCTACCTTTGGTGCAATTCCAATACGGCTACATAGCTAACGGTTTCTTGAGGAGGCATGTACCCGACCGCCTTATTCGTTTGATGGAAGACGTTTCCGCTAACCGGTTCTAAGGTTTCGGACCGCCGTAAGCGAGTATCCTAGTACATGGTTTCGCTGCCGGCGGAATGCTATTCCTGCTAACGAGTTTCATGAGATGGATCGCATACTACGAAGCGTCATGTTTGGCCGGTACTCTCCAACCGTTCGAGCTTCTTCCTTTCGAGTCGGAGGCGAGTACTTACGATGATAGAGCTTTCAGACACGCGATCGGGACCACCGCAGCTAACAGTAATAGGACCGACCGACCGTTCGATTCAGGGAGATTGCCCTACACTATGCGGCAGCTGGCATAGACTCCTAAGGAGATGCGTACTTGTTAAATAGGACTCTTTCATCGAAATCCGCTTGGTAACCGCTAGGTTACGTTTGTTGTTCATCATAACACTGAACGTAACTATGTCTTTCGAGTTATCCGACACATCAAACTAAGTATGAGCTTTCCGAGCGCCACAGCCCAACTTCTAGCTAATCTCTTTGGCCGGTACTCTCCAACCGTTCTATTATGCCTGTTTGGCTGCCGGCGGAATGCTATTCTCCTGCTACACTTTGTAGTCACGGACGAACGATAACTGCTTAGAACCTTTGCAGTGGCCGCAGTTATTCCACTTAACGCGGAGTTTGGACGCACTCGATAACCGGGAGAACATTAGCTCTTTCATGAGATGGATCGCATACTAGAAGACATTAGGTTTCATGACGCCGTGATTATATCAGAAGTGTTACGGTTTCGGCGTCCAAGTTGCCGTCAAGAATTGATGAGCTTTCGTGCCAATGTTAACCGGTCAGACATTATAGCCTTTGGTGCAATTCCAATACGGCTAGACTCACACGACTTTGGAGTCGGAGGCGAGTACTTAGAGTTGTTGACGTTTCCAGCTAACGGTCGCCTAATAGATGATAGAACGTTTGCAGTCTGAGTTCATGCGAGAGCAATAAGCTACTTTCGGACCGCCGTAAGCGAGTATGCGATTAAGTAGTTTGATGTCGGGAGCCGCTTTGTAGGATACTCGACGTTTCAGACACGCGATCGGGACCACGGTAATATGACGTTTCTTGAGGAGGCATGTACCCGAGGTGTTGCATGGTTTCATGGAAGACGTTTCCGCTAAGTAACGGTATTCTTTGATGAGTGGCCTCTCGAATCAGTACCAGACTTGTTTCATGGCGGTAACGGTATCCAAGTAGCATAATCGTTTGGTGCATGCAGTGTCGGACTAGTATGCTATCGCAGGAGGATGGGGCAGGACAGGACGCGGCCACCCCAGGCCTTTCCAGAGCAAACCTGGAGAAGATTCACAATAGACAGGCCAAGAAACCCGGTGCTTCCTCCAGAGCCGTTTAAAGCTGATATGAGGAAATAAAGAGTGAACTGGAAGGATCCGATATCGCCACCGTGGCTGAATGAGACTGGTGTCGACCTGTGCCT").as_bytes().to_owned();
        let read = String::from("GTATTGCTCATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCAAATGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCATAGAAGAAGCTTAGCTCTTGTGCGAGCTACAGGAACGATGTTTGATTAGAGTAAGCAGAGGACAAGGGCTCGCGTGCAGCCGAAGTTTGGCCGGTACTCTCCATACAGTGTGCTACCTTTGGTGCAATTCCAATACGGCTACATAGCTAACGGTTTCTTGAGGAGGCATGTACCCGACCGCCTTATTCGTTTGATGGAAGACGTTTCCGCTAACCGGTTCTAAGGTTTCGGGACCGCCGTAAGCGATTGATGAGCTTTCGTGCCAATGTTAACCGGTCAGACATTATAGCCTTTGGTGCAATTCCAATACGGTAATATGACGTTTCTTGAGGAGGCATGTACCCGAGGTGTTGCATGGTTTCATGGAAGACGTTTCCGCTAAGTAACGGTATTCTTTGATGAGTGGCCTCTCGAATCAGTACCAGACTTGTTTCATGGCGGTAACGGTATCCAAGTAGCATAATCGTTTGGTGCATGCAGTGTCGGACTAGTATGCTATCGCAGGAGGATGGGGCAGGACAGGACGCGGCCACCCCAGGCCTTTCCAGAGCAAACCTGGAGAAGATTCACAATAGACAGGCCAAGAAACCCGGTGCTTCCTCCAGAGCCGTTTAAAGCTGATATGAGGAAATAAAGAGTGAACTGGAAGGATCCCATATCGACAATACGTAACTGAACGAAGTACACCAGTATT").as_bytes().to_owned();
        let reference = ReferenceManager::find_seeds(&refseq, 20);

        let hits = find_greedy_non_overlapping_segments(&read, &refseq, &reference);

        for hit in hits.alignment_segments {
            trace!("Overlapping seeds ref: {} search: {}, length: {}, endref: {}, endsearch: {}\n", hit.ref_start, hit.search_start, hit.length, hit.ref_start + hit.length, hit.search_start + hit.length);
        }
    }

    #[test]
    fn find_greedy_simple_ref_test() {
        let refseq = String::from("GTGGAAAGGACGAAACACCGGTACTTTCGAAAGTACGCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTCGCATTCTACCGTGACTTTAGCAAGGTGATCATTCGCAACAGTATCGACCTGCGTACTTTCGAAAGTACGCCGGTCGATGTTTGAATTCGAATTTAAATCGGATCCGCGGCCAA").as_bytes().to_owned();
        let read = String::from("GTGGAAAGGACGAAACACCGGTACTTTCGAAAGTACGCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTCGCATTCTACCGTGACTTTAGCAAGGTGATCATTCGCAACAGTATCGACGGGCGTACTTTCGAAAGTACGCCCGTCGATGTTTGAATTCGAATTTAAATCGGATCCGCGGCCAA").as_bytes().to_owned();
        let reference = ReferenceManager::find_seeds(&refseq, 20);

        let hits = find_greedy_non_overlapping_segments(&read, &refseq, &reference);

        for hit in hits.alignment_segments {
            println!("Overlapping seeds ref: {} search: {}, length: {}, endref: {}, endsearch: {}\n", hit.ref_start, hit.search_start, hit.length, hit.ref_start + hit.length, hit.search_start + hit.length);
        }
    }

    /*
    #[test]
    fn inversion_alignment_big() {
        let reference = String::from("NNNNNNNNCATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCATATGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCGTAGAAGAAAGTTAGCTCTTGTGCGAGCTACAGGAACGATGTTTGATTAGAGTAAGCAGAGGACAAGGGCTCGCGTGCAGCCGAAGTTTGGCCGGTACTCTCCAACCGTTAACAACAACACCTTTCATCGAAATCCGCTTGGTAACAACACTAGGATGTTTCGACGCACTCGATAACCGGGAAACCAAGAGAAGTTTCCGAGCGCCACAGCCCAACTTAACTTCGCCATGTTTGAGACACGCGATCGGGACCACAAGACGTTACTCTTTGGGACCGCCGTAAGCGAGTATAAGCAGATTGTGTTTCGGCGTCCAAGTTGCCGTCAAAAGCTTACTGAGTTTCGCTGCCGGCGGAATGCTATTAATCGCGCCTACTTTCATGGAAGACGTTTCCGCTAAAATGGACTGGTGTTTCATGTCGGGAGCCGCTTTGTAAGATGGAGCTACTTTCCAGTCTGAGTTCATGCGAGAACACCACGAGAGTTTCATGAGTGGCCTCTCGAATCAACAGTCTACAAGTTTGGAGTTATCCGACACATCAAAACCAGCCATTCGTTTCATGAGATGGATCGCATACTAACCTTAACGGAGTTTGTAGTCACGGACGAACGATAAACGAGCATAATCTTTCGAGTCGGAGGCGAGTACTTAACGGATATAACGTTTCGTGCCAATGTTAACCGGTCAACTACCACTCAGTTTCTTGTTCATCATAACACTGAAACTGAGATCGTCTTTGGTGCAATTCCAATACGGCTAACTTACGCATACTTTGATGACGCCGTGATTATATCAAGAACCTACCGCTTTCATGGCGGTAACGGTATCCAAAGAATTGGTGTGTTTCGTGCATGCAGTGTCGGACTAAGACTAGGAATGTTTGCAGTGGCCGCAGTTATTCCAAGAGGATGCTTCTTTCCAGCTAACGGTCGCCTAATAAGATGTAACTGGTTTCTTGAGGAGGCATGTACCCGAAGCTTAGAGTAGTCTCCTCTATGAATACTGAAGGACTTGCGTAGTTATGTACAAGCTCACCAACGGACGGGTGCTTCCACATATAACGTTAGCATCTCGTGTGCTATTCGTAAGAGTTTCTAGTCACGGACGAACGATAAAGTACCAACGCCTTTCATGAGTGGCCTCTCGAATCAAGTGATCGGACCTTTGGACGCACTCGATAACCGGGAAGTTATCCAGACTTTCGTGCCAATGTTAACCGGTCAATAAGAGCTACCTTTGATGACGCCGTGATTATATCAATACGCTTCTGGTTTGGGCGTCCAAGTTGCCGTCAAATAGTAGTGACCTTTGCAGTCTGAGTTCATGCGAGAATCACCGCGAAGTTTGTTGTTCATCATAACACTGAAATCCGCAATTAGTTTCCAGCTAACGGTCGCCTAATAATCGGTAGCACGTTTCCGAGCGCCACAGCCCAACTTATCTTACACACGTTTCATCGAAATCCGCTTGGTAACATGCAAGTGTAGTTTGCAGTGGCCGCAGTTATTCCAATGTGTGTGAGCTTTCGAGTTATCCGACACATCAAACAACCGATTAACTTTCGTGCATGCAGTGTCGGACTACAAGAATAGTGCTTTGATGGCGGTAACGGTATCCAACACACTATTACCTTTCATGTCGGGAGCCGCTTTGTACAGTGTGCTACCTTTGGTGCAATTCCAATACGGCTACATAGCTAACGGTTTCTTGAGGAGGCATGTACCCGACCGCCTTATTCGTTTGATGGAAGACGTTTCCGCTAACCGGTTCTAAGGTTTCGGACCGCCGTAAGCGAGTATCCTAGTACATGGTTTCGCTGCCGGCGGAATGCTATTCCTGCTAACGAGTTTCATGAGATGGATCGCATACTACGAAGCGTCATGTTTGGCCGGTACTCTCCAACCGTTCGAGCTTCTTCCTTTCGAGTCGGAGGCGAGTACTTACGATGATAGAGCTTTCAGACACGCGATCGGGACCACCGCAGCTAACAGTAATAGGACCGACCGACCGTTCGATTCAGGGAGATTGCCCTACACTATGCGGCAGCTGGCATAGACTCCTAAGGAGATGCGTACTTGTTAAATAGGACTCTTTCATCGAAATCCGCTTGGTAACCGCTAGGTTACGTTTGTTGTTCATCATAACACTGAACGTAACTATGTCTTTCGAGTTATCCGACACATCAAACTAAGTATGAGCTTTCCGAGCGCCACAGCCCAACTTCTAGCTAATCTCTTTGGCCGGTACTCTCCAACCGTTCTATTATGCCTGTTTGGCTGCCGGCGGAATGCTATTCTCCTGCTACACTTTGTAGTCACGGACGAACGATAACTGCTTAGAACCTTTGCAGTGGCCGCAGTTATTCCACTTAACGCGGAGTTTGGACGCACTCGATAACCGGGAGAACATTAGCTCTTTCATGAGATGGATCGCATACTAGAAGACATTAGGTTTCATGACGCCGTGATTATATCAGAAGTGTTACGGTTTCGGCGTCCAAGTTGCCGTCAAGAATTGATGAGCTTTCGTGCCAATGTTAACCGGTCAGACATTATAGCCTTTGGTGCAATTCCAATACGGCTAGACTCACACGACTTTGGAGTCGGAGGCGAGTACTTAGAGTTGTTGACGTTTCCAGCTAACGGTCGCCTAATAGATGATAGAACGTTTGCAGTCTGAGTTCATGCGAGAGCAATAAGCTACTTTCGGACCGCCGTAAGCGAGTATGCGATTAAGTAGTTTGATGTCGGGAGCCGCTTTGTAGGATACTCGACGTTTCAGACACGCGATCGGGACCACGGTAATATGACGTTTCTTGAGGAGGCATGTACCCGAGGTGTTGCATGGTTTCATGGAAGACGTTTCCGCTAAGTAACGGTATTCTTTGATGAGTGGCCTCTCGAATCAGTACCAGACTTGTTTCATGGCGGTAACGGTATCCAAGTAGCATAATCGTTTGGTGCATGCAGTGTCGGACTAGTATGCTATCGCAGGAGGATGGGGCAGGACAGGACGCGGCCACCCCAGGCCTTTCCAGAGCAAACCTGGAGAAGATTCACAATAGACAGGCCAAGAAACCCGGTGCTTCCTCCAGAGCCGTTTAAAGCTGATATGAGGAAATAAAGAGTGAACTGGAAGGATCCGATATCGCCACCGTGGCTGAATGAGACTGGTGTCGACCTGTGCCT").as_bytes().to_owned();
        let test_read = String::from("GTATTGCTCATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCAAATGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCATAGAAGAAGCTTAGCTCTTGTGCGAGCTACAGGAACGATGTTTGATTAGAGTAAGCAGAGGACAAGGGCTCGCGTGCAGCCGAAGTTTGGCCGGTACTCTCCATACAGTGTGCTACCTTTGGTGCAATTCCAATACGGCTACATAGCTAACGGTTTCTTGAGGAGGCATGTACCCGACCGCCTTATTCGTTTGATGGAAGACGTTTCCGCTAACCGGTTCTAAGGTTTCGGGACCGCCGTAAGCGATTGATGAGCTTTCGTGCCAATGTTAACCGGTCAGACATTATAGCCTTTGGTGCAATTCCAATACGGTAATATGACGTTTCTTGAGGAGGCATGTACCCGAGGTGTTGCATGGTTTCATGGAAGACGTTTCCGCTAAGTAACGGTATTCTTTGATGAGTGGCCTCTCGAATCAGTACCAGACTTGTTTCATGGCGGTAACGGTATCCAAGTAGCATAATCGTTTGGTGCATGCAGTGTCGGACTAGTATGCTATCGCAGGAGGATGGGGCAGGACAGGACGCGGCCACCCCAGGCCTTTCCAGAGCAAACCTGGAGAAGATTCACAATAGACAGGCCAAGAAACCCGGTGCTTCCTCCAGAGCCGTTTAAAGCTGATATGAGGAAATAAAGAGTGAACTGGAAGGATCCCATATCGACAATACGTAACTGAACGAAGTACACCAGTATT").as_bytes().to_owned();
        let ref_fasta = str_to_fasta_vec(str::from_utf8(reference.as_slice()).unwrap());
        let test_read_fasta = str_to_fasta_vec(str::from_utf8(test_read.as_slice()).unwrap());

        let my_score = InversionScoring {
            match_score: 10.0,
            mismatch_score: -11.0,
            //special_character_score: 9.0,
            gap_open: -15.0,
            gap_extend: -5.0,
            inversion_penalty: -2.0,
            min_inversion_length: 20
        };


        let my_aff_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -11.0,
            special_character_score: 8.0,
            gap_open: -15.0,
            gap_extend: -5.0,
            final_gap_multiplier: 1.0,

        };
        let reference_lookup = ReferenceManager::find_seeds(&reference, 20);
        let fwd_score_mp = find_greedy_non_overlapping_segments(&test_read, &reference, &reference_lookup);
        let mut alignment_mat: Alignment<Ix3> = create_scoring_record_3d((reference.len() + 1) * 2, (test_read.len() + 1) * 2, AlignmentType::Affine, false);

        let results = align_string_with_anchors(&test_read_fasta, &ref_fasta,  &fwd_score_mp, Some(&my_score), &my_aff_score,  &mut alignment_mat);

        trace!("CIGAR: {:?}",results.read_aligned);
    }

    #[test]
    fn inversion_alignment_big_inversion() {
        let reference = String::from("CATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCATATGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCGTAGAAGAAAGTTAGCTCTTGTGCGAGCTACAGGAACGATGTTTGATTAGAGTAAGCAGAGGACAAGGGCTCGCGTGCAGCCGAAGTTTGGCCGGTACTCTCCAACCGTTAACAACAACACCTTTCATCGAAATCCGCTTGGTAACAACACTAGGATGTTTCGACGCACTCGATAACCGGGAAACCAAGAGAAGTTTCCGAGCGCCACAGCCCAACTTAACTTCGCCATGTTTGAGACACGCGATCGGGACCACAAGACGTTACTCTTTGGGACCGCCGTAAGCGAGTATAAGCAGATTGTGTTTCGGCGTCCAAGTTGCCGTCAAAAGCTTACTGAGTTTCGCTGCCGGCGGAATGCTATTAATCGCGCCTACTTTCATGGAAGACGTTTCCGCTAAAATGGACTGGTGTTTCATGTCGGGAGCCGCTTTGTAAGATGGAGCTACTTTCCAGTCTGAGTTCATGCGAGAACACCACGAGAGTTTCATGAGTGGCCTCTCGAATCAACAGTCTACAAGTTTGGAGTTATCCGACACATCAAAACCAGCCATTCGTTTCATGAGATGGATCGCATACTAACCTTAACGGAGTTTGTAGTCACGGACGAACGATAAACGAGCATAATCTTTCGAGTCGGAGGCGAGTACTTAACGGATATAACGTTTCGTGCCAATGTTAACCGGTCAACTACCACTCAGTTTCTTGTTCATCATAACACTGAAACTGAGATCGTCTTTGGTGCAATTCCAATACGGCTAACTTACGCATACTTTGATGACGCCGTGATTATATCAAGAACCTACCGCTTTCATGGCGGTAACGGTATCCAAAGAATTGGTGTGTTTCGTGCATGCAGTGTCGGACTAAGACTAGGAATGTTTGCAGTGGCCGCAGTTATTCCAAGAGGATGCTTCTTTCCAGCTAACGGTCGCCTAATAAGATGTAACTGGTTTCTTGAGGAGGCATGTACCCGAAGCTTAGAGTAGTCTCCTCTATGAATACTGAAGGACTTGCGTAGTTATGTACAAGCTCACCAACGGACGGGTGCTTCCACATATAACGTTAGCATCTCGTGTGCTATTCGTAAGAGTTTCTAGTCACGGACGAACGATAAAGTACCAACGCCTTTCATGAGTGGCCTCTCGAATCAAGTGATCGGACCTTTGGACGCACTCGATAACCGGGAAGTTATCCAGACTTTCGTGCCAATGTTAACCGGTCAATAAGAGCTACCTTTGATGACGCCGTGATTATATCAATACGCTTCTGGTTTGGGCGTCCAAGTTGCCGTCAAATAGTAGTGACCTTTGCAGTCTGAGTTCATGCGAGAATCACCGCGAAGTTTGTTGTTCATCATAACACTGAAATCCGCAATTAGTTTCCAGCTAACGGTCGCCTAATAATCGGTAGCACGTTTCCGAGCGCCACAGCCCAACTTATCTTACACACGTTTCATCGAAATCCGCTTGGTAACATGCAAGTGTAGTTTGCAGTGGCCGCAGTTATTCCAATGTGTGTGAGCTTTCGAGTTATCCGACACATCAAACAACCGATTAACTTTCGTGCATGCAGTGTCGGACTACAAGAATAGTGCTTTGATGGCGGTAACGGTATCCAACACACTATTACCTTTCATGTCGGGAGCCGCTTTGTACAGTGTGCTACCTTTGGTGCAATTCCAATACGGCTACATAGCTAACGGTTTCTTGAGGAGGCATGTACCCGACCGCCTTATTCGTTTGATGGAAGACGTTTCCGCTAACCGGTTCTAAGGTTTCGGACCGCCGTAAGCGAGTATCCTAGTACATGGTTTCGCTGCCGGCGGAATGCTATTCCTGCTAACGAGTTTCATGAGATGGATCGCATACTACGAAGCGTCATGTTTGGCCGGTACTCTCCAACCGTTCGAGCTTCTTCCTTTCGAGTCGGAGGCGAGTACTTACGATGATAGAGCTTTCAGACACGCGATCGGGACCACCGCAGCTAACAGTAATAGGACCGACCGACCGTTCGATTCAGGGAGATTGCCCTACACTATGCGGCAGCTGGCATAGACTCCTAAGGAGATGCGTACTTGTTAAATAGGACTCTTTCATCGAAATCCGCTTGGTAACCGCTAGGTTACGTTTGTTGTTCATCATAACACTGAACGTAACTATGTCTTTCGAGTTATCCGACACATCAAACTAAGTATGAGCTTTCCGAGCGCCACAGCCCAACTTCTAGCTAATCTCTTTGGCCGGTACTCTCCAACCGTTCTATTATGCCTGTTTGGCTGCCGGCGGAATGCTATTCTCCTGCTACACTTTGTAGTCACGGACGAACGATAACTGCTTAGAACCTTTGCAGTGGCCGCAGTTATTCCACTTAACGCGGAGTTTGGACGCACTCGATAACCGGGAGAACATTAGCTCTTTCATGAGATGGATCGCATACTAGAAGACATTAGGTTTCATGACGCCGTGATTATATCAGAAGTGTTACGGTTTCGGCGTCCAAGTTGCCGTCAAGAATTGATGAGCTTTCGTGCCAATGTTAACCGGTCAGACATTATAGCCTTTGGTGCAATTCCAATACGGCTAGACTCACACGACTTTGGAGTCGGAGGCGAGTACTTAGAGTTGTTGACGTTTCCAGCTAACGGTCGCCTAATAGATGATAGAACGTTTGCAGTCTGAGTTCATGCGAGAGCAATAAGCTACTTTCGGACCGCCGTAAGCGAGTATGCGATTAAGTAGTTTGATGTCGGGAGCCGCTTTGTAGGATACTCGACGTTTCAGACACGCGATCGGGACCACGGTAATATGACGTTTCTTGAGGAGGCATGTACCCGAGGTGTTGCATGGTTTCATGGAAGACGTTTCCGCTAAGTAACGGTATTCTTTGATGAGTGGCCTCTCGAATCAGTACCAGACTTGTTTCATGGCGGTAACGGTATCCAAGTAGCATAATCGTTTGGTGCATGCAGTGTCGGACTAGTATGCTATCGCAGGAGGATGGGGCAGGACAGGACGCGGCCACCCCAGGCCTTTCCAGAGCAAACCTGGAGAAGATTCACAATAGACAGGCCAAGAAACCCGGTGCTTCCTCCAGAGCCGTTTAAAGCTGATATGAGGAAATAAAGAGTGAACTGGAAGGATCCGATATCGCCACCGTGGCTGAATGAGACTGGTGTCGACCTGTGCCT").as_bytes().to_owned();
        let test_read = String::from("CATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGACGTTCCTGTAGCTCGCACAAGAGCTAACTTTCTTCTACGTGGACCATTCGAAAGGATAGGGTACCTTTCTCTAGACATATGGTCTCATGTTTGATTAGAGTAAGCAGAGGACAAGGGCTCGCGTGCAGCCGAAGTTTGGCCGGTACTCTCCAACCGTTAACAACAACACCTTTCATCGAAATCCGCTTGGTAACAACACTAGGATGTTTCGACGCACTCGATAACCGGGAAACCAAGAGAAGTTTCCGAGCGCCACAGCCCAACTTAACTTCGCCATGTTTGAGACACGCGATCGGGACCACAAGACGTTACTCTTTGGGACCGCCGTAAGCGAGTATAAGCAGATTGTGTTTCGGCGTCCAAGTTGCCGTCAAAAGCTTACTGAGTTTCGCTGCCGGCGGAATGCTATTAATCGCGCCTACTTTCATGGAAGACGTTTCCGCTAAAATGGACTGGTGTTTCATGTCGGGAGCCGCTTTGTAAGATGGAGCTACTTTCCAGTCTGAGTTCATGCGAGAACACCACGAGAGTTTCATGAGTGGCCTCTCGAATCAACAGTCTACAAGTTTGGAGTTATCCGACACATCAAAACCAGCCATTCGTTTCATGAGATGGATCGCATACTAACCTTAACGGAGTTTGTAGTCACGGACGAACGATAAACGAGCATAATCTTTCGAGTCGGAGGCGAGTACTTAACGGATATAACGTTTCGTGCCAATGTTAACCGGTCAACTACCACTCAGTTTCTTGTTCATCATAACACTGAAACTGAGATCGTCTTTGGTGCAATTCCAATACGGCTAACTTACGCATACTTTGATGACGCCGTGATTATATCAAGAACCTACCGCTTTCATGGCGGTAACGGTATCCAAAGAATTGGTGTGTTTCGTGCATGCAGTGTCGGACTAAGACTAGGAATGTTTGCAGTGGCCGCAGTTATTCCAAGAGGATGCTTCTTTCCAGCTAACGGTCGCCTAATAAGATGTAACTGGTTTCTTGAGGAGGCATGTACCCGAAGCTTAGAGTAGTCTCCTCTATGAATACTGAAGGACTTGCGTAGTTATGTACAAGCTCACCAACGGACGGGTGCTTCCACATATAACGTTAGCATCTCGTGTGCTATTCGTAAGAGTTTCTAGTCACGGACGAACGATAAAGTACCAACGCCTTTCATGAGTGGCCTCTCGAATCAAGTGATCGGACCTTTGGACGCACTCGATAACCGGGAAGTTATCCAGACTTTCGTGCCAATGTTAACCGGTCAATAAGAGCTACCTTTGATGACGCCGTGATTATATCAATACGCTTCTGGTTTGGGCGTCCAAGTTGCCGTCAAATAGTAGTGACCTTTGCAGTCTGAGTTCATGCGAGAATCACCGCGAAGTTTGTTGTTCATCATAACACTGAAATCCGCAATTAGTTTCCAGCTAACGGTCGCCTAATAATCGGTAGCACGTTTCCGAGCGCCACAGCCCAACTTATCTTACACACGTTTCATCGAAATCCGCTTGGTAACATGCAAGTGTAGTTTGCAGTGGCCGCAGTTATTCCAATGTGTGTGAGCTTTCGAGTTATCCGACACATCAAACAACCGATTAACTTTCGTGCATGCAGTGTCGGACTACAAGAATAGTGCTTTGATGGCGGTAACGGTATCCAACACACTATTACCTTTCATGTCGGGAGCCGCTTTGTACAGTGTGCTACCTTTGGTGCAATTCCAATACGGCTACATAGCTAACGGTTTCTTGAGGAGGCATGTACCCGACCGCCTTATTCGTTTGATGGAAGACGTTTCCGCTAACCGGTTCTAAGGTTTCGGACCGCCGTAAGCGAGTATCCTAGTACATGGTTTCGCTGCCGGCGGAATGCTATTCCTGCTAACGAGTTTCATGAGATGGATCGCATACTACGAAGCGTCATGTTTGGCCGGTACTCTCCAACCGTTCGAGCTTCTTCCTTTCGAGTCGGAGGCGAGTACTTACGATGATAGAGCTTTCAGACACGCGATCGGGACCACCGCAGCTAACAGTAATAGGACCGACCGACCGTTCGATTCAGGGAGATTGCCCTACACTATGCGGCAGCTGGCATAGACTCCTAAGGAGATGCGTACTTGTTAAATAGGACTCTTTCATCGAAATCCGCTTGGTAACCGCTAGGTTACGTTTGTTGTTCATCATAACACTGAACGTAACTATGTCTTTCGAGTTATCCGACACATCAAACTAAGTATGAGCTTTCCGAGCGCCACAGCCCAACTTCTAGCTAATCTCTTTGGCCGGTACTCTCCAACCGTTCTATTATGCCTGTTTGGCTGCCGGCGGAATGCTATTCTCCTGCTACACTTTGTAGTCACGGACGAACGATAACTGCTTAGAACCTTTGCAGTGGCCGCAGTTATTCCACTTAACGCGGAGTTTGGACGCACTCGATAACCGGGAGAACATTAGCTCTTTCATGAGATGGATCGCATACTAGAAGACATTAGGTTTCATGACGCCGTGATTATATCAGAAGTGTTACGGTTTCGGCGTCCAAGTTGCCGTCAAGAATTGATGAGCTTTCGTGCCAATGTTAACCGGTCAGACATTATAGCCTTTGGTGCAATTCCAATACGGCTAGACTCACACGACTTTGGAGTCGGAGGCGAGTACTTAGAGTTGTTGACGTTTCCAGCTAACGGTCGCCTAATAGATGATAGAACGTTTGCAGTCTGAGTTCATGCGAGAGCAATAAGCTACTTTCGGACCGCCGTAAGCGAGTATGCGATTAAGTAGTTTGATGTCGGGAGCCGCTTTGTAGGATACTCGACGTTTCAGACACGCGATCGGGACCACGGTAATATGACGTTTCTTGAGGAGGCATGTACCCGAGGTGTTGCATGGTTTCATGGAAGACGTTTCCGCTAAGTAACGGTATTCTTTGATGAGTGGCCTCTCGAATCAGTACCAGACTTGTTTCATGGCGGTAACGGTATCCAAGTAGCATAATCGTTTGGTGCATGCAGTGTCGGACTAGTATGCTATCGCAGGAGGATGGGGCAGGACAGGACGCGGCCACCCCAGGCCTTTCCAGAGCAAACCTGGAGAAGATTCACAATAGACAGGCCAAGAAACCCGGTGCTTCCTCCAGAGCCGTTTAAAGCTGATATGAGGAAATAAAGAGTGAACTGGAAGGATCCGATATCGCCACCGTGGCTGAATGAGACTGGTGTCGACCTGTGCCT").as_bytes().to_owned();
        let ref_fasta = str_to_fasta_vec(str::from_utf8(reference.as_slice()).unwrap());
        let test_read_fasta = str_to_fasta_vec(str::from_utf8(test_read.as_slice()).unwrap());

        let my_score = InversionScoring {
            match_score: 10.0,
            mismatch_score: -11.0,
            //special_character_score: 9.0,
            gap_open: -15.0,
            gap_extend: -5.0,
            inversion_penalty: -2.0,
            min_inversion_length: 20
        };


        let my_aff_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -11.0,
            special_character_score: 8.0,
            gap_open: -15.0,
            gap_extend: -5.0,
            final_gap_multiplier: 1.0,

        };
        let reference_lookup = ReferenceManager::find_seeds(&reference, 20);
        let fwd_score_mp = find_greedy_non_overlapping_segments(&test_read, &reference, &reference_lookup);
        let mut alignment_mat: Alignment<Ix3> = create_scoring_record_3d((reference.len() + 1) * 2, (test_read.len() + 1) * 2, AlignmentType::Affine, false);

        let results = align_string_with_anchors(&test_read_fasta, &ref_fasta, &fwd_score_mp, Some(&my_score), &my_aff_score, &mut alignment_mat);

        println!("{} from {} reference {}", FastaBase::string(&results.read_aligned), String::from_utf8(test_read).unwrap(), String::from_utf8(reference).unwrap());
    }

    #[test]
    fn test_anchor_alignment() {


        let my_aff_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -11.0,
            special_character_score: 8.0,
            gap_open: -15.0,
            gap_extend: -5.0,
            final_gap_multiplier: 1.0,
        };

        let reference = String::from("CATGGTNNNNNNNNNNNNNNNNNNCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCATATGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCGTAGAAGAAAGTTAGCTCTTGTGCGA").as_bytes().to_owned();
        let ref_fasta = str_to_fasta_vec(str::from_utf8(reference.as_slice()).unwrap());
        let reference_lookup = ReferenceManager::find_seeds(&reference, 20);

        let test_read = String::from("CATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCAACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTGAAAGGTACCCTATCCTTTCGAATGGTCCACGCGTAGAAGAAAGTTAGTGCGA").as_bytes().to_owned();
        let test_read_fasta = str_to_fasta_vec(str::from_utf8(test_read.as_slice()).unwrap());

        let fwd_score_mp = find_greedy_non_overlapping_segments(&test_read, &reference, &reference_lookup);
        let mut alignment_mat: Alignment<Ix3> = create_scoring_record_3d((reference.len() + 1) * 2, (test_read.len() + 1) * 2, AlignmentType::Affine, false);

        let results = align_string_with_anchors(&test_read_fasta, &ref_fasta, &fwd_score_mp, None, &my_aff_score, &mut alignment_mat);

        println!("{} from {} reference {}", FastaBase::string(&results.read_aligned), String::from_utf8(test_read).unwrap(), String::from_utf8(reference.clone()).unwrap());
        assert_eq!(String::from("CATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCA------------ACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCT---------------------GAAAGGTACCCTATCCTTTCGAATGGTCCACGCGTAGAAGAAAGTTAG------TGCGA"),FastaBase::string(&results.read_aligned));


        // bigger dup
        let test_read = String::from("CATGGTAAAAAAAAAAAAAAAAAACGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCATAGTAACCGTTGCTAGGAGAGACCATATGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCGTAGAAGAAAGTTAGCTCTTGTGCGA").as_bytes().to_owned();
        let test_read_fasta = str_to_fasta_vec(str::from_utf8(test_read.as_slice()).unwrap());

        let fwd_score_mp = find_greedy_non_overlapping_segments(&test_read, &reference, &reference_lookup);
        let mut alignment_mat: Alignment<Ix3> = create_scoring_record_3d((reference.len() + 1) * 2, (test_read.len() + 1) * 2, AlignmentType::Affine, false);

        let results = align_string_with_anchors(&test_read_fasta, &ref_fasta, &fwd_score_mp, None, &my_aff_score, &mut alignment_mat);

        println!("{} from {} reference {}", FastaBase::string(&results.read_aligned), String::from_utf8(test_read).unwrap(), String::from_utf8(reference.clone()).unwrap());
        assert_eq!(String::from("CATGGTAAAAAAAAAAAAAAAAAACGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCATAGTAACCGTTGCTAGGAGAGACCATATGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCGTAGAAGAAAGTTAGCTCTTGTGCGA"),FastaBase::string(&results.read_aligned));
        assert_eq!(String::from("CATGGTNNNNNNNNNNNNNNNNNNCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCATA-------------------------TGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCGTAGAAGAAAGTTAGCTCTTGTGCGA"),FastaBase::string(&results.reference_aligned));

        // trailing gap
        let test_read = String::from("CATGGTAAAAAAAAAAAAAAAAAACGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCATAGTAACCGTTGCTAGGAGAGACCATATGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCGTAG").as_bytes().to_owned();
        let test_read_fasta = str_to_fasta_vec(str::from_utf8(test_read.as_slice()).unwrap());

        let fwd_score_mp = find_greedy_non_overlapping_segments(&test_read, &reference, &reference_lookup);
        let mut alignment_mat: Alignment<Ix3> = create_scoring_record_3d((reference.len() + 1) * 2, (test_read.len() + 1) * 2, AlignmentType::Affine, false);

        let results = align_string_with_anchors(&test_read_fasta, &ref_fasta, &fwd_score_mp, None, &my_aff_score, &mut alignment_mat);

        println!("{} from {} reference {}", FastaBase::string(&results.read_aligned), String::from_utf8(test_read).unwrap(), String::from_utf8(reference.clone()).unwrap());
        assert_eq!(String::from("CATGGTAAAAAAAAAAAAAAAAAACGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCATAGTAACCGTTGCTAGGAGAGACCATATGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCGTAG----------------------"),FastaBase::string(&results.read_aligned));
        assert_eq!(String::from("CATGGTNNNNNNNNNNNNNNNNNNCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCATA-------------------------TGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCGTAGAAGAAAGTTAGCTCTTGTGCGA"),FastaBase::string(&results.reference_aligned));


        // real world
        let reference = String::from("GTGAATTCGCCACCATGCTCGAGAGCTAGCGAATTCGAATTTNNNNNNNNNNNNNNNNAAATCGAGTAACCGTTGCTAGGAGAGACCCTAGCAGATCACCGTAAGGACTACCTTTCGTCCGGATAAGTTCATGTCGTTTCGCCCGTCGTAACCGGCGTAGTTTCCTACGAGTAATCTCTCGGGCTTTCGATCGAAGCAACCTGGCACCTTTCGTACTAATATACAGGATAGCTTTCCACCCAGGAGTAGCCCACGATTTCATAACCGCAGACGCTCTTATTTTCCCGGCACAGAGACTGTACGCTTTCTAAGCCACACCATCGAGTACTTTCGCACAGTAAGGCAATCGGATGTAATCCAGAGGTTGATTGTCGACCTTGCACTGTACTCTACGCGACTCTTTCACTGGCGCATAACTAGAGGGTTTCAGGGGCGAACAGTATATCGCTTTCAGCCCGTACGCATACGAGCGTTTCGCGTAGACAAGATAAGAGACTTTCATGCGAGATAGGTATACGCTTTTCCGCGGCCCGAACGCAGAAGATTTCCATACATACAAGGTTTGGGTTTTCGTGACGGACATGGGTGTGCCTTTCATCTGGTACGGATCCGAGTCTTTCACGACTCCCAGAGGCCGAGCTTTCGTACGAACAGAGGTCAATCTTTTCTGTGCGTACATCAGGACCAATTTCGTACCACAACGGAATATGACTTTCGTCCGCCACTATCCCTTCATTTTCGAAGTAACATACAGGTGGCGTTTCGCGGGAGCAGGAACCCGCTGTTTCATCCACTAATAGACACCGACTTTCGAGAGTGCACAGCTTAAGCTTTTCGCCTCACACTCAGCTCATTCTTTCGCGAGGCCAGACAAGCAGCGTTTCATCGGGGACGTACCCATGTTTTTCAGTCCCCTACATAATAGGGCTTTCGGCTAGCACGATCCCTAGGGTTTCGCGTGTTAGCTACCTATACCTTTCGTCCGGGATTATAGCCTCGGTTTCGTATGTTGTACAACCCGCGATTTCGCACCAGAGACGTCCAATCGTTTCACTATGACAGAAGTCGTCGCTTTCGGTAAAACAACCTCGGGATCTTTCGGAGCGGTAGAACGCGAGCGTTTCGCGATTCTAACTCTAATCCTTTTCTCGCAAAACAGTCCTAGACCTTTCTCGTCCCGACATCCATGGCTTTTCACCCCAGATCAGGTTGCTGCTTTCGCCATATACACCCGCGTGTCTTTCGATGTCGAACCTCCGCGAAGTTTCGCTGACACAACACTTAGAGCTTTCGCGAGGTAGTAGTCGGTCCTTTTCGTTTAAACACACCGGGTTAATTTCGGAGACGCATATTGTACCACGATTCTGTGGATAACCGTATTACCGCCTTTGAGTGAGCTGATACCGCTCGGCNNNNNNNNNNNNNNNNGCCACGTGGCTTTAAGGCCGGTCCTAGC").as_bytes().to_owned();
        let ref_fasta = str_to_fasta_vec(str::from_utf8(reference.as_slice()).unwrap());
        let reference_lookup = ReferenceManager::find_seeds(&reference, 20);

        let test_read = String::from("TTTTTGTTTTATATTACTTGGACCGTTGCGTATTGCTGCTAGGACCGGCCTTAAAGCCACCAGATAACCAAGCAAGGTGGCGGCCGCCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCGTCGTACGATGCATCTGCAAATTAACCCGGTGTGTTTAAACGAAAAGGACCGACTACTACCTCGCGAAAGCTCTAAGTGTTGTGTCAGCGAAACTTCGCGGAGGTTCGACATCGAAAGACACGCGGGTGTATATGGCGAAAGCAGCAACCTGATCTGGGGTGAAAAGCCATGGATGTCGGGACGAGAAAGGTCTAGGACTGTTTTGCGAGAAAAGGATTAGGTTAGAATCGCGAAACGCTCGCGTTCTACCACTCCGAAAGATCCCGAGGTTGTTTTACCGAAAGCGACGACTTCTGTCATAGTGAAACGATTGGACGTCTCTGGTGCGAAATCGCGGGTTGTACAACATACGAAACCGAGGCTATAATCCCGGACGAAAGGTATAGGTAGCTAACACGCGAAACCCTAGGGATCGTGCTTGTCGAAGGCCTCTATGTAGGGGACTGAAAAACATGGGTACGTCCCCGATGAAACGCTGCTTGTCTCACACGCGAAAGAATGAGCTGAGTGTGAGGCGAAAAGCTGAGGCTGTTAGCTCCTCGAAAGTCGGATCTATTAGTGGATGAAACAGCGGGTTCCTGCTCCCGCGAAACGCCACCTCCTGTTACTTCGAAAATGAAGGGATAGTGGCGGACGAAGAAAGTCATATTCCGTTGTGGTACGAAATTGGTCCTGATGTACGCACAGAAAAGATTGACCTCTGTTCGTACGAAAGCTCGGCCTCTGGGAGTCGTGAAAGACTGGATCCATACAGATGAAAGGCACACCCATGTCCGTCACGAAAACCCAAACCTTGTATGTATGGAAATCTTGCTCGTTCGGGCCGCGGAAAAGCGTATACCTATCTCGCATGAAAGTCTCTTATCTTGTCTACGCGAAACGCTCGTATGCGTACGGGCTGAAAGCGATATACTGTTCGCCCCTGAAACCCTCTAGTTATGCGCCAGTGAAAGAGTCGCGTAGAGTACAGTGCAAGGTCGACAATCAACCTCTGGATTACCCATCCGATTGCCTTACTGTGCGAAAGTACTCGATGGTGTGGCTTAGAAAGCGTACAGTCTCTGTGCGGGAAAATAAGAGCGTCTGCGGTTATGAAATCGTGGGCTACTCCTCGGTAGAAACTATGCTGTATATTAGTACGAAAGGTGCCAGGTTGCTTCGATCGAAAGCCCGAGAGATTCTTGCGTAGGAAACTACGCCGGTTACGACGGGCGAAACGACATGAACTTATCCGGACGAAAGGTAGTCCTTACGGTGATCTGCTAGGGTCTCTCCTAGCAACGGGGTTCTGGATTTTAAATTCGAATTCGCTAGCTCTCGAGCATGGTGGCAAATTCACAGCAATACGTTGC").as_bytes().to_owned();
        let test_read_fasta = str_to_fasta_vec(str::from_utf8(test_read.as_slice()).unwrap());

        let fwd_score_mp = find_greedy_non_overlapping_segments(&test_read, &reference, &reference_lookup);
        let mut alignment_mat: Alignment<Ix3> = create_scoring_record_3d((reference.len() + 1) * 2, (test_read.len() + 1) * 2, AlignmentType::Affine, false);

        let results = align_string_with_anchors(&test_read_fasta, &ref_fasta, &fwd_score_mp, None, &my_aff_score, &mut alignment_mat);

        println!("{} from {} reference {}", FastaBase::string(&results.read_aligned), String::from_utf8(test_read).unwrap(), FastaBase::string(&results.reference_aligned));
        assert_eq!(String::from("-TTTTTGTTTT----ATATTACTTG-GACCGTTGCGTATTGC-----TGCTAGGACCGGCCTTAAAGCCA-CCAGATAACCAAGCAAGGTGGC--GGCCGCCGAG-CGGTATCAGCTCACTCAAAGGCGGTAATACGGTT-ATCCACAGAA--TC--GTCGTACGATGCATCTGCAAAT--TAACCCGGTGT-GTTTAAACGAAAAGGACCGACTACTACCTCGCGAAAGCTCTAAGTGTTGTGTC--AGCGAAACTTCGCGGAGGTTCGACATCGAAAGACACGCGGGTGTATATGGCGAAAGC-AGCAACCTGATCTGGGGTGAAAAGCCATGGATGTCGGGA-CG-AGAAAGGTCTAGGACTGTTTTGCGAGAAAAG-GA-T-TAGG--TTAGAATCGCGAAACGCTCGCGTTCTACCACTCCGAAAG-------------ATCCCGAGGTTG-TT-------TT--ACCGAAAGCGA--CGACTTCTGTCA-TAGTG-A-AACGATTGGACGTCTCTGGTGCGAA------ATCGCGGGTTGTACAACATACGAAACCGAG-GCTATAATCCCG--GACGAAAGGTATAGGTAG-CTAACACGCGAAACCCTAGGGAT-CGTGCTTGTCGAAGGCCTCTATGTAGGGGA----CTGAAAAACA-----TGGG--TACGTCCCCG--ATGAAACGCTGCTTGTCTCACA-CGCGAAAGAATGAGCTGAGTGTGAGGCGAAAAGCTGAGGCTGTTAGC--TCCT-CGAA---A-GTCGGATCTATT-AGTGGATGAAACAGCGGG------TTCCTGCTCCCGCGAAA-------------CGCCAC-CTCCTGTTA-CTTCGAAAATG-AAGGGATA---GTGGCG----GACGAAGAAAGTCATATTCCGTTGTGGTA-CGA--AATTGGTCCTGATGTACGCACAGA----AAAGATTGACCTCTGTTCG--T-ACGAA---AGCTCGGCCTCT-GGGAGTCGTGA-AAGACTG-GATCCATACAGATGAAAGGCAC--ACCCATGTCCGTCACGAAAACCCAAACCTTGTATGTATGGAAATCTTGCTCGTTCGGGCCGCGGAAAAGCGTATACCTA---TCTCGCATGAAAGTCT--CT-TATCTT--GTCTACGCGA----AACGCTCGTATGCGTACGGGCTGAAA---GCGATATACTGTTCGCCCCTGAAAC--CCTCTAGTTATGCGCCAGTGAAAG-AGTCGCGTAGAGTACAGTGCAAGGTCGACAATCAACCTCTGGATTACCCATCCGATTGCCTTACTGTGCGAAAGTACTCGATGGTGTGGC--TTAGAAAGCGTACAGTCTCTGTGC-GGGAAAATAAGAGCGTCTGCGGTTATGAAATCGTGGGCTACTCCTCGGTAGA-AACTAT-GCTGTATATTAGT-ACGAAAGGTGCCAGGTTGC----TTCG--ATCGAAAGCCCGAGAGATTCTTGCGTAG--GAAACTACGCCG--GTTACGACGGGCGAAACGACAT-GAACTTATCCGGACGAAAGGTAGTCCTTACGGTGATCTGC---TA--------GGGT--CTCTCCTA---GCA---------ACGGGGTTC--TGGAT-----T-TTA---AATTCGAATTCGC---TAGCTCTCGAGCATGGTGGCAAATTCA---CA----GCAATA----CG----TTGC"),FastaBase::string(&results.read_aligned));


    }*/
}


