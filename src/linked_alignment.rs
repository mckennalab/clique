use std::{str, fmt, cmp};

use bio::alignment::AlignmentOperation;


use crate::extractor::*;
use std::convert::TryFrom;
use crate::alignment::alignment_matrix;
use crate::AlignmentType;
use crate::SuffixTableLookup;
use crate::*;

/// find a read's orientation using the greatest total number of matching bases
///
/// # Arguments
///
/// * `search_string` - a u8 Vec representing the search string
/// * `reference` - a u8 Vec representing the reference string
/// * `seeds` - a suffix array lookup object
pub fn orient_by_longest_segment(search_string: &Vec<u8>, reference: &Vec<u8>, seeds: &SuffixTableLookup) -> (bool, SharedSegments, SharedSegments) {
    let fwd_score_mp = find_greedy_non_overlapping_segments(search_string, reference, seeds);
    let fwd_score: usize = fwd_score_mp.alignment_cigar.clone().into_iter().map(|p| p.length).sum();

    let rev_score_mp =  find_greedy_non_overlapping_segments(&bio::alphabets::dna::revcomp(search_string), reference, seeds);
    let rev_score: usize = rev_score_mp.alignment_cigar.clone().into_iter().map(|p| p.length).sum();

    (fwd_score > rev_score,fwd_score_mp,rev_score_mp)
}


/// create an alignment string, with matching lengths, from a read, reference, and their CIGAR string
///
/// # Arguments
///
/// * `search_string` - a u8 Vec representing the search string
/// * `reference` - a u8 Vec representing the reference string
/// * `alignment` - contains the starting position and the alignment CIGAR strings
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
            AlignmentTags::MatchMismatch(size) => {
                read_align.push_str(str::from_utf8(&read[current_read_pos..(current_read_pos + size)]).unwrap());
                ref_align.push_str(str::from_utf8(&reference[current_ref_pos..(current_ref_pos + size)]).unwrap());
                current_read_pos += size;
                current_ref_pos += size;
            }
            AlignmentTags::Del(size) => {
                read_align.push_str(String::from_utf8(vec![b'-'; size]).unwrap().as_str());
                ref_align.push_str(str::from_utf8(&reference[current_ref_pos..(current_ref_pos + size)]).unwrap());
                current_ref_pos += size;
            }
            AlignmentTags::Ins(size) => {
                ref_align.push_str(String::from_utf8(vec![b'-'; size]).unwrap().as_str());
                read_align.push_str(str::from_utf8(&read[current_read_pos..(current_read_pos + size)]).unwrap());
                current_read_pos += size;
            }
            AlignmentTags::Inversion(_size) => {
                assert_eq!(0, 1);
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
pub fn find_greedy_non_overlapping_segments(search_string: &Vec<u8>, reference: &Vec<u8>, seeds: &SuffixTableLookup) -> SharedSegments {
    let mut return_hits: Vec<MatchedPosition> = Vec::new();
    let mut position = 0;
    let mut least_ref_pos = reference.len() as usize;
    let mut greatest_ref_pos = 0;

    while position < search_string.len() - seeds.seed_size {
        let ref_positions = seeds.suffix_table.positions(str::from_utf8(&search_string[position..(position + seeds.seed_size)]).unwrap());
        let mut longest_hit = 0;
        for ref_position in ref_positions {
            if ref_position >= &greatest_ref_pos {
                let extended_hit_size = extend_hit(search_string, position, reference, *ref_position as usize);
                if extended_hit_size > longest_hit {
                    return_hits.push(MatchedPosition { search_start: position, ref_start: *ref_position as usize, length: extended_hit_size });
                    position += extended_hit_size;
                    least_ref_pos = cmp::min(usize::try_from(*ref_position).unwrap(),least_ref_pos);
                    greatest_ref_pos = cmp::max(ref_position + &(extended_hit_size as u32),greatest_ref_pos);
                    longest_hit = extended_hit_size;
                }
            }
        }
        position += 1;
    }
    SharedSegments { start_position: least_ref_pos as usize, alignment_cigar: return_hits}
}

/// find a series of exact matches between the search string and the reference
///
/// # Arguments
///
/// * `search_string` - a u8 Vec representing the search string
/// * `reference` - a u8 Vec representing the reference string
/// * `seeds` - a suffix array lookup object
pub fn align_with_anchors(search_string: &Vec<u8>, reference: &Vec<u8>, min_alignment_seg_length: usize, overlaps: &SharedSegments) -> AlignmentCigar {
    let mut alignment_tags: Vec<AlignmentTags> = Vec::new();
    let mut read_alignment_last_position: usize = 0;
    let mut ref_alignment_last_position: usize = 0;

    for overlap in &overlaps.alignment_cigar {
        //println!("read_alignment_last_position : {},{}({}) ref_alignment_last_position : {},{}({}), length {}",read_alignment_last_position,overlap.search_start,search_string.len(),ref_alignment_last_position,overlap.ref_start,reference.len(),overlap.length);
        assert!(read_alignment_last_position <= overlap.search_start,"READ START FAILURE: {} and {}",read_alignment_last_position,overlap.search_start);
        assert!(ref_alignment_last_position <= overlap.ref_start,"REF START FAILURE: {} and {} from {}",ref_alignment_last_position,overlap.ref_start,overlap.length);

        // look back to see what segment we haven't aligned in the read
        let read_slice = slice_for_alignment(&search_string, read_alignment_last_position, overlap.search_start);
        let ref_slice = slice_for_alignment(&reference, ref_alignment_last_position, overlap.ref_start);
        //println!("sizes {} and {} ",read_slice.len(), ref_slice.len());

        let alignment = unaligned_segment_to_alignment(&read_slice, &ref_slice, min_alignment_seg_length);

        let read_ref_aligned_length = read_ref_alignment_lengths(&alignment);
        read_alignment_last_position += read_ref_aligned_length.0;
        ref_alignment_last_position += read_ref_aligned_length.1;
        //println!("22 read_alignment_last_position : {} ref_alignment_last_position : {} OVERLAP {}",read_alignment_last_position,ref_alignment_last_position,overlap.length);
        alignment_tags.extend(alignment);

        // now add the matching segment
        alignment_tags.push(AlignmentTags::MatchMismatch(overlap.length));
        read_alignment_last_position += overlap.length;
        ref_alignment_last_position += overlap.length;
    }

    if overlaps.alignment_cigar.len() > 0 {
        let read_stop = overlaps.alignment_cigar[overlaps.alignment_cigar.len() - 1].search_start + overlaps.alignment_cigar[overlaps.alignment_cigar.len() - 1].length;
        if read_stop < search_string.len() {

            // look back to see what segment we haven't aligned in the read
            let read_slice = slice_for_alignment(&search_string, read_alignment_last_position, search_string.len());
            let ref_slice = slice_for_alignment(&reference, ref_alignment_last_position, reference.len());
            let alignment = unaligned_segment_to_alignment(&read_slice, &ref_slice, min_alignment_seg_length);

            alignment_tags.extend(alignment);
        }
    } else {
        let alignment = unaligned_segment_to_alignment(&search_string, &reference, min_alignment_seg_length);

        alignment_tags.extend(alignment);
    }

    AlignmentCigar{alignment_start: 0, alignment_tags}
}

pub fn read_ref_alignment_lengths(alignment_tags: &Vec<AlignmentTags>) -> (usize, usize) {
    let mut read_len = 0;
    let mut ref_len = 0;
    for tag in alignment_tags {
        //println!("TAG {}",tag);
        match tag {
            AlignmentTags::MatchMismatch(s) => {
                read_len += s;
                ref_len += s;
            }
            AlignmentTags::Del(s) => {
                ref_len += s;
            }
            AlignmentTags::Ins(s) => {
                read_len += s;
            }
            AlignmentTags::Inversion(_size) => {
                assert_eq!(0, 1);
            }
        }
    }
    (read_len, ref_len)
}


pub fn slice_for_alignment(read: &Vec<u8>, read_start: usize, read_stop: usize) -> Vec<u8> {
    assert!(read_stop <= read.len(),"Read position requested {} when our length is only {}",read_stop,read.len());
    read[read_start..read_stop].to_vec()
}

pub fn unaligned_segment_to_alignment(read_segment: &Vec<u8>, reference_segment: &Vec<u8>, min_alignment_size: usize) -> Vec<AlignmentTags> {
    if read_segment.len() == reference_segment.len() && read_segment.len() <= min_alignment_size {
        vec!(AlignmentTags::MatchMismatch(read_segment.len()))
    } else if read_segment.len() < min_alignment_size && reference_segment.len() > 1 { // TODO: THIS IS VERY LAX AND WRONG BUT OK FOR TODAY
        let mut vec = Vec::new();
        vec.push(AlignmentTags::MatchMismatch(read_segment.len()));
        vec.push(AlignmentTags::Del(reference_segment.len()-read_segment.len()));
        vec
    } else {
        let alignment = align_forward_read( reference_segment, read_segment );
        convert_alignments(&alignment.0.operations)
    }
}


pub fn convert_alignments(bio: &Vec<AlignmentOperation>) -> Vec<AlignmentTags> {
    let mut new_tags = Vec::new();
    if bio.len() == 0 {
        return new_tags
    }
    let mut last_tag = AlignmentTags::MatchMismatch(0);
    for tag in bio {
        match tag {
            AlignmentOperation::Match | AlignmentOperation::Subst => {
                match last_tag {
                    AlignmentTags::MatchMismatch(size) => {
                        last_tag = AlignmentTags::MatchMismatch(size + 1);
                    }
                    _ => {
                        if last_tag != AlignmentTags::MatchMismatch(0) {new_tags.push(last_tag);}
                        last_tag = AlignmentTags::MatchMismatch(1);
                    }
                }
            }
            AlignmentOperation::Del => {
                match last_tag {
                    AlignmentTags::Del(size) => {
                        last_tag = AlignmentTags::Del(size + 1);
                    }
                    _ => {
                        if last_tag != AlignmentTags::MatchMismatch(0) {new_tags.push(last_tag);}
                        last_tag = AlignmentTags::Del(1);
                    }
                }
            }
            AlignmentOperation::Ins => {
                match last_tag {
                    AlignmentTags::Ins(size) => {
                        last_tag = AlignmentTags::Ins(size + 1);
                    }
                    _ => {
                        if last_tag != AlignmentTags::MatchMismatch(0) {new_tags.push(last_tag);}
                        last_tag = AlignmentTags::Ins(1);
                    }
                }
            }
            AlignmentOperation::Xclip(xclip_size) => {
                match last_tag {
                    AlignmentTags::Ins(size) => {
                        last_tag = AlignmentTags::Ins(size + xclip_size);
                    }
                    _ => {
                        if last_tag != AlignmentTags::MatchMismatch(0) {new_tags.push(last_tag);}
                        last_tag = AlignmentTags::Ins(*xclip_size);
                    }
                }
            }
            AlignmentOperation::Yclip(yclip_size) => {
                match last_tag {
                    AlignmentTags::Del(size) => {
                        last_tag = AlignmentTags::Del(size + yclip_size);
                    }
                    _ => {
                        if last_tag != AlignmentTags::MatchMismatch(0) {new_tags.push(last_tag);}
                        last_tag = AlignmentTags::Del(*yclip_size);
                    }
                }
            }
        }
    }
    new_tags.push(last_tag);
    new_tags
}


/// Extend a seed hit within the reference to its maximum length, using degenerate base matching
pub fn extend_hit(search_string: &Vec<u8>, search_location: usize, reference: &Vec<u8>, reference_location: usize) -> usize {
    let mut current_length = 0;
    while current_length + search_location < search_string.len() && current_length + reference_location < reference.len() {
        let search_loc = current_length + search_location;
        let ref_loc = current_length + reference_location;

        assert!(DEGENERATEBASES.contains_key(&search_string[search_loc]));
        assert!(DEGENERATEBASES.contains_key(&reference[ref_loc]));

        if DEGENERATEBASES.get(&search_string[search_loc]).unwrap().contains_key(&reference[ref_loc]) ||
            DEGENERATEBASES.get(&reference[ref_loc]).unwrap().contains_key(&search_string[search_loc]) {
            current_length += 1;
        } else {
            return current_length;
        }
    }
    current_length
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let reference = find_seeds(&refseq, 20);
        assert!(reference.suffix_table.contains("AAT"));
        assert!(!reference.suffix_table.contains("TAAT"));
    }

    #[test]
    fn find_greedy_non_overlapping_segments_test() {
        let refseq = String::from("NNNNNNNNCATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCATATGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCGTAGAAGAAAGTTAGCTCTTGTGCGAGCTACAGGAACGATGTTTGATTAGAGTAAGCAGAGGACAAGGGCTCGCGTGCAGCCGAAGTTTGGCCGGTACTCTCCAACCGTTAACAACAACACCTTTCATCGAAATCCGCTTGGTAACAACACTAGGATGTTTCGACGCACTCGATAACCGGGAAACCAAGAGAAGTTTCCGAGCGCCACAGCCCAACTTAACTTCGCCATGTTTGAGACACGCGATCGGGACCACAAGACGTTACTCTTTGGGACCGCCGTAAGCGAGTATAAGCAGATTGTGTTTCGGCGTCCAAGTTGCCGTCAAAAGCTTACTGAGTTTCGCTGCCGGCGGAATGCTATTAATCGCGCCTACTTTCATGGAAGACGTTTCCGCTAAAATGGACTGGTGTTTCATGTCGGGAGCCGCTTTGTAAGATGGAGCTACTTTCCAGTCTGAGTTCATGCGAGAACACCACGAGAGTTTCATGAGTGGCCTCTCGAATCAACAGTCTACAAGTTTGGAGTTATCCGACACATCAAAACCAGCCATTCGTTTCATGAGATGGATCGCATACTAACCTTAACGGAGTTTGTAGTCACGGACGAACGATAAACGAGCATAATCTTTCGAGTCGGAGGCGAGTACTTAACGGATATAACGTTTCGTGCCAATGTTAACCGGTCAACTACCACTCAGTTTCTTGTTCATCATAACACTGAAACTGAGATCGTCTTTGGTGCAATTCCAATACGGCTAACTTACGCATACTTTGATGACGCCGTGATTATATCAAGAACCTACCGCTTTCATGGCGGTAACGGTATCCAAAGAATTGGTGTGTTTCGTGCATGCAGTGTCGGACTAAGACTAGGAATGTTTGCAGTGGCCGCAGTTATTCCAAGAGGATGCTTCTTTCCAGCTAACGGTCGCCTAATAAGATGTAACTGGTTTCTTGAGGAGGCATGTACCCGAAGCTTAGAGTAGTCTCCTCTATGAATACTGAAGGACTTGCGTAGTTATGTACAAGCTCACCAACGGACGGGTGCTTCCACATATAACGTTAGCATCTCGTGTGCTATTCGTAAGAGTTTCTAGTCACGGACGAACGATAAAGTACCAACGCCTTTCATGAGTGGCCTCTCGAATCAAGTGATCGGACCTTTGGACGCACTCGATAACCGGGAAGTTATCCAGACTTTCGTGCCAATGTTAACCGGTCAATAAGAGCTACCTTTGATGACGCCGTGATTATATCAATACGCTTCTGGTTTGGGCGTCCAAGTTGCCGTCAAATAGTAGTGACCTTTGCAGTCTGAGTTCATGCGAGAATCACCGCGAAGTTTGTTGTTCATCATAACACTGAAATCCGCAATTAGTTTCCAGCTAACGGTCGCCTAATAATCGGTAGCACGTTTCCGAGCGCCACAGCCCAACTTATCTTACACACGTTTCATCGAAATCCGCTTGGTAACATGCAAGTGTAGTTTGCAGTGGCCGCAGTTATTCCAATGTGTGTGAGCTTTCGAGTTATCCGACACATCAAACAACCGATTAACTTTCGTGCATGCAGTGTCGGACTACAAGAATAGTGCTTTGATGGCGGTAACGGTATCCAACACACTATTACCTTTCATGTCGGGAGCCGCTTTGTACAGTGTGCTACCTTTGGTGCAATTCCAATACGGCTACATAGCTAACGGTTTCTTGAGGAGGCATGTACCCGACCGCCTTATTCGTTTGATGGAAGACGTTTCCGCTAACCGGTTCTAAGGTTTCGGACCGCCGTAAGCGAGTATCCTAGTACATGGTTTCGCTGCCGGCGGAATGCTATTCCTGCTAACGAGTTTCATGAGATGGATCGCATACTACGAAGCGTCATGTTTGGCCGGTACTCTCCAACCGTTCGAGCTTCTTCCTTTCGAGTCGGAGGCGAGTACTTACGATGATAGAGCTTTCAGACACGCGATCGGGACCACCGCAGCTAACAGTAATAGGACCGACCGACCGTTCGATTCAGGGAGATTGCCCTACACTATGCGGCAGCTGGCATAGACTCCTAAGGAGATGCGTACTTGTTAAATAGGACTCTTTCATCGAAATCCGCTTGGTAACCGCTAGGTTACGTTTGTTGTTCATCATAACACTGAACGTAACTATGTCTTTCGAGTTATCCGACACATCAAACTAAGTATGAGCTTTCCGAGCGCCACAGCCCAACTTCTAGCTAATCTCTTTGGCCGGTACTCTCCAACCGTTCTATTATGCCTGTTTGGCTGCCGGCGGAATGCTATTCTCCTGCTACACTTTGTAGTCACGGACGAACGATAACTGCTTAGAACCTTTGCAGTGGCCGCAGTTATTCCACTTAACGCGGAGTTTGGACGCACTCGATAACCGGGAGAACATTAGCTCTTTCATGAGATGGATCGCATACTAGAAGACATTAGGTTTCATGACGCCGTGATTATATCAGAAGTGTTACGGTTTCGGCGTCCAAGTTGCCGTCAAGAATTGATGAGCTTTCGTGCCAATGTTAACCGGTCAGACATTATAGCCTTTGGTGCAATTCCAATACGGCTAGACTCACACGACTTTGGAGTCGGAGGCGAGTACTTAGAGTTGTTGACGTTTCCAGCTAACGGTCGCCTAATAGATGATAGAACGTTTGCAGTCTGAGTTCATGCGAGAGCAATAAGCTACTTTCGGACCGCCGTAAGCGAGTATGCGATTAAGTAGTTTGATGTCGGGAGCCGCTTTGTAGGATACTCGACGTTTCAGACACGCGATCGGGACCACGGTAATATGACGTTTCTTGAGGAGGCATGTACCCGAGGTGTTGCATGGTTTCATGGAAGACGTTTCCGCTAAGTAACGGTATTCTTTGATGAGTGGCCTCTCGAATCAGTACCAGACTTGTTTCATGGCGGTAACGGTATCCAAGTAGCATAATCGTTTGGTGCATGCAGTGTCGGACTAGTATGCTATCGCAGGAGGATGGGGCAGGACAGGACGCGGCCACCCCAGGCCTTTCCAGAGCAAACCTGGAGAAGATTCACAATAGACAGGCCAAGAAACCCGGTGCTTCCTCCAGAGCCGTTTAAAGCTGATATGAGGAAATAAAGAGTGAACTGGAAGGATCCGATATCGCCACCGTGGCTGAATGAGACTGGTGTCGACCTGTGCCT").as_bytes().to_owned();
        let read = String::from("GTATTGCTCATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCAAATGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCATAGAAGAAGCTTAGCTCTTGTGCGAGCTACAGGAACGATGTTTGATTAGAGTAAGCAGAGGACAAGGGCTCGCGTGCAGCCGAAGTTTGGCCGGTACTCTCCATACAGTGTGCTACCTTTGGTGCAATTCCAATACGGCTACATAGCTAACGGTTTCTTGAGGAGGCATGTACCCGACCGCCTTATTCGTTTGATGGAAGACGTTTCCGCTAACCGGTTCTAAGGTTTCGGGACCGCCGTAAGCGATTGATGAGCTTTCGTGCCAATGTTAACCGGTCAGACATTATAGCCTTTGGTGCAATTCCAATACGGTAATATGACGTTTCTTGAGGAGGCATGTACCCGAGGTGTTGCATGGTTTCATGGAAGACGTTTCCGCTAAGTAACGGTATTCTTTGATGAGTGGCCTCTCGAATCAGTACCAGACTTGTTTCATGGCGGTAACGGTATCCAAGTAGCATAATCGTTTGGTGCATGCAGTGTCGGACTAGTATGCTATCGCAGGAGGATGGGGCAGGACAGGACGCGGCCACCCCAGGCCTTTCCAGAGCAAACCTGGAGAAGATTCACAATAGACAGGCCAAGAAACCCGGTGCTTCCTCCAGAGCCGTTTAAAGCTGATATGAGGAAATAAAGAGTGAACTGGAAGGATCCCATATCGACAATACGTAACTGAACGAAGTACACCAGTATT").as_bytes().to_owned();
        let reference = find_seeds(&refseq, 20);

        let hits = find_greedy_non_overlapping_segments(&read, &refseq, &reference);

        for hit in hits.alignment_cigar {
            println!("SEEEEEDS ref: {} search: {}, length: {}, endref: {}, endsearch: {}\n", hit.ref_start, hit.search_start, hit.length, hit.ref_start + hit.length, hit.search_start + hit.length);
        }
    }
    /* align_with_anchors(search_string: &Vec<u8>, reference: &Vec<u8>, seeds: &SuffixTableLookup, min_alignment_seg_length: usize)
    #[test]
    fn test_basic_align_with_anchors() {
        let refseq = String::from("NNNNNNNNCATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCATATGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCGTAGAAGAAAGTTAGCTCTTGTGCGAGCTACAGGAACGATGTTTGATTAGAGTAAGCAGAGGACAAGGGCTCGCGTGCAGCCGAAGTTTGGCCGGTACTCTCCAACCGTTAACAACAACACCTTTCATCGAAATCCGCTTGGTAACAACACTAGGATGTTTCGACGCACTCGATAACCGGGAAACCAAGAGAAGTTTCCGAGCGCCACAGCCCAACTTAACTTCGCCATGTTTGAGACACGCGATCGGGACCACAAGACGTTACTCTTTGGGACCGCCGTAAGCGAGTATAAGCAGATTGTGTTTCGGCGTCCAAGTTGCCGTCAAAAGCTTACTGAGTTTCGCTGCCGGCGGAATGCTATTAATCGCGCCTACTTTCATGGAAGACGTTTCCGCTAAAATGGACTGGTGTTTCATGTCGGGAGCCGCTTTGTAAGATGGAGCTACTTTCCAGTCTGAGTTCATGCGAGAACACCACGAGAGTTTCATGAGTGGCCTCTCGAATCAACAGTCTACAAGTTTGGAGTTATCCGACACATCAAAACCAGCCATTCGTTTCATGAGATGGATCGCATACTAACCTTAACGGAGTTTGTAGTCACGGACGAACGATAAACGAGCATAATCTTTCGAGTCGGAGGCGAGTACTTAACGGATATAACGTTTCGTGCCAATGTTAACCGGTCAACTACCACTCAGTTTCTTGTTCATCATAACACTGAAACTGAGATCGTCTTTGGTGCAATTCCAATACGGCTAACTTACGCATACTTTGATGACGCCGTGATTATATCAAGAACCTACCGCTTTCATGGCGGTAACGGTATCCAAAGAATTGGTGTGTTTCGTGCATGCAGTGTCGGACTAAGACTAGGAATGTTTGCAGTGGCCGCAGTTATTCCAAGAGGATGCTTCTTTCCAGCTAACGGTCGCCTAATAAGATGTAACTGGTTTCTTGAGGAGGCATGTACCCGAAGCTTAGAGTAGTCTCCTCTATGAATACTGAAGGACTTGCGTAGTTATGTACAAGCTCACCAACGGACGGGTGCTTCCACATATAACGTTAGCATCTCGTGTGCTATTCGTAAGAGTTTCTAGTCACGGACGAACGATAAAGTACCAACGCCTTTCATGAGTGGCCTCTCGAATCAAGTGATCGGACCTTTGGACGCACTCGATAACCGGGAAGTTATCCAGACTTTCGTGCCAATGTTAACCGGTCAATAAGAGCTACCTTTGATGACGCCGTGATTATATCAATACGCTTCTGGTTTGGGCGTCCAAGTTGCCGTCAAATAGTAGTGACCTTTGCAGTCTGAGTTCATGCGAGAATCACCGCGAAGTTTGTTGTTCATCATAACACTGAAATCCGCAATTAGTTTCCAGCTAACGGTCGCCTAATAATCGGTAGCACGTTTCCGAGCGCCACAGCCCAACTTATCTTACACACGTTTCATCGAAATCCGCTTGGTAACATGCAAGTGTAGTTTGCAGTGGCCGCAGTTATTCCAATGTGTGTGAGCTTTCGAGTTATCCGACACATCAAACAACCGATTAACTTTCGTGCATGCAGTGTCGGACTACAAGAATAGTGCTTTGATGGCGGTAACGGTATCCAACACACTATTACCTTTCATGTCGGGAGCCGCTTTGTACAGTGTGCTACCTTTGGTGCAATTCCAATACGGCTACATAGCTAACGGTTTCTTGAGGAGGCATGTACCCGACCGCCTTATTCGTTTGATGGAAGACGTTTCCGCTAACCGGTTCTAAGGTTTCGGACCGCCGTAAGCGAGTATCCTAGTACATGGTTTCGCTGCCGGCGGAATGCTATTCCTGCTAACGAGTTTCATGAGATGGATCGCATACTACGAAGCGTCATGTTTGGCCGGTACTCTCCAACCGTTCGAGCTTCTTCCTTTCGAGTCGGAGGCGAGTACTTACGATGATAGAGCTTTCAGACACGCGATCGGGACCACCGCAGCTAACAGTAATAGGACCGACCGACCGTTCGATTCAGGGAGATTGCCCTACACTATGCGGCAGCTGGCATAGACTCCTAAGGAGATGCGTACTTGTTAAATAGGACTCTTTCATCGAAATCCGCTTGGTAACCGCTAGGTTACGTTTGTTGTTCATCATAACACTGAACGTAACTATGTCTTTCGAGTTATCCGACACATCAAACTAAGTATGAGCTTTCCGAGCGCCACAGCCCAACTTCTAGCTAATCTCTTTGGCCGGTACTCTCCAACCGTTCTATTATGCCTGTTTGGCTGCCGGCGGAATGCTATTCTCCTGCTACACTTTGTAGTCACGGACGAACGATAACTGCTTAGAACCTTTGCAGTGGCCGCAGTTATTCCACTTAACGCGGAGTTTGGACGCACTCGATAACCGGGAGAACATTAGCTCTTTCATGAGATGGATCGCATACTAGAAGACATTAGGTTTCATGACGCCGTGATTATATCAGAAGTGTTACGGTTTCGGCGTCCAAGTTGCCGTCAAGAATTGATGAGCTTTCGTGCCAATGTTAACCGGTCAGACATTATAGCCTTTGGTGCAATTCCAATACGGCTAGACTCACACGACTTTGGAGTCGGAGGCGAGTACTTAGAGTTGTTGACGTTTCCAGCTAACGGTCGCCTAATAGATGATAGAACGTTTGCAGTCTGAGTTCATGCGAGAGCAATAAGCTACTTTCGGACCGCCGTAAGCGAGTATGCGATTAAGTAGTTTGATGTCGGGAGCCGCTTTGTAGGATACTCGACGTTTCAGACACGCGATCGGGACCACGGTAATATGACGTTTCTTGAGGAGGCATGTACCCGAGGTGTTGCATGGTTTCATGGAAGACGTTTCCGCTAAGTAACGGTATTCTTTGATGAGTGGCCTCTCGAATCAGTACCAGACTTGTTTCATGGCGGTAACGGTATCCAAGTAGCATAATCGTTTGGTGCATGCAGTGTCGGACTAGTATGCTATCGCAGGAGGATGGGGCAGGACAGGACGCGGCCACCCCAGGCCTTTCCAGAGCAAACCTGGAGAAGATTCACAATAGACAGGCCAAGAAACCCGGTGCTTCCTCCAGAGCCGTTTAAAGCTGATATGAGGAAATAAAGAGTGAACTGGAAGGATCCGATATCGCCACCGTGGCTGAATGAGACTGGTGTCGACCTGTGCCT").as_bytes().to_owned();
        let read = String::from("GTATTGCTCATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCAAATGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCATAGAAGAAGCTTAGCTCTTGTGCGAGCTACAGGAACGATGTTTGATTAGAGTAAGCAGAGGACAAGGGCTCGCGTGCAGCCGAAGTTTGGCCGGTACTCTCCATACAGTGTGCTACCTTTGGTGCAATTCCAATACGGCTACATAGCTAACGGTTTCTTGAGGAGGCATGTACCCGACCGCCTTATTCGTTTGATGGAAGACGTTTCCGCTAACCGGTTCTAAGGTTTCGGGACCGCCGTAAGCGATTGATGAGCTTTCGTGCCAATGTTAACCGGTCAGACATTATAGCCTTTGGTGCAATTCCAATACGGTAATATGACGTTTCTTGAGGAGGCATGTACCCGAGGTGTTGCATGGTTTCATGGAAGACGTTTCCGCTAAGTAACGGTATTCTTTGATGAGTGGCCTCTCGAATCAGTACCAGACTTGTTTCATGGCGGTAACGGTATCCAAGTAGCATAATCGTTTGGTGCATGCAGTGTCGGACTAGTATGCTATCGCAGGAGGATGGGGCAGGACAGGACGCGGCCACCCCAGGCCTTTCCAGAGCAAACCTGGAGAAGATTCACAATAGACAGGCCAAGAAACCCGGTGCTTCCTCCAGAGCCGTTTAAAGCTGATATGAGGAAATAAAGAGTGAACTGGAAGGATCCCATATCGACAATACGTAACTGAACGAAGTACACCAGTATT").as_bytes().to_owned();
        let reference = find_seeds(&refseq, 20);

        let fwd_score_mp = find_greedy_non_overlapping_segments(&read, &refseq, &reference);

        let hits = align_with_anchors(&read, &refseq,  10, &fwd_score_mp);

        for hit in hits.alignment_tags {
            print!("{}",hit);
        }
    }

    // align_with_anchors(search_string: &Vec<u8>, reference: &Vec<u8>, seeds: &SuffixTableLookup, min_alignment_seg_length: usize)
    #[test]
    fn test_to_string() {
        let refseq = String::from("NNNNNNNNCATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCATATGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCGTAGAAGAAAGTTAGCTCTTGTGCGAGCTACAGGAACGATGTTTGATTAGAGTAAGCAGAGGACAAGGGCTCGCGTGCAGCCGAAGTTTGGCCGGTACTCTCCAACCGTTAACAACAACACCTTTCATCGAAATCCGCTTGGTAACAACACTAGGATGTTTCGACGCACTCGATAACCGGGAAACCAAGAGAAGTTTCCGAGCGCCACAGCCCAACTTAACTTCGCCATGTTTGAGACACGCGATCGGGACCACAAGACGTTACTCTTTGGGACCGCCGTAAGCGAGTATAAGCAGATTGTGTTTCGGCGTCCAAGTTGCCGTCAAAAGCTTACTGAGTTTCGCTGCCGGCGGAATGCTATTAATCGCGCCTACTTTCATGGAAGACGTTTCCGCTAAAATGGACTGGTGTTTCATGTCGGGAGCCGCTTTGTAAGATGGAGCTACTTTCCAGTCTGAGTTCATGCGAGAACACCACGAGAGTTTCATGAGTGGCCTCTCGAATCAACAGTCTACAAGTTTGGAGTTATCCGACACATCAAAACCAGCCATTCGTTTCATGAGATGGATCGCATACTAACCTTAACGGAGTTTGTAGTCACGGACGAACGATAAACGAGCATAATCTTTCGAGTCGGAGGCGAGTACTTAACGGATATAACGTTTCGTGCCAATGTTAACCGGTCAACTACCACTCAGTTTCTTGTTCATCATAACACTGAAACTGAGATCGTCTTTGGTGCAATTCCAATACGGCTAACTTACGCATACTTTGATGACGCCGTGATTATATCAAGAACCTACCGCTTTCATGGCGGTAACGGTATCCAAAGAATTGGTGTGTTTCGTGCATGCAGTGTCGGACTAAGACTAGGAATGTTTGCAGTGGCCGCAGTTATTCCAAGAGGATGCTTCTTTCCAGCTAACGGTCGCCTAATAAGATGTAACTGGTTTCTTGAGGAGGCATGTACCCGAAGCTTAGAGTAGTCTCCTCTATGAATACTGAAGGACTTGCGTAGTTATGTACAAGCTCACCAACGGACGGGTGCTTCCACATATAACGTTAGCATCTCGTGTGCTATTCGTAAGAGTTTCTAGTCACGGACGAACGATAAAGTACCAACGCCTTTCATGAGTGGCCTCTCGAATCAAGTGATCGGACCTTTGGACGCACTCGATAACCGGGAAGTTATCCAGACTTTCGTGCCAATGTTAACCGGTCAATAAGAGCTACCTTTGATGACGCCGTGATTATATCAATACGCTTCTGGTTTGGGCGTCCAAGTTGCCGTCAAATAGTAGTGACCTTTGCAGTCTGAGTTCATGCGAGAATCACCGCGAAGTTTGTTGTTCATCATAACACTGAAATCCGCAATTAGTTTCCAGCTAACGGTCGCCTAATAATCGGTAGCACGTTTCCGAGCGCCACAGCCCAACTTATCTTACACACGTTTCATCGAAATCCGCTTGGTAACATGCAAGTGTAGTTTGCAGTGGCCGCAGTTATTCCAATGTGTGTGAGCTTTCGAGTTATCCGACACATCAAACAACCGATTAACTTTCGTGCATGCAGTGTCGGACTACAAGAATAGTGCTTTGATGGCGGTAACGGTATCCAACACACTATTACCTTTCATGTCGGGAGCCGCTTTGTACAGTGTGCTACCTTTGGTGCAATTCCAATACGGCTACATAGCTAACGGTTTCTTGAGGAGGCATGTACCCGACCGCCTTATTCGTTTGATGGAAGACGTTTCCGCTAACCGGTTCTAAGGTTTCGGACCGCCGTAAGCGAGTATCCTAGTACATGGTTTCGCTGCCGGCGGAATGCTATTCCTGCTAACGAGTTTCATGAGATGGATCGCATACTACGAAGCGTCATGTTTGGCCGGTACTCTCCAACCGTTCGAGCTTCTTCCTTTCGAGTCGGAGGCGAGTACTTACGATGATAGAGCTTTCAGACACGCGATCGGGACCACCGCAGCTAACAGTAATAGGACCGACCGACCGTTCGATTCAGGGAGATTGCCCTACACTATGCGGCAGCTGGCATAGACTCCTAAGGAGATGCGTACTTGTTAAATAGGACTCTTTCATCGAAATCCGCTTGGTAACCGCTAGGTTACGTTTGTTGTTCATCATAACACTGAACGTAACTATGTCTTTCGAGTTATCCGACACATCAAACTAAGTATGAGCTTTCCGAGCGCCACAGCCCAACTTCTAGCTAATCTCTTTGGCCGGTACTCTCCAACCGTTCTATTATGCCTGTTTGGCTGCCGGCGGAATGCTATTCTCCTGCTACACTTTGTAGTCACGGACGAACGATAACTGCTTAGAACCTTTGCAGTGGCCGCAGTTATTCCACTTAACGCGGAGTTTGGACGCACTCGATAACCGGGAGAACATTAGCTCTTTCATGAGATGGATCGCATACTAGAAGACATTAGGTTTCATGACGCCGTGATTATATCAGAAGTGTTACGGTTTCGGCGTCCAAGTTGCCGTCAAGAATTGATGAGCTTTCGTGCCAATGTTAACCGGTCAGACATTATAGCCTTTGGTGCAATTCCAATACGGCTAGACTCACACGACTTTGGAGTCGGAGGCGAGTACTTAGAGTTGTTGACGTTTCCAGCTAACGGTCGCCTAATAGATGATAGAACGTTTGCAGTCTGAGTTCATGCGAGAGCAATAAGCTACTTTCGGACCGCCGTAAGCGAGTATGCGATTAAGTAGTTTGATGTCGGGAGCCGCTTTGTAGGATACTCGACGTTTCAGACACGCGATCGGGACCACGGTAATATGACGTTTCTTGAGGAGGCATGTACCCGAGGTGTTGCATGGTTTCATGGAAGACGTTTCCGCTAAGTAACGGTATTCTTTGATGAGTGGCCTCTCGAATCAGTACCAGACTTGTTTCATGGCGGTAACGGTATCCAAGTAGCATAATCGTTTGGTGCATGCAGTGTCGGACTAGTATGCTATCGCAGGAGGATGGGGCAGGACAGGACGCGGCCACCCCAGGCCTTTCCAGAGCAAACCTGGAGAAGATTCACAATAGACAGGCCAAGAAACCCGGTGCTTCCTCCAGAGCCGTTTAAAGCTGATATGAGGAAATAAAGAGTGAACTGGAAGGATCCGATATCGCCACCGTGGCTGAATGAGACTGGTGTCGACCTGTGCCT").as_bytes().to_owned();
        let read = String::from("TCTATAACCTACTGGTTGGTTATGTATTCATCTACTTCGTTCAGTTACGTATTGCTCATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCATATGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCGTAGAAGAAAGTTAGCTCTTGTGCGAGCTACAGGAACGATGTTTGATTAGAGTAAGCGGACAAGGGCTCGCGTGCAGCCGAAGTTTGGCCGGTACTCTCCAACAACACCTTTATCGAAATCCGCTTGGTAACAACACTAGGATGTTTCGACGCACTCGATAACCGGGAAACCAAGAGAAGTTTCCGAGCGCCACAGCCCAACTTAACTTCGCCATGTTTGACGCGCAATAGGACCACAACGTTACTACTGGGACCGCCGTAAGCGAGTATCAGATTGTGTTTCGGCGTCCAAGTTGCCGTCAAAAGCTTACTGAGTTTCGCTGCCGGCGGAATGCTATTAATCGCGCCTACTTTCATGGAAGACGTTTCCGCTAAAATGGACTGGTGTTTCATGTCGGGAGCCGCTTTGTAAGATGGAGCTACTTTCCAGTCTGAGTTCATGCGAGAACACCACGAGAGTTTCATGAGTGGCCTCTCGAATCAACAGTCTACAAGTTTGGAGTTATCCGACACATCAAAACCAGCCATTCGTTTCATGAGATGGATCGCATACTAACCTTAACGGAGTTTGTAGTCACGGACGAACAATAACGACCATAATCTTTCGAGTCGGAGGCGAGTACTTAACGGATATAACGTTTCGTGCCAATGTTAACCGGTCAACTACCACTCAGTTTCTTGTTCATCATAACACTGAAACTGAGATCGTCTTTGGTGCAATTCCAATACGGCTAACTTACGCATACTTTGATGACGCCGTGATTATATCAAGAACCTACCGCTTTCATGGCGGTAACGGTATCCAAAGAATTGGTGTGTTTCGTGCATGCAGTGTCGGACTAAGACTAGGAATGTTTGCAGTGGCCGCAGTTATTCCAAGAGGATGCTTCTTTCCAGCTAACGGTCGCCTAATAAGATGTAACTGGTTTTTATTGAGGAGGCATGTACCCGAAGCTTAGAGTAGTCTCCTCTATGAATACTCAAGGACTTGCGTAGTTATGTACAAGCTCACCAACGGACGGGTGCTTCCACATATAACGTTAGCATCTCGTGTGCTATTCGTAAGAGTTTCTAGTCACGGACGAACGATAAAGTACCAACGCCTTTCATGAGTGGCCTCTCGAATCAAGTGATCGGACCTTTGGACGCACTACGGTTGGAGAGTACCGGCCAAACATGACGCTTCGTAGTATGCGATCCATCTCATGAAACTCGTTAGCGGATCGCATTCCGCCGGCAGCGAAATTATAGGATACTCGCTTACGGCGGTCCGAAACCTTAGAACCGGTTCCGGAAACGTCTTCCATCAAACGAATAAGGCGGTCGGGTACATGCCTCCTCAAGAAACCGTTAGCTATGTAGCCGTATTGGAATTGCACCAAAGGTAGCACACTGTACAAAGCGGCTCCCGACATGAAAGGTAATAGTGTGTTGGATACCGTTACCGCCATCAAAGCACTATTCTTGTAGTCCGACACTGCATGCACGAAAGTTAATCGGTTGTTTGATGTGTCGGATAACTCGAAAGCTCACACACATTGGAATAACTGCGGCCACTGCAAACTACAAGCGGATTTCGATGAAACGTGTGTATGAACTCAGACTGCAAACGTTCTATCATCTATTAGGCGACCGTTAGCTGGAAACGTCAACAACTCTAAGTACTCGCCTCCGACTCCAAAGTCGTGTGAGTCTAGCCGTATTGGAATTGCACCAAAGGCTATAATGTCTGACCGGTTAACAGCAATACGT").as_bytes().to_owned();
        let reference = find_seeds(&refseq, 20);

        let fwd_score_mp = find_greedy_non_overlapping_segments(&read, &refseq, &reference);

        let hits = align_with_anchors(&read, &refseq,  10, &fwd_score_mp);

        let alignment_strings = cigar_alignment_to_full_string(&read, &refseq, &hits);

        print!("ref:  {}\nread: {}\n",alignment_strings.0,alignment_strings.1);
    }*/


    //unaligned_segment_to_alignment(read_segment: &Vec<u8>, reference_segment: &Vec<u8>, min_alignment_size: usize)
    #[test]
    fn test_basic_unaligned_segment() {
        let refseq = String::from("AATGATACGG").as_bytes().to_owned();
        let read = String::from("AATGATACGG").as_bytes().to_owned();

        let hits = unaligned_segment_to_alignment(&read, &refseq, 10);

        assert_eq!(hits.len(),1);
        assert_eq!(hits[0],AlignmentTags::MatchMismatch(10));
    }

    //unaligned_segment_to_alignment(read_segment: &Vec<u8>, reference_segment: &Vec<u8>, min_alignment_size: usize)
    #[test]
    fn test_basic_indel_segment() {
        let refseq = String::from("AATGATACGGTTTTT").as_bytes().to_owned();
        let read = String::from("AATGAGGTTTTT").as_bytes().to_owned();

        let hits = unaligned_segment_to_alignment(&read, &refseq, 10);

        assert_eq!(hits.len(),3);
        assert_eq!(hits[0],AlignmentTags::MatchMismatch(5));
        assert_eq!(hits[1],AlignmentTags::Del(3));
        assert_eq!(hits[2],AlignmentTags::MatchMismatch(7));
    }
}


