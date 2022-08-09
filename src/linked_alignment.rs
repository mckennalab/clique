use std::{str, fmt};

use bio::alignment::{Alignment, AlignmentOperation};

use suffix::SuffixTable;

use crate::extractor::*;

pub struct SuffixTableLookup<'s, 't> {
    suffixTable: SuffixTable<'s, 't>,
    seed_size: usize,
}

#[derive(Copy, Clone, Debug)]
pub struct MatchedPosition {
    search_start: usize,
    ref_start: usize,
    length: usize,
}

pub struct MatchedPositions {
    positions: Vec<MatchedPosition>
}

#[derive(Eq, PartialEq, Debug, Copy, Clone, Hash)]
pub enum AlignmentTags {
    MatchMismatch(usize),
    Del(usize),
    Ins(usize),
}

impl fmt::Display for AlignmentTags {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            AlignmentTags::MatchMismatch(size) => write!(f, "{}M",size),
            AlignmentTags::Del(size) => write!(f, "{}D",size),
            AlignmentTags::Ins(size) => write!(f, "{}I",size),
        }
    }
}


pub fn find_seeds(reference: &Vec<u8>, seed_size: usize) -> SuffixTableLookup {
    return SuffixTableLookup { suffixTable: SuffixTable::new(String::from_utf8(reference.clone()).unwrap()), seed_size };
}


pub fn is_forward_orientation(search_string: &Vec<u8>, reference: &Vec<u8>, seeds: &SuffixTableLookup) -> (bool,MatchedPositions,MatchedPositions) {
    let fwd_score_mp = find_greedy_non_overlapping_segments(search_string, reference, seeds);
    let fwd_score: usize = fwd_score_mp.positions.clone().into_iter().map(|p| p.length).sum();

    let rev_score_mp =  find_greedy_non_overlapping_segments(&bio::alphabets::dna::revcomp(search_string), reference, seeds);
    let rev_score: usize = rev_score_mp.positions.clone().into_iter().map(|p| p.length).sum();

    (fwd_score > rev_score,fwd_score_mp,rev_score_mp)
}

pub fn find_greedy_non_overlapping_segments(search_string: &Vec<u8>, reference: &Vec<u8>, seeds: &SuffixTableLookup) -> MatchedPositions {
    let mut return_hits: Vec<MatchedPosition> = Vec::new();
    let mut position = 0;
    let mut highest_ref_pos = 0;

    while position < search_string.len() - seeds.seed_size {
        let ref_positions = seeds.suffixTable.positions(str::from_utf8(&search_string[position..(position + seeds.seed_size)]).unwrap());
        let mut longest_hit = 0;
        for ref_position in ref_positions {
            if ref_position >= &highest_ref_pos {
                let extended_hit_size = extend_hit(search_string, position, reference, *ref_position as usize);
                if extended_hit_size > longest_hit {
                    return_hits.push(MatchedPosition { search_start: position, ref_start: *ref_position as usize, length: extended_hit_size });
                    //println!("adding {},{},{}",position, *ref_position as usize, extended_hit_size);
                    position += extended_hit_size;
                    highest_ref_pos = ref_position + &(extended_hit_size as u32);
                }
            }
        }
        position += 1;
    }
    MatchedPositions { positions: return_hits }
}

pub fn align_with_anchors(search_string: &Vec<u8>, reference: &Vec<u8>, seeds: &SuffixTableLookup, min_alignment_seg_length: usize, overlaps: &MatchedPositions) -> Vec<AlignmentTags> {
    let mut alignmentTags: Vec<AlignmentTags> = Vec::new();
    let mut read_alignment_last_position: usize = 0;
    let mut ref_alignment_last_position: usize = 0;

    for overlap in &overlaps.positions {
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
        alignmentTags.extend(alignment);

        // now add the matching segment
        alignmentTags.push(AlignmentTags::MatchMismatch(overlap.length));
        read_alignment_last_position += overlap.length;
        ref_alignment_last_position += overlap.length;
    }

    if overlaps.positions.len() > 0 {
        let read_stop = overlaps.positions[overlaps.positions.len() - 1].search_start + overlaps.positions[overlaps.positions.len() - 1].length;
        let ref_stop = overlaps.positions[overlaps.positions.len() - 1].ref_start + overlaps.positions[overlaps.positions.len() - 1].length;
        //println!("LASTBIT read_alignment_last_position : {} ref_alignment_last_position : {}",read_stop,ref_stop);
        if read_stop < search_string.len() {

            // look back to see what segment we haven't aligned in the read
            let read_slice = slice_for_alignment(&search_string, read_alignment_last_position, search_string.len());
            let ref_slice = slice_for_alignment(&reference, ref_alignment_last_position, reference.len());
            let alignment = unaligned_segment_to_alignment(&read_slice, &ref_slice, min_alignment_seg_length);

            let read_ref_aligned_length = read_ref_alignment_lengths(&alignment);
            read_alignment_last_position += read_ref_aligned_length.0;
            ref_alignment_last_position += read_ref_aligned_length.1;

            alignmentTags.extend(alignment);
        }
    }

    alignmentTags
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
        assert!(reference.suffixTable.contains("AAT"));
        assert!(!reference.suffixTable.contains("TAAT"));
    }

    #[test]
    fn find_greedy_non_overlapping_segments_test() {
        let refseq = String::from("NNNNNNNNCATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCATATGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCGTAGAAGAAAGTTAGCTCTTGTGCGAGCTACAGGAACGATGTTTGATTAGAGTAAGCAGAGGACAAGGGCTCGCGTGCAGCCGAAGTTTGGCCGGTACTCTCCAACCGTTAACAACAACACCTTTCATCGAAATCCGCTTGGTAACAACACTAGGATGTTTCGACGCACTCGATAACCGGGAAACCAAGAGAAGTTTCCGAGCGCCACAGCCCAACTTAACTTCGCCATGTTTGAGACACGCGATCGGGACCACAAGACGTTACTCTTTGGGACCGCCGTAAGCGAGTATAAGCAGATTGTGTTTCGGCGTCCAAGTTGCCGTCAAAAGCTTACTGAGTTTCGCTGCCGGCGGAATGCTATTAATCGCGCCTACTTTCATGGAAGACGTTTCCGCTAAAATGGACTGGTGTTTCATGTCGGGAGCCGCTTTGTAAGATGGAGCTACTTTCCAGTCTGAGTTCATGCGAGAACACCACGAGAGTTTCATGAGTGGCCTCTCGAATCAACAGTCTACAAGTTTGGAGTTATCCGACACATCAAAACCAGCCATTCGTTTCATGAGATGGATCGCATACTAACCTTAACGGAGTTTGTAGTCACGGACGAACGATAAACGAGCATAATCTTTCGAGTCGGAGGCGAGTACTTAACGGATATAACGTTTCGTGCCAATGTTAACCGGTCAACTACCACTCAGTTTCTTGTTCATCATAACACTGAAACTGAGATCGTCTTTGGTGCAATTCCAATACGGCTAACTTACGCATACTTTGATGACGCCGTGATTATATCAAGAACCTACCGCTTTCATGGCGGTAACGGTATCCAAAGAATTGGTGTGTTTCGTGCATGCAGTGTCGGACTAAGACTAGGAATGTTTGCAGTGGCCGCAGTTATTCCAAGAGGATGCTTCTTTCCAGCTAACGGTCGCCTAATAAGATGTAACTGGTTTCTTGAGGAGGCATGTACCCGAAGCTTAGAGTAGTCTCCTCTATGAATACTGAAGGACTTGCGTAGTTATGTACAAGCTCACCAACGGACGGGTGCTTCCACATATAACGTTAGCATCTCGTGTGCTATTCGTAAGAGTTTCTAGTCACGGACGAACGATAAAGTACCAACGCCTTTCATGAGTGGCCTCTCGAATCAAGTGATCGGACCTTTGGACGCACTCGATAACCGGGAAGTTATCCAGACTTTCGTGCCAATGTTAACCGGTCAATAAGAGCTACCTTTGATGACGCCGTGATTATATCAATACGCTTCTGGTTTGGGCGTCCAAGTTGCCGTCAAATAGTAGTGACCTTTGCAGTCTGAGTTCATGCGAGAATCACCGCGAAGTTTGTTGTTCATCATAACACTGAAATCCGCAATTAGTTTCCAGCTAACGGTCGCCTAATAATCGGTAGCACGTTTCCGAGCGCCACAGCCCAACTTATCTTACACACGTTTCATCGAAATCCGCTTGGTAACATGCAAGTGTAGTTTGCAGTGGCCGCAGTTATTCCAATGTGTGTGAGCTTTCGAGTTATCCGACACATCAAACAACCGATTAACTTTCGTGCATGCAGTGTCGGACTACAAGAATAGTGCTTTGATGGCGGTAACGGTATCCAACACACTATTACCTTTCATGTCGGGAGCCGCTTTGTACAGTGTGCTACCTTTGGTGCAATTCCAATACGGCTACATAGCTAACGGTTTCTTGAGGAGGCATGTACCCGACCGCCTTATTCGTTTGATGGAAGACGTTTCCGCTAACCGGTTCTAAGGTTTCGGACCGCCGTAAGCGAGTATCCTAGTACATGGTTTCGCTGCCGGCGGAATGCTATTCCTGCTAACGAGTTTCATGAGATGGATCGCATACTACGAAGCGTCATGTTTGGCCGGTACTCTCCAACCGTTCGAGCTTCTTCCTTTCGAGTCGGAGGCGAGTACTTACGATGATAGAGCTTTCAGACACGCGATCGGGACCACCGCAGCTAACAGTAATAGGACCGACCGACCGTTCGATTCAGGGAGATTGCCCTACACTATGCGGCAGCTGGCATAGACTCCTAAGGAGATGCGTACTTGTTAAATAGGACTCTTTCATCGAAATCCGCTTGGTAACCGCTAGGTTACGTTTGTTGTTCATCATAACACTGAACGTAACTATGTCTTTCGAGTTATCCGACACATCAAACTAAGTATGAGCTTTCCGAGCGCCACAGCCCAACTTCTAGCTAATCTCTTTGGCCGGTACTCTCCAACCGTTCTATTATGCCTGTTTGGCTGCCGGCGGAATGCTATTCTCCTGCTACACTTTGTAGTCACGGACGAACGATAACTGCTTAGAACCTTTGCAGTGGCCGCAGTTATTCCACTTAACGCGGAGTTTGGACGCACTCGATAACCGGGAGAACATTAGCTCTTTCATGAGATGGATCGCATACTAGAAGACATTAGGTTTCATGACGCCGTGATTATATCAGAAGTGTTACGGTTTCGGCGTCCAAGTTGCCGTCAAGAATTGATGAGCTTTCGTGCCAATGTTAACCGGTCAGACATTATAGCCTTTGGTGCAATTCCAATACGGCTAGACTCACACGACTTTGGAGTCGGAGGCGAGTACTTAGAGTTGTTGACGTTTCCAGCTAACGGTCGCCTAATAGATGATAGAACGTTTGCAGTCTGAGTTCATGCGAGAGCAATAAGCTACTTTCGGACCGCCGTAAGCGAGTATGCGATTAAGTAGTTTGATGTCGGGAGCCGCTTTGTAGGATACTCGACGTTTCAGACACGCGATCGGGACCACGGTAATATGACGTTTCTTGAGGAGGCATGTACCCGAGGTGTTGCATGGTTTCATGGAAGACGTTTCCGCTAAGTAACGGTATTCTTTGATGAGTGGCCTCTCGAATCAGTACCAGACTTGTTTCATGGCGGTAACGGTATCCAAGTAGCATAATCGTTTGGTGCATGCAGTGTCGGACTAGTATGCTATCGCAGGAGGATGGGGCAGGACAGGACGCGGCCACCCCAGGCCTTTCCAGAGCAAACCTGGAGAAGATTCACAATAGACAGGCCAAGAAACCCGGTGCTTCCTCCAGAGCCGTTTAAAGCTGATATGAGGAAATAAAGAGTGAACTGGAAGGATCCGATATCGCCACCGTGGCTGAATGAGACTGGTGTCGACCTGTGCCT").as_bytes().to_owned();
        let read = String::from("GTATTGCTCATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCAAATGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCATAGAAGAAGCTTAGCTCTTGTGCGAGCTACAGGAACGATGTTTGATTAGAGTAAGCAGAGGACAAGGGCTCGCGTGCAGCCGAAGTTTGGCCGGTACTCTCCATACAGTGTGCTACCTTTGGTGCAATTCCAATACGGCTACATAGCTAACGGTTTCTTGAGGAGGCATGTACCCGACCGCCTTATTCGTTTGATGGAAGACGTTTCCGCTAACCGGTTCTAAGGTTTCGGGACCGCCGTAAGCGATTGATGAGCTTTCGTGCCAATGTTAACCGGTCAGACATTATAGCCTTTGGTGCAATTCCAATACGGTAATATGACGTTTCTTGAGGAGGCATGTACCCGAGGTGTTGCATGGTTTCATGGAAGACGTTTCCGCTAAGTAACGGTATTCTTTGATGAGTGGCCTCTCGAATCAGTACCAGACTTGTTTCATGGCGGTAACGGTATCCAAGTAGCATAATCGTTTGGTGCATGCAGTGTCGGACTAGTATGCTATCGCAGGAGGATGGGGCAGGACAGGACGCGGCCACCCCAGGCCTTTCCAGAGCAAACCTGGAGAAGATTCACAATAGACAGGCCAAGAAACCCGGTGCTTCCTCCAGAGCCGTTTAAAGCTGATATGAGGAAATAAAGAGTGAACTGGAAGGATCCCATATCGACAATACGTAACTGAACGAAGTACACCAGTATT").as_bytes().to_owned();
        let reference = find_seeds(&refseq, 20);

        let hits = find_greedy_non_overlapping_segments(&read, &refseq, &reference);

        for hit in hits.positions {
            println!("SEEEEEDS ref: {} search: {}, length: {}, endref: {}, endsearch: {}\n", hit.ref_start, hit.search_start, hit.length, hit.ref_start + hit.length, hit.search_start + hit.length);
        }
    }
    // align_with_anchors(search_string: &Vec<u8>, reference: &Vec<u8>, seeds: &SuffixTableLookup, min_alignment_seg_length: usize)
    #[test]
    fn test_basic_align_with_anchors() {
        let refseq = String::from("NNNNNNNNCATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCATATGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCGTAGAAGAAAGTTAGCTCTTGTGCGAGCTACAGGAACGATGTTTGATTAGAGTAAGCAGAGGACAAGGGCTCGCGTGCAGCCGAAGTTTGGCCGGTACTCTCCAACCGTTAACAACAACACCTTTCATCGAAATCCGCTTGGTAACAACACTAGGATGTTTCGACGCACTCGATAACCGGGAAACCAAGAGAAGTTTCCGAGCGCCACAGCCCAACTTAACTTCGCCATGTTTGAGACACGCGATCGGGACCACAAGACGTTACTCTTTGGGACCGCCGTAAGCGAGTATAAGCAGATTGTGTTTCGGCGTCCAAGTTGCCGTCAAAAGCTTACTGAGTTTCGCTGCCGGCGGAATGCTATTAATCGCGCCTACTTTCATGGAAGACGTTTCCGCTAAAATGGACTGGTGTTTCATGTCGGGAGCCGCTTTGTAAGATGGAGCTACTTTCCAGTCTGAGTTCATGCGAGAACACCACGAGAGTTTCATGAGTGGCCTCTCGAATCAACAGTCTACAAGTTTGGAGTTATCCGACACATCAAAACCAGCCATTCGTTTCATGAGATGGATCGCATACTAACCTTAACGGAGTTTGTAGTCACGGACGAACGATAAACGAGCATAATCTTTCGAGTCGGAGGCGAGTACTTAACGGATATAACGTTTCGTGCCAATGTTAACCGGTCAACTACCACTCAGTTTCTTGTTCATCATAACACTGAAACTGAGATCGTCTTTGGTGCAATTCCAATACGGCTAACTTACGCATACTTTGATGACGCCGTGATTATATCAAGAACCTACCGCTTTCATGGCGGTAACGGTATCCAAAGAATTGGTGTGTTTCGTGCATGCAGTGTCGGACTAAGACTAGGAATGTTTGCAGTGGCCGCAGTTATTCCAAGAGGATGCTTCTTTCCAGCTAACGGTCGCCTAATAAGATGTAACTGGTTTCTTGAGGAGGCATGTACCCGAAGCTTAGAGTAGTCTCCTCTATGAATACTGAAGGACTTGCGTAGTTATGTACAAGCTCACCAACGGACGGGTGCTTCCACATATAACGTTAGCATCTCGTGTGCTATTCGTAAGAGTTTCTAGTCACGGACGAACGATAAAGTACCAACGCCTTTCATGAGTGGCCTCTCGAATCAAGTGATCGGACCTTTGGACGCACTCGATAACCGGGAAGTTATCCAGACTTTCGTGCCAATGTTAACCGGTCAATAAGAGCTACCTTTGATGACGCCGTGATTATATCAATACGCTTCTGGTTTGGGCGTCCAAGTTGCCGTCAAATAGTAGTGACCTTTGCAGTCTGAGTTCATGCGAGAATCACCGCGAAGTTTGTTGTTCATCATAACACTGAAATCCGCAATTAGTTTCCAGCTAACGGTCGCCTAATAATCGGTAGCACGTTTCCGAGCGCCACAGCCCAACTTATCTTACACACGTTTCATCGAAATCCGCTTGGTAACATGCAAGTGTAGTTTGCAGTGGCCGCAGTTATTCCAATGTGTGTGAGCTTTCGAGTTATCCGACACATCAAACAACCGATTAACTTTCGTGCATGCAGTGTCGGACTACAAGAATAGTGCTTTGATGGCGGTAACGGTATCCAACACACTATTACCTTTCATGTCGGGAGCCGCTTTGTACAGTGTGCTACCTTTGGTGCAATTCCAATACGGCTACATAGCTAACGGTTTCTTGAGGAGGCATGTACCCGACCGCCTTATTCGTTTGATGGAAGACGTTTCCGCTAACCGGTTCTAAGGTTTCGGACCGCCGTAAGCGAGTATCCTAGTACATGGTTTCGCTGCCGGCGGAATGCTATTCCTGCTAACGAGTTTCATGAGATGGATCGCATACTACGAAGCGTCATGTTTGGCCGGTACTCTCCAACCGTTCGAGCTTCTTCCTTTCGAGTCGGAGGCGAGTACTTACGATGATAGAGCTTTCAGACACGCGATCGGGACCACCGCAGCTAACAGTAATAGGACCGACCGACCGTTCGATTCAGGGAGATTGCCCTACACTATGCGGCAGCTGGCATAGACTCCTAAGGAGATGCGTACTTGTTAAATAGGACTCTTTCATCGAAATCCGCTTGGTAACCGCTAGGTTACGTTTGTTGTTCATCATAACACTGAACGTAACTATGTCTTTCGAGTTATCCGACACATCAAACTAAGTATGAGCTTTCCGAGCGCCACAGCCCAACTTCTAGCTAATCTCTTTGGCCGGTACTCTCCAACCGTTCTATTATGCCTGTTTGGCTGCCGGCGGAATGCTATTCTCCTGCTACACTTTGTAGTCACGGACGAACGATAACTGCTTAGAACCTTTGCAGTGGCCGCAGTTATTCCACTTAACGCGGAGTTTGGACGCACTCGATAACCGGGAGAACATTAGCTCTTTCATGAGATGGATCGCATACTAGAAGACATTAGGTTTCATGACGCCGTGATTATATCAGAAGTGTTACGGTTTCGGCGTCCAAGTTGCCGTCAAGAATTGATGAGCTTTCGTGCCAATGTTAACCGGTCAGACATTATAGCCTTTGGTGCAATTCCAATACGGCTAGACTCACACGACTTTGGAGTCGGAGGCGAGTACTTAGAGTTGTTGACGTTTCCAGCTAACGGTCGCCTAATAGATGATAGAACGTTTGCAGTCTGAGTTCATGCGAGAGCAATAAGCTACTTTCGGACCGCCGTAAGCGAGTATGCGATTAAGTAGTTTGATGTCGGGAGCCGCTTTGTAGGATACTCGACGTTTCAGACACGCGATCGGGACCACGGTAATATGACGTTTCTTGAGGAGGCATGTACCCGAGGTGTTGCATGGTTTCATGGAAGACGTTTCCGCTAAGTAACGGTATTCTTTGATGAGTGGCCTCTCGAATCAGTACCAGACTTGTTTCATGGCGGTAACGGTATCCAAGTAGCATAATCGTTTGGTGCATGCAGTGTCGGACTAGTATGCTATCGCAGGAGGATGGGGCAGGACAGGACGCGGCCACCCCAGGCCTTTCCAGAGCAAACCTGGAGAAGATTCACAATAGACAGGCCAAGAAACCCGGTGCTTCCTCCAGAGCCGTTTAAAGCTGATATGAGGAAATAAAGAGTGAACTGGAAGGATCCGATATCGCCACCGTGGCTGAATGAGACTGGTGTCGACCTGTGCCT").as_bytes().to_owned();
        let read = String::from("GTATTGCTCATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAACGAAGAGTAACCGTTGCTAGGAGAGACCAAATGTCTAGAGAAAGGTACCCTATCCTTTCGAATGGTCCACGCATAGAAGAAGCTTAGCTCTTGTGCGAGCTACAGGAACGATGTTTGATTAGAGTAAGCAGAGGACAAGGGCTCGCGTGCAGCCGAAGTTTGGCCGGTACTCTCCATACAGTGTGCTACCTTTGGTGCAATTCCAATACGGCTACATAGCTAACGGTTTCTTGAGGAGGCATGTACCCGACCGCCTTATTCGTTTGATGGAAGACGTTTCCGCTAACCGGTTCTAAGGTTTCGGGACCGCCGTAAGCGATTGATGAGCTTTCGTGCCAATGTTAACCGGTCAGACATTATAGCCTTTGGTGCAATTCCAATACGGTAATATGACGTTTCTTGAGGAGGCATGTACCCGAGGTGTTGCATGGTTTCATGGAAGACGTTTCCGCTAAGTAACGGTATTCTTTGATGAGTGGCCTCTCGAATCAGTACCAGACTTGTTTCATGGCGGTAACGGTATCCAAGTAGCATAATCGTTTGGTGCATGCAGTGTCGGACTAGTATGCTATCGCAGGAGGATGGGGCAGGACAGGACGCGGCCACCCCAGGCCTTTCCAGAGCAAACCTGGAGAAGATTCACAATAGACAGGCCAAGAAACCCGGTGCTTCCTCCAGAGCCGTTTAAAGCTGATATGAGGAAATAAAGAGTGAACTGGAAGGATCCCATATCGACAATACGTAACTGAACGAAGTACACCAGTATT").as_bytes().to_owned();
        let reference = find_seeds(&refseq, 20);

        let fwd_score_mp = find_greedy_non_overlapping_segments(&read, &refseq, &reference);

        let hits = align_with_anchors(&read, &refseq, &reference, 10, &fwd_score_mp);

        for hit in hits {
            print!("{}",hit);
        }
    }
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


