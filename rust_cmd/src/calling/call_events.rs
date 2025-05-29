//! # Call Events Module
//!
//! This module contains functionality for processing BAM files to extract and analyze
//! genetic events from aligned sequencing reads. It provides tools for:
//!
//! - Extracting alignment information from BAM records
//! - Parsing CIGAR strings to identify insertions, deletions, and mismatches
//! - Mapping events to target regions within reference sequences
//! - Extracting custom BAM tags for read metadata
//! - Generating event summaries for downstream analysis
//!
//! The primary entry point is the `BamCallingParser` struct, which processes
//! BAM files according to a sequence layout configuration and outputs
//! tabular data describing the events found in each read.
//!
//! ## Key Components
//!
//! - `ExtractorTags`: Enum for custom BAM tags containing read metadata
//! - `BaseModifications`: Enum representing nucleotide substitutions
//! - `FullAlignment`: Enum representing alignment events (matches, mismatches, indels)
//! - `TargetRange`: Struct representing genomic regions of interest
//! - `BamCallingParser`: Main parser for processing BAM files

use std::cmp::{min};
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use bio::bio_types::genome::AbstractInterval;
use bio::bio_types::sequence::{SequenceRead};
use itertools::Itertools;
use rust_htslib::bam;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{Read, Record};
use rust_htslib::bam::record::{Aux, Cigar, CigarStringView};
use std::io::Write;
use crate::read_strategies::sequence_layout::{ReferenceRecord, SequenceLayout};
use crate::reference::fasta_reference::ReferenceManager;


/// Custom BAM tags used to store metadata about sequencing reads.
///
/// These tags contain information about read counts, alignment scores,
/// barcode sequences, and target regions. The tags are stored in BAM
/// auxiliary fields and extracted during processing.
#[derive(Clone, PartialEq)]
#[allow(dead_code)]
pub enum ExtractorTags {
    // AC represents the read count with a value of type u64, coded as 'a', 'c'
    // Old value: [b'r', b'c']
    AC { value: u64 },

    // AD represents the downsampled read count with a value of type u64, coded as 'a', 'd'
    // Old value: [b'd', b'c']
    AD { value: u64 },

    // AN represents the read names with a value of type String, coded as 'a', 'n'
    // Old value: [b'a', b'r']
    AN { value: String },

    // AR represents the alignment rate with a value of type f64, coded as 'a', 'r'
    // Old value: [b'r', b'm']
    AR { value: f64 },

    // AS represents the alignment score with a value of type f64, coded as 'a', 's'
    // Old value: [b'r', b's']
    AS { value: f64 },

    // BARCODE represents a extracted barcode sequence with a reference name and nucleotide value, both of type String
    // Reference name starts with 'b' and ends in 0-9
    BARCODE { reference_name: String, value: String },

    // TARGET represents a CRISPR / editing target with a reference name, wild-type sequence, and editor type
    // Reference name starts with 'e-z' and ends in 0-9
    // 'a-z' are extracted tags
    TARGET { reference_name: String, target_offset_in_reference: usize},
}

/// A pair linking an extractor tag with its string value.
#[derive(Clone, PartialEq)]
struct ExtractorPair {
    extractor: ExtractorTags,
    value: String,
}

/// State tracking for target extraction operations.
#[allow(dead_code)]
pub struct TargetExtractorState {
    current_target_id: usize,
}

impl ExtractorTags {
    /// Extracts and parses a custom BAM tag based on its two-byte identifier.
    ///
    /// This function matches two-byte tag identifiers against known patterns and
    /// constructs the appropriate `ExtractorTags` variant. It handles read counts,
    /// alignment scores, barcode sequences, and target region tags.
    ///
    /// # Arguments
    /// * `tag` - Two-byte tag identifier from BAM auxiliary field
    /// * `value` - String value associated with the tag
    /// * `reference_name` - Name of the reference sequence
    /// * `sequence_layout` - Layout configuration for validation
    ///
    /// # Returns
    /// `Some(ExtractorTags)` if the tag is recognized, `None` otherwise
    pub fn extract_matching_tag(tag: &[u8; 2], value: &String, reference_name: &String, sequence_layout: &SequenceLayout) -> Option<ExtractorTags> {
        //println!("tag {} {} value {}",char::from(tag[0]).as_ascii().unwrap(),char::from(tag[1]).as_ascii().unwrap(),value);
        match tag {
            &[b'r', b'c'] => { Some(ExtractorTags::AC{value: value.parse::<u64>().expect("Unable to parse integer from rc tag")})}, // TODO fix this,
            &[b'r', b'm'] => { Some(ExtractorTags::AR{value: value.parse::<f64>().expect("Unable to parse integer from rc tag")})}, // TODO fix this,
            &[b'r', b's'] => { Some(ExtractorTags::AS{value: value.parse::<f64>().expect("Unable to parse integer from rc tag")})}, // TODO fix this

            &[b'a', b'c'] => { Some(ExtractorTags::AC{value: value.parse::<u64>().expect("Unable to parse integer from rc tag")})},
            &[b'a', b'd'] => { Some(ExtractorTags::AC{value: value.parse::<u64>().expect("Unable to parse integer from rc tag")})},
            &[b'a', b'n'] => { Some(ExtractorTags::AC{value: value.parse::<u64>().expect("Unable to parse integer from rc tag")})},
            &[b'a', b'r'] => { None}, // TODO process read names
            &[b'a', b's'] => { Some(ExtractorTags::AS{value: value.parse::<f64>().expect("Unable to parse integer from rc tag")})},
            &[x, _y] if x == b'b' => {
                // make sure the reference is in the lookup table
                assert!(sequence_layout.references.contains_key(reference_name));
                Some(ExtractorTags::BARCODE{
                    reference_name: reference_name.clone(),
                    value: value.clone(),
                })
            },
            &[x, y] if x >= b'e' && x <= b'z' => {

                assert!(y >= b'0' && y <= b'9');
                let second_offset : usize = (y as usize) - 48; // 48 = '0' in ASCII
                let target_offset = second_offset * ((x as usize) - 101); // 101 = 'e' in ASCII

                // make sure the reference is in the lookup table
                assert!(sequence_layout.references.contains_key(reference_name));
                let _ref_obj = sequence_layout.references.get(reference_name).unwrap();

                Some(ExtractorTags::TARGET{
                    reference_name: reference_name.clone(),
                    target_offset_in_reference: target_offset,
                })
            },
            _ => {
                // we don't know what this is
                None
            },
        }
    }

}

/// Represents all possible nucleotide substitutions (from -> to).
///
/// Each variant represents a specific base change, where the first letter
/// is the reference base and the second letter is the observed base.
/// For example, `AC` represents A in the reference changed to C in the read.
#[derive(Copy, Clone, Eq, Hash, PartialEq, Debug)]
enum BaseModifications {
    AC,
    AG,
    AT,
    AN,
    CA,
    CG,
    CT,
    CN,
    GA,
    GC,
    GT,
    GN,
    TA,
    TC,
    TG,
    TN,
    NA,
    NC,
    NG,
    NT,
}

impl BaseModifications {
    /// Returns the target nucleotide for this base modification.
    ///
    /// For a modification like `AC` (A->C), this returns `C` (the target base).
    pub fn modified_base(&self) -> u8 {
        match self {
            BaseModifications::AC => b'C',
            BaseModifications::AG => b'G',
            BaseModifications::AT => b'T',
            BaseModifications::AN => b'N',
            BaseModifications::CA => b'A',
            BaseModifications::CG => b'G',
            BaseModifications::CT => b'T',
            BaseModifications::CN => b'N',
            BaseModifications::GA => b'A',
            BaseModifications::GC => b'C',
            BaseModifications::GT => b'T',
            BaseModifications::GN => b'N',
            BaseModifications::TA => b'A',
            BaseModifications::TC => b'C',
            BaseModifications::TG => b'G',
            BaseModifications::TN => b'N',
            BaseModifications::NA => b'A',
            BaseModifications::NC => b'C',
            BaseModifications::NG => b'G',
            BaseModifications::NT => b'T',
        }
    }

    /// Creates a `BaseModifications` variant from reference and observed bases.
    ///
    /// # Arguments
    /// * `base1` - The reference nucleotide
    /// * `base2` - The observed nucleotide in the read
    ///
    /// # Returns
    /// The corresponding `BaseModifications` variant
    ///
    /// # Panics
    /// Panics if either base is not a valid nucleotide (A, C, G, T, N)
    pub fn from_modified_bases(base1: u8, base2: u8) -> BaseModifications {
        match (base1.to_ascii_uppercase(), base2.to_ascii_uppercase()) {
            (b'A', b'C') => BaseModifications::AC,
            (b'A', b'G') => BaseModifications::AG,
            (b'A', b'T') => BaseModifications::AT,
            (b'A', b'N') => BaseModifications::AN,
            (b'C', b'A') => BaseModifications::CA,
            (b'C', b'G') => BaseModifications::CG,
            (b'C', b'T') => BaseModifications::CT,
            (b'C', b'N') => BaseModifications::CN,
            (b'G', b'A') => BaseModifications::GA,
            (b'G', b'C') => BaseModifications::GC,
            (b'G', b'T') => BaseModifications::GT,
            (b'G', b'N') => BaseModifications::GN,
            (b'T', b'A') => BaseModifications::TA,
            (b'T', b'C') => BaseModifications::TC,
            (b'T', b'G') => BaseModifications::TG,
            (b'T', b'N') => BaseModifications::TN,
            (b'N', b'A') => BaseModifications::NA,
            (b'N', b'C') => BaseModifications::NC,
            (b'N', b'G') => BaseModifications::NG,
            (b'N', b'T') => BaseModifications::NT,
            (b1, b2) => panic!("Unable to convert bases {} and {}", b1 as char, b2 as char)
        }
    }
}


/// Represents different types of alignment events between a read and reference.
///
/// Each variant contains the reference position and additional information
/// about the event type:
/// - `Insertion`: Extra bases in the read not present in reference
/// - `Deletion`: Bases missing from the read that are present in reference  
/// - `Match`: Exact sequence match between read and reference
/// - `Mismatch`: Sequence differences between read and reference
#[derive(PartialEq, Debug, Eq, Hash, Clone)]
enum FullAlignment {
    Insertion(usize, Vec<u8>),
    Deletion(usize, u32),
    Match(usize, u32),
    Mismatch(usize, Vec<BaseModifications>),
}

impl FullAlignment {
    /// Returns the genomic range (start, end) covered by this alignment event.
    pub fn to_range(&self) -> (usize, usize) {
        match self {
            FullAlignment::Insertion(x, y) => (*x as usize, *x as usize + y.len()),
            FullAlignment::Deletion(x, y) => (*x as usize, *x as usize + *y as usize),
            FullAlignment::Match(x, y) => (*x as usize, *x as usize + *y as usize),
            FullAlignment::Mismatch(x, y) => (*x as usize, *x as usize + y.len())
        }
    }


    /// Encodes the alignment event as a string for output.
    ///
    /// Returns a formatted string describing the event, or `None` for matches
    /// which are typically not reported. Format examples:
    /// - Insertion: "4I+10+ACGT" (4 bases inserted at position 10: ACGT)
    /// - Deletion: "5D+15" (5 bases deleted starting at position 15)
    /// - Mismatch: "2S+20+CG" (2 substitutions at position 20: to CG)
    pub fn to_encoding(&self) -> Option<String> {
        match self {
            FullAlignment::Insertion(x, y) => {
                Some(y.len().to_string() + "I+" + x.to_string().as_str() + "+" + std::str::from_utf8(y).unwrap())
            }
            FullAlignment::Deletion(x, y) => {
                Some(y.to_string() + "D+" + x.to_string().as_str())
            }
            FullAlignment::Match(_, _) => {
                None
            }
            FullAlignment::Mismatch(x, y) => {
                Some(y.len().to_string() + "S+" + x.to_string().as_str() + "+" + y.iter().map(|bm| bm.modified_base() as char).join("").as_str())
            }
        }
    }
}

/// A reference sequence with its name and nucleotide sequence.
#[derive(Clone, Eq, Hash, PartialEq)]
struct Reference {
    name: String,
    sequence: String,
}

/// Compares two nucleotide sequences and identifies contiguous regions of matches and mismatches.
///
/// This function performs a base-by-base comparison between reference and read sequences,
/// grouping consecutive matches and mismatches into separate alignment events.
///
/// # Arguments
/// * `reference` - Reference sequence as bytes
/// * `sequence` - Read sequence as bytes to compare against reference
/// * `reference_offset` - Starting position in the reference coordinate system
///
/// # Returns
/// Vector of `FullAlignment` events representing matches and mismatches
fn breakup_nucleotide_sequences(reference: &[u8], sequence: &[u8], reference_offset: &usize) -> Vec<FullAlignment> {
    let mut return_sections = Vec::new();
    let mut current_section = Vec::new();
    let mut in_match = false;
    let mut segment_length = 0;

    println!("{} {} ", String::from_utf8(reference.to_vec()).unwrap(), String::from_utf8(sequence.to_vec()).unwrap());

    for (position, (reference_base, read_base)) in reference.iter().zip(sequence.iter()).enumerate() {
        let section_start = reference_offset + (position as usize) - (current_section.len() as usize);

        if reference_base.to_ascii_uppercase() == read_base.to_ascii_uppercase() {
            if in_match {
                segment_length += 1;
            } else if !current_section.is_empty() {
                return_sections.push(FullAlignment::Mismatch(section_start, current_section.clone()));
                segment_length = 1;
                current_section.clear();
            }
            in_match = true;
        } else {
            if !in_match {
                segment_length += 1;
            } else {
                if segment_length > 0 {
                    return_sections.push(FullAlignment::Match(section_start, segment_length));
                }
                segment_length = 1;
                current_section.clear();
            }
            current_section.push(BaseModifications::from_modified_bases(reference_base.clone(), read_base.clone()));
            assert_eq!(current_section.len(),segment_length as usize);
            in_match = false;
        }
    }

    // Add the last section
    let section_start = reference_offset + reference.len()- (current_section.len());
    if in_match {
        return_sections.push(FullAlignment::Match(section_start, segment_length));
    } else {
        return_sections.push(FullAlignment::Mismatch(section_start, current_section.clone()));
    }

    return_sections
}

/// Extracts alignment events from a CIGAR string and sequence data.
///
/// This function processes CIGAR operations to generate detailed alignment events
/// including matches, mismatches, insertions, and deletions. It handles soft clipping
/// at read ends by optionally realigning those regions.
///
/// # Arguments
/// * `ref_start` - Starting position on the reference sequence
/// * `reference_sequence` - Reference DNA sequence as bytes
/// * `read_seq` - Read DNA sequence as bytes
/// * `cigar` - CIGAR string describing the alignment
/// * `realign_soft_clipped_ends` - Whether to realign soft-clipped regions
///
/// # Returns
/// Vector of `FullAlignment` events describing the complete alignment
///
/// # Panics
/// Panics on unsupported CIGAR operations (RefSkip, HardClip, Pad, Equal, Diff)
fn extract_read_cigar_elements(ref_start: &usize, reference_sequence: &Vec<u8>, read_seq: &Vec<u8>, cigar: &CigarStringView, realign_soft_clipped_ends: &bool) -> Vec<FullAlignment> {
    let mut ref_pos: usize = if *ref_start == 0 {0} else {*ref_start - 1}; // positions are 1-based, offsets needed here are 0-based
    let mut read_pos: usize = 0;
    let mut alignments = Vec::new();

    println!("red read {} {} ref_pos {} cigar {:?}", String::from_utf8(reference_sequence.to_vec()).unwrap(), String::from_utf8(read_seq.to_vec()).unwrap(),ref_pos,cigar);

    cigar.iter().enumerate().for_each(|(index,x)| {
        println!("x {} read pos {} ref pos {}",x.to_string(),read_pos, ref_pos);
        match x {
            Cigar::Match(ln) => {
                let ref_seq = &reference_sequence[(ref_pos as usize)..((ref_pos + *ln as usize) as usize)];
                let read_seq = &read_seq[(read_pos as usize)..((read_pos + *ln as usize) as usize)];
                alignments.extend(breakup_nucleotide_sequences(ref_seq, read_seq, &ref_pos));
                ref_pos += *ln as usize;
                read_pos += *ln as usize;
            }
            Cigar::Ins(ln) => {
                alignments.push(FullAlignment::Insertion(ref_pos, read_seq[(read_pos as usize)..((read_pos + *ln as usize) as usize)].to_vec()));
                read_pos += *ln as usize;
            }
            Cigar::Del(ln) => {
                alignments.push(FullAlignment::Deletion(ref_pos, *ln));
                println!("del {} {} {} {}", String::from_utf8(reference_sequence.to_vec()).unwrap(), String::from_utf8(read_seq.to_vec()).unwrap(),ref_pos, ref_pos + *ln as usize);
                ref_pos += *ln as usize;
            }
            Cigar::SoftClip(ln) if index == 0 => {
                // we're looking backward -- we don't want walk off the end of the reference
                let backward_length = min(ref_pos, *ln as usize);

                let ref_seq = &reference_sequence[((ref_pos - backward_length) as usize)..((ref_pos) as usize)];
                let read_seq = &read_seq[(((read_pos + *ln as usize) - backward_length))..((read_pos + (*ln as usize)) as usize)];

                if !*realign_soft_clipped_ends {
                    // TODO: this right?
                    alignments.extend(breakup_nucleotide_sequences(ref_seq, read_seq, &ref_pos));
                } else {

                }
                read_pos += *ln as usize;
            }
            Cigar::SoftClip(ln) if index == cigar.len() -1  => {
                // we're looking forward
                let ref_forward_len = min(*ln as usize,reference_sequence.len() - ref_pos);

                let ref_seq = &reference_sequence[(ref_pos as usize)..((ref_pos + ref_forward_len) as usize)];
                let read_seq = &read_seq[(read_pos as usize)..((read_pos + ref_forward_len) as usize)];
                alignments.extend(breakup_nucleotide_sequences(ref_seq, read_seq, &ref_pos));
                read_pos += ref_forward_len;
                ref_pos += ref_forward_len;
            }
            Cigar::RefSkip(_ln) => {
                panic!("RefSkip Unsupported")
            }
            Cigar::HardClip(_ln) => {
                panic!("HardClip Unsupported")
            }
            Cigar::Pad(_ln) => {
                panic!("Pad Unsupported")
            }
            Cigar::Equal(_ln) => {
                panic!("Equal Unsupported")
            }
            Cigar::Diff(_ln) => {
                panic!("Diff Unsupported")
            }
            _ => {
                panic!("Unprocessed cigar element {}", x);
            }
        }
    });
    assert_eq!(reference_sequence.len(), ref_pos); // make sure we used up the whole reference

    alignments
}
/// Computes the reverse complement of a DNA sequence.
///
/// # Arguments
/// * `dna` - Input DNA sequence string
///
/// # Returns
/// Reverse complement sequence with A<->T and C<->G swaps
fn reverse_complement(dna: &str) -> String {
    dna.chars()
        .rev() // Reverse the sequence
        .map(|nucleotide| match nucleotide {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => 'N'
        })
        .collect()
}

/// Represents a genomic range of interest for event analysis.
///
/// Target ranges define regions where we want to detect and report
/// alignment events. The orientation indicates whether the target
/// is on the forward (true) or reverse (false) strand.
#[derive(PartialEq, PartialOrd, Ord, Debug, Eq, Hash, Clone, Copy)]
struct TargetRange {
    start: usize,
    end: usize,
    orientation: bool,
}

impl TargetRange {
    /// Creates a new target range with the specified coordinates and orientation.
    pub fn new(start: &usize, end: &usize, orientation: &bool) -> TargetRange {
        TargetRange { start: *start, end: *end, orientation: *orientation }
    }
    /// Checks if this target range intersects with the given coordinates.
    ///
    /// # Arguments
    /// * `start` - Start position to check for intersection
    /// * `end` - End position to check for intersection
    ///
    /// # Returns
    /// `true` if the ranges overlap, `false` otherwise
    pub fn intersect_position(&self, start: usize, end: usize) -> bool {
        if start < self.start {
            end > self.start && start < self.end
        } else {
            self.end > start && self.start < end
        }
    }

    /// Checks if this target range intersects with another target range.
    #[allow(dead_code)]
    pub fn intersect(&self, other: &TargetRange) -> bool {
        self.intersect_position(other.start, other.end)
    }
}

/// Maps target sequences to their positions within a reference sequence.
///
/// This structure contains the reference sequence, target sequences,
/// and a mapping from target indices to their genomic positions.
#[allow(dead_code)]
struct TargetPositions {
    reference: String,
    targets: Vec<String>,
    positions: HashMap<usize, Vec<TargetRange>>,
}

impl TargetPositions {
    /// Creates target positions from a sequence layout record.
    ///
    /// This function finds all occurrences of target sequences within the
    /// reference, including both forward and reverse complement matches.
    pub fn from_sequence_layout_entry(record: &ReferenceRecord) -> TargetPositions {
        //println!("Targets {}", record.targets.join(","));
        TargetPositions {
            reference: record.sequence.clone(),
            targets: record.targets.clone(),
            positions: TargetPositions::validate_target_positions(record.sequence.as_str(), &record.targets.clone()),
        }
    }

    /// Validates and maps target sequences to their positions in the reference.
    ///
    /// # Arguments
    /// * `reference` - The reference sequence to search within
    /// * `targets` - Vector of target sequences to locate
    ///
    /// # Returns
    /// HashMap mapping target indices to their genomic positions
    fn validate_target_positions(reference: &str, targets: &Vec<String>) -> HashMap<usize, Vec<TargetRange>> {
        let mut target_locations: HashMap<usize, Vec<TargetRange>> = HashMap::new();

        targets.iter().enumerate().for_each(|(index, target)| {
            let mut target_indices: Vec<TargetRange> = reference.match_indices(target).into_iter().map(|(pos, _str)| TargetRange::new(&pos, &(pos + target.len()), &true)).collect();
            let mut target_rev: Vec<TargetRange> = reference.match_indices(&reverse_complement(target.as_str())).into_iter().map(|(pos, _str)| TargetRange::new(&pos, &(pos + target.len()), &true)).collect();
            target_indices.append(&mut target_rev);
            target_indices.sort();
            target_locations.insert(index.clone(), target_indices);
        });

        target_locations
    }
}

/// Main parser for processing BAM files and extracting alignment events.
///
/// This parser reads BAM files, extracts alignment information, and maps
/// events to target regions defined in the sequence layout. It outputs
/// tabular data describing the events found in each read.
#[allow(dead_code)]
pub struct BamCallingParser<'a, 's, 't> {
    sequence_layout: SequenceLayout,
    reference_manager: ReferenceManager<'a, 's, 't>,
    target_positions: HashMap<String, TargetPositions>,
    ordered_target_ranges: HashMap<String, Vec<(TargetRange, String)>>, // keep a ordered list of target ranges for simplicity
}

impl BamCallingParser<'_, '_, '_> {
    /// Creates a new BAM calling parser from a sequence layout configuration.
    ///
    /// # Arguments
    /// * `sequence_layout_design` - Configuration defining references and targets
    ///
    /// # Returns
    /// Configured parser ready to process BAM files
    pub fn new(sequence_layout_design: &SequenceLayout) -> BamCallingParser {
        let rm = ReferenceManager::from_yaml_input(sequence_layout_design, 12, 6);

        let mut ordered_target_ranges = HashMap::new();
        let target_positions = sequence_layout_design.references.iter().map(|(name, reference)| {
            let mut ordered_targets: Vec<(TargetRange, String)> = Vec::new();

            let positions = TargetPositions::from_sequence_layout_entry(reference);
            positions.targets.iter().enumerate().for_each(|(index, target_name)| {
                positions.positions.get(&index).unwrap().iter().for_each(|target_pos| ordered_targets.push((target_pos.clone(), target_name.clone())));
            });

            ordered_targets.sort();
            ordered_target_ranges.insert(name.clone(), ordered_targets);

            (name.clone(), positions)
        }).collect();


        BamCallingParser {
            sequence_layout: sequence_layout_design.clone(),
            reference_manager: rm,
            target_positions,
            ordered_target_ranges,
        }
    }

    /// Extracts custom tags from a BAM record.
    ///
    /// # Arguments
    /// * `bam_entry` - BAM record to extract tags from
    /// * `reference` - Reference sequence name
    ///
    /// # Returns
    /// Vector of extracted tags found in the BAM record
    fn extraction_tags(&self, bam_entry: &Record, reference: &String) -> Vec<ExtractorTags> {
        let mut extractors = Vec::new();

        bam_entry.aux_iter().for_each(|c| {
            match c {
                Ok((aux_tag, aux_enum)) => {
                    // match the aux tag against our known list of tags
                    match aux_enum {
                        Aux::String(x) => {
                            match ExtractorTags::extract_matching_tag(&[aux_tag[0],aux_tag[1]],&x.to_string(),reference, &self.sequence_layout) {
                                None => {}
                                Some(x) => {
                                    extractors.push(x);
                                }
                            }
                        }
                        _ => {}
                    }

                }
                Err(_) => { panic!("Unable to read auxillary tags from BAM file") }
            }
        });
        extractors
    }

    /// Computes overlaps between reference target ranges and full alignment tokens.
    ///
    /// This function checks each `TargetRange` from `reference_targets` against
    /// each token in `full_alignment_tokens` to determine overlaps. Overlaps are
    /// determined based on the position intersection of the `TargetRange` with the
    /// range derived from each token. It then encodes these overlaps into strings.
    ///
    /// # Parameters
    /// - `reference_targets`: A reference to a vector of tuples, each containing a `TargetRange`
    ///   and an associated string identifier. The `TargetRange` typically specifies a range
    ///   within a sequence or dataset.
    /// - `full_alignment_tokens`: A reference to a vector of `FullAlignment` tokens, which are
    ///   used to check for positional overlaps with the `TargetRange`. Each token can be converted
    ///   to a range and optionally holds alignment information.
    ///
    /// # Returns
    /// A `Vec<String>` where each element corresponds to a comma-separated string of encoded
    /// overlaps for each `TargetRange`. If there are no overlaps for a particular range, the
    /// corresponding string will be "NONE".
    ///
    /// # Examples
    /// Suppose we have the following `TargetRange` and `FullAlignment` data structures:
    /// ```
    /// struct TargetRange {
    ///     start: usize,
    ///     end: usize,
    /// }
    ///
    /// enum FullAlignment {
    ///     Match(usize, usize),
    ///     Mismatch(usize, usize),
    /// }
    /// ```
    ///
    /// Example usage might look like this:
    /// ```
    /// let reference_targets = vec![(TargetRange{start: 1, end: 5}, "target1")];
    /// let tokens = vec![FullAlignment::Mismatch(2, 4)];
    /// let overlaps = target_overlaps(&reference_targets, &tokens);
    /// assert_eq!(overlaps, vec!["encoded_value"]);
    /// ```
    ///
    /// This function utilizes `target_pos.intersect_position` to determine if a `TargetRange`
    /// intersects with a range derived from each `FullAlignment` token. It filters out tokens
    /// based on the intersection results, then converts the intersecting tokens to a string
    /// encoding using `to_encoding`. If no valid encoding is found, an uppercase empty string
    /// is used, which is converted to "NONE" if the entire string for a target is empty.
    fn target_overlaps(reference_targets: &Vec<(TargetRange, String)>, full_alignment_tokens: &Vec<FullAlignment>) -> Vec<String> {
        let target_overlap = reference_targets.iter().map(|(target_pos, _name)| {
            full_alignment_tokens.iter().filter(|tk| {
                let range = tk.to_range();
                match tk {
                    FullAlignment::Match(_, _) => false,
                    _ => target_pos.intersect_position(range.0, range.1),
                }
            }).map(|tk| { tk.to_encoding().unwrap_or("".to_ascii_uppercase()) }).into_iter().collect::<Vec<String>>().join(",")
        }).collect::<Vec<String>>();

        target_overlap.into_iter().map(|st| if st == "" { "NONE".to_ascii_uppercase() } else { st }).collect()
    }


    fn convert_reference_sequence(seq: &[u8])-> Vec<u8> {
        seq.iter().map(|x| {match *x {
            b'A' | b'a' | b'C' | b'c' | b'G' | b'g' |b'T' | b't' => {*x}
            _ => b'N'
        }}).collect()
    }

    /// Processes a BAM file and outputs event analysis to a TSV file.
    ///
    /// This is the main entry point for BAM processing. It reads each record,
    /// extracts alignment events, maps them to target regions, and outputs
    /// a summary table.
    ///
    /// # Arguments
    /// * `bam_file` - Path to input BAM file
    /// * `output_file` - Path to output TSV file
    ///
    /// # Returns
    /// `Result<(), std::io::Error>` indicating success or failure
    pub fn output_bam_file_entries(&self, bam_file: &str, output_file: &str) -> std::io::Result<()> {
        let mut bam = bam::Reader::from_path(bam_file).unwrap();
        let header = bam::Header::from_template(bam.header());

        for (key, records) in header.to_hashmap() {
            for record in records {
                //println!("@{}\tSN:{}\tLN:{}", key, record["SN"], record["LN"]);
                println!("@{}\tSN:{:?}", key, record);
            }
        }

        let mut file = File::create(output_file)?;
        write!(file, "read\tref\talignment_start\ttarget_outcomes\n").expect("Unable to write to output file");

        let mut read_count = 0;

        for r in bam.records() {
            read_count += 1;
            if read_count % 100000 == 0 { println!("Read count {}", read_count); }

            match r {
                Ok(record) => {
                    let ref_name = record.contig().as_bytes();
                    let layout_record = self.sequence_layout.references.get(record.contig()).unwrap();

                    // BTreeMap to sort the keys
                    let extractor_id_to_name : BTreeMap<usize,String> = layout_record.umi_configurations.iter().map(|(k,v)| (v.order.clone(),k.clone())).collect();

                    let reference = self.reference_manager.reference_name_to_ref.get(ref_name).expect("Unable to find reference");
                    let ref_name = String::from_utf8(self.reference_manager.references.get(reference).unwrap().name.clone()).unwrap();
                    let reference_sequence = BamCallingParser::convert_reference_sequence(self.reference_manager.references.get(reference).unwrap().sequence.as_slice());
                    let alignment_start = record.reference_start() as usize;

                    println!("read name {} start {}",String::from_utf8(record.name().to_vec()).unwrap(),alignment_start);

                    // TODO output extractor tags
                    let extractor_tokens = self.extraction_tags(&record, &ref_name);
                    let _token_output : Vec<(String,Option<&ExtractorTags>)> = extractor_id_to_name.iter().map(|(id,name)| (name.clone(),extractor_tokens.get(*id))).into_iter().collect::<Vec<(String,Option<&ExtractorTags>)>>();

                    println!("read name {} start {}",String::from_utf8(record.name().to_vec()).unwrap(),alignment_start);

                    let cigar_tokens = extract_read_cigar_elements(&alignment_start, &reference_sequence, &record.seq().as_bytes(), &record.cigar(), &true);

                    // now we need to intersect any events with the known target lists
                    let reference_targets = self.ordered_target_ranges.get(record.contig()).unwrap();

                    let target_output = BamCallingParser::target_overlaps(&reference_targets, &cigar_tokens);

                    write!(file, "{}\t{}\t{}\t{}\n", String::from_utf8(record.name().to_vec()).unwrap(), ref_name, alignment_start, target_output.join(",")).expect("Unable to write to output file");
                }
                Err(_x) => {
                    panic!("Underlying BAM file error!")
                }
            };
        }
        Ok(())
    }
}


#[cfg(test)]
mod tests {
    use super::*;


    #[test]
    fn event_target_intersections() {
        let target_positions: Vec<(TargetRange,String)> = vec![(TargetRange::new(&25, &50, &true),"Range1".to_ascii_uppercase()),(TargetRange::new(&75, &100, &true),"Range2".to_ascii_uppercase())];
        let full_alignment_tokens: Vec<FullAlignment> = vec![FullAlignment::Mismatch(28, vec![BaseModifications::AG,BaseModifications::AC])];
        let overlap = BamCallingParser::target_overlaps(&target_positions,&full_alignment_tokens);

        println!("{:?} {:?}",overlap.get(0).unwrap(),overlap.get(1).unwrap());
        assert_eq!(overlap.len(),2);
        assert_eq!(overlap.get(0).unwrap(),&"2S+28+GC".to_ascii_uppercase());
        assert_eq!(overlap.get(1).unwrap(),&"NONE".to_ascii_uppercase());

        let full_alignment_tokens: Vec<FullAlignment> = vec![FullAlignment::Mismatch(10, vec![BaseModifications::AG,BaseModifications::AC])];
        let overlap = BamCallingParser::target_overlaps(&target_positions,&full_alignment_tokens);

        assert_eq!(overlap.len(),2);
        assert_eq!(overlap.get(0).unwrap(),&"NONE".to_ascii_uppercase());
        assert_eq!(overlap.get(1).unwrap(),&"NONE".to_ascii_uppercase());

    }

    #[test]
    fn target_intersections() {
        let target_pos1 = TargetRange::new(&0, &20, &true);
        let target_pos2 = TargetRange::new(&19, &39, &true);
        assert!(target_pos1.intersect(&target_pos2));
        assert!(target_pos2.intersect(&target_pos1));

        let target_pos2 = TargetRange::new(&20, &40, &true);
        assert!(!target_pos1.intersect(&target_pos2));
        assert!(!target_pos2.intersect(&target_pos1));

        let target_pos1 = TargetRange::new(&80, &100, &true);
        assert!(!target_pos1.intersect(&target_pos2));
        assert!(!target_pos2.intersect(&target_pos1));

        let target_pos1 = TargetRange::new(&10, &100, &true);
        assert!(target_pos1.intersect(&target_pos2));
        assert!(target_pos2.intersect(&target_pos1));
    }

    #[test]
    fn target_positions() {
        let reference = "TTTTACGTAACGTAACGTAACGTACGGTTTT".to_ascii_uppercase();
        let targets = vec!["ACGTAACGTAACGTAACGTACGG".to_ascii_uppercase()];
        let _reference_offset = 10;

        let positions = TargetPositions::validate_target_positions(&reference, &targets);
        assert!(positions.contains_key(&0));
        assert_eq!(positions.get(&0).unwrap().len(), 1);
        assert_eq!(positions.get(&0).unwrap()[0].start, 4);

        // two positions in serial
        let reference = "TTTTACGTAACGTAACGTAACGTACGGTTTTTTTTACGTAACGTAACGTAACGTACGGTTTT".to_ascii_uppercase();
        let _reference_offset = 10;

        let positions = TargetPositions::validate_target_positions(&reference, &targets);
        assert!(positions.contains_key(&0));
        assert_eq!(positions.get(&0).unwrap().len(), 2);
        assert_eq!(positions.get(&0).unwrap()[0].start, 4);
        assert_eq!(positions.get(&0).unwrap()[1].start, 35);

        // two targets with a reverse comp target
        let targets = vec!["AAAAATTTTTAAAAATTTTTCGG".to_ascii_uppercase(), "ACGTAACGTAACGTAACGTACGG".to_ascii_uppercase()];
        let reference = "TTTTACGTAACGTAACGTAACGTACGGTTTTTTTTACGTAACGTAACGTAACGTACGGTTTTAAAAATTTTTAAAAATTTTTCGGAAAACCGTACGTTACGTTACGTTACGT".to_ascii_uppercase();
        let _reference_offset = 10;

        let positions = TargetPositions::validate_target_positions(&reference, &targets);

        assert!(positions.contains_key(&0));
        assert_eq!(positions.get(&0).unwrap().len(), 1);
        assert_eq!(positions.get(&0).unwrap()[0].start, 62);

        assert!(positions.contains_key(&1));
        assert_eq!(positions.get(&1).unwrap().len(), 3);
        assert_eq!(positions.get(&1).unwrap()[0].start, 4);
        assert_eq!(positions.get(&1).unwrap()[1].start, 35);
        assert_eq!(positions.get(&1).unwrap()[2].start, 89);


        // two targets with a reverse comp target
        let targets = vec!["AAAAATTTTTAAAAATTTTTCGG".to_ascii_uppercase(), "ACGTAACGTAACGTAACGTACGG".to_ascii_uppercase()];
        let reference = "TTTTACGTAACGTAACGTAACGTACGGTTTTTTTTACGTAACGTAACGTAACGTACGGTTTTAAAAATTTTTAAAAATTTTTCGGAAAACCGTACGTTACGTTACGTTACGTCCGAAAAATTTTTAAAAATTTTTCAGAAAAATTTTTAAAAATTTTT".to_ascii_uppercase();
        let _reference_offset = 10;

        let positions = TargetPositions::validate_target_positions(&reference, &targets);

        assert!(positions.contains_key(&0));
        assert_eq!(positions.get(&0).unwrap().len(), 2);
        assert_eq!(positions.get(&0).unwrap()[0].start, 62);
        assert_eq!(positions.get(&0).unwrap()[1].start, 112);

        assert!(positions.contains_key(&1));
        assert_eq!(positions.get(&1).unwrap().len(), 3);
        assert_eq!(positions.get(&1).unwrap()[0].start, 4);
        assert_eq!(positions.get(&1).unwrap()[1].start, 35);
        assert_eq!(positions.get(&1).unwrap()[2].start, 89);
    }

    #[test]
    fn full_alignment_encoding() {
        assert_eq!(FullAlignment::Match(10,10).to_encoding(),None);
        assert_eq!(FullAlignment::Mismatch(10,vec![BaseModifications::AC,BaseModifications::AG]).to_encoding().unwrap().as_str(),"2S+10+CG");
        assert_eq!(FullAlignment::Insertion(1,vec![b'A',b'C',b'T',b'T']).to_encoding().unwrap().as_str(),"4I+1+ACTT");
        assert_eq!(FullAlignment::Deletion(1,400).to_encoding().unwrap().as_str(),"400D+1");



    }
}