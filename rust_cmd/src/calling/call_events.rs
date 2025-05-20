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


// Derive the Clone and PartialEq traits for the ExtractorTags enum
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

#[derive(Clone, PartialEq)]
struct ExtractorPair {
    extractor: ExtractorTags,
    value: String,
}

#[allow(dead_code)]
pub struct TargetExtractorState {
    current_target_id: usize,
}

impl ExtractorTags {
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
            (b1, b2) => panic!("Unable to convert bases {} and {}", b1, b2)
        }
    }
}


#[derive(PartialEq, Debug, Eq, Hash, Clone)]
enum FullAlignment {
    Insertion(usize, Vec<u8>),
    Deletion(usize, u32),
    Match(usize, u32),
    Mismatch(usize, Vec<BaseModifications>),
}

impl FullAlignment {
    pub fn to_range(&self) -> (usize, usize) {
        match self {
            FullAlignment::Insertion(x, y) => (*x as usize, *x as usize + y.len()),
            FullAlignment::Deletion(x, y) => (*x as usize, *x as usize + *y as usize),
            FullAlignment::Match(x, y) => (*x as usize, *x as usize + *y as usize),
            FullAlignment::Mismatch(x, y) => (*x as usize, *x as usize + y.len())
        }
    }


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

#[derive(Clone, Eq, Hash, PartialEq)]
struct Reference {
    name: String,
    sequence: String,
}

/// Compares two nucleotide sequences and breaks them up into sections of matches and mismatches.
///
/// This function iterates over two nucleotide sequences represented as byte vectors, comparing
/// each base of the `reference` sequence with the corresponding base in the `sequence`. It
/// categorizes each section into either a match or a mismatch, based on the comparison. The
/// function takes into account an offset for the reference sequence to properly align the
/// sequences for comparison.
///
/// # Parameters
/// - `reference`: A `Vec<u8>` representing the reference sequence of nucleotides.
/// - `sequence`: A slice of u8 (`&[u8]`) representing the sequence of nucleotides to compare
///   against the reference.
/// - `reference_offset`: A reference to a u32 value representing the offset at which the comparison
///   of the `reference` sequence starts.
///
/// # Returns
/// A `Vec<FullAlignment>` where each `FullAlignment` is either a `Match` or `Mismatch`.
/// `Match` contains the starting position and length of the match, while `Mismatch` contains
/// the starting position and the mismatched sequence as a `Vec<u8>`.
///
/// # Examples
/// ```
/// enum FullAlignment {
///     Match(u32, u32),       // Start position and length of the match
///     Mismatch(u32, Vec<u8>) // Start position and mismatched sequence
/// }
///
/// fn main() {
///     let reference = vec![65, 67, 71, 84, 65, 67, 71, 84]; // ACGTACGT
///     let sequence = [65, 67, 71, 84, 84, 71, 67, 116];     // ACGTTGCt
///     let reference_offset = 0;
///     let alignments = breakup_nucleotide_sequences(&reference, &sequence, &reference_offset);
///
///     for alignment in alignments {
///         match alignment {
///             FullAlignment::Match(start, len) => println!("Match at {}: Length {}", start, len),
///             FullAlignment::Mismatch(start, seq) => println!("Mismatch at {}: {:?}", start, seq),
///         }
///     }
/// }
/// ```
/// Note: Ensure the `FullAlignment` enum is defined in your code as it is used in the function's return type.
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

/// Extracts and returns a vector of `FullAlignment` elements from given reference and read sequences based on the CIGAR string.
///
/// This function processes a CIGAR string to align a read sequence against a reference sequence starting from a specified position. It supports matches, insertions, and deletions but panics on unsupported CIGAR operations like reference skips, soft clips, hard clips, pads, equals, and diffs.
///
/// # Parameters
/// - `ref_start`: A reference to a `u32` indicating the start position on the reference sequence where the alignment should begin.
/// - `reference_sequence`: A reference to a vector of `u8` bytes representing the reference DNA sequence.
/// - `read_seq`: A reference to a vector of `u8` bytes representing the read DNA sequence to be aligned.
/// - `cigar`: A reference to a `CigarStringView` which contains the CIGAR operations that describe how the read sequence aligns to the reference sequence.
///
/// # Returns
/// Returns a vector of `FullAlignment` enumerations that describe the full alignment between the reference and read sequences. This vector may contain elements representing matches, insertions, and deletions according to the CIGAR string.
///
/// # Panics
/// The function panics if it encounters unsupported CIGAR operations (`RefSkip`, `SoftClip`, `HardClip`, `Pad`, `Equal`, `Diff`). It is designed to work with a subset of CIGAR operations that represent simple alignments.
///
/// # Examples
/// ```
/// // Example usage of `extract_read_cigar_elements`
/// let ref_start = 0;
/// let reference_sequence = b"ACGTACGT".to_vec();
/// let read_seq = b"ACGTTACGT".to_vec();
/// let cigar = CigarStringView::from_string("8M1I".to_string()).unwrap(); // Simplified example; actual initialization may vary
///
/// let alignments = extract_read_cigar_elements(&ref_start, &reference_sequence, &read_seq, &cigar);
/// // Process `alignments` as needed
/// ```
///
/// # Notes
/// - The function asserts that the entire reference sequence is aligned by the end of the process, which might not always be the case in real-world scenarios. This assertion should be adjusted according to the specific requirements of the alignment algorithm being implemented.
fn extract_read_cigar_elements(ref_start: &usize, reference_sequence: &Vec<u8>, read_seq: &Vec<u8>, cigar: &CigarStringView, realign_soft_clipped_ends: &bool) -> Vec<FullAlignment> {
    let mut ref_pos: usize = *ref_start;
    let mut read_pos: usize = 0;
    let mut alignments = Vec::new();

    println!("red read {} {} ref_pos {}", String::from_utf8(reference_sequence.to_vec()).unwrap(), String::from_utf8(read_seq.to_vec()).unwrap(),ref_pos);

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
    assert_eq!(reference_sequence.len(), ref_pos);

    alignments
}

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

#[derive(PartialEq, PartialOrd, Ord, Debug, Eq, Hash, Clone, Copy)]
struct TargetRange {
    start: usize,
    end: usize,
    orientation: bool,
}

impl TargetRange {
    pub fn new(start: &usize, end: &usize, orientation: &bool) -> TargetRange {
        TargetRange { start: *start, end: *end, orientation: *orientation }
    }
    pub fn intersect_position(&self, start: usize, end: usize) -> bool {
        if start < self.start {
            end > self.start && start < self.end
        } else {
            self.end > start && self.start < end
        }
    }

    #[allow(dead_code)]
    pub fn intersect(&self, other: &TargetRange) -> bool {
        self.intersect_position(other.start, other.end)
    }
}

#[allow(dead_code)]
struct TargetPositions {
    reference: String,
    targets: Vec<String>,
    positions: HashMap<usize, Vec<TargetRange>>,
}

impl TargetPositions {
    pub fn from_sequence_layout_entry(record: &ReferenceRecord) -> TargetPositions {
        //println!("Targets {}", record.targets.join(","));
        TargetPositions {
            reference: record.sequence.clone(),
            targets: record.targets.clone(),
            positions: TargetPositions::validate_target_positions(record.sequence.as_str(), &record.targets.clone()),
        }
    }

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

#[allow(dead_code)]
pub struct BamCallingParser<'a, 's, 't> {
    sequence_layout: SequenceLayout,
    reference_manager: ReferenceManager<'a, 's, 't>,
    target_positions: HashMap<String, TargetPositions>,
    ordered_target_ranges: HashMap<String, Vec<(TargetRange, String)>>, // keep a ordered list of target ranges for simplicity
}

impl BamCallingParser<'_, '_, '_> {
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


    //fn process_alignment(reference: &str, query: &str, cigar: &str) -> Vec<String> {}
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
                    let reference_sequence = self.reference_manager.references.get(reference).unwrap().sequence.clone();
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