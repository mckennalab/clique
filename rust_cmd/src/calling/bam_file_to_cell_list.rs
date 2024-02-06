use std::collections::HashMap;
use bio::bio_types::genome::AbstractInterval;
use bio::bio_types::sequence::SequenceRead;
use rust_htslib::bam;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{Read, Record};
use rust_htslib::bam::record::{Aux, Cigar, CigarStringView};
use crate::reference::fasta_reference::ReferenceManager;
use phf::phf_map;

#[derive(Copy, Clone, Eq, Hash, PartialEq)]
pub enum ExtractorTag {
    E0,
    E1,
    E2,
    E3,
    E4,
    E5,
    E6,
    E7,
    E8,
    E9,
}

#[derive(Copy, Clone, Eq, Hash, PartialEq)]
struct ExtractorPair {
    extractor: ExtractorTag,
    value: u32,
}

static EXTRACTORTYPE: phf::Map<&'static str, ExtractorTag> = phf_map! {
        "e0" => ExtractorTag::E0,
        "e1" => ExtractorTag::E1,
        "e2" => ExtractorTag::E2,
        "e3" => ExtractorTag::E3,
        "e4" => ExtractorTag::E4,
        "e5" => ExtractorTag::E5,
        "e6" => ExtractorTag::E6,
        "e7" => ExtractorTag::E7,
        "e8" => ExtractorTag::E8,
        "e9" => ExtractorTag::E9,
};

#[derive(PartialEq, Debug, Eq, Hash, Clone)]
enum FullAlignment {
    Insertion(u32, Vec<u8>),
    Deletion(u32, u32),
    Match(u32, u32),
    Mismatch(u32, Vec<u8>),
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
fn breakup_nucleotide_sequences(reference: &[u8], sequence: &[u8], reference_offset: &u32) -> Vec<FullAlignment> {
    let mut sections = Vec::new();
    let mut current_section = Vec::new();
    let mut last_was_match = None;

    for (position, (reference_base, read_base)) in reference.iter().zip(sequence.iter()).enumerate() {
        let is_match = reference_base == read_base;
        let section_start = reference_offset + (position as u32) - (current_section.len() as u32);
        match last_was_match {
            Some(last) if last != is_match => {
                // End the current section and start a new one
                if last {
                    sections.push(FullAlignment::Match(section_start, current_section.len() as u32));
                } else {
                    sections.push(FullAlignment::Mismatch(section_start, current_section));
                }
                current_section = Vec::new();
            }
            _ => {}
        }

        current_section.push(read_base.clone());
        last_was_match = Some(is_match);
    }

    let section_start = reference_offset + (sequence.len() as u32) - (current_section.len() as u32);
    // Add the last section
    if let Some(last) = last_was_match {
        if last {
            sections.push(FullAlignment::Match(section_start, current_section.len() as u32));
        } else {
            sections.push(FullAlignment::Mismatch(section_start, current_section));
        }
    }

    sections
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
fn extract_read_cigar_elements(ref_start: &u32, reference_sequence: &Vec<u8>, read_seq: &Vec<u8>, cigar: &CigarStringView) -> Vec<FullAlignment> {
    let mut ref_pos = *ref_start;
    let mut read_pos = 0;
    let mut alignments = Vec::new();

    cigar.iter().for_each(|x| {
        match x {
            Cigar::Match(ln) => {
                alignments.extend(breakup_nucleotide_sequences(&reference_sequence[(ref_pos as usize)..((ref_pos + ln) as usize)], &read_seq[(read_pos as usize)..((read_pos + ln) as usize)], &ref_pos));
                ref_pos += ln;
                read_pos += ln;
            }
            Cigar::Ins(ln) => {
                alignments.push(FullAlignment::Insertion(ref_pos, read_seq[(read_pos as usize)..((read_pos + ln) as usize)].to_vec()));
                read_pos += ln;
            }
            Cigar::Del(ln) => {
                alignments.push(FullAlignment::Deletion(ref_pos, ln.clone()));
                ref_pos += ln;
            }
            Cigar::RefSkip(_ln) => {
                panic!("RefSkip Unsupported")
            }
            Cigar::SoftClip(_ln) => {
                panic!("SoftClip Unsupported")
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
        }
    });
    assert_eq!(reference_sequence.len(), ref_pos.try_into().unwrap());

    alignments
}




struct ReadTableEntry {
    name: String,
    extractions: Vec<ExtractorPair>,
    alignments: Vec<u32>,
    reference: u16,
}

struct ExtractionTable {
    extractions: HashMap<String, u32>,
    extractions_reverse: HashMap<u32, String>,
    extraction_entries: Vec<String>,
    extraction_current_index: u32,
}

impl ExtractionTable {
    pub fn new() -> ExtractionTable {
        ExtractionTable {
            extractions: HashMap::new(),
            extractions_reverse: HashMap::new(),
            extraction_entries: vec![],
            extraction_current_index: 0,
        }
    }

    pub fn add(&mut self, value: &str) -> u32 {
        let vl = String::from(value);
        match self.extractions.get(&vl) {
            None => {
                self.extractions.insert(vl.clone(), self.extraction_current_index);
                self.extractions_reverse.insert(self.extraction_current_index, vl.clone());
                self.extraction_entries.push(vl);
                let ret = self.extraction_current_index;
                self.extraction_current_index += 1;
                ret
            }
            Some(hit) => {
                *hit
            }
        }
    }
}

struct BAMStream<'a, 's, 't> {
    reference_manager: ReferenceManager<'a, 's, 't>,

}
impl BAMStream<'_, '_, '_>  {
    pub fn new(fasta_file: &String) -> Box<BAMStream> {

        // read in the fasta file, creating a new BAM registry from the file
        let rm = ReferenceManager::from(&fasta_file, 12, 6);
        Box::new(BAMStream {
            reference_manager: rm,
        })
    }
}

// the goal here is to make a large lookup table for all aspects we want to store for a read; then each
// read can be (relatively) lightweight
struct FullBAMAlignmentRegistry<'a, 's, 't> {
    extractions: HashMap<ExtractorTag, ExtractionTable>,

    alignments: HashMap<FullAlignment, u32>,
    alignments_reverse: HashMap<u32, FullAlignment>,
    alignment_entries: Vec<FullAlignment>,
    alignment_current_index: u32,

    reference_manager: ReferenceManager<'a, 's, 't>,

    read_table: Vec<ReadTableEntry>,
}

impl FullBAMAlignmentRegistry<'_, '_, '_> {
    pub fn new(fasta_file: &String) -> Box<FullBAMAlignmentRegistry> {

        // read in the fasta file, creating a new BAM registry from the file
        let rm = ReferenceManager::from(&fasta_file, 12, 6);
        Box::new(FullBAMAlignmentRegistry {
            extractions: HashMap::new(),
            alignments: HashMap::new(),
            alignments_reverse: HashMap::new(),
            alignment_entries: vec![],
            alignment_current_index: 0,
            reference_manager: rm,
            read_table: vec![],
        })
    }

    pub fn push_extractor(&mut self, extractor: &ExtractorTag, value: &str) -> u32 {
        if self.extractions.get(extractor).is_none() {
            self.extractions.insert(extractor.clone(), ExtractionTable::new());
        }
        let rt = self.extractions.get_mut(extractor).unwrap().add(value);
        rt
    }

    pub fn extraction_tags(&mut self, bam_entry: &Record ) -> Vec<ExtractorPair> {
        let mut extractors = Vec::new();

        bam_entry.aux_iter().for_each(|c| {
            match c {
                Ok((aux_tag, aux_enum)) => {
                    if aux_tag[0] == b'e' && aux_tag[1] >= b'0' && aux_tag[1] <= b'9' {
                        let extractor = EXTRACTORTYPE.get(String::from_utf8(aux_tag.to_vec()).unwrap().as_str()).unwrap().clone();
                        match aux_enum {
                            Aux::String(x) => {
                                extractors.push(ExtractorPair { extractor, value: self.push_extractor(&extractor, x) })//;
                            }
                            _ => {}
                        }
                    }
                }
                Err(_) => { panic!("Unable to read auxillary tags") }
            }
        });
        extractors
    }

    pub fn alignment_tags(&mut self, alignment_start: &i64, reference_sequence: &Vec<u8>, read_cigar: &CigarStringView, sequence: &Vec<u8>) -> Vec<u32> {
        let mut alignments = Vec::new();

        match alignment_start {
            x if x > &0 => {
                alignments.push(FullAlignment::Deletion(0, x.clone() as u32));
                alignments.extend(extract_read_cigar_elements(&(*x as u32), reference_sequence, sequence, read_cigar))
            }
            _ => {
                alignments.extend(extract_read_cigar_elements(&0, reference_sequence, sequence, read_cigar))
            }
        }

        alignments.extend(extract_read_cigar_elements(&(*alignment_start as u32),
                                                      reference_sequence,
                                                      sequence,
                                                      read_cigar));

        alignments.into_iter().map(|x| {
            match self.alignments.get(&x) {
                None => {
                    self.alignments.insert(x.clone(), self.alignment_current_index);
                    self.alignments_reverse.insert(self.alignment_current_index, x.clone());
                    self.alignment_entries.push(x.clone());
                    let index = self.alignment_current_index;
                    self.alignment_current_index += 1;
                    index
                }
                Some(tag) => { tag.clone() }
            }
        }).collect()
    }

    pub fn add_bam_file_entries(&mut self, bam_file: &str) {
        let mut bam = bam::Reader::from_path(bam_file).unwrap();
        let header = bam::Header::from_template(bam.header());

        // print header records to the terminal, akin to samtool
        for (key, records) in header.to_hashmap() {
            for record in records {
                println!("@{}\tSN:{}\tLN:{}", key, record["SN"], record["LN"]);
            }
        }

        let mut read_count = 0;
        // copy reverse reads to new BAM file
        for r in bam.records() {
            read_count += 1;
            if read_count % 100000 == 0 {
                println!("Read count {}",read_count);
            }
            match r {
                Ok(record) => {
                    let ref_name = record.contig().as_bytes();
                    let reference = self.reference_manager.reference_name_to_ref.get(ref_name).expect("Unable to find reference");
                    let reference_id = reference.clone() as u16;
                    let reference_sequence = &self.reference_manager.references.get(reference).unwrap().sequence_u8.to_vec();
                    let alignment_start = record.reference_start();

                    let extractor_tokens = self.extraction_tags(&record);
                    let cigar_tokens = self.alignment_tags(&alignment_start, reference_sequence, &record.cigar(), &record.seq().as_bytes());

                    self.read_table.push(ReadTableEntry {
                        name: String::from_utf8(record.name().to_vec()).unwrap(),
                        extractions: extractor_tokens,
                        alignments: cigar_tokens,
                        reference: reference_id,
                    });
                }
                Err(_) => {
                    panic!("Unable to process read")
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use rust_htslib::bam::record::CigarString;
    use super::*;

    #[test]
    fn test_breakup_nucleotide_sequences() {
        let reference = b"ACGTACGT"; // Example reference sequence
        let sequence = b"ACGTTGCt"; // Example sequence to compare
        let reference_offset = 10;

        let alignments = breakup_nucleotide_sequences(&reference.to_vec(), sequence, &reference_offset);

        let expected_alignments = vec![
            FullAlignment::Match(10, 4), // "ACG" matches
            FullAlignment::Mismatch(14, vec![b'T', b'G', b'C', b't']), // "TGCt" mismatch
        ];

        assert_eq!(alignments, expected_alignments);
    }

    #[test]
    fn test_bam_file() {
        let reference = "./test_data/bam/v10_ref_plus_castag.fa".to_string();
        let mut reader = FullBAMAlignmentRegistry::new(&reference);
        reader.add_bam_file_entries("./test_data/bam/test.bam")
    }

    #[test]
    fn test_extract_read_cigar_elements() {
        let ref_start = 0;          // ____XXX_
        let reference_sequence = b"ACGTACGT".to_vec(); // Example reference sequence
        let read_seq = b"ACGTTGCT".to_vec(); // Example read sequence

        // Constructing a simple CigarStringView
        // Assuming CigarStringView and Cigar are defined and can be created like this
        let cigar = CigarStringView::new(CigarString::from(vec![
            Cigar::Match(8),
        ]), 0);

        let alignments = extract_read_cigar_elements(&ref_start, &reference_sequence, &read_seq, &cigar);

        let expected_alignments = vec![
            FullAlignment::Match(0, 4), // "ACG" matches
            FullAlignment::Mismatch(4, vec![b'T', b'G', b'C']), // "TGC" matches (note the shifted position due to insertion)
            FullAlignment::Match(7, 1),
        ];

        assert_eq!(alignments, expected_alignments);
    }
}