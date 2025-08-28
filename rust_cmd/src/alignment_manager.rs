use crate::alignment::scoring_functions::{AffineScoring, InversionScoring};
use crate::read_strategies::read_disk_sorter::SortingReadSetContainer;
use crate::read_strategies::sequence_layout::SequenceLayout;
use crate::reference::fasta_reference::ReferenceManager;
use bstr::BString;

use noodles_sam::header::record::value::map::ReferenceSequence;
use noodles_sam::header::record::value::Map;
use noodles_sam::Header;
use noodles_util::alignment;
use std::collections::HashMap;
use std::convert::TryFrom;

use std::io::Result;
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use indexmap::IndexMap;

use crate::alignment::alignment_matrix::{
    create_scoring_record_3d, perform_3d_global_traceback, perform_affine_alignment,
    AlignmentResult, AlignmentTag, AlignmentType,
};
use noodles_sam;
use utils::read_utils::u8s;


/// A trait for writing aligned reads to various output formats.
/// The output may or may not respect all the fields of the SortingReadSetContainer.
pub trait OutputAlignmentWriter: Sync + Send {
    /// Writes a read to the output format.
    ///
    /// # Arguments
    ///
    /// * `read_set_container` - The container holding the aligned read and its metadata
    /// * `additional_tags` - Additional tags to be included in the output
    ///
    /// # Returns
    ///
    /// A `Result<()>` indicating success or failure of the write operation.
    fn write_read(
        &mut self,
        read_set_container: &SortingReadSetContainer,
        additional_tags: &HashMap<[u8; 2], String>,
    ) -> Result<()>;

    
    fn close(&mut self) -> Result<()>;
}

/// Implementation of `OutputAlignmentWriter` for BAM files.
///
/// This struct provides functionality to write aligned reads to BAM format files
/// using the noodles library for SAM/BAM operations.
pub struct BamFileAlignmentWriter<'a> {
    underlying_bam_file: Arc<Mutex<noodles_util::alignment::io::Writer>>,
    header: noodles_sam::Header,
    reference_manager: ReferenceManager<'a, 'a, 'a>,
}

unsafe impl<'a> Send for BamFileAlignmentWriter<'a> {}
unsafe impl<'a> Sync for BamFileAlignmentWriter<'a> {}

impl<'a> BamFileAlignmentWriter<'a> {
    /// Creates a new BAM file writer.
    ///
    /// # Arguments
    ///
    /// * `path` - The path where the BAM file will be written
    /// * `reference_manager` - Manager containing reference sequences and metadata
    ///
    /// # Returns
    ///
    /// A new `BamFileAlignmentWriter` instance ready to write aligned reads.
    pub fn new(
        path: &PathBuf,
        reference_manager: &ReferenceManager<'a, 'a, 'a>,
    ) -> BamFileAlignmentWriter<'a> {
        // Collect and sort by key
        let mut sorted: Vec<_> = reference_manager.references.iter().collect();
        sorted.sort_by_key(|&(k, _)| k);
        
        let reference_sequences : IndexMap<BString,Map<ReferenceSequence>> = sorted.iter().map(|(_k, v)| {
            (BString::from(v.name.clone()), 
            Map::<ReferenceSequence>::new(NonZeroUsize::try_from(v.sequence.len()).unwrap()))
        }).into_iter().collect();
        
        let header = Header::builder()
            .set_header(Default::default())
            .add_comment("Clique processed")
            .set_reference_sequences(reference_sequences)
            .build();

        let mut writer = alignment::io::writer::builder::Builder::default()
            .set_format(alignment::io::Format::Bam)
            .build_from_path(path)
            .unwrap();
        
        writer
            .write_header(&header)
            .expect("Unable to write header to output bam file");

        BamFileAlignmentWriter {
            underlying_bam_file: Arc::new(Mutex::new(writer)),
            header,
            reference_manager: reference_manager.clone(),
        }
    }
}

impl<'a> OutputAlignmentWriter for BamFileAlignmentWriter<'a> {
    /// Writes a read to the BAM file.
    ///
    /// This method converts the read container to a SAM record and writes it to the BAM file.
    /// It includes additional tags and sorting keys from the read container.
    ///
    /// # Arguments
    ///
    /// * `read_set_container` - The container holding the aligned read and its metadata
    /// * `additional_tags` - Additional tags to be included in the BAM record
    ///
    /// # Returns
    ///
    /// A `Result<()>` indicating success or failure of the write operation.
    ///
    /// # Panics
    ///
    /// Panics if the reference name cannot be found in the reference manager or
    /// if writing to the BAM file fails.
    fn write_read(
        &mut self,
        read_set_container: &SortingReadSetContainer,
        additional_tags: &HashMap<[u8; 2], String>
    ) -> Result<()> {

        let reference_record = match self
            .reference_manager
            .reference_name_to_ref
            .get(
                &read_set_container
                    .aligned_read
                    .reference_name
                    .as_bytes()
                    .to_vec(),
            ) {
            None => {

                let ref_name = u8s(&read_set_container.aligned_read.reference_name.as_bytes().to_vec());
                println!("Unable to find reference name {}" , ref_name);
                panic!("unable to find reference");
            }
            Some(x) => {x}
        };

        let mut extra_annotations = additional_tags.clone();
        read_set_container.ordered_sorting_keys.iter().for_each(|(key, value)| {
            extra_annotations.insert(
                [b'e', *key as u8],
                u8s(&value.corrected),
            );
            extra_annotations.insert(
                [b'o', *key as u8],
                u8s(&value.original),
            );
        });
            
        let samrecord =
            read_set_container
                .aligned_read
                .to_sam_record(&(*reference_record as i32), &extra_annotations, None);

        let writer = Arc::clone(&self.underlying_bam_file);

        let mut output = writer.lock().unwrap();
        
        match output
            .write_record(&self.header, &samrecord) {
                Ok(_x) => {
                    //println!("Added read!")
                },
                Err(e) => {
                    println!("samrecord {:?}",samrecord);
                    println!("Sequence: {} {:?}", u8s(&samrecord.name().unwrap().to_vec()),samrecord);
                    println!("Sequence: {} {}", u8s(&samrecord.quality_scores().clone().into()),u8s(&samrecord.sequence().clone().into()));
                    println!("error kind {} {}", e.kind().to_string(),e.to_string());
                    panic!("Unable to write record to bam file; error {}", e);
                },
            }
            
        Ok(())
    }

    fn close(&mut self) -> Result<()> {

        let writer = Arc::clone(&self.underlying_bam_file);

        let mut output = writer.lock().unwrap();

        match output.finish(&self.header) {
            Ok(_x) => {
                //println!("Added read!")
            },
            Err(e) => {
                panic!("error, unable to close output BAM writer {} {}", e.kind().to_string(),e.to_string());
            },
        }
        Ok(())
    }
}

/// Performs pairwise alignment between a reference sequence and a read sequence.
///
/// This function creates a 3D scoring matrix and performs affine gap alignment
/// followed by traceback to generate the alignment result.
///
/// # Arguments
///
/// * `reference_sequence` - The reference sequence as a vector of bytes
/// * `read_sequence` - The read sequence as a vector of bytes
/// * `read_qual` - Optional quality scores for the read sequence
/// * `scoring_function` - The affine scoring parameters for alignment
/// * `local` - Whether to perform local alignment (true) or global alignment (false)
/// * `ref_name` - Name of the reference sequence
/// * `read_name` - Name of the read sequence
/// * `_reference_manager` - Optional reference manager (currently unused)
///
/// # Returns
///
/// An `AlignmentResult` containing the alignment information including CIGAR string,
/// scores, and positions.
pub fn align_two_strings(
    reference_sequence: &Vec<u8>,
    read_sequence: &Vec<u8>,
    read_qual: Option<Vec<u8>>,
    scoring_function: &AffineScoring,
    local: bool,
    ref_name: &String,
    read_name: &String,
    _reference_manager: Option<&ReferenceManager>,
) -> AlignmentResult {
    let mut alignment_mat = create_scoring_record_3d(
        reference_sequence.len() + 1,
        read_sequence.len() + 1,
        AlignmentType::Affine,
        local,
    );

    /*match (reference_manager, ref_name) {
        (Some(x), Some(y)) => {
            panic!("for now");
        }

        _ => {*/
            perform_affine_alignment(
                &mut alignment_mat,
                reference_sequence,
                read_sequence,
                scoring_function,
            );

            perform_3d_global_traceback(
                &mut alignment_mat,
                None,
                reference_sequence,
                read_sequence,
                ref_name,
                read_name,
                read_qual,
                None,
            )
       // }
    //}
}

/// Configuration parameters for sequence alignment operations.
///
/// This struct holds all the necessary parameters and settings required
/// for performing sequence alignments, including file paths, scoring functions,
/// and filtering criteria.
#[allow(dead_code)]
struct AlignmentParameters<'a> {
    read_structure: &'a SequenceLayout,
    rm: &'a ReferenceManager<'a, 'a, 'a>,
    minimum_read_to_reference_length_ratio: f64,
    maximum_read_to_reference_length_ratio: f64,
    minimum_read_length: usize,
    find_inversions: bool,
    read1: PathBuf,
    read2: Option<PathBuf>,
    index1: Option<PathBuf>,
    index2: Option<PathBuf>,
    affine_alignment_scores: AffineScoring,
    inversion_alignment_scores: InversionScoring,
}

impl<'a> AlignmentParameters<'a> {

    /// Creates default alignment parameters for affine gap alignment.
    ///
    /// This constructor sets up standard parameters suitable for most DNA sequence
    /// alignment tasks using affine gap penalties.
    ///
    /// # Arguments
    ///
    /// * `read1` - Path to the first read file
    /// * `read2` - Optional path to the second read file (for paired-end)
    /// * `index1` - Optional path to the first index file
    /// * `index2` - Optional path to the second index file
    /// * `read_structure` - Layout describing the structure of the sequencing reads
    /// * `reference_manager` - Manager containing reference sequences
    ///
    /// # Returns
    ///
    /// An `AlignmentParameters` instance with default settings for affine alignment.
    #[allow(dead_code)]
    pub fn default_affine_alignment(
        read1: PathBuf,
        read2: Option<PathBuf>,
        index1: Option<PathBuf>,
        index2: Option<PathBuf>,
        read_structure: &'a SequenceLayout,
        reference_manager: &'a ReferenceManager,
    ) -> AlignmentParameters<'a> {
        let affine_scores = AffineScoring::default_dna();
        let inversion_scores = InversionScoring::default();

        AlignmentParameters {
            read_structure,
            rm: reference_manager,
            minimum_read_to_reference_length_ratio: 0.1,
            maximum_read_to_reference_length_ratio: 1.5,
            minimum_read_length: 20,
            find_inversions: false,
            read1,
            read2,
            index1,
            index2,
            affine_alignment_scores: affine_scores,
            inversion_alignment_scores: inversion_scores,
        }
    }
}

/// Simplifies a CIGAR string represented as a vector of `AlignmentTag`.
///
/// This function consolidates consecutive alignment tags of the same type by summing their counts.
/// It is particularly useful for compressing a CIGAR string by merging adjacent tags that
/// represent the same kind of operation.
///
/// # Arguments
///
/// * `cigar_tokens` - A reference to a vector of `AlignmentTag` enum items representing the CIGAR operations.
///
/// # Returns
///
/// A new vector of `AlignmentTag` items where consecutive tags of the same type have been combined.
///
/// # Panics
///
/// The function panics if it encounters two consecutive inversion open tags or two consecutive inversion close tags,
/// as this would represent an invalid state in a CIGAR string.
///
/// # Examples
///
/// ```
/// let cigar_tokens = vec![
///     AlignmentTag::MatchMismatch(3),
///     AlignmentTag::MatchMismatch(5),
///     AlignmentTag::Ins(2),
///     AlignmentTag::Del(4),
///     AlignmentTag::Del(1),
/// ];
/// let simplified_cigar = simplify_cigar_string(&cigar_tokens);
/// assert_eq!(simplified_cigar, vec![
///     AlignmentTag::MatchMismatch(8),
///     AlignmentTag::Ins(2),
///     AlignmentTag::Del(5),
/// ]);
/// ```
///
/// # Notes
///
/// The `AlignmentTag` is assumed to be an enum that can represent different CIGAR operations such as
/// matches/mismatches (`MatchMismatch`), insertions (`Ins`), deletions (`Del`), and inversions (`InversionOpen` and `InversionClose`).
///
pub fn simplify_cigar_string(cigar_tokens: &Vec<AlignmentTag>) -> Vec<AlignmentTag> {
    let mut new_cigar = Vec::new();

    let mut last_token: Option<AlignmentTag> = None; // zero length, so combining won't affect the final cigar string

    cigar_tokens
        .iter()
        .for_each(|token| match (&last_token, token) {
            (None, _) => last_token = Some(token.clone()),
            (Some(AlignmentTag::InversionOpen), AlignmentTag::InversionOpen) => {
                panic!("Cannot have two inversion open tags in a row");
            }
            (Some(AlignmentTag::InversionClose), AlignmentTag::InversionClose) => {
                panic!("Cannot have two inversion closed tags in a row");
            }
            (
                Some(AlignmentTag::MatchMismatch(last_count)),
                AlignmentTag::MatchMismatch(this_count),
            ) => {
                last_token = Some(AlignmentTag::MatchMismatch(last_count + this_count));
            }
            (Some(AlignmentTag::Del(last_count)), AlignmentTag::Del(this_count)) => {
                last_token = Some(AlignmentTag::Del(last_count + this_count));
            }
            (Some(AlignmentTag::Ins(last_count)), AlignmentTag::Ins(this_count)) => {
                last_token = Some(AlignmentTag::Ins(last_count + this_count));
            }
            (Some(x), y) => {
                new_cigar.push(x.clone());
                last_token = Some(y.clone());
            }
        });

    if let Some(x) = last_token {
        new_cigar.push(x);
    }
    new_cigar
}
