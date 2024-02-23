use crate::alignment::scoring_functions::{AffineScoring, InversionScoring};
use crate::read_strategies::read_disk_sorter::SortingReadSetContainer;
use crate::read_strategies::sequence_layout::SequenceLayout;
use crate::reference::fasta_reference::ReferenceManager;
use bstr::BString;
use ndarray::Ix3;
use noodles_sam::header::record::value::map::ReferenceSequence;
use noodles_sam::header::record::value::Map;
use noodles_sam::Header;
use noodles_util::alignment;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::Result;
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use crate::alignment::alignment_matrix::{
    create_scoring_record_3d, perform_3d_global_traceback, perform_affine_alignment,
    AlignmentResult, AlignmentTag, AlignmentType,
};
use crate::alignment::fasta_bit_encoding::FastaBase;
use noodles_bam;
use noodles_bgzf::Writer;
use noodles_sam;

/// something that writes aligned reads. The output may or may not respect all the fields
/// of the SortingReadSetContainer
pub trait OutputAlignmentWriter {
    fn write_read(
        &mut self,
        read_set_container: &SortingReadSetContainer,
        additional_tags: &HashMap<(u8, u8), String>,
    ) -> Result<()>;
}

/// implement a OutputAlignmentWriter for BAM files
pub struct BamFileAlignmentWriter<'a> {
    underlying_bam_file: Arc<Mutex<noodles_bam::io::Writer<noodles_bgzf::Writer<Writer<File>>>>>,
    header: noodles_sam::Header,
    reference_manager: ReferenceManager<'a, 'a, 'a>,
}

impl<'a> BamFileAlignmentWriter<'a> {
    pub fn new(path: &PathBuf, reference_manager: &ReferenceManager<'a,'a,'a>) -> BamFileAlignmentWriter<'a> {
        
           let reference_sequences = [(
    BString::from("sq0"),
    Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13).unwrap()),
    )].into_iter().collect();
    /*
        let reference_sequences:  IndexMap<BString, Map<ReferenceSequence>> = IndexMap::new();
        reference_manager
            .references
            .iter()
            .for_each(|(id, reference)| {
                reference_sequences.insert(BString::from(String::from_utf8(reference.name.clone()).unwrap()),
                    Map::<ReferenceSequence>::new(NonZeroUsize::try_from(*id).unwrap()));
            });
*/
        let header = Header::builder()
            .set_header(Default::default())
            .add_comment("an example BAM written by noodles-bam")
            .set_reference_sequences(reference_sequences)
            .build();

        let bgzf = Writer::new(File::open(path).expect("Unable to create output bam file"));

        let mut writer = noodles_bam::io::Writer::new(bgzf);

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
    fn write_read(
        &mut self,
        read_set_container: &SortingReadSetContainer,
        additional_tags: &HashMap<(u8, u8), String>,
    ) -> Result<()> {
        let writer = Arc::clone(&self.underlying_bam_file);

        let reference_record = self
            .reference_manager
            .reference_name_to_ref
            .get(
                &read_set_container
                    .aligned_read
                    .reference_name
                    .as_bytes()
                    .to_vec(),
            )
            .unwrap();
        //let reference_record = self.reference_manager.references.get(reference_record).unwrap();
        let samrecord =
            read_set_container
                .aligned_read
                .to_sam_record(&(*reference_record as i32), &true, None);

        let mut output = writer.lock().unwrap();
        let mut bam_writer = alignment::io::writer::builder::Builder::default()
            .set_format(alignment::io::Format::Bam)
            .build_from_writer(std::io::stdout().lock())
            .unwrap();
        bam_writer
            .write_record(&self.header, &samrecord)
            .expect("Unable to write read to output bam file");
        Ok(())
    }
}

pub fn align_two_strings(
    read1_seq: &Vec<FastaBase>,
    rev_comp_read2: &Vec<FastaBase>,
    scoring_function: &AffineScoring,
    local: bool,
    ref_name: Option<&Vec<u8>>,
    reference_manager: Option<&ReferenceManager>,
) -> AlignmentResult {
    let mut alignment_mat = create_scoring_record_3d(
        read1_seq.len() + 1,
        rev_comp_read2.len() + 1,
        AlignmentType::Affine,
        local,
    );

    match (reference_manager, ref_name) {
        (Some(x), Some(y)) => {
            panic!("for now");
        }

        _ => {
            perform_affine_alignment(
                &mut alignment_mat,
                read1_seq,
                rev_comp_read2,
                scoring_function,
            );

            perform_3d_global_traceback(
                &mut alignment_mat,
                None,
                read1_seq,
                rev_comp_read2,
                &String::from_utf8(ref_name.unwrap_or(&"reference".as_bytes().to_vec()).clone())
                    .unwrap(),
                &String::from("read_name"),
                None,
            )
        }
    }
}

/*
impl BamFileAlignmentWriter {
    pub fn new(output_path: &PathBuf) -> BamFileAlignmentWriter {
        let mut writer = bam::io::Writer::new(output_path.to_str().expect("Unable to access output bam file path"));

        let header = sam::Header::builder()
            .set_header(Default::default())
            .add_program("noodles-bam", Map::<Program>::default())
            .add_comment("a BAM written by clique")
            .build();

        writer.write_header(&header).expect("Unable to write bam file header");

        BamFileAlignmentWriter{
            underlying_bam_file: Arc::new(Mutex::new(writer)),
            header,
        }
    }
}*/

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
    pub fn default_affine_alignment(
        read1: PathBuf,
        read2: Option<PathBuf>,
        index1: Option<PathBuf>,
        index2: Option<PathBuf>,
        read_structure: &'a SequenceLayout,
        reference_manager: &'a ReferenceManager,
    ) -> AlignmentParameters<'a> {
        let affine_scores = AffineScoring::default();
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
/*
struct AlignmentResults {

}

trait AlignmentFilter {

}

/// Takes a recipe (alignment parameters) and an output sink and uses threads to
/// align reads
struct AlignmentEngine {

}

impl AlignmentEngine {
    pub fn align_sequences(parameters: &AlignmentParameters, output_writer: &OutputAlignmentWriter, alignment_approach: &AlignmentType, threads: &usize) -> AlignmentResults {

        let read_iterator = MergedReadSequence::new(ReadIterator::new(PathBuf::from(&read1),
                                                                      Some(PathBuf::from(&read2)),
                                                                      Some(PathBuf::from(&index1)),
                                                                      Some(PathBuf::from(&index2))), read_structure);
        // setup the alignment matrix with a thread pool
        let alignment_mat: Alignment<Ix3> = create_scoring_record_3d((parameters.rm.longest_ref + 1),
                                                                     (parameters.maximum_read_to_reference_length_ratio * parameters.rm.longest_ref).ceil() as usize + 1,
                                                                     *alignment_approach,
                                                                     false);

        lazy_static! {static ref STORE_CLONES: Mutex<Vec<SharedStore>> = Mutex::new(Vec::new());}
        thread_local!(static STORE: SharedStore = Arc::new(Mutex::new(None)));



        AlignmentResults{}
    }
}
*/

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
