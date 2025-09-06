use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::slice;

use crate::alignment::alignment_matrix::{
    create_scoring_record_3d, perform_3d_global_traceback, perform_affine_alignment,
    perform_affine_alignment_bandwidth, Alignment, AlignmentResult, AlignmentTag, AlignmentType,
};
use crate::alignment::scoring_functions::{AffineScoring, InversionScoring};
use crate::linked_alignment::{
    align_string_with_anchors, find_greedy_non_overlapping_segments, orient_by_longest_segment,
};
use crate::rayon::iter::ParallelBridge;
use crate::rayon::iter::ParallelIterator;
use crate::read_strategies::read_set::ReadIterator;
use crate::reference::fasta_reference::ReferenceManager;
use ndarray::Ix3;
use std::time::Instant;
use bio::alignment::AlignmentOperation;
use crate::merger::{MergedReadSequence, UnifiedRead};

use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation::*;
use bio::scores::blosum62;

use crate::read_strategies::sequence_layout::SequenceLayout;

use crate::alignment_manager::BamFileAlignmentWriter;
use crate::alignment_manager::OutputAlignmentWriter;
use crate::read_strategies::read_disk_sorter::SortingReadSetContainer;
use consensus::consensus_builders::get_reference_alignment_rate;
use extractor::extract_tagged_sequences;
use itertools::Itertools;
use sigalign::results::AlignmentOperations;
use reference::fasta_reference::Reference;
use utils::read_utils::{reverse_complement, u8s};
use ::{Aligner as RustAligner, FASTA_UNSET};

/// Include the generated bindings into a separate module.
#[allow(non_upper_case_globals)]
#[allow(non_snake_case)]
#[allow(non_camel_case_types)]
#[allow(unused)]
mod wfa {
    include!(concat!(env!("OUT_DIR"), "/bindings_wfa.rs"));
}

/// Compute the affine alignment score between `a` and `b` with the given substitution,
/// gap-open, and gap-extend penalties.
#[allow(dead_code)]
fn wfa_alignment(
    read: &[u8],
    reference: &[u8],
    mismatch: &i32,
    gap_open: &i32,
    gap_extend: &i32
) -> Vec<AlignmentOperation> {
    let score = |a: u8, b: u8| if a == b || a == b'N' { 1i32 } else { -1i32 };
    // gap open score: -5, gap extension score: -1
    let mut aligner = Aligner::with_capacity(read.len(), reference.len(), -5, -1, &score);
    let alignment = aligner.global(reference,read);
    // x is global (target sequence) and y is local (reference sequence)
    alignment.operations
}

pub fn align_reads(
    read_structure: &SequenceLayout,
    rm: &ReferenceManager,
    output: &Path,
    max_reference_multiplier: &usize,
    min_read_length: &usize,
    read1: &String,
    read2: &String,
    index1: &String,
    index2: &String,
    threads: &usize,
    _aligner: &RustAligner,
) {
    let read_iterator = ReadIterator::new(
        PathBuf::from(&read1),
        Some(PathBuf::from(&read2)),
        Some(PathBuf::from(&index1)),
        Some(PathBuf::from(&index2)),
    );

    let read_iterator = MergedReadSequence::new(read_iterator, read_structure);

    let writer = BamFileAlignmentWriter::new(&PathBuf::from(output), &rm);

    let output = Arc::new(Mutex::new(writer));

    // setup our thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(*threads)
        .build_global()
        .unwrap();

    let my_score = InversionScoring {
        match_score: 9.0,
        mismatch_score: -21.0,
        gap_open: -25.0,
        gap_extend: -1.0,
        inversion_penalty: -40.0,
        min_inversion_length: 20,
    };

    let my_aff_score = AffineScoring {
        match_score: 10.0,
        mismatch_score: -9.0,
        special_character_score: 9.0,
        gap_open: -20.0,
        gap_extend: -2.0,
        final_gap_multiplier: 1.0,
    };
    let start = Instant::now();
    let read_count = Arc::new(Mutex::new(0)); // we rely on this Arc for output file access control

    type SharedStore = Arc<Mutex<Option<Alignment<Ix3>>>>;

    lazy_static! {
        static ref STORE_CLONES: Mutex<Vec<SharedStore>> = Mutex::new(Vec::new());
    }
    thread_local!(static STORE: SharedStore = Arc::new(Mutex::new(None)));

    let max_read_size = (rm.longest_ref + 1) * max_reference_multiplier;
    info!(
        "Longest reference found: {}, max read size set at {}",
        rm.longest_ref, max_read_size
    );

    let alignment_mat: Alignment<Ix3> = create_scoring_record_3d(
        rm.longest_ref + 1,
        max_read_size,
        AlignmentType::Affine,
        false,
    );

    read_iterator.par_bridge().for_each(|mut xx: UnifiedRead| {
        STORE.with(|arc_mtx| {
            let mut local_alignment = arc_mtx.lock().unwrap();
            if local_alignment.is_none() {
                *local_alignment = Some(alignment_mat.clone());
                STORE_CLONES.lock().unwrap().push(arc_mtx.clone());
            }

            let name = &String::from_utf8(xx.name().clone()).unwrap();

            let seq_len = &xx.seq().len();
            let qual = Some(xx.quals.as_ref().unwrap().clone());
            if seq_len < &max_read_size {
                let aligned = align_to_reference_choices(
                    name,
                    xx.seq(),
                    qual,
                    rm,
                    &true,
                    read_structure,
                    local_alignment.as_mut().unwrap(),
                    &my_aff_score,
                    &my_score,
                    &false,
                    *max_reference_multiplier as f64,
                    *min_read_length,
                    seq_len,
                );

                match aligned {
                    None => {
                        // TODO: we should track this and provide a final summary
                        debug!("Unable to create alignment for read {}", name);
                    }
                    Some(alignment_obj) => {
                        let results = alignment_obj.alignment;

                        let _orig_ref_seq = alignment_obj.ref_sequence;
                        match results {
                            None => {
                                // TODO: we should track this and provide a final summary

                                debug!("Unable to create alignment for read {}", name);
                            }
                            Some(aln) => {
                                let mut read_count = read_count.lock().unwrap();
                                *read_count += 1;
                                if *read_count % 1000000 == 0 {
                                    let duration = start.elapsed();
                                    info!(
                                        "Time elapsed in aligning reads ({:?}) is: {:?}",
                                        read_count, duration
                                    );
                                }
                                assert_eq!(aln.reference_aligned.len(), aln.read_aligned.len());

                                let read = SortingReadSetContainer::empty_tags(aln);

                                let extracted_tags = extract_tagged_sequences(
                                    &read.aligned_read.read_aligned,
                                    &read.aligned_read.reference_aligned,
                                );

                                let mut added_tags: HashMap<[u8; 2], String> = HashMap::new();

                                let structure = read_structure
                                    .references
                                    .get(&read.aligned_read.reference_name)
                                    .unwrap();

                                extracted_tags.iter().for_each(|(x, y)| {
                                    structure.umi_configurations.iter().for_each(|xi| {
                                        if xi.1.symbol as u8 == *x {
                                            added_tags.insert([b'e', xi.1.symbol as u8], y.clone());
                                        }
                                    })
                                });

                                added_tags.insert([b'r', b'c'], 1.to_string());

                                added_tags
                                    .insert([b'a', b'r'], read.aligned_read.read_name.clone());
                                added_tags.insert(
                                    [b'r', b'm'],
                                    get_reference_alignment_rate(
                                        &read.aligned_read.reference_aligned,
                                        &read.aligned_read.read_aligned,
                                    )
                                    .to_string(),
                                );
                                added_tags
                                    .insert([b'a', b's'], read.aligned_read.score.to_string());

                                let output = Arc::clone(&output);
                                let arc_writer = output.clone();
                                let mut arc_writer = arc_writer
                                    .lock()
                                    .expect("Unable to access multi-threaded writer");
                                arc_writer
                                    .write_read(&read, &added_tags)
                                    .expect("Unable to write a read to the arc writer (LOC1)");
                            }
                        }
                    }
                }
            } else {
                warn!(
                    "Dropped read {} is it's length {} exceeds 2x the reference length {}",
                    String::from_utf8(xx.name().clone()).unwrap(),
                    xx.seq().len(),
                    max_read_size
                );
            }
        });
    });

    let output = Arc::clone(&output);
    let arc_writer = output.clone();
    let mut arc_writer = arc_writer
        .lock()
        .expect("Unable to access multi-threaded writer");
    arc_writer.close().unwrap();
}

#[allow(dead_code)]
pub fn align_two_strings(
    read1_name: &String,
    read1_seq: &[u8],
    sequence_2_seq: &[u8],
    scoring_function: &AffineScoring,
    local: bool,
    ref_name: &String,
    reference_manager: Option<&ReferenceManager>,
) -> AlignmentResult {
    let mut alignment_mat = create_scoring_record_3d(
        read1_seq.len() + 1,
        sequence_2_seq.len() + 1,
        AlignmentType::Affine,
        local,
    );

    match (reference_manager, ref_name) {
        (Some(x), y) => {
            let ref_id = x.reference_name_to_ref.get(y.as_bytes()).unwrap();
            let shared_segments = &x.references.get(ref_id).unwrap().suffix_table;
            let ref_name =
                String::from_utf8(x.references.get(ref_id).unwrap().name.clone()).unwrap();

            let ref_seq = read1_seq;
            let read_seq = sequence_2_seq;

            let shared_segs =
                find_greedy_non_overlapping_segments(&ref_seq, &read_seq, shared_segments);

            align_string_with_anchors(
                read1_name,
                &ref_name,
                read1_seq,
                sequence_2_seq,
                &shared_segs,
                None,
                scoring_function,
                &mut alignment_mat,
            )
        }

        _ => {
            perform_affine_alignment(
                &mut alignment_mat,
                read1_seq,
                sequence_2_seq,
                scoring_function,
            );

            perform_3d_global_traceback(
                &mut alignment_mat,
                None,
                read1_seq,
                sequence_2_seq,
                read1_name,
                ref_name,
                None,
                None,
            )
        }
    }
}

/// Aligns two sequences based on a provided matrix, using either specialized or affine alignment methods.
///
/// This function aligns two DNA or RNA sequences (`read1_seq` and `read2_seq`) and utilizes an
/// affine scoring function. It can optionally use a reference manager and a reference name to
/// perform specialized alignment using shared segments. If these are not provided, it defaults
/// to affine alignment.
///
/// # Arguments
///
/// * `read1_seq` - A reference to the vector containing the first sequence (as `FastaBase` elements) to be aligned.
/// * `read2_seq` - A reference to the vector containing the second sequence to be aligned.
/// * `scoring_function` - A reference to an object implementing the `AffineScoringFunction` trait for scoring alignments.
/// * `ref_name` - An optional reference to a vector of bytes representing the name of the reference sequence.
/// * `reference_manager` - An optional reference to a `ReferenceManager` object managing reference sequences.
/// * `alignment_mat` - A mutable reference to an `Alignment` object, a 3D matrix (`Ix3`) used for storing alignment scores.
///
/// # Returns
///
/// An `AlignmentResult` that may contain details such as the alignment score, the aligned sequences, and other relevant information.
///
/// # Behavior
///
/// The function operates in two modes depending on the availability of `reference_manager` and `ref_name`:
///
/// 1. **With Reference Manager and Reference Name**:
///    - Retrieves a reference ID and shared segments from the reference manager using the provided `ref_name`.
///    - Converts both sequences to `Vec<u8>`.
///    - Finds greedy non-overlapping segments between the sequences using the shared segments.
///    - Performs alignment of `read2_seq` with `read1_seq` using these segments and the scoring function.
///
/// 2. **Without Reference Manager and Reference Name**:
///    - Performs an affine alignment between `read1_seq` and `read2_seq` using the provided scoring function.
///    - Executes a 3D global traceback on the alignment matrix to produce the final alignment result.
///
/// # Panics
///
/// The function may panic if invalid references or indices are encountered, particularly when unwrapping
/// `Option` types or retrieving references from the reference manager.
///
/// # Example
///
/// ```
/// let read1 = vec![/* ... FastaBase elements ... */];
/// let read2 = vec![/* ... FastaBase elements ... */];
/// let scoring_function = /* ... implementation of AffineScoringFunction ... */;
/// let ref_name = Some(&vec![/* reference name as Vec<u8> */]);
/// let reference_manager = Some(&/* ReferenceManager instance */);
/// let mut alignment_matrix = Alignment::new(/* ... dimensions ... */);
///
/// let alignment_result = align_two_strings_passed_matrix(
///     &read1,
///     &read2,
///     &scoring_function,
///     ref_name,
///     reference_manager,
///     &mut alignment_matrix
/// );
///
/// // Use `alignment_result` here
/// ```
pub fn align_two_strings_passed_matrix(
    read1_name: &String,
    read2_name: &String,
    read1_seq: &[u8],
    read2_seq: &[u8],
    qual_sequence: Option<Vec<u8>>,
    scoring_function: &AffineScoring,
    alignment_mat: &mut Alignment<Ix3>,
    max_indel: &usize,
) -> AlignmentResult {
    /*match (reference_manager, ref_name) {
    (Some(x), Some(y)) => {
        //let ref_id = x.reference_name_to_ref.get(y).unwrap();
        //let shared_segments = &x.references.get(ref_id).unwrap().suffix_table;

        //let ref_seq = FastaBase::to_vec_u8(read1_seq);
        //let read_seq = FastaBase::to_vec_u8(read2_seq);

        // TODO fix this
        //let shared_segs = find_greedy_non_overlapping_segments(
        //    &read_seq,
        //    &ref_seq,
        //    shared_segments);

        perform_affine_alignment(
            alignment_mat,
            read1_seq,
            read2_seq,
            scoring_function);

        perform_3d_global_traceback(
            alignment_mat,
            None,
            read1_seq,
            read2_seq,
            None)

        /*align_string_with_anchors(read2_seq,
                                  read1_seq,
                                  &shared_segs,
                                  None,
                                  scoring_function,
                                  alignment_mat)*/
    }

    _ => {*/
    perform_affine_alignment_bandwidth(
        alignment_mat,
        read1_seq,
        read2_seq,
        scoring_function,
        &max_indel,
    );

    perform_3d_global_traceback(
        alignment_mat,
        None,
        read1_seq,
        read2_seq,
        read1_name,
        read2_name,
        qual_sequence,
        None,
    )
    //}
    //}
}

#[allow(dead_code)]
#[derive(Clone)]
pub struct AlignmentWithRef {
    alignment: Option<AlignmentResult>,
    ref_name: Vec<u8>,
    ref_sequence: Vec<u8>,
}

/// Aligns two DNA or RNA sequences using affine alignment or a specialized alignment with anchors.
///
/// This function aligns two sequences represented by `Vec<FastaBase>` using either affine alignment
/// or a specialized alignment with anchors, depending on the availability of a reference manager
/// and a reference name. It employs a scoring function to evaluate the alignments.
///
/// # Arguments
///
/// * `read1_seq` - A reference to the first sequence to be aligned, each element being of type `FastaBase`.
/// * `read2_seq` - A reference to the second sequence to be aligned.
/// * `scoring_function` - A reference to an object implementing the `AffineScoringFunction` trait,
///    used for scoring alignments.
/// * `ref_name` - An optional reference to a byte vector representing the name of the reference sequence.
/// * `reference_manager` - An optional reference to a `ReferenceManager` object, which manages reference sequences.
/// * `alignment_mat` - A mutable reference to an `Alignment` object (3D matrix) used for storing alignment scores.
///
/// # Returns
///
/// An `AlignmentResult` containing the details of the alignment, which may include the alignment score,
/// aligned sequences, and other relevant information.
///
/// # Behavior
///
/// If both `reference_manager` and `ref_name` are provided, the function:
/// - Retrieves a reference ID and shared segments from the reference manager.
/// - Converts input sequences to `Vec<u8>`.
/// - Finds greedy non-overlapping segments between the read and reference sequences.
/// - Aligns `read2_seq` with `read1_seq` using these segments, the scoring function, and the alignment matrix.
///
/// If either `reference_manager` or `ref_name` is `None`, the function defaults to:
/// - Performing an affine alignment between `read1_seq` and `read2_seq`.
/// - Conducting a 3D global traceback on the alignment matrix for the final alignment result.
///
/// # Panics
///
/// The function may panic if invalid references or indices are encountered, particularly when unwrapping
/// `Option` types or retrieving references from the reference manager. Proper exception handling should be
/// considered in the calling code.
///
/// # Example
///
/// ```
/// let read1 = vec![/* ... FastaBase elements ... */];
/// let read2 = vec![/* ... FastaBase elements ... */];
/// let scoring_function = /* ... implementation of AffineScoringFunction ... */;
/// let ref_name = Some(&vec![/* reference name as Vec<u8> */]);
/// let reference_manager = Some(&/* ReferenceManager instance */);
/// let mut alignment_matrix = Alignment::new(/* ... dimensions ... */);
///
/// let result = align_two_strings_passed_matrix(
///     &read1,
///     &read2,
///     &scoring_function,
///     ref_name,
///     reference_manager,
///     &mut alignment_matrix
/// );
///
/// // Use `result` here
/// ```
///
pub fn align_to_reference_choices(
    read_name: &String,
    read: &Vec<u8>,
    qual_sequence: Option<Vec<u8>>,
    rm: &ReferenceManager,
    fast_lookup: &bool,
    read_structure: &SequenceLayout,
    alignment_mat: &mut Alignment<Ix3>,
    my_aff_score: &AffineScoring,
    _my_score: &InversionScoring,
    _use_inversions: &bool,
    _max_reference_multiplier: f64,
    _min_read_length: usize,
    _max_indel: &usize,
) -> Option<AlignmentWithRef> {
    match rm.references.len() {
        0 => {
            // TODO: we should track this and provide a final summary
            warn!(
                "Unable to align read {} as it has no candidate references",
                u8s(read)
            );
            None
        }
        1 => {
            let ref_base = &rm.references.get(&0).unwrap();
            let ref_name = String::from_utf8(ref_base.name.clone()).unwrap();
            let forward_oriented_seq = if !read_structure.known_strand {
                let orientation =
                    orient_by_longest_segment(&read, &ref_base.sequence, &ref_base.suffix_table).0;
                if orientation {
                    read.clone()
                } else {
                    reverse_complement(&read)
                }
            } else {
                read.clone()
            };

            let alignment = wfa_alignment(&ref_base.sequence, &forward_oriented_seq, &4, &10, &1);
            //println!("{}", alignment);

            let alignment = cigar_to_alignment(
                &ref_base.sequence,
                &forward_oriented_seq,
                &alignment,
            );

            let result = AlignmentResult {
                reference_name: ref_name.clone(),
                read_name: read_name.clone(),
                reference_aligned: alignment.0,
                read_aligned: alignment.1,
                read_quals: qual_sequence,
                cigar_string: alignment.2,
                path: vec!(),
                score: 0.0,
                reference_start: 0,
                read_start: 0,
                bounding_box: None,
            };


            /*
            let aln = align_two_strings_passed_matrix(
                &ref_name,
                read_name,
                &ref_base.sequence,
                &forward_oriented_seq,
                qual_sequence,
                my_aff_score,
                alignment_mat,
                &100,
            );
*/
            Some(AlignmentWithRef {
                alignment: Some(result),
                ref_name: ref_base.name.clone(),
                ref_sequence: ref_base.sequence.clone(),
            })
        }
        x if x > 1 => {
            if *fast_lookup {
                quick_alignment_search(
                    read_name,
                    read,
                    qual_sequence,
                    &rm,
                    alignment_mat,
                    my_aff_score,
                    &0.90,
                ) // TODO: parameterize this
            } else {
                exhaustive_alignment_search(
                    read_name,
                    read,
                    qual_sequence,
                    &rm,
                    alignment_mat,
                    my_aff_score,
                    None,
                )
            }
        }
        x => {
            panic!("we dont know what to do with a reference count of {}", x)
        }
    }
}

/// Performs a quick alignment search for a given read against a reference manager.
///
/// This function takes a read sequence and aligns it to the most likely reference sequence
/// in a reference manager using k-mers. It supports both affine scoring and inversion scoring,
/// and can optionally include inversions in the alignment process.
///
/// # Arguments
///
/// * `read` - A reference to the sequence to be aligned, represented as `Vec<FastaBase>`.
/// * `rm` - A reference to the `ReferenceManager` containing reference sequences and k-mer information.
/// * `read_structure` - A reference to the `SequenceLayoutDesign` that describes the layout of the read sequence.
/// * `alignment_mat` - A mutable reference to an `Alignment` object (3D matrix) used for storing alignment scores.
/// * `my_aff_score` - A reference to an `AffineScoring` object for scoring alignments.
/// * `my_score` - A reference to an `InversionScoring` object for scoring inversions, if used.
/// * `use_inversions` - A boolean reference indicating whether inversions should be considered in the alignment.
/// * `max_reference_multiplier` - A floating-point value representing a multiplier to determine the maximum reference size.
/// * `min_read_length` - The minimum length of the read sequence for performing the alignment.
///
/// # Returns
///
/// Returns `Option<(Option<AlignmentResult>, Vec<u8>, Vec<u8>)>`. The tuple inside the `Option` contains:
/// - An `Option<AlignmentResult>` which is `None` if no alignment is found, or contains the alignment result.
/// - A `Vec<u8>` representing the sequence of the matched reference.
/// - A `Vec<u8>` representing the name of the matched reference.
///
/// # Behavior
///
/// The function converts the read sequence into a `Vec<u8>` and generates k-mers. It then finds
/// the reference sequence in the reference manager with the highest count of matching k-mers.
/// If such a reference is found, it performs an alignment using `align_two_strings_passed_matrix`
/// and returns the alignment result along with the sequence and name of the matched reference.
/// If no matching reference is found, it returns `None`.
///
/// # Example
///
/// ```
/// let read = vec![/* ... FastaBase elements ... */];
/// let reference_manager = /* ... ReferenceManager instance ... */;
/// let read_structure = /* ... SequenceLayoutDesign instance ... */;
/// let mut alignment_matrix = Alignment::new(/* ... dimensions ... */);
/// let affine_scoring = /* ... AffineScoring instance ... */;
/// let inversion_scoring = /* ... InversionScoring instance ... */;
/// let use_inversions = &true;
/// let max_reference_multiplier = 1.5;
/// let min_read_length = 100;
///
/// let result = quick_alignment_search(
///     &read,
///     &reference_manager,
///     &read_structure,
///     &mut alignment_matrix,
///     &affine_scoring,
///     &inversion_scoring,
///     use_inversions,
///     max_reference_multiplier,
///     min_read_length
/// );
///
/// // Process `result` here
/// ```
fn quick_alignment_search(
    read_name: &String,
    read: &Vec<u8>,
    qual_sequence: Option<Vec<u8>>,
    rm: &ReferenceManager,
    alignment_mat: &mut Alignment<Ix3>,
    my_aff_score: &AffineScoring,
    match_threshold: &f64,
) -> Option<AlignmentWithRef> {
    let read_u8 = read;
    let read_kmers = ReferenceManager::sequence_to_kmers(&read_u8, &rm.kmer_size, &rm.kmer_skip);

    let max_ref = read_kmers
        .iter()
        .map(|(kmer, _c)| rm.unique_kmers.kmer_to_reference.get(kmer))
        .flatten()
        .counts();

    let count: f64 = max_ref.iter().map(|(_x, y)| y).sum::<usize>() as f64;
    let proportions = max_ref
        .iter()
        .map(|(x, y)| (x, *y as f64 / count))
        .collect::<Vec<(&&Reference, f64)>>();
    let max_ref = proportions.iter().max_by(|x, y| x.1.total_cmp(&y.1));

    match max_ref {
        None => {
            info!("No reference found; moving to exhaustive_alignment_search");
            exhaustive_alignment_search(
                read_name,
                read,
                qual_sequence,
                rm,
                alignment_mat,
                my_aff_score,
                None,
            )
        }
        Some(x) => match x.1 {
            prop if prop > *match_threshold => {
                let ref_name = String::from_utf8(x.0.name.clone()).unwrap();
                Some(AlignmentWithRef {
                    alignment: Some(align_two_strings_passed_matrix(
                        &ref_name,
                        read_name,
                        &x.0.sequence,
                        read,
                        qual_sequence,
                        my_aff_score,
                        alignment_mat,
                        &read.len(),
                    )),
                    ref_name: x.0.name.clone(),
                    ref_sequence: x.0.sequence.clone(),
                })
            }
            _ => {
                let ref_names: HashSet<Vec<u8>> = proportions
                    .iter()
                    .map(|(x, _y)| x.name.clone())
                    .into_iter()
                    .collect();
                exhaustive_alignment_search(
                    read_name,
                    read,
                    qual_sequence,
                    rm,
                    alignment_mat,
                    my_aff_score,
                    Some(ref_names),
                )
            }
        },
    }
}

fn exhaustive_alignment_search(
    read_name: &String,
    read: &Vec<u8>,
    qual_sequence: Option<Vec<u8>>,
    rm: &ReferenceManager,
    alignment_mat: &mut Alignment<Ix3>,
    my_aff_score: &AffineScoring,
    reference_subset: Option<HashSet<Vec<u8>>>,
) -> Option<AlignmentWithRef> {
    let references = &rm.references;

    let ranked_alignments = references
        .iter()
        .map(|reference| {
            if reference_subset.is_none()
                || reference_subset
                    .as_ref()
                    .unwrap()
                    .contains(&reference.1.name)
            {
                let qual = qual_sequence.clone();
                let lt = align_two_strings_passed_matrix(
                    &String::from_utf8(reference.1.name.clone()).unwrap(),
                    read_name,
                    &reference.1.sequence,
                    read,
                    qual,
                    my_aff_score,
                    alignment_mat,
                    &read.len(),
                );

                Some((lt, reference.1.sequence.clone(), reference.1.name.clone()))
            } else {
                None
            }
        })
        .filter(|x| x.is_some())
        .map(|c| c.unwrap());

    let ranked_alignments = ranked_alignments.into_iter().enumerate().max_by(|al, al2| {
        let score1 = al.1 .0.score; // / al.1.0.reference_aligned.len() as f64;
        let score2 = al2.1 .0.score; // / al2.1.0.reference_aligned.len() as f64;
        score1.partial_cmp(&score2).unwrap()
    });

    match ranked_alignments.iter().next() {
        None => None,
        Some((_x, y)) => {
            //Some((Some(y.0.clone()), , y.2.clone()))
            //println!("---- {} {}",String::from_utf8(y.1.clone()).unwrap(),String::from_utf8(y.2.clone()).unwrap());
            Some(AlignmentWithRef {
                alignment: Some(y.0.clone()),
                ref_name: y.2.clone(),
                ref_sequence: y.1.clone(),
            })
        }
    }
}

#[allow(dead_code)]
fn cigar_to_alignment(reference: &Vec<u8>, read: &Vec<u8>, cigar: &Vec<AlignmentOperation>) -> (Vec<u8>, Vec<u8>, Vec<AlignmentTag>) {
    let mut alignment_string1 = Vec::new();
    let mut alignment_string2 = Vec::new();
    let mut cigar_vec = Vec::new();

    let mut seq1_index = 0;
    let mut seq2_index = 0;

    for c in cigar {
        match c {
            Match => {
                alignment_string1.push(reference.get(seq1_index).unwrap().clone());
                alignment_string2.push(read.get(seq2_index).unwrap().clone());
                cigar_vec.push(AlignmentTag::MatchMismatch(1));
                seq1_index += 1;
                seq2_index += 1;
            }
            Del => {
                alignment_string1.push(reference.get(seq1_index).unwrap().clone());
                alignment_string2.push(FASTA_UNSET);
                cigar_vec.push(AlignmentTag::Del(1));
                seq1_index += 1;
            }
           Ins => {
                alignment_string1.push(FASTA_UNSET);
                alignment_string2.push(read.get(seq2_index).unwrap().clone());
                cigar_vec.push(AlignmentTag::Ins(1));
                seq2_index += 1;
            },
            Subst => {
                alignment_string1.push(reference.get(seq1_index).unwrap().clone());
                alignment_string2.push(read.get(seq2_index).unwrap().clone());
                cigar_vec.push(AlignmentTag::MatchMismatch(1));
                seq1_index += 1;
                seq2_index += 1;
            }
            _ => {
                panic!("Unknown cigar operation {:?}",c)
            }
        }
    }
    (alignment_string1, alignment_string2, simplify_cigar_string(&cigar_vec))
}

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

#[cfg(test)]
mod tests {
    use sigalign::algorithms::SemiGlobal;
    use sigalign::{Aligner, ReferenceBuilder};
    use std::collections::{BTreeMap, HashMap};
    use std::fs::File;
    use std::io::BufReader;

    use crate::alignment::alignment_matrix::{
        create_scoring_record_3d, AlignmentTag, AlignmentType,
    };
    use crate::alignment::scoring_functions::{AffineScoring, InversionScoring};
    use crate::alignment_functions::{exhaustive_alignment_search, simplify_cigar_string};
    use crate::read_strategies::sequence_layout::{
        AlignedReadOrientation, ReadPosition, SequenceLayout,
    };
    use crate::reference::fasta_reference::ReferenceManager;

    #[test]
    fn test_twist() {
        let ref_location = &"./test_data/18guide1.fa".to_string();
        let rm = ReferenceManager::from_fa_file(&ref_location, 8, 8);

        rm.references.iter().for_each(|reference| {
            let read_one = reference.1.sequence.clone().to_ascii_uppercase();
            println!("Reference: {}", u8s(&reference.1.name));
            let _read_structure = SequenceLayout {
                aligner: None,
                merge: None,
                reads: vec![ReadPosition::Read1 {
                    orientation: AlignedReadOrientation::Forward,
                }],
                known_strand: true,
                references: BTreeMap::new(),
            };

            let mut read_mat = create_scoring_record_3d(300, 300, AlignmentType::Affine, false);

            let _my_score = InversionScoring {
                match_score: 9.0,
                mismatch_score: -21.0,
                gap_open: -25.0,
                gap_extend: -1.0,
                inversion_penalty: -40.0,
                min_inversion_length: 20,
            };

            let my_aff_score = AffineScoring {
                match_score: 10.0,
                mismatch_score: -9.0,
                special_character_score: 9.0,
                gap_open: -20.0,
                gap_extend: -1.0,
                final_gap_multiplier: 1.0,
            };

            let best_ref = quick_alignment_search(
                &"testread".to_string(),
                &read_one,
                None,
                &&rm,
                &mut read_mat,
                &my_aff_score,
                &0.9,
            );
            match best_ref {
                None => {
                    println!(
                        "No reference found {}",
                        String::from_utf8(reference.1.name.clone()).unwrap()
                    );
                }
                Some(x) => {
                    assert_eq!(
                        String::from_utf8(x.ref_name).unwrap(),
                        String::from_utf8(reference.1.name.clone()).unwrap()
                    );
                }
            }
        });
    }
    #[test]
    fn test_twist_read() {
        let ref_location = &"./test_data/18guide1.fa".to_string();
        let rm = ReferenceManager::from_fa_file(&ref_location, 8, 8);
        let read_one = "GTGGAAAGGACGAAACACCGGTACTTTCGAAAGTACGCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTCGCATTCTACCGTGACTTTAGCAAGGTGATCATTCGCAACAGTATCGACGGGCGTACTTTCGAAAGTACGCCCGTCGATGTTTGAATTCGAATTTAAATCGGATCCGCGGCCAA".to_string().to_ascii_uppercase().into_bytes();

        let _read_structure = SequenceLayout {
            aligner: None,
            merge: None,
            reads: vec![ReadPosition::Read1 {
                orientation: AlignedReadOrientation::Forward,
            }],
            known_strand: true,
            references: BTreeMap::new(),
        };

        let mut read_mat = create_scoring_record_3d(
            read_one.len() + 100,
            read_one.len() + 100,
            AlignmentType::Affine,
            false,
        );

        let _my_score = InversionScoring {
            match_score: 9.0,
            mismatch_score: -21.0,
            gap_open: -25.0,
            gap_extend: -1.0,
            inversion_penalty: -40.0,
            min_inversion_length: 20,
        };

        let my_aff_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -9.0,
            special_character_score: 9.0,
            gap_open: -20.0,
            gap_extend: -1.0,
            final_gap_multiplier: 1.0,
        };
        println!("here");

        let best_ref = quick_alignment_search(
            &"testread".to_string(),
            &read_one,
            None,
            &&rm,
            &mut read_mat,
            &my_aff_score,
            &0.9,
        );
        assert_eq!(
            String::from_utf8(best_ref.unwrap().ref_name).unwrap(),
            String::from_utf8(
                "A5_GCGTACTTTCGAAAGTACGCCGG_18GUIDE_1_2_1"
                    .to_string()
                    .into_bytes()
            )
            .unwrap()
        );
    }

    #[test]
    fn test_twist_read2() {
        let ref_location = &"./test_data/18guide1.fa".to_string();
        let rm = ReferenceManager::from_fa_file(&ref_location, 8, 8);
        let read_one = "GTGGAAAGGACGAAACACCGACTCCCGCGCGGGAGTACGTTGTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTCGCATTCTACCGTGACTTTAGCAAGGTGATCATTCGCAACAGTATCGACCTCGTAAATTCGCGAATTTACGAGGTCGATGTTTGAATTCGAATTTAAATCGGATCCGCGGCCAA".to_string().to_ascii_uppercase().into_bytes();

        let _read_structure = SequenceLayout {
            aligner: None,
            merge: None,
            reads: vec![ReadPosition::Read1 {
                orientation: AlignedReadOrientation::Forward,
            }],
            known_strand: true,
            references: BTreeMap::new(),
        };

        let mut read_mat = create_scoring_record_3d(
            read_one.len() + 100,
            read_one.len() + 100,
            AlignmentType::Affine,
            false,
        );

        let _my_score = InversionScoring {
            match_score: 9.0,
            mismatch_score: -21.0,
            gap_open: -25.0,
            gap_extend: -1.0,
            inversion_penalty: -40.0,
            min_inversion_length: 20,
        };

        let my_aff_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -9.0,
            special_character_score: 9.0,
            gap_open: -20.0,
            gap_extend: -1.0,
            final_gap_multiplier: 1.0,
        };
        println!("here");

        let best_ref = quick_alignment_search(
            &"testread".to_string(),
            &read_one,
            None,
            &&rm,
            &mut read_mat,
            &my_aff_score,
            &0.9,
        );
        assert_eq!(
            String::from_utf8(best_ref.unwrap().ref_name).unwrap(),
            String::from_utf8(
                "A3_GTACTCCCGCGCGGGAGTACGGG_18GUIDE_1_2_1"
                    .to_string()
                    .into_bytes()
            )
            .unwrap()
        );
    }

    #[test]
    fn test_twist_read3() {
        let ref_location = &"./test_data/18guide1.fa".to_string();
        let rm = ReferenceManager::from_fa_file(&ref_location, 8, 8);
        let read_one = "GTGGAAAGGACGAAACACCGACAAGCGGCCGCTTGTAGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTCGCATTCTACCGTGACTTTAGCAAGGTGATCATTCGCAACAGTATCGACCGCTACAAGCGGCCGCTTGTAGCGGTCGATGTTTGAATTCGAATTTAAATCGGATCCGCGGCCAA".to_string().to_ascii_uppercase().into_bytes();

        let _read_structure = SequenceLayout {
            aligner: None,
            merge: None,
            reads: vec![ReadPosition::Read1 {
                orientation: AlignedReadOrientation::Forward,
            }],
            known_strand: true,
            references: BTreeMap::new(),
        };

        let mut read_mat = create_scoring_record_3d(
            read_one.len() + 100,
            read_one.len() + 100,
            AlignmentType::Affine,
            false,
        );

        let _my_score = InversionScoring {
            match_score: 9.0,
            mismatch_score: -21.0,
            gap_open: -25.0,
            gap_extend: -1.0,
            inversion_penalty: -40.0,
            min_inversion_length: 20,
        };

        let my_aff_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -9.0,
            special_character_score: 9.0,
            gap_open: -20.0,
            gap_extend: -1.0,
            final_gap_multiplier: 1.0,
        };
        println!("here");

        let best_ref = quick_alignment_search(
            &"testread".to_string(),
            &read_one,
            None,
            &&rm,
            &mut read_mat,
            &my_aff_score,
            &0.9,
        );
        assert_eq!(
            String::from_utf8(best_ref.unwrap().ref_name).unwrap(),
            String::from_utf8(
                "A356_CTACAAGCGGCCGCTTGTAGCGG_18GUIDE_1_2_1"
                    .to_string()
                    .into_bytes()
            )
            .unwrap()
        );
    }

    #[test]
    fn test_find_best_reference() {
        let ref_location = &"test_data/test_best_alignment.fasta".to_string();
        let rm = ReferenceManager::from_fa_file(&ref_location, 8, 8);

        let read_one = "atggactatcatatgcttaccgtaacttgaaagtatttcgatttcttggctttatatatcttgtggaaaggacgaaacaccgGGTAGCAAACGTTTGGACGTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTTTCCTGCAGGAAACCCCGGGgaat".to_string().to_ascii_uppercase().into_bytes();

        let _read_structure = SequenceLayout {
            aligner: None,
            merge: None,
            reads: vec![ReadPosition::Read1 {
                orientation: AlignedReadOrientation::Forward,
            }],
            known_strand: true,
            references: BTreeMap::new(),
        };

        let mut read_mat = create_scoring_record_3d(
            read_one.len() + 100,
            read_one.len() + 100,
            AlignmentType::Affine,
            false,
        );

        let _my_score = InversionScoring {
            match_score: 9.0,
            mismatch_score: -21.0,
            gap_open: -25.0,
            gap_extend: -1.0,
            inversion_penalty: -40.0,
            min_inversion_length: 20,
        };

        let my_aff_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -9.0,
            special_character_score: 9.0,
            gap_open: -20.0,
            gap_extend: -1.0,
            final_gap_multiplier: 1.0,
        };

        let best_ref = exhaustive_alignment_search(
            &"testread".to_string(),
            &read_one,
            None,
            &&rm,
            &mut read_mat,
            &my_aff_score,
            None,
        );
        assert_eq!(
            String::from_utf8(best_ref.unwrap().ref_name).unwrap(),
            String::from_utf8(
                "1_AAACCCCGGG_GGTAGCAAACGTTTGGACGTG"
                    .to_string()
                    .into_bytes()
            )
            .unwrap()
        );

        let read_one = "atggactatcatatgcttaccgtaacttgaaagtatttcgatttcttggctttatatatcttgtggaaaggacgaaacaccgGGTGCCCTTACTCTCACCTGATTACTTAATCCGTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTTTCCTGCAGGAACGCCCTACgaattcgggcccattggtatggc".to_string().to_ascii_uppercase().into_bytes();
        let best_ref = exhaustive_alignment_search(
            &"testread".to_string(),
            &read_one,
            None,
            &&rm,
            &mut read_mat,
            &my_aff_score,
            None,
        );

        assert_eq!(
            String::from_utf8(best_ref.unwrap().ref_name).unwrap(),
            String::from_utf8(
                "2_AACGCCCTAC_GGTGCCCTTACTCTCACCTGATTACTTAATCCGTG"
                    .to_string()
                    .into_bytes()
            )
            .unwrap()
        );
    }

    #[test]
    fn test_find_best_reference2() {
        let ref_location = &"test_data/test_ref_alignment.fasta".to_string();
        let rm = ReferenceManager::from_fa_file(&ref_location, 8, 8);

        let read_one = "ATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCGGTAAATTTGAGGCTCCGGCATGCAGGAGGCCGTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCTAGTCCGTTATCAACTTG".to_string().to_ascii_uppercase().into_bytes();

        let _read_structure = SequenceLayout {
            aligner: None,
            merge: None,
            reads: vec![ReadPosition::Read1 {
                orientation: AlignedReadOrientation::Forward,
            }],
            known_strand: true,
            references: BTreeMap::new(),
        };

        let mut read_mat = create_scoring_record_3d(
            read_one.len() + 100,
            read_one.len() + 100,
            AlignmentType::Affine,
            false,
        );

        let _my_score = InversionScoring {
            match_score: 9.0,
            mismatch_score: -21.0,
            gap_open: -25.0,
            gap_extend: -1.0,
            inversion_penalty: -40.0,
            min_inversion_length: 20,
        };

        let my_aff_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -9.0,
            special_character_score: 9.0,
            gap_open: -20.0,
            gap_extend: -1.0,
            final_gap_multiplier: 1.0,
        };

        let best_ref = exhaustive_alignment_search(
            &"testread".to_string(),
            &read_one,
            None,
            &&rm,
            &mut read_mat,
            &my_aff_score,
            None,
        );
        assert_eq!(
            String::from_utf8(best_ref.unwrap().ref_name).unwrap(),
            String::from_utf8(
                "ref_48_GGTAAATTTGAGGCTCCGGCATGCAGGAGGCCGTG"
                    .to_string()
                    .into_bytes()
            )
            .unwrap()
        );
    }

    #[test]
    fn simplify_cigar_test() {
        let input_cigar = vec![
            AlignmentTag::MatchMismatch(1),
            AlignmentTag::MatchMismatch(1),
            AlignmentTag::MatchMismatch(1),
        ];
        let merged_cigar = vec![AlignmentTag::MatchMismatch(3)];
        let resulting_cigar = simplify_cigar_string(&input_cigar);
        assert_eq!(resulting_cigar, merged_cigar);

        let input_cigar = vec![
            AlignmentTag::MatchMismatch(1),
            AlignmentTag::Ins(1),
            AlignmentTag::MatchMismatch(1),
            AlignmentTag::MatchMismatch(1),
        ];
        let merged_cigar = vec![
            AlignmentTag::MatchMismatch(1),
            AlignmentTag::Ins(1),
            AlignmentTag::MatchMismatch(2),
        ];
        let resulting_cigar = simplify_cigar_string(&input_cigar);
        assert_eq!(resulting_cigar, merged_cigar);
    }

    #[test]
    fn test_alignment_speed_sigaln() {
        let test_read = b"TTCCGATCTGTCATAACACCACACTAGAATCACGCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAGCGATGCAATTTCCTCATTTTATTAGGAAAGGACAGTGGGAGTGGCACCTTCCAGGGTCAAGGAAGGCACGGGGGAGGGGCAAACAACAGATGGCTGGCAACTAGAAGGCACAGTGAGCTTGTACATAACTACGCAAGTCCTTGCTAGGACCGGCCTTAAAGCCACGTGGCGGCCGCCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCGTGGTACAATATGCGTCTCCGAAATTAACCCGGTGCGTTTAAACGAAAAGGACCGACTACTACCTCGCGAAAGCTCTAAGCGTCGTGTCAGCGAAACTTCGCGGAGGTTCGACATCGAAAGACACGCGGGTGTATGTGGCGAAAGCAGCAACCTGATCTGGGGTGAAAAGCCATGGACGCCGGGACGAGAAAGGTCTAGGACTGTTTTGCGAGAAAAGGATTAGAGTTAGAATCGCGAAACGCTCGCGTTCCACCGCTCCGAAAGATCCCGAGGTCGTTTTACCGAAAGCGACGACTTCTGTCATAGTGAAACGATTGGACGTCTCTGGTGCGAAATCGCGGGTTGTACAACATACGAAACCGAGGCTACAACCCCGGACGAAAAGGTATAGGTAGCTAACACGCGAAACCCTAGGGATCGTGCTAGCCGAAAGCCCTATCACGCAGGGGACTGAAAAACATGGGCACGCCCCCGATGAAACGCTGCTTGTCTGGCCTCGCGAAAGAATGAGCAGAGCGTGAGGCGAAAAGCTTAAGCTGTGCACTCTCGAAAGTCGGTGTCCATCAGTGGATGAAACAGCGGGTTCCTGCTCCCGCGAAACGCCACCTGTACGTTACTTCGAAAATGAAGGGACAGCGGCGGACGAAAGTCATATTCCGTTGTGGTACGAAATTGGTCCTGATGCACGCACAGAAAAGATTGACCTCCGTTCGTACGAAAGCTCGGCCTCTGGGAGTCGTGAAAGACTCGGATCCGCACCAGATGAAAGGCACACCCACGCCCGTCACGAAAACCCAAACCTTGTATGTATGGAAATCTTCTGCGTCCGGGCCGCGGAAAAGCGTATACCTATCTCGCATGAAAGTCTCTCACCTCGTCTACGCGAAACGCTCGTACGCGTACGGGCTGAAAGCGATACACCGCTCGCCCCTGAAACCCTCTAGTTACGCGCCAGTGAAAGAGTCGCGTAGAGTACAGTGCAAGGTCGACAATCAACCTCTGGATTACATCCGATTGCCTCACTGTGCGAAAGTACTCGATGGCGTGGCTTAGAAAGCGTACAGTCTCCGTGCCGGGAAAATAAGAGCGCCTGCGGTTATGAAATCGTGGGCTACTCCTGGGTGGAAAGCTATCCTGCACATTAGTACGAAAGGTGCCAGGTTGCTTCGATCGAAAGCCCGAGAGATCACTCGTAGGAAACTACGCCGGTCACGACGGGCGAAACGACATGAACTCATCCGGACGAAAGGTAGTCCTTACGGTGATCTGCTAGGGTCTCTCCTAGCAACGGTTACTCCATCTGGTACACCCCCTGCTCGGGGCAAGTACCTGATGCGGCACAATGTCTAGCAGGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTTAACCTTAGCAATACGTCTGTCGAAGCAGCTACAA";

        /* (1) Build `Reference`
        let fasta =
            br#">target_1
CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTAGGAAAGGACAGTGGGAGTGGCACCTTCCAGGGTCAAGGAAGGCACGGGGGAGGGGCAAACAACAGATGGCTGGCAACTAGAAGGCACAGTGAGCTTGTACATAACTACGCAAGTCCTTGCTAGGACCGGCCTTAAAGCCACGTGGCGGCCGCCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCGTGGTACAATATGCGTCTCCGAAATTAACCCGGTGTGTTTAAACGAAAAGGACCGACTACTACCTCGCGAAAGCTCTAAGTGTTGTGTCAGCGAAACTTCGCGGAGGTTCGACATCGAAAGACACGCGGGTGTATATGGCGAAAGCAGCAACCTGATCTGGGGTGAAAAGCCATGGATGTCGGGACGAGAAAGGTCTAGGACTGTTTTGCGAGAAAAGGATTAGAGTTAGAATCGCGAAACGCTCGCGTTCTACCGCTCCGAAAGATCCCGAGGTTGTTTTACCGAAAGCGACGACTTCTGTCATAGTGAAACGATTGGACGTCTCTGGTGCGAAATCGCGGGTTGTACAACATACGAAACCGAGGCTATAATCCCGGACGAAAGGTATAGGTAGCTAACACGCGAAACCCTAGGGATCGTGCTAGCCGAAAGCCCTATTATGTAGGGGACTGAAAAACATGGGTACGTCCCCGATGAAACGCTGCTTGTCTGGCCTCGCGAAAGAATGAGCTGAGTGTGAGGCGAAAAGCTTAAGCTGTGCACTCTCGAAAGTCGGTGTCTATTAGTGGATGAAACAGCGGGTTCCTGCTCCCGCGAAACGCCACCTGTATGTTACTTCGAAAATGAAGGGATAGTGGCGGACGAAAGTCATATTCCGTTGTGGTACGAAATTGGTCCTGATGTACGCACAGAAAAGATTGACCTCTGTTCGTACGAAAGCTCGGCCTCTGGGAGTCGTGAAAGACTCGGATCCGTACCAGATGAAAGGCACACCCATGTCCGTCACGAAAACCCAAACCTTGTATGTATGGAAATCTTCTGCGTTCGGGCCGCGGAAAAGCGTATACCTATCTCGCATGAAAGTCTCTTATCTTGTCTACGCGAAACGCTCGTATGCGTACGGGCTGAAAGCGATATACTGTTCGCCCCTGAAACCCTCTAGTTATGCGCCAGTGAAAGAGTCGCGTAGAGTACAGTGCAAGGTCGACAATCAACCTCTGGATTACATCCGATTGCCTTACTGTGCGAAAGTACTCGATGGTGTGGCTTAGAAAGCGTACAGTCTCTGTGCCGGGAAAATAAGAGCGTCTGCGGTTATGAAATCGTGGGCTACTCCTGGGTGGAAAGCTATCCTGTATATTAGTACGAAAGGTGCCAGGTTGCTTCGATCGAAAGCCCGAGAGATTACTCGTAGGAAACTACGCCGGTTACGACGGGCGAAACGACATGAACTTATCCGGACGAAAGGTAGTCCTTACGGTGATCTGCTAGGGTCTCTCCTAGCAACGGTTACTCGATTTGGTACNNNNNNNNNNNNNNNNNNGTACCTGATGCGGCACAATGTCTAGC
>target_2
CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTAGGAAAGGACAGTGGGAGTGGCACCTTCCAGGGTCAAGGAAGGCACGGGGGAGGGGCAAACAACAGATGGCTGGCAACTAGAAGGCACAGTGAGCTTGTACATAACTACGCAAGTCCTTGCTAGGACCGGCCTTAAAGCCACGTGGCGGCCGCCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCGTGGTACAATATGCGTCTCCGAAATTAACCCGGTGTGTTTAAACGAAAAGGACCGACTACTACCTCGCGAAAGCTCTAAGTGTTGTGTCAGCGAAACTTCGCGGAGGTTCGACATCGAAAGACACGCGGGTGTATATGGCGAAAGCAGCAACCTGATCTGGGGTGAAAAGCCATGGATGTCGGGACGAGAAAGGTCTAGGACTGTTTTGCGAGAAAAGGATTAGAGTTAGAATCGCGAAACGCTCGCGTTCTACCGCTCCGAAAGATCCCGAGGTTGTTTTACCGAAAGCGACGACTTCTGTCATAGTGAAACGATTGGACGTCTCTGGTGCGAAATCGCGGGTTGTACAACATACGAAACCGAGGCTATAATCCCGGACGAAAGGTATAGGTAGCTAACACGCGAAACCCTAGGGATCGTGCTAGCCGAAAGCCCTATTATGTAGGGGACTGAAAAACATGGGTACGTCCCCGATGAAACGCTGCTTGTCTGGCCTCGCGAAAGAATGAGCTGAGTGTGAGGCGAAAAGCTTAAGCTGTGCACTCTCGAAAGTCGGTGTCTATTAGTGGATGAAACAGCGGGTTCCTGCTCCCGCGAAACGCCACCTGTATGTTACTTCGAAAATGAAGGGATAGTGGCGGACGAAAGTCATATTCCGTTGTGGTACGAAATTGGTCCTGATGTACGCACAGAAAAGATTGACCTCTGTTCGTACGAAAGCTCGGCCTCTGGGAGTCGTGAAAGACTCGGATCCGTACCAGATGAAAGGCACACCCATGTCCGTCACGAAAACCCAAACCTTGTATGTATGGAAATCTTCTGCGTTCGGGCCGCGGAAAAGCGTATACCTATCTCGCATGAAAGTCTCTTATCTTGTCTACGCGAAACGCTCGTATGCGTACGGGCTGAAAGCGATATACTGTTCGCCCCTGAAACCCTCTAGTTATGCGCCAGTGAAAGAGTCGCGTAGAGTACAGTGCAAGGTCGACAATCAACCTCTGGATTACATCCGATTGCCTTACTGTGCGAAAGTACTCGATGGTGTGGCTTAGAAAGCGTACAGTCTCTGTGCCGGGAAAATAAGAGCGTCTGCGGTTATGAAATCGTGGGCTACTCCTGGGTGGAAAGCTATCCTGTATATTAGTACGAAAGGTGCCAGGTTGCTTCGATCGAAAGCCCGAGAGATTACTCGTAGGAAACTACGCCGGTTACGACGGGCACGTACGTACGTTTTTGGGGGTGTGTGTGTTTGTCCTTACGGTGATCTGCTAGGGTCTCTCCTAGCAACGGTTACTCGATTTGGTACNNNNNNNNNNNNNNNNNNGTACCTGATGCGGCACAATGTCTAGC"#;
*/
        let fasta =
            br#">target_1
CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTAGGAAAGGACAGTGGGAGTGGCACCTTCCAGGGTCAAGGAAGGCACGGGGGAGGGGCAAACAACAGATGGCTGGCAACTAGAAGGCACAGTGAGCTTGTACATAACTACGCAAGTCCTTGCTAGGACCGGCCTTAAAGCCACGTGGCGGCCGCCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCGTGGTACAATATGCGTCTCCGAAATTAACCCGGTGTGTTTAAACGAAAAGGACCGACTACTACCTCGCGAAAGCTCTAAGTGTTGTGTCAGCGAAACTTCGCGGAGGTTCGACATCGAAAGACACGCGGGTGTATATGGCGAAAGCAGCAACCTGATCTGGGGTGAAAAGCCATGGATGTCGGGACGAGAAAGGTCTAGGACTGTTTTGCGAGAAAAGGATTAGAGTTAGAATCGCGAAACGCTCGCGTTCTACCGCTCCGAAAGATCCCGAGGTTGTTTTACCGAAAGCGACGACTTCTGTCATAGTGAAACGATTGGACGTCTCTGGTGCGAAATCGCGGGTTGTACAACATACGAAACCGAGGCTATAATCCCGGACGAAAGGTATAGGTAGCTAACACGCGAAACCCTAGGGATCGTGCTAGCCGAAAGCCCTATTATGTAGGGGACTGAAAAACATGGGTACGTCCCCGATGAAACGCTGCTTGTCTGGCCTCGCGAAAGAATGAGCTGAGTGTGAGGCGAAAAGCTTAAGCTGTGCACTCTCGAAAGTCGGTGTCTATTAGTGGATGAAACAGCGGGTTCCTGCTCCCGCGAAACGCCACCTGTATGTTACTTCGAAAATGAAGGGATAGTGGCGGACGAAAGTCATATTCCGTTGTGGTACGAAATTGGTCCTGATGTACGCACAGAAAAGATTGACCTCTGTTCGTACGAAAGCTCGGCCTCTGGGAGTCGTGAAAGACTCGGATCCGTACCAGATGAAAGGCACACCCATGTCCGTCACGAAAACCCAAACCTTGTATGTATGGAAATCTTCTGCGTTCGGGCCGCGGAAAAGCGTATACCTATCTCGCATGAAAGTCTCTTATCTTGTCTACGCGAAACGCTCGTATGCGTACGGGCTGAAAGCGATATACTGTTCGCCCCTGAAACCCTCTAGTTATGCGCCAGTGAAAGAGTCGCGTAGAGTACAGTGCAAGGTCGACAATCAACCTCTGGATTACATCCGATTGCCTTACTGTGCGAAAGTACTCGATGGTGTGGCTTAGAAAGCGTACAGTCTCTGTGCCGGGAAAATAAGAGCGTCTGCGGTTATGAAATCGTGGGCTACTCCTGGGTGGAAAGCTATCCTGTATATTAGTACGAAAGGTGCCAGGTTGCTTCGATCGAAAGCCCGAGAGATTACTCGTAGGAAACTACGCCGGTTACGACGGGCGAAACGACATGAACTTATCCGGACGAAAGGTAGTCCTTACGGTGATCTGCTAGGGTCTCTCCTAGCAACGGTTACTCGATTTGGTACNNNNNNNNNNNNNNNNNNGTACCTGATGCGGCACAATGTCTAGC"#;

        let reference = ReferenceBuilder::new()
            .set_uppercase(true) // Ignore case
            .ignore_base(b'N') // 'N' is never matched
            .add_fasta(&fasta[..])
            .unwrap() // Add sequences from FASTA
            .add_target("target_3", b"AAAAAAAAAAA") // Add sequence manually
            .build()
            .unwrap();

        // (2) Initialize `Aligner`
        let algorithm = SemiGlobal::new(
            4,   // Mismatch penalty
            6,   // Gap-open penalty
            2,   // Gap-extend penalty
            50,  // Minimum length
            0.2, // Maximum penalty per length
        )
        .unwrap();
        let mut aligner = Aligner::new(algorithm);

        // (3) Align query to reference
        for _i in 0..10000 {
            let _ = aligner.align(test_read, &reference);
        }
        //println!("{:#?}", result);
    }

    #[test]
    fn test_alignment_speed_seqam() {
        let test_read = b"TTCCGATCTGTCATAACACCACACTAGAATCACGCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAGCGATGCAATTTCCTCATTTTATTAGGAAAGGACAGTGGGAGTGGCACCTTCCAGGGTCAAGGAAGGCACGGGGGAGGGGCAAACAACAGATGGCTGGCAACTAGAAGGCACAGTGAGCTTGTACATAACTACGCAAGTCCTTGCTAGGACCGGCCTTAAAGCCACGTGGCGGCCGCCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCGTGGTACAATATGCGTCTCCGAAATTAACCCGGTGCGTTTAAACGAAAAGGACCGACTACTACCTCGCGAAAGCTCTAAGCGTCGTGTCAGCGAAACTTCGCGGAGGTTCGACATCGAAAGACACGCGGGTGTATGTGGCGAAAGCAGCAACCTGATCTGGGGTGAAAAGCCATGGACGCCGGGACGAGAAAGGTCTAGGACTGTTTTGCGAGAAAAGGATTAGAGTTAGAATCGCGAAACGCTCGCGTTCCACCGCTCCGAAAGATCCCGAGGTCGTTTTACCGAAAGCGACGACTTCTGTCATAGTGAAACGATTGGACGTCTCTGGTGCGAAATCGCGGGTTGTACAACATACGAAACCGAGGCTACAACCCCGGACGAAAAGGTATAGGTAGCTAACACGCGAAACCCTAGGGATCGTGCTAGCCGAAAGCCCTATCACGCAGGGGACTGAAAAACATGGGCACGCCCCCGATGAAACGCTGCTTGTCTGGCCTCGCGAAAGAATGAGCAGAGCGTGAGGCGAAAAGCTTAAGCTGTGCACTCTCGAAAGTCGGTGTCCATCAGTGGATGAAACAGCGGGTTCCTGCTCCCGCGAAACGCCACCTGTACGTTACTTCGAAAATGAAGGGACAGCGGCGGACGAAAGTCATATTCCGTTGTGGTACGAAATTGGTCCTGATGCACGCACAGAAAAGATTGACCTCCGTTCGTACGAAAGCTCGGCCTCTGGGAGTCGTGAAAGACTCGGATCCGCACCAGATGAAAGGCACACCCACGCCCGTCACGAAAACCCAAACCTTGTATGTATGGAAATCTTCTGCGTCCGGGCCGCGGAAAAGCGTATACCTATCTCGCATGAAAGTCTCTCACCTCGTCTACGCGAAACGCTCGTACGCGTACGGGCTGAAAGCGATACACCGCTCGCCCCTGAAACCCTCTAGTTACGCGCCAGTGAAAGAGTCGCGTAGAGTACAGTGCAAGGTCGACAATCAACCTCTGGATTACATCCGATTGCCTCACTGTGCGAAAGTACTCGATGGCGTGGCTTAGAAAGCGTACAGTCTCCGTGCCGGGAAAATAAGAGCGCCTGCGGTTATGAAATCGTGGGCTACTCCTGGGTGGAAAGCTATCCTGCACATTAGTACGAAAGGTGCCAGGTTGCTTCGATCGAAAGCCCGAGAGATCACTCGTAGGAAACTACGCCGGTCACGACGGGCGAAACGACATGAACTCATCCGGACGAAAGGTAGTCCTTACGGTGATCTGCTAGGGTCTCTCCTAGCAACGGTTACTCCATCTGGTACACCCCCTGCTCGGGGCAAGTACCTGATGCGGCACAATGTCTAGCAGGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTTAACCTTAGCAATACGTCTGTCGAAGCAGCTACAA";

        // (1) Build `Reference`
        let fasta =
            br#">target_1
CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTAGGAAAGGACAGTGGGAGTGGCACCTTCCAGGGTCAAGGAAGGCACGGGGGAGGGGCAAACAACAGATGGCTGGCAACTAGAAGGCACAGTGAGCTTGTACATAACTACGCAAGTCCTTGCTAGGACCGGCCTTAAAGCCACGTGGCGGCCGCCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCGTGGTACAATATGCGTCTCCGAAATTAACCCGGTGTGTTTAAACGAAAAGGACCGACTACTACCTCGCGAAAGCTCTAAGTGTTGTGTCAGCGAAACTTCGCGGAGGTTCGACATCGAAAGACACGCGGGTGTATATGGCGAAAGCAGCAACCTGATCTGGGGTGAAAAGCCATGGATGTCGGGACGAGAAAGGTCTAGGACTGTTTTGCGAGAAAAGGATTAGAGTTAGAATCGCGAAACGCTCGCGTTCTACCGCTCCGAAAGATCCCGAGGTTGTTTTACCGAAAGCGACGACTTCTGTCATAGTGAAACGATTGGACGTCTCTGGTGCGAAATCGCGGGTTGTACAACATACGAAACCGAGGCTATAATCCCGGACGAAAGGTATAGGTAGCTAACACGCGAAACCCTAGGGATCGTGCTAGCCGAAAGCCCTATTATGTAGGGGACTGAAAAACATGGGTACGTCCCCGATGAAACGCTGCTTGTCTGGCCTCGCGAAAGAATGAGCTGAGTGTGAGGCGAAAAGCTTAAGCTGTGCACTCTCGAAAGTCGGTGTCTATTAGTGGATGAAACAGCGGGTTCCTGCTCCCGCGAAACGCCACCTGTATGTTACTTCGAAAATGAAGGGATAGTGGCGGACGAAAGTCATATTCCGTTGTGGTACGAAATTGGTCCTGATGTACGCACAGAAAAGATTGACCTCTGTTCGTACGAAAGCTCGGCCTCTGGGAGTCGTGAAAGACTCGGATCCGTACCAGATGAAAGGCACACCCATGTCCGTCACGAAAACCCAAACCTTGTATGTATGGAAATCTTCTGCGTTCGGGCCGCGGAAAAGCGTATACCTATCTCGCATGAAAGTCTCTTATCTTGTCTACGCGAAACGCTCGTATGCGTACGGGCTGAAAGCGATATACTGTTCGCCCCTGAAACCCTCTAGTTATGCGCCAGTGAAAGAGTCGCGTAGAGTACAGTGCAAGGTCGACAATCAACCTCTGGATTACATCCGATTGCCTTACTGTGCGAAAGTACTCGATGGTGTGGCTTAGAAAGCGTACAGTCTCTGTGCCGGGAAAATAAGAGCGTCTGCGGTTATGAAATCGTGGGCTACTCCTGGGTGGAAAGCTATCCTGTATATTAGTACGAAAGGTGCCAGGTTGCTTCGATCGAAAGCCCGAGAGATTACTCGTAGGAAACTACGCCGGTTACGACGGGCGAAACGACATGAACTTATCCGGACGAAAGGTAGTCCTTACGGTGATCTGCTAGGGTCTCTCCTAGCAACGGTTACTCGATTTGGTACNNNNNNNNNNNNNNNNNNGTACCTGATGCGGCACAATGTCTAGC
>target_2
CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTAGGAAAGGACAGTGGGAGTGGCACCTTCCAGGGTCAAGGAAGGCACGGGGGAGGGGCAAACAACAGATGGCTGGCAACTAGAAGGCACAGTGAGCTTGTACATAACTACGCAAGTCCTTGCTAGGACCGGCCTTAAAGCCACGTGGCGGCCGCCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCGTGGTACAATATGCGTCTCCGAAATTAACCCGGTGTGTTTAAACGAAAAGGACCGACTACTACCTCGCGAAAGCTCTAAGTGTTGTGTCAGCGAAACTTCGCGGAGGTTCGACATCGAAAGACACGCGGGTGTATATGGCGAAAGCAGCAACCTGATCTGGGGTGAAAAGCCATGGATGTCGGGACGAGAAAGGTCTAGGACTGTTTTGCGAGAAAAGGATTAGAGTTAGAATCGCGAAACGCTCGCGTTCTACCGCTCCGAAAGATCCCGAGGTTGTTTTACCGAAAGCGACGACTTCTGTCATAGTGAAACGATTGGACGTCTCTGGTGCGAAATCGCGGGTTGTACAACATACGAAACCGAGGCTATAATCCCGGACGAAAGGTATAGGTAGCTAACACGCGAAACCCTAGGGATCGTGCTAGCCGAAAGCCCTATTATGTAGGGGACTGAAAAACATGGGTACGTCCCCGATGAAACGCTGCTTGTCTGGCCTCGCGAAAGAATGAGCTGAGTGTGAGGCGAAAAGCTTAAGCTGTGCACTCTCGAAAGTCGGTGTCTATTAGTGGATGAAACAGCGGGTTCCTGCTCCCGCGAAACGCCACCTGTATGTTACTTCGAAAATGAAGGGATAGTGGCGGACGAAAGTCATATTCCGTTGTGGTACGAAATTGGTCCTGATGTACGCACAGAAAAGATTGACCTCTGTTCGTACGAAAGCTCGGCCTCTGGGAGTCGTGAAAGACTCGGATCCGTACCAGATGAAAGGCACACCCATGTCCGTCACGAAAACCCAAACCTTGTATGTATGGAAATCTTCTGCGTTCGGGCCGCGGAAAAGCGTATACCTATCTCGCATGAAAGTCTCTTATCTTGTCTACGCGAAACGCTCGTATGCGTACGGGCTGAAAGCGATATACTGTTCGCCCCTGAAACCCTCTAGTTATGCGCCAGTGAAAGAGTCGCGTAGAGTACAGTGCAAGGTCGACAATCAACCTCTGGATTACATCCGATTGCCTTACTGTGCGAAAGTACTCGATGGTGTGGCTTAGAAAGCGTACAGTCTCTGTGCCGGGAAAATAAGAGCGTCTGCGGTTATGAAATCGTGGGCTACTCCTGGGTGGAAAGCTATCCTGTATATTAGTACGAAAGGTGCCAGGTTGCTTCGATCGAAAGCCCGAGAGATTACTCGTAGGAAACTACGCCGGTTACGACGGGCACGTACGTACGTTTTTGGGGGTGTGTGTGTTTGTCCTTACGGTGATCTGCTAGGGTCTCTCCTAGCAACGGTTACTCGATTTGGTACNNNNNNNNNNNNNNNNNNGTACCTGATGCGGCACAATGTCTAGC"#;
        let reference = ReferenceBuilder::new()
            .set_uppercase(true) // Ignore case
            .ignore_base(b'N') // 'N' is never matched
            .add_fasta(&fasta[..])
            .unwrap() // Add sequences from FASTA
            .add_target("target_3", b"AAAAAAAAAAA") // Add sequence manually
            .build()
            .unwrap();

        // (2) Initialize `Aligner`
        let algorithm = SemiGlobal::new(
            4,   // Mismatch penalty
            6,   // Gap-open penalty
            2,   // Gap-extend penalty
            50,  // Minimum length
            0.2, // Maximum penalty per length
        )
        .unwrap();
        let mut aligner = Aligner::new(algorithm);

        // (3) Align query to reference
        for _i in 0..10 {
            let _ = aligner.align(test_read, &reference);
        }
        //println!("{:#?}", result);
    }

    #[test]
    fn test_alignment_speed_wfa() {
        let reference = "CTACACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATTAGGAAAGGACAGTGGGAGTGGCACCTTCCAGGGTCAAGGAAGGCACGGGGGAGGGGCAAACAACAGATGGCTGGCAACTAGAAGGCACAGTGAGCTTGTACATAACTACGCAAGTCCTTGCTAGGACCGGCCTTAAAGCCACGTGGCGGCCGCCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCGTGGTACAATATGCGTCTCCGAAATTAACCCGGTGTGTTTAAACGAAAAGGACCGACTACTACCTCGCGAAAGCTCTAAGTGTTGTGTCAGCGAAACTTCGCGGAGGTTCGACATCGAAAGACACGCGGGTGTATATGGCGAAAGCAGCAACCTGATCTGGGGTGAAAAGCCATGGATGTCGGGACGAGAAAGGTCTAGGACTGTTTTGCGAGAAAAGGATTAGAGTTAGAATCGCGAAACGCTCGCGTTCTACCGCTCCGAAAGATCCCGAGGTTGTTTTACCGAAAGCGACGACTTCTGTCATAGTGAAACGATTGGACGTCTCTGGTGCGAAATCGCGGGTTGTACAACATACGAAACCGAGGCTATAATCCCGGACGAAAGGTATAGGTAGCTAACACGCGAAACCCTAGGGATCGTGCTAGCCGAAAGCCCTATTATGTAGGGGACTGAAAAACATGGGTACGTCCCCGATGAAACGCTGCTTGTCTGGCCTCGCGAAAGAATGAGCTGAGTGTGAGGCGAAAAGCTTAAGCTGTGCACTCTCGAAAGTCGGTGTCTATTAGTGGATGAAACAGCGGGTTCCTGCTCCCGCGAAACGCCACCTGTATGTTACTTCGAAAATGAAGGGATAGTGGCGGACGAAAGTCATATTCCGTTGTGGTACGAAATTGGTCCTGATGTACGCACAGAAAAGATTGACCTCTGTTCGTACGAAAGCTCGGCCTCTGGGAGTCGTGAAAGACTCGGATCCGTACCAGATGAAAGGCACACCCATGTCCGTCACGAAAACCCAAACCTTGTATGTATGGAAATCTTCTGCGTTCGGGCCGCGGAAAAGCGTATACCTATCTCGCATGAAAGTCTCTTATCTTGTCTACGCGAAACGCTCGTATGCGTACGGGCTGAAAGCGATATACTGTTCGCCCCTGAAACCCTCTAGTTATGCGCCAGTGAAAGAGTCGCGTAGAGTACAGTGCAAGGTCGACAATCAACCTCTGGATTACATCCGATTGCCTTACTGTGCGAAAGTACTCGATGGTGTGGCTTAGAAAGCGTACAGTCTCTGTGCCGGGAAAATAAGAGCGTCTGCGGTTATGAAATCGTGGGCTACTCCTGGGTGGAAAGCTATCCTGTATATTAGTACGAAAGGTGCCAGGTTGCTTCGATCGAAAGCCCGAGAGATTACTCGTAGGAAACTACGCCGGTTACGACGGGCGAAACGACATGAACTTATCCGGACGAAAGGTAGTCCTTACGGTGATCTGCTAGGGTCTCTCCTAGCAACGGTTACTCGATTTGGTACNNNNNNNNNNNNNNNNNNGTACCTGATGCGGCACAATGTCTAGC";
        let test_read = "TTCCGATCTGTCATAACACCACACTAGAATCACGCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAGCGATGCAATTTCCTCATTTTATTAGGAAAGGACAGTGGGAGTGGCACCTTCCAGGGTCAAGGAAGGCACGGGGGAGGGGCAAACAACAGATGGCTGGCAACTAGAAGGCACAGTGAGCTTGTACATAACTACGCAAGTCCTTGCTAGGACCGGCCTTAAAGCCACGTGGCGGCCGCCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCGTGGTACAATATGCGTCTCCGAAATTAACCCGGTGCGTTTAAACGAAAAGGACCGACTACTACCTCGCGAAAGCTCTAAGCGTCGTGTCAGCGAAACTTCGCGGAGGTTCGACATCGAAAGACACGCGGGTGTATGTGGCGAAAGCAGCAACCTGATCTGGGGTGAAAAGCCATGGACGCCGGGACGAGAAAGGTCTAGGACTGTTTTGCGAGAAAAGGATTAGAGTTAGAATCGCGAAACGCTCGCGTTCCACCGCTCCGAAAGATCCCGAGGTCGTTTTACCGAAAGCGACGACTTCTGTCATAGTGAAACGATTGGACGTCTCTGGTGCGAAATCGCGGGTTGTACAACATACGAAACCGAGGCTACAACCCCGGACGAAAAGGTATAGGTAGCTAACACGCGAAACCCTAGGGATCGTGCTAGCCGAAAGCCCTATCACGCAGGGGACTGAAAAACATGGGCACGCCCCCGATGAAACGCTGCTTGTCTGGCCTCGCGAAAGAATGAGCAGAGCGTGAGGCGAAAAGCTTAAGCTGTGCACTCTCGAAAGTCGGTGTCCATCAGTGGATGAAACAGCGGGTTCCTGCTCCCGCGAAACGCCACCTGTACGTTACTTCGAAAATGAAGGGACAGCGGCGGACGAAAGTCATATTCCGTTGTGGTACGAAATTGGTCCTGATGCACGCACAGAAAAGATTGACCTCCGTTCGTACGAAAGCTCGGCCTCTGGGAGTCGTGAAAGACTCGGATCCGCACCAGATGAAAGGCACACCCACGCCCGTCACGAAAACCCAAACCTTGTATGTATGGAAATCTTCTGCGTCCGGGCCGCGGAAAAGCGTATACCTATCTCGCATGAAAGTCTCTCACCTCGTCTACGCGAAACGCTCGTACGCGTACGGGCTGAAAGCGATACACCGCTCGCCCCTGAAACCCTCTAGTTACGCGCCAGTGAAAGAGTCGCGTAGAGTACAGTGCAAGGTCGACAATCAACCTCTGGATTACATCCGATTGCCTCACTGTGCGAAAGTACTCGATGGCGTGGCTTAGAAAGCGTACAGTCTCCGTGCCGGGAAAATAAGAGCGCCTGCGGTTATGAAATCGTGGGCTACTCCTGGGTGGAAAGCTATCCTGCACATTAGTACGAAAGGTGCCAGGTTGCTTCGATCGAAAGCCCGAGAGATCACTCGTAGGAAACTACGCCGGTCACGACGGGCGAAACGACATGAACTCATCCGGACGAAAGGTAGTCCTTACGGTGATCTGCTAGGGTCTCTCCTAGCAACGGTTACTCCATCTGGTACACCCCCTGCTCGGGGCAAGTACCTGATGCGGCACAATGTCTAGCAGGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTTAACCTTAGCAATACGTCTGTCGAAGCAGCTACAA";
        let mut cnt = 0;
        (0..10000).for_each(|_x| {
            cnt += 1;
            let alignment = wfa_alignment(reference.as_bytes(), test_read.as_bytes(),&4, &10, &1);
            //println!("{:?}", alignment);
            let _alignment = cigar_to_alignment(
                &reference.as_bytes().to_vec(),
                &test_read.as_bytes().to_vec(),
                &alignment.into_bytes(),
            );
        });
        println!("{:?}", cnt);
        //println!("{:#?}", result);
    }

    #[test]
    fn test_alignment_wfa() {
        let reference = "CCGTAGATTAACTGCTTGCCCGTAGATTAACTGCTTGCCCGTAGATTAACTGCTTGCCCGTAGATTAACTGCTTGCCCGTAGATTAACTGCTTGCCCGTAGATTAACTGCTTGCCCGTAGATTAACTGCTTGCCCGTAGATTAACTGCTTGCCCGTAGATTAACTGCTTGC";
        let test_read = "CCGTAGATTTATGCTTGCCCGTAGATTTATGCTTGCCCGTAGATTTATGCTTGCCCGTAGATTTATGCTTGCCCGTAGATTTATGCTTGCCCGTAGATTTATGCTTGCCCGTAGATTTATGCTTGCCCGTAGATTTATGCTTGCCCGTAGATTTATGCTTGC";
        let alignment = wfa_alignment(reference.as_bytes(), test_read.as_bytes(), &4, &10, &1);
        println!("{}", alignment);

        let alignment = cigar_to_alignment(
            &reference.as_bytes().to_vec(),
            &test_read.as_bytes().to_vec(),
            &alignment.into_bytes(),
        );

        println!("{}\n{}", String::from_utf8(alignment.0).unwrap(), String::from_utf8(alignment.1).unwrap());
    }

    #[test]
    fn test_how_fast_do_we_do_100_long_reads() {
        let ref_location = &"test_data/larger_two_reference_example.fa".to_string();
        let rm = ReferenceManager::from_fa_file(&ref_location, 8, 8);

        let read_one = "TTCCGATCTGTCATAACACCACACTAGAATCACGCGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAGCGATGCAATTTCCTCATTTTATTAGGAAAGGACAGTGGGAGTGGCACCTTCCAGGGTCAAGGAAGGCACGGGGGAGGGGCAAACAACAGATGGCTGGCAACTAGAAGGCACAGTGAGCTTGTACATAACTACGCAAGTCCTTGCTAGGACCGGCCTTAAAGCCACGTGGCGGCCGCCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCGTGGTACAATATGCGTCTCCGAAATTAACCCGGTGCGTTTAAACGAAAAGGACCGACTACTACCTCGCGAAAGCTCTAAGCGTCGTGTCAGCGAAACTTCGCGGAGGTTCGACATCGAAAGACACGCGGGTGTATGTGGCGAAAGCAGCAACCTGATCTGGGGTGAAAAGCCATGGACGCCGGGACGAGAAAGGTCTAGGACTGTTTTGCGAGAAAAGGATTAGAGTTAGAATCGCGAAACGCTCGCGTTCCACCGCTCCGAAAGATCCCGAGGTCGTTTTACCGAAAGCGACGACTTCTGTCATAGTGAAACGATTGGACGTCTCTGGTGCGAAATCGCGGGTTGTACAACATACGAAACCGAGGCTACAACCCCGGACGAAAAGGTATAGGTAGCTAACACGCGAAACCCTAGGGATCGTGCTAGCCGAAAGCCCTATCACGCAGGGGACTGAAAAACATGGGCACGCCCCCGATGAAACGCTGCTTGTCTGGCCTCGCGAAAGAATGAGCAGAGCGTGAGGCGAAAAGCTTAAGCTGTGCACTCTCGAAAGTCGGTGTCCATCAGTGGATGAAACAGCGGGTTCCTGCTCCCGCGAAACGCCACCTGTACGTTACTTCGAAAATGAAGGGACAGCGGCGGACGAAAGTCATATTCCGTTGTGGTACGAAATTGGTCCTGATGCACGCACAGAAAAGATTGACCTCCGTTCGTACGAAAGCTCGGCCTCTGGGAGTCGTGAAAGACTCGGATCCGCACCAGATGAAAGGCACACCCACGCCCGTCACGAAAACCCAAACCTTGTATGTATGGAAATCTTCTGCGTCCGGGCCGCGGAAAAGCGTATACCTATCTCGCATGAAAGTCTCTCACCTCGTCTACGCGAAACGCTCGTACGCGTACGGGCTGAAAGCGATACACCGCTCGCCCCTGAAACCCTCTAGTTACGCGCCAGTGAAAGAGTCGCGTAGAGTACAGTGCAAGGTCGACAATCAACCTCTGGATTACATCCGATTGCCTCACTGTGCGAAAGTACTCGATGGCGTGGCTTAGAAAGCGTACAGTCTCCGTGCCGGGAAAATAAGAGCGCCTGCGGTTATGAAATCGTGGGCTACTCCTGGGTGGAAAGCTATCCTGCACATTAGTACGAAAGGTGCCAGGTTGCTTCGATCGAAAGCCCGAGAGATCACTCGTAGGAAACTACGCCGGTCACGACGGGCGAAACGACATGAACTCATCCGGACGAAAGGTAGTCCTTACGGTGATCTGCTAGGGTCTCTCCTAGCAACGGTTACTCCATCTGGTACACCCCCTGCTCGGGGCAAGTACCTGATGCGGCACAATGTCTAGCAGGTGCTGAAGAAAGTTGTCGGTGTCTTTGTGTTAACCTTAGCAATACGTCTGTCGAAGCAGCTACAA".to_string().to_ascii_uppercase().into_bytes();

        let _read_structure = SequenceLayout {
            aligner: None,
            merge: None,
            reads: vec![ReadPosition::Read1 {
                orientation: AlignedReadOrientation::Forward,
            }],
            known_strand: true,
            references: BTreeMap::new(),
        };

        let mut read_mat = create_scoring_record_3d(
            read_one.len() + 100,
            read_one.len() + 100,
            AlignmentType::Affine,
            false,
        );

        let _my_score = InversionScoring {
            match_score: 9.0,
            mismatch_score: -21.0,
            gap_open: -25.0,
            gap_extend: -1.0,
            inversion_penalty: -40.0,
            min_inversion_length: 20,
        };

        let my_aff_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -9.0,
            special_character_score: 9.0,
            gap_open: -20.0,
            gap_extend: -1.0,
            final_gap_multiplier: 1.0,
        };

        for _i in 0..10 {
            let _ = exhaustive_alignment_search(
                &"testread".to_string(),
                &read_one,
                None,
                &&rm,
                &mut read_mat,
                &my_aff_score,
                None,
            );
        }
    }

    use alignment_functions::{cigar_to_alignment, quick_alignment_search, wfa_alignment};
    use bio::io::{fasta, fastq};
    use flate2::read::GzDecoder;
    use utils::read_utils::u8s;

    #[test]
    fn test_alignment_marc1_data() {
        // Open the FASTA file
        let file = File::open("test_data/all_MARC1_references.fa").unwrap();
        let reader = BufReader::new(file);

        // Create a FASTA reader
        let fasta_reader = fasta::Reader::new(reader);

        let mut indexes = HashMap::new();
        // Iterate over the records in the FASTA file
        for (index, result) in fasta_reader.records().enumerate() {
            let record = result.unwrap();

            // Extract the sequence ID and sequence
            let id = record.id().to_owned();
            let seq = record.seq();
            println!("id {} seq {}", id, String::from_utf8(seq.to_vec()).unwrap());
            indexes.insert(index, id);
        }

        let reference = ReferenceBuilder::new()
            .set_uppercase(true) // Ignore case
            .ignore_base(b'N') // 'N' is never matched
            .add_fasta_file("test_data/all_MARC1_references.fa")
            .unwrap() // Add sequences from FASTA
            .build()
            .unwrap();

        // (2) Initialize `Aligner`
        let algorithm = SemiGlobal::new(
            4,   // Mismatch penalty
            6,   // Gap-open penalty
            2,   // Gap-extend penalty
            50,  // Minimum length
            0.2, // Maximum penalty per length
        )
        .unwrap();

        let mut aligner = Aligner::new(algorithm);

        let input_fasta = "test_data/mouse_24c_known_barcode_with_tags_r1.fq.gz";

        // Create a GzDecoder to decompress the file
        let gz = GzDecoder::new(File::open(input_fasta).unwrap());

        // Use a BufReader to efficiently read the decompressed data
        let reader = BufReader::new(gz);

        // Create a FASTQ reader from the decompressed data
        let fastq_reader = fastq::Reader::new(reader);

        // Iterate over the records in the FASTQ file
        let mut count = 0;
        let mut correct_count = 0;
        for result in fastq_reader.records() {
            count += 1;
            if count < 1000 {
                let record = result.unwrap();

                // Extract the sequence ID, sequence, and quality scores
                let id = record.id();
                let seq = record.seq();

                let correct = id.split('_').next().unwrap_or("failed");

                let result = aligner.align(seq, &reference);
                //println!("{:#?}", result);

                let best_alignment = result
                    .0
                    .iter()
                    .map(|aln| {
                        let mut min_score = u32::MAX;
                        let mut min_pos = 0;
                        let _ = aln.alignments.iter().enumerate().map(|(index, ln)| {
                            if ln.penalty < min_score {
                                min_score = ln.penalty;
                                min_pos = index;
                            }
                        });
                        (min_score, min_pos, aln.index)
                    })
                    .into_iter()
                    .min()
                    .unwrap_or((0, 0, 0));

                let align_name = indexes.get(&(best_alignment.2 as usize)).unwrap();
                let candidate_split = align_name.split('_').collect::<Vec<&str>>()[2];
                if correct == candidate_split {
                    correct_count += 1;
                }
            } else {
                println!("count {} correct {}", count, correct_count);
                return;
            }
        }
    }
}
