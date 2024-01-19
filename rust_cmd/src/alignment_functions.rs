use std::collections::{BTreeMap, HashMap, VecDeque};
use std::sync::{Arc, Mutex};
use std::path::{Path, PathBuf};
use std::str;
use crate::rayon::iter::ParallelBridge;
use crate::rayon::iter::ParallelIterator;
use crate::alignment::alignment_matrix::{Alignment, AlignmentResult, AlignmentTag, AlignmentType, create_scoring_record_3d, perform_3d_global_traceback, perform_affine_alignment};
use crate::alignment::scoring_functions::{AffineScoring, AffineScoringFunction, InversionScoring};
use crate::extractor::{extract_tagged_sequences, gap_proportion_per_tag, stretch_sequence_to_alignment};
use crate::linked_alignment::{align_string_with_anchors, find_greedy_non_overlapping_segments, orient_by_longest_segment};
use crate::read_strategies::read_set::{ReadIterator};
use crate::reference::fasta_reference::{Reference, ReferenceManager};
use std::time::{Instant};
use bio::alignment::AlignmentOperation;
use ndarray::Ix3;
use shardio::{ShardReader, ShardSender, ShardWriter};
use crate::alignment::fasta_bit_encoding::{FASTA_UNSET, FastaBase, reverse_complement};
use crate::merger::{MergedReadSequence, UnifiedRead};
use crate::read_strategies::sequence_layout::{SequenceLayoutDesign, UMIConfiguration};
use rust_htslib::bam;
use rust_htslib::bam::{Record};
use rust_htslib::bam::record::{Aux, CigarString};
use crate::read_strategies::read_disk_sorter::{SortedAlignment, SortingReadSetContainer};
use bio::alignment::pairwise::{Aligner, Scoring};
use itertools::Itertools;


pub fn align_reads(read_structure: &SequenceLayoutDesign,
                   rm: &ReferenceManager,
                   output: &Path,
                   max_reference_multiplier: &usize,
                   min_read_length: &usize,
                   read1: &String,
                   read2: &String,
                   index1: &String,
                   index2: &String,
                   threads: &usize,
                   inversions: &bool) {
    let (_reference_mapper, output_file) = setup_sam_writer(&output.to_str().unwrap().to_string(), rm);
    let output_file = output_file.expect("Unable to create output bam file");

    let read_iterator = ReadIterator::new(PathBuf::from(&read1),
                                          Some(PathBuf::from(&read2)),
                                          Some(PathBuf::from(&index1)),
                                          Some(PathBuf::from(&index2)));

    let read_iterator = MergedReadSequence::new(read_iterator, read_structure);

    let output = Arc::new(Mutex::new(output_file));

    // setup our thread pool
    rayon::ThreadPoolBuilder::new().num_threads(*threads).build_global().unwrap();

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
        gap_extend: -1.0,
        final_gap_multiplier: 1.0,

    };
    let start = Instant::now();
    let read_count = Arc::new(Mutex::new(0)); // we rely on this Arc for output file access control

    type SharedStore = Arc<Mutex<Option<Alignment<Ix3>>>>;

    lazy_static! {static ref STORE_CLONES: Mutex<Vec<SharedStore>> = Mutex::new(Vec::new());}
    thread_local!(static STORE: SharedStore = Arc::new(Mutex::new(None)));

    let max_read_size = (rm.longest_ref + 1) * 2;
    let alignment_mat: Alignment<Ix3> = create_scoring_record_3d(rm.longest_ref + 1, max_read_size, AlignmentType::AFFINE, false);

    read_iterator.par_bridge().for_each(|mut xx: UnifiedRead| {
        STORE.with(|arc_mtx| {
            let mut local_alignment = arc_mtx.lock().unwrap();
            if local_alignment.is_none() {
                *local_alignment = Some(alignment_mat.clone());
                STORE_CLONES.lock().unwrap().push(arc_mtx.clone());
            }

            let name = &String::from_utf8(xx.name().clone()).unwrap();

            if xx.seq().len() < max_read_size {
                let aligned = align_to_reference_choices(&xx.seq(),
                                                         &rm,
                                                         &true,
                                                         read_structure,
                                                         &mut local_alignment.as_mut().unwrap(),
                                                         &my_aff_score,
                                                         &my_score,
                                                         inversions,
                                                         *max_reference_multiplier as f64,
                                                         *min_read_length);

                match aligned {
                    None => {
                        // TODO: we should track this and provide a final summary
                        debug!("Unable to create alignment for read {}",name);
                    }
                    Some((results, orig_ref_seq, _ref_name)) => {
                        match results {
                            None => {
                                // TODO: we should track this and provide a final summary

                                debug!("Unable to create alignment for read {}",name);
                            }
                            Some(aln) => {
                                let output = Arc::clone(&output);
                                let mut read_count = read_count.lock().unwrap();
                                *read_count += 1;
                                if *read_count % 1000000 == 0 {
                                    let duration = start.elapsed();
                                    info!("Time elapsed in aligning reads ({:?}) is: {:?}", read_count, duration);
                                }
                                assert_eq!(aln.reference_aligned.len(), aln.read_aligned.len());

                                let samrecord = aln.to_sam_record(&name, &orig_ref_seq, &true);

                                let mut output = output.lock().unwrap();
                                output.write(&samrecord).expect("Unable to write read to output bam file");
                            }
                        }
                    }
                }
            } else {
                warn!("Dropped read {} is it's length {} exceeds 2x the reference length {}", String::from_utf8(xx.name().clone()).unwrap(), xx.seq().len(), max_read_size);
            }
        });
    });
}

type SharedStore = Arc<Mutex<Option<Alignment<Ix3>>>>;

pub fn fast_align_reads(_use_capture_sequences: &bool,
                        read_structure: &SequenceLayoutDesign,
                        _only_output_captured_ref: &bool,
                        fast_lookup: &bool,
                        rm: &ReferenceManager,
                        output: &Path,
                        max_reference_multiplier: &f64,
                        min_read_length: &usize,
                        read1: &String,
                        read2: &String,
                        index1: &String,
                        index2: &String,
                        max_gaps_proportion: &f64,
                        threads: &usize,
                        inversions: &bool) -> (usize, HashMap<String, ShardReader<SortingReadSetContainer>>) {
    let read_iterator = MergedReadSequence::new(ReadIterator::new(PathBuf::from(&read1),
                                                                  Some(PathBuf::from(&read2)),
                                                                  Some(PathBuf::from(&index1)),
                                                                  Some(PathBuf::from(&index2))), read_structure);


    let mut output_files: HashMap<String, String> = HashMap::new();

    let sharded_outputs = read_structure.references.iter().map(|(ref_name, _ref_obj)| {
        let mut output_path = output.to_str().unwrap().to_string().clone();
        output_path.push_str(ref_name.as_str());
        output_files.insert(ref_name.clone(), output_path.clone());
        (ref_name.clone(), ShardWriter::new(Path::new(output_path.as_str()), 32,
                                            256,
                                            1 << 16).unwrap())
    }).collect::<HashMap<String, ShardWriter<SortingReadSetContainer>>>();

    let alignment_mat: Alignment<Ix3> = create_scoring_record_3d((rm.longest_ref + 1) * 2, (rm.longest_ref + 1) * 2, AlignmentType::AFFINE, false);

    // keep track of timing for progress output
    let start = Instant::now();

    let read_count = Arc::new(Mutex::new(0 as usize));
    let skipped_count = Arc::new(Mutex::new(0 as usize));
    let gap_rejected = Arc::new(Mutex::new(0 as usize));


    // setup local thread-safe storage for our alignment matricies
    lazy_static! {
        static ref STORE_CLONES: Mutex<Vec<SharedStore>> = Mutex::new(Vec::new());
        }
    thread_local!(static STORE: SharedStore = Arc::new(Mutex::new(None)));

    {
        let sender = Arc::new(
            Mutex::new(
                sharded_outputs.iter().map(|(name, sharder)| { (name.clone(), sharder.get_sender()) }).collect::<HashMap<String, ShardSender<SortingReadSetContainer>>>()
            ));

        // setup our thread pool
        rayon::ThreadPoolBuilder::new().num_threads(*threads).build_global().unwrap();


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
            gap_extend: -1.0,
            final_gap_multiplier: 1.0,

        };


        read_iterator.par_bridge().for_each(|mut xx| {
            STORE.with(|arc_mtx| {
                let mut local_alignment = arc_mtx.lock().unwrap();
                if local_alignment.is_none() {
                    *local_alignment = Some(alignment_mat.clone());
                    STORE_CLONES.lock().unwrap().push(arc_mtx.clone());
                }

                let aligned = align_to_reference_choices(&xx.seq(),
                                                         &rm,
                                                         &fast_lookup,
                                                         read_structure,
                                                         &mut local_alignment.as_mut().unwrap(),
                                                         &my_aff_score,
                                                         &my_score,
                                                         inversions,
                                                         *max_reference_multiplier,
                                                         *min_read_length);

                match aligned {
                    None | Some((None, _, _)) => {
                        // TODO: track and output summary
                        debug!("Unable to create alignment for read {}",String::from_utf8(xx.name().clone()).unwrap());
                    }
                    Some((Some(alignment), ref_seq, ref_name)) => {
                        let ref_al = FastaBase::to_vec_u8(&alignment.reference_aligned);
                        let read_al = FastaBase::to_vec_u8(&alignment.read_aligned);

                        let full_ref = stretch_sequence_to_alignment(&ref_al, &ref_seq);
                        let ets = extract_tagged_sequences(&read_al, &full_ref);

                        let gap_proportion = gap_proportion_per_tag(&ets);

                        let sorted_tags = get_sorting_order(read_structure, &String::from_utf8(ref_name.clone()).unwrap());

                        if gap_proportion.len() == 0 || gap_proportion.iter().max_by(|a, b| a.total_cmp(b)).unwrap() <= max_gaps_proportion {
                            let (invalid_read, read_tags_ordered) = extract_tag_sequences(&sorted_tags, ets);

                            if !invalid_read {
                                //println!("READ {} Ref {} ",String::from_utf8(xx.name().clone()).unwrap(), String::from_utf8(ref_name.clone()).unwrap());
                                let new_sorted_read_container = SortingReadSetContainer {
                                    ordered_sorting_keys: vec![], // for future use
                                    ordered_unsorted_keys: read_tags_ordered, // the current unsorted tag collection
                                    aligned_read: SortedAlignment {
                                        aligned_read: alignment.read_aligned.clone(),
                                        aligned_ref: alignment.reference_aligned,
                                        ref_name: String::from_utf8(ref_name.clone()).unwrap(),
                                        read_name: String::from_utf8(xx.name().clone()).unwrap(),
                                        cigar_string: alignment.cigar_string.clone(),
                                        score: alignment.score,
                                    },
                                };
                                let subsender = Arc::clone(&sender);
                                let mut sender_unwrapped = subsender.lock().unwrap();
                                let dispatch: &mut ShardSender<SortingReadSetContainer> = sender_unwrapped.get_mut(new_sorted_read_container.aligned_read.ref_name.as_str()).unwrap();
                                dispatch.send(new_sorted_read_container).unwrap();
                            } else {
                                *skipped_count.lock().unwrap() += 1;
                            }
                        } else {
                            *gap_rejected.lock().unwrap() += 1;
                        }
                    }
                };

                *read_count.lock().unwrap() += 1;
                if *read_count.lock().unwrap() % 1000000 == 0 {
                    let duration = start.elapsed();
                    info!("Time elapsed in aligning reads ({:?}) is: {:?}", read_count.lock().unwrap(), duration);
                }
            });
        });

        Arc::clone(&sender).lock().unwrap().iter_mut().for_each(|(_x, y)| { y.finished().unwrap() })
    }

    sharded_outputs.into_iter().for_each(|(_name, mut sender_obj)| {
        sender_obj.finish().unwrap();
    });

    STORE_CLONES.lock().unwrap().iter_mut().for_each(|x| std::mem::drop(x.lock().unwrap()));

    let final_count = (read_count.lock().unwrap().clone() - skipped_count.lock().unwrap().clone()) - gap_rejected.lock().unwrap().clone();
    info!("Aligned {} reads; {} gap-rejected and {} skipped for being longer than {} their reference size", final_count, gap_rejected.lock().unwrap(), skipped_count.lock().unwrap(),*max_reference_multiplier);

    let outputs = output_files.into_iter().filter(|(ref_name, out)| {
        let exists = Path::new(out).exists();
        if !exists {
            warn!("Dropping reference {} as no reads were aligned to it, file {} exists {}", ref_name, out, exists);
        } else {
            warn!("KEEP reference {} as no reads were aligned to it, file {} exists {}", ref_name, out, exists);
        }
        exists
    }).map(|(ref_name, out)| {
        (ref_name, ShardReader::open(out).unwrap())
    }).collect::<HashMap<String, ShardReader<SortingReadSetContainer>>>();
    (final_count, outputs)
}

fn extract_tag_sequences(sorted_tags: &Vec<char>, ets: BTreeMap<u8, String>) -> (bool, VecDeque<(char, Vec<FastaBase>)>) {
    let mut invalid_read = false;
    let queue = VecDeque::from(sorted_tags.iter().
        map(|x| {
            let ets_hit = ets.get(&x.to_string().as_bytes()[0]);
            match ets_hit {
                Some(e) => {
                    Some((x.clone(), e.as_bytes().iter().map(|f| FastaBase::from(f.clone())).collect::<Vec<FastaBase>>()))
                }
                None => {
                    invalid_read = true;
                    None
                }
            }
        }).filter(|x| x.is_some()).map(|x| x.unwrap()).collect::<Vec<(char, Vec<FastaBase>)>>());
    (invalid_read, queue)
}

fn get_sorting_order(read_structure: &SequenceLayoutDesign, reference_name: &String) -> Vec<char> {
    match read_structure.references.get(reference_name) {
        None => {
            panic!("Unable to find reference {} ", reference_name);
        }
        Some(reference) => {
            let mut sorting_order = reference.umi_configurations.iter().map(|x| x.1.clone()).collect::<Vec<UMIConfiguration>>();
            sorting_order.sort_by(|a, b| a.order.cmp(&b.order));
            let sorted_tags = sorting_order.iter().map(|x| x.symbol).collect::<Vec<char>>();
            sorted_tags
        }
    }
}

#[allow(dead_code)]
pub fn align_using_selected_aligner(read_structure: &SequenceLayoutDesign,
                                    reference_manager: &ReferenceManager,
                                    reference_name: &Vec<u8>,
                                    reference_bases: &Vec<FastaBase>,
                                    xx: &Vec<FastaBase>,
                                    alignment_mat: &mut Alignment<Ix3>) -> AlignmentResult {
    match &read_structure.aligner {
        None => {
            let score = AffineScoring::default();
            align_two_strings_passed_matrix(&reference_bases, &xx, &score, Some(reference_name), Some(reference_manager), alignment_mat)
            //perform_rust_bio_alignment(&reference_bases, &xx) // current default
        }
        Some(x) => {
            match x.as_str() {
                "rustbio" => { perform_rust_bio_alignment(&reference_bases, &xx) }
                //"inversion_aware" => {perform_inversion_aware_alignment(&reference_bases, &xx.seq)}
                //"degenerate" => {perform_degenerate_alignment(&reference_bases, &xx.seq)}
                _ => { panic!("Unknown alignment method {}", x) }
            }
        }
    }
}

pub fn align_two_strings(read1_seq: &Vec<FastaBase>,
                         rev_comp_read2: &Vec<FastaBase>,
                         scoring_function: &dyn AffineScoringFunction,
                         local: bool,
                         ref_name: Option<&Vec<u8>>,
                         reference_manager: Option<&ReferenceManager>) -> AlignmentResult {
    let mut alignment_mat = create_scoring_record_3d(
        read1_seq.len() + 1,
        rev_comp_read2.len() + 1,
        AlignmentType::AFFINE,
        local);


    match (reference_manager, ref_name) {
        (Some(x), Some(y)) => {
            let ref_id = x.reference_name_to_ref.get(y).unwrap();
            let shared_segments = &x.references.get(ref_id).unwrap().suffix_table;

            let ref_seq = FastaBase::to_vec_u8(read1_seq);
            let read_seq = FastaBase::to_vec_u8(rev_comp_read2);

            let shared_segs = find_greedy_non_overlapping_segments(
                &ref_seq,
                &read_seq,
                shared_segments);

            align_string_with_anchors(read1_seq,
                                      rev_comp_read2,
                                      &shared_segs,
                                      None,
                                      scoring_function,
                                      &mut alignment_mat)
        }

        _ => {
            perform_affine_alignment(
                &mut alignment_mat,
                read1_seq,
                rev_comp_read2,
                scoring_function);

            perform_3d_global_traceback(
                &mut alignment_mat,
                None,
                read1_seq,
                rev_comp_read2,
                None)
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
pub fn align_two_strings_passed_matrix(read1_seq: &Vec<FastaBase>,
                                       read2_seq: &Vec<FastaBase>,
                                       scoring_function: &dyn AffineScoringFunction,
                                       ref_name: Option<&Vec<u8>>,
                                       reference_manager: Option<&ReferenceManager>,
                                       alignment_mat: &mut Alignment<Ix3>) -> AlignmentResult {
    match (reference_manager, ref_name) {
        (Some(x), Some(y)) => {
            let ref_id = x.reference_name_to_ref.get(y).unwrap();
            let shared_segments = &x.references.get(ref_id).unwrap().suffix_table;

            let ref_seq = FastaBase::to_vec_u8(read1_seq);
            let read_seq = FastaBase::to_vec_u8(read2_seq);

            let shared_segs = find_greedy_non_overlapping_segments(
                &read_seq,
                &ref_seq,
                shared_segments);

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

        _ => {
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
        }
    }
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
pub fn align_to_reference_choices(read: &Vec<FastaBase>,
                                  rm: &ReferenceManager,
                                  fast_lookup: &bool,
                                  read_structure: &SequenceLayoutDesign,
                                  alignment_mat: &mut Alignment<Ix3>,
                                  my_aff_score: &AffineScoring,
                                  my_score: &InversionScoring,
                                  use_inversions: &bool,
                                  max_reference_multiplier: f64,
                                  min_read_length: usize) -> Option<(Option<AlignmentResult>, Vec<u8>, Vec<u8>)> {
    match rm.references.len() {
        0 => {
            // TODO: we should track this and provide a final summary
            warn!("Unable to align read {} as it has no candidate references",FastaBase::to_string(read));
            None
        }
        1 => {
            let ref_base = &rm.references.get(&0).unwrap();
            let forward_oriented_seq = if !read_structure.known_strand {
                let orientation = orient_by_longest_segment(&read, &ref_base.sequence_u8, &ref_base.suffix_table).0;
                if orientation {
                    read.clone()
                } else {
                    reverse_complement(&read)
                }
            } else {
                read.clone()
            };

            let aln = align_two_strings_passed_matrix(
                &ref_base.sequence,
                &forward_oriented_seq,
                my_aff_score,
                Some(&ref_base.name),
                Some(rm),
                alignment_mat);

            Some(
                (Some(aln),
                 ref_base.sequence_u8.clone(),
                 ref_base.name.clone()))
        }
        x if x > 1 => {
            if *fast_lookup {
                quick_alignment_search(read, &rm, read_structure, alignment_mat, my_aff_score, my_score, use_inversions, max_reference_multiplier, min_read_length)
            } else {
                exhaustive_alignment_search(read, &rm, alignment_mat, my_aff_score)
            }
        }
        x => { panic!("we dont know what to do with a reference count of {}", x) }
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
fn quick_alignment_search(read: &Vec<FastaBase>,
                          rm: &ReferenceManager,
                          _read_structure: &SequenceLayoutDesign,
                          alignment_mat: &mut Alignment<Ix3>,
                          my_aff_score: &AffineScoring,
                          _my_score: &InversionScoring,
                          _use_inversions: &bool,
                          _max_reference_multiplier: f64,
                          _min_read_length: usize) -> Option<(Option<AlignmentResult>, Vec<u8>, Vec<u8>)> {
    let read_u8 = FastaBase::to_vec_u8(read);
    let read_kmers = ReferenceManager::sequence_to_kmers(&read_u8, &rm.kmer_size, &rm.kmer_skip);

    let max_ref = read_kmers.iter().map(|(kmer, _c)| {
        rm.unique_kmers.kmer_to_reference.get(kmer)
    }).flatten().counts();
    let max_ref = max_ref.iter().max_by(|x, y| x.1.cmp(y.1));

    match max_ref {
        None => {
            None
        }
        Some(x) => {
            Some((Some(align_two_strings_passed_matrix(&x.0.sequence,
                                                       read,
                                                       my_aff_score,
                                                       Some(&x.0.name),
                                                       Some(rm),
                                                       alignment_mat)), x.0.sequence_u8.clone(), x.0.name.clone()))
        }
    }
}

fn exhaustive_alignment_search(read: &Vec<FastaBase>,
                               rm: &ReferenceManager,
                               alignment_mat: &mut Alignment<Ix3>,
                               my_aff_score: &AffineScoring) -> Option<(Option<AlignmentResult>, Vec<u8>, Vec<u8>)> {

    let references = &rm.references;

    let ranked_alignments = references.iter().map(|reference| {
        let lt = align_two_strings_passed_matrix(&reference.1.sequence, read, my_aff_score, Some(&reference.1.name), Some(rm), alignment_mat);
        //println!(">{}\n{}\n{}\n{}\n",String::from_utf8(reference.1.sequence_u8.clone()).unwrap(),String::from_utf8(reference.1.name.clone()).unwrap(),FastaBase::to_string(&lt.reference_aligned),FastaBase::to_string(&lt.read_aligned));

        Some((lt, reference.1.sequence_u8.clone(), reference.1.name.clone()))
    }).filter(|x| x.is_some()).map(|c| c.unwrap());


    //println!("read {} ", FastaBase::to_string(read));
    let ranked_alignments = ranked_alignments.into_iter().enumerate().max_by(|al, al2| {
        let score1 = al.1.0.score;// / al.1.0.reference_aligned.len() as f64;
        let score2 = al2.1.0.score;// / al2.1.0.reference_aligned.len() as f64;
        //println!("Score 1 {}, al {} score 2 {} al2 {}", score1, String::from_utf8(al.1.2.clone()).unwrap(), score2, String::from_utf8(al2.1.2.clone()).unwrap());
        score1.partial_cmp(&score2).unwrap()
    });

    match ranked_alignments.iter().next() {
        None => { None }
        Some((_x, y)) => {
            Some((Some(y.0.clone()), y.1.clone(), y.2.clone()))
        }
    }
}

#[allow(dead_code)]
fn cigar_to_alignment(reference: &Vec<FastaBase>,
                      read: &Vec<FastaBase>,
                      cigar: &Vec<u8>) -> (Vec<FastaBase>, Vec<FastaBase>) {
    let mut alignment_string1 = Vec::new();
    let mut alignment_string2 = Vec::new();
    let mut seq1_index = 0;
    let mut seq2_index = 0;

    for c in cigar {
        match c {
            b'M' | b'X' => {
                alignment_string1.push(reference.get(seq1_index).unwrap().clone());
                alignment_string2.push(read.get(seq2_index).unwrap().clone());
                seq1_index += 1;
                seq2_index += 1;
            }
            b'D' => {
                alignment_string1.push(reference.get(seq1_index).unwrap().clone());
                alignment_string2.push(FASTA_UNSET);
                seq1_index += 1;
            }
            b'I' => {
                alignment_string1.push(FASTA_UNSET);
                alignment_string2.push(read.get(seq2_index).unwrap().clone());
                seq2_index += 1;
            }
            _ => { panic!("Unknown cigar operation {}", c) }
        }
    };
    (alignment_string1, alignment_string2)
}

pub fn perform_rust_bio_alignment(reference: &Vec<FastaBase>, read: &Vec<FastaBase>) -> AlignmentResult {
    // TODO: do a better look at scoring here!
    let score = |a: u8, b: u8| if a == b { 2i32 } else if a == 'N' as u8 || b == 'N' as u8 { 2i32 } else { -2i32 };

    let alignment_scoring = Scoring::new(-20, -2, &score).xclip_prefix(0).xclip_suffix(0).yclip_suffix(0).yclip_suffix(0);
    let mut aligner = Aligner::with_capacity_and_scoring(reference.len(), read.len(), alignment_scoring);
    let x = FastaBase::to_vec_u8(&reference);
    let y = FastaBase::to_vec_u8(&read);
    let alignment = aligner.global(y.as_slice(), x.as_slice());
    bio_to_alignment_result(alignment, reference, read)
}


pub fn bio_to_alignment_result(alignment: bio::alignment::Alignment, reference: &Vec<FastaBase>, read: &Vec<FastaBase>) -> AlignmentResult {
    let mut aligned_ref = Vec::new();
    let mut aligned_read = Vec::new();
    let mut ref_pos = alignment.ystart;
    let mut read_pos = alignment.xstart;

    let mut resulting_cigar = Vec::new();
    for al in alignment.operations {
        match al {
            AlignmentOperation::Match => {
                aligned_ref.push(reference.get(ref_pos).unwrap().clone());
                aligned_read.push(read.get(read_pos).unwrap().clone());
                ref_pos += 1;
                read_pos += 1;
                resulting_cigar.push(AlignmentTag::MatchMismatch(1));
            }
            AlignmentOperation::Subst => {
                aligned_ref.push(reference.get(ref_pos).unwrap().clone());
                aligned_read.push(read.get(read_pos).unwrap().clone());
                ref_pos += 1;
                read_pos += 1;
                resulting_cigar.push(AlignmentTag::MatchMismatch(1));
            }
            AlignmentOperation::Del => {
                aligned_ref.push(reference.get(ref_pos).unwrap().clone());
                aligned_read.push(FASTA_UNSET);
                ref_pos += 1;
                resulting_cigar.push(AlignmentTag::Del(1));
            }
            AlignmentOperation::Ins => {
                aligned_ref.push(FASTA_UNSET);
                aligned_read.push(read.get(read_pos).unwrap().clone());
                read_pos += 1;
                resulting_cigar.push(AlignmentTag::Ins(1));
            }
            AlignmentOperation::Xclip(x) => {
                aligned_ref.extend(reference[ref_pos..ref_pos + x].iter());
                aligned_read.extend(FastaBase::from_vec_u8(&vec![b'-'; x]).iter());
                ref_pos += x.clone();
                resulting_cigar.push(AlignmentTag::Ins(x));
            }
            AlignmentOperation::Yclip(y) => {
                aligned_read.extend(read[read_pos..read_pos + y].iter());
                aligned_ref.extend(FastaBase::from_vec_u8(&vec![b'-'; y]).iter());
                read_pos += y.clone();
                resulting_cigar.push(AlignmentTag::Del(y));
            }
        }
    }
    AlignmentResult {
        reference_aligned: aligned_ref,
        read_aligned: aligned_read,
        cigar_string: simplify_cigar_string(&resulting_cigar),
        path: vec![],
        score: alignment.score as f64,
        reference_start: alignment.xstart,
        read_start: alignment.ystart,
        bounding_box: None,
    }
}


pub fn matching_read_bases_prop(read: &Vec<FastaBase>, reference: &Vec<FastaBase>) -> f32 {
    assert_eq!(read.len(), reference.len());
    let mut total_read_bases = 0;
    let mut matched = 0;
    read.iter().zip(reference.iter()).for_each(|(readb, refb)| {
        if *readb != FASTA_UNSET {
            total_read_bases += 1;
        }
        if *readb != FASTA_UNSET && readb == refb {
            matched += 1;
        }
    });
    if total_read_bases == 0 {
        0.0
    } else {
        matched as f32 / total_read_bases as f32
    }
}

pub fn setup_sam_writer(filename: &String, reference_manger: &ReferenceManager) -> (HashMap<Vec<u8>, u16>, Result<bam::Writer, rust_htslib::errors::Error>) {
    let mut header = bam::Header::new();

    let mut reference_to_bin = HashMap::new();

    reference_manger.references.iter().enumerate().for_each(|(index, reference)| {
        let mut header_record = bam::header::HeaderRecord::new(b"SQ");
        header_record.push_tag(b"SN", String::from_utf8(reference.1.name.clone()).unwrap());
        header_record.push_tag(b"LN", &reference.1.sequence.len());
        header.push_record(&header_record);
        reference_to_bin.insert(reference.1.name.clone(), index as u16);
    });

    (reference_to_bin, bam::Writer::from_path(&filename, &header, bam::Format::Bam))
}


pub fn create_sam_record(
    reference_bin: &u16,
    read_name: &str,
    read_seq: &Vec<FastaBase>,
    reference_aligned_seq: &Vec<FastaBase>,
    reference_original_seq: &Vec<u8>,
    cigar_string: &CigarString,
    extract_capture_tags: &bool,
    additional_tags: HashMap<(u8, u8), String>) -> Record {
    let mut record = Record::new();

    let seq = FastaBase::to_vec_u8_strip_gaps(&read_seq);

    // we don't currently calculate real quality scores
    let quals = vec![255 as u8; seq.len()];

    record.set(read_name.as_bytes(), Some(&cigar_string), &seq.as_slice(), &quals.as_slice());
    record.set_bin(*reference_bin);
    record.set_mate_unmapped();
    additional_tags.iter().for_each(|(x, y)| {
        record.push_aux(vec![x.0, x.1].as_slice(), Aux::String(y)).unwrap();
    });

    if *extract_capture_tags {
        let full_ref = stretch_sequence_to_alignment(&FastaBase::to_vec_u8(&reference_aligned_seq), reference_original_seq);
        let ets = extract_tagged_sequences(&FastaBase::to_vec_u8(&read_seq), &full_ref);

        ets.iter().for_each(|(tag, seq)| {
            record.push_aux(vec![b'e', (*tag)].as_slice(), Aux::String(seq)).unwrap();
        });
    }
    record
}

pub fn simplify_cigar_string(cigar_tokens: &Vec<AlignmentTag>) -> Vec<AlignmentTag> {
    let mut new_cigar = Vec::new();

    let mut last_token: Option<AlignmentTag> = None; // zero length, so combining won't affect the final cigar string

    cigar_tokens.iter().for_each(|token| {
        match (&last_token, token) {
            (None, _) => { last_token = Some(token.clone()) }
            (Some(AlignmentTag::InversionOpen), AlignmentTag::InversionOpen) => {
                panic!("Cannot have two inversion open tags in a row");
            }
            (Some(AlignmentTag::InversionClose), AlignmentTag::InversionClose) => {
                panic!("Cannot have two inversion closed tags in a row");
            }
            (Some(AlignmentTag::MatchMismatch(last_count)), AlignmentTag::MatchMismatch(this_count)) => {
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
        }
    });

    if let Some(x) = last_token {
        new_cigar.push(x);
    }
    new_cigar
}


#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;
    use rust_htslib::bam;
    use crate::alignment::alignment_matrix::{AlignmentResult, AlignmentTag, AlignmentType, create_scoring_record_3d};
    use crate::alignment::fasta_bit_encoding::FastaBase;
    use crate::alignment::scoring_functions::{AffineScoring, InversionScoring};
    use crate::alignment_functions::{exhaustive_alignment_search, simplify_cigar_string};
    use crate::read_strategies::sequence_layout::{AlignedReadOrientation, ReadPosition, SequenceLayoutDesign};
    use crate::reference::fasta_reference::{Reference, reference_sequences_to_structs, ReferenceManager, SuffixTableLookup};

    #[test]
    fn test_find_best_reference() {
        let ref_location = &"test_data/test_best_alignment.fasta".to_string();
        let rm = ReferenceManager::from_fa_file(&ref_location, 8, 8);

        let read_one = FastaBase::from_string(&"atggactatcatatgcttaccgtaacttgaaagtatttcgatttcttggctttatatatcttgtggaaaggacgaaacaccgGGTAGCAAACGTTTGGACGTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTTTCCTGCAGGAAACCCCGGGgaat".to_string().to_ascii_uppercase());

        let read_structure = SequenceLayoutDesign {
            aligner: None,
            merge: None,
            reads: vec![ReadPosition::READ1 { chain_align: None, orientation: AlignedReadOrientation::Forward }],
            known_strand: true,
            references: BTreeMap::new(),
        };

        let mut read_mat = create_scoring_record_3d(read_one.len() + 100, read_one.len() + 100, AlignmentType::AFFINE, false);

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
            gap_extend: -1.0,
            final_gap_multiplier: 1.0,

        };

        let best_ref = exhaustive_alignment_search(&read_one, &&rm, &mut read_mat, &my_aff_score);
        assert_eq!(String::from_utf8(best_ref.unwrap().2).unwrap(),
                   String::from_utf8("1_AAACCCCGGG_GGTAGCAAACGTTTGGACGTG".to_string().into_bytes()).unwrap());

        let read_one = FastaBase::from_string(&"atggactatcatatgcttaccgtaacttgaaagtatttcgatttcttggctttatatatcttgtggaaaggacgaaacaccgGGTGCCCTTACTCTCACCTGATTACTTAATCCGTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTTTCCTGCAGGAACGCCCTACgaattcgggcccattggtatggc".to_string().to_ascii_uppercase());
        let best_ref = exhaustive_alignment_search(&read_one, &&rm, &mut read_mat, &my_aff_score);

        assert_eq!(String::from_utf8(best_ref.unwrap().2).unwrap(),
                   String::from_utf8("2_AACGCCCTAC_GGTGCCCTTACTCTCACCTGATTACTTAATCCGTG".to_string().into_bytes()).unwrap());
    }

    fn str_to_Reference<'a>(st: &'a str, name: &'a str) -> Reference<'a, 'a> {
        Reference{
            sequence: FastaBase::from_str(st),
            sequence_u8: st.as_bytes().to_vec(),
            name: name.to_ascii_uppercase().into_bytes(),
            suffix_table: ReferenceManager::find_seeds(&st.as_bytes().to_vec(), 8),
        }
    }

    #[test]
    fn simplify_cigar_test() {
        let input_cigar = vec![AlignmentTag::MatchMismatch(1), AlignmentTag::MatchMismatch(1), AlignmentTag::MatchMismatch(1)];
        let merged_cigar = vec![AlignmentTag::MatchMismatch(3)];
        let resulting_cigar = simplify_cigar_string(&input_cigar);
        assert_eq!(resulting_cigar, merged_cigar);

        let input_cigar = vec![AlignmentTag::MatchMismatch(1), AlignmentTag::Ins(1), AlignmentTag::MatchMismatch(1), AlignmentTag::MatchMismatch(1)];
        let merged_cigar = vec![AlignmentTag::MatchMismatch(1), AlignmentTag::Ins(1), AlignmentTag::MatchMismatch(2)];
        let resulting_cigar = simplify_cigar_string(&input_cigar);
        assert_eq!(resulting_cigar, merged_cigar);
    }

    #[test]
    fn writing_sam_file() {
        let mut header_record = bam::header::HeaderRecord::new(b"SQ");
        header_record.push_tag(b"SN", &"chr1");
        header_record.push_tag(b"LN", &"150");
        let mut header = bam::Header::new();
        header.push_record(&header_record);

        let mut out = bam::Writer::from_path(&"test_data/out.bam", &header, bam::Format::Bam).unwrap();
        let alignment_record = AlignmentResult {
            reference_aligned: FastaBase::from_string(&"CCAATCTACTACTGCTTGCA".to_string()),
            read_aligned: FastaBase::from_string(&"CCAATCTACTACTGCTTGCA".to_string()),
            cigar_string: vec![AlignmentTag::MatchMismatch(20)],
            path: vec![],
            score: 0.0,
            reference_start: 0,
            read_start: 0,
            bounding_box: None,
        };

        for _i in 0..1000 {
            out.write(&alignment_record.to_sam_record("test",
                                                      &"CCAATCTACTACTGCTTGCA".to_string().into_bytes(),
                                                      &false)).unwrap();
        }
    }
}