use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use std::str;
use crate::itertools::Itertools;
use crate::rayon::iter::ParallelBridge;
use crate::rayon::iter::ParallelIterator;
use crate::alignment::alignment_matrix::{Alignment, AlignmentResult, AlignmentTag, AlignmentType, create_scoring_record_3d, perform_3d_global_traceback, perform_affine_alignment, reverse_complement};
use crate::alignment::scoring_functions::{AffineScoring, AffineScoringFunction, InversionScoring};
use crate::extractor::{extract_tagged_sequences, READ_CHAR, REFERENCE_CHAR};
use crate::linked_alignment::{align_string_with_anchors, AlignmentResults, find_greedy_non_overlapping_segments, orient_by_longest_segment};
use crate::read_strategies::read_set::{ReadIterator, ReadSetContainer};
use crate::reference::fasta_reference::{Reference, ReferenceManager};
use std::time::{Instant};
use ndarray::Ix3;


pub fn align_reads(use_capture_sequences: &bool,
                   only_output_captured_ref: &bool,
                   to_fake_fastq: &bool,
                   reference: &String,
                   output: &String,
                   max_reference_multiplier: &usize,
                   min_read_length: &usize,
                   read1: &String,
                   read2: &String,
                   index1: &String,
                   index2: &String,
                   threads: &usize,
                   inversions: &bool) {

    //let reference = reference_file_to_structs(reference, 20);
    let rm = ReferenceManager::from(&reference, 8);

    let output_file = File::create(&output).unwrap();

    let read_iterator = ReadIterator::new(PathBuf::from(&read1),
                                          Some(PathBuf::from(&read2)),
                                          Some(PathBuf::from(&index1)),
                                          Some(PathBuf::from(&index2)));

    let output = Arc::new(Mutex::new(output_file)); // Mutex::new(gz));

    // setup our thread pool
    rayon::ThreadPoolBuilder::new().num_threads(*threads).build_global().unwrap();

    let my_score = InversionScoring {
        match_score: 9.0,
        mismatch_score: -21.0,
        special_character_score: 8.0,
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

    lazy_static! {
        static ref STORE_CLONES: Mutex<Vec<SharedStore>> = Mutex::new(Vec::new());
    }
    thread_local!(static STORE: SharedStore = Arc::new(Mutex::new(None)));

    let alignment_mat: Alignment<Ix3> = create_scoring_record_3d((rm.longest_ref + 1) * 2, (rm.longest_ref + 1) * 2, AlignmentType::AFFINE, false);

    read_iterator.par_bridge().for_each(|xx| {
        STORE.with(|arc_mtx| {
            let mut local_alignment = arc_mtx.lock().unwrap();
            if local_alignment.is_none() {
                *local_alignment = Some(alignment_mat.clone());
                STORE_CLONES.lock().unwrap().push(arc_mtx.clone());
            }

            let name = &String::from(xx.read_one.id()).to_string();

            let aligned = best_reference(&xx,
                                         &rm.references,
                                         &mut local_alignment.as_mut().unwrap(),
                                         &my_aff_score,
                                         &my_score,
                                         inversions,
                                         *max_reference_multiplier,
                                         *min_read_length);

            match aligned {
                None => {
                    warn!("Unable to find alignment for read {}",xx.read_one.id());
                }
                Some((results, ref_name)) => {
                    match results {
                        None => { warn!("Unable to find alignment for read {}",xx.read_one.id()); }
                        Some(aln) => {
                            let output = Arc::clone(&output);
                            let mut read_count = read_count.lock().unwrap();
                            *read_count += 1;
                            if *read_count % 10000 == 0 {
                                let duration = start.elapsed();
                                println!("Time elapsed in aligning reads ({:?}) is: {:?}", read_count, duration);
                            }

                            let cigar_string = format!("{}", "");

                            output_alignment(&aln.alignment_string2,
                                             name,
                                             &aln.alignment_string1,
                                             &ref_name,
                                             use_capture_sequences,
                                             only_output_captured_ref,
                                             to_fake_fastq,
                                             output,
                                             &cigar_string,
                                             min_read_length);
                        }
                    }
                }
            }
        });
    });
}

pub fn align_two_strings(read1_seq: &Vec<u8>, rev_comp_read2: &Vec<u8>, scoring_function: &dyn AffineScoringFunction, local: bool) -> AlignmentResult {
    let mut alignment_mat = create_scoring_record_3d(
        read1_seq.len() + 1,
        rev_comp_read2.len() + 1,
        AlignmentType::AFFINE,
        local);

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

pub fn best_reference(read: &ReadSetContainer,
                      references: &Vec<Reference>,
                      alignment_mat: &mut Alignment<Ix3>,
                      my_aff_score: &AffineScoring,
                      my_score: &InversionScoring,
                      use_inversions: &bool,
                      max_reference_multiplier: usize,
                      min_read_length: usize) -> Option<(Option<AlignmentResult>, Vec<u8>)> {
    match references.len() {
        0 => {
            warn!("Unable to align read {} as it has no candidate references",read.read_one.id());
            None
        }
        1 => {
            let aln = alignment(read, &references[0], alignment_mat, my_aff_score, my_score, use_inversions, max_reference_multiplier, min_read_length);
            Some((aln, references[0].name.clone()))
        }
        x if x > 1 => {
            let ranked_alignments = references.iter().map(|reference| {
                match alignment(read, reference, alignment_mat, my_aff_score, my_score, use_inversions, max_reference_multiplier, min_read_length) {
                    None => None,
                    Some(aln) => Some((aln, reference.name.clone())),
                }
            }).filter(|x| x.is_some()).map(|c| c.unwrap());

            let ranked_alignments = ranked_alignments.into_iter().max_by(|al, al2|
                matching_read_bases_prop(&al.0.alignment_string2, &al.0.alignment_string1).
                    partial_cmp(&matching_read_bases_prop(&al2.0.alignment_string2, &al2.0.alignment_string1)).unwrap());

            match ranked_alignments.iter().next() {
                None => { None }
                Some((x, y)) => {
                    Some((Some(x.clone()), y.clone()))
                }
            }
        }
        x => { panic!("we dont know what to do with a reference count of {}", x) }
    }
}

pub fn alignment(x: &ReadSetContainer,
                 reference: &Reference,
                 alignment_mat: &mut Alignment<Ix3>,
                 my_aff_score: &AffineScoring,
                 my_score: &InversionScoring,
                 use_inversions: &bool,
                 max_reference_multiplier: usize,
                 min_read_length: usize) -> Option<(AlignmentResult)> {

    // find the best reference(s)
    let orientation = orient_by_longest_segment(&x.read_one.seq().to_vec(), &reference.sequence, &reference.suffix_table).0;
    let forward_oriented_seq = if orientation {
        x.read_one.seq().to_vec()
    } else {
        reverse_complement(&x.read_one.seq().to_vec())
    };

    if forward_oriented_seq.len() > reference.sequence.len() * max_reference_multiplier || forward_oriented_seq.len() < min_read_length {
        warn!("Dropping read of length {}",forward_oriented_seq.len());
        None
    } else {
        perform_affine_alignment(
            alignment_mat,
            &reference.sequence,
            &forward_oriented_seq,
            my_aff_score);

        let results = perform_3d_global_traceback(
            alignment_mat,
            None,
            &reference.sequence,
            &forward_oriented_seq,
            None);

        Some(results)
    }
}

pub fn matching_read_bases_prop(read: &Vec<u8>, reference: &Vec<u8>) -> f32 {
    assert_eq!(read.len(), reference.len());
    let mut total_read_bases = 0;
    let mut matched = 0;
    vec_to_uppercase(read).iter().zip(vec_to_uppercase(reference).iter()).for_each(|(readb, refb)| {
        if *readb != b'-' {
            total_read_bases += 1;
        }
        if *readb != b'-' && readb == refb {
            matched += 1;
        }
    });
    if total_read_bases == 0 {
        0.0
    } else {
        matched as f32 / total_read_bases as f32
    }
}

pub fn vec_to_uppercase(inp: &Vec<u8>) -> Vec<u8> {
    inp.iter().map(|x| x.to_ascii_uppercase()).collect::<Vec<u8>>()
}

pub fn tags_to_output(tags: &BTreeMap<u8, String>) -> String {
    tags.iter().filter(|(k, v)| **k != READ_CHAR && **k != REFERENCE_CHAR).map(|(k, v)| format!("key={}:{}", k, v)).join(";")
}

/// Handle output of the alignment results based on the requested output formats
///
pub fn output_alignment(aligned_read: &Vec<u8>,
                        aligned_name: &String,
                        aligned_ref: &Vec<u8>,
                        reference_name: &Vec<u8>,
                        use_capture_sequences: &bool,
                        only_output_captured_ref: &bool,
                        to_fake_fastq: &bool,
                        output: Arc<Mutex<File>>,
                        cigar_string: &String,
                        min_read_length: &usize) {
    let mut output = output.lock().unwrap();

    let (ref_seq, read_seq, read_tags) =
        if *use_capture_sequences {
            let ets = extract_tagged_sequences(aligned_read, aligned_ref);
            let read_tags = tags_to_output(&ets);
            if *only_output_captured_ref {
                (ets.get(&REFERENCE_CHAR).unwrap().clone(), ets.get(&READ_CHAR).unwrap().clone(), read_tags)
            } else {
                (String::from_utf8(aligned_ref.clone()).unwrap(), String::from_utf8(aligned_read.clone()).unwrap(), read_tags)
            }
        } else {
            (String::from_utf8(aligned_ref.clone()).unwrap(), String::from_utf8(aligned_read.clone()).unwrap(), String::from(""))
        };

    if *to_fake_fastq {
        let replaced = read_seq.replace("-", "");

        if replaced.len() >= *min_read_length {
            let fake_qual = (0..replaced.len()).map(|_| "H").collect::<String>();
            write!(output, "@{}_{}_ref_{}\n{}\n+\n{}\n",
                   str::replace(aligned_name, " ", "_"),
                   read_tags,
                   String::from_utf8(reference_name.clone()).unwrap(),
                   replaced,
                   fake_qual,
            ).expect("Unable to write to output file");
        } else {
            warn!("Final read product too short after trimming: {}, dropping read", replaced.len());
        }
    } else {
        if read_seq.len() >= *min_read_length {
            write!(output, ">{}\n{}\n>{}_{}\n{}\n",
                   String::from_utf8(reference_name.clone()).unwrap(),
                   ref_seq,
                   str::replace(aligned_name, " ", "_"),
                   read_tags,
                   read_seq,
            ).expect("Unable to write to output file");
        } else {
            warn!("Final read product too short after trimming: {}, dropping read", read_seq.len());
        }
    }
}

pub fn simplify_cigar_string(cigar_tokens: &Vec<AlignmentTag>) -> Vec<AlignmentTag> {
    let mut new_cigar = Vec::new();

    let mut last_token: Option<AlignmentTag> = None; // zero length, so combining won't affect the final cigar string

    for token in cigar_tokens.into_iter() {
        last_token = match token {
            AlignmentTag::MatchMismatch(size) => {
                match last_token {
                    Some(AlignmentTag::MatchMismatch(size_old)) => Some(AlignmentTag::MatchMismatch(size + size_old)),
                    _ => {
                        if last_token != None { new_cigar.push(last_token.unwrap()); }

                        Some(AlignmentTag::MatchMismatch(*size))
                    }
                }
            }
            AlignmentTag::Del(size) => {
                match last_token {
                    Some(AlignmentTag::Del(size_old)) => Some(AlignmentTag::Del(size + size_old)),
                    _ => {
                        if last_token != None { new_cigar.push(last_token.unwrap()); }
                        Some(AlignmentTag::Del(*size))
                    }
                }
            }
            AlignmentTag::Ins(size) => {
                match last_token {
                    Some(AlignmentTag::Ins(size_old)) => Some(AlignmentTag::Ins(size + size_old)),
                    _ => {
                        if last_token != None { new_cigar.push(last_token.unwrap()); }
                        Some(AlignmentTag::Ins(*size))
                    }
                }
            }
            AlignmentTag::InversionOpen => {
                match last_token {
                    _ => {
                        // we don't combine inversions
                        if last_token != None { new_cigar.push(last_token.unwrap()); }
                        Some(AlignmentTag::InversionOpen)
                    }
                }
            }
            AlignmentTag::InversionClose => {
                match last_token {
                    _ => {
                        // we don't combine inversions
                        if last_token != None { new_cigar.push(last_token.unwrap()); }
                        Some(AlignmentTag::InversionClose)
                    }
                }
            }
        }
    }
    if let Some(x) = last_token {
        new_cigar.push(x);
    }
    new_cigar.reverse();
    new_cigar
}
