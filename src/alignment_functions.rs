use std::sync::{Arc, Mutex};
use ::std::io::Result;
use std::borrow::Borrow;
use std::collections::{BTreeMap, HashMap};
use std::fs::{File, read};
use std::io::{BufRead, BufReader, Error, Write};
use std::io;
use std::path::{Path, PathBuf};
use std::rc::Rc;
use std::str;
use std::str::FromStr;
use crate::itertools::Itertools;
use crate::rayon::iter::ParallelBridge;
use crate::rayon::iter::ParallelIterator;
use crate::alignment::alignment_matrix::{find_seeds, reverse_complement};
use crate::alignment::scoring_functions::{AffineScoring, InversionScoring};
use crate::extractor::extract_tagged_sequences;
use crate::linked_alignment::{align_string_with_anchors, find_greedy_non_overlapping_segments, orient_by_longest_segment};
use crate::read_strategies::read_set::ReadIterator;
use crate::reference::fasta_reference::reference_file_to_struct;

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

    let reference = reference_file_to_struct(reference);

    let output_file = File::create(&output).unwrap();

    let read_iterator = ReadIterator::new(PathBuf::from(&read1),
                                          Some(PathBuf::from(&read2)),
                                          Some(PathBuf::from(&index1)),
                                          Some(PathBuf::from(&index2)));

    let reference_lookup = find_seeds(&reference.sequence, 20);
    let output = Arc::new(Mutex::new(output_file)); // Mutex::new(gz));

    // setup our thread pool
    rayon::ThreadPoolBuilder::new().num_threads(*threads).build_global().unwrap();

    let my_score = InversionScoring {
        match_score: 9.0,
        mismatch_score: -21.0,
        special_character_score: 9.0,
        gap_open: -25.0,
        gap_extend: -5.0,
        inversion_penalty: -40.0,
        min_inversion_length: 20,
    };

    let my_aff_score = AffineScoring {
        match_score: 10.0,
        mismatch_score: -11.0,
        special_character_score: 10.0,
        gap_open: -20.0,
        gap_extend: -5.0,
    };

    read_iterator.par_bridge().for_each(|xx| {
        let x = xx.read_one;
        let name = &String::from(x.id()).to_string();

        let read_length = x.seq().to_vec().len();
        if read_length > reference.sequence.len() * max_reference_multiplier || read_length < *min_read_length {
            info!("Dropping read of length {}",read_length);
        } else {
            let orientation = orient_by_longest_segment(&x.seq().to_vec(), &reference.sequence, &reference_lookup).0;
            let forward_oriented_seq = if orientation {
                x.seq().to_vec()
            } else {
                reverse_complement(&x.seq().to_vec())
            };

            let fwd_score_mp = find_greedy_non_overlapping_segments(&forward_oriented_seq, &reference.sequence, &reference_lookup);

            let results = align_string_with_anchors(&forward_oriented_seq, &reference.sequence, &fwd_score_mp, &my_score, &my_aff_score, inversions);

            let cigar_string = results.cigar_tags.iter().map(|tag| format!("{}", tag)).collect::<Vec<String>>().join(",");

            let output = Arc::clone(&output);

            output_alignment(&results.aligned_read,
                             name,
                             &results.aligned_ref,
                             use_capture_sequences,
                             only_output_captured_ref,
                             to_fake_fastq,
                             output,
                             &cigar_string);
        }
    });
}


/// Handle output of the alignment results based on the requested output formats
///
pub fn output_alignment(aligned_read: &Vec<u8>,
                        aligned_name: &String,
                        aligned_ref: &Vec<u8>,
                        use_capture_sequences: &bool,
                        only_output_captured_ref: &bool,
                        to_fake_fastq: &bool,
                        output: Arc<Mutex<File>>,
                        cigar_string: &String) {
    let mut output = output.lock().unwrap();

    let extracted_seqs = if *use_capture_sequences {
        Some(extract_tagged_sequences(aligned_read, aligned_ref))
    } else {
        None
    };

    if *only_output_captured_ref {
        match extracted_seqs {
            None => {
                warn!("unable to extract sequences from read {}",aligned_name)
            }
            Some(btree) => {
                let read_seq = btree.iter().map(|kv| {
                    if kv.0.starts_with('r') {
                        kv.1.clone()
                    } else {
                        String::from("")
                    }
                }).join("");
                let ref_seq = btree.iter().map(|kv| {
                    if kv.0.starts_with('e') {
                        kv.1.clone()
                    } else {
                        String::from("")
                    }
                }).join("");
                let others = btree.iter().
                    filter(|k| !k.0.starts_with('e') && !k.0.starts_with('r')).
                    map(|k| format!("key={}:{}", &k.0, &k.1)).join(";");

                if *to_fake_fastq {
                    let replaced = read_seq.replace("-", "");
                    let fake_qual = (0..replaced.len()).map(|_| "H").collect::<String>();
                    write!(output, "@{}_{}\n{}\n+\n{}\n",
                           str::replace(aligned_name, " ", "_"),
                           others,
                           replaced,
                           fake_qual,
                    ).expect("Unable to write to output file");
                } else {
                    write!(output, ">ref\n{}\n>{}_{}\n{}\n",
                           ref_seq,
                           str::replace(aligned_name, " ", "_"),
                           others,
                           read_seq,
                    ).expect("Unable to write to output file");
                }
            }
        };
    } else {
        let collapsed_tags = match extracted_seqs {
            None => { String::from("NONE") }
            Some(x) => { x.iter().map(|k| format!("key={}:{}", &k.0, &k.1)).join(";") }
        };
        write!(output, ">ref\n{}\n>{};{};{}\n{}\n",
               str::from_utf8(&aligned_ref).unwrap(),
               str::replace(aligned_name, " ", "_"),
               collapsed_tags,
               cigar_string,
               str::from_utf8(&aligned_read).unwrap(),
        ).expect("Unable to write to output file");
    }
}
