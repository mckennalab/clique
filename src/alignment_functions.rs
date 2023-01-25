use std::sync::{Arc, Mutex};
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use std::str;
use crate::itertools::Itertools;
use crate::rayon::iter::ParallelBridge;
use crate::rayon::iter::ParallelIterator;
use crate::alignment::alignment_matrix::{AlignmentTag, reverse_complement};
use crate::alignment::scoring_functions::{AffineScoring, InversionScoring};
use crate::extractor::extract_tagged_sequences;
use crate::linked_alignment::{align_string_with_anchors, AlignmentResults, find_greedy_non_overlapping_segments, orient_by_longest_segment};
use crate::read_strategies::read_set::{ReadIterator, ReadSetContainer};
use crate::reference::fasta_reference::{Reference, ReferenceManager};
use std::time::{Duration, Instant};


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
    let start = Instant::now();
    let mut read_count = Arc::new(Mutex::new(0)); // we rely on the Arc for output as access control

    read_iterator.par_bridge().for_each(|xx| {
        let name = &String::from(xx.read_one.id()).to_string();
        let aligned = best_reference(&xx,
                                     &rm.references,
                                     &my_aff_score,
                                     &my_score,
                                     inversions,
                                     *max_reference_multiplier,
                                     *min_read_length);

        match aligned {
            None => {
                warn!("Unable to find alignment for read {}",xx.read_one.id());
            }
            Some((results,cigar_string, ref_name)) => {

                let output = Arc::clone(&output);
                let mut read_count = read_count.lock().unwrap();
                *read_count += 1;
                if *read_count % 100 == 0 {
                    let duration = start.elapsed();
                    println!("Time elapsed in aligning reads ({:?}) is: {:?}", read_count, duration);
                }


                output_alignment(&results.aligned_read,
                                 name,
                                 &results.aligned_ref,
                                 &ref_name,
                                 use_capture_sequences,
                                 only_output_captured_ref,
                                 to_fake_fastq,
                                 output,
                                 &cigar_string,
                                 min_read_length);
            }
        }
    });
}


pub fn best_reference(read: &ReadSetContainer,
                      references: &Vec<Reference>,
                      my_aff_score: &AffineScoring,
                      my_score: &InversionScoring,
                      use_inversions: &bool,
                      max_reference_multiplier: usize,
                      min_read_length: usize) -> Option<(AlignmentResults, String, Vec<u8>)> {

    match references.len() {
        0 => {
            warn!("Unable to align read {} as it has no candidate references",read.read_one.id());
            None
        }
        1 => {
            let aln = alignment(read, &references[0], my_aff_score, my_score, use_inversions, max_reference_multiplier, min_read_length);
            Some((aln.clone().unwrap().0, aln.clone().unwrap().1, references[0].name.clone()))
        }
        x if x > 1 => {
            let ranked_alignments = references.iter().map(|reference| {
                let aln = alignment(read, reference, my_aff_score, my_score, use_inversions, max_reference_multiplier, min_read_length);
                (aln.clone().unwrap().0, aln.clone().unwrap().1, reference.name.clone())
            });

            let ranked_alignments = ranked_alignments.into_iter().max_by(|al, al2|
                matching_read_bases_prop(&al.0.aligned_read, &al.0.aligned_ref).
                    partial_cmp(&matching_read_bases_prop(&al2.0.aligned_read, &al2.0.aligned_ref)).unwrap());

            ranked_alignments.iter().next().cloned()
        }
        x => { panic!("we dont know what to do with a reference count of {}", x) }
    }
}

pub fn alignment(x: &ReadSetContainer,
                 reference: &Reference,
                 my_aff_score: &AffineScoring,
                 my_score: &InversionScoring,
                 use_inversions: &bool,
                 max_reference_multiplier: usize,
min_read_length: usize ) -> Option<(AlignmentResults, String)> {

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
        let fwd_score_mp = find_greedy_non_overlapping_segments(
            &forward_oriented_seq,
            &reference.sequence,
            &reference.suffix_table);

        let results = align_string_with_anchors(
            &forward_oriented_seq,
            &reference.sequence,
            &fwd_score_mp,
            my_score,
            my_aff_score,
            use_inversions);

        let cigar_string = results.cigar_tags.iter().map(
            |tag| format!("{}", tag)).collect::<Vec<String>>().join(",");

        Some((results, cigar_string))
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

                    if replaced.len() >= *min_read_length {
                        let fake_qual = (0..replaced.len()).map(|_| "H").collect::<String>();
                        write!(output, "@{}_{}_ref_{}\n{}\n+\n{}\n",
                               str::replace(aligned_name, " ", "_"),
                               others,
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
                               others,
                               read_seq,
                        ).expect("Unable to write to output file");
                    } else {
                        warn!("Final read product too short after trimming: {}, dropping read", read_seq.len());
                    }
                }
            }
        };
    } else {
        let collapsed_tags = match extracted_seqs {
            None => { String::from("NONE") }
            Some(x) => { x.iter().map(|k| format!("key={}:{}", &k.0, &k.1)).join(";") }
        };
        write!(output, ">{}\n{}\n>{};{};{}\n{}\n",
               String::from_utf8(reference_name.clone()).unwrap(),
               str::from_utf8(&aligned_ref).unwrap(),
               str::replace(aligned_name, " ", "_"),
               collapsed_tags,
               cigar_string,
               str::from_utf8(&aligned_read).unwrap(),
        ).expect("Unable to write to output file");
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
