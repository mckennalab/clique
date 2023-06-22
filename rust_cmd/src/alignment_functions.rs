use std::collections::{BTreeMap, VecDeque};
use std::sync::{Arc, Mutex};
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Output;
use std::str;
use crate::itertools::Itertools;
use crate::rayon::iter::ParallelBridge;
use crate::rayon::iter::ParallelIterator;
use crate::alignment::alignment_matrix::{Alignment, AlignmentResult, AlignmentTag, AlignmentType, create_scoring_record_3d, perform_3d_global_traceback, perform_affine_alignment};
use crate::alignment::scoring_functions::{AffineScoring, AffineScoringFunction, InversionScoring};
use crate::extractor::{extract_tagged_sequences, READ_CHAR, REFERENCE_CHAR, stretch_sequence_to_alignment};
use crate::linked_alignment::{orient_by_longest_segment};
use crate::read_strategies::read_set::{ReadIterator};
use crate::reference::fasta_reference::{Reference, ReferenceManager};
use std::time::{Instant};
use ndarray::Ix3;
use shardio::{ShardReader, ShardWriter};
use crate::alignment::fasta_bit_encoding::{FASTA_UNSET, FastaBase, reverse_complement};
use crate::merger::MergedReadSequence;
use crate::read_strategies::sequence_layout::{SequenceLayoutDesign, UMIConfiguration};
use libwfa::{affine_wavefront::*, bindings::*, mm_allocator::*, penalties::*};
use crate::read_strategies::read_disk_sorter::{SortedAlignment, SortingReadSetContainer};


pub fn align_reads(use_capture_sequences: &bool,
                   read_structure: &SequenceLayoutDesign,
                   only_output_captured_ref: &bool,
                   to_fake_fastq: &bool,
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
    let output_file = File::create(&output).unwrap();

    let read_iterator = ReadIterator::new(PathBuf::from(&read1),
                                          Some(PathBuf::from(&read2)),
                                          Some(PathBuf::from(&index1)),
                                          Some(PathBuf::from(&index2)));

    let read_iterator = MergedReadSequence::new(read_iterator, read_structure);

    let output = Arc::new(Mutex::new(output_file)); // Mutex::new(gz));

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

            let name = &String::from_utf8(xx.name).unwrap();

            let aligned = best_reference(&xx.seq,
                                         &rm.references,
                                         read_structure,
                                         &mut local_alignment.as_mut().unwrap(),
                                         &my_aff_score,
                                         &my_score,
                                         inversions,
                                         *max_reference_multiplier,
                                         *min_read_length);
            //println!("Aligning read {}", FastaBase::to_string(&xx.seq));
            match aligned {
                None => {
                    warn!("Unable to find alignment for read {}",name);
                }
                Some((results, orig_ref_seq, ref_name)) => {
                    match results {
                        None => { warn!("Unable to find alignment for read {}",name); }
                        Some(aln) => {
                            let output = Arc::clone(&output);
                            let mut read_count = read_count.lock().unwrap();
                            *read_count += 1;
                            if *read_count % 1000000 == 0 {
                                let duration = start.elapsed();
                                //println!("Aligning read {} ({})", FastaBase::to_string(&aln.alignment_string1), FastaBase::to_string(&aln.alignment_string2));
                                println!("Time elapsed in aligning reads ({:?}) is: {:?}", read_count, duration);
                            }
                            assert_eq!(aln.alignment_string1.len(), aln.alignment_string2.len());
                            let cigar_string = format!("{}", "");

                            output_alignment(&FastaBase::to_vec_u8(&aln.alignment_string2),
                                             name,
                                             &FastaBase::to_vec_u8(&aln.alignment_string1),
                                             &orig_ref_seq,
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

pub fn fast_align_reads(use_capture_sequences: &bool,
                        read_structure: &SequenceLayoutDesign,
                        only_output_captured_ref: &bool,
                        to_fake_fastq: &bool,
                        rm: &ReferenceManager,
                        output: &Path,
                        max_reference_multiplier: &usize,
                        min_read_length: &usize,
                        read1: &String,
                        read2: &String,
                        index1: &String,
                        index2: &String,
                        threads: &usize,
                        inversions: &bool) -> (usize,ShardReader<SortingReadSetContainer>) {
    let output_file = File::create(&output).unwrap();

    let read_iterator = ReadIterator::new(PathBuf::from(&read1),
                                          Some(PathBuf::from(&read2)),
                                          Some(PathBuf::from(&index1)),
                                          Some(PathBuf::from(&index2)));

    let read_iterator = MergedReadSequence::new(read_iterator, read_structure);

    let mut sharded_output: ShardWriter<SortingReadSetContainer> = ShardWriter::new(&output, 32,
                                                                                    256,
                                                                                    1 << 16).unwrap();

    let mut sender = Arc::new(Mutex::new(sharded_output.get_sender()));

    let mut sorting_order = read_structure.umi_configurations.iter().map(|x| x.1.clone()).collect::<Vec<UMIConfiguration>>();
    sorting_order.sort_by(|a, b| a.order.cmp(&b.order));
    let sorted_tags = sorting_order.iter().map(|x| x.symbol).collect::<Vec<char>>();

    // setup our thread pool
    rayon::ThreadPoolBuilder::new().num_threads(*threads).build_global().unwrap();

    let start = Instant::now();

    assert_eq!(rm.references.len(), 1);
    let reference_seq = rm.references.iter().next().unwrap();
    let reference_bases = FastaBase::from_vec_u8_default_ns(&reference_seq.sequence_u8);
    let reference_name = String::from_utf8(reference_seq.name.clone()).unwrap();
    let mut read_count = Arc::new(Mutex::new(0 as usize));
    read_iterator.par_bridge().for_each(|xx| {
        let name = &String::from_utf8(xx.name.clone()).unwrap();
        let alignment = perform_wavefront_alignment(&reference_bases, &xx.seq);
        let ref_al = FastaBase::to_vec_u8(&alignment.alignment_string1);
        let read_al = FastaBase::to_vec_u8(&alignment.alignment_string2);
        let full_ref = stretch_sequence_to_alignment(&ref_al, &reference_seq.sequence_u8);
        let ets = extract_tagged_sequences(&read_al, &full_ref);


        let read_tags_ordered = VecDeque::from(sorted_tags.iter().
            map(|x| (x.clone(), ets.get(&x.to_string().as_bytes()[0]).unwrap().as_bytes().iter().map(|f| FastaBase::from(f.clone())).collect::<Vec<FastaBase>>())).collect::<Vec<(char, Vec<FastaBase>)>>());

        let new_sorted_read_container = SortingReadSetContainer {
            ordered_sorting_keys: vec![],
            ordered_unsorted_keys: read_tags_ordered,

            aligned_read: SortedAlignment {
                aligned_read: alignment.alignment_string1,
                aligned_ref: alignment.alignment_string2,
                ref_name: reference_name.clone(),
                read_name: String::from_utf8(xx.name.clone()).unwrap(),
            },
        };
        assert_eq!(new_sorted_read_container.ordered_unsorted_keys.len(), read_structure.umi_configurations.len());

        let sender = Arc::clone(&sender);

        sender.lock().unwrap().send(new_sorted_read_container).unwrap();
        *read_count.lock().unwrap() += 1;
        if *read_count.lock().unwrap() % 1000000 == 0 {
            let duration = start.elapsed();
            println!("Time elapsed in aligning reads ({:?}) is: {:?}", read_count.lock().unwrap(), duration);
        }
        if *read_count.lock().unwrap() < 100 {
            println!("Read: {:?}", xx);
            println!("ref: {:?}", ref_al);
            println!("read: {:?}", read_al);
            println!("Extracted tags: {:?}", ets);
        }
    });

    sender.lock().unwrap().finished().unwrap();
    sharded_output.finish().unwrap();

    let final_count = read_count.lock().unwrap();
    (final_count.clone(),ShardReader::open(output).unwrap())
}

#[allow(dead_code)]
pub fn align_two_strings(read1_seq: &Vec<FastaBase>, rev_comp_read2: &Vec<FastaBase>, scoring_function: &dyn AffineScoringFunction, local: bool) -> AlignmentResult {
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

pub fn best_reference(read: &Vec<FastaBase>,
                      references: &Vec<Reference>,
                      read_structure: &SequenceLayoutDesign,
                      alignment_mat: &mut Alignment<Ix3>,
                      my_aff_score: &AffineScoring,
                      my_score: &InversionScoring,
                      use_inversions: &bool,
                      max_reference_multiplier: usize,
                      min_read_length: usize) -> Option<(Option<AlignmentResult>, Vec<u8>, Vec<u8>)> {
    match references.len() {
        0 => {
            warn!("Unable to align read {} as it has no candidate references",FastaBase::to_string(read));
            None
        }
        1 => {
            let aln = alignment(read, &references[0], read_structure, alignment_mat, my_aff_score, my_score, use_inversions, max_reference_multiplier, min_read_length);
            Some((aln, references[0].sequence_u8.clone(), references[0].name.clone()))
        }
        x if x > 1 => {
            let ranked_alignments = references.iter().map(|reference| {
                match alignment(read, reference, read_structure, alignment_mat, my_aff_score, my_score, use_inversions, max_reference_multiplier, min_read_length) {
                    None => None,
                    Some(aln) => Some((aln, reference.sequence_u8.clone(), reference.name.clone())),
                }
            }).filter(|x| x.is_some()).map(|c| c.unwrap());

            let ranked_alignments = ranked_alignments.into_iter().enumerate().max_by(|al, al2|
                matching_read_bases_prop(&al.1.0.alignment_string2, &al.1.0.alignment_string1).
                    partial_cmp(&matching_read_bases_prop(&al2.1.0.alignment_string2, &al2.1.0.alignment_string1)).unwrap());

            match ranked_alignments.iter().next() {
                None => { None }
                Some((x, y)) => {
                    Some((Some(y.0.clone()), y.1.clone(), y.2.clone()))
                }
            }
        }
        x => { panic!("we dont know what to do with a reference count of {}", x) }
    }
}

// TODO bring back inversions
pub fn alignment(x: &Vec<FastaBase>,
                 reference: &Reference,
                 read_structure: &SequenceLayoutDesign,
                 alignment_mat: &mut Alignment<Ix3>,
                 my_aff_score: &AffineScoring,
                 _my_score: &InversionScoring,
                 _use_inversions: &bool,
                 max_reference_multiplier: usize,
                 min_read_length: usize) -> Option<AlignmentResult> {

    // find the best reference(s)
    let forward_oriented_seq = //if !read_structure.known_orientation {
        /*let orientation = orient_by_longest_segment(x, &reference.sequence_u8, &reference.suffix_table).0;
        if orientation {*/
        x.clone();
    /*} else {
        reverse_complement(&x)
    }
} else {
    x.clone()
};*/

    if forward_oriented_seq.len() > reference.sequence.len() * max_reference_multiplier || forward_oriented_seq.len() < min_read_length {
        warn!("Dropping read of length {}",forward_oriented_seq.len());
        None
    } else {
        /*perform_affine_alignment(
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
        Some(results)*/
        Some(perform_wavefront_alignment(&reference.sequence,
                                         &forward_oriented_seq))
    }
}

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

fn perform_wavefront_alignment(reference: &Vec<FastaBase>,
                               read: &Vec<FastaBase>) -> AlignmentResult {
    let alloc = MMAllocator::new(BUFFER_SIZE_8M as u64);

    let mut penalties = AffinePenalties {
        match_: 0,
        mismatch: 4,
        gap_opening: 6,
        gap_extension: 2,
    };

    let pat_len = read.len();
    let text_len = reference.len();

    let mut wavefronts = AffineWavefronts::new_complete(
        pat_len,
        text_len,
        &mut penalties,
        &alloc,
    );

    wavefronts
        .align(FastaBase::to_vec_u8(reference).as_slice(), FastaBase::to_vec_u8(read).as_slice())
        .unwrap();

    let cigar = wavefronts.cigar_bytes_raw();
    let score = wavefronts.edit_cigar_score(&mut penalties);
    let ref_read = cigar_to_alignment(reference, read, &cigar);


    AlignmentResult {
        alignment_string1: ref_read.0,
        alignment_string2: ref_read.1,
        cigar_string: cigar.iter().map(|x| AlignmentTag::from(x.clone())).collect::<Vec<AlignmentTag>>(),
        path: vec![],
        score: score as f64,
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

pub fn tags_to_output(tags: &BTreeMap<u8, String>) -> String {
    tags.iter().filter(|(k, _v)| **k != READ_CHAR && **k != REFERENCE_CHAR).map(|(k, v)| format!("key={}:{}", String::from_utf8(vec![*k]).unwrap(), v)).join(";")
}


/// Handle output of the alignment results based on the requested output formats
///
pub fn output_alignment(aligned_read: &Vec<u8>,
                        aligned_name: &String,
                        aligned_ref: &Vec<u8>,
                        original_ref: &Vec<u8>,
                        reference_name: &Vec<u8>,
                        use_capture_sequences: &bool,
                        only_output_captured_ref: &bool,
                        to_fake_fastq: &bool,
                        output: Arc<Mutex<File>>,
                        _cigar_string: &String,
                        min_read_length: &usize) {
    let mut output = output.lock().unwrap();

    let (ref_seq, read_seq, read_tags) =
        if *use_capture_sequences {
            let full_ref = stretch_sequence_to_alignment(aligned_ref, original_ref);
            let ets = extract_tagged_sequences(aligned_read, &full_ref);
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
                   str::replace(aligned_name, " ", "."),
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


#[cfg(test)]
mod tests {
    use libwfa::{affine_wavefront::*, bindings::*, mm_allocator::*, penalties::*};

    #[test]
    fn simple_libwfa() {
        let alloc = MMAllocator::new(BUFFER_SIZE_8M as u64);

        let pattern = String::from("AAACGCTTCTGCACTTCGCGTGATATCATTACGTTTTGTCCCTCTACGATGTACCTGTCATCTTAGCTAAGACGACAGGACATGTCCAGGAAGTGCTCGAGTACTTCCTGGCCCATGTACTCTGCGTTGATACCACTGCTT");
        let text = String::from("NNNNNNNNNNNNNNNNNNNNNNNNNNNNTTACGTTNNNNNNNNNNNNGATGTACCTGTCATCTTAGCTAAGATGACAGGACATGTCCAGGAAGTACTCGAGTACTTCCTGGCCCATGTACTCTGCGTTGATACCACTGCTT");

        let mut penalties = AffinePenalties {
            match_: 0,
            mismatch: 4,
            gap_opening: 6,
            gap_extension: 2,
        };

        let pat_len = pattern.as_bytes().len();
        let text_len = text.as_bytes().len();


        for i in 0..1000000 {
            let mut wavefronts = AffineWavefronts::new_complete(
                pat_len,
                text_len,
                &mut penalties,
                &alloc,
            );

            wavefronts
                .align(pattern.as_bytes(), text.as_bytes())
                .unwrap();

            let score = wavefronts.edit_cigar_score(&mut penalties);

            if i == 0 {
                println!("score: {}", score);
                wavefronts.print_cigar(pattern.as_bytes(), text.as_bytes());

                // The cigar can also be extracted as a byte vector
                let cigar = wavefronts.cigar_bytes_raw();
                let cg_str = std::str::from_utf8(&cigar).unwrap();
                println!("cigar: {}", cg_str);

                // Or as a prettier byte vector

                let cigar = wavefronts.cigar_bytes();
                let cg_str = std::str::from_utf8(&cigar).unwrap();
                println!("cigar: {}", cg_str);
            }
        }
    }
}