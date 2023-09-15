extern crate rust_spoa;

use std::borrow::Borrow;
use std::cmp::Ordering;
use std::collections::{HashMap, VecDeque};
use rayon::prelude::*;
use rust_spoa::poa_consensus;
use shardio::{Range, ShardReader, ShardWriter};
use crate::alignment::fasta_bit_encoding::{FASTA_UNSET, FastaBase};
use crate::read_strategies::read_disk_sorter::{SortedAlignment, SortingReadSetContainer};
use std::io::Write;
use counter::Counter;
use indicatif::ProgressBar;
use rust_htslib::bam::record::{Cigar, CigarString};
use crate::alignment::alignment_matrix::AlignmentTag;
use crate::alignment_functions::{create_sam_record, perform_wavefront_alignment, setup_sam_writer, simplify_cigar_string};
use crate::reference::fasta_reference::ReferenceManager;
use rand::prelude::*;
use rand::seq::SliceRandom;


pub struct ConsensusCandidate {
    pub reads: Vec<SortingReadSetContainer>,
    pub global_umi: Vec<u8>,
}

pub struct ConsensusResult {
    pub read_one: Vec<u8>,
    pub read_two: Option<Vec<u8>>,
    pub global_umi: Vec<u8>,
}

pub fn null_cap(strs: &[Vec<u8>]) -> Vec<Vec<u8>> {
    strs.iter().map(|r| {
        let mut ret = r.clone();
        ret.push(b'\0');
        ret
    }).collect::<Vec<Vec<u8>>>()
}




pub fn write_consensus_reads(reader: &ShardReader<SortingReadSetContainer>,
                             output_file: &String,
                             threads: &usize,
                             levels: usize,
                             read_counts: &usize,
                             reference_manager: &ReferenceManager,
                             maximum_reads_before_downsampling: &usize) {

    info!("Writing consensus reads to {}", output_file);
    let mut writer = setup_sam_writer(output_file, reference_manager).unwrap();

    let mut last_read: Option<SortingReadSetContainer> = None;
    let mut buffered_reads = VecDeque::new();
    let bar = ProgressBar::new(read_counts.clone() as u64);

    reader.iter_range(&Range::all()).unwrap().enumerate().for_each(|(i, x)| {
        bar.inc(1);
        let x = x.unwrap();
        assert_eq!(x.ordered_sorting_keys.len(), levels);
        if last_read.is_some() && &x.cmp(last_read.as_ref().unwrap()) == &Ordering::Equal {
            buffered_reads.push_back(x.clone());
            last_read = Some(x);
        } else {
            if buffered_reads.len() > 0 {
                let mut added_tags = HashMap::new();
                added_tags.insert((b'r', b'c'), buffered_reads.len().to_string());

                // downsample if needed
                if buffered_reads.len() > *maximum_reads_before_downsampling {
                    let mut rng = rand::thread_rng();
                    let target_proportion = *maximum_reads_before_downsampling as f64 / buffered_reads.len() as f64;
                    println!("Downsampling from {} to {} reads, target prop {} random 1 {} random 2 {}", buffered_reads.len(), *maximum_reads_before_downsampling,target_proportion, rng.gen::<f64>(), rng.gen::<f64>());
                    //buffered_reads = buffered_reads.iter().enumerate().filter(|(_, _)| rng.gen::<f64>() < target_proportion).map(|(_, x)| x.clone()).collect::<VecDeque<SortingReadSetContainer>>();
                    buffered_reads = VecDeque::from(buffered_reads.clone().into_iter().choose_multiple(&mut rng, *maximum_reads_before_downsampling));
                }
                added_tags.insert((b'd', b'c'), buffered_reads.len().to_string());

                let mut consensus_reads = create_poa_consensus(&buffered_reads);
                let consensus_reference = Counter::<Vec<u8>, usize>::init(buffered_reads.iter().map(|x| x.aligned_read.ref_name.clone().as_bytes().to_vec()).collect::<Vec<Vec<u8>>>()).most_common_ordered();

                let top_ref = &consensus_reference.get(0).unwrap().clone().0;
                let reference_pointer = reference_manager.references.get(reference_manager.reference_name_to_ref.get(top_ref).unwrap()).unwrap();
                let new_alignment = perform_wavefront_alignment(&reference_pointer.sequence, &FastaBase::from_vec_u8(&consensus_reads));
                let read_names = buffered_reads.iter().map(|x| x.aligned_read.read_name.clone()).collect::<Vec<String>>();

                added_tags.insert((b'a', b'r'), read_names.join(","));

                let sam_read = create_sam_record(read_names.get(0).clone().unwrap(),
                                                 &new_alignment.read_aligned,
                                                 &new_alignment.reference_aligned,
                                                 &reference_pointer.sequence_u8,
                                                 &reference_read_to_cigar_string(&new_alignment.reference_aligned, &new_alignment.read_aligned),
                                                 &true,
                                                 added_tags);

                writer.write(&sam_read).unwrap();

                buffered_reads.clear();
            }
            buffered_reads.push_back(x.clone());
            last_read = Some(x)
        }
    });
}

pub fn reference_read_to_cigar_string(reference_seq: &Vec<FastaBase>, read_seq: &Vec<FastaBase>) -> CigarString {
    let mut result: Vec<AlignmentTag> = Vec::new();

    for index in 0..reference_seq.len() {
        if *reference_seq.get(index).unwrap() == FASTA_UNSET {
            result.push(AlignmentTag::Ins(1))
        } else if *read_seq.get(index).unwrap() == FASTA_UNSET {
            result.push(AlignmentTag::Del(1))
        } else {
            result.push(AlignmentTag::MatchMismatch(1))
        }
    }

    let simplied_cigar = simplify_cigar_string(&result);
    CigarString::try_from(
        simplied_cigar.iter().map(|m| format!("{}", m)).collect::<Vec<String>>().join("").as_bytes()).
        expect("Unable to parse cigar string.")
}


pub fn create_poa_consensus(sequences: &VecDeque<SortingReadSetContainer>) -> Vec<u8> {
    let max_length = sequences.iter().map(|n| n.aligned_read.aligned_read.len()).collect::<Vec<usize>>();
    let max_length = max_length.iter().max().unwrap();
    let first_name = match sequences.iter().next() {
        Some(x) => x.aligned_read.read_name.clone(),
        None => "UNKNOWN".to_string(),
    };

    let base_sequences = sequences.iter().map(|n| {
        let mut y = FastaBase::to_vec_u8(&n.aligned_read.aligned_read);
        y.push(b'\0');
        y
    }).collect::<Vec<Vec<u8>>>();
    poa_consensus(&base_sequences, max_length.clone(), 1, 5, -4, -3, -1).
        iter().filter(|x| *x != &b'-').map(|x| *x).collect::<Vec<u8>>()
}