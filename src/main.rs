extern crate backtrace;
extern crate bgzip;
extern crate bio;
extern crate fastq;
extern crate flate2;
extern crate indicatif;
extern crate itertools;
#[macro_use]
extern crate lazy_static;
extern crate log;
extern crate ndarray;
extern crate needletail;
extern crate noodles_fastq;
extern crate num_traits;
extern crate petgraph;
extern crate rand;
extern crate rayon;
extern crate rust_htslib;
extern crate rust_spoa;
extern crate seq_io;
extern crate suffix;
extern crate tempfile;

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, Error, Write};
use std::io;
use std::path::{Path, PathBuf};
use std::str;
use std::str::FromStr;
use std::sync::{Arc, Mutex};

use backtrace::Backtrace;
use bio::alignment::Alignment;
use log::{LevelFilter, SetLoggerError};
use noodles_fastq as Fastq;
use petgraph::algo::connected_components;
use rand::Rng;
use rayon::prelude::*;
use seq_io::fasta::{OwnedRecord, Reader, Record};
use tempfile::{NamedTempFile, TempDir};

use alignment::alignment_matrix::*;
use alignment::scoring_functions::*;
use clap::Parser;
use clap::StructOpt;

use crate::consensus::consensus_builders::threaded_write_consensus_reads;
use crate::extractor::extract_tagged_sequences;
use crate::linked_alignment::*;
use crate::read_strategies::sequence_file_containers::*;
use crate::read_strategies::sequence_layout::*;
use crate::reference::fasta_reference::reference_file_to_struct;
use crate::sorters::known_list::KnownList;
use crate::sorters::known_list::KnownListConsensus;
use crate::sorters::sorter::{Sorter, SortStructure};
use crate::umis::sequence_clustering::*;

mod linked_alignment;
pub mod extractor;
mod simple_umi_clustering;

mod umis {
    pub mod bronkerbosch;
    pub mod sequence_clustering;
}

mod alignment {
    pub mod alignment_matrix;
    pub mod scoring_functions;
}

mod consensus {
    pub mod consensus_builders;
}

pub mod fasta_comparisons;

mod read_strategies {
    pub mod sequence_layout;
    pub mod sequence_file_containers;
    pub mod ten_x;
}

mod utils {
    pub mod file_utils;
    pub mod base_utils;
}

mod sorters {
    pub mod known_list;
    pub mod sorter;
    pub mod sort_streams;
    pub mod balanced_split_output;
}

mod reference {
    pub mod fasta_reference;
}

#[derive(Debug, StructOpt)]
enum Cmd {
    Collapse,
    Align,
}

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    #[structopt(subcommand)]
    cmd: Cmd,

    #[clap(long)]
    reference: String,

    #[clap(long)]
    output_base: String,

    #[clap(long)]
    output: String,

    #[clap(long)]
    read1: String,

    #[structopt(long, default_value = "NONE")]
    read2: String,

    #[structopt(long, default_value = "NONE")]
    index1: String,

    #[structopt(long, default_value = "NONE")]
    index2: String,

    #[clap(long, default_value_t = 1)]
    threads: usize,

    #[clap(long)]
    outputupper: bool,

    #[clap(long)]
    read_template: String,

    #[structopt(long, default_value = "NONE")]
    known_list: String,

    #[structopt(long, default_value = "250")]
    max_bins: usize,

    #[structopt(long, default_value = "1")]
    sorting_threads: usize,

    #[structopt(long, default_value = "1")]
    processing_threads: usize,

    #[structopt(long, default_value = "NONE")]
    temp_dir: String,
}

fn main() {
    let parameters = Args::parse();
    println!("{:?}",&parameters.cmd);
    match &parameters.cmd {
        Cmd::Collapse => {
            let read_layout = LayoutType::from_str(&parameters.read_template).expect("Unable to parse read template type");

            let read_bundle = ReadFileContainer::new(&parameters.read1, &parameters.read2, &parameters.index1, &parameters.index2);

            let est_read_count = estimate_read_count(&parameters.read1).unwrap_or(0);

            println!("estimated read set/pair count = {}", est_read_count);

            let no_temp_file: String = String::from("NONE");
            let temp_dir: Option<PathBuf> = match &parameters.temp_dir {
                no_temp_file => None,
                _ => Some(PathBuf::from(parameters.temp_dir)),
            };

            let tmp_location = TempDir::new().unwrap();

            let mut run_specs = RunSpecifications {
                estimated_reads: est_read_count,
                sorting_file_count: parameters.max_bins,
                sorting_threads: parameters.sorting_threads,
                processing_threads: parameters.processing_threads,
                tmp_location,
                file_count: 0,
            };


            // setup our thread pool
            rayon::ThreadPoolBuilder::new().num_threads(parameters.threads).build_global().unwrap();

            let mut known_list = KnownList::new(&parameters.known_list, 7);
            let mut known_list_hash = HashMap::new();
            known_list_hash.insert(read_layout.clone(), Arc::new(Mutex::new(known_list)));

            println!("Sorting");
            let sort_structure = SortStructure::from_layout(&read_layout, known_list_hash);
            let read_piles = Sorter::sort(sort_structure, &read_bundle, &"./tmp/".to_string(), &"test_sorted.txt.gz".to_string(), &read_layout, &mut run_specs);

            println!("Consensus");
            threaded_write_consensus_reads(read_piles,
                                           &parameters.output_base,
                                           &ReadPattern::from_read_file_container(&read_bundle),
                                           &run_specs);
        }
        Cmd::Align => {
            let reference = reference_file_to_struct(&parameters.reference);

            let output_file = File::create(parameters.output).unwrap();

            let read_iterator = ReadIterator::new(PathBuf::from(&parameters.read1),
                                                  Some(PathBuf::from(&parameters.read2)),
                                                  Some(PathBuf::from(&parameters.index1)),
                                                  Some(PathBuf::from(&parameters.index2)));
            ;

            let reference_lookup = find_seeds(&reference.sequence, 20);
            let output = Arc::new(Mutex::new(output_file)); // Mutex::new(gz));

            // setup our thread pool
            rayon::ThreadPoolBuilder::new().num_threads(parameters.threads).build_global().unwrap();


            let my_score = InversionScoring {
                match_score: 9.0,
                mismatch_score: -21.0,
                gap_open: -16.0,
                gap_extend: -6.0,
                inversion_penalty: -20.0,
                min_inversion_length: 20,
            };


            let my_aff_score = AffineScoring {
                match_score: 10.0,
                mismatch_score: -11.0,
                special_character_score: 7.0,
                gap_open: -15.0,
                gap_extend: -5.0,
            };

            // This is a little ugly since we're wrapping in the para_bridge, which introduces some scope issues


            read_iterator.par_bridge().for_each(|xx| {
                let x = xx.read_one;
                let name = &String::from(x.id()).to_string();

                let is_forward = orient_by_longest_segment(&x.seq().to_vec(), &reference.sequence, &reference_lookup);

                if is_forward.0 {
                    let fwd_score_mp = find_greedy_non_overlapping_segments(&x.seq().to_vec(), &reference.sequence, &reference_lookup);
                    let results = align_string_with_anchors(&x.seq().to_vec(), &reference.sequence, &fwd_score_mp, &my_score, &my_aff_score);

                    let output = Arc::clone(&output);
                    let mut output = output.lock().unwrap();
                    write!(output, ">ref\n{}\n>{}\n{}\n{}\n{}\n",
                           str::from_utf8(&results.1).unwrap(),
                           str::replace(name, " ", "_"),
                           str::from_utf8(&results.0).unwrap(),
                           str::from_utf8(&x.seq().to_vec()).unwrap(),
                            results.2.iter().map(|tag| format!("{}",tag)).collect::<Vec<String>>().join(",")).
                        expect("Unable to write to output file");

                    output.flush().expect("Unable to flush output");
                } else {
                    let fwd_score_mp = find_greedy_non_overlapping_segments(&x.seq().to_vec(), &reference.sequence, &reference_lookup);

                    let results = align_string_with_anchors(&reverse_complement(&x.seq().to_vec()), &reference.sequence, &fwd_score_mp, &my_score, &my_aff_score);

                    let output = Arc::clone(&output);
                    let mut output = output.lock().unwrap();
                    write!(output, ">ref\n{}\n>{}\n{}\n{}\n{}\n",
                           str::from_utf8(&results.1).unwrap(),
                           str::replace(name, " ", "_"),
                           str::from_utf8(&results.0).unwrap(),
                           str::from_utf8(&reverse_complement(&x.seq().to_vec())).unwrap(),
                           results.2.iter().map(|tag| format!("{}",tag)).collect::<Vec<String>>().join(",")).
                        expect("Unable to write to output file");

                    output.flush().expect("Unable to flush output");
                }
            });
        }
    }
}

pub struct RunSpecifications {
    pub estimated_reads: usize,
    pub sorting_file_count: usize,
    pub sorting_threads: usize,
    pub processing_threads: usize,
    pub tmp_location: TempDir,
    pub file_count: u64,
}

impl RunSpecifications {
    pub fn create_temp_file(&mut self) -> PathBuf {
        let file_path = PathBuf::from("/analysis/tmp/").join(self.file_count.to_string());
        self.file_count += 1;
        file_path
    }
}

#[allow(dead_code)]
struct AlignedWithFeatures {
    alignment: Alignment,
    read_id: String,
    read: Vec<u8>,
    reference: Vec<u8>,
    features: BTreeMap<String, String>,
}

fn to_two_line_fasta(align_features: AlignedWithFeatures, output_upper: bool) -> String {
    if output_upper {
        format!(">@{}_{}\n{}\n>ref\n{}\n", align_features.read_id,
                align_features.features.iter().filter(|(s, _t)| **s != "r".to_string() && **s != "e".to_string()).map(|(s, t)| format!("{}{}", &**s, &**t)).collect::<Vec<_>>().join(","),
                align_features.features.get(&"r".to_string()).unwrap(),
                align_features.features.get(&"e".to_string()).unwrap())
    } else {
        format!(">@{}_{}\n{}\n>ref\n{}\n",
                align_features.read_id,
                align_features.features.iter().filter(|(s, _t)| **s != "r".to_string() && **s != "e".to_string()).map(|(s, t)| format!("{}{}", &**s, &**t)).collect::<Vec<_>>().join(","),
                &format!("{}", String::from_utf8_lossy(align_features.read.as_slice())),
                &format!("{}", String::from_utf8_lossy(align_features.reference.as_slice())))
    }
}


