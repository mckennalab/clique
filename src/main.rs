extern crate backtrace;
extern crate bgzip;
extern crate bio;
extern crate chrono;
extern crate fastq;
extern crate flate2;
extern crate indicatif;
extern crate itertools;
#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate log;
extern crate ndarray;
extern crate needletail;
extern crate noodles_fastq;
extern crate num_traits;
extern crate petgraph;
extern crate pretty_env_logger;
extern crate rand;
extern crate rayon;
extern crate rust_htslib;
extern crate rust_spoa;
extern crate seq_io;
extern crate serde;
extern crate suffix;
extern crate tempfile;
extern crate serde_yaml;
extern crate symspell;

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
use std::sync::{Arc, Mutex};

use backtrace::Backtrace;
use bio::alignment::Alignment;
use chrono::Local;
use log::{debug, error, info, Level, LevelFilter, log_enabled};
use noodles_fastq as Fastq;
use petgraph::algo::connected_components;
use rand::Rng;
use rayon::prelude::*;
use seq_io::fasta::{OwnedRecord, Reader, Record};
use tempfile::{Builder, NamedTempFile};
use tempfile::TempDir as ActualTempDir;

use alignment::alignment_matrix::*;
use alignment::scoring_functions::*;
use clap::Parser;
use clap::Subcommand;
use itertools::Itertools;
use nanoid::nanoid;

use crate::extractor::extract_tagged_sequences;
use crate::linked_alignment::*;
use crate::reference::fasta_reference::reference_file_to_struct;

use pretty_trace::*;
use crate::read_strategies::read_set::ReadIterator;

mod linked_alignment;
pub mod extractor;
pub mod sequence_lookup;

mod read_strategies {
    pub mod read_set;
    pub mod sequence_layout;
}
mod alignment {
    pub mod alignment_matrix;
    pub mod scoring_functions;
    pub mod fasta_bit_encoding;
}

pub mod fasta_comparisons;

mod utils {
    pub mod file_utils;
    pub mod base_utils;
    pub mod read_utils;
}


mod reference {
    pub mod fasta_reference;
}

#[derive(Subcommand, Debug)]
enum Cmd {
    Collapse{
        #[clap(long)]
        /// Name of the package to search
        package_name: String,

        #[clap(long)]
        output_base: String,

        #[clap(long)]
        output: String,

        #[clap(long)]
        read_template: String,

        #[clap(long, default_value = "NONE")]
        known_list: String,

        #[clap(long, default_value = "250")]
        max_bins: usize,

        #[clap(long, default_value = "1")]
        sorting_threads: usize,

        #[clap(long, default_value = "NONE")]
        temp_dir: String,

        #[clap(long)]
        read1: String,

        #[clap(long, default_value = "NONE")]
        read2: String,

        #[clap(long, default_value = "NONE")]
        index1: String,

        #[clap(long, default_value = "NONE")]
        index2: String,

        #[clap(long, default_value_t = 1)]
        threads: usize,

        #[clap(long, default_value = "1")]
        processing_threads: usize,
    },
    Align{
        #[clap(long)]
        use_capture_sequences: bool,

        #[clap(long)]
        only_output_captured_ref: bool,

        #[clap(long)]
        reference: String,

        #[clap(long)]
        output: String,

        #[clap(long, default_value = "2")]
        max_reference_multiplier: usize,

        #[clap(long)]
        read1: String,

        #[clap(long, default_value = "NONE")]
        read2: String,

        #[clap(long, default_value = "NONE")]
        index1: String,

        #[clap(long, default_value = "NONE")]
        index2: String,

        #[clap(long, default_value_t = 1)]
        threads: usize,

        #[clap(long, default_value = "1")]
        processing_threads: usize,
    },
}

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    #[clap(subcommand)]
    cmd: Cmd,



}


fn main() {
    PrettyTrace::new().ctrlc().on();

    if let Err(_) = std::env::var("RUST_LOG") {
        std::env::set_var("RUST_LOG", "warn");
    }

    pretty_env_logger::init();

    let parameters = Args::parse();
    trace!("{:?}", &parameters.cmd);

    match &parameters.cmd {
        Cmd::Collapse{ .. } => {
            //merger(&parameters);
        }

        Cmd::Align{ use_capture_sequences,
            only_output_captured_ref,
            reference,
            output, max_reference_multiplier,
            read1,
            read2,
            index1,
            index2,
            threads,
            processing_threads
        } => {
            align_reads(use_capture_sequences,
                        only_output_captured_ref,
                        reference,
                        output, max_reference_multiplier,
                        read1,
                        read2,
                        index1,
                        index2,
                        threads,
                        processing_threads);
        }
    }
}

fn align_reads(use_capture_sequences : &bool,
               only_output_captured_ref: &bool,
               reference: &String,
               output :&String,
               max_reference_multiplier: &usize,
               read1: &String,
               read2: &String,
               index1: &String,
               index2: &String,
               threads: &usize,
               processing_threads: &usize) {

    let reference = reference_file_to_struct(reference);

    let output_file = File::create(&output).unwrap();

    let read_iterator = ReadIterator::new(PathBuf::from(&read1),
                                          Some(PathBuf::from(&read2)),
                                          Some(PathBuf::from(&index1)),
                                          Some(PathBuf::from(&index2)));
    ;

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

    //let mut too_long : Arc<Mutex<usize>> = Arc::new(Mutex::new(0));

    read_iterator.par_bridge().for_each(|xx| {
        let x = xx.read_one;
        let name = &String::from(x.id()).to_string();

        let read_length = x.seq().to_vec().len();
        if read_length > reference.sequence.len() * max_reference_multiplier {
            info!("Dropping read of length {}",read_length);
            //let mut cnt = *too_long.clone().lock().unwrap();
            //cnt += 1;
        } else {
            let orientation = orient_by_longest_segment(&x.seq().to_vec(), &reference.sequence, &reference_lookup).0;
            let forward_oriented_seq = if orientation {
                x.seq().to_vec()
            } else {
                reverse_complement(&x.seq().to_vec())
            };

            let fwd_score_mp = find_greedy_non_overlapping_segments(&forward_oriented_seq, &reference.sequence, &reference_lookup);
            let first_hit = fwd_score_mp.alignment_segments.get(0).unwrap_or(&MatchedPosition {
                search_start: 0,
                ref_start: 0,
                length: 0,
            });
            let results = align_string_with_anchors(&forward_oriented_seq, &reference.sequence, &fwd_score_mp, &my_score, &my_aff_score);


            let extracted_seqs = if *use_capture_sequences {
                Some(extract_tagged_sequences(&results.aligned_read, &results.aligned_ref))
            } else {
                None
            };

            let cigar_string = results.cigar_tags.iter().map(|tag| format!("{}", tag)).collect::<Vec<String>>().join(",");

            let output = Arc::clone(&output);
            let mut output = output.lock().unwrap();

            if *only_output_captured_ref {
                match extracted_seqs {
                    None => {
                        warn!("unable to extract sequences from read {}",x.id())
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

                        write!(output, ">ref\n{}\n>{}\n{}\n",
                               ref_seq,
                               str::replace(name, " ", "_"),
                               read_seq,
                        ).expect("Unable to write to output file");
                    }
                };
            } else {
                let collapsed_tags = match extracted_seqs {
                    None => { String::from("NONE") }
                    Some(x) => { x.iter().map(|k| format!("key={}:{}", &k.0, &k.1)).join(";") }
                };
                write!(output, ">ref\n{}\n>{};{};{}\n{}\n",
                       str::from_utf8(&results.aligned_ref).unwrap(),
                       str::replace(name, " ", "_"),
                       collapsed_tags,
                       cigar_string,
                       str::from_utf8(&results.aligned_read).unwrap(),
                ).expect("Unable to write to output file");
            }
        }
    });

}
/*
fn merger(parameters: &Args) {
    let read_layout = LayoutType::from_str(&parameters.read_template).expect("Unable to parse read template type");

    let read_bundle = ReadFileContainer::new(&parameters.read1, &parameters.read2, &parameters.index1, &parameters.index2);

    let est_read_count = estimate_read_count(&parameters.read1).unwrap_or(0);

    info!("estimated read set/pair count = {}", est_read_count);

    let no_temp_file: String = String::from("NONE");
    let temp_dir: Option<PathBuf> = match &parameters.temp_dir {
        no_temp_file => None,
        _ => Some(PathBuf::from(&parameters.temp_dir)),
    };

    let tmp_location = TempDir::new().unwrap();

    let mut run_specs = RunSpecifications {
        estimated_reads: est_read_count,
        sorting_file_count: parameters.max_bins,
        sorting_threads: parameters.sorting_threads,
        processing_threads: parameters.processing_threads,
        tmp_location: tmp_location.into(),
    };

    // setup our thread pool
    rayon::ThreadPoolBuilder::new().num_threads(parameters.threads).build_global().unwrap();

    let mut known_list = KnownList::new(&parameters.known_list);

    info!("------------------------- Sorting ------------------------- ");
    let sort_structure = SortStructure::from_umi_type(&read_layout, &known_list);

    let read_pattern = ReadPattern::from_read_file_container(&read_bundle);
    info!("Running with read pattern {},{},{},{}",read_pattern,read_bundle.read_two.is_some(),read_bundle.index_one.is_some(),read_bundle.index_two.is_some());

    let read_piles = Sorter::sort(
        UMIType::TENXRT,
        sort_structure,
        &read_bundle,
        &"./tmp/".to_string(),
        &"test_sorted.txt.gz".to_string(),
        &read_layout,
        &read_pattern,
        &mut run_specs);

    info!("------------------------- Building Consensus ------------------------- ");
    threaded_write_consensus_reads(read_piles,
                                   &parameters.output_base,
                                   &ReadPattern::from_read_file_container(&read_bundle),
                                   &run_specs);
}
*/

pub struct RunSpecifications {
    pub estimated_reads: usize,
    pub sorting_file_count: usize,
    pub sorting_threads: usize,
    pub processing_threads: usize,
    pub tmp_location: Arc<TempDir>,
}


#[derive(Debug)]
pub struct TempDir(Option<ActualTempDir>);

// Forward inherent methods to the tempdir crate.
impl TempDir {
    pub fn new() -> Result<TempDir>
    { ActualTempDir::new().map(Some).map(TempDir) }

    pub fn path(&self) -> &Path
    { self.0.as_ref().unwrap().path() }
}

/// Leaks the inner TempDir if we are unwinding.
impl Drop for TempDir {
    fn drop(&mut self) {
        if ::std::thread::panicking() {
            ::std::mem::forget(self.0.take())
        }
    }
}

impl RunSpecifications {
    pub fn create_temp_file(&self) -> PathBuf {
        let file_path = PathBuf::from(&self.tmp_location.clone().path()).join(nanoid!());
        file_path
    }
}

impl Clone for RunSpecifications {
    fn clone(&self) -> RunSpecifications {
        RunSpecifications {
            estimated_reads: self.estimated_reads,
            sorting_file_count: self.sorting_file_count,
            sorting_threads: self.sorting_threads,
            processing_threads: self.processing_threads,
            tmp_location: Arc::clone(&self.tmp_location),
        }
    }
}
