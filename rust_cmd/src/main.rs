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
extern crate seq_io;
extern crate serde;
extern crate suffix;
extern crate tempfile;
extern crate serde_yaml;
extern crate symspell;
extern crate derive_more;
extern crate shardio;
extern crate anyhow;
extern crate sift4;


use ::std::io::Result;
use std::path::{Path, PathBuf};
use std::str;
use std::sync::{Arc};

use tempfile::{TempDir as ActualTempDir};

use alignment::alignment_matrix::*;
use alignment::scoring_functions::*;
use clap::Parser;
use clap::Subcommand;
use nanoid::nanoid;

use pretty_trace::*;
use crate::alignment_functions::align_reads;
use crate::collapse::collapse;
use crate::read_strategies::sequence_layout::SequenceLayoutDesign;
use crate::reference::fasta_reference::ReferenceManager;

mod linked_alignment;
pub mod extractor;
pub mod sequence_lookup;

mod read_strategies {
    pub mod read_set;
    pub mod sequence_layout;
    pub mod read_disk_sorter;
}

mod alignment {
    pub mod alignment_matrix;
    pub mod scoring_functions;
    pub mod fasta_bit_encoding;
}

mod umis {
    pub mod sequence_clustering;
    pub mod bronkerbosch;
    pub mod known_list;
}
mod consensus {
    pub mod consensus_builders;
}
pub mod fasta_comparisons;

mod utils {
    pub mod base_utils;
    pub mod read_utils;
}

mod alignment_functions;
mod sorter;
pub mod merger;
mod collapse;

mod reference {
    pub mod fasta_reference;
}


#[derive(Subcommand, Debug)]
enum Cmd {
    Collapse {
        #[clap(long)]
        reference: String,

        #[clap(long)]
        outbam: String,

        #[clap(long)]
        read_structure: String,

        #[clap(long, default_value = "1")]
        threads: usize,

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

    },
    Align {
        #[clap(long)]
        use_capture_sequences: bool,

        #[clap(long)]
        read_structure: String,

        #[clap(long)]
        only_output_captured_ref: bool,

        #[clap(long)]
        to_fake_fastq: bool,

        #[clap(long)]
        reference: String,

        #[clap(long)]
        output: String,

        #[clap(long, default_value = "2")]
        max_reference_multiplier: usize,

        #[clap(long, default_value = "50")]
        min_read_length: usize,

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

        #[clap(long)]
        find_inversions: bool,

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
        std::env::set_var("RUST_LOG", "info");
    }

    pretty_env_logger::init_timed();

    let parameters = Args::parse();
    trace!("{:?}", &parameters.cmd);

    match &parameters.cmd {
        Cmd::Collapse {
            reference,
            outbam,
            read_structure,
            threads,
            temp_dir: _,
            read1,
            read2,
            index1,
            index2,
        } => {
            let my_yaml = SequenceLayoutDesign::from_yaml(read_structure).unwrap();

            let mut tmp = InstanceLivedTempDir::new().unwrap();
            collapse(reference,
                     outbam,
                     &mut tmp,
                     &my_yaml,
                     &1.2,
                     &50,
                     read1,
                     read2,
                     index1,
                     index2,
                     threads);
        }

        Cmd::Align {
            use_capture_sequences,
            read_structure,
            only_output_captured_ref,
            to_fake_fastq,
            reference,
            output,
            max_reference_multiplier,
            min_read_length,
            read1,
            read2,
            index1,
            index2,
            threads,
            find_inversions,
        } => {

            let my_yaml = SequenceLayoutDesign::from_yaml(read_structure).unwrap();

            // load up the reference files
            let rm = ReferenceManager::from(&reference, 8);

            let output_path = Path::new(&output);

            align_reads(use_capture_sequences,
                        &my_yaml,
                        only_output_captured_ref,
                        to_fake_fastq,
                        &rm,
                        &output_path,
                        max_reference_multiplier,
                        min_read_length,
                        read1,
                        read2,
                        index1,
                        index2,
                        threads,
                        find_inversions);
        }
    }
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
    pub tmp_location: Arc<InstanceLivedTempDir>,
}
/*
#[derive(Debug)]
pub struct InstanceLivedTempDir{
    temp_dir: Option<ActualTempDir>,
    next_file_id: usize,
}



//Forward inherent methods to the tempdir crate.
impl InstanceLivedTempDir {
    pub fn new() -> Result<InstanceLivedTempDir>
    { ActualTempDir::new().map(Some).map(InstanceLivedTempDir) }

    pub fn temp_file(&mut self, name: &str) -> PathBuf
    {
        self.next_file_id += 1; // indexing will then start on one
        self.temp_dir.as_ref().unwrap().path().join(self.next_file_id.to_string()).join(name.clone()).clone()
    }

    pub fn path(&self) -> &Path
    { self.temp_dir.as_ref().unwrap().path() }

}

/// Leaks the inner TempDir if we are unwinding.
impl Drop for InstanceLivedTempDir {
    fn drop(&mut self) {
        if ::std::thread::panicking() {
            ::std::mem::forget(self.temp_dir.take())
        }
    }
}*/

#[derive(Debug)]
pub struct InstanceLivedTempDir(Option<ActualTempDir>);

// Forward inherent methods to the tempdir crate.
impl InstanceLivedTempDir {
    pub fn new() -> Result<InstanceLivedTempDir>
    { ActualTempDir::new().map(Some).map(InstanceLivedTempDir) }

    pub fn temp_file(&mut self, name: &str) -> PathBuf
    {
        self.0.as_ref().unwrap().path().join(name.clone()).clone()
    }

    pub fn path(&self) -> &Path
    { self.0.as_ref().unwrap().path() }
}

/// Leaks the inner TempDir if we are unwinding.
impl Drop for InstanceLivedTempDir {
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
