#![feature(extern_types)]


#![feature(ascii_char)]
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
extern crate num_traits;
extern crate petgraph;
extern crate pretty_env_logger;
extern crate rand;
extern crate rayon;
extern crate seq_io;
extern crate serde;
extern crate suffix;
extern crate tempfile;
extern crate serde_yaml;
extern crate symspell;
extern crate shardio;
extern crate anyhow;
extern crate phf;
extern crate rust_htslib;
extern crate noodles_bam;
extern crate noodles_util;
extern crate bstr;
extern crate rust_star;
extern crate libc;
extern crate clap;
extern crate indexmap;
extern crate sigalign;
extern crate counter;
extern crate rustc_hash;
extern crate noodles_sam;
extern crate nohash_hasher;
extern crate vpsearch;
extern crate nanoid;
extern crate spoa;

use ::std::io::Result;
use std::path::{Path, PathBuf};
use std::str;
use std::sync::{Arc};

use tempfile::{TempDir as ActualTempDir};

use clap::Parser;
use clap::Subcommand;
use nanoid::nanoid;
use consensus::consensus_builders::{MergeStrategy, ReadOutputApproach};
use crate::alignment_functions::align_reads;
use crate::collapse::collapse;
use crate::read_strategies::sequence_layout::SequenceLayout;
use crate::reference::fasta_reference::ReferenceManager;

mod linked_alignment;
pub mod extractor;
pub mod sequence_lookup;

#[allow(dead_code)]
const FASTA_UNSET: u8 = b'-';
#[allow(dead_code)]
const FASTA_N: u8 = b'N';
#[allow(dead_code)]
const FASTA_A: u8 = b'A';
#[allow(dead_code)]
const FASTA_G: u8 = b'G';
#[allow(dead_code)]
const FASTA_C: u8 = b'C';
#[allow(dead_code)]
const FASTA_T: u8 = b'T';


mod read_strategies {
    pub mod read_set;
    pub mod sequence_layout;
    pub mod read_disk_sorter;
}

mod alignment {
    pub mod alignment_matrix;
    pub mod scoring_functions;
    //pub mod fasta_bit_encoding;
}

mod umis {
    pub mod sequence_clustering;
    pub mod bronkerbosch;
    pub mod correct_tags;
    
    pub mod known_list;
}

mod consensus {
    pub mod consensus_builders;

    pub mod stretcher;
}

pub mod fasta_comparisons;

mod utils {
    pub mod base_utils;
    pub mod read_utils;
}

// mod alignment_functions;
mod sorter;
pub mod merger;
mod collapse;
mod alignment_manager;
mod alignment_functions;

mod reference {
    pub mod fasta_reference;
}

#[derive(Subcommand, Debug)]
enum Cmd {
    Collapse {
        #[clap(long)]
        output_bam_file: String,

        #[clap(long)]
        read_structure: String,

        #[clap(long, default_value = "1")]
        threads: usize,

        #[clap(long, default_value = "NONE")]
        temp_dir: String,

        #[clap(long)]
        input_bam_file: String,

        #[clap(long)]
        find_inversions: bool,

        #[clap(long)]
        fast_reference_lookup: bool,

        #[clap(long, default_value = "0")]
        max_deletion: usize,

        #[clap(long, action=clap::ArgAction::SetTrue)]
        correct_only: bool,

    },
    Align {
        #[clap(long)]
        read_structure: String,

        #[clap(long)]
        output_bam_file: String,

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

    if let Err(_) = std::env::var("RUST_LOG") {
        std::env::set_var("RUST_LOG", "info");
    }

    pretty_env_logger::init_timed();

    let parameters = Args::parse();
    trace!("{:?}", &parameters.cmd);

    match &parameters.cmd {
        Cmd::Collapse {
            output_bam_file: outbam,
            read_structure,
            threads: _,
            temp_dir: _,
            input_bam_file: inbam,
            find_inversions: _,
            fast_reference_lookup: _,
            max_deletion: _,
            correct_only: correction_only,

        } => {
            let my_yaml = SequenceLayout::from_yaml(read_structure);

            let mut tmp = InstanceLivedTempDir::new().unwrap();

            let correction = match *correction_only {
                true => {
                    info!("Reads will be corrected and not collapsed");
                    ReadOutputApproach::Correct
                },
                false => {
                    info!("Reads will be collapsed by the combination of tags");
                    ReadOutputApproach::Collapse
                }
            };

            collapse(outbam,
                     &mut tmp,
                     &my_yaml,
                     inbam,
                     &MergeStrategy::Stretcher, // TODO parameterize,
                     &correction,
            );
        },

        Cmd::Align {
            read_structure,
            output_bam_file: output,
            max_reference_multiplier,
            min_read_length,
            read1,
            read2,
            index1,
            index2,
            threads,
            find_inversions,
        } => {
            let my_yaml = SequenceLayout::from_yaml(read_structure);
            let rm = ReferenceManager::from_yaml_input(&my_yaml, 8, 4);

            let output_path = Path::new(&output);

            align_reads(&my_yaml,
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

pub struct RunSpecifications {
    pub estimated_reads: usize,
    pub sorting_file_count: usize,
    pub sorting_threads: usize,
    pub processing_threads: usize,
    pub tmp_location: Arc<InstanceLivedTempDir>,
}

#[derive(Debug)]
pub struct InstanceLivedTempDir(Option<ActualTempDir>);

// Forward inherent methods to the tempdir crate.
impl InstanceLivedTempDir {
    pub fn new() -> Result<InstanceLivedTempDir>
    { ActualTempDir::new().map(Some).map(InstanceLivedTempDir) }

    pub fn temp_file(&mut self, name: &str) -> PathBuf
    {
        self.0.as_ref().unwrap().path().join(name).clone()
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
