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
extern crate derive_more;
extern crate shardio;
extern crate anyhow;
extern crate sift4;
extern crate phf;
extern crate rust_htslib;
extern crate noodles_bam;
#[feature(alignment)]
extern crate noodles_util;
extern crate bstr; 

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
use crate::calling::call_events::BamCallingParser;
use crate::collapse::collapse;
use crate::read_strategies::sequence_layout::SequenceLayout;
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
mod calling {
//    pub mod bam_file_to_cell_list;
    pub mod call_events;
}

mod consensus {
    pub mod consensus_builders;
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

mod reference {
    pub mod fasta_reference;
}


#[derive(Subcommand, Debug)]
enum Cmd {
    Collapse {
        #[clap(long)]
        outbam: String,

        #[clap(long)]
        read_structure: String,

        #[clap(long, default_value = "1")]
        threads: usize,

        #[clap(long, default_value = "NONE")]
        temp_dir: String,

        #[clap(long)]
        inbam: String,

        #[clap(long)]
        find_inversions: bool,

        #[clap(long)]
        fast_reference_lookup: bool,

        #[clap(long, default_value = "0")]
        max_deletion: usize,
    }/*,
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

    }*/,

    Call {
        #[clap(long)]
        read_structure: String,

        #[clap(long)]
        bam: String,

        #[clap(long)]
        output: String,
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
            outbam,
            read_structure,
            threads,
            temp_dir: _,
            inbam,
            find_inversions,
            fast_reference_lookup,
            max_deletion,

        } => {
            let my_yaml = SequenceLayout::from_yaml(read_structure);

            let mut tmp = InstanceLivedTempDir::new().unwrap();

            collapse(outbam,
                     &mut tmp,
                     &my_yaml,
                     inbam);
        },
/*
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
            let my_yaml = SequenceLayoutDesign::from_yaml(read_structure);
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
        },*/
        Cmd::Call {
            read_structure,
            bam,
            output,
        } => {
            let my_yaml = SequenceLayout::from_yaml(read_structure);

            let mut parser = BamCallingParser::new(&my_yaml);

            parser.output_bam_file_entries(bam.as_str(), output.as_str() ).expect("Unable to process events");

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
