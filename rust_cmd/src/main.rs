#![feature(extern_types)]


#![feature(ascii_char)]
//! # clique
//!
//! A high-performance CLI for amplicon / lineage-tracing sequencing data from
//! both Illumina and long-read (Nanopore / PacBio) platforms. It aligns reads
//! to one or more reference amplicons, collapses the reads of each molecule by
//! their unique molecular identifiers (UMIs) and static identifiers into a
//! consensus, and calls the CRISPR edits ("events") each molecule carries.
//!
//! ## Pipeline
//!
//! 1. **`align`** — FASTQ(s) + a read-structure YAML → a BAM. Each read is
//!    oriented, matched to its best reference, aligned, and annotated; UMI/tag
//!    subsequences and edit events are written as BAM aux tags.
//! 2. **`collapse`** — an aligned BAM + the YAML → a consensus BAM. Reads are
//!    grouped down a hierarchy of UMIs (each level corrected to a known list or
//!    clustered), then merged into one consensus (or corrected) read per
//!    molecule.
//! 3. **`genbank-to-yaml`** — an annotated GenBank file → a read-structure YAML
//!    (see [`genbank`]).
//!
//! ## Read structure
//!
//! A [`read_strategies::sequence_layout::SequenceLayout`] (loaded from YAML)
//! describes each reference amplicon: the UMI/barcode slots (marked by symbol
//! characters embedded in the reference sequence) and the CRISPR target sites
//! with their editing chemistry.
//!
//! ## BAM tag contract
//!
//! Downstream tools (e.g. the `cliqueR` R package) read these aux tags:
//! `e<symbol>` corrected extracted tag, `o<symbol>` original, `rc` read count,
//! `ar` read name(s), `as`/`rs` alignment score, `rm` alignment rate, and `ce`
//! the called-edit string (see [`events`]).
//!
//! ## Module map
//!
//! - [`read_strategies`] — FASTQ reading, the read-structure model, the on-disk
//!   sort container.
//! - [`reference`](mod@reference) — the multi-reference manager and its k-mer index.
//! - [`alignment`] / `alignment_functions` / `alignment_manager` — the aligner,
//!   scoring, reference selection, and BAM output.
//! - [`merger`] — merging paired / overlapping reads.
//! - [`extractor`] — pulling UMI/tag subsequences out of an alignment.
//! - [`umis`] — tag correction, clustering, known-list lookup, clique finding.
//! - [`consensus`] — per-molecule consensus building and writing.
//! - [`events`] — CRISPR edit calling; [`genbank`] — GenBank → YAML.

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
extern crate counter;
extern crate rustc_hash;
extern crate noodles_sam;
extern crate nohash_hasher;
extern crate vpsearch;
extern crate nanoid;
extern crate gb_io;

use ::std::io::Result;
use std::path::{Path, PathBuf};
use std::str;
use std::sync::{Arc};

use tempfile::{TempDir as ActualTempDir};

use clap::Parser;
use clap::Subcommand;
use clap::ValueEnum;
use nanoid::nanoid;
use consensus::consensus_builders::{MergeStrategy, ReadOutputApproach};
use crate::alignment_functions::align_reads;
use crate::collapse::{collapse, AlignmentFilterConfig};
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
pub mod events;
pub mod genbank;
mod collapse;
mod alignment_manager;
mod alignment_functions;

mod reference {
    pub mod fasta_reference;
    pub mod discriminating;
}

/// Aligner selection. Currently informational: the affine-gap aligner is used
/// regardless of this choice.
#[derive(Debug, Default, Clone, ValueEnum)]
enum Aligner {
    #[default]
    WFA,
    Degenerate,
    Inversion,
}

#[derive(Subcommand, Debug)]
enum Cmd {
    /// Collapse an aligned BAM into one consensus (or corrected) read per
    /// molecule, grouping reads by the UMI hierarchy in the read structure.
    Collapse {
        /// Output BAM path for the collapsed reads.
        #[clap(long)]
        output_bam_file: String,

        /// Read-structure YAML describing references, UMIs, and targets.
        #[clap(long)]
        read_structure: String,

        /// Number of worker threads.
        #[clap(long, default_value = "1")]
        threads: usize,

        /// Directory for temporary sort files ("NONE" uses the system temp dir).
        #[clap(long, default_value = "NONE")]
        temp_dir: String,

        /// Input aligned BAM. A `<BAM>.bai` index is used when available.
        #[clap(long)]
        input_bam_file: String,

        /// Detect inversions while collapsing.
        #[clap(long)]
        find_inversions: bool,

        /// Use the fast k-mer reference lookup instead of exhaustive search.
        #[clap(long)]
        fast_reference_lookup: bool,

        /// Maximum deletion length to tolerate.
        #[clap(long, default_value = "0")]
        max_deletion: usize,

        /// Minimum number of aligned non-UMI bases (capped to the reference's available bases).
        #[clap(long, default_value = "45")]
        min_aligned_bases: usize,

        /// Minimum identity among aligned non-UMI bases, from 0.0 to 1.0.
        #[clap(long, default_value = "0.8")]
        min_aligned_identity: f64,

        /// Only correct UMI/tag sequences; do not build consensus reads.
        #[clap(long, action=clap::ArgAction::SetTrue)]
        correct_only: bool,

    },
    /// Align FASTQ reads to their best-matching reference and write an annotated
    /// BAM (extracted tags plus called edits).
    Align {
        /// Read-structure YAML describing references, UMIs, and targets.
        #[clap(long)]
        read_structure: String,

        /// Output BAM path for the aligned reads.
        #[clap(long)]
        output_bam_file: String,

        /// Drop reads longer than this multiple of the longest reference.
        #[clap(long, default_value = "2")]
        max_reference_multiplier: usize,

        /// Skip reads shorter than this many bases.
        #[clap(long, default_value = "50")]
        min_read_length: usize,

        /// Read 1 FASTQ (required).
        #[clap(long)]
        read1: String,

        /// Read 2 FASTQ ("NONE" for single-end / long reads).
        #[clap(long, default_value = "NONE")]
        read2: String,

        /// Index 1 FASTQ ("NONE" if absent).
        #[clap(long, default_value = "NONE")]
        index1: String,

        /// Index 2 FASTQ ("NONE" if absent).
        #[clap(long, default_value = "NONE")]
        index2: String,

        /// Number of worker threads.
        #[clap(long, default_value_t = 1)]
        threads: usize,

        /// Aligner selection (currently informational).
        #[clap(long, arg_enum, default_value_t = Aligner::WFA)]
        aligner: Aligner,

        /// Select the reference for near-identical panels by comparing only the
        /// columns where the panel references differ (discriminating positions),
        /// emitting a top-2 margin confidence (dm/di/da BAM tags).
        #[clap(long, action=clap::ArgAction::SetTrue)]
        discriminating_classifier: bool,

        /// Minimum top-2 margin at the discriminating positions for a confident
        /// (non-ambiguous) reference call.
        #[clap(long, default_value = "1")]
        discriminating_min_margin: usize,

    },
    /// Generate a read-structure YAML from an annotated GenBank file.
    GenbankToYaml {
        /// Path to the annotated GenBank file.
        #[clap(long)]
        genbank: String,

        /// Path to write the generated read-structure YAML.
        #[clap(long)]
        output: String,

        /// Text that must appear in a feature's name for it to be included.
        #[clap(long, default_value = "lineage_target")]
        tag: String,

        /// Reference name for the emitted layout (defaults to the GenBank LOCUS).
        #[clap(long, default_value = "NONE")]
        reference_name: String,
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
            min_aligned_bases,
            min_aligned_identity,
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
            let alignment_filter = AlignmentFilterConfig::new(
                *min_aligned_bases,
                *min_aligned_identity,
            );

            collapse(outbam,
                     &mut tmp,
                     &my_yaml,
                     inbam,
                     &MergeStrategy::Stretcher, // TODO parameterize,
                     &correction,
                     &alignment_filter,
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
            aligner,
            discriminating_classifier,
            discriminating_min_margin,
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
                        aligner,
                        *discriminating_classifier,
                        *discriminating_min_margin);
        }

        Cmd::GenbankToYaml {
            genbank,
            output,
            tag,
            reference_name,
        } => {
            let records = gb_io::reader::parse_file(genbank)
                .unwrap_or_else(|e| panic!("Unable to parse GenBank file {}: {:?}", genbank, e));
            if records.is_empty() {
                panic!("No records found in GenBank file {}", genbank);
            }
            if records.len() > 1 {
                warn!(
                    "GenBank file {} has {} records; using the first ({:?})",
                    genbank, records.len(), records[0].name
                );
            }
            let ref_name = if reference_name == "NONE" {
                None
            } else {
                Some(reference_name.clone())
            };
            let opts = genbank::GenbankToYamlOptions { tag: tag.clone(), reference_name: ref_name };
            let layout = genbank::genbank_to_layout(&records[0], &opts)
                .unwrap_or_else(|e| panic!("Unable to build read structure from GenBank: {}", e));

            let yaml = serde_yaml::to_string(&layout)
                .unwrap_or_else(|e| panic!("Unable to serialize YAML: {}", e));
            std::fs::write(output, &yaml)
                .unwrap_or_else(|e| panic!("Unable to write {}: {}", output, e));

            // Validate the emitted file by re-parsing it through the same loader
            // the align/collapse commands use (panics on an invalid layout).
            let reloaded = SequenceLayout::from_yaml(output);
            let reference = reloaded.references.values().next().unwrap();
            info!(
                "Wrote read-structure YAML to {} ({} reference(s), {} target(s), {} UMI(s))",
                output,
                reloaded.references.len(),
                reference.targets.len(),
                reference.umi_configurations.len()
            );
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
