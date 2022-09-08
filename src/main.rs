#[macro_use]
extern crate lazy_static;
extern crate needletail;
extern crate seq_io;
extern crate petgraph;
extern crate rand;
extern crate bio;
extern crate flate2;
extern crate suffix;
extern crate ndarray;
extern crate num_traits;
extern crate fastq;
extern crate noodles_fastq;
extern crate rust_spoa;

use std::fs::File;
use std::str;
use std::io::{Write, BufReader, BufRead};
use seq_io::fasta::{Reader, Record, OwnedRecord};
use noodles_fastq as Fastq;


use std::io;

use clap::Parser;
use extractor::extract_tagged_sequences;
use rayon::prelude::*;
use std::sync::{Arc, Mutex};
use std::collections::{BTreeMap, HashMap};
use bio::alignment::Alignment;
//use flate2::GzBuilder;
//use flate2::Compression;

use alignment::alignment_matrix::*;
use alignment::scoring_functions::*;

mod linked_alignment;
use crate::linked_alignment::*;
use umis::sequenceclustering::*;
use petgraph::algo::connected_components;
use std::path::Path;
use reference::fasta_reference::reference_file_to_struct;
use read_strategies::sequence_layout::{ReadIterator, LayoutType, transform};
use std::str::FromStr;

pub mod extractor;
mod simple_umi_clustering;

mod umis {
    pub mod bronkerbosch;
    pub mod sequenceclustering;
}

mod alignment {
    pub mod alignment_matrix;
    pub mod scoring_functions;
}


mod consensus {
    pub mod serial_passage_read_corrector;
    pub mod consensus_builders;
}

pub mod fasta_comparisons;

mod read_strategies {
    pub mod sequence_layout;
    pub mod ten_x;
}

mod reference {
    pub mod fasta_reference;
}


#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    #[clap(long)]
    reference: String,

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
}
/*
fn old_main() {
    let mut test_set = Vec::new();

    if let Ok(lines) = read_lines("test_data/just_sequences_20000_16s.txt") {
        // Consumes the iterator, returns an (Optional) String
        for line in lines {
            if let Ok(ip) = line {
                test_set.push(ip.as_bytes().to_vec());
            }
        }
    }
    println!("TEST SIZE {}, making graph",test_set.len());
    let graph = input_list_to_graph(&InputList{strings: test_set, max_dist: 1},string_distance, true);
    println!("Making clique");
    println!("Connected size: {}", connected_components(&graph.graph));
    process_cliques(&graph);
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
    where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())


    let my_score = InversionScoring {
        match_score: 10.0,
        mismatch_score: -11.0,
        gap_open: -15.0,
        gap_extend: -5.0,
        inversion_penalty: -2.0,
    };


    let my_aff_score = AffineScoring {
        match_score: 10.0,
        mismatch_score: -11.0,
        special_character_score: 7.0,
        gap_open: -15.0,
        gap_extend: -5.0,
    };
}*/

fn main() {
    let parameters = Args::parse();

    let read_layout = LayoutType::from_str(&parameters.read_template).expect("Unable to parse read template type");

    let reference = reference_file_to_struct(&parameters.reference);

    let output_file = File::create(parameters.output).unwrap();

    let reference_lookup = find_seeds(&reference.name,20);
    let output = Arc::new(Mutex::new(output_file)); // Mutex::new(gz));

    // setup our thread pool
    rayon::ThreadPoolBuilder::new().num_threads(parameters.threads).build_global().unwrap();

    let mut known_list = load_knownlist(&parameters.known_list);

    let read_iterator = ReadIterator::new(parameters.read1, parameters.read2,parameters.index1,parameters.index2);

    if read_layout.has_umi() {
        let mut cnt = 0;
        let mut read_mapping = HashMap::new();
        for rd in read_iterator {

            if cnt % 100 == 0 {
                println!("Count {}, hash size: {}",cnt, read_mapping.len());
            }
            cnt += 1;
            let transformed_reads = transform(rd, &read_layout);
            let first_hit = transformed_reads.get_unique_sequences().unwrap()[0].clone();
            //if known_list.as_ref().is_some() {
                let corrected_hits = correct_to_known_list(&first_hit, &mut known_list, 1);
                read_mapping.insert(first_hit.clone(),corrected_hits);
            //} else {
            //    read_mapping.insert(first_hit.clone(),BestHits{ hits: vec![first_hit.clone()], distance: 0 });
            //}
        }
        println!("Read mapping size {}",read_mapping.len());
    }

    //read_iterator.par_bridge().for_each(|xx|
}
/*let x = xx.unwrap().clone();
            let name = &String::from_utf8_lossy(&x.name()).to_string();

            let is_forward = orient_by_longest_segment(&x.sequence().to_vec(), &reference.sequence, &reference_lookup);

            if is_forward.0 {
                let fwd_score_mp = find_greedy_non_overlapping_segments(&x.sequence().to_vec(), &reference.sequence, &reference_lookup);
                let results = align_string_with_anchors(&x.sequence().to_vec(), &reference.sequence, &fwd_score_mp, &my_score,&my_aff_score);

                let output = Arc::clone(&output);
                let mut output = output.lock().unwrap();
                write!(output,">ref{}\n{}\n>{}\n{}\n",str::from_utf8(&reference.sequence).unwrap(),str::from_utf8(&results.1).unwrap(),str::replace(name," ","_"),str::from_utf8(&results.0).unwrap()).expect("Unable to write to output file");
                output.flush().expect("Unable to flush output");
            } else {
                let fwd_score_mp = find_greedy_non_overlapping_segments(&x.sequence().to_vec(), &reference.sequence, &reference_lookup);

                let results = align_string_with_anchors(&reverse_complement(&x.sequence().to_vec()), &reference.sequence, &fwd_score_mp, &my_score,&my_aff_score);

                let output = Arc::clone(&output);
                let mut output = output.lock().unwrap();
                write!(output,">ref{}\n{}\n>{}\n{}\n",str::from_utf8(&reference.sequence).unwrap(),str::from_utf8(&results.0).unwrap(),str::replace(name," ","_"),str::from_utf8(&results.1).unwrap()).expect("Unable to write to output file");
                output.flush().expect("Unable to flush output");
            }*/

#[allow(dead_code)]
struct AlignedWithFeatures {
    alignment: Alignment,
    read_id: String,
    read: Vec<u8>,
    reference: Vec<u8>,
    features: BTreeMap<String,String>
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


