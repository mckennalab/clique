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

use std::fs::File;
use std::str;
use std::io::{Write, BufReader};
use seq_io::fasta::{Reader, Record, OwnedRecord};
use noodles_fastq as Fastq;


use std::io;

use clap::Parser;
use extractor::{extract_tagged_sequences, align_unknown_orientation_read};
use rayon::prelude::*;
use std::sync::{Arc, Mutex};
use std::collections::BTreeMap;
use bio::alignment::Alignment;
//use flate2::GzBuilder;
//use flate2::Compression;

use alignment::alignment_matrix::*;
use alignment::scoring_functions::*;

mod linked_alignment;
use crate::linked_alignment::*;
pub mod extractor;
mod simple_umi_clustering;


mod alignment {
    pub mod alignment_matrix;
    pub mod fasta_bit_encoding;
    pub mod scoring_functions;
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

    #[clap(long, default_value_t = 1)]
    threads: usize,

    #[clap(long)]
    outputupper: bool,
}


fn main() {
    let parameters = Args::parse();

    // reference loading / checking
    let mut reader = Reader::from_path(parameters.reference).unwrap();
    let fasta_entries: Vec<OwnedRecord> = reader.records().map(|f| f.unwrap()).collect();
    assert_eq!(fasta_entries.len(), 1, "We can only run with single entry FASTA files");
    let ref_string = fasta_entries.get(0).unwrap().seq.clone();
    let ref_name = fasta_entries.get(0).unwrap().id().unwrap();

    let output_file = File::create(parameters.output).unwrap();
    //let mut gz = GzBuilder::new()
     //   .comment("aligned fasta file")
     //   .write(output_file, Compression::best());
    //write!(gz,"@HD\tVN:1.6\n@SQ\tSN:{}\tLN:{}\n",ref_name,ref_string.len()).expect("Unable to write to output file");


    // open read one
    //let f1 = File::open(parameters.read1).unwrap();

    let f1gz = File::open(parameters.read1)
        .map(BufReader::new)
        .map(Fastq::Reader::new).unwrap();
    let mut readers = Readers{first: Some(f1gz), second: None};

    // check for a second read file
    let f2 = File::open(parameters.read2.clone());
    if f2.is_ok() {
        let f2gz = File::open(parameters.read2)
            .map(BufReader::new)
            .map(Fastq::Reader::new).unwrap();
        readers.second = Some(f2gz);
    }

    let reference_lookup = find_seeds(&ref_string,20);
    let output = Arc::new(Mutex::new(output_file)); // Mutex::new(gz));

    // setup our thread pool
    rayon::ThreadPoolBuilder::new().num_threads(parameters.threads).build_global().unwrap();


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
        gap_open: -15.0,
        gap_extend: -5.0,
    };

    // This is a little ugly since we're wrapping in the para_bridge, which introduces some scope issues
    if readers.second.is_some() {
        readers.first.unwrap().records().into_iter().zip(readers.second.unwrap().records()).par_bridge().for_each(|(xx, yy)| {
            let x = xx.unwrap().clone();
            let y = yy.unwrap().clone();
            let alignment1 = process_read_into_enriched_obj(&String::from_utf8_lossy(&x.name()).to_string(), &x.sequence().to_vec(), &ref_string);
            let alignment2 = process_read_into_enriched_obj(&String::from_utf8_lossy(&y.name()).to_string(), &y.sequence().to_vec(), &ref_string);

            let output = Arc::clone(&output);
            let mut output = output.lock().unwrap();

            write!(output,"{}",to_two_line_fasta(alignment1,parameters.outputupper)).expect("Unable to write to output file");
            write!(output,"{}",to_two_line_fasta(alignment2,parameters.outputupper)).expect("Unable to write to output file");
        });
    } else {
        readers.first.unwrap().records().par_bridge().for_each(|xx| {
            let x = xx.unwrap().clone();
            let name = &String::from_utf8_lossy(&x.name()).to_string();

            let is_forward = orient_by_longest_segment(&x.sequence().to_vec(), &ref_string, &reference_lookup);

            if is_forward.0 {
                let fwd_score_mp = find_greedy_non_overlapping_segments(&x.sequence().to_vec(), &ref_string, &reference_lookup);
                let results = align_string_with_anchors(&x.sequence().to_vec(), &ref_string, &fwd_score_mp, &my_score,&my_aff_score);
                //let alignment_string: String = alignment.alignment_tags.into_iter().map(|m| m.to_string()).collect();

                let output = Arc::clone(&output);
                let mut output = output.lock().unwrap();
                //write!(output,"{}\t0\t{}\t1\t250\t{}\t*\t0\t{}\t{}\t{}\n",str::replace(name," ","_"),ref_name,alignment_string,&x.sequence().len(),seq,qual).expect("Unable to write to output file");
                write!(output,">ref{}\n{}\n>{}\n{}\n",ref_name,str::from_utf8(&results.1).unwrap(),str::replace(name," ","_"),str::from_utf8(&results.0).unwrap()).expect("Unable to write to output file");
                output.flush();
            } else {
                let fwd_score_mp = find_greedy_non_overlapping_segments(&x.sequence().to_vec(), &ref_string, &reference_lookup);

                let results = align_string_with_anchors(&reverse_complement(&x.sequence().to_vec()), &ref_string, &fwd_score_mp, &my_score,&my_aff_score);

                let output = Arc::clone(&output);
                let mut output = output.lock().unwrap();
                //write!(output,"{}\t0\t{}\t1\t250\t{}\t*\t0\t{}\t{}\t{}\n",str::replace(name," ","_"),ref_name,alignment_string,&x.sequence().len(),seq,qual).expect("Unable to write to output file");
                write!(output,">ref{}\n{}\n>{}\n{}\n",ref_name,str::from_utf8(&results.0).unwrap(),str::replace(name," ","_"),str::from_utf8(&results.1).unwrap()).expect("Unable to write to output file");
                output.flush();
            }

        });
    }
}

struct Readers<R: io::Read> {
    first: Option<Fastq::Reader<R>>,
    second:Option<Fastq::Reader<R>>,
}

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


fn process_read_into_enriched_obj(read_id:&String, read_rec: &Vec<u8>, reference: &Vec<u8>) -> AlignedWithFeatures {
    let aligned_read = align_unknown_orientation_read(&read_rec, &reference);
    let features = extract_tagged_sequences(&aligned_read.2, &aligned_read.1);
    AlignedWithFeatures{
        alignment: aligned_read.0,
        read_id: read_id.clone(),
        read: aligned_read.2,
        reference: aligned_read.1,
        features: features.clone()
    }
}
