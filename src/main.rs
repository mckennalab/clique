#[macro_use]
extern crate lazy_static;
extern crate needletail;
extern crate seq_io;
extern crate petgraph;
extern crate rand;
extern crate bio;

use std::fs::File;
use std::str;
use std::io::Write;

use bio::io::fasta;

use seq_io::fastq::{Reader,Record};
use clap::Parser;
use extractor::{extract_tagged_sequences, align_unknown_orientation_read};
use rayon::prelude::*;
use std::sync::{Arc, Mutex};

pub mod extractor;
pub mod knownlist;
pub mod sequencelayout;
mod sequenceclustering;
mod bronkerbosch;


#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    #[clap(long)]
    reference: String,

    #[clap(long)]
    output: String,

    #[clap(long)]
    read: String,

    #[clap(long, default_value_t = 1)]
    threads: usize,

    #[clap(long)]
    outputupper: bool,
}

fn main() {
    let parameters = Args::parse();
    let reader = fasta::Reader::from_file(parameters.reference).unwrap();
    let _reference = reader.records().next().unwrap().unwrap().clone();
    let reference = _reference.seq();
    let ref_string = reference.to_vec();;

    let output = Arc::new(Mutex::new(File::create(parameters.output).unwrap()));

    //let mut reader1 = Reader::new(parameters.read);
    let f = File::open(parameters.read).unwrap();
    let mut reader = Reader::new(f);

    reader.records().par_bridge().for_each(|x| {
        let rec = x.unwrap();
        let norm_seq1 = rec.seq().to_vec();
        let aligned_read1 = align_unknown_orientation_read(&norm_seq1, &ref_string);
        let features = extract_tagged_sequences(&aligned_read1.2, &aligned_read1.1);

        let output = Arc::clone(&output);
        let mut output = output.lock().unwrap();

        if parameters.outputupper {
            write!(output, ">read1_{}\n{}\n", features.iter().filter(|(s,_t)| **s != "r".to_string() && **s != "e".to_string()).map(|(s, t)| format!("{}{}", &**s, &**t)).collect::<Vec<_>>().join(","), features.get(&"r".to_string()).unwrap()).unwrap();
            write!(output, ">ref\n{}\n", features.get(&"e".to_string()).unwrap()).unwrap();
        } else {
            write!(output, ">read1_{}\n{}\n", features.iter().filter(|(s,_t)| **s != "r".to_string() && **s != "e".to_string()).map(|(s, t)| format!("{}{}", &**s, &**t)).collect::<Vec<_>>().join(","),&format!("{}", String::from_utf8_lossy(aligned_read1.2.as_slice()))).unwrap();
            write!(output, ">ref\n{}\n", &format!("{}", String::from_utf8_lossy(aligned_read1.1.as_slice())));
        }
    });
}
