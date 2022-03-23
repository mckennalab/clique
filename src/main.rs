#[macro_use]
extern crate lazy_static;
extern crate needletail;
extern crate seq_io;

use std::borrow::Borrow;
use std::collections::HashMap;
use std::fs::File;
use std::str;
use std::sync::mpsc;
use std::time::{Duration, SystemTime};
use std::io::Write;
use std::sync::mpsc::sync_channel;
use std::thread;

use needletail::{parse_fastx_file, Sequence};

use bio::io::fasta;

use seq_io::fastq::{Reader,Record};
use std::io;
use clap::Parser;
use extractor::{align_unknown_orientation_read_u8_ref, extract_tagged_sequences, align_unknown_orientation_read};
use rayon::prelude::*;
use std::sync::{Arc, Mutex};

pub mod extractor;
pub mod knownlist;
pub mod sequencelayout;


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
}

fn main() {
    let parameters = Args::parse();
    let reader = fasta::Reader::from_file(parameters.reference).unwrap();
    let _reference = reader.records().next().unwrap().unwrap().clone();
    let reference = _reference.seq();
    let ref_string = str::from_utf8(&reference).unwrap().to_string();
    let counter = Arc::new(Mutex::new(0));


    let output = Arc::new(Mutex::new(File::create(parameters.output).unwrap()));

    //let mut reader1 = Reader::new(parameters.read);
    let mut f = File::open(parameters.read).unwrap();
    let mut reader = Reader::new(f);

    reader.records().par_bridge().for_each(|x| {
        let rec = x.unwrap();
        let norm_seq1 = str::from_utf8(rec.seq()).unwrap().to_string();
        let aligned_read1 = align_unknown_orientation_read(&norm_seq1, &ref_string);
        let features = extract_tagged_sequences(&aligned_read1.2, &aligned_read1.1);

        let output = Arc::clone(&output);
        let mut output = output.lock().unwrap();

        write!(output, ">read1_{}\n{}\n", features.iter().map(|(s,t)| format!("{}{}",&**s,&**t)).collect::<Vec<_>>().join(","), aligned_read1.2).unwrap();
        write!(output, ">ref\n{}\n", aligned_read1.1).unwrap();
    });
}
