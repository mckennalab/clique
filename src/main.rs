#[macro_use]
extern crate lazy_static;
extern crate needletail;
extern crate seq_io;
extern crate petgraph;
extern crate rand;
extern crate bio;
extern crate flate2;
extern crate rust_htslib;

use std::fs::File;
use std::str;
use std::io::Write;
use rust_htslib::bgzf::Reader as HtslibReader;
use seq_io::fasta::{Reader, Record, OwnedRecord};
use seq_io::fastq::Reader as Fastq;
use seq_io::fastq::Record as FastqRecord;

use std::io;

use clap::Parser;
use extractor::{extract_tagged_sequences, align_unknown_orientation_read};
use rayon::prelude::*;
use std::sync::{Arc, Mutex};
use std::collections::BTreeMap;
use bio::alignment::Alignment;
use bio::io::fasta::Records;
use std::borrow::BorrowMut;
use std::io::prelude::*;
use flate2::read::GzDecoder;
use flate2::GzBuilder;
use flate2::Compression;


pub mod extractor;
pub mod knownlist;
pub mod sequencelayout;
mod sequenceclustering;
mod bronkerbosch;
mod simple_umi_clustering;


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


    let output_file = File::create(parameters.output).unwrap();
    let mut gz = GzBuilder::new()
        .comment("aligned fasta file")
        .write(output_file, Compression::best());
    let output = Arc::new(Mutex::new(gz));

    // open read one
    //let f1 = File::open(parameters.read1).unwrap();
    let f1gz = HtslibReader::from_path(parameters.read1).unwrap();
    let mut readers = Readers{first: Some(Fastq::new(f1gz)), second: None};

    // check for a second read file
    let f2 = File::open(parameters.read2.clone());
    if f2.is_ok() {
        let f2gz = HtslibReader::from_path(parameters.read2).unwrap();
        readers.second = Some(Fastq::new(f2gz));
    }

    // setup our thread pool
    rayon::ThreadPoolBuilder::new().num_threads(parameters.threads).build_global().unwrap();

    // This is a little ugly since we're wrapping in the para_bridge, which introduces some scope issues
    if readers.second.is_some() {
        readers.first.unwrap().records().zip(readers.second.unwrap().records()).for_each(|(xx, yy)| {
            let x = xx.unwrap().clone();
            let y = yy.unwrap().clone();
            let alignment1 = process_read_into_enriched_obj(&x.id().unwrap().to_string(), &x.seq().to_vec(), &ref_string);
            let alignment2 = process_read_into_enriched_obj(&y.id().unwrap().to_string(), &y.seq().to_vec(), &ref_string);

            let output = Arc::clone(&output);
            let mut output = output.lock().unwrap();

            write!(output,"{}",to_two_line_fasta(alignment1,parameters.outputupper));
            write!(output,"{}",to_two_line_fasta(alignment2,parameters.outputupper));
        });
    } else {
        readers.first.unwrap().records().for_each(|(xx)| {
            let x = xx.unwrap().clone();
            let alignment1 = process_read_into_enriched_obj(&x.id().unwrap().to_string(),&x.seq().to_vec(), &ref_string);

            let output = Arc::clone(&output);
            let mut output = output.lock().unwrap();

            write!(output,"{}",to_two_line_fasta(alignment1,parameters.outputupper));
        });
    }
}

struct Readers<R: io::Read> {
    first: Option<Fastq<R>>,
    second:Option<Fastq<R>>,
}

struct AlignedWithFeatures {
    alignment: Alignment,
    read_id: String,
    read: Vec<u8>,
    reference: Vec<u8>,
    features: BTreeMap<String,String>
}

fn to_two_line_fasta(alignFeatures: AlignedWithFeatures, output_upper: bool) -> String {
    if output_upper {
        format!(">@{}_{}\n{}\n>ref\n{}\n",alignFeatures.read_id,
                alignFeatures.features.iter().filter(|(s, _t)| **s != "r".to_string() && **s != "e".to_string()).map(|(s, t)| format!("{}{}", &**s, &**t)).collect::<Vec<_>>().join(","),
                alignFeatures.features.get(&"r".to_string()).unwrap(),
                alignFeatures.features.get(&"e".to_string()).unwrap())

    } else {
        format!(">@{}_{}\n{}\n>ref\n{}\n",
                alignFeatures.read_id,
                alignFeatures.features.iter().filter(|(s, _t)| **s != "r".to_string() && **s != "e".to_string()).map(|(s, t)| format!("{}{}", &**s, &**t)).collect::<Vec<_>>().join(","),
                &format!("{}", String::from_utf8_lossy(alignFeatures.read.as_slice())),
                &format!("{}", String::from_utf8_lossy(alignFeatures.reference.as_slice())))
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
