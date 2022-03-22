extern crate bio;
extern crate clap;
extern crate rand;
extern crate needletail;

#[macro_use]
extern crate lazy_static;

use clap::{App, Arg};
use bio::io::fasta;
use std::fs::File;
use needletail::{parse_fastx_file, Sequence};
use extractor::{align_unknown_orientation_read_u8_ref, extract_tagged_sequences};
use std::borrow::Borrow;
use std::str;

use std::io::Write;
use std::time::{SystemTime,Duration};

pub mod extractor;
pub mod knownlist;
pub mod sequencelayout;

fn main() -> std::io::Result<()> {
    let matches = App::new("Cabal")
        .version("1.0")
        .author("Aaron M. <aaronatwpi@gmail.com>")
        .about("Collect and collapse UMIs")
        .arg(Arg::with_name("reference")
            .short("r")
            .long("ref")
            .value_name("FILE")
            .help("The reference we will align to")
            .required(true)
            .takes_value(true))
        .arg(Arg::with_name("output")
            .short("o")
            .long("output")
            .value_name("FILE")
            .help("output reads")
            .required(true)
            .takes_value(true))
        .arg(Arg::with_name("read")
            .short("r")
            .long("read")
            .value_name("FILE")
            .help("the reads")
            .required(true)
            .takes_value(true))
        .get_matches();


    let reference_file = matches.value_of("reference");
    let reader = fasta::Reader::from_file(reference_file.unwrap()).unwrap();
    let _reference = reader.records().next().unwrap().unwrap().clone();
    let reference = _reference.seq();
    let ref_string = str::from_utf8(&reference).unwrap().to_string();


    let mut output = File::create(matches.value_of("output").unwrap()).unwrap();


    let mut reader1 = parse_fastx_file(matches.value_of("read").unwrap()).expect("valid path/file");

    let mut counter = 0;
    let now = SystemTime::now();

    while let Some(record) = reader1.next() {
        counter += 1;
        let norm_seq1 = record.unwrap().normalize(false).into_owned();
        let aligned_read1 = align_unknown_orientation_read_u8_ref(norm_seq1.borrow(), reference);
        let features = extract_tagged_sequences(&aligned_read1.2, &aligned_read1.1);

        write!(output, ">read1\n{}\n", aligned_read1.2).unwrap();

        write!(output, ">ref\n{}\n", aligned_read1.1).unwrap();
        write!(output, ">special\n{}\n", features.iter().map(|(s,t)| format!("{}{}",&**s,&**t)).collect::<Vec<_>>().join(", ")).unwrap();
        if counter % 100 == 0 {
            println!("Counter {}, elapsed {}",counter,now.elapsed().unwrap().as_secs());
        }
    }

    Ok(())
}
