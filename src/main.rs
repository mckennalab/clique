extern crate actix;
#[macro_use]
extern crate lazy_static;
extern crate needletail;

use std::borrow::Borrow;
use std::collections::HashMap;
use std::fs::File;
use std::str;
use std::sync::mpsc;
use std::time::{Duration, SystemTime};

use actix::{Actor, Context, Handler, SyncArbiter, System};
use actix::prelude::*;
use needletail::{parse_fastx_file, Sequence};
use tokio::signal;

use alignment_actor::{Aligner, ReadRefAligner};
use bio::io::fasta;
use clap::Parser;
use extractor::{align_unknown_orientation_read_u8_ref, extract_tagged_sequences};

pub mod extractor;
pub mod knownlist;
pub mod sequencelayout;
pub mod alignment_actor;


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

// https://github.com/actix/actix/blob/HEAD/actix/examples/fibonacci.rs

fn main() {
    let mut system = System::new();
    let parameters = Args::parse();
    let reader = fasta::Reader::from_file(parameters.reference).unwrap();
    let _reference = reader.records().next().unwrap().unwrap().clone();
    let reference = _reference.seq();
    let ref_string = str::from_utf8(&reference).unwrap().to_string();


    let mut output = File::create(parameters.output).unwrap();


    let mut reader1 = parse_fastx_file(parameters.read).expect("valid path/file");

    system.block_on(async {
        let output_actor_addr = alignment_actor::AlignmentWriter { output }.start().recipient();

        let addr = SyncArbiter::start(parameters.threads, || Aligner);

        while let Some(record) = reader1.next() {
            println!("STUFF!");
            let norm_seq1 = str::from_utf8(record.unwrap().normalize(false).borrow()).unwrap().to_string();
            addr.do_send(ReadRefAligner {
                read: norm_seq1,
                reference: ref_string.clone(),
                output_actor: output_actor_addr.clone(),
            });
        }

        tokio::signal::ctrl_c().await.unwrap();
        println!("Ctrl-C received, shutting down");
        System::current().stop();
    });
}
