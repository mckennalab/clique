extern crate rust_spoa;

use std::borrow::Borrow;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use bio::io::fastq::Record;
use flate2::{Compression, GzBuilder};
use rayon::prelude::{IntoParallelIterator, IntoParallelRefIterator};
use rayon::prelude::*;
use rust_htslib::bgzf::Writer;
use rust_spoa::poa_consensus;

use crate::read_strategies::sequence_file_containers::*;
use crate::read_strategies::sequence_layout::*;
use crate::RunSpecifications;
use crate::sorters::known_list::KnownListBinSplit;
use crate::umis::sequence_clustering::*;

pub struct ConsensusCandidate {
    pub reads: Vec<ReadSetContainer>,
    pub global_umi: Vec<u8>,
}

pub struct ConsensusResult {
    pub read_one: Vec<u8>,
    pub read_two: Option<Vec<u8>>,
    pub global_umi: Vec<u8>,
}

pub fn null_cap(strs: &[Vec<u8>]) -> Vec<Vec<u8>> {
    strs.iter().map(|r| {
        let mut ret = r.clone();
        ret.push(b'\0');
        ret
    }).collect::<Vec<Vec<u8>>>()
}

pub fn threaded_write_consensus_reads(read_iterators: Vec<SuperClusterOnDiskIterator>, output_file_base: &String, pattern: &ReadPattern, run_spec: &RunSpecifications) {
    let mut output_writer = OutputReadSetWriter::from_pattern(&PathBuf::from(output_file_base), pattern);

    let output = Arc::new(Mutex::new(output_writer));
    let pool = rayon::ThreadPoolBuilder::new().num_threads(run_spec.processing_threads).build().unwrap();
    let pooled_install = pool.install(|| {
        //let iters = read_iterators.into_iter();
        for read_iterator in read_iterators {
            read_iterator.into_iter().par_bridge().for_each(|xx| { //.par_bridge()
                let output = Arc::clone(&output);
                output_poa_consensus(Box::new(xx.reads), output);
            });
        }
    });
}

pub fn create_poa_consensus(sequences: &Vec<Vec<u8>>) -> Vec<u8> {
    let max_length_vec = &sequences.iter().map(|n| n.len()).collect::<Vec<usize>>();
    let max_length = max_length_vec.iter().max().unwrap();
    let consensus = poa_consensus(&sequences, *max_length, 1, 5, -4, -3, -1);
    consensus
}


/// Create a consensus sequence from an input set of reads. We merge reads by their underlying
/// sequence, not the ReadLayout we get from the transformed reads
pub fn output_poa_consensus(reads: Box<dyn Iterator<Item=ReadSetContainer>>, output: Arc<Mutex<OutputReadSetWriter>>) {
    let mut read1_agg = Vec::new();
    let mut read2_agg = Vec::new();
    let mut read3_agg = Vec::new();
    let mut read4_agg = Vec::new();

    let mut read_name: Option<Vec<u8>> = None;
    let mut read_count: usize = 0;
    &reads.for_each(|n| {
        if !read_name.is_some() { read_name = Some(n.read_one.id().clone().as_bytes().to_vec()) }
        read1_agg.push(n.read_one.seq().clone().to_vec());
        n.read_two.as_ref().map(|r| read2_agg.push(r.seq().clone().to_vec()));
        n.index_one.as_ref().map(|r| read3_agg.push(r.seq().clone().to_vec()));
        n.index_two.as_ref().map(|r| read4_agg.push(r.seq().clone().to_vec()));
        read_count += 1;
    });

    let mut read1_agg= if read1_agg.len() > 100 {
        &read1_agg[read1_agg.len()-100..read1_agg.len()]
    } else {
        &read1_agg
    };
    let mut read2_agg= if read2_agg.len() > 100 {
        &read2_agg[read2_agg.len()-100..read2_agg.len()]
    } else {
        &read2_agg
    };
    let mut read3_agg= if read3_agg.len() > 100 {
        &read3_agg[read3_agg.len()-100..read3_agg.len()]
    } else {
        &read3_agg
    };
    let mut read4_agg= if read4_agg.len() > 100 {
        &read4_agg[read4_agg.len()-100..read4_agg.len()]
    } else {
        &read4_agg
    };

    let read_one_conc: Record = to_read(read_name.as_ref().unwrap(),
                                        &create_poa_consensus(&null_cap(read1_agg)).into_iter().collect::<Vec<u8>>(),
                                        &1, &read_count);

    let read_two_conc: Option<Record> = if read2_agg.len() > 0 {
        Some(to_read(read_name.as_ref().unwrap(),
                     &create_poa_consensus(&null_cap(read2_agg)).into_iter().collect::<Vec<u8>>(),
                     &2, &read_count))
    } else {
        None
    };
    let read_three_conc: Option<Record> = if read3_agg.len() > 0 {
        Some(to_read(read_name.as_ref().unwrap(),
                     &create_poa_consensus(&null_cap(read3_agg)).into_iter().collect::<Vec<u8>>(),
                     &3, &read_count))
    } else {
        None
    };
    let read_four_conc: Option<Record> = if read4_agg.len() > 0 {
        Some(to_read(read_name.as_ref().unwrap(),
                     &create_poa_consensus(&null_cap(read4_agg)).into_iter().collect::<Vec<u8>>(),
                     &4, &read_count))
    } else {
        None
    };

    let conc = ReadSetContainer {
        read_one: read_one_conc,
        read_two: read_two_conc,
        index_one: read_three_conc,
        index_two: read_four_conc,
    };

    let mut output_unwrapped = output.lock().unwrap();
    output_unwrapped.write(&conc);
}

pub fn to_read(base_name: &Vec<u8>, sequence: &Vec<u8>, read_position: &usize, read_count: &usize) -> Record {
    let qual = sequence.iter().map(|n| b'H').collect::<Vec<u8>>();
    let mut name = String::from_utf8(base_name.clone()).unwrap();
    name.push_str(&"-read_");
    name.push_str(&read_position.to_string());
    name.push_str(&"-count_");
    name.push_str(&read_count.to_string());
    Record::with_attrs(&name, None, sequence, qual.as_slice())
}


#[cfg(test)]
mod tests {
    use bio::io::fastq::Record;

    use super::*;

    #[test]
    fn simple_poa() {
        let reference = String::from("AATGATACGG\0").as_bytes().to_owned();
        let test_read = String::from("TATGATAAGG\0").as_bytes().to_owned();
        let test_read1 = String::from("TATGAAGG\0").as_bytes().to_owned();
        let test_read2 = String::from("TATGATAAGG\0").as_bytes().to_owned();
        let mut all = Vec::new();
        all.push(reference);
        all.push(test_read);
        all.push(test_read1);
        all.push(test_read2);


        assert_eq!(create_poa_consensus(&all), "TATGATAAGG".as_bytes().to_owned());
    }
}