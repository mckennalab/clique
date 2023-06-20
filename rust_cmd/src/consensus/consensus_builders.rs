extern crate rust_spoa;

use std::borrow::Borrow;
use std::cmp::Ordering;
use std::collections::VecDeque;
use rayon::prelude::*;
use rust_spoa::poa_consensus;
use shardio::{Range, ShardReader, ShardWriter};
use crate::alignment::fasta_bit_encoding::FastaBase;
use crate::read_strategies::read_disk_sorter::{SortedAlignment, SortingReadSetContainer};
use std::io::Write;


pub struct ConsensusCandidate {
    pub reads: Vec<SortingReadSetContainer>,
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

pub fn write_consensus_reads(reader: &ShardReader<SortingReadSetContainer>, output_file: &String, threads: &usize, levels: usize) {
    //let mut sharded_output: ShardWriter<SortingReadSetContainer> = ShardWriter::new(&output_file, 32,
    //                                                                                256,
    //                                                                                1 << 16).unwrap();
    //let mut sender = sharded_output.get_sender();
    warn!("Writing consensus reads to {}", output_file);
    let mut output_fl = std::fs::File::create(output_file).expect("Unable to create file");

    let mut last_read: Option<SortingReadSetContainer> = None;
    let mut buffered_reads = VecDeque::new();


    reader.iter_range(&Range::all()).unwrap().enumerate().for_each(|(i, x)| {
        if i % 100000 == 0 { println!("processing read {}", i); }
        let x = x.unwrap();
        assert_eq!(x.ordered_sorting_keys.len(), levels);
        if last_read.is_some() && &x.cmp(last_read.as_ref().unwrap()) == &Ordering::Equal {
            if i < 10000 { println!("processing read {}", FastaBase::to_string(&x.aligned_read.aligned_read)); }
            buffered_reads.push_back(x.clone());
            last_read = Some(x);
        } else {
            if i < 50 && last_read.is_some() {
                let tags = x.ordered_sorting_keys.iter().map(|x|FastaBase::to_string(&x.1)).collect::<Vec<String>>().join(",");
                let tags2 = last_read.as_ref().unwrap().ordered_sorting_keys.iter().map(|x|FastaBase::to_string(&x.1)).collect::<Vec<String>>().join(",");

                println!("processing read {} tags {} tags2 {} ", FastaBase::to_string(&x.aligned_read.aligned_read),tags,tags2);
            }
            if buffered_reads.len() > 0 {
                let mut consensus_reads = create_poa_consensus(&buffered_reads);
                let ref_seq = FastaBase::to_string(&consensus_reads.aligned_read.aligned_ref);
                let read_seq = FastaBase::to_string(&consensus_reads.aligned_read.aligned_read);
                let count = buffered_reads.len();
                let tags = buffered_reads.pop_front().unwrap().ordered_sorting_keys.iter().map(|x|FastaBase::to_string(&x.1)).collect::<Vec<String>>().join(",");

                write!(output_fl, ">{};{}\n{}\n{}\n", tags, count, ref_seq, read_seq).unwrap(); ;
                buffered_reads.clear();
            }
            buffered_reads.push_back(x.clone());
            last_read = Some(x)
        }
    });
    output_fl.flush();


//                output_poa_consensus(Box::new(xx.reads), output);
}

pub fn create_poa_consensus(sequences: &VecDeque<SortingReadSetContainer>) -> SortingReadSetContainer {
    let max_length = sequences.iter().map(|n| n.aligned_read.aligned_read.len()).collect::<Vec<usize>>();
    let max_length = max_length.iter().max().unwrap();
    let base_sequences = sequences.iter().map(|n| {
        let mut y = FastaBase::to_vec_u8(&n.aligned_read.aligned_read);
        y.push(b'\0');
        y
    }).collect::<Vec<Vec<u8>>>();
    let consensus = poa_consensus(&base_sequences, max_length.clone(), 1, 5, -4, -3, -1);
    SortingReadSetContainer {
        ordered_sorting_keys: sequences[0].ordered_sorting_keys.clone(),
        ordered_unsorted_keys: sequences[0].ordered_unsorted_keys.clone(),
        aligned_read: SortedAlignment {
            aligned_read: FastaBase::from_vec_u8(&consensus),
            aligned_ref: sequences[0].aligned_read.aligned_ref.clone(),
            ref_name: sequences[0].aligned_read.ref_name.clone(),
        },
    }
}
/*

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

    let mut read1_agg = if read1_agg.len() > 100 {
        &read1_agg[read1_agg.len() - 100..read1_agg.len()]
    } else {
        &read1_agg
    };
    let mut read2_agg = if read2_agg.len() > 100 {
        &read2_agg[read2_agg.len() - 100..read2_agg.len()]
    } else {
        &read2_agg
    };
    let mut read3_agg = if read3_agg.len() > 100 {
        &read3_agg[read3_agg.len() - 100..read3_agg.len()]
    } else {
        &read3_agg
    };
    let mut read4_agg = if read4_agg.len() > 100 {
        &read4_agg[read4_agg.len() - 100..read4_agg.len()]
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
}*/