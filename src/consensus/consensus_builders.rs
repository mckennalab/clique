extern crate rust_spoa;

use rust_spoa::poa_consensus;
use crate::read_strategies::sequence_layout::{SequenceLayout, ReadLayout};

pub struct ConsensusCandidate {
    pub reads: Vec<Box<dyn SequenceLayout>>,
    pub global_umi: Vec<u8>,
}



pub struct ConsensusResult {
    pub read_one: Vec<u8>,
    pub read_two: Option<Vec<u8>>,
    pub global_umi: Vec<u8>,
}

pub fn create_seq_layout_poa_consensus(candidates: &ConsensusCandidate) -> Vec<u8> {
    let read_one_conc: Vec<u8> = create_poa_consensus(&candidates.reads.iter().map(|n| {
        let mut temp = n.read_one().clone();
        temp.push(b'\0');
        temp
    }).into_iter().collect());


    let has_read_two = candidates.reads.iter().map(|r| r.read_two().is_some()).fold(true, |acc, mk| acc && mk);

    let read_two_conc: Option<Vec<u8>> = if has_read_two {Some(
        create_poa_consensus(&candidates.reads.iter().map(|n| {
        let mut temp = n.read_two().unwrap().clone();
        temp.push(b'\0');
        temp
    }).into_iter().collect()))} else {
        None
    };
    read_one_conc // TODO: yeah this isn't final..
}


pub fn create_poa_consensus(sequences: &Vec<Vec<u8>>) -> Vec<u8> {
    let consensus = poa_consensus(&sequences, 20, 1, 5, -4, -3, -1);
    consensus
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simple_poa_test() {
        let reference = String::from(  "AATGATACGG\0").as_bytes().to_owned();
        let test_read = String::from(  "TATGATAAGG\0").as_bytes().to_owned();
        let test_read1 = String::from( "TATGAAGG\0").as_bytes().to_owned();
        let test_read2 = String::from( "TATGATAAGG\0").as_bytes().to_owned();
        let mut all = Vec::new();
        all.push(reference);
        all.push(test_read);
        all.push(test_read1);
        all.push(test_read2);

        //let conc = ConsensusCandidate{ reads: : vec!["read1".as_bytes().to_vec(), "read2".as_bytes().to_vec(), "read3".as_bytes().to_vec()], read_sequences: all, global_umi: vec![] };

        //assert_eq!(create_poa_consensus(conc),"TATGATAAGG".as_bytes().to_owned());
    }
}