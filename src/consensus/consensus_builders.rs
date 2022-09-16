extern crate rust_spoa;

use rust_spoa::poa_consensus;
use crate::read_strategies::sequence_structures::*;
use crate::read_strategies::sequence_layout::*;


pub struct ConsensusCandidate {
    pub reads: Vec<ReadSetContainer>,
    pub global_umi: Vec<u8>,
}

pub struct ConsensusResult {
    pub read_one: Vec<u8>,
    pub read_two: Option<Vec<u8>>,
    pub global_umi: Vec<u8>,
}

pub fn null_cap(strs: Vec<Vec<u8>>) -> Vec<Vec<u8>> {
    strs.iter().map(|r| {
        let mut ret = r.clone();
        ret.push(b'\0');
        ret
    }).collect::<Vec<Vec<u8>>>()
}


pub fn create_seq_layout_poa_consensus(reads: &Vec<ReadSetContainer>) -> SequenceSetContainer {
    let mut read1_agg = Vec::new();
    let mut read2_agg = Vec::new();
    let mut read3_agg = Vec::new();
    let mut read4_agg = Vec::new();

    &reads.iter().for_each(|n| {
        read1_agg.push(n.read_one.seq().clone().to_vec());
        n.read_two.as_ref().map(|r| read2_agg.push(r.seq().clone().to_vec()));
        n.index_one.as_ref().map(|r| read3_agg.push(r.seq().clone().to_vec()));
        n.index_two.as_ref().map(|r| read4_agg.push(r.seq().clone().to_vec()));
    });


    let read_one_conc: Vec<u8> = create_poa_consensus(&null_cap(read1_agg)).into_iter().collect::<Vec<u8>>();
    let read_two_conc: Option<Vec<u8>> = if read2_agg.len() > 0 {
        Some(create_poa_consensus(&null_cap(read2_agg)).into_iter().collect::<Vec<u8>>())
    } else {
        None
    };
    let read_three_conc: Option<Vec<u8>> = if read3_agg.len() > 0 {
        Some(create_poa_consensus(&null_cap(read3_agg)).into_iter().collect::<Vec<u8>>())
    } else {
        None
    };
    let read_four_conc: Option<Vec<u8>> = if read4_agg.len() > 0 {
        Some(create_poa_consensus(&null_cap(read4_agg)).into_iter().collect::<Vec<u8>>())
    } else {
        None
    };
    SequenceSetContainer {
        read_one: read_one_conc,
        read_two: read_two_conc,
        index_one: read_three_conc,
        index_two: read_four_conc
    }
}


pub fn create_poa_consensus(sequences: &Vec<Vec<u8>>) -> Vec<u8> {
    let max_length_vec = &sequences.iter().map(|n|n.len()).collect::<Vec<usize>>();
    let max_length = max_length_vec.iter().max().unwrap();
    let consensus = poa_consensus(&sequences, *max_length, 1, 5, -4, -3, -1);
    consensus
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::io::fastq::Record;

    #[test]
    fn simple_poa() {
        let reference = String::from(  "AATGATACGG\0").as_bytes().to_owned();
        let test_read = String::from(  "TATGATAAGG\0").as_bytes().to_owned();
        let test_read1 = String::from( "TATGAAGG\0").as_bytes().to_owned();
        let test_read2 = String::from( "TATGATAAGG\0").as_bytes().to_owned();
        let mut all = Vec::new();
        all.push(reference);
        all.push(test_read);
        all.push(test_read1);
        all.push(test_read2);


        assert_eq!(create_poa_consensus(&all),"TATGATAAGG".as_bytes().to_owned());
    }

    #[test]
    fn basic_read_consensus() {
        let read1 = ReadSetContainer::new_from_read1(Record::with_attrs("id_str", Some("desc"), b"ATGCGGG", b"QQQQQQQ"));
        let read2 = ReadSetContainer::new_from_read1(Record::with_attrs("id_str", Some("desc"), b"ATGCGGG", b"QQQQQQQ"));
        let read3 = ReadSetContainer::new_from_read1(Record::with_attrs("id_str", Some("desc"), b"ATACGGG", b"QQQQQQQ"));
        let read_vec = vec![&read1,&read2,&read3];
        let conc = create_seq_layout_poa_consensus(read_vec);
        assert_eq!(b"ATGCGGG".to_vec(),conc.read_one);
    }
}