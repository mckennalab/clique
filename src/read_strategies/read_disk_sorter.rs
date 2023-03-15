use std::cmp::Ordering;
use std::fs;
use std::fs::File;
use flate2::bufread::GzEncoder;
use crate::read_strategies::read_set::ReadSetContainer;
use serde::ser::{SerializeSeq, Serializer};
use serde::{Serialize, Deserialize};
use tempfile::TempPath;
use crate::alignment::fasta_bit_encoding::FastaBase;

#[derive(Serialize, Deserialize, Debug)]
pub struct SortingReadSetContainer {
    sorting_keys: Vec<Vec<FastaBase>>,
    reads: ReadSetContainer,
}

impl Eq for SortingReadSetContainer {}

impl PartialEq<Self> for SortingReadSetContainer {
    fn eq(&self, other: &Self) -> bool {
        assert_eq!(self.sorting_keys.len(), other.sorting_keys.len());
        self.sorting_keys.iter().zip(other.sorting_keys.iter()).map(|(a, b)| a == b).count() == 0
    }
}

impl PartialOrd<Self> for SortingReadSetContainer {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        for (a, b) in self.sorting_keys.iter().zip(other.sorting_keys.iter()) {
            if a > b {
                return Some(Ordering::Greater);
            }
            else if a < b {
                return Some(Ordering::Less);
            }
        }
        Some(Ordering::Equal)
    }
}

pub struct SplitWriteAndSort {
    split_files: Vec<Option<GzEncoder<File>>>,

}

#[derive(Serialize, Deserialize, Debug)]
pub struct ReadSetNestingContainer {
    /// A structure that contains zero or more subcontainers, with an optional set of unsorted reads.
    /// We use this as a multilayer hierarchy of reads to be arranged on disk and then read back in for
    /// sorting / merging.
    nest: Option<Box<ReadSetNestingContainer>>,
    unsorted: Option<Vec<ReadSetContainer>>,
}

impl ReadSetNestingContainer {
    #[allow(dead_code)]
    pub fn from(input_file: &str) -> Option<ReadSetNestingContainer> {
        let file_contents =
            fs::read(&input_file).
                expect(&format!("Unable to open encoded ReadSetNestingContainer binary file: {}", &input_file));

        bincode::deserialize(&file_contents).unwrap()
    }

    #[allow(dead_code)]
    pub fn to(filename: &str, input_container: ReadSetNestingContainer) -> std::io::Result<()> {
        let encoded: Vec<u8> = bincode::serialize(&Some(input_container)).unwrap();

        std::fs::write(filename, encoded)
    }
}


pub trait ReadSetNestingContainerSortingTransform {
    /// sort the retained bin of reads
    fn sort_bin(container: &mut ReadSetNestingContainer) -> Box<ReadSetNestingContainer>;
}

/// all
pub trait ReadSetNestingContainerStream {
    /// Use this method to correct reads to known barcodes, and to build
    /// up a database of observed subsequences for UMIs
    fn transform_stream(read_set: &ReadSetContainer) -> ReadSetContainer;
}


#[cfg(test)]
mod tests {
    use crate::alignment::fasta_bit_encoding::{FASTA_A, FASTA_N};
    use super::*;
    use crate::utils::read_utils::fake_reads;

    #[test]
    fn simple_read_set_container() {
        let fake_reads = fake_reads(80, 8);
        let _srsc_reads = ReadSetNestingContainer { nest: None, unsorted: Some(fake_reads) };
    }

    #[test]
    fn nested_read_set_container() {
        let fake_reads = fake_reads(80, 8);
        let srsc_reads1 = ReadSetNestingContainer { nest: None, unsorted: Some(fake_reads) };

        let _srsc2 = ReadSetNestingContainer { nest: Some(Box::new(srsc_reads1)), unsorted: None };
    }

    #[test]
    fn simple_serialize_deserialize() {
        let fake_reads = fake_reads(80, 8);
        let srsc_reads = ReadSetNestingContainer { nest: None, unsorted: Some(fake_reads.clone()) };

        let encoded: Vec<u8> = bincode::serialize(&Some(srsc_reads)).unwrap();

        let decoded: Option<ReadSetNestingContainer> = bincode::deserialize(&encoded[..]).unwrap();

        match decoded {
            None => { panic!("Unable to unwrap the decoded string") }
            Some(x) => {
                assert_eq!(x.unsorted.unwrap().len(), fake_reads.len());
            }
        }
    }

    #[test]
    fn enbedded_serialize_deserialize() {
        let fake_reads = fake_reads(80, 8);
        let fake_read_size = fake_reads.len();
        let srsc_reads = ReadSetNestingContainer { nest: None, unsorted: Some(fake_reads.clone()) };

        // now another layer!
        let srsc_reads2 = ReadSetNestingContainer { nest: Some(Box::new(srsc_reads)), unsorted: None };

        let encoded: Vec<u8> = bincode::serialize(&Some(srsc_reads2)).unwrap();

        let decoded: Option<ReadSetNestingContainer> = bincode::deserialize(&encoded[..]).unwrap();

        match decoded {
            None => { panic!("Unable to unwrap the decoded string") }
            Some(x) => {
                assert_eq!(x.nest.unwrap().unsorted.unwrap().len(), fake_read_size);
            }
        }
    }

    #[test]
    fn to_disk_serialize_deserialize() {
        let fake_reads = fake_reads(80, 8);
        let fake_read_size = fake_reads.len();
        let srsc_reads = ReadSetNestingContainer { nest: None, unsorted: Some(fake_reads.clone()) };

        // now another layer!
        let srsc_reads2 = ReadSetNestingContainer { nest: Some(Box::new(srsc_reads)), unsorted: None };

        let filename = "test_data/test_rsnc.bin";
        ReadSetNestingContainer::to(filename, srsc_reads2).expect("Unable to write read set container");

        let decoded = ReadSetNestingContainer::from(&filename).unwrap();

        assert_eq!(decoded.nest.unwrap().unsorted.unwrap().len(), fake_read_size);
    }

    #[test]
    fn test_sorting_read_container() {
        let key1 = vec![FASTA_N,FASTA_A];
        let key2 = vec![FASTA_N,FASTA_N];

        let reads = fake_reads(10, 1);

        let st1 = SortingReadSetContainer{ sorting_keys: vec![key1.clone()], reads: reads.get(0).unwrap().clone() };
        let st2 = SortingReadSetContainer{ sorting_keys: vec![key2.clone()], reads: reads.get(1).unwrap().clone() };
        assert!(st1 < st2);

        let st1 = SortingReadSetContainer{ sorting_keys: vec![key2.clone()], reads: reads.get(0).unwrap().clone() };
        let st2 = SortingReadSetContainer{ sorting_keys: vec![key1.clone()], reads: reads.get(1).unwrap().clone() };
        assert!(st1 > st2);

        let st1 = SortingReadSetContainer{ sorting_keys: vec![key1.clone()], reads: reads.get(0).unwrap().clone() };
        let st2 = SortingReadSetContainer{ sorting_keys: vec![key1.clone()], reads: reads.get(1).unwrap().clone() };
        println!("{} {}",st1> st2, st1 < st2);
        assert!(!(st1 > st2) & !(st2 > st1));

        let st1 = SortingReadSetContainer{ sorting_keys: vec![key1.clone(),key2.clone()], reads: reads.get(0).unwrap().clone() };
        let st2 = SortingReadSetContainer{ sorting_keys: vec![key1.clone(),key1.clone()], reads: reads.get(1).unwrap().clone() };
        assert!(st1 > st2);

    }
}