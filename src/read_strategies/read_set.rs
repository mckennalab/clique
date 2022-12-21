use std::fs::{File, OpenOptions};
use bio::io::fastq::{Error, Record, Records};
use bio::io::fastq::Reader as FqReader;
use serde::{Serialize, Deserialize};
use std::fs;
use bincode::*;
use std::io::Write;

/// holds a set of reads for reading and writing to disk
#[derive(Serialize, Deserialize, Debug)]
pub struct ReadSetContainer {
    pub read_one: Record,
    pub read_two: Option<Record>,
    pub index_one: Option<Record>,
    pub index_two: Option<Record>,
}

impl Clone for ReadSetContainer {
    fn clone(&self) -> ReadSetContainer {
        ReadSetContainer {
            read_one: self.read_one.clone(),
            read_two: if self.read_two.as_ref().is_some() { Some(self.read_two.as_ref().unwrap().clone()) } else { None },
            index_one: if self.index_one.as_ref().is_some() { Some(self.index_one.as_ref().unwrap().clone()) } else { None },
            index_two: if self.index_two.as_ref().is_some() { Some(self.index_two.as_ref().unwrap().clone()) } else { None },
        }
    }
}

impl ReadSetContainer {
    pub fn new_from_read1(rec: Record) -> ReadSetContainer {
        ReadSetContainer {
            read_one: rec,
            read_two: None,
            index_one: None,
            index_two: None,
        }
    }
    pub fn new_from_read2(rec: Record, old: &ReadSetContainer) -> ReadSetContainer {
        ReadSetContainer {
            read_one: old.read_one.clone(),
            read_two: Some(rec),
            index_one: old.index_one.clone(),
            index_two: old.index_two.clone(),
        }
    }
}

impl std::fmt::Display for ReadSetContainer {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let res = write!(f, "{}", &self.read_one);
        if let Some(x) = &self.read_two {
            write!(f, "{}", x);
        }
        if let Some(x) = &self.index_one {
            write!(f, "{}", x);
        }
        if let Some(x) = &self.index_two {
            write!(f, "{}", x);
        }
        res
    }
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
    pub fn from(input_file: &str) -> Option<ReadSetNestingContainer> {
        let mut file_contents =
            fs::read(&input_file).
                expect(&format!("Unable to open encoded ReadSetNestingContainer binary file: {}", &input_file));


        // here's where we get a little trickier -- we want to test if the file is, in order:
        // 1) a serialized ReadSetNestingContainer written
        // 2) a flat iterator of serialized ReadSetContainer
        // 3) error out
        let decoded_try: Result<Option<ReadSetNestingContainer>> = bincode::deserialize(&file_contents);
        match decoded_try {
            Ok(x) => { return (x); }
            Err(_) => {None}
        }
    }

    pub fn to(filename: &str, input_container: ReadSetNestingContainer) -> std::io::Result<()> {
        let encoded: Vec<u8> = bincode::serialize(&Some(input_container)).unwrap();

        std::fs::write(filename, encoded)
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str;
    use crate::utils::read_utils::fake_reads;
// find_max_value_2d_array(matrix: Array::<f64, Ix2>) -> Option<(AlignmentLocation,f64)>

    #[test]
    fn simple_read_set_container() {
        let fake_reads = fake_reads(80, 8);
        let srsc_reads = ReadSetNestingContainer { nest: None, unsorted: Some(fake_reads) };
    }

    #[test]
    fn nested_read_set_container() {
        let fake_reads = fake_reads(80, 8);
        let srsc_reads1 = ReadSetNestingContainer { nest: None, unsorted: Some(fake_reads) };

        let srsc2 = ReadSetNestingContainer { nest: Some(Box::new(srsc_reads1)), unsorted: None };
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
}