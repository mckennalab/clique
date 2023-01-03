use std::fs::{File, OpenOptions};
use bio::io::fastq::{Error, Record, Records};
use bio::io::fastq::Reader as FqReader;
use serde::{Serialize, Deserialize};
use std::fs;
use bincode::*;
use std::io::{BufReader, Write};
use std::path::PathBuf;
use rust_htslib::bgzf::{Reader, Writer};

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

unsafe impl Send for ReadIterator {}

unsafe impl Sync for ReadIterator {}

pub struct ReadIterator {
    read_one: Option<Records<BufReader<Reader>>>,
    read_two: Option<Records<BufReader<Reader>>>,
    index_one: Option<Records<BufReader<Reader>>>,
    index_two: Option<Records<BufReader<Reader>>>,

    pub reads_processed: usize,
    pub broken_reads: usize,
}


impl ReadIterator
{
    pub fn new(read_1: PathBuf,
               read_2: Option<PathBuf>,
               index_1: Option<PathBuf>,
               index_2: Option<PathBuf>
               ,
    ) -> ReadIterator {
        let r_one = ReadIterator::open_reader(&Some(&read_1));
        if !r_one.is_some() { panic!("Unable to open input file"); }

        let read2 = if read_2.is_some() { ReadIterator::open_reader(&Some(&read_2.clone().unwrap())) } else { None };
        let index1 = if index_1.is_some() { ReadIterator::open_reader(&Some(&index_1.clone().unwrap())) } else { None };
        let index2 = if index_2.is_some() { ReadIterator::open_reader(&Some(&index_2.clone().unwrap())) } else { None };

        let read2_active = read2.is_some();
        let index1_active = index1.is_some();
        let index2_active = index2.is_some();

        ReadIterator {
            read_one: r_one,
            read_two: read2,
            index_one: index1,
            index_two: index2,
            reads_processed: 0,
            broken_reads: 0,
        }
    }

    fn open_reader(check_path: &Option<&PathBuf>) -> Option<Records<BufReader<Reader>>> {
        if check_path.is_some() && check_path.as_ref().unwrap().exists() {
            let mut bgr = FqReader::new(Reader::from_path(&check_path.unwrap()).unwrap());
            let records = bgr.records();
            Some(records)
        } else {
            None
        }
    }

}
/// A helper function to manage unwrapping reads, which has a lot of layers
pub fn unwrap_reader(read: &mut Option<Records<BufReader<Reader>>>) -> Option<Record> {
    match read
    {
        Some(ref mut read_pointer) => {
            let next = read_pointer.next();
            match next {
                Some(rp) => {
                    match rp {
                        Ok(x) => Some(x),
                        Err(x) => None,
                    }
                }
                None => None,
            }
        }
        None => None,
    }
}

impl Iterator for ReadIterator {
    type Item = ReadSetContainer;

    fn next(&mut self) -> Option<ReadSetContainer> {
        let next_read_one = self.read_one.as_mut().unwrap().next();
        if next_read_one.is_some() {
            match next_read_one.unwrap() {
                Ok(v) => {
                    Some(ReadSetContainer {
                        read_one: v,
                        read_two: unwrap_reader(&mut self.read_two),
                        index_one: unwrap_reader(&mut self.index_one),
                        index_two: unwrap_reader(&mut self.index_two),
                    })
                }
                Err(e) => {
                    println!("error parsing header: {:?}", e);
                    None
                }
            }
        } else {
            None
        }
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