use bio::io::fastq::{Record, Records};
use bio::io::fastq::Reader as FqReader;
use serde::{Serialize, Deserialize};
use std::io::{BufReader};
use std::path::PathBuf;
use rust_htslib::bgzf::Reader;

/// holds a set of reads for reading and writing to disk
#[derive(Serialize, Deserialize, Debug, PartialEq)]
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
    #[allow(dead_code)]
    pub fn new_from_read1(rec: Record) -> ReadSetContainer {
        ReadSetContainer {
            read_one: rec,
            read_two: None,
            index_one: None,
            index_two: None,
        }
    }
}

impl std::fmt::Display for ReadSetContainer {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let res = write!(f, "{}", &self.read_one);
        if let Some(x) = &self.read_two {
            write!(f, "{}", x).expect("Unable to write ReadSetContainer for read 2");
        }
        if let Some(x) = &self.index_one {
            write!(f, "{}", x).expect("Unable to write ReadSetContainer for index 1");
        }
        if let Some(x) = &self.index_two {
            write!(f, "{}", x).expect("Unable to write ReadSetContainer for index 2");
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
            info!("Opening file: {:?}", check_path);
            Some(FqReader::new(Reader::from_path(&check_path.unwrap()).unwrap()).records())
        } else {
            info!("Unable to open file: {:?}", check_path);
            None
        }
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
            info!("Done processing reads");
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
                        Err(_) => None,
                    }
                }
                None => None,
            }
        }
        None => None,
    }
}
