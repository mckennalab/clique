
use bio::io::fastq::*;
use core::clone::Clone;
use core::option::Option;
use core::option::Option::{None, Some};
use std::fs::File;
use std::path::{Path, PathBuf};
use crate::sorters::sorter::ReadSortingOnDiskContainer;
use std::io::BufReader;
use std::io;
use flate2::read::GzEncoder;
use std::borrow::BorrowMut;
use bgzip::BGZFReader;
use std::collections::HashMap;

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
            read_two: if self.read_two.as_ref().is_some() {Some(self.read_two.as_ref().unwrap().clone())} else {None},
            index_one: if self.index_one.as_ref().is_some() {Some(self.index_one.as_ref().unwrap().clone())} else {None},
            index_two: if self.index_two.as_ref().is_some() {Some(self.index_two.as_ref().unwrap().clone())} else {None}
        }
    }
}

pub struct SequenceSetContainer {
    pub read_one: Vec<u8>,
    pub read_two: Option<Vec<u8>>,
    pub index_one: Option<Vec<u8>>,
    pub index_two: Option<Vec<u8>>,
}

impl Clone for SequenceSetContainer {
    fn clone(&self) -> SequenceSetContainer {
        SequenceSetContainer {
            read_one: self.read_one.clone(),
            read_two: if self.read_two.as_ref().is_some() {Some(self.read_two.as_ref().unwrap().clone())} else {None},
            index_one: if self.index_one.as_ref().is_some() {Some(self.index_one.as_ref().unwrap().clone())} else {None},
            index_two: if self.index_two.as_ref().is_some() {Some(self.index_two.as_ref().unwrap().clone())} else {None}
        }
    }
}

impl ReadSetContainer {
    pub fn new_from_read1(rec: Record) -> ReadSetContainer {
        ReadSetContainer{
            read_one: rec,
            read_two: None,
            index_one: None,
            index_two: None
        }
    }
}

impl std::fmt::Display for ReadSetContainer {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let res = write!(f,"{}",&self.read_one);
        if let Some(x) = &self.read_two {
            write!(f,"{}",x);
        }
        if let Some(x) = &self.index_one {
            write!(f,"{}",x);
        }
        if let Some(x) = &self.index_two {
            write!(f,"{}",x);
        }
        res
    }
}

pub struct ReadFileContainer {
    pub read_one: PathBuf,
    pub read_two: PathBuf,
    pub index_one: PathBuf,
    pub index_two: PathBuf,
}

pub enum ReadPattern {
    READ1, READ2, INDEX1, INDEX2
}


pub struct ReadIterator {
    pub read_one: Records<BufReader<BGZFReader<BufReader<File>>>>,
    pub read_two: Option<Records<BufReader<BGZFReader<BufReader<File>>>>>,
    pub index_one: Option<Records<BufReader<BGZFReader<BufReader<File>>>>>,
    pub index_two: Option<Records<BufReader<BGZFReader<BufReader<File>>>>>,

    pub reads_processed: usize,
    pub read2: bool,
    pub index1: bool,
    pub index2: bool,
    pub broken_reads: usize,
}

pub fn unwrap_read(read: &mut Option<Records<BufReader<BGZFReader<BufReader<File>>>>>) -> Option<Record> {
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
                },
                None => None,
            }
        },
        None => None,
    }
}

impl Iterator for ReadIterator {
    type Item = ReadSetContainer;


    fn next(&mut self) -> Option<ReadSetContainer> {
        let next_read_one = self.read_one.next();
        if next_read_one.is_some() {
            let ret = match next_read_one.unwrap() {
                Ok(v) => {
                    Some(ReadSetContainer {
                        read_one: v,
                        read_two:  unwrap_read(&mut self.read_two),
                        index_one: unwrap_read(&mut self.index_one),
                        index_two: unwrap_read(&mut self.index_two),
                    })
                }
                Err(e) => {
                    println!("error parsing header: {e:?}");
                    None
                }
            };
            match ret {
                Some(x) => {
                    if self.read2 && !self.read_two.is_some() {
                        self.broken_reads += 1;
                        None
                    } else if self.index1 && !self.index_one.is_some() {
                        self.broken_reads += 1;
                        None
                    } else if self.index2 && !self.index_two.is_some() {
                        self.broken_reads += 1;
                        None
                    } else {
                        self.reads_processed += 1;
                        Some(x)
                    }
                }
                None => None
            }
        } else {
            None
        }
    }


}

impl ReadIterator
{
    pub fn new(read_1: &Path,
               read_2: &Path,
               index_1: &Path,
               index_2: &Path,
    ) -> ReadIterator  {
        let r_one = ReadIterator::open_reader(&Some(read_1));
        if !r_one.is_some() {panic!("Unable to open input file");}

        let read2 = ReadIterator::open_reader(&Some(read_2));
        let index1 = ReadIterator::open_reader(&Some(index_1));
        let index2 = ReadIterator::open_reader(&Some(index_2));

        let read2_active = read2.is_some();
        let index1_active = index1.is_some();
        let index2_active = index2.is_some();

        ReadIterator {
            read_one: r_one.unwrap(),
            read_two: read2,
            index_one: index1,
            index_two: index2,
            reads_processed: 0,
            read2: read2_active,
            index1: index1_active,
            index2: index2_active,
            broken_reads: 0
        }
    }

    pub fn new_from_bundle(read_files: &ReadFileContainer) -> ReadIterator  {
        println!("opening read 1 {:?}",&read_files.read_one);
        ReadIterator::new(&read_files.read_one, &read_files.read_two, &read_files.index_one, &read_files.index_two)
    }

    pub fn new_from_on_disk_sorter(read_sorter: &ReadSortingOnDiskContainer) -> ReadIterator {
        let path1 = read_sorter.file_1.as_path();
        let path2 = if read_sorter.file_2.is_some() {Some(read_sorter.file_2.as_ref().unwrap().as_path())} else {None};
        let path3 = if read_sorter.file_3.is_some() {Some(read_sorter.file_3.as_ref().unwrap().as_path())} else {None};
        let path4 = if read_sorter.file_4.is_some() {Some(read_sorter.file_4.as_ref().unwrap().as_path())} else {None};

        let read2 = ReadIterator::open_reader(&path2);
        let index1 = ReadIterator::open_reader(&path3);
        let index2 = ReadIterator::open_reader(&path4);

        let read2_active = read2.is_some();
        let index1_active = index1.is_some();
        let index2_active = index2.is_some();

        ReadIterator {
            read_one: ReadIterator::open_reader(&Some(&path1)).unwrap(),
            read_two: read2,
            index_one: index1,
            index_two: index2,
            reads_processed: 0,
            read2: read2_active,
            index1: index1_active,
            index2: index2_active,
            broken_reads: 0
        }
    }

    fn open_reader(check_path: &Option<&Path>) -> Option<Records<BufReader<BGZFReader<BufReader<File>>>>> {
        if check_path.is_some() && check_path.as_ref().unwrap().exists() {
            println!("Opening {}",check_path.as_ref().unwrap().to_str().unwrap());
            let mut bgr = BGZFReader::new(BufReader::new(File::open(check_path.unwrap().to_str().unwrap()).unwrap()));
            let mut f2gz = Reader::new(bgr);
            let records = f2gz.records();
            Some(records)
        } else {
            None
        }
    }
}
