use std::borrow::BorrowMut;
use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::io;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::slice::Chunks;

use bgzip::BGZFReader;
use bio::io::fastq::{Error, Record, Records};
use bio::io::fastq::Reader as FqReader;
use flate2::*;
use flate2::GzBuilder;
use flate2::read::GzEncoder;
use num_traits::FromPrimitive;
use rust_htslib::bgzf::{Reader, Writer};

use core::clone::Clone;
use core::option::Option;
use core::option::Option::{None, Some};
use serde::{Deserialize, Serialize};

use crate::itertools::Itertools;
use tempfile::TempPath;

#[derive(Clone)]
pub enum ReadPattern {
    ONE,
    ONETWO,
    ONETWOI1,
    ONETWOI2,
    ONETWOI1I2,
    ONEI1,
    ONEI2,
    ONEI1I2,
}

impl ReadPattern {
    pub fn contains_r2(&self) -> bool {
        match self {
            ReadPattern::ONE => { false }
            ReadPattern::ONETWO => { true }
            ReadPattern::ONETWOI1 => { true }
            ReadPattern::ONETWOI2 => { true }
            ReadPattern::ONETWOI1I2 => { true }
            ReadPattern::ONEI1 => { false }
            ReadPattern::ONEI2 => { false }
            ReadPattern::ONEI1I2 => { false }
        }
    }
    pub fn contains_i1(&self) -> bool {
        match self {
            ReadPattern::ONE => { false }
            ReadPattern::ONETWO => { false }
            ReadPattern::ONETWOI1 => { true }
            ReadPattern::ONETWOI2 => { false }
            ReadPattern::ONETWOI1I2 => { true }
            ReadPattern::ONEI1 => { true }
            ReadPattern::ONEI2 => { false }
            ReadPattern::ONEI1I2 => { true }
        }
    }
    pub fn contains_i2(&self) -> bool {
        match self {
            ReadPattern::ONE => { false }
            ReadPattern::ONETWO => { false }
            ReadPattern::ONETWOI1 => { false }
            ReadPattern::ONETWOI2 => { true }
            ReadPattern::ONETWOI1I2 => { true }
            ReadPattern::ONEI1 => { false }
            ReadPattern::ONEI2 => { true }
            ReadPattern::ONEI1I2 => { true }
        }
    }

    pub fn from_read_iterator(iter: &ReadIterator) -> ReadPattern {
        ReadPattern::pattern(iter.path_two.is_some(), iter.index_one.is_some(), iter.index_two.is_some())
    }

    pub fn from_read_file_container(set: &ReadFileContainer) -> ReadPattern {
        ReadPattern::pattern(set.read_two.is_some(), set.index_one.is_some(), set.index_two.is_some())
    }

    pub fn from_read_set_container(set: &ReadSetContainer) -> ReadPattern {
        ReadPattern::pattern(set.read_two.is_some(), set.index_one.is_some(), set.index_two.is_some())
    }

    fn pattern(r2: bool, i1: bool, i2: bool) -> ReadPattern {
        match (r2, i2, i2) {
            (true, true, true) => ReadPattern::ONETWOI1I2,
            (false, true, true) => ReadPattern::ONEI1I2,
            (false, false, true) => ReadPattern::ONEI2,
            (false, true, false) => ReadPattern::ONEI1,
            (false, false, false) => ReadPattern::ONE,
            (true, false, true) => ReadPattern::ONETWOI2,
            (true, false, false) => ReadPattern::ONETWO,
            (true, true, false) => ReadPattern::ONETWOI1,
        }
    }
}

/// holds a set of reads for reading and writing to disk
#[derive(Serialize, Deserialize)]
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
}

pub struct NamedReadSetContainer {
    pub sort_string: Vec<u8>,
    pub reads: ReadSetContainer,
}

impl Clone for NamedReadSetContainer {
    fn clone(&self) -> NamedReadSetContainer {
        NamedReadSetContainer {
            sort_string: self.sort_string.clone(),
            reads: self.reads.clone(),
        }
    }
}

// with_attrs(id: &str, desc: Option<&str>, seq: TextSlice<'_>, qual: &[u8])
//#[derive(Serialize, Deserialize)]
//#[serde(remote = "Record")]
struct RecordDef {
    //#[serde(getter = "Record::id")]
    id: String,
    //#[serde(getter = "Record::seq")]
    seq: Vec<u8>,
    //#[serde(getter = "Record::qual")]
    qual: Vec<u8>,
}

// Provide a conversion to construct the remote type.
impl From<RecordDef> for Record {
    fn from(def: RecordDef) -> Record {
        Record::with_attrs(def.id.as_str(), None, def.seq.as_slice(), def.qual.as_slice())
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

pub fn estimate_read_count(location: &String) -> Option<usize> {
    let path = PathBuf::from(location);
    let mut reader = ReadIterator::open_reader(&Some(&path));
    let file_size = path.as_path().metadata().unwrap().len();

    if reader.is_some() {
        let mut first_100_reads_sizes = 0;
        for i in 0..100 {
            let read = unwrap_reader(&mut reader);
            match read {
                None => { return None; }
                Some(rd) => {
                    first_100_reads_sizes += rd.seq().len() * 2;
                    first_100_reads_sizes += rd.id().len();
                }
            }
        }
        let average = first_100_reads_sizes as f64 / 100.0;
        Some(((file_size as f64) / ((average as f64) * 0.15)) as usize) //
    } else {
        None
    }

}


/// Where we can find either input or output fastq files
pub struct ReadFileContainer {
    pub read_one: PathBuf,
    pub read_two: Option<PathBuf>,
    pub index_one: Option<PathBuf>,
    pub index_two: Option<PathBuf>,
}

impl Clone for ReadFileContainer {
    fn clone(&self) -> ReadFileContainer {
        ReadFileContainer {
            read_one: self.read_one.clone(),
            read_two: self.read_two.clone(),
            index_one: self.index_one.clone(),
            index_two: self.index_two.clone(),
        }
    }
}

impl ReadFileContainer {
    pub fn new(read1: &String, read2: &String, index1: &String, index2: &String) -> ReadFileContainer {
        let path1 = PathBuf::from(read1);
        let path2 = PathBuf::from(read2);
        let path3 = PathBuf::from(index1);
        let path4 = PathBuf::from(index2);

        ReadFileContainer {
            read_one: path1,
            read_two: if path2.exists() { Some(path2) } else { None },
            index_one: if path3.exists() { Some(path3) } else { None },
            index_two: if path4.exists() { Some(path4) } else { None },
        }
    }
    pub fn new_from_paths(path1: &PathBuf, path2: &PathBuf, path3: &PathBuf, path4: &PathBuf) -> ReadFileContainer {
        ReadFileContainer {
            read_one: path1.clone(),
            read_two: if path2.exists() { Some(path2.clone()) } else { None },
            index_one: if path3.exists() { Some(path3.clone()) } else { None },
            index_two: if path4.exists() { Some(path4.clone()) } else { None },
        }
    }

    fn format_temp_read_name(index: usize, prefix: String) -> String { format!("{}_{}.fastq.gz", prefix, index).to_string() }
/*
    pub fn new_temp_files(tmp_generator: fn() -> TempPath) -> ReadFileContainer {
        let mut read1 = tmp_generator();
        let mut read2 = tmp_generator();
        let mut read3 = tmp_generator();
        let mut read4 = tmp_generator();

        ReadFileContainer {
            read_one: read1,
            read_two: if rd.read_two.is_some() { Some(temp_dir.join(ReadFileContainer::format_temp_read_name(id, read2))) } else { None },
            index_one: if rd.index_one.is_some() { Some(temp_dir.join(ReadFileContainer::format_temp_read_name(id, read3))) } else { None },
            index_two: if rd.index_two.is_some() { Some(temp_dir.join(ReadFileContainer::format_temp_read_name(id, read4))) } else { None },
        }
*/
    pub fn new_from_temp_dir(rd: &ReadIterator, prefix: &String, id: usize, temp_dir: &Path) -> ReadFileContainer {
        let mut read1 = prefix.clone();
        read1.push_str("_read1");
        let mut read2 = prefix.clone();
        read2.push_str("_read2");
        let mut read3 = prefix.clone();
        read3.push_str("_read3");
        let mut read4 = prefix.clone();
        read4.push_str("_read4");

        ReadFileContainer {
            read_one: temp_dir.join(ReadFileContainer::format_temp_read_name(id, read1)),
            read_two: if rd.read_two.is_some() { Some(temp_dir.join(ReadFileContainer::format_temp_read_name(id, read2))) } else { None },
            index_one: if rd.index_one.is_some() { Some(temp_dir.join(ReadFileContainer::format_temp_read_name(id, read3))) } else { None },
            index_two: if rd.index_two.is_some() { Some(temp_dir.join(ReadFileContainer::format_temp_read_name(id, read4))) } else { None },
        }
    }
}


/// manage writing read sets to disk
pub struct OutputReadSetWriter {
    file_1: flate2::write::GzEncoder<File>,
    file_2: Option<flate2::write::GzEncoder<File>>,
    file_3: Option<flate2::write::GzEncoder<File>>,
    file_4: Option<flate2::write::GzEncoder<File>>,

    written_read1: usize,
    written_read2: usize,
    written_read3: usize,
    written_read4: usize,
}

impl OutputReadSetWriter {
    pub fn format_read_name(index: usize, prefix: String) -> String { format!("{}_{}.txt", prefix, index).to_string() }

    pub fn from_read_file_container(sc: &ReadFileContainer) -> OutputReadSetWriter {
        OutputReadSetWriter {
            file_1: OutputReadSetWriter::create_writer(&sc.read_one),//.unwrap(),
            file_2: if let Some(x) = &sc.read_two { Some(OutputReadSetWriter::create_writer(&x.clone())) } else { None },
            file_3: if let Some(x) = &sc.index_one { Some(OutputReadSetWriter::create_writer(&x.clone())) } else { None },
            file_4: if let Some(x) = &sc.index_two { Some(OutputReadSetWriter::create_writer(&x.clone())) } else { None },
            written_read1: 0,
            written_read2: 0,
            written_read3: 0,
            written_read4: 0,
        }
    }


    pub fn from_pattern(base: &PathBuf, pt: &ReadPattern) -> OutputReadSetWriter {
        let base_path = base.as_path().to_str().unwrap();
        let read1: PathBuf = [base_path, "read1.fq.gz"].iter().collect();
        println!("READ 1 = {:?}", read1);
        OutputReadSetWriter {
            file_1: OutputReadSetWriter::create_writer(&read1),
            file_2: if pt.contains_r2() { Some(OutputReadSetWriter::create_writer(&[base_path, "read2.fq.gz"].iter().collect())) } else { None },
            file_3: if pt.contains_i1() { Some(OutputReadSetWriter::create_writer(&[base_path, "index1.fq.gz"].iter().collect())) } else { None },
            file_4: if pt.contains_i2() { Some(OutputReadSetWriter::create_writer(&[base_path, "index2.fq.gz"].iter().collect())) } else { None },
            written_read1: 0,
            written_read2: 0,
            written_read3: 0,
            written_read4: 0,
        }
    }

    fn create_writer(filename: &PathBuf) -> flate2::write::GzEncoder<File> {
        let f = File::create(filename).unwrap();
        let mut gz = GzBuilder::new()
            .write(f, Compression::best());
        gz
    }

    pub fn write(&mut self, rl: &ReadSetContainer) {
        //self.file_1.borrow_mut().write(&rl.read_one.to_string().as_bytes());
        write!(self.file_1.borrow_mut(), "{}", rl.read_one);
        self.written_read1 += 1;
        if let Some(x) = self.file_2.as_mut() {
            if let Some(rd) = &rl.read_two {
                write!(x, "{}", rd);
                //x.write(&rd.to_string().as_bytes());
                self.written_read2 += 1;
            };
        };
        if let Some(x) = self.file_3.as_mut() {
            if let Some(rd) = &rl.index_one {
                write!(x, "{}", rd);
                //x.write(&rd.to_string().as_bytes());
                self.written_read3 += 1;
            };
        };
        if let Some(x) = self.file_4.as_mut() {
            if let Some(rd) = &rl.index_two {
                write!(x, "{}", rd);
                //x.write(&rd.to_string().as_bytes());
                self.written_read4 += 1;
            };
        };
    }
    pub fn print_read_count(&self) {
        println!("Read 1 {}, read 2 {} read 3 {} read 4 {}", self.written_read1, self.written_read2, self.written_read3, self.written_read4);
    }

    pub fn create_x_bins(rd: &ReadIterator, prefix: &String, x_bins: usize, temp_dir: &Path) -> Vec<(usize, ReadFileContainer)> {
        (0..x_bins).map(|id|
            (id, ReadFileContainer::new_from_temp_dir(rd, prefix, id, temp_dir))).collect::<Vec<(usize, ReadFileContainer)>>()
    }
}

unsafe impl Send for ReadIterator {}

unsafe impl Sync for ReadIterator {}

pub struct ReadIterator {
    read_one: Option<Records<BufReader<Reader>>>,
    read_two: Option<Records<BufReader<Reader>>>,
    index_one: Option<Records<BufReader<Reader>>>,
    index_two: Option<Records<BufReader<Reader>>>,

    read_collection: Option<ClusteredReads>,

    path_one: PathBuf,
    path_two: Option<PathBuf>,
    path_i1: Option<PathBuf>,
    path_i2: Option<PathBuf>,

    pub reads_processed: usize,
    pub read2: bool,
    pub index1: bool,
    pub index2: bool,
    pub broken_reads: usize,
}

/// A helper function to manage unwrapping reads, which has a lot of conditions
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
        match self.read_collection.as_mut() {
            Some(mut x) => {
                x.reads.pop_front()
            }
            None => {
                let next_read_one = self.read_one.as_mut().unwrap().next();
                if next_read_one.is_some() {
                    let ret = match next_read_one.unwrap() {
                        Ok(v) => {
                            Some(ReadSetContainer {
                                read_one: v,
                                read_two: unwrap_reader(&mut self.read_two),
                                index_one: unwrap_reader(&mut self.index_one),
                                index_two: unwrap_reader(&mut self.index_two),
                            })
                        }
                        Err(e) => {
                            println!("error parsing header: {e:?}");
                            None
                        }
                    };
                    match ret {
                        Some(x) => {
                            if self.read2 && !x.read_two.is_some() {
                                self.broken_reads += 1;
                                None
                            } else if self.index1 && !x.index_one.is_some() {
                                self.broken_reads += 1;
                                None
                            } else if self.index2 && !x.index_two.is_some() {
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
    }
}

impl ReadIterator
{
    pub fn from_collection(read_collection: ClusteredReads) -> ReadIterator {
        ReadIterator {
            read_one: None,
            read_two: None,
            index_one: None,
            index_two: None,
            read_collection: Some(read_collection),
            path_one: Default::default(),
            path_two: None,
            path_i1: None,
            path_i2: None,
            reads_processed: 0,
            read2: false,
            index1: false,
            index2: false,
            broken_reads: 0,
        }
    }
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
            read_collection: None,
            path_one: read_1.clone(),
            path_two: read_2.clone(),
            path_i1: index_1.clone(),
            path_i2: index_2.clone(),
            reads_processed: 0,
            read2: read2_active,
            index1: index1_active,
            index2: index2_active,
            broken_reads: 0,
        }
    }

    pub fn new_from_bundle(read_files: &ReadFileContainer) -> ReadIterator {
        println!("opening read 1 {:?}", &read_files.read_one);
        ReadIterator::new(read_files.read_one.clone(), read_files.read_two.clone(), read_files.index_one.clone(), read_files.index_two.clone())
    }

    fn open_reader(check_path: &Option<&PathBuf>) -> Option<Records<BufReader<Reader>>> {
        if check_path.is_some() && check_path.as_ref().unwrap().exists() {
            println!("Opening {}", check_path.as_ref().unwrap().to_str().unwrap());
            let mut bgr = FqReader::new(Reader::from_path(&check_path.unwrap()).unwrap());
            let records = bgr.records();
            Some(records)
        } else {
            None
        }
    }

    pub fn new_reset(&self) -> Self {
        ReadIterator::new(self.path_one.clone(), self.path_two.clone(), self.path_i1.clone(), self.path_i2.clone())
    }
}

/// This is ugly, but I don't have the energy for some dyn nightmare right now
pub struct ClusteredReadIterator {
    reads: VecDeque<ClusteredReads>,
    read_files: VecDeque<ReadFileContainer>,
    read_pattern: ReadPattern,
}

impl Iterator for ClusteredReadIterator {
    type Item = ClusteredReads;

    fn next(&mut self) -> Option<ClusteredReads> {
        match self.reads.len() {
            0 => {
                while self.read_files.len() > 0 && self.reads.len() == 0 {
                    match self.read_files.pop_front() {
                        None => {}
                        Some(x) => {
                            let mut readset = VecDeque::new();
                            for r in ReadIterator::new_from_bundle(&x) {
                                readset.push_back(r);
                            }
                            self.reads.push_back(ClusteredReads { reads: readset, pattern: self.read_pattern.clone(), average_distance: None })
                        }
                    }
                }
                self.reads.pop_front()
            }
            _ => self.reads.pop_front(),
        }
    }
}

impl ClusteredReadIterator {
    pub fn new_from_vec(reads: Vec<ClusteredReads>, read_pattern: ReadPattern) -> ClusteredReadIterator {
        ClusteredReadIterator {
            reads: VecDeque::from(reads),
            read_files: VecDeque::new(),
            read_pattern,
        }
    }

    pub fn new_from_files(reads: VecDeque<ReadFileContainer>, read_pattern: ReadPattern) -> ClusteredReadIterator {
        ClusteredReadIterator {
            reads: VecDeque::new(),
            read_files: reads,
            read_pattern,
        }
    }
}

pub struct ClusteredReads {
    reads: VecDeque<ReadSetContainer>,
    pattern: ReadPattern,
    average_distance: Option<f64>,
}

impl ClusteredReads {
    pub fn new(reads: VecDeque<ReadSetContainer>, pattern: ReadPattern, average_dist: f64) -> ClusteredReads {
        ClusteredReads { reads, pattern, average_distance: Some(average_dist) }
    }
}
