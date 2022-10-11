use std::{fmt, fs, io};
use std::borrow::BorrowMut;
use std::collections::{HashMap, VecDeque};
use std::fs::*;
use std::io::{BufRead, BufReader, BufWriter};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::slice::Chunks;
use std::str::FromStr;

use bgzip::BGZFReader;
use bio::io::fastq::{Error, Record, Records};
use bio::io::fastq::Reader as FqReader;
use flate2::*;
use flate2::GzBuilder;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use num_traits::FromPrimitive;
use rust_htslib::bgzf::{Reader, Writer};
use serde::{Deserialize, Serialize};
use serde::ser::{Serializer, SerializeSeq};
use tempfile::TempPath;

use core::clone::Clone;
use core::option::Option;
use core::option::Option::{None, Some};

use crate::itertools::Itertools;
use crate::RunSpecifications;

#[derive(Serialize, Deserialize, Debug, Clone)]
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

impl FromStr for ReadPattern {
    type Err = ();

    fn from_str(input: &str) -> Result<ReadPattern, Self::Err> {
        match input {
            "ONE" => Ok(ReadPattern::ONE),
            "ONETWO" => Ok(ReadPattern::ONETWO),
            "ONETWOI1" => Ok(ReadPattern::ONETWOI1),
            "ONETWOI2" => Ok(ReadPattern::ONETWOI2),
            "ONETWOI1I2" => Ok(ReadPattern::ONETWOI1I2),
            "ONEI1" => Ok(ReadPattern::ONEI1),
            "ONEI2" => Ok(ReadPattern::ONEI2),
            "ONEI1I2" => Ok(ReadPattern::ONEI1I2),
            _ => Err(()),
        }
    }
}

impl fmt::Display for ReadPattern {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl ReadPattern {
    pub fn read_count(&self) -> usize {
        match self {
            ReadPattern::ONE => { 1 }
            ReadPattern::ONETWO => { 2 }
            ReadPattern::ONETWOI1 => { 3 }
            ReadPattern::ONETWOI2 => { 3 }
            ReadPattern::ONETWOI1I2 => { 4 }
            ReadPattern::ONEI1 => { 2 }
            ReadPattern::ONEI2 => { 2 }
            ReadPattern::ONEI1I2 => { 3 }
        }
    }

    pub fn to_read_collection(&self, reads: &mut VecDeque<Record>) -> Option<ReadSetContainer> {
        if reads.len() != self.read_count() {
            None
        } else {
            Some(match self {
                ReadPattern::ONE => {
                    ReadSetContainer {
                        read_one: reads.pop_front().unwrap(),
                        read_two: None,
                        index_one: None,
                        index_two: None,
                    }
                }
                ReadPattern::ONETWO => {
                    ReadSetContainer {
                        read_one: reads.pop_front().unwrap(),
                        read_two: Some(reads.pop_front().unwrap()),
                        index_one: None,
                        index_two: None,
                    }
                }
                ReadPattern::ONETWOI1 => {
                    ReadSetContainer {
                        read_one: reads.pop_front().unwrap(),
                        read_two: Some(reads.pop_front().unwrap()),
                        index_one: Some(reads.pop_front().unwrap()),
                        index_two: None,
                    }
                }
                ReadPattern::ONETWOI2 => {
                    ReadSetContainer {
                        read_one: reads.pop_front().unwrap(),
                        read_two: Some(reads.pop_front().unwrap()),
                        index_one: None,
                        index_two: Some(reads.pop_front().unwrap()),
                    }
                }
                ReadPattern::ONETWOI1I2 => {
                    ReadSetContainer {
                        read_one: reads.pop_front().unwrap(),
                        read_two: Some(reads.pop_front().unwrap()),
                        index_one: Some(reads.pop_front().unwrap()),
                        index_two: Some(reads.pop_front().unwrap()),
                    }
                }
                ReadPattern::ONEI1 => {
                    ReadSetContainer {
                        read_one: reads.pop_front().unwrap(),
                        read_two: None,
                        index_one: Some(reads.pop_front().unwrap()),
                        index_two: None,
                    }
                }
                ReadPattern::ONEI2 => {
                    ReadSetContainer {
                        read_one: reads.pop_front().unwrap(),
                        read_two: None,
                        index_one: None,
                        index_two: Some(reads.pop_front().unwrap()),
                    }
                }
                ReadPattern::ONEI1I2 => {
                    ReadSetContainer {
                        read_one: reads.pop_front().unwrap(),
                        read_two: None,
                        index_one: Some(reads.pop_front().unwrap()),
                        index_two: Some(reads.pop_front().unwrap()),
                    }
                }
            })
        }
    }

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
    fn reverse_pattern(pat: &ReadPattern) -> (bool, bool, bool) {
        match pat {
            ReadPattern::ONETWOI1I2 => (true, true, true),
            ReadPattern::ONEI1I2 => (false, true, true),
            ReadPattern::ONEI2 => (false, false, true),
            ReadPattern::ONEI1 => (false, true, false),
            ReadPattern::ONE => (false, false, false),
            ReadPattern::ONETWOI2 => (true, false, true),
            ReadPattern::ONETWO => (true, false, false),
            ReadPattern::ONETWOI1 => (true, true, false),
        }
    }
}

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
        Some(((file_size as f64) / ((average as f64) * 0.14)) as usize) //
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

    pub fn temporary(rd: &ReadIterator, run_specs: &mut RunSpecifications) -> ReadFileContainer {
        ReadFileContainer {
            read_one: run_specs.create_temp_file(),
            read_two: if rd.read_two.is_some() { Some(run_specs.create_temp_file()) } else { None },
            index_one: if rd.index_one.is_some() { Some(run_specs.create_temp_file()) } else { None },
            index_two: if rd.index_two.is_some() { Some(run_specs.create_temp_file()) } else { None },
        }
    }

    pub fn temporary_from_pattern(pattern: &ReadPattern, run_specs: &mut RunSpecifications) -> ReadFileContainer {
        let file_layout = ReadPattern::reverse_pattern(pattern);
        ReadFileContainer {
            read_one: run_specs.create_temp_file(),
            read_two: if file_layout.0 { Some(run_specs.create_temp_file()) } else { None },
            index_one: if file_layout.1 { Some(run_specs.create_temp_file()) } else { None },
            index_two: if file_layout.2 { Some(run_specs.create_temp_file()) } else { None },
        }
    }
}


/// manage writing read sets to disk
pub struct OutputReadSetWriter {
    file_1: GzEncoder<File>,
    file_2: Option<GzEncoder<File>>,
    file_3: Option<GzEncoder<File>>,
    file_4: Option<GzEncoder<File>>,

    files: ReadFileContainer,

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
            files: sc.clone(),
            written_read1: 0,
            written_read2: 0,
            written_read3: 0,
            written_read4: 0,
        }
    }

    pub fn from_pattern(base: &PathBuf, pt: &ReadPattern) -> OutputReadSetWriter {
        let base_path = base.as_path().to_str().unwrap();
        let rfc = ReadFileContainer {
            read_one: [base_path, "read1.fq.gz"].iter().collect(),
            read_two: if pt.contains_r2() { Some(PathBuf::from(format!("{}{}", base_path, "/read2.fq.gz"))) } else { None },
            index_one: if pt.contains_i1() { Some(PathBuf::from(format!("{}{}", base_path, "/index1.fq.gz"))) } else { None },
            index_two: if pt.contains_i2() { Some(PathBuf::from(format!("{}{}", base_path, "/index2.fq.gz"))) } else { None },
        };

        OutputReadSetWriter::from_read_file_container(&rfc)
    }

    pub fn temp(pt: &ReadPattern, run_specs: &mut RunSpecifications) -> OutputReadSetWriter {
        //println!("r2 {} i1 {} i2 {}",pt.contains_r2(),pt.contains_i1(),pt.contains_i2());
        let rfc = ReadFileContainer {
            read_one: run_specs.create_temp_file(),
            read_two: if pt.contains_r2() { Some(run_specs.create_temp_file()) } else { None },
            index_one: if pt.contains_i1() { Some(run_specs.create_temp_file()) } else { None },
            index_two: if pt.contains_i2() { Some(run_specs.create_temp_file()) } else { None },
        };

        OutputReadSetWriter::from_read_file_container(&rfc)
    }

    pub fn create_writer(filename: &PathBuf) -> GzEncoder<File> {
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
    pub fn files(&mut self) -> ReadFileContainer {
        self.files.clone()
    }

    pub fn print_read_count(&self) {
        println!("Read 1 {}, read 2 {} read 3 {} read 4 {}", self.written_read1, self.written_read2, self.written_read3, self.written_read4);
    }

    pub fn create_x_bins(rd: &ReadIterator, x_bins: usize, run_specs: &mut RunSpecifications) -> Vec<(usize, ReadFileContainer)> {
        (0..x_bins).map(|id|
            (id, ReadFileContainer::temporary(rd, run_specs))).collect::<Vec<(usize, ReadFileContainer)>>()
    }
}

impl Clone for OutputReadSetWriter {
    fn clone(&self) -> OutputReadSetWriter {
        OutputReadSetWriter::from_read_file_container(&self.files)
    }
}

unsafe impl Send for ReadIterator {}

unsafe impl Sync for ReadIterator {}

pub struct ReadIterator {
    read_one: Option<Records<BufReader<Reader>>>,
    read_two: Option<Records<BufReader<Reader>>>,
    index_one: Option<Records<BufReader<Reader>>>,
    index_two: Option<Records<BufReader<Reader>>>,

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
                    println!("error parsing header: {:?}", e);
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
        ReadIterator::new(read_files.read_one.clone(), read_files.read_two.clone(), read_files.index_one.clone(), read_files.index_two.clone())
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

    pub fn new_reset(&self) -> Self {
        ReadIterator::new(self.path_one.clone(), self.path_two.clone(), self.path_i1.clone(), self.path_i2.clone())
    }
}

/// A light wrapper around a read set iterator with some underlying details of it's layout
//#[derive(Serialize, Deserialize)]
pub struct ClusteredReads {
    pub reads: Box<dyn Iterator<Item=ReadSetContainer>>,
    pub pattern: ReadPattern,
    pub known_size: Option<i64>,
}

impl fmt::Debug for ClusteredReads {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ClusteredReads")
            .field("pattern", &self.pattern)
            .finish()
    }
}

impl Iterator for ClusteredReads {
    type Item = ReadSetContainer;

    fn next(&mut self) -> Option<ReadSetContainer> {
        self.reads.next()
    }
}

unsafe impl Send for ClusteredReads {}

unsafe impl Sync for ClusteredReads {}

impl ClusteredReads {
    pub fn new(reads: Box<dyn Iterator<Item=ReadSetContainer>>, pattern: ReadPattern) -> ClusteredReads {
        ClusteredReads { reads, pattern, known_size: None }
    }
    pub fn new_sized(reads: Box<dyn Iterator<Item=ReadSetContainer>>, pattern: ReadPattern, size: i64) -> ClusteredReads {
        ClusteredReads { reads, pattern, known_size: Some(size) }
    }

    pub fn write_records(self, output: &mut GzEncoder<File>) {
        let length: i64 = match self.known_size {
            None => { -1 }
            Some(x) => { x }
        };
        ClusteredReads::to_disk(output, &self.pattern, length, self.reads.into_iter());
    }

    pub fn try_to_read_line_to_string(reader: &mut BufReader<GzDecoder<File>>) -> Option<String> {
        let mut header = String::new();
        let len = reader.read_line(&mut header);
        match len {
            Ok(0) => {
                None
            }
            Ok(_) => {
                //println!("Header {}",&header);
                header.pop();
                Some(header)
            }
            Err(_) => { None }
        }
    }

    pub fn from_disk(reader: &mut BufReader<GzDecoder<File>>) -> Option<ClusteredReads> {
        let mut line = String::new();
        let len = reader.read_line(&mut line).unwrap();
        println!("line ---aaaa--{}--aaaa--", &line);
        line.pop();
        let pt_read = ReadPattern::from_str(line.as_str());
        match pt_read {
            Ok(pattern) => {
                let mut line = String::new();
                let len = reader.read_line(&mut line).unwrap();
                line.pop();
                let read_count = i64::from_str(line.as_str()).unwrap();
                let mut return_vec = Vec::new();

                if read_count >= 0 {
                    for x in 0..read_count {
                        let reads_per_pattern = pattern.read_count();
                        let mut collected_reads = Vec::new();

                        for i in 0..reads_per_pattern {
                            let read = ClusteredReads::read_record(reader);
                            if let Some(x) = read { collected_reads.push(x) }
                        }
                        let converted_read = pattern.to_read_collection(&mut VecDeque::from(collected_reads));
                        match converted_read {
                            None => { panic!("Unable to read from underlying sequence stream") }
                            Some(x) => { return_vec.push(x); }
                        }
                    }
                    let ln = return_vec.len();
                    Some(ClusteredReads { reads: Box::new(return_vec.into_iter()), pattern, known_size: Some(ln as i64) })
                } else {
                    let mut cnn = 0;
                    loop {
                        cnn += 1;
                        let reads_per_pattern = pattern.read_count();
                        let mut collected_reads = Vec::new();

                        for i in 0..reads_per_pattern {
                            let read = ClusteredReads::read_record(reader);
                            if let Some(x) = read { collected_reads.push(x) }
                        }
                        let converted_read = pattern.to_read_collection(&mut VecDeque::from(collected_reads));
                        match converted_read {
                            None => { break; }
                            Some(x) => { return_vec.push(x); }
                        }
                    }
                    println!("Yup! {}",cnn);
                    let ln = return_vec.len();
                    Some(ClusteredReads { reads: Box::new(return_vec.into_iter()), pattern, known_size: Some(ln as i64) })
                }
            }
            Err(_) => {
                println!("Errored out reading!");
                panic!("errrrrr");
                None
            }
        }
    }
    pub fn read_record(reader: &mut BufReader<GzDecoder<File>>) -> Option<Record> {
        //println!("reading record...");
        let line1 = ClusteredReads::try_to_read_line_to_string(reader);
        let line2 = ClusteredReads::try_to_read_line_to_string(reader);
        let line3 = ClusteredReads::try_to_read_line_to_string(reader);
        let line4 = ClusteredReads::try_to_read_line_to_string(reader);
        match (line1, line2, line3, line4) {
            (Some(l1), Some(l2), Some(l3), Some(l4)) => {
                //let mut header = l1[1..].trim_end().splitn(2, ' ');
                //let hd = header.next().unwrap_or_default().to_owned();

                Some(Record::with_attrs(l1.as_str(), None, l2.as_bytes(), l4.as_bytes()))
            }
            (_, _, _, _) => None
        }
    }

    pub fn write_header(output: &mut GzEncoder<File>, pattern: &ReadPattern, length: i64) {
        write!(output, "{}\n", pattern);
        write!(output, "{}\n", length);
    }

    pub fn to_disk(output: &mut GzEncoder<File>, pattern: &ReadPattern, length: i64, reads: impl Iterator<Item=ReadSetContainer>) {
        ClusteredReads::write_header(output, pattern, length);

        for rl in reads {
            write!(output, "{}", rl.read_one);
            if let Some(rd) = &rl.read_two {
                write!(output, "{}", rd);
            };
            if let Some(rd) = &rl.index_one {
                write!(output, "{}", rd);
            };
            if let Some(rd) = &rl.index_two {
                write!(output, "{}", rd);
            };
        }
        output.flush();
    }

}


#[derive(Debug)]
pub struct SuperCluster {
    pub clusters: Vec<ClusteredReads>
}

impl SuperCluster {
    pub fn from_clusters(clusters: Vec<ClusteredReads>) -> SuperCluster {
        SuperCluster { clusters }
    }
}

unsafe impl Send for SuperCluster {}

unsafe impl Sync for SuperCluster {}

impl Iterator for SuperCluster {
    type Item = ClusteredReads;

    fn next(&mut self) -> Option<ClusteredReads> {
        self.clusters.pop()
    }
}

unsafe impl Send for SuperClusterOnDiskIterator {}

unsafe impl Sync for SuperClusterOnDiskIterator {}

pub struct SuperClusterOnDiskIterator {
    cluster_counts: VecDeque<i64>,
    read_files: VecDeque<PathBuf>,
    current_cluster_count: Option<i64>,
    current_reader: Option<BufReader<GzDecoder<File>>>,
    override_clusters: Option<VecDeque<ClusteredReads>>,
    pattern: ReadPattern,
}

impl Iterator for SuperClusterOnDiskIterator {
    type Item = ClusteredReads;

    fn next(&mut self) -> Option<ClusteredReads> {
        match &mut self.override_clusters {
            None => {
                let ret = if self.current_reader.is_some() &&
                    self.current_cluster_count.is_some() &&
                    (self.current_cluster_count.unwrap() > 0 ||
                        self.current_cluster_count.unwrap() < 0) {
                    println!("555 : No current cluster, trying to reload");

                    ClusteredReads::from_disk(&mut self.current_reader.as_mut().unwrap())

                } else {
                    println!("111 : No cluster");
                    None
                };

                match self.current_cluster_count {
                    Some(x) if x <= 0 => {
                        match self.read_files.pop_front() {
                            None => {
                                // we're done
                                println!("444 : No cluster");
                                self.current_cluster_count = None;
                                self.current_reader = None;
                            }
                            Some(x) => {
                                self.current_cluster_count = self.cluster_counts.pop_front();
                                println!("opening the file {:?}", &x.as_path());
                                self.current_reader = Some(BufReader::new(GzDecoder::new(File::open(&x.as_path()).unwrap())));
                            }
                        }
                    }
                    _ => {}
                }
                ret
            }
            Some(x) => {
                println!("299 : No cluster");
                x.pop_front()
            }
        }
    }
}

impl SuperClusterOnDiskIterator {
    fn to_disk(clusters: Vec<ClusteredReads>, read_pattern: ReadPattern, run_specs: &mut RunSpecifications) -> SuperClusterOnDiskIterator {
        let temp_file = run_specs.create_temp_file();
        let mut temp_file_writer: GzEncoder<File> = OutputReadSetWriter::create_writer(&temp_file);
        let cluster_count = clusters.len();
        for cluster in clusters {
            ClusteredReads::to_disk(&mut temp_file_writer, &read_pattern, cluster_count as i64, cluster.reads);
        }

        let input = BufReader::new(GzDecoder::new(File::open(&temp_file.as_path()).unwrap()));
        let files: Vec<PathBuf> = Vec::new();
        let counts: Vec<i64> = Vec::new();
        SuperClusterOnDiskIterator {
            cluster_counts: VecDeque::from(counts),
            read_files: VecDeque::from(files),
            current_cluster_count: Some(cluster_count as i64),
            current_reader: Some(BufReader::new(GzDecoder::new(File::open(temp_file.as_path()).unwrap()))),
            override_clusters: None,
            pattern: read_pattern,
        }
    }


    pub fn new_from_vec(clusters: Vec<ClusteredReads>, read_pattern: ReadPattern, run_specs: &mut RunSpecifications) -> SuperClusterOnDiskIterator {
        SuperClusterOnDiskIterator::to_disk(clusters, read_pattern, run_specs)
    }

    pub fn new_from_flat_iterator(iter: ReadIterator, read_pattern: ReadPattern, run_specs: &mut RunSpecifications) -> SuperClusterOnDiskIterator {
        let cr = vec![ClusteredReads {
            reads: Box::new(iter),
            pattern: read_pattern.clone(),
            known_size: None,
        }];

        SuperClusterOnDiskIterator{
            cluster_counts: VecDeque::new(),
            read_files: VecDeque::new(),
            current_cluster_count: None,
            current_reader: None,
            override_clusters: Some(VecDeque::from(cr)),
            pattern: read_pattern.clone(),
        }
    }

    pub fn new_from_read_file_container(writer_files: Vec<PathBuf>, read_pattern: ReadPattern, run_specs: RunSpecifications) -> SuperClusterOnDiskIterator {
        let ln = writer_files.len();
        let mut queued_files = VecDeque::from(writer_files);
        let next = queued_files.pop_front().unwrap();
        println!("First file opening {:?}",&next);

        SuperClusterOnDiskIterator {
            read_files: queued_files,
            cluster_counts: VecDeque::from(vec![-1; ln]),
            current_cluster_count: Some(-1),
            current_reader: Some(BufReader::new(GzDecoder::new(File::open(next.as_path()).unwrap()))),
            override_clusters: None,
            pattern: read_pattern,
        }
    }
}


