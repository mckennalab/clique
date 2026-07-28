//! Reading raw sequencing input: [`ReadSetContainer`] bundles the one-to-four
//! FASTQ records of a spot, and [`ReadIterator`] streams them from the input
//! files.

use bio::io::fastq::{Record, Records};
use bio::io::fastq::Reader as FqReader;
use serde::{Serialize, Deserialize};
use std::io::{BufReader};
use std::path::PathBuf;
use rust_htslib::bgzf::Reader;

/// Return the identifier portion of a FASTQ header, excluding optional
/// whitespace-delimited metadata that is not valid in a SAM/BAM QNAME.
pub fn fastq_record_id(id: &str) -> &str {
    id.split_ascii_whitespace().next().unwrap_or(id)
}

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
    read_one: Records<BufReader<Reader>>,
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
               index_2: Option<PathBuf>,
    ) -> ReadIterator {
        let r_one = ReadIterator::open_reader(&read_1, "read 1");
        let read2 = ReadIterator::open_optional_reader(read_2, "read 2");
        let index1 = ReadIterator::open_optional_reader(index_1, "index 1");
        let index2 = ReadIterator::open_optional_reader(index_2, "index 2");

        ReadIterator {
            read_one: r_one,
            read_two: read2,
            index_one: index1,
            index_two: index2,
            reads_processed: 0,
            broken_reads: 0,
        }
    }

    fn open_reader(path: &PathBuf, label: &str) -> Records<BufReader<Reader>> {
        if !path.exists() {
            panic!("Unable to open {} FASTQ file: {}", label, path.display());
        }

        info!("Opening {} file: {}", label, path.display());
        let reader = Reader::from_path(path).unwrap_or_else(|error| {
            panic!(
                "Unable to open {} FASTQ file {}: {:?}",
                label,
                path.display(),
                error
            )
        });
        FqReader::new(reader).records()
    }

    fn open_optional_reader(
        path: Option<PathBuf>,
        label: &str,
    ) -> Option<Records<BufReader<Reader>>> {
        match path {
            None => None,
            Some(path) if path == PathBuf::from("NONE") => None,
            Some(path) => Some(ReadIterator::open_reader(&path, label)),
        }
    }

    fn canonical_read_id(id: &str) -> &str {
        let id = fastq_record_id(id);
        id.strip_suffix("/1")
            .or_else(|| id.strip_suffix("/2"))
            .or_else(|| id.strip_suffix("/3"))
            .or_else(|| id.strip_suffix("/4"))
            .unwrap_or(id)
    }

    fn validate_record(record: &Record, label: &str, record_number: usize) {
        if let Err(error) = record.check() {
            panic!(
                "Invalid {} FASTQ record {}: {}",
                label, record_number, error
            );
        }
    }

    fn next_companion(
        reader: &mut Option<Records<BufReader<Reader>>>,
        label: &str,
        expected_id: &str,
        record_number: usize,
    ) -> Option<Record> {
        let records = match reader.as_mut() {
            None => return None,
            Some(records) => records,
        };
        let record = match records.next() {
            Some(Ok(record)) => record,
            Some(Err(error)) => panic!(
                "Unable to parse {} FASTQ record {}: {:?}",
                label, record_number, error
            ),
            None => panic!(
                "{} FASTQ ended before read 1 at record {}",
                label, record_number
            ),
        };
        ReadIterator::validate_record(&record, label, record_number);

        let actual_id = ReadIterator::canonical_read_id(record.id());
        if actual_id != expected_id {
            panic!(
                "FASTQ identifiers differ at record {}: read 1 is '{}', {} is '{}'",
                record_number,
                expected_id,
                label,
                actual_id
            );
        }

        Some(record)
    }

    fn assert_companion_exhausted(
        reader: &mut Option<Records<BufReader<Reader>>>,
        label: &str,
        reads_processed: usize,
    ) {
        if let Some(records) = reader.as_mut() {
            match records.next() {
                None => {}
                Some(Ok(record)) => panic!(
                    "{} FASTQ contains extra record '{}' after read 1 ended at {} records",
                    label,
                    record.id(),
                    reads_processed
                ),
                Some(Err(error)) => panic!(
                    "Unable to parse {} FASTQ after read 1 ended at {} records: {:?}",
                    label, reads_processed, error
                ),
            }
        }
    }

}
impl Iterator for ReadIterator {
    type Item = ReadSetContainer;

    fn next(&mut self) -> Option<ReadSetContainer> {
        match self.read_one.next() {
            Some(Ok(read_one)) => {
                let record_number = self.reads_processed + 1;
                ReadIterator::validate_record(&read_one, "read 1", record_number);
                let expected_id = ReadIterator::canonical_read_id(read_one.id()).to_string();
                let read_two = ReadIterator::next_companion(
                    &mut self.read_two,
                    "read 2",
                    &expected_id,
                    record_number,
                );
                let index_one = ReadIterator::next_companion(
                    &mut self.index_one,
                    "index 1",
                    &expected_id,
                    record_number,
                );
                let index_two = ReadIterator::next_companion(
                    &mut self.index_two,
                    "index 2",
                    &expected_id,
                    record_number,
                );
                self.reads_processed += 1;

                Some(ReadSetContainer {
                    read_one,
                    read_two,
                    index_one,
                    index_two,
                })
            }
            Some(Err(error)) => {
                self.broken_reads += 1;
                panic!(
                    "Unable to parse read 1 FASTQ record {}: {:?}",
                    self.reads_processed + 1,
                    error
                );
            }
            None => {
                ReadIterator::assert_companion_exhausted(
                    &mut self.read_two,
                    "read 2",
                    self.reads_processed,
                );
                ReadIterator::assert_companion_exhausted(
                    &mut self.index_one,
                    "index 1",
                    self.reads_processed,
                );
                ReadIterator::assert_companion_exhausted(
                    &mut self.index_two,
                    "index 2",
                    self.reads_processed,
                );
                info!("Done processing {} reads", self.reads_processed);
                None
            }
        }
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::io::fastq::Record;
    use serde_yaml;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_record(id: &str, seq: &[u8], qual: &[u8]) -> Record {
        Record::with_attrs(id, None, seq, qual)
    }

    fn write_fastq(contents: &str) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(contents.as_bytes()).unwrap();
        file.flush().unwrap();
        file
    }

    #[test]
    fn test_read_iterator_validates_paired_records() {
        let read_one = write_fastq("@spot/1\nACGT\n+\nHHHH\n");
        let read_two = write_fastq("@spot/2\nTGCA\n+\nIIII\n");
        let mut iterator = ReadIterator::new(
            read_one.path().to_path_buf(),
            Some(read_two.path().to_path_buf()),
            None,
            None,
        );

        let record = iterator.next().unwrap();
        assert_eq!(record.read_one.id(), "spot/1");
        assert_eq!(record.read_two.unwrap().id(), "spot/2");
        assert!(iterator.next().is_none());
        assert_eq!(iterator.reads_processed, 1);
    }

    #[test]
    fn test_fastq_record_id_removes_nanopore_metadata() {
        assert_eq!(
            fastq_record_id("12d75f4f-f926-49fa-afdc-116c3b382e16\tqs:f:23.5\tch:i:2366"),
            "12d75f4f-f926-49fa-afdc-116c3b382e16"
        );
        assert_eq!(fastq_record_id("spot/1 comment"), "spot/1");
        assert_eq!(fastq_record_id("plain_id"), "plain_id");
    }

    #[test]
    fn test_read_iterator_accepts_none_sentinel() {
        let read_one = write_fastq("@spot\nACGT\n+\nHHHH\n");
        let iterator = ReadIterator::new(
            read_one.path().to_path_buf(),
            Some(PathBuf::from("NONE")),
            Some(PathBuf::from("NONE")),
            Some(PathBuf::from("NONE")),
        );

        assert_eq!(iterator.count(), 1);
    }

    #[test]
    #[should_panic(expected = "FASTQ identifiers differ")]
    fn test_read_iterator_rejects_mismatched_ids() {
        let read_one = write_fastq("@spot_a\nACGT\n+\nHHHH\n");
        let read_two = write_fastq("@spot_b\nTGCA\n+\nIIII\n");
        let mut iterator = ReadIterator::new(
            read_one.path().to_path_buf(),
            Some(read_two.path().to_path_buf()),
            None,
            None,
        );

        iterator.next();
    }

    #[test]
    #[should_panic(expected = "ended before read 1")]
    fn test_read_iterator_rejects_shorter_companion() {
        let read_one = write_fastq(
            "@spot_1\nACGT\n+\nHHHH\n@spot_2\nACGT\n+\nHHHH\n",
        );
        let read_two = write_fastq("@spot_1\nTGCA\n+\nIIII\n");
        let mut iterator = ReadIterator::new(
            read_one.path().to_path_buf(),
            Some(read_two.path().to_path_buf()),
            None,
            None,
        );

        assert!(iterator.next().is_some());
        iterator.next();
    }

    #[test]
    #[should_panic(expected = "contains extra record")]
    fn test_read_iterator_rejects_longer_companion() {
        let read_one = write_fastq("@spot_1\nACGT\n+\nHHHH\n");
        let read_two = write_fastq(
            "@spot_1\nTGCA\n+\nIIII\n@spot_2\nTGCA\n+\nIIII\n",
        );
        let mut iterator = ReadIterator::new(
            read_one.path().to_path_buf(),
            Some(read_two.path().to_path_buf()),
            None,
            None,
        );

        assert!(iterator.next().is_some());
        iterator.next();
    }

    #[test]
    #[should_panic(expected = "Invalid read 1 FASTQ record")]
    fn test_read_iterator_rejects_malformed_read_one() {
        let read_one = write_fastq("@spot\nACGT\n+\nHHH\n");
        let mut iterator = ReadIterator::new(read_one.path().to_path_buf(), None, None, None);

        iterator.next();
    }

    #[test]
    #[should_panic(expected = "Invalid read 2 FASTQ record")]
    fn test_read_iterator_rejects_malformed_companion() {
        let read_one = write_fastq("@spot\nACGT\n+\nHHHH\n");
        let read_two = write_fastq("@spot\nTGCA\n+\nIII\n");
        let mut iterator = ReadIterator::new(
            read_one.path().to_path_buf(),
            Some(read_two.path().to_path_buf()),
            None,
            None,
        );

        iterator.next();
    }

    #[test]
    #[should_panic(expected = "Unable to open read 2 FASTQ file")]
    fn test_read_iterator_rejects_missing_companion_path() {
        let read_one = write_fastq("@spot\nACGT\n+\nHHHH\n");
        let missing_read_two = NamedTempFile::new().unwrap();
        let missing_read_two_path = missing_read_two.path().to_path_buf();
        drop(missing_read_two);

        ReadIterator::new(
            read_one.path().to_path_buf(),
            Some(missing_read_two_path),
            None,
            None,
        );
    }

    #[test]
    fn test_read_set_container_new_from_read1() {
        let rec = make_record("read1", b"ACGT", b"HHHH");
        let rsc = ReadSetContainer::new_from_read1(rec.clone());
        assert_eq!(rsc.read_one.id(), "read1");
        assert!(rsc.read_two.is_none());
        assert!(rsc.index_one.is_none());
        assert!(rsc.index_two.is_none());
    }

    #[test]
    fn test_read_set_container_clone_read_only() {
        let rec = make_record("read1", b"ACGT", b"HHHH");
        let rsc = ReadSetContainer::new_from_read1(rec);
        let cloned = rsc.clone();
        assert_eq!(cloned.read_one.id(), rsc.read_one.id());
        assert_eq!(cloned.read_one.seq(), rsc.read_one.seq());
        assert!(cloned.read_two.is_none());
        assert!(cloned.index_one.is_none());
        assert!(cloned.index_two.is_none());
    }

    #[test]
    fn test_read_set_container_clone_all_fields() {
        let rsc = ReadSetContainer {
            read_one: make_record("r1", b"ACGT", b"HHHH"),
            read_two: Some(make_record("r2", b"TGCA", b"IIII")),
            index_one: Some(make_record("i1", b"AA", b"HH")),
            index_two: Some(make_record("i2", b"CC", b"HH")),
        };
        let cloned = rsc.clone();
        assert_eq!(cloned.read_one.id(), "r1");
        assert_eq!(cloned.read_two.as_ref().unwrap().id(), "r2");
        assert_eq!(cloned.index_one.as_ref().unwrap().id(), "i1");
        assert_eq!(cloned.index_two.as_ref().unwrap().id(), "i2");
    }

    #[test]
    fn test_read_set_container_display() {
        let rsc = ReadSetContainer::new_from_read1(
            make_record("r1", b"ACGT", b"HHHH"),
        );
        let display = format!("{}", rsc);
        assert!(display.contains("r1"));
    }

    #[test]
    fn test_read_set_container_equality() {
        let rsc1 = ReadSetContainer::new_from_read1(make_record("r1", b"ACGT", b"HHHH"));
        let rsc2 = ReadSetContainer::new_from_read1(make_record("r1", b"ACGT", b"HHHH"));
        assert_eq!(rsc1, rsc2);
    }

    #[test]
    fn test_read_set_container_inequality() {
        let rsc1 = ReadSetContainer::new_from_read1(make_record("r1", b"ACGT", b"HHHH"));
        let rsc2 = ReadSetContainer::new_from_read1(make_record("r2", b"TGCA", b"IIII"));
        assert_ne!(rsc1, rsc2);
    }

    #[test]
    fn test_read_set_container_serialize_deserialize() {
        let rsc = ReadSetContainer {
            read_one: make_record("r1", b"ACGT", b"HHHH"),
            read_two: Some(make_record("r2", b"TGCA", b"IIII")),
            index_one: None,
            index_two: None,
        };
        let serialized = serde_yaml::to_string(&rsc).unwrap();
        let deserialized: ReadSetContainer = serde_yaml::from_str(&serialized).unwrap();
        assert_eq!(deserialized.read_one.id(), "r1");
        assert_eq!(deserialized.read_two.as_ref().unwrap().id(), "r2");
        assert!(deserialized.index_one.is_none());
    }
}
