use std::collections::{HashMap, VecDeque};
use std::path::PathBuf;

use crate::read_strategies::sequence_file_containers::{ReadPattern, ReadSetContainer};
use crate::read_strategies::sequence_file_containers::OutputReadSetWriter;
use crate::RunSpecifications;
use indicatif::style::ProgressTracker;
use std::borrow::{BorrowMut, Borrow};
use crate::read_strategies::sequence_file_containers::ReadFileContainer;
use crate::read_strategies::sequence_file_containers::SuperClusterOnDiskIterator;
use crate::read_strategies::sequence_file_containers::ReadIterator;
use crate::read_strategies::sequence_file_containers::ClusteredReads;
use flate2::write::GzEncoder;
use std::fs::File;
use std::io::Write;


pub struct RoundRobinDiskWriter {
    writers: Vec<Option<GzEncoder<File>>>,
    underlying_files: Vec<PathBuf>,
    read_buffer_size: usize,
    buffered_reads: Vec<VecDeque<ReadSetContainer>>,
    assigned_sequences: HashMap<Vec<u8>, usize>,
    assigned_sequences_count: HashMap<Vec<u8>, usize>,
    output_counts: HashMap<usize, usize>,
    current_bin: usize,
    read_pattern: ReadPattern,
    run_specs: RunSpecifications,
    to_disk_count: usize,
}

impl RoundRobinDiskWriter {
    pub fn from(run_specs: &mut RunSpecifications, read_pattern: &ReadPattern, read_buffer_size: usize) -> RoundRobinDiskWriter {
        let mut writers = Vec::new();
        let mut buffered_reads: Vec<VecDeque<ReadSetContainer>> = Vec::new();
        let mut underlying_files = Vec::new();
        let mut assigned_sequences: HashMap<Vec<u8>, usize> = HashMap::new();
        let mut assigned_sequences_count: HashMap<Vec<u8>, usize> = HashMap::new();

        let mut output_counts: HashMap<usize, usize> = HashMap::new();

        for i in 0..run_specs.sorting_file_count {
            let temp_file = run_specs.create_temp_file();
            underlying_files.push(temp_file.clone());
            buffered_reads.push(VecDeque::with_capacity(read_buffer_size));

            writers.push(None);
            output_counts.insert(i, 0);
        };
        println!("len : {}",buffered_reads.len());

        RoundRobinDiskWriter {
            writers,
            underlying_files,
            read_buffer_size,
            buffered_reads,
            assigned_sequences,
            assigned_sequences_count,
            output_counts,
            current_bin: 0,
            read_pattern: read_pattern.clone(),
            run_specs: run_specs.clone(),
            to_disk_count: 0,
        }
    }

    pub fn write_all_reads(writer: &mut GzEncoder<File>, read: &ReadSetContainer) {
        write!(writer, "{}", read.read_one.to_string());
        if let Some(x) = &read.read_two { write!(writer, "{}", x.to_string()); }
        if let Some(x) = &read.index_one { write!(writer, "{}", x.to_string()); }
        if let Some(x) = &read.index_two { write!(writer, "{}", x.to_string()); }
    }

    pub fn dump_buffer(&mut self, bin: usize) {
        let reads = &self.buffered_reads[bin];

        let mut writer: GzEncoder<File> = OutputReadSetWriter::create_writer(&self.underlying_files.get(bin).unwrap());

        ClusteredReads::write_header(&mut writer, &self.read_pattern, -1);

        for read in reads.into_iter() {
            self.to_disk_count += 1;
            RoundRobinDiskWriter::write_all_reads(&mut writer, &read);
        }

        self.writers.insert(bin,Some(writer));
    }

    pub fn store_or_open(&mut self, bin: usize, read: &ReadSetContainer) {
        match self.writers.get(bin).unwrap().is_some() {
            true => {
                let mut writer = self.writers.get_mut(bin).unwrap().as_mut().unwrap();
                self.to_disk_count += 1;
                RoundRobinDiskWriter::write_all_reads(writer, read);
            }
            false => {
                self.buffered_reads[bin].push_back(read.clone());
                if self.buffered_reads[bin].len() > self.read_buffer_size {
                    self.dump_buffer(bin);
                }
            }
        }
    }

    pub fn assign_new_bin(&mut self, sort_string: &Vec<u8>) {
        assert!(self.assigned_sequences.get(sort_string).is_none());
        self.assigned_sequences.insert(sort_string.clone(), self.current_bin);

        self.current_bin += 1;
        if self.current_bin >= self.output_counts.len() {
            self.current_bin = 0;
        };
    }

    pub fn write(&mut self, sort_string: &Vec<u8>, read: &ReadSetContainer) {
        if self.assigned_sequences.get(sort_string).is_none() {
            self.assign_new_bin(sort_string);
        }

        self.assigned_sequences_count.insert(
            sort_string.clone(),
            self.assigned_sequences_count.get(sort_string).unwrap_or(&(0 as usize)) + 1,
        );

        self.store_or_open(*self.assigned_sequences.get(sort_string).unwrap(),read);
    }

    /// Destructively get a link to the underlying writers
    pub fn close_and_recover_iterator(self) -> SuperClusterOnDiskIterator {
        let full_output_file = self.run_specs.create_temp_file().clone();
        let mut output_file = OutputReadSetWriter::create_writer(&full_output_file);
        let mut written_reads = 0;
        for i in (0..self.buffered_reads.len()) {
            println!("testing bin {} {} {}",i,self.writers.get(i).unwrap().is_none(),self.buffered_reads.get(i).unwrap().len());
            if self.writers.get(i).unwrap().is_none() && self.buffered_reads.get(i).unwrap().len() > 0 {
                ClusteredReads::write_header(&mut output_file,
                                             &self.read_pattern,
                                             self.buffered_reads.get(i).unwrap().len() as i64);

                written_reads += self.buffered_reads.get(i).unwrap().len();

                for read in self.buffered_reads.get(i).unwrap().into_iter() {
                    RoundRobinDiskWriter::write_all_reads(&mut output_file, &read);
                }
            }
        }

        for brb in self.buffered_reads {
            println!("bin size = {}", brb.len());
        }

        println!("written {} post {}",self.to_disk_count, written_reads);
        let mut written_files : Vec<PathBuf> = Vec::new();

        if written_reads > 0 {
            written_files.push(full_output_file);
        }

        for i in 0..self.underlying_files.len() {
            if self.writers.get(i).unwrap().is_some() {
                println!("adding");
                written_files.push(self.underlying_files.get(i).unwrap().clone());
            }
            println!("bin = {} size = {}, file {:?}", i, self.output_counts.get(&i).unwrap(), self.underlying_files.get(i).unwrap());
        }
        println!("lookup set size = {}",self.assigned_sequences.len());

        SuperClusterOnDiskIterator::new_from_read_file_container(
            written_files,
            self.read_pattern.clone(),
            self.run_specs.clone())
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::str;
    use std::sync::Arc;
    use bio::io::fastq::Record;
    use rand::{seq::IteratorRandom, thread_rng};
    use itertools::Itertools;
    use crate::TempDir;

    pub fn random_sequence(length: usize) -> String {
        let bases = vec![b'A', b'C', b'G', b'T'];
        let mut rng = thread_rng();

        String::from_utf8(bases.iter().choose_multiple(&mut rng, length).iter().map(|c| **c).collect::<Vec<u8>>()).unwrap()
    }

    pub fn all_combinations(n: usize) -> Vec<String> {
        let characters = vec!["A", "C", "G", "T"];// vec![vec!["A"],vec!["C"],vec!["G"],vec!["T"]];
        let mut combs: Vec<String> = Vec::new();

        (2..n).fold(
            characters.iter().map(|c| characters.iter().map(move |&d| d.to_owned() + *c)).flatten().collect(),
            |acc, _| acc.into_iter().map(|c| characters.iter().map(move |&d| d.to_owned() + &*c)).flatten().collect(),
        )
    }

    pub fn create_fake_quality_scores(length: usize) -> Vec<u8> {
        vec![b'H'; length]
    }

    pub fn fake_reads(full_length: usize, permutation_leader_size: usize) -> Vec<ReadSetContainer> {
        //sort_string: &Vec<u8>, read: &ReadSetContainer
        let mut fake_reads = Vec::new();
        let all_perm = all_combinations(permutation_leader_size);
        println!("All perm size {} for size {}", all_perm.len(), permutation_leader_size);
        for sequence_leader in all_perm {
            let mut read_seq = sequence_leader.clone();
            read_seq.push_str(random_sequence(full_length - permutation_leader_size).as_str());
            let qual = create_fake_quality_scores(full_length);
            let record = bio::io::fastq::Record::with_attrs("fakeRead", None, read_seq.as_bytes(), qual.as_slice());
            let fake_rsc = ReadSetContainer::new_from_read1(record);
            fake_reads.push(fake_rsc);
        }
        fake_reads
    }

    #[test]
    fn fill_balancer_with_reads() {
        let mut run_specs = RunSpecifications {
            estimated_reads: 0,
            sorting_file_count: 20,
            sorting_threads: 0,
            processing_threads: 0,
            tmp_location: Arc::new(TempDir::new().unwrap()),
        };

        let mut balancer = RoundRobinDiskWriter::from(&mut run_specs,
                                                      &ReadPattern::ONE,
                                                      100);

        let fake_reads = fake_reads(50, 6);
        let fake_read_len = fake_reads.len();
        println!("adding {} reads", fake_read_len);

        for read in fake_reads {
            let sort_string = read.read_one.seq()[0..6].to_vec();
            balancer.write(&sort_string, &read);
        }

        let recovered = balancer.close_and_recover_iterator();
        let mut cnt = 0;
        for r in recovered {
            for read in r {
                cnt += 1;
            }
        }
        assert_eq!(cnt, fake_read_len);
    }

    #[test]
    fn light_balancer_with_reads() {
        let mut run_specs = RunSpecifications {
            estimated_reads: 0,
            sorting_file_count: 20,
            sorting_threads: 0,
            processing_threads: 0,
            tmp_location: Arc::new(TempDir::new().unwrap()),
        };

        let mut balancer = RoundRobinDiskWriter::from(&mut run_specs,
                                                      &ReadPattern::ONE,
                                                      100);

        let fake_reads = fake_reads(50, 3);
        let fake_read_len = fake_reads.len();
        println!("adding {} reads", fake_read_len);

        for read in fake_reads {
            let sort_string = read.read_one.seq()[0..6].to_vec();
            balancer.write(&sort_string, &read);
        }

        let recovered = balancer.close_and_recover_iterator();
        let mut cnt = 0;
        for r in recovered {
            for read in r {
                cnt += 1;
            }
        }
        assert_eq!(cnt, fake_read_len);
    }
}