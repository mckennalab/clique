use std::collections::{HashMap, VecDeque};
use std::path::PathBuf;

use crate::read_strategies::sequence_file_containers::{ReadPattern, ReadSetContainer};
use crate::read_strategies::sequence_file_containers::OutputReadSetWriter;
use crate::RunSpecifications;
use indicatif::style::ProgressTracker;
use std::borrow::{BorrowMut, Borrow};
use crate::read_strategies::sequence_file_containers::ReadFileContainer;

pub struct RoundRobinDiskWriter {
    writers: Vec<OutputReadSetWriter>,
    assigned_sequences: HashMap<Vec<u8>, usize>,
    output_counts: HashMap<usize, usize>,
    current_bin: usize,
}

impl RoundRobinDiskWriter {
    pub fn from(run_specs: &mut RunSpecifications, read_pattern: &ReadPattern) -> RoundRobinDiskWriter {
        let mut writers = Vec::new();
        let mut assigned_sequences: HashMap<Vec<u8>, usize> = HashMap::new();
        let mut output_counts: HashMap<usize, usize> = HashMap::new();

        for i in 0..run_specs.sorting_file_count {
            let writer = OutputReadSetWriter::temp(read_pattern, run_specs);
            writers.push(writer);
            output_counts.insert(i, 0);
        }
        RoundRobinDiskWriter {
            writers,
            assigned_sequences,
            output_counts,
            current_bin: 0
        }
    }

    pub fn write(&mut self, sort_string: &Vec<u8>, read: &ReadSetContainer) {
        match self.assigned_sequences.get(sort_string) {
            Some(x) => {
                let mut writer = self.writers.get_mut(*x).unwrap();
                writer.write(read);
            },
            _ => {
                self.assigned_sequences.insert(sort_string.clone(),self.current_bin);
                self.output_counts.insert(self.current_bin,
                                          self.output_counts.get(&self.current_bin).unwrap_or(&(0 as usize)) + 1);
                self.writers.get_mut(self.current_bin).unwrap().write(read);
                self.current_bin += 1;
                if self.current_bin >= self.output_counts.len() {
                    self.current_bin = 0;
                }
            }
        }
    }

    pub fn get_writers(&mut self) -> VecDeque<ReadFileContainer> {
        let mut writer_files = Vec::new();
        for writer in &mut self.writers {
            writer_files.push(writer.files().clone());
        }
        VecDeque::from(writer_files)
    }
}