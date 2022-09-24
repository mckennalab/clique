use std::collections::HashMap;
use std::path::PathBuf;

use crate::read_strategies::sequence_file_containers::{ReadPattern, ReadSetContainer};
use crate::read_strategies::sequence_file_containers::OutputReadSetWriter;
use crate::RunSpecifications;
use indicatif::style::ProgressTracker;
use std::borrow::BorrowMut;

pub struct RoundRobinDiskWriter {
    writers: Vec<OutputReadSetWriter>,
    assigned_sequences: HashMap<Vec<u8>, usize>,
    output_counts: HashMap<usize, usize>,
    current_bin: usize,
}
/*
impl RoundRobinDiskWriter {
    pub fn from(run_specs: &RunSpecifications, read_pattern: &ReadPattern) -> RoundRobinDiskWriter {
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
        }
    }

    pub fn write(&mut self, sort_string: &Vec<u8>, read: &ReadSetContainer) {
        match self.assigned_sequences.get(sort_string) {
            Some(x) => {
                self.writers.get(x).unwrap().write(read);
            },
            _ => {
                self.assigned_sequences.insert(sort_string.clone(),self.current_bin);
                self.output_counts.insert(self.current_bin,self.output_counts.get(&self.current_bin).unwrap_or(&(0 as usize)) + 1);
                self.writers.get(x).unwrap().write(read);
            }
        }
    }
}*/