use crate::read_strategies::sequence_file_containers::OutputReadSetWriter;
use std::collections::HashMap;
use crate::RunSpecifications;
use tempfile::NamedTempFile;

pub struct StaticCorrectionBalancer {
    files: Vec<NamedTempFile>,
    writers: Vec<OutputReadSetWriter>,
    assigned_sequences: HashMap<Vec<u8>,usize>,
    output_counts: HashMap<usize,usize>,
}
/*
impl StaticCorrectionBalancer {
pub fn from(run_specs: &RunSpecifications) {
    let mut files = Vec::new();
    let mut assigned_sequences : HashMap<Vec<u8>,usize> = HashMap::new();
    let mut output_counts : HashMap<usize,usize> = HashMap::new();

    for i in 0..run_specs.sorting_file_count {
        tempfile
        files.push();
        output_counts.insert(i,0);
    }
    StaticCorrectionBalancer{
        files: vec![],
        assigned_sequences,
        output_counts
    }
}

}*/