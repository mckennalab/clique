use std::collections::{HashMap, VecDeque};
use std::slice::Iter;

use crate::consensus::consensus_builders::create_poa_consensus;
use crate::read_strategies::sequence_file_containers::{ReadIterator, ReadSetContainer};
use crate::read_strategies::sequence_file_containers::ClusteredReads;
use crate::read_strategies::sequence_file_containers::SuperClusterOnDiskIterator;
use crate::read_strategies::sequence_file_containers::ReadPattern;
use crate::read_strategies::sequence_layout::LayoutType;
use crate::read_strategies::sequence_layout::SequenceLayout;
use crate::read_strategies::sequence_layout::transform;
use crate::sorters::sorter::SortStructure;
use crate::umis::sequence_clustering::{split_subgroup, string_distance};
use crate::umis::sequence_clustering::get_connected_components;
use crate::umis::sequence_clustering::input_list_to_graph;
use crate::umis::sequence_clustering::InputList;
use crate::umis::sequence_clustering::average_dist;
use crate::RunSpecifications;
use crate::read_strategies::sequence_file_containers::{ReadFileContainer, OutputReadSetWriter};
use std::cell::RefCell;
use std::io::Read;
use std::sync::Mutex;


pub trait SortStream<'z> {
    fn from_read_iterator(read_iter: Box<dyn Iterator<Item=ReadSetContainer>>, read_pattern: &ReadPattern, sort_structure: &SortStructure, layout: &LayoutType, run_specs: &'z mut RunSpecifications) -> Self;
    fn from_read_collection(read_collection: ClusteredReads, sort_structure: &SortStructure, layout: &LayoutType, pattern: ReadPattern, run_specs: &'z mut RunSpecifications) -> Self;
    fn sorted_read_set(self) -> Option<SuperClusterOnDiskIterator>;
}

pub struct ClusteredDiskSortStream<'z> {
    reads: HashMap<Vec<u8>, Vec<ReadSetContainer>>,
    sort_structure: SortStructure,
    layout: LayoutType,
    pattern: ReadPattern,
    run_specs: &'z mut RunSpecifications,
}

impl<'z> ClusteredDiskSortStream<'z> {
    fn push(&mut self, read: ReadSetContainer) {
        let seq_struct = transform(read, &self.layout);
        let sequence = self.sort_structure.get_field(&seq_struct).unwrap();

        let elements = self.reads.entry(sequence).or_insert(vec![]);
        elements.push(seq_struct.original_reads().unwrap());
    }
}

impl<'z> SortStream<'z> for ClusteredDiskSortStream<'z> {
    fn from_read_iterator(read_iter: Box<dyn Iterator<Item=ReadSetContainer>>, read_pattern: &ReadPattern, sort_structure: &SortStructure, layout: &LayoutType, run_specs: &'z mut RunSpecifications) -> Self {
        let mut reads: HashMap<Vec<u8>, Vec<ReadSetContainer>> = HashMap::new();
        let mut mem_sort = ClusteredDiskSortStream {
            reads,
            sort_structure: sort_structure.clone(),
            layout: layout.clone(),
            pattern: read_pattern.clone(),
            run_specs: run_specs,
        };
        for read in read_iter {
            mem_sort.push(read.clone());
        }
        mem_sort
    }

    fn from_read_collection(read_collection: ClusteredReads, sort_structure: &SortStructure, layout: &LayoutType, pattern: ReadPattern, run_specs: &'z mut RunSpecifications) -> Self {
        let mut reads: HashMap<Vec<u8>, Vec<ReadSetContainer>> = HashMap::new();
        let mut mem_sort = ClusteredDiskSortStream {
            reads,
            sort_structure: sort_structure.clone(),
            layout: layout.clone(),
            pattern: pattern.clone(),
            run_specs: run_specs,
        };
        for read in read_collection {
            mem_sort.push(read.clone());
        }
        mem_sort
    }

    fn sorted_read_set(self) -> Option<SuperClusterOnDiskIterator> {
        let max_dist = match &self.sort_structure {
            SortStructure::KNOWN_LIST { layout_type, max_distance: maximum_distance, on_disk, known_list } => { maximum_distance }
            SortStructure::HD_UMI { layout_type, max_distance, on_disk } => { max_distance }
            SortStructure::LD_UMI { layout_type, max_distance, on_disk } => { max_distance }
        };

        let collection = InputList { strings: self.reads.keys().map(|x| x.clone()).collect::<Vec<Vec<u8>>>(), max_dist: *max_dist as u64 };
        let mut graph = input_list_to_graph(&collection, string_distance, false);

        let cc = get_connected_components(&graph);

        let mut final_vec: Vec<ClusteredReads> = Vec::new();

        for group in cc {
            let minilist = InputList { strings: group, max_dist: *max_dist as u64 };
            let mut minigraph = input_list_to_graph(&minilist, string_distance, false);

            let issubgroups = split_subgroup(&mut minigraph);
            match issubgroups {
                None => {
                    let mut rc = Vec::new();
                    for s in &minilist.strings {
                        if let Some(i) = self.reads.get(s) {
                            rc.extend(i.clone());
                        }
                    }
                    final_vec.push(ClusteredReads::new(Box::new(VecDeque::from(rc).into_iter()), self.pattern.clone()));
                }
                Some(x) => {
                    for sgroup in &x {
                        let average_distance = average_dist(&sgroup, string_distance);
                        let mut rc = Vec::new();
                        for s in sgroup {
                            if let Some(i) = self.reads.get(s) {
                                rc.extend(i.clone());
                            }
                        }
                        let iter  = Box::new(VecDeque::from(rc).into_iter());
                        final_vec.push(ClusteredReads::new(iter, self.pattern.clone()));
                    }
                }
            }
        }
        Some(SuperClusterOnDiskIterator::new_from_vec(final_vec, self.pattern.clone(), &mut self.run_specs.clone()))
    }
}
