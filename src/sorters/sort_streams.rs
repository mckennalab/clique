use std::collections::{HashMap, VecDeque};
use std::slice::Iter;

use crate::consensus::consensus_builders::create_poa_consensus;
use crate::read_strategies::sequence_file_containers::{ReadIterator, ReadSetContainer};
use crate::read_strategies::sequence_file_containers::NamedReadSetContainer;
use crate::read_strategies::sequence_file_containers::ClusteredReads;
use crate::read_strategies::sequence_file_containers::ClusteredReadIterator;
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

pub trait SortStream {
    fn from_read_iterator(read_iter: ReadIterator, sort_structure: &SortStructure, layout: &LayoutType) -> Self;
    fn from_read_collection(read_collection: ClusteredReads, sort_structure: &SortStructure, layout: &LayoutType, pattern: ReadPattern) -> Self;
    fn sorted_read_set(&mut self) -> Option<ClusteredReadIterator>;
}

pub struct ClusteredMemorySortStream {
    reads: HashMap<Vec<u8>, Vec<ReadSetContainer>>,
    sort_structure: SortStructure,
    layout: LayoutType,
    pattern: ReadPattern,
}

impl ClusteredMemorySortStream {
    fn push(&mut self, read: ReadSetContainer) {
        let seq_struct = transform(read, &self.layout);
        let sequence = self.sort_structure.get_field(&seq_struct).unwrap();

        let elements = self.reads.entry(sequence).or_insert(vec![]);
        elements.push(seq_struct.original_reads().unwrap());
    }
}

impl SortStream for ClusteredMemorySortStream {
    fn from_read_iterator(read_iter: ReadIterator, sort_structure: &SortStructure, layout: &LayoutType) -> Self {
        let mut reads: HashMap<Vec<u8>, Vec<ReadSetContainer>> = HashMap::new();
        let mut mem_sort = ClusteredMemorySortStream {
            reads,
            sort_structure: sort_structure.clone(),
            layout: layout.clone(),
            pattern: ReadPattern::from_read_iterator(&read_iter),
        };
        for read in read_iter {
            mem_sort.push(read.clone());
        }
        mem_sort
    }

    fn from_read_collection(read_collection: ClusteredReads, sort_structure: &SortStructure, layout: &LayoutType, pattern: ReadPattern) -> Self {
        let mut reads: HashMap<Vec<u8>, Vec<ReadSetContainer>> = HashMap::new();
        let mut mem_sort = ClusteredMemorySortStream {
            reads,
            sort_structure: sort_structure.clone(),
            layout: layout.clone(),
            pattern: pattern.clone(),
        };
        for read in ReadIterator::from_collection(read_collection) {
            mem_sort.push(read.clone());
        }
        mem_sort
    }

    fn sorted_read_set(&mut self) -> Option<ClusteredReadIterator> {
        let max_dist = match &self.sort_structure {
            SortStructure::KNOWN_LIST { layout_type, max_distance: maximum_distance, on_disk, known_list } => { maximum_distance }
            SortStructure::HD_UMI { layout_type, max_distance, on_disk } => { max_distance }
            SortStructure::LD_UMI { layout_type, max_distance, on_disk } => { max_distance }
        };
        let collection = InputList { strings: self.reads.keys().map(|x| x.clone()).collect::<Vec<Vec<u8>>>(), max_dist: *max_dist as u64 };
        let mut graph = input_list_to_graph(&collection, string_distance, false);

        let cc = get_connected_components(&graph);

        //println!("UMI Read count {}", &cc.len());
        let mut final_vec = Vec::new();

        for group in cc {
            let minilist = InputList { strings: group, max_dist: *max_dist as u64 };
            let mut minigraph = input_list_to_graph(&minilist, string_distance, false);

            let issubgroups = split_subgroup(&mut minigraph);
            match issubgroups {
                None => {
                    let average_distance = average_dist(&minilist.strings, string_distance);
                    let mut rc = Vec::new();
                    for s in &minilist.strings {
                        if let Some(i) = self.reads.get(s) {
                            rc.extend(i.clone());
                        }
                    }
                    final_vec.push(ClusteredReads::new(VecDeque::from(rc), self.pattern.clone(), average_distance));
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
                        final_vec.push(ClusteredReads::new(VecDeque::from(rc), self.pattern.clone(), average_distance));
                    }
                }
            }
        }

        //println!("And final size = {}", &final_vec.len());
        Some(ClusteredReadIterator::new_from_vec(final_vec, self.pattern.clone()))
    }
}
