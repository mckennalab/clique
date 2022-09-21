use std::collections::{HashMap, VecDeque};
use std::slice::Iter;

use crate::consensus::consensus_builders::create_poa_consensus;
use crate::read_strategies::sequence_file_containers::{ReadIterator, ReadSetContainer};
use crate::read_strategies::sequence_file_containers::NamedReadSetContainer;
use crate::read_strategies::sequence_file_containers::ReadCollection;
use crate::read_strategies::sequence_file_containers::ReadCollectionIterator;
use crate::read_strategies::sequence_layout::LayoutType;
use crate::read_strategies::sequence_layout::SequenceLayout;
use crate::read_strategies::sequence_layout::transform;
use crate::sorters::sorter::SortStructure;
use crate::umis::sequence_clustering::{split_subgroup, string_distance};
use crate::umis::sequence_clustering::input_list_to_graph;
use crate::umis::sequence_clustering::InputList;
use crate::read_strategies::sequence_file_containers::ReadPattern;


pub trait SortStream {
    fn from_read_iterator(read_iter: ReadIterator, sort_structure: &SortStructure, layout: &LayoutType) -> Self;
    fn from_read_collection(read_collection: ReadCollection, sort_structure: &SortStructure, layout: &LayoutType, pattern: ReadPattern) -> Self;
    fn sorted_read_set(&mut self) -> Option<ReadCollectionIterator>;
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

    fn from_read_collection(read_collection: ReadCollection, sort_structure: &SortStructure, layout: &LayoutType, pattern: ReadPattern) -> Self {
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

    fn sorted_read_set(&mut self) -> Option<ReadCollectionIterator> {
        let collection = InputList { strings: self.reads.keys().map(|x| x.clone()).collect::<Vec<Vec<u8>>>(), max_dist: 6 };
        let mut graph = input_list_to_graph(&collection, string_distance, true);
        let groups = split_subgroup(&mut graph);

        let mut final_vec = Vec::new();

        match groups {
            None => {}
            Some(grouping) => {
                for group in grouping {
                    let mut rc = Vec::new();

                    for tag in &group {
                        if let Some(i) = self.reads.get(tag) {
                            rc.extend(i.clone());
                        }
                    }
                    let consensus_umi = create_poa_consensus(&group);
                    final_vec.push(ReadCollection::new( VecDeque::from(rc),  ReadPattern::ONE ));
                }
            }
        }


        Some(ReadCollectionIterator::new_from_vec(final_vec, self.pattern.clone()))
    }
}
