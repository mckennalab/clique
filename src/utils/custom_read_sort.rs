use std::borrow::{Borrow, BorrowMut};
use std::fs::{File, remove_file};
use std::io::{BufReader, BufWriter, Lines, prelude::*};
use std::mem;
use std::path::{Path, PathBuf};

use bio::io::fastq::{Reader, Record, Records};

use crate::read_strategies::sequence_layout::*;
use crate::read_strategies::sequence_structures::{ReadIterator, ReadSetContainer};
use crate::read_strategies::sequence_structures::*;
use crate::sorters::known_list::KnownListSort;
use crate::sorters::sorter::SortStructure;

const BUFFER_CAPACITY: usize = 4_000_000_000;
const MAX_MEM_USE: usize = 4_000_000_000;

/*
fn main() {
    let filename = "logistic-chaos.txt";
    create_large_file(8_000_000_000, filename);
    extern_sort(filename, MAX_MEM_USE);
}
*/
/*
impl SortedInputContainer {
    pub fn from_sorting_container(sc: &ReadSortingFileContainer, layout: &LayoutType, cms: &KnownListSort) -> SortedInputContainer {
        let read_iterator = ReadIterator::new_from_on_disk_sorter(sc);

        let mut to_be_sorted = Vec::new();

        for rd in read_iterator {
            let transformed_reads = transform(rd, layout);
            assert!(transformed_reads.has_original_reads());

            let first_barcode = transformed_reads.umi_sorting_path().unwrap()[0].clone();

            let mapped_to_id = cms.original_to_best_match.get(&first_barcode).unwrap();
            to_be_sorted.push((String::from_utf8(mapped_to_id.clone()).unwrap(),transformed_reads));
        }
        to_be_sorted.sort_by(|a, b|a.0.partial_cmp(&b.0).unwrap());

        SortedInputContainer{ sorted_records: to_be_sorted }
    }
}
*/