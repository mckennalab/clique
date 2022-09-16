use std::borrow::BorrowMut;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};

use crate::read_strategies::sequence_structures::{ReadSetContainer, ReadFileContainer, ReadIterator};

use crate::read_strategies::sequence_layout::LayoutType;
use crate::read_strategies::sequence_layout::SequenceLayout;
use crate::read_strategies::sequence_layout::transform;
use crate::sorters::known_list::KnownListSort;
use std::io::Write;
use flate2::write::GzEncoder;
use flate2::Compression;
use rust_htslib::bgzf::Writer;

// passed_fn: impl FnOnce(i32) -> ()

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum SortStructure {
    KNOWN_LIST { layout_type: LayoutType, maximum_distance: usize },
    HD_UMI { layout_type: LayoutType },
    LD_UMI { layout_type: LayoutType },
}

impl SortStructure {
    pub fn from_layout(layout: &LayoutType) -> Vec<SortStructure> {
        match layout {
            LayoutType::TENXV3 => {
                let mut ret = Vec::new();
                ret.push(SortStructure::KNOWN_LIST{ layout_type: LayoutType::TENXV3, maximum_distance: 1 });
                ret.push(SortStructure::LD_UMI{ layout_type: LayoutType::TENXV3});
                ret
            }
            LayoutType::TENXV2 => {
                let mut ret = Vec::new();
                ret.push(SortStructure::KNOWN_LIST{ layout_type: LayoutType::TENXV2, maximum_distance: 1 });
                ret.push(SortStructure::LD_UMI{ layout_type: LayoutType::TENXV2});
                ret
            }
            LayoutType::PAIREDUMI => {unimplemented!()}
            LayoutType::PAIRED => {unimplemented!()}
            LayoutType::SCI => {unimplemented!()}
            _ => {unimplemented!()}
        }
    }

    /// I couldn't get my sort structure to take a generic 'function parameter' to call for the correct sequences, so this
    /// method has to do that work in a less elegant way. Maybe I'll figure it out later
    ///
    pub fn get_known_list_value(&self, layout_type: &LayoutType, seq_layout: &dyn SequenceLayout, known_lists: HashMap<LayoutType,PathBuf>) -> PathBuf {
        match self {
            SortStructure::KNOWN_LIST{layout_type, maximum_distance} => {
                match layout_type {
                    LayoutType::TENXV3 => { known_lists.get(layout_type).unwrap().clone() }
                    LayoutType::TENXV2 => { unimplemented!() }
                    LayoutType::PAIREDUMI => { unimplemented!() }
                    LayoutType::PAIRED => { unimplemented!() }
                    LayoutType::SCI => { unimplemented!() }
                }
            }
            SortStructure::HD_UMI{layout_type} => {
                unimplemented!()
            }
            SortStructure::HD_UMI{layout_type} => {
                unimplemented!()
            }
            _ => {unimplemented!()}
        }
    }

    pub fn get_field<'a>(&self, seq_layout: &'a dyn SequenceLayout) -> Option<&'a Vec<u8>> {
        match self {
            SortStructure::KNOWN_LIST{layout_type, maximum_distance} => {
                match layout_type {
                    LayoutType::TENXV3 => { seq_layout.cell_id().to_owned() }
                    LayoutType::TENXV2 => {  seq_layout.cell_id().to_owned() }
                    LayoutType::PAIREDUMI => { unimplemented!() }
                    LayoutType::PAIRED => { unimplemented!() }
                    LayoutType::SCI => { unimplemented!() }
                }
            }
            SortStructure::HD_UMI{layout_type} => {
                match layout_type {
                    LayoutType::TENXV3 => {  seq_layout.umi().to_owned() }
                    LayoutType::TENXV2 => { seq_layout.umi().to_owned() }
                    LayoutType::PAIREDUMI => { unimplemented!() }
                    LayoutType::PAIRED => { unimplemented!() }
                    LayoutType::SCI => { unimplemented!() }
                }
            }
            SortStructure::LD_UMI{layout_type} => {
                match layout_type {
                    LayoutType::TENXV3 => { seq_layout.umi().to_owned() }
                    LayoutType::TENXV2 => { seq_layout.umi().to_owned() }
                    LayoutType::PAIREDUMI => { unimplemented!() }
                    LayoutType::PAIRED => { unimplemented!() }
                    LayoutType::SCI => { unimplemented!() }
                }
            }
        }
    }
}

pub struct ReadSortingOnDiskContainer {
    pub file_1: PathBuf,
    pub file_2: Option<PathBuf>,
    pub file_3: Option<PathBuf>,
    pub file_4: Option<PathBuf>
}


impl ReadSortingOnDiskContainer {

    // a bit ugly
    pub fn format_read_one(index: usize) -> String {format!("bucket_read1_{}.txt",index).to_string()}
    pub fn format_read_two(index: usize) -> String {format!("bucket_read2_{}.txt",index).to_string()}
    pub fn format_read_three(index: usize) -> String {format!("bucket_read3_{}.txt",index).to_string()}
    pub fn format_read_four(index: usize) -> String {format!("bucket_read4_{}.txt",index).to_string()}

    pub fn read_count(&self) -> usize {
        1 +
            if self.file_2.is_some() {1} else {0} +
            if self.file_3.is_some() {1} else {0} +
            if self.file_4.is_some() {1} else {0}
    }

    pub fn new_from_temp_dir(rd: &ReadIterator, id: usize, temp_dir: &Path) -> ReadSortingOnDiskContainer {
        ReadSortingOnDiskContainer {
            file_1: temp_dir.join(ReadSortingOnDiskContainer::format_read_one(id)),
            file_2: if rd.read_two.is_some() {Some(temp_dir.join(ReadSortingOnDiskContainer::format_read_two(id)))} else {println!("No read 2!!!"); None},
            file_3: if rd.index_one.is_some() {Some(temp_dir.join(ReadSortingOnDiskContainer::format_read_three(id)))} else {None},
            file_4: if rd.index_two.is_some() {Some(temp_dir.join(ReadSortingOnDiskContainer::format_read_four(id)))} else {None},
        }
    }

    pub fn create_x_bins(rd: &ReadIterator, x_bins: usize, temp_dir: &Path) -> Vec<(usize,ReadSortingOnDiskContainer)> {
        (0..x_bins).map(|id| (id,ReadSortingOnDiskContainer::new_from_temp_dir(rd, id, temp_dir))).collect::<Vec<(usize,ReadSortingOnDiskContainer)>>()
    }
}

pub struct SortingOutputContainer {
    file_1: BufWriter<Writer>,
    file_2: Option<BufWriter<Writer>>,
    file_3: Option<BufWriter<Writer>>,
    file_4: Option<BufWriter<Writer>>,
}

impl SortingOutputContainer {
    pub fn from_sorting_container(sc: &ReadSortingOnDiskContainer) -> SortingOutputContainer {
        SortingOutputContainer {
            file_1: BufWriter::new(Writer::from_path(&sc.file_1).unwrap()),
            file_2: if let Some(x) = &sc.file_2 {Some(BufWriter::new(Writer::from_path(&sc.file_2.as_ref().unwrap()).unwrap()))} else {None},
            file_3: if let Some(x) = &sc.file_3 {Some(BufWriter::new(Writer::from_path(&sc.file_3.as_ref().unwrap()).unwrap()))} else {None},
            file_4: if let Some(x) = &sc.file_4 {Some(BufWriter::new(Writer::from_path(&sc.file_4.as_ref().unwrap()).unwrap()))} else {None},
        }
    }

    pub fn close(&mut self)  {
        self.file_1.flush();
        if let Some(x) = &mut self.file_2 {x.flush();};
        if let Some(x) = &mut self.file_3 {x.flush();};
        if let Some(x) = &mut self.file_4 {x.flush();};
    }

    pub fn write_reads(&mut self, rl: ReadSetContainer) {
        write!(self.file_1.borrow_mut(),"{}",rl.read_one);
        if let Some(x) = self.file_2.as_mut() {
            if let Some(rd) = rl.read_two {
                write!(x, "{}", rd);
            };
        };
        if let Some(x) = self.file_3.as_mut() {
            if let Some(rd) = rl.index_one {
                write!(x, "{}", rd);
            };
        };
        if let Some(x) = self.file_4.as_mut() {
            if let Some(rd) = rl.index_two {
                write!(x, "{}", rd);
            };
        };
    }
}


pub struct SortedInputContainer {
    pub sorted_records: Vec<(String,Box<dyn SequenceLayout>)>
}

impl SortedInputContainer {
    pub fn from_sorting_container(bin: &ReadSortingOnDiskContainer, id: usize, splits: &KnownListSort, layout: &LayoutType, temp_dir: &PathBuf) -> SortedInputContainer {
        let file1 = temp_dir.join(ReadSortingOnDiskContainer::format_read_one(id));
        let file2 = temp_dir.join(ReadSortingOnDiskContainer::format_read_two(id));
        let file3 = temp_dir.join(ReadSortingOnDiskContainer::format_read_three(id));
        let file4 = temp_dir.join(ReadSortingOnDiskContainer::format_read_four(id));
        let container = ReadFileContainer{
            read_one: file1,
            read_two: file2,
            index_one: file3,
            index_two: file4
        };

        let read_iterator = ReadIterator::new_from_bundle(&container);
        println!("Opening {:?}",container.read_one);
        let mut return_reads = Vec::new();
        for rd in read_iterator {
            let transformed_reads = transform(rd, layout);
            assert!(transformed_reads.has_original_reads());

            let first_barcode = String::from_utf8(transformed_reads.cell_id().unwrap().clone()).unwrap();
            return_reads.push((first_barcode,transformed_reads));
        }
        println!("SORTING");
        return_reads.sort_by(|a,b| b.0.cmp(&a.0));
        println!("SORTING DONE");
        SortedInputContainer{ sorted_records: return_reads }
    }
}