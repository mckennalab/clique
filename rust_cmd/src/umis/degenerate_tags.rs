use std::collections::{HashMap, VecDeque};
use std::hash::BuildHasherDefault;
use std::path::PathBuf;

use rustc_hash::{FxHasher, FxHashMap};
use shardio::{ShardSender, ShardWriter};
use crate::alignment::fasta_bit_encoding::{FastaBase};
use crate::read_strategies::read_disk_sorter::SortingReadSetContainer;
use crate::read_strategies::sequence_layout::UMIConfiguration;

use rust_starcode::StarcodeAlignment;

pub struct DegenerateBuffer {
    buffer: VecDeque<SortingReadSetContainer>,
    max_buffer_size: usize,
    writen_reads: usize,
    shard_writer: Option<Box<ShardWriter<SortingReadSetContainer>>>,
    shard_sender: Option<Box<ShardSender<SortingReadSetContainer>>>,
    output_file: PathBuf,
    tag: UMIConfiguration,
    hash_map: FxHashMap<Vec<u8>, usize>,
}

impl DegenerateBuffer {
    pub fn new(output_file: PathBuf, max_size: &usize, tag: UMIConfiguration) -> DegenerateBuffer {
        DegenerateBuffer {
            buffer: VecDeque::new(),
            max_buffer_size: *max_size,
            writen_reads: 0,
            shard_writer: None,
            shard_sender: None,
            output_file,
            tag,
            hash_map: FxHashMap::default(),
        }
    }

    pub fn bin_size_exceeded(&self) -> bool {
        self.buffer.len() >= self.max_buffer_size
    }

    /// Pushes a new item onto the buffer
    /// we buffer writing to prevent the costly disk writes when we have smaller sets of barcodes.
    /// When we overflow we write the whole buffer to disk and continue to write additional reads to
    /// disk until we've finished.
    ///
    pub fn push(&mut self, mut item: SortingReadSetContainer) {
        let key_value = item.ordered_unsorted_keys.pop_front().unwrap();
        item.ordered_unsorted_keys.push_front(key_value.clone()); // we want to keep the key in the list for now, we'll remove it later
        assert_eq!(key_value.0, self.tag.symbol);


        let gaplengthgapless = FastaBase::strip_gaps(&key_value.1);

        *self
            .hash_map
            .entry(FastaBase::vec_u8(&gaplengthgapless))
            .or_insert(0) += 1;

        match (
            &self.shard_writer,
            self.writen_reads >= self.max_buffer_size,
        ) {
            (None, true) => {
                self.buffer.push_back(item);
                self.dump_buffer();
            }
            (None, false) => {
                self.buffer.push_back(item);
            }
            (Some(_), _) => {
                self.shard_sender.as_mut().unwrap().send(item).unwrap();
            }
        }
    }

    /// Dumps the buffer to disk
    pub fn dump_buffer(&mut self) {
        self.shard_writer = Some(Box::new(
            ShardWriter::new(&self.output_file, 32, 256, 1 << 16).unwrap(),
        ));
        self.shard_sender = Some(Box::new(self.shard_writer.as_mut().unwrap().get_sender()));
        self.buffer.iter().for_each(|x| {
            self.shard_sender.as_mut().unwrap().send(x.clone()).unwrap();
        });
    }

    /// This function 'corrects' a list of barcodes using starcode
    pub fn correct_list(&self) -> FxHashMap<Vec<u8>, Vec<u8>> {
        println!("Correcting list of length {}",self.hash_map.len());
        self.hash_map.iter().for_each(|(k,v)| {
            for x in k {
                match x {
                    &b'a' | &b'A' | &b'c' | &b'C' | &b'g' | &b'G' | &b't' | &b'T' => {

                    },
                    _  => {
                        println!("Invalid character {} in {}",x,String::from_utf8(k.clone()).unwrap());
                    }
                }
            }
        });
        let correction = StarcodeAlignment::align_sequences(&self.hash_map, &(i32::try_from(self.tag.max_distance).unwrap()), &3.0);
        let mut knowns: FxHashMap<Vec<u8>,Vec<u8>> = FxHashMap::default();
        for i in 0..correction.cluster_centers.len() {
            let center = correction.cluster_centers.get(i).unwrap();
            correction.cluster_members.get(i).unwrap().iter().for_each(|v| {
                assert!(!knowns.contains_key(v));
                //println!("COrrection center {} for {}, count {}",String::from_utf8(center.clone()).unwrap(),String::from_utf8(v.clone()).unwrap(),correction.cluster_count.get(i).unwrap());
                knowns.insert(v.clone(),center.clone());
            });
        };
        knowns
    }

    pub fn close(
        &mut self,
    ) -> (
        VecDeque<SortingReadSetContainer>,
        FxHashMap<Vec<u8>, Vec<u8>>,
    ) {
        let final_correction = self.correct_list();
        let ret = (
            std::mem::replace(&mut self.buffer, VecDeque::new()),
            final_correction,
        );
        self.clean();
        ret
    }

    pub fn clean(&mut self) {
        self.buffer.clear();
        self.hash_map.clear();
    }
}

trait ListCorrector {
    fn correct_list(&self, counts: &FxHashMap<String, usize>) -> HashMap<Vec<u8>, Vec<u8>, BuildHasherDefault<FxHasher>>;
    fn get_max_distance(&self) -> u32;
}
