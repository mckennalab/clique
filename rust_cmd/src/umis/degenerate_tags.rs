use std::collections::{HashMap, VecDeque};
use std::hash::BuildHasherDefault;
use std::path::PathBuf;

use rustc_hash::{FxHasher, FxHashMap};
use shardio::{Range, ShardReader, ShardSender, ShardWriter};
use crate::alignment::fasta_bit_encoding::{FastaBase};
use crate::read_strategies::read_disk_sorter::SortingReadSetContainer;
use crate::read_strategies::sequence_layout::UMIConfiguration;

use rust_starcode::StarcodeAlignment;

pub struct DegenerateBuffer {
    buffer: VecDeque<SortingReadSetContainer>,
    max_buffer_size: usize,
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
            shard_writer: None,
            shard_sender: None,
            output_file,
            tag,
            hash_map: FxHashMap::default(),
        }
    }

    /// Pushes a new item onto the buffer
    /// we buffer writing to prevent the costly disk writes when we have smaller sets of barcodes.
    /// When we overflow we write the whole buffer to disk and continue to write additional reads to
    /// disk until we've finished.
    ///
    pub fn push(&mut self, mut item: SortingReadSetContainer) {
        assert!(self.tag.length >= self.tag.max_distance);

        let key_value = item.ordered_unsorted_keys.pop_front().unwrap();
        item.ordered_unsorted_keys.push_front(key_value.clone()); // we want to keep the key in the list for now, we'll remove it later
        assert_eq!(key_value.0, self.tag.symbol);


        let gapless = FastaBase::strip_gaps(&key_value.1);

        if gapless.len() >= self.tag.length - self.tag.max_distance && gapless.len() <= self.tag.length + self.tag.max_distance {
            *self
                .hash_map
                .entry(FastaBase::vec_u8(&gapless))
                .or_insert(0) += 1;

            //println!("{} {} {} {}",self.shard_writer.is_some(), self.buffer.len(), self.max_buffer_size,self.buffer.len() >= self.max_buffer_size);
            match (
                &self.shard_sender.is_some(),
                self.buffer.len() >= self.max_buffer_size,
            ) {
                (false, true) => {
                    self.buffer.push_back(item);
                    self.dump_buffer_to_disk();
                }
                (false, false) => {
                    self.buffer.push_back(item);
                }
                (true, _) => {
                    self.shard_sender.as_mut().unwrap().send(item).unwrap();
                }
            }
        } else {
            // TODO record missing reads here
        }
    }
    /// Dumps the buffer to disk
    pub fn dump_buffer_to_disk(&mut self) {
        self.shard_writer = Some(Box::new(
            ShardWriter::new(&self.output_file, 32, 256, 1 << 16).unwrap(),
        ));
        self.shard_sender = Some(Box::new(self.shard_writer.as_mut().unwrap().get_sender()));
        self.buffer.iter().for_each(|x| {
            self.shard_sender.as_mut().unwrap().send(x.clone()).unwrap();
        });
        self.buffer.clear();
    }

    /// This function 'corrects' a list of barcodes using starcode
    pub fn correct_list(&self) -> FxHashMap<Vec<u8>, Vec<u8>> {
        //println!("Correcting list of length {}",self.hash_map.len());
        self.hash_map.iter().for_each(|(k, _v)| {
            for x in k {
                match x {
                    &b'a' | &b'A' | &b'c' | &b'C' | &b'g' | &b'G' | &b't' | &b'T' => {}
                    _ => {
                        println!("Invalid character {} in {}", x, String::from_utf8(k.clone()).unwrap());
                    }
                }
            }
        });

        let mut knowns: FxHashMap<Vec<u8>, Vec<u8>> = FxHashMap::default();

        if self.hash_map.len() > 0 {
            let correction = StarcodeAlignment::align_sequences(&self.hash_map, &(i32::try_from(self.tag.max_distance + 1).unwrap()), &3.0); // TODO: the plus 1 here is a patch until it's clear why starcode distance metrics are weird
            for i in 0..correction.cluster_centers.len() {
                let center = correction.cluster_centers.get(i).unwrap();
                correction.cluster_members.get(i).unwrap().iter().for_each(|v| {
                    assert!(!knowns.contains_key(v));
                    //println!("COrrection center {} for {}, count {}",String::from_utf8(center.clone()).unwrap(),String::from_utf8(v.clone()).unwrap(),correction.cluster_count.get(i).unwrap());
                    knowns.insert(v.clone(), center.clone());
                });
            };
        }
        knowns

    }

    pub fn close_and_write_to_shard_writer(
        &mut self,
        sender: &mut ShardSender<SortingReadSetContainer>,
    ) -> usize {
        let mut read_count: usize = 0;
        let mut buffered_reads = 0;
        let mut unbuffered_reads = 0;
        let final_correction = self.correct_list();

        self.buffer.iter().for_each(|y| {
            let mut y = y.clone();
            let key_value = y.ordered_unsorted_keys.pop_front().unwrap();
            let corrected = match final_correction.get(&FastaBase::vec_u8_strip_gaps(&key_value.1)) {
                None => { panic!("Unable to find match for key {} in corrected values", FastaBase::string(&key_value.1)); }
                Some(x) => { x }
            };
            y.ordered_sorting_keys
                .push((key_value.0, FastaBase::from_vec_u8(corrected)));
            sender.send(y).unwrap();
            read_count += 1;
            buffered_reads += 1;
        });

        if self.shard_writer.is_some() {
            self.shard_sender.as_mut().unwrap().finished().unwrap();
            self.shard_writer.as_mut().unwrap().as_mut().finish().unwrap();
            let reader: ShardReader<SortingReadSetContainer> = ShardReader::open(&self.output_file).unwrap();
            reader.iter_range(&Range::all()).unwrap().for_each(|current_read| {
                let mut current_read: SortingReadSetContainer = current_read.unwrap();
                let key_value = current_read.ordered_unsorted_keys.pop_front().unwrap();
                let corrected = match final_correction.get(&FastaBase::vec_u8_strip_gaps(&key_value.1)) {
                    None => { panic!("Unable to find match for key {} in corrected values", FastaBase::string(&key_value.1)); }
                    Some(x) => { x }
                };
                current_read.ordered_sorting_keys
                    .push((key_value.0, FastaBase::from_vec_u8(corrected)));
                sender.send(current_read).unwrap();
                read_count += 1;
                unbuffered_reads += 1;
            });
        }

        // clear everything out
        self.shard_sender = None;
        self.shard_writer = None;
        self.buffer.clear();
        self.hash_map.clear();

        //println!("COUNTS {} {} {}",read_count,buffered_reads,unbuffered_reads);
        read_count
    }
}

trait ListCorrector {
    fn correct_list(&self, counts: &FxHashMap<String, usize>) -> HashMap<Vec<u8>, Vec<u8>, BuildHasherDefault<FxHasher>>;
    fn get_max_distance(&self) -> u32;
}
