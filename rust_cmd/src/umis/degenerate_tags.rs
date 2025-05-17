use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::hash::BuildHasherDefault;
use std::io::BufWriter;
use std::path::PathBuf;

use rustc_hash::{FxHasher, FxHashMap};
use shardio::{Range, ShardReader, ShardSender, ShardWriter};
use crate::read_strategies::read_disk_sorter::SortingReadSetContainer;
use crate::read_strategies::sequence_layout::UMIConfiguration;
use std::io::Write;
use std::ops::Deref;
use rand::distr::uniform::SampleBorrow;
use rust_star::{DistanceGraphNode, LinkedDistances};
use alignment::alignment_matrix::AlignmentResult;
use utils::read_utils::{strip_gaps, u8s};

pub struct DegenerateBuffer {
    buffer: VecDeque<SortingReadSetContainer>,
    max_buffer_size: usize,
    collapse_ratio: f64,
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
            collapse_ratio: *tag.minimum_collapsing_difference.as_ref().unwrap_or(&5.0),
            shard_writer: None,
            shard_sender: None,
            output_file,
            tag,
            hash_map: FxHashMap::default(),
        }
    }

    /// Pushes a new item onto the buffer
    /// we buffer writing to disk, to prevent the costly disk writes when we have smaller sets of barcodes.
    /// When we overflow we write the whole buffer to disk and continue to write additional reads to
    /// disk until we've finished.
    ///
    pub fn push(&mut self, mut item: SortingReadSetContainer) {
        assert!(self.tag.length >= self.tag.max_distance);

        let key_value = item.ordered_unsorted_keys.pop_front().unwrap();
        item.ordered_unsorted_keys.push_front(key_value.clone()); // we want to keep the key in the list for now, we'll remove it later
        assert_eq!(key_value.0, self.tag.symbol);


        let gapless = strip_gaps(&key_value.1);

        if gapless.len() >= self.tag.length - self.tag.max_distance && gapless.len() <= self.tag.length + self.tag.max_distance {
            *self
                .hash_map
                .entry(gapless)
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
        let sender = self.shard_sender.as_mut().unwrap();
        self.buffer.iter().for_each(|x| {
            sender.send(x.clone()).unwrap();
        });
        self.buffer.clear();
    }

    /// This function 'corrects' a list of barcodes using starcode
    pub fn correct_list(&self) -> FxHashMap<Vec<u8>, Vec<u8>> {

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

        match self.hash_map.len() {
            0 => {
                knowns
            }
            1 => {
                let kn = self.hash_map.iter().next().unwrap().0;
                knowns.insert(kn.clone(),kn.clone());
                knowns
            }
            _ => {
                let tags = self.hash_map.iter().map(|x| (x.0.clone(),*x.1)).collect::<Vec<(Vec<u8>,usize)>>();

                tags.iter().for_each(|x| println!("tag {} size {}",u8s(&x.0),&x.1));

                let correction = LinkedDistances::cluster_string_vector_list(tags, &self.tag.max_distance, &self.collapse_ratio);

                correction.into_iter().for_each(|(mut center,mut dist_graph)| {
                    let mut connected_node = dist_graph.borrow_mut();
                    let mut connected_node: &DistanceGraphNode = connected_node.deref().to_owned();
                    let string_name = connected_node.string.clone();
                    knowns.insert(string_name.clone(),string_name.clone());
                    println!("Known {} to known {}",u8s(&string_name),u8s(&string_name));
                    connected_node.swallowed_links.iter().for_each(|(x,y)| {
                        knowns.insert(x.clone(), string_name.clone());
                        println!("Unknown {} to known {}",u8s(&x),u8s(&string_name));
                    });
                });
                knowns
            }
        }

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
            let corrected = match final_correction.get(&strip_gaps(&key_value.1)) {
                None => { panic!("Unable to find match for key {} in corrected values", u8s(&key_value.1)); }
                Some(x) => { x.clone() }
            };
            y.ordered_sorting_keys
                .push((key_value.0, corrected));
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
                let corrected = match final_correction.get(&strip_gaps(&key_value.1)) {
                    None => { panic!("Unable to find match for key {} in corrected values", u8s(&key_value.1)); }
                    Some(x) => { x }
                };
                current_read.ordered_sorting_keys
                    .push((key_value.0, corrected.clone()));
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



#[cfg(test)]
mod tests {
    use bio::io::fastq::Record;
    use tempfile::tempfile;
    use ::{FASTA_A, FASTA_T};
    use read_strategies::read_set::ReadSetContainer;
    use read_strategies::sequence_layout::UMISortType;
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_tag_buffer_corrects() {
        let tempfile = NamedTempFile::new().unwrap(); // creates a new temp file
        let path: PathBuf = tempfile.path().to_path_buf(); // clone the path
        //let path_str = path.to_string_lossy(); // convert to a string-like type
        //
        //
        let config = UMIConfiguration{
            symbol: '1',
            file: None,
            reverse_complement_sequences: None,
            sort_type: UMISortType::DegenerateTag,
            length: 10,
            order: 0,
            pad: None,
            max_distance: 2,
            maximum_subsequences: Some(5000),
            minimum_collapsing_difference: Some(5.0),
            max_gaps: None,
        };


        // above the threshold for merging
        let mut tag_buffer = create_tag_buffer_with_set_anchor_seq_count(&10,&path, &config);
        let list = tag_buffer.correct_list();
        let minimum_thresh = 5.0;

        assert!(list.contains_key("AAAAATTTTT".as_bytes()));
        assert_eq!(list.get("AAAAATTTTT".as_bytes()).unwrap().clone(),"AAAAATTTTT".as_bytes().to_vec());

        assert!(list.contains_key("AAAAATTTGT".as_bytes()));
        assert_eq!(list.get("AAAAATTTGT".as_bytes()).unwrap().clone(),"AAAAATTTTT".as_bytes().to_vec());

        assert!(list.contains_key("GGGGGCCCCC".as_bytes()));
        assert_eq!(list.get("GGGGGCCCCC".as_bytes()).unwrap().clone(),"GGGGGCCCCC".as_bytes().to_vec());

        assert!(list.contains_key("GCGGGCCCCC".as_bytes()));
        assert_eq!(list.get("GCGGGCCCCC".as_bytes()).unwrap().clone(),"GGGGGCCCCC".as_bytes().to_vec());

        // below the threshold for merging
        let mut tag_buffer = create_tag_buffer_with_set_anchor_seq_count(&3,&path, &config);
        let list = tag_buffer.correct_list();

        assert!(list.contains_key("AAAAATTTTT".as_bytes()));
        assert_eq!(list.get("AAAAATTTTT".as_bytes()).unwrap().clone(),"AAAAATTTTT".as_bytes().to_vec());

        assert!(list.contains_key("AAAAATTTGT".as_bytes()));
        assert_eq!(list.get("AAAAATTTGT".as_bytes()).unwrap().clone(),"AAAAATTTGT".as_bytes().to_vec()); // not corrected

        assert!(list.contains_key("GGGGGCCCCC".as_bytes()));
        assert_eq!(list.get("GGGGGCCCCC".as_bytes()).unwrap().clone(),"GGGGGCCCCC".as_bytes().to_vec());

        assert!(list.contains_key("GCGGGCCCCC".as_bytes()));
        assert_eq!(list.get("GCGGGCCCCC".as_bytes()).unwrap().clone(),"GCGGGCCCCC".as_bytes().to_vec()); // not corrected

        let mut tag_buffer = create_tag_buffer_with_set_anchor_seq_count(&10,&path, &config);
        tag_buffer.push(create_fake_read_set_container(&"read1".to_string(),&"GGGGGCCCC-".to_string(),&config));
        tag_buffer.push(create_fake_read_set_container(&"read1".to_string(),&"GGGGCCCC--".to_string(),&config));
        tag_buffer.push(create_fake_read_set_container(&"read1".to_string(),&"GGGGGCCC--".to_string(),&config));
        let list = tag_buffer.correct_list();

        assert!(list.contains_key("GGGGGCCCC".as_bytes()));
        assert!(list.contains_key("GGGGGCCC".as_bytes()));
        assert!(list.contains_key("GGGGCCCC".as_bytes()));
        assert_eq!(list.get("GGGGGCCCC".as_bytes()).unwrap().clone(),"GGGGGCCCCC".as_bytes().to_vec()); // not corrected
        assert_eq!(list.get("GGGGGCCC".as_bytes()).unwrap().clone(),"GGGGGCCCCC".as_bytes().to_vec()); // not corrected
        assert_eq!(list.get("GGGGCCCC".as_bytes()).unwrap().clone(),"GGGGGCCCCC".as_bytes().to_vec()); // not corrected

    }

    fn create_tag_buffer_with_set_anchor_seq_count(count: &usize, path: &PathBuf, config: &UMIConfiguration) -> DegenerateBuffer {
        let mut tag_buffer = crate::umis::degenerate_tags::DegenerateBuffer::new(path.clone(),&5000, config.clone());

        for i in 0..*count {
            tag_buffer.push(create_fake_read_set_container(&"read1".to_string(), &"AAAAATTTTT".to_string(), &config));
        }
        tag_buffer.push(create_fake_read_set_container(&"read1".to_string(),&"AAAAATTTGT".to_string(),&config));

        for i in 0..*count {
            tag_buffer.push(create_fake_read_set_container(&"read1".to_string(), &"GGGGGCCCCC".to_string(), &config));
        }
        tag_buffer.push(create_fake_read_set_container(&"read1".to_string(),&"GCGGGCCCCC".to_string(),&config));

        tag_buffer
    }

    fn create_fake_read_set_container(name: &String, read_one_seq: &String, config : &UMIConfiguration) -> SortingReadSetContainer {
        let read_1_qual = vec![b'H'; read_one_seq.len()];

        let mut unsorted_keys = VecDeque::default();
        unsorted_keys.push_front((config.symbol.clone(), read_one_seq.clone().into_bytes().to_vec()));

        SortingReadSetContainer{
            ordered_sorting_keys: Vec::new(),
            ordered_unsorted_keys: unsorted_keys,
            aligned_read: AlignmentResult {
                reference_name: "default".to_string(),
                read_name: name.clone(),
                reference_aligned: vec![], // can be empty -- only the tags are used in this case
                read_aligned: vec![], // can be empty -- only the tags are used in this case
                read_quals: None, // can be empty -- only the tags are used in this case
                cigar_string: vec![], // can be empty -- only the tags are used in this case
                path: vec![], // can be empty -- only the tags are used in this case
                score: 0.0,
                reference_start: 0,
                read_start: 0,
                bounding_box: None,
            }
        }
    }
}