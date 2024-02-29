use std::collections::{HashMap, VecDeque};
use std::path::PathBuf;
use rustc_hash::FxHashMap;
use shardio::{ShardSender, ShardWriter};
use crate::alignment::fasta_bit_encoding::FastaBase;
use crate::read_strategies::read_disk_sorter::SortingReadSetContainer;
use crate::read_strategies::sequence_layout::UMIConfiguration;
use crate::umis::known_list::FastaString;
use crate::umis::sequence_clustering::{get_connected_components, InputList, vantage_point_string_graph};

pub struct DegenerateBuffer {
    buffer: VecDeque<SortingReadSetContainer>,
    max_buffer_size: usize,
    writen_reads: usize,
    shard_writer: Option<Box<ShardWriter<SortingReadSetContainer>>>,
    shard_sender: Option<Box<ShardSender<SortingReadSetContainer>>>,
    output_file: PathBuf,
    tag: UMIConfiguration,
    hash_map: FxHashMap<String, usize>,
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
        *self
            .hash_map
            .entry(FastaBase::string(&key_value.1))
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

    /// This function 'corrects' a list of barcodes using a connected components algorithm
    pub fn correct_list(&self) -> FxHashMap<Vec<u8>, Vec<u8>> {
        let string_set = Vec::from_iter(self.hash_map.keys())
            .iter()
            .map(|s| s.as_bytes().to_vec())
            .collect::<Vec<Vec<u8>>>();
        let collection = InputList {
            strings: string_set,
            max_dist: self.tag.max_distance.clone(),
        };
        let graph = vantage_point_string_graph(&collection, self.hash_map.len() > 50000);

        let cc = get_connected_components(&graph);
        let mut final_correction: FxHashMap<Vec<u8>, Vec<u8>> = FxHashMap::default();

        for group in cc {
            let minilist = InputList {
                strings: group,
                max_dist: self.tag.max_distance.clone(),
            };

            // TODO enable this for a future improvement -- spliting connected components -- but a lot of work needs to go into validation here
            /*let mut minigraph = input_list_to_graph(&minilist, string_distance, false);

            let is_subgroups = split_subgroup(&mut minigraph);

            match is_subgroups {
                None => {*/
            let group = minilist.strings.clone();
            let conc = crate::collapse::consensus(&group); // we need to do this beforehand -- we can't do it in the loop below because we're consuming the list below which changes the results
            for s in minilist.strings {
                //info!("correct_list from {} -> {}", String::from_utf8(s.clone()).unwrap(), String::from_utf8(conc.clone()).unwrap());
                final_correction.insert(s.clone(), conc.clone());
            }
            /*
                }
                Some(x) => {
                    for sgroup in &x {
                        for s in sgroup {
                            final_correction.insert(s.clone(), consensus(&sgroup));
                        }
                    }
                }
            }*/
        }
        final_correction
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
    fn correct_list(counts: &HashMap<String, usize>) -> HashMap<Vec<u8>, Vec<u8>>;
}

struct MLCorrector {}

impl ListCorrector for MLCorrector {
    fn correct_list(counts: &HashMap<String, usize>) -> HashMap<Vec<u8>, Vec<u8>> {
        let keys:Vec<FastaString> = counts.iter().map(|(key,value)|
            FastaString::new(FastaBase::from_str(key.as_str()))
        ).collect();

        let vantage_tree = vpsearch::Tree::new(&keys);


        let return_map = HashMap::default();
        return_map

    }
}


