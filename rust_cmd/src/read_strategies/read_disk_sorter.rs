use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs;
use std::fs::File;
use std::hash::BuildHasherDefault;
use std::path::{Path, PathBuf};
use flate2::bufread::GzEncoder;
use nohash_hasher::NoHashHasher;
use crate::read_strategies::read_set::ReadSetContainer;
use serde::ser::{SerializeSeq, Serializer};
use serde::{Serialize, Deserialize};
use tempfile::TempPath;
use crate::alignment::fasta_bit_encoding::FastaBase;
use shardio::*;
use crate::LivedTempDir;
use crate::read_strategies::sequence_layout::{SequenceLayoutDesign, UMIConfiguration, UMISortType};
use crate::umis::sequence_clustering::StringGraph;


/// a sortable read set container that sorts on a set of keys -- which we populate with
/// extracted barcode sequences. These sorting sequences could have been corrected to a known list
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SortingReadSetContainer {
    sorting_keys: Vec<Vec<FastaBase>>,
    reads: ReadSetContainer,
}

impl Eq for SortingReadSetContainer {}

impl PartialEq<Self> for SortingReadSetContainer {
    fn eq(&self, other: &Self) -> bool {
        assert_eq!(self.sorting_keys.len(), other.sorting_keys.len());
        self.sorting_keys.iter().zip(other.sorting_keys.iter()).map(|(a, b)| a == b).count() == 0
    }
}

impl PartialOrd for SortingReadSetContainer {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Right now we simply offer a stable sort for bases -- there's no natural ordering to nucleotides;
/// you could argue alphabetical, but we simply sort on their underlying bit encoding
impl Ord for SortingReadSetContainer {
    fn cmp(&self, other: &Self) -> Ordering {
        for (a, b) in self.sorting_keys.iter().zip(other.sorting_keys.iter()) {
            if a > b {
                return Ordering::Greater;
            }
            else if a < b {
                return Ordering::Less;
            }
        }
        Ordering::Equal
    }
}

pub struct StreamSorter {
    shard_writer: ShardWriter<SortingReadSetContainer>,
    sender: ShardSender<SortingReadSetContainer>,
    number_of_splits: usize,
    filename: PathBuf,
}

impl StreamSorter {
    pub fn from(filedir: &LivedTempDir, splits: &usize) -> StreamSorter {
        let filename = filedir.path().join("my-temporary-shard-sort.txt");
        let mut writer: ShardWriter<SortingReadSetContainer> = ShardWriter::new(filename.clone(), 64, 256, 1<<16).unwrap();

        // Get a handle to send data to the file
        let mut sender = writer.get_sender();

        StreamSorter{
            shard_writer: writer,
            sender,
            number_of_splits: 0,
            filename: filename,
        }
    }


    pub fn preprocess_reads(reads: &dyn Iterator<Item=ReadSetContainer>, sld: &SequenceLayoutDesign) { // -> HashMap<char,StringGraph> {
        let hash_characters = sld.umi_configurations.iter().map(|(k,v)| (v.symbol,v.sort_type == UMISortType::DegenerateTag)).filter(|(c,isd)| *isd).collect::<Vec<(char,bool)>>();
        let mut hashedvalues: HashMap::<char, StringGraph, BuildHasherDefault<NoHashHasher<u8>>> = HashMap::with_capacity_and_hasher(hash_characters.len(), BuildHasherDefault::default());



    }


    /// sort the iterator of ReadSetContainers by a specific capture sequence
    pub fn sort_level(reads: &dyn Iterator<Item=ReadSetContainer>, sort_on: UMIConfiguration) {


    }

}

#[cfg(test)]
mod tests {
    use crate::alignment::fasta_bit_encoding::{FASTA_A, FASTA_N};
    use super::*;
    use crate::utils::read_utils::fake_reads;

    #[test]
    fn test_sorting_read_container() {
        let key1 = vec![FASTA_N,FASTA_A];
        let key2 = vec![FASTA_N,FASTA_N];

        let reads = fake_reads(10, 1);

        let st1 = SortingReadSetContainer{ sorting_keys: vec![key1.clone()], reads: reads.get(0).unwrap().clone() };
        let st2 = SortingReadSetContainer{ sorting_keys: vec![key2.clone()], reads: reads.get(1).unwrap().clone() };
        assert!(st1 < st2);

        let st1 = SortingReadSetContainer{ sorting_keys: vec![key2.clone()], reads: reads.get(0).unwrap().clone() };
        let st2 = SortingReadSetContainer{ sorting_keys: vec![key1.clone()], reads: reads.get(1).unwrap().clone() };
        assert!(st1 > st2);

        let st1 = SortingReadSetContainer{ sorting_keys: vec![key1.clone()], reads: reads.get(0).unwrap().clone() };
        let st2 = SortingReadSetContainer{ sorting_keys: vec![key1.clone()], reads: reads.get(1).unwrap().clone() };
        println!("{} {}",st1> st2, st1 < st2);
        assert!(!(st1 > st2) & !(st2 > st1));

        let st1 = SortingReadSetContainer{ sorting_keys: vec![key1.clone(),key2.clone()], reads: reads.get(0).unwrap().clone() };
        let st2 = SortingReadSetContainer{ sorting_keys: vec![key1.clone(),key1.clone()], reads: reads.get(1).unwrap().clone() };
        assert!(st1 > st2);
    }
}