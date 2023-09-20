use std::cmp::Ordering;
use std::collections::{VecDeque};
use serde::{Serialize, Deserialize};
use crate::alignment::fasta_bit_encoding::FastaBase;

/// a sortable read set container that sorts on a set of keys -- which we populate with
/// extracted barcode sequences. These sorting sequences could have been corrected to a known list
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SortingReadSetContainer {
    pub ordered_sorting_keys: Vec<(char,Vec<FastaBase>)>,   // grow in order
    pub ordered_unsorted_keys: VecDeque<(char, Vec<FastaBase>)>, // use the default behavior to push back, pop front
    pub aligned_read: SortedAlignment,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SortedAlignment {
    pub aligned_read: Vec<FastaBase>,
    pub aligned_ref: Vec<FastaBase>,
    pub ref_name: String,
    pub read_name: String,
}

impl Eq for SortingReadSetContainer {}

impl PartialEq<Self> for SortingReadSetContainer {
    fn eq(&self, other: &Self) -> bool {
        assert_eq!(self.ordered_sorting_keys.len(), other.ordered_sorting_keys.len(), "SortingReadSetContainer: mismatched number of sorting keys from {} and {}",self.ordered_sorting_keys.len(), other.ordered_sorting_keys.len());
        self.ordered_sorting_keys.iter().zip(other.ordered_sorting_keys.iter()).map(|(a, b)| a.0 == b.0 && a.1 == b.1).count() == 0
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
        for (a, b) in self.ordered_sorting_keys.iter().zip(other.ordered_sorting_keys.iter()) {
            assert_eq!(a.0, b.0, "SortingReadSetContainer: mismatched sorting keys");
            if a.1 > b.1 {
                return Ordering::Greater;
            } else if a.1 < b.1 {
                return Ordering::Less;
            }
        }
        Ordering::Equal
    }
}

/*
impl StreamSorter {

    pub fn from_fasta(filedir: &InstanceLivedTempDir,
                reads: &dyn Iterator<Item=ReadSetContainer>,
                sort_on: UMIConfiguration) -> StreamSorter {

        let filename = filedir.path().join("my-temporary-shard-sort.txt");

        let mut writer: ShardWriter<SortingReadSetContainer> =
            ShardWriter::new(filename.clone(),
                             32,
                             256,
                             1 << 16).unwrap();

        // Get a handle to send data to the file
        let mut sender = writer.get_sender();

        // now transform individual reads into sorting read set containers by the UMIconfiguration we're sorting on

        StreamSorter {
            shard_reader: (),
            shard_writer: writer,
            sender,
            filename: filename,
        }
    }


    /// sort the iterator of ReadSetContainers by a specific capture sequence
    pub fn sort_level(reads: &dyn Iterator<Item=ReadSetContainer>, sort_on: UMIConfiguration) {}
}
*/
#[cfg(test)]
mod tests {
    use crate::alignment::fasta_bit_encoding::{FASTA_A, FASTA_N, FASTA_T};
    use super::*;
    use crate::utils::read_utils::fake_reads;


    #[test]
    fn test_ordinal_nature() {
        let srsc1 = SortingReadSetContainer{
            ordered_sorting_keys: vec![('*', vec![FASTA_A, FASTA_A])],
            ordered_unsorted_keys: Default::default(),
            aligned_read: SortedAlignment {
                aligned_read: vec![],
                aligned_ref: vec![],
                ref_name: "".to_string(),
                read_name: "".to_string(),
            },
        };
        let srsc2 = SortingReadSetContainer{
            ordered_sorting_keys: vec![('*', vec![FASTA_A, FASTA_A])],
            ordered_unsorted_keys: Default::default(),
            aligned_read: SortedAlignment {
                aligned_read: vec![],
                aligned_ref: vec![],
                ref_name: "".to_string(),
                read_name: "".to_string(),
            },
        };

        assert_eq!(srsc1.cmp(&srsc2), Ordering::Equal);

        let srsc1 = SortingReadSetContainer{
            ordered_sorting_keys: vec![('*', vec![FASTA_A, FASTA_A]),('*', vec![FASTA_A, FASTA_A])],
            ordered_unsorted_keys: Default::default(),
            aligned_read: SortedAlignment {
                aligned_read: vec![],
                aligned_ref: vec![],
                ref_name: "".to_string(),
                read_name: "".to_string(),
            },
        };
        let srsc2 = SortingReadSetContainer{
            ordered_sorting_keys: vec![('*', vec![FASTA_A, FASTA_A]),('*', vec![FASTA_A, FASTA_T])],
            ordered_unsorted_keys: Default::default(),
            aligned_read: SortedAlignment {
                aligned_read: vec![],
                aligned_ref: vec![],
                ref_name: "".to_string(),
                read_name: "".to_string(),
            },
        };

        assert_eq!(srsc1.cmp(&srsc2), Ordering::Less);
        assert_eq!(srsc2.cmp(&srsc1), Ordering::Greater);

        let srsc1 = SortingReadSetContainer{
            ordered_sorting_keys: vec![('*', vec![FASTA_A, FASTA_A]),
                                       ('*', vec![FASTA_A, FASTA_A]),
                                       ('*', vec![FASTA_A, FASTA_A])],
            ordered_unsorted_keys: Default::default(),
            aligned_read: SortedAlignment {
                aligned_read: vec![],
                aligned_ref: vec![],
                ref_name: "".to_string(),
                read_name: "".to_string(),
            },
        };
        let srsc2 = SortingReadSetContainer{
            ordered_sorting_keys: vec![('*', vec![FASTA_A, FASTA_A]),
                                       ('*', vec![FASTA_A, FASTA_T]),
                                       ('*', vec![FASTA_A, FASTA_A])],
            ordered_unsorted_keys: Default::default(),
            aligned_read: SortedAlignment {
                aligned_read: vec![],
                aligned_ref: vec![],
                ref_name: "".to_string(),
                read_name: "".to_string(),
            },
        };

        assert_eq!(srsc1.cmp(&srsc2), Ordering::Less);
        assert_eq!(srsc2.cmp(&srsc1), Ordering::Greater);


    }

    #[test]
    fn test_sorting_read_container() {
        let key1 = ('$',vec![FASTA_N, FASTA_A]);
        let key2 = ('$',vec![FASTA_N, FASTA_N]);

        let reads = fake_reads(10, 1);
        let read_seq = reads.get(0).unwrap().read_one.seq().iter().map(|x| FastaBase::from(x.clone())).collect::<Vec<FastaBase>>();
        let fake_read = SortedAlignment{
            aligned_read: read_seq.clone(),
            aligned_ref:  read_seq.clone(),
            ref_name: "".to_string(),
            read_name: "".to_string(),
        };

        let st1 = SortingReadSetContainer { ordered_sorting_keys: vec![key1.clone()], ordered_unsorted_keys: VecDeque::new(), aligned_read: fake_read.clone() };
        let st2 = SortingReadSetContainer { ordered_sorting_keys: vec![key2.clone()], ordered_unsorted_keys: VecDeque::new(), aligned_read: fake_read.clone() };
        assert!(st1 < st2);

        let st1 = SortingReadSetContainer { ordered_sorting_keys: vec![key2.clone()], ordered_unsorted_keys: VecDeque::new(),aligned_read: fake_read.clone() };
        let st2 = SortingReadSetContainer { ordered_sorting_keys: vec![key1.clone()], ordered_unsorted_keys: VecDeque::new(),aligned_read: fake_read.clone() };
        assert!(st1 > st2);

        let st1 = SortingReadSetContainer { ordered_sorting_keys: vec![key1.clone()], ordered_unsorted_keys: VecDeque::new(),aligned_read: fake_read.clone() };
        let st2 = SortingReadSetContainer { ordered_sorting_keys: vec![key1.clone()], ordered_unsorted_keys: VecDeque::new(),aligned_read: fake_read.clone() };
        println!("{} {}", st1 > st2, st1 < st2);
        assert!(!(st1 > st2) & !(st2 > st1));

        let st1 = SortingReadSetContainer { ordered_sorting_keys: vec![key1.clone(), key2.clone()], ordered_unsorted_keys: VecDeque::new(),aligned_read: fake_read.clone() };
        let st2 = SortingReadSetContainer { ordered_sorting_keys: vec![key1.clone(), key1.clone()], ordered_unsorted_keys: VecDeque::new(),aligned_read: fake_read.clone() };
        assert!(st1 > st2);
    }
}