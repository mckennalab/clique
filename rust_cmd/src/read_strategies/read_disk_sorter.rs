use std::cmp::Ordering;
use std::collections::{VecDeque};
use serde::{Serialize, Deserialize};
use crate::alignment::alignment_matrix::{AlignmentResult};

/// a sortable read set container that sorts on a set of keys -- which we populate with
/// extracted barcode sequences. These sorting sequences could have been corrected to a known list
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct CorrectedKey {
    pub key: char,
    pub original: Vec<u8>,
    pub corrected: Vec<u8>,
}

impl CorrectedKey {
    #[allow(dead_code)]
    pub fn new(key: char, original: Vec<u8>, corrected: Vec<u8>) -> Self {
        Self { key, original, corrected }
    }
}
impl Eq for CorrectedKey {}

impl PartialEq<Self> for CorrectedKey {
    fn eq(&self, other: &Self) -> bool {
        self.corrected.eq(&other.corrected)
    }
}

impl PartialOrd<Self> for CorrectedKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.corrected.partial_cmp(&other.corrected)
    }
}

impl Ord for CorrectedKey {
    fn cmp(&self, other: &Self) -> Ordering {
        self.corrected.cmp(&other.corrected)
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SortingReadSetContainer {
    pub ordered_sorting_keys: Vec<(char,CorrectedKey)>,   // contains an ordered list of how this SortingReadSetContainer has been sorted so far
    pub ordered_unsorted_keys: VecDeque<(char, Vec<u8>)>, // a list of keys in this SortingReadSetContainer that have yet to be used to sort the container
    pub aligned_read: AlignmentResult,
}
impl SortingReadSetContainer {
    pub fn with_new_alignment(&self,new_alignment: AlignmentResult) -> SortingReadSetContainer {
        SortingReadSetContainer{
            ordered_sorting_keys: self.ordered_sorting_keys.clone(),
            ordered_unsorted_keys: self.ordered_unsorted_keys.clone(),
            aligned_read: new_alignment,
        }
    }
    pub fn empty_tags(new_alignment: AlignmentResult) -> SortingReadSetContainer {
        SortingReadSetContainer{
            ordered_sorting_keys: vec![],
            ordered_unsorted_keys: VecDeque::new(),
            aligned_read: new_alignment,
        }
    }


}
impl Eq for SortingReadSetContainer {}

impl PartialEq<Self> for SortingReadSetContainer {
    fn eq(&self, other: &Self) -> bool {
        assert_eq!(self.ordered_sorting_keys.len(), other.ordered_sorting_keys.len(), "SortingReadSetContainer: mismatched number of sorting keys from {} and {}",self.ordered_sorting_keys.len(), other.ordered_sorting_keys.len());
        for (a,b) in self.ordered_sorting_keys.iter().zip(other.ordered_sorting_keys.iter()) {
            if !(a.0 == b.0 && a.1 == b.1) {
                return false
            }
        }
        true
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
        match self.aligned_read.reference_name.cmp(&other.aligned_read.reference_name) {
            Ordering::Less => Ordering::Less,
            Ordering::Equal => {
                for (a, b) in self.ordered_sorting_keys.iter().zip(other.ordered_sorting_keys.iter()) {
                    assert_eq!(a.0, b.0, "SortingReadSetContainer: mismatched sorting keys");
                    match a.1.cmp(&b.1) {
                        Ordering::Less => {return Ordering::Less}
                        Ordering::Greater => {return Ordering::Greater}
                        Ordering::Equal => {/* do nothing, someone else will solve our problem*/}
                    }
                }
                Ordering::Equal
            }
            Ordering::Greater => Ordering::Greater,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::cmp::Ordering;
    use std::collections::VecDeque;
    use ::{FASTA_A, FASTA_N};
    use FASTA_T;
    use read_strategies::read_disk_sorter::CorrectedKey;
    use crate::alignment::alignment_matrix::AlignmentResult;
    use crate::read_strategies::read_disk_sorter::{SortingReadSetContainer};
    use crate::utils::read_utils::fake_reads;

    #[test]
    fn test_ordinal_nature() {
        let srsc1 = SortingReadSetContainer{
            ordered_sorting_keys: vec![('*', CorrectedKey::new('*', vec![FASTA_A, FASTA_A],vec![FASTA_A, FASTA_A]))],
            ordered_unsorted_keys: Default::default(),
            aligned_read: AlignmentResult {
                reference_name: "".to_string(),
                read_name: "".to_string(),
                reference_aligned: vec![],
                read_aligned: vec![],
                read_quals: None,
                cigar_string: vec![],
                path: vec![],
                score: 0.0,
                reference_start: 0,
                read_start: 0,
                bounding_box: None,
            }
        };
        let srsc2 = SortingReadSetContainer{
            ordered_sorting_keys: vec![('*', CorrectedKey::new('*', vec![FASTA_A, FASTA_A],vec![FASTA_A, FASTA_A]))],
            ordered_unsorted_keys: Default::default(),
            aligned_read: AlignmentResult {
                reference_name: "".to_string(),
                read_name: "".to_string(),
                reference_aligned: vec![],
                read_aligned: vec![],
                read_quals: None,
                cigar_string: vec![],
                path: vec![],
                score: 0.0,
                reference_start: 0,
                read_start: 0,
                bounding_box: None,
            }
        };

        assert_eq!(srsc1.cmp(&srsc2), Ordering::Equal);

        let srsc1 = SortingReadSetContainer{
            ordered_sorting_keys: vec![('*', CorrectedKey::new('*', vec![FASTA_A, FASTA_A],vec![FASTA_A, FASTA_A])),('*', CorrectedKey::new('*', vec![FASTA_A, FASTA_A],vec![FASTA_A, FASTA_A]))],
            ordered_unsorted_keys: Default::default(),
            aligned_read: AlignmentResult {
                reference_name: "".to_string(),
                read_name: "".to_string(),
                reference_aligned: vec![],
                read_aligned: vec![],
                read_quals: None,
                cigar_string: vec![],
                path: vec![],
                score: 0.0,
                reference_start: 0,
                read_start: 0,
                bounding_box: None,
            }
        };
        let srsc2 = SortingReadSetContainer{
            ordered_sorting_keys: vec![('*', CorrectedKey::new('*', vec![FASTA_A, FASTA_A],vec![FASTA_A, FASTA_A])),('*', CorrectedKey::new('*', vec![FASTA_A, FASTA_A],vec![FASTA_A, FASTA_T]))],
            ordered_unsorted_keys: Default::default(),
            aligned_read: AlignmentResult {
                reference_name: "".to_string(),
                read_name: "".to_string(),
                reference_aligned: vec![],
                read_aligned: vec![],
                read_quals: None,
                cigar_string: vec![],
                path: vec![],
                score: 0.0,
                reference_start: 0,
                read_start: 0,
                bounding_box: None,
            }
        };

        assert_eq!(srsc1.cmp(&srsc2), Ordering::Less);
        assert_eq!(srsc2.cmp(&srsc1), Ordering::Greater);

        let srsc1 = SortingReadSetContainer{
            ordered_sorting_keys: vec![('*', CorrectedKey::new('*', vec![FASTA_A, FASTA_A],vec![FASTA_A, FASTA_A])),
                                       ('*', CorrectedKey::new('*', vec![FASTA_A, FASTA_A],vec![FASTA_A, FASTA_A])),
                                       ('*', CorrectedKey::new('*', vec![FASTA_A, FASTA_A],vec![FASTA_A, FASTA_A]))],
            ordered_unsorted_keys: Default::default(),
            aligned_read: AlignmentResult {
                reference_name: "".to_string(),
                read_name: "".to_string(),
                reference_aligned: vec![],
                read_aligned: vec![],
                read_quals: None,
                cigar_string: vec![],
                path: vec![],
                score: 0.0,
                reference_start: 0,
                read_start: 0,
                bounding_box: None,
            }
        };
        let srsc2 = SortingReadSetContainer{
            ordered_sorting_keys: vec![('*', CorrectedKey::new('*', vec![FASTA_A, FASTA_A],vec![FASTA_A, FASTA_A])),
                                       ('*', CorrectedKey::new('*', vec![FASTA_A, FASTA_A],vec![FASTA_A, FASTA_T])),
                                       ('*', CorrectedKey::new('*', vec![FASTA_A, FASTA_A],vec![FASTA_A, FASTA_A]))],
            ordered_unsorted_keys: Default::default(),
            aligned_read: AlignmentResult {
                reference_name: "".to_string(),
                read_name: "".to_string(),
                reference_aligned: vec![],
                read_aligned: vec![],
                read_quals: None,
                cigar_string: vec![],
                path: vec![],
                score: 0.0,
                reference_start: 0,
                read_start: 0,
                bounding_box: None,
            }
        };

        assert_eq!(srsc1.cmp(&srsc2), Ordering::Less);
        assert_eq!(srsc2.cmp(&srsc1), Ordering::Greater);


    }

    #[test]
    fn test_sorting_read_container() {
        let key1 = ('*', CorrectedKey::new('*', vec![FASTA_N, FASTA_A],vec![FASTA_N, FASTA_A]));
        let key2 = ('*', CorrectedKey::new('*', vec![FASTA_N, FASTA_N],vec![FASTA_N, FASTA_N]));

        let reads = fake_reads(10, 1);
        let read_seq = reads.get(0).unwrap().read_one.seq().iter().map(|x| x.clone()).collect::<Vec<u8>>();
        let fake_read = AlignmentResult{
            reference_name: "".to_string(),
            read_name: "".to_string(),
            reference_aligned: read_seq.clone(),
            cigar_string: vec![],
            path: vec![],
            score: 0.0,
            reference_start: 0,
            read_start: 0,
            read_aligned: read_seq.clone(),
            bounding_box: None,
            read_quals: None,
        };

        let st1 = SortingReadSetContainer { ordered_sorting_keys: vec![key1.clone()], ordered_unsorted_keys: VecDeque::new(), aligned_read: fake_read.clone() };
        let st2 = SortingReadSetContainer { ordered_sorting_keys: vec![key2.clone()], ordered_unsorted_keys: VecDeque::new(), aligned_read: fake_read.clone() };
        assert!(st1 < st2);

        let st1 = SortingReadSetContainer { ordered_sorting_keys: vec![key2.clone()], ordered_unsorted_keys: VecDeque::new(),aligned_read: fake_read.clone() };
        let st2 = SortingReadSetContainer { ordered_sorting_keys: vec![key1.clone()], ordered_unsorted_keys: VecDeque::new(),aligned_read: fake_read.clone() };
        assert!(st1 > st2);

        let st1 = SortingReadSetContainer { ordered_sorting_keys: vec![key1.clone()], ordered_unsorted_keys: VecDeque::new(),aligned_read: fake_read.clone() };
        let st2 = SortingReadSetContainer { ordered_sorting_keys: vec![key1.clone()], ordered_unsorted_keys: VecDeque::new(),aligned_read: fake_read.clone() };
        assert!(!(st1 > st2) & !(st2 > st1));

        let st1 = SortingReadSetContainer { ordered_sorting_keys: vec![key1.clone(), key2.clone()], ordered_unsorted_keys: VecDeque::new(),aligned_read: fake_read.clone() };
        let st2 = SortingReadSetContainer { ordered_sorting_keys: vec![key1.clone(), key1.clone()], ordered_unsorted_keys: VecDeque::new(),aligned_read: fake_read.clone() };
        assert!(st1 > st2);

        // real example we hit
        let st1 = SortingReadSetContainer { ordered_sorting_keys: vec![('a',CorrectedKey::new('a',"AAACCCATCAGCATTA".as_bytes().to_vec(),"AAACCCATCAGCATTA".as_bytes().to_vec())),
                                                                       ('a',CorrectedKey::new('a',"TATTGACAACCT".as_bytes().to_vec(),"TATTGACAACCT".as_bytes().to_vec()))],
            ordered_unsorted_keys: VecDeque::new(),aligned_read: fake_read.clone() };
        let st2 = st1.clone();
        assert_eq!(st1.cmp(&st2) == Ordering::Equal, true);


    }
}