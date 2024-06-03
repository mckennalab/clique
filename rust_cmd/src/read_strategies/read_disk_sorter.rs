use std::cmp::Ordering;
use std::collections::{VecDeque};
use serde::{Serialize, Deserialize};
use crate::alignment::alignment_matrix::{AlignmentResult};
use crate::alignment::fasta_bit_encoding::{FastaBase};

/// a sortable read set container that sorts on a set of keys -- which we populate with
/// extracted barcode sequences. These sorting sequences could have been corrected to a known list
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SortingReadSetContainer {
    pub ordered_sorting_keys: Vec<(char,Vec<FastaBase>)>,   // grow in order
    pub ordered_unsorted_keys: VecDeque<(char, Vec<FastaBase>)>, // use the default behavior to push back, pop front
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
    use crate::alignment::alignment_matrix::AlignmentResult;
    use crate::alignment::fasta_bit_encoding::{FASTA_A, FASTA_N, FASTA_T, FastaBase};
    use crate::read_strategies::read_disk_sorter::{SortingReadSetContainer};
    use crate::utils::read_utils::fake_reads;

    #[test]
    fn test_ordinal_nature() {
        let srsc1 = SortingReadSetContainer{
            ordered_sorting_keys: vec![('*', vec![FASTA_A, FASTA_A])],
            ordered_unsorted_keys: Default::default(),
            aligned_read: AlignmentResult {
                reference_name: "".to_string(),
                read_name: "".to_string(),
                reference_aligned: vec![],
                read_aligned: vec![],
                cigar_string: vec![],
                path: vec![],
                score: 0.0,
                reference_start: 0,
                read_start: 0,
                bounding_box: None,
            }
        };
        let srsc2 = SortingReadSetContainer{
            ordered_sorting_keys: vec![('*', vec![FASTA_A, FASTA_A])],
            ordered_unsorted_keys: Default::default(),
            aligned_read: AlignmentResult {
                reference_name: "".to_string(),
                read_name: "".to_string(),
                reference_aligned: vec![],
                read_aligned: vec![],
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
            ordered_sorting_keys: vec![('*', vec![FASTA_A, FASTA_A]),('*', vec![FASTA_A, FASTA_A])],
            ordered_unsorted_keys: Default::default(),
            aligned_read: AlignmentResult {
                reference_name: "".to_string(),
                read_name: "".to_string(),
                reference_aligned: vec![],
                read_aligned: vec![],
                cigar_string: vec![],
                path: vec![],
                score: 0.0,
                reference_start: 0,
                read_start: 0,
                bounding_box: None,
            }
        };
        let srsc2 = SortingReadSetContainer{
            ordered_sorting_keys: vec![('*', vec![FASTA_A, FASTA_A]),('*', vec![FASTA_A, FASTA_T])],
            ordered_unsorted_keys: Default::default(),
            aligned_read: AlignmentResult {
                reference_name: "".to_string(),
                read_name: "".to_string(),
                reference_aligned: vec![],
                read_aligned: vec![],
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
            ordered_sorting_keys: vec![('*', vec![FASTA_A, FASTA_A]),
                                       ('*', vec![FASTA_A, FASTA_A]),
                                       ('*', vec![FASTA_A, FASTA_A])],
            ordered_unsorted_keys: Default::default(),
            aligned_read: AlignmentResult {
                reference_name: "".to_string(),
                read_name: "".to_string(),
                reference_aligned: vec![],
                read_aligned: vec![],
                cigar_string: vec![],
                path: vec![],
                score: 0.0,
                reference_start: 0,
                read_start: 0,
                bounding_box: None,
            }
        };
        let srsc2 = SortingReadSetContainer{
            ordered_sorting_keys: vec![('*', vec![FASTA_A, FASTA_A]),
                                       ('*', vec![FASTA_A, FASTA_T]),
                                       ('*', vec![FASTA_A, FASTA_A])],
            ordered_unsorted_keys: Default::default(),
            aligned_read: AlignmentResult {
                reference_name: "".to_string(),
                read_name: "".to_string(),
                reference_aligned: vec![],
                read_aligned: vec![],
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
        let key1 = ('$',vec![FASTA_N, FASTA_A]);
        let key2 = ('$',vec![FASTA_N, FASTA_N]);

        let reads = fake_reads(10, 1);
        let read_seq = reads.get(0).unwrap().read_one.seq().iter().map(|x| FastaBase::from(x.clone())).collect::<Vec<FastaBase>>();
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
        let st1 = SortingReadSetContainer { ordered_sorting_keys: vec![('a',FastaBase::from_str("AAACCCATCAGCATTA")),
                                                                        ('a',FastaBase::from_str("TATTGACAACCT"))],
            ordered_unsorted_keys: VecDeque::new(),aligned_read: fake_read.clone() };
        let st2 = st1.clone();
        assert_eq!(st1.cmp(&st2) == Ordering::Equal, true);


    }
}