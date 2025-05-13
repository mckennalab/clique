use std::fmt;
use itertools::enumerate;
use serde::{Deserialize, Serialize};
use crate::alignment::alignment_matrix::AlignmentResult;
use crate::consensus::consensus_builders::{calculate_qual_scores, combine_qual_scores, prob_to_phred};

// take a series of reads aligned to a common reference and create a consensus using their alignments
struct StretchConsensus {
    max_alignment_distance: f64,
}

struct AlignmentPile {
    //Vec<Alignment>
}


#[derive(Clone, Serialize, Deserialize, Hash)]
struct NucCounts {
    pub a: usize,
    pub c: usize,
    pub g: usize,
    pub t: usize,
    pub n: usize,
    pub gap: usize,

    pub a_qual: Vec<u8>,
    pub c_qual: Vec<u8>,
    pub g_qual: Vec<u8>,
    pub t_qual: Vec<u8>,
    pub n_qual: Vec<u8>,
}

impl fmt::Debug for NucCounts {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("NucCount")
            .field("a", &self.a)
            .field("c", &self.c)
            .field("g", &self.g)
            .field("t", &self.t)
            .field("n", &self.n)
            .field("gap", &self.gap)
            .field("a_qual", &self.a)
            .field("c_qual", &self.c)
            .field("g_qual", &self.g)
            .field("t_qual", &self.t)
            .field("n_qual", &self.n)
            .finish()
    }
}

impl fmt::Display for NucCounts {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "a: {} c {} g{} t{} n {} gap {}",self.a,self.c,self.g,self.t,self.n,self.gap)
    }
}
impl NucCounts {
    pub fn new() -> NucCounts {
        NucCounts {
            a: 0,
            c: 0,
            g: 0,
            t: 0,
            n: 0,
            gap: 0,
            a_qual: Vec::new(),
            c_qual: Vec::new(),
            g_qual: Vec::new(),
            t_qual: Vec::new(),
            n_qual: Vec::new(),
        }
    }

    pub fn new_from(base: u8, qual: u8) -> NucCounts {
        let mut new_nc = NucCounts::new();
        new_nc.update(base,Some(qual));
        new_nc
    }

    pub fn update(&mut self, base: u8, qual: Option<u8>) {
        match base {
            b'a' | b'A' => {
                self.a += 1;
                self.a_qual.push(qual.unwrap());
            }
            b'c' | b'C' => {
                self.c += 1;
                self.c_qual.push(qual.unwrap());
            }
            b'g' | b'G' => {
                self.g += 1;
                self.g_qual.push(qual.unwrap());
            }
            b't' | b'T' => {
                self.t += 1;
                self.t_qual.push(qual.unwrap());
            }
            b'-' => {
                self.gap += 1;
            }
            _ => {
                self.n += 1;
                self.a_qual.push(qual.unwrap());
            }
        }
    }

    pub fn total_count(&self) -> usize {
        self.a + self.c + self.g + self.t + self.n
    }

    pub fn consensus_base(&self, gap_proportion_to_call: &f64) -> (u8,u8) {

        if (self.gap as f64) / (self.total_count() as f64) >= *gap_proportion_to_call {
            let mut bases = vec![b'A'; self.a];
            bases.extend(vec![b'C'; self.c]);
            bases.extend(vec![b'G'; self.g]);
            bases.extend(vec![b'T'; self.t]);
            bases.extend(vec![b'N'; self.n]);

            let mut quals = self.a_qual.clone();
            quals.extend(self.c_qual.clone());
            quals.extend(self.g_qual.clone());
            quals.extend(self.t_qual.clone());
            quals.extend(self.n_qual.clone());

            let mut allele_props = combine_qual_scores(&bases, &quals, &0.75, &true);
            let qual_normalized = calculate_qual_scores(&mut allele_props);
            let index_of_max: usize = qual_normalized
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.total_cmp(b))
                .map(|(index, _)| index).unwrap();

            let prob = prob_to_phred(&qual_normalized[index_of_max]);

            match index_of_max {
                0 => (b'A', prob),
                1 => (b'A', prob),
                2 => (b'A', prob),
                3 => (b'A', prob),
                4 => (b'A', prob),
                _ => panic!("Unknown index")
            }
        } else {
            (b'-', b'H')
        }
    }
}

/// The idea is that we need to keep track of bases that came from the reference originally and bases
/// that were inserted into the reference. Bases that are inserted in the reference only get matched
/// against other reads that have inserted bases at that position
#[derive(Clone, Serialize, Deserialize, Hash)]
enum ReferenceStatus {
    Original { base: u8, original_position: usize, counts: NucCounts },
    Insertion { base: u8, counts: NucCounts },
}

impl PartialEq<u8> for ReferenceStatus {
    fn eq(&self, other: &u8) -> bool {
        match self {
            ReferenceStatus::Original { base, .. } => base == other,
            ReferenceStatus::Insertion { base, .. } => base == other,
        }
    }
}
impl fmt::Debug for ReferenceStatus {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ReferenceStatus::Original {
                base,
                original_position,
                counts,
            } => f
                .debug_struct("Original")
                .field("base", &format_args!("{:?}", *base as char))
                .field("original_position", original_position)
                .field("counts", counts)
                .finish(),
            ReferenceStatus::Insertion { base, counts } => f
                .debug_struct("Insertion")
                .field("base", &format_args!("{:?}", *base as char))
                .field("counts", counts)
                .finish(),
        }
    }
}


impl fmt::Display for ReferenceStatus {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ReferenceStatus::Original { base, original_position, counts } => {
                write!(
                    f,
                    "Original: base = {}, position = {}, counts = {}",
                    *base as char, original_position, counts
                )
            }
            ReferenceStatus::Insertion { base, counts } => {
                write!(f, "Insertion: base = {}, counts = {}", *base as char, counts)
            }
        }
    }
}
struct AlignmentCandidate {
    reference: Vec<ReferenceStatus>,
    read_names : Vec<String>,
    reference_name: Option<String>,
}

impl fmt::Debug for AlignmentCandidate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("AlignmentCandidate")
            .field("ref", &self.reference)
            .field("alignments", &self.read_names)
            .finish()
    }
}

impl AlignmentCandidate {
    pub fn new(alignment: &AlignmentResult) -> AlignmentCandidate {
        AlignmentCandidate {
            reference: alignment.reference_aligned.
                iter().enumerate().
                map(|(index,x)| {
                    println!("adding ref base {:?}",*x as char);
                    ReferenceStatus::Original {base: *x, original_position: index, counts: NucCounts::new()}
                }).collect(),
            read_names: Vec::new(),
            reference_name: Some(alignment.reference_name.clone()),

        }
    }

    pub fn add_alignment(&mut self, alignment: &AlignmentResult) {
        let mut existing_index = 0;
        let mut incoming_index = 0;
        let existing_ref_size =  self.reference.len();
        self.read_names.push(alignment.read_name.clone());
        let read_qual = alignment.read_quals.clone().unwrap_or(alignment.read_aligned.iter().map(|x| b'h').collect());


        while existing_index < existing_ref_size && incoming_index < alignment.reference_aligned.len() {
            let incoming_ref_base = &alignment.reference_aligned[incoming_index];
            let incoming_read_base = &alignment.read_aligned[incoming_index];
            let incoming_read_qual = read_qual[incoming_index];

            let existing_ref_base = &mut self.reference.get_mut(existing_index).unwrap();
            //println!("Loop: {:?} and {}, pos {} and {}, max is {} and {}", *existing_ref_base, incoming_ref_base as char, existing_index, incoming_index, existing_ref_size, alignment.read_aligned.len());

            match (existing_ref_base, incoming_ref_base) {
                // we're in an insertion on both references -- we're not going to concern ourselves with multiple paths and just record insertion bases
                (ReferenceStatus::Insertion { base, counts }, b'-') => {
                    //println!("step1");
                    counts.update(*incoming_read_base, Some(incoming_read_qual));
                    incoming_index += 1;
                    existing_index += 1;
                }
                // an existing insertion to the reference but the new read isn't an insertion -- just skip over it
                (ReferenceStatus::Insertion { base, counts}, new_ref) => {
                    //println!("step2b");
                    // we're going to add an insertion -- we have to choose if we're going to left or right align them (we choose right)
                    existing_index += 1;
                }
                // we have a new insertion in the reference we haven't seen before
                (ReferenceStatus::Original { base, original_position, counts}, b'-') => {
                    //println!("step2");
                    // we're going to add an insertion -- we have to choose if we're going to left or right align them (we choose right)
                    self.reference.insert(existing_index, ReferenceStatus::Insertion{base: *incoming_read_base, counts: NucCounts::new_from(*incoming_read_base, incoming_read_qual) });
                    incoming_index += 1;
                    existing_index += 1;
                }
                // we have a new insertion in the reference we haven't seen before

                // this should not happen, where we're at the same position but we have different reference nucleotides; panic
                (ReferenceStatus::Original { base, original_position, counts}, new_ref)
                if base != new_ref && *base != b'-' && *new_ref != b'-' => {
                    panic!("Two mismatched reference nucleotides that are not gaps: {} and {}, pos {} and {}", *base as char, *new_ref as char, existing_index, incoming_index);
                }
                // easy -- two reference aligned bases -- make sure they agree and then add to the counts
                (ReferenceStatus::Original { base, original_position, counts}, new_ref)
                if base == new_ref && *base != b'-' && *new_ref != b'-' => {
                    //println!("step3 {}",incoming_read_base as char);
                    counts.update(*incoming_read_base, Some(incoming_read_qual));;
                    incoming_index += 1;
                    existing_index += 1;
                }
                (x, y) => {
                    panic!("Unmanaged alignment merging issue for {} and {}", x, y);
                }
            }
        };
    }

    pub fn to_consensus(&mut self, gap_call_threshold: &f64) -> AlignmentResult {
        assert!(self.read_names.len() > 0);
        let mut resulting_alignmented_read = Vec::new();
        let mut resulting_alignmented_ref = Vec::new();
        let mut resulting_alignmented_qual = Vec::new();
        let cigar_tokens = Vec::new();
        panic!("implement cigars");
        self.reference.iter().for_each(|rb| {
            match rb {
                ReferenceStatus::Original { base, original_position, counts} => {
                    let base_qual = counts.consensus_base(gap_call_threshold);
                    resulting_alignmented_ref.push(*base);
                    resulting_alignmented_read.push(base_qual.0);
                    resulting_alignmented_qual.push(base_qual.1);
                }
                ReferenceStatus::Insertion { base, counts} => {
                    let base_qual = counts.consensus_base(gap_call_threshold);
                    resulting_alignmented_ref.push(*base);
                    resulting_alignmented_read.push(base_qual.0);
                    resulting_alignmented_qual.push(base_qual.1);
                }
            }
        });

        AlignmentResult{
            reference_name: self.reference_name.as_ref().unwrap().clone(),
            read_name: self.read_names.get(0).unwrap_or(&"UnnamedRead".to_string()).clone(),
            reference_aligned: resulting_alignmented_ref.clone(),
            read_aligned: resulting_alignmented_read.clone(),
            read_quals: Some(resulting_alignmented_qual),
            cigar_string: cigar_tokens,
            path: vec![],
            score: 0.0,
            reference_start: 0,
            read_start: 0,
            bounding_box: None,
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    fn create_alignment_result(read_bases: &str, ref_bases: &str) -> AlignmentResult {
        AlignmentResult {
            reference_name: "testref".to_string(),
            read_name: "testread".to_string(),
            reference_aligned: ref_bases.as_bytes().to_vec(),
            read_aligned: read_bases.as_bytes().to_vec(),
            read_quals: Some(read_bases.as_bytes().iter().filter(|x| **x != b'-').map(|x| b'H').collect::<Vec<u8>>()), // not right
            cigar_string: vec![],
            path: vec![],
            score: 0.0,
            reference_start: 0,
            read_start: 0,
            bounding_box: None,
        }
    }

    #[test]
    fn test_merge_two_references() {
        let ref_bases = "ACGTACGT";
        let read_bases = "ACG--CGT";
        let mut candidate = AlignmentCandidate::new(&create_alignment_result(read_bases, ref_bases));

        //     t ref_bases =   "ACGTACGT";
        let ref_bases =  "ACGT-ACGT";
        let read_bases = "ACGTAACGT";

        candidate.add_alignment(&create_alignment_result(read_bases, ref_bases));
        println!("Alignment {:?}", candidate);

        let ref_bases =  "ACGTACGT";
        let read_bases = "ACGTACGT";

        candidate.add_alignment(&create_alignment_result(read_bases, ref_bases));

        let ref_bases =  "ACGTACGT";
        let read_bases = "--------";
        println!("Alignment {:?}", candidate);
        candidate.add_alignment(&create_alignment_result(read_bases, ref_bases));

        println!("Alignment {:?}", candidate);
    }
}