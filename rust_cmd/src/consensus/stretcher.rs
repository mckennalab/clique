use std::fmt;
use serde::{Deserialize, Serialize};
use alignment::alignment_matrix::AlignmentTag;
use alignment_functions::simplify_cigar_string;
use crate::alignment::alignment_matrix::AlignmentResult;
use crate::consensus::consensus_builders::{calculate_qual_scores, combine_qual_scores, prob_to_phred};

#[allow(dead_code)]
const DEFAULT_QUAL_FOR_UNKNOWN_QUAL: u8 = 32u8;

#[derive(Clone, Serialize, Deserialize, Hash)]
struct NucCounts {
    pub ref_base: u8,

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
        write!(f, "a: {} c {} g{} t{} n {} gap {}", self.a, self.c, self.g, self.t, self.n, self.gap)
    }
}

impl NucCounts {
    pub fn new(ref_base: &u8) -> NucCounts {
        NucCounts {
            ref_base: *ref_base,
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

    pub fn new_from(ref_base: &u8, base: u8, qual: u8) -> NucCounts {
        let mut new_nc = NucCounts::new(ref_base);
        new_nc.update(base, Some(qual));
        new_nc
    }

    pub fn proportion(&self, base: &u8, read_count: &usize) -> f64 {
        let cnt = match base {
            b'a' | b'A' => {
                self.a
            }
            b'c' | b'C' => {
                self.c
            }
            b'g' | b'G' => {
                self.g
            }
            b't' | b'T' => {
                self.t
            }
            b'-' => {
                self.gap
            }
            _ => {
                self.n
            }
        };
        //println!("Cnt {} {} {}",*base as char, cnt,self.total_count());
        cnt as f64 / *read_count as f64
    }


    pub fn update(&mut self, base: u8, qual: Option<u8>) {
        //println!("Base {}", base as char);
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
                self.n_qual.push(qual.unwrap());
            }
        }
    }

    pub fn total_count(&self) -> usize {
        self.a + self.c + self.g + self.t + self.n + self.gap
    }

    pub fn consensus_base(&self, gap_proportion_to_call: &f64) -> (u8, Option<u8>) {
        if (self.gap as f64) / (self.total_count() as f64) < *gap_proportion_to_call {
            let a_slice = vec![b'A'; self.a];
            let c_slice = vec![b'C'; self.c];
            let g_slice = vec![b'G'; self.g];
            let t_slice = vec![b'T'; self.t];
            let n_slice = vec![b'N'; self.n];
            println!("A {} C {} G {} T {}",self.a, self.c, self.g, self.t);

            let bases = vec![a_slice.as_slice(), c_slice.as_slice(), g_slice.as_slice(), t_slice.as_slice(), n_slice.as_slice()];

            let quals = vec![self.a_qual.as_slice(), self.c_qual.as_slice(), self.g_qual.as_slice(), self.t_qual.as_slice(), self.n_qual.as_slice()];

            let mut allele_props = combine_qual_scores(bases.as_slice(), quals.as_slice(), &self.ref_base, &0.75);
            let qual_normalized = calculate_qual_scores(&mut allele_props);
            println!("{:?}",qual_normalized);
            let index_of_max: usize = qual_normalized[0..4]
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.total_cmp(b))
                .map(|(index, _)| index).unwrap();

            let prob = prob_to_phred(&qual_normalized[index_of_max]);
            //println!("quals {:?} {:?} {:?} prop {} qual {}",allele_props,qual_normalized,quals,&qual_normalized[index_of_max],prob);
            match index_of_max {
                0 => (b'A', Some(prob)),
                1 => (b'C', Some(prob)),
                2 => (b'G', Some(prob)),
                3 => (b'T', Some(prob)),
                4 => (b'N', Some(prob)),
                _ => panic!("Unknown index")
            }
        } else {
            (b'-', None) // call a high-quality gap
        }
    }
}

/// The idea is that we need to keep track of bases that came from the reference originally and bases
/// that were inserted into the reference. Bases that are inserted in the reference only get matched
/// against other reads that have inserted bases at that position. Read deletions are captured within
/// bases in the reference sequence ('-' at that position)
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

pub struct AlignmentCandidate {
    reference: Vec<ReferenceStatus>,
    read_names: Vec<String>,
    reference_name: String,
}

impl fmt::Debug for AlignmentCandidate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let reference_string: Vec<u8> = self.reference.iter().map(|x| {
            match x {
                ReferenceStatus::Original { base, original_position: _, counts: _ } => { *base }
                ReferenceStatus::Insertion { base, counts: _ } => { *base }
            }
        }).collect();
        f.debug_struct("AlignmentCandidate")
            .field("ref string version", &String::from_utf8(reference_string).unwrap())
            .field("ref", &self.reference)
            .field("alignments", &self.read_names)
            .finish()
    }
}

impl AlignmentCandidate {
    pub fn new(reference: &[u8], reference_name: &[u8]) -> AlignmentCandidate {

        AlignmentCandidate {
            reference: reference.iter().enumerate().
                map(|(index, x)| {
                    //println!("adding ref base {:?}", *x as char);
                    ReferenceStatus::Original { base: *x, original_position: index, counts: NucCounts::new(x) }
                }).collect(),
            read_names: Vec::new(),
            reference_name: String::from_utf8(reference_name.to_vec()).unwrap(),

        }
    }

    // TODO this is broken: quals and gaps aren't handled right
    pub fn add_alignment(&mut self, alignment: &AlignmentResult) -> Result<(), String> {
        let mut existing_index = 0;
        let mut incoming_ref_index = 0;
        let mut incoming_read_qual_index = 0;


        let existing_ref_size = self.reference.len();
        self.read_names.push(alignment.read_name.clone());
        let read_qual = alignment.read_quals.clone().unwrap_or(alignment.read_aligned.iter().map(|_x: _| b'h').collect());
        //println!("*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*- entry with {} {} ", u8s(&alignment.reference_aligned), u8s(&alignment.read_aligned));

        while existing_index < existing_ref_size && incoming_ref_index < alignment.reference_aligned.len() {
            let incoming_ref_base = &alignment.reference_aligned[incoming_ref_index];
            let incoming_read_base = &alignment.read_aligned[incoming_ref_index];

            let incoming_read_qual = if incoming_read_base == &b'-' { &b'+' } else { &(read_qual[incoming_read_qual_index])};
            let existing_ref_base = &mut self.reference.get_mut(existing_index).unwrap();
            //println!("Loop: {:?} and {}, pos {} and {}, max is {} and {}", *existing_ref_base, *incoming_ref_base as char, existing_index, incoming_ref_index, existing_ref_size, alignment.read_aligned.len());

            match (existing_ref_base, incoming_ref_base) {
                // we're in an insertion on both references -- we're not going to concern ourselves with resolving multiple paths and just record insertion bases
                (ReferenceStatus::Insertion { base: _, counts  }, b'-') => {
                    counts.update(*incoming_read_base, Some(*incoming_read_qual));
                    //println!("step1 {} {:?}", *incoming_read_base as char, counts);

                    incoming_ref_index += 1;
                    existing_index += 1;
                }
                // an existing insertion to the reference but the new read isn't an insertion -- just skip over it
                (ReferenceStatus::Insertion { base: _, counts: _ }, _new_ref) => {
                    //println!("step2b");
                    // just move past it
                    existing_index += 1;
                }
                // we have a new insertion in the reference we haven't seen before
                (ReferenceStatus::Original { base: _, original_position: _, counts: _ }, b'-') => {
                    //println!("step2");
                    // we're going to add an insertion -- we have to choose if we're going to left or right align them (we choose right)
                    self.reference.insert(existing_index, ReferenceStatus::Insertion { base: *incoming_read_base, counts: NucCounts::new_from(&b'-', *incoming_read_base, *incoming_read_qual) });
                    incoming_ref_index += 1;
                    existing_index += 1;
                    if *incoming_read_base != b'-' {
                        incoming_read_qual_index += 1;
                    }
                }
                // we have a new insertion in the reference we haven't seen before

                // this should not happen, where we're at the same position but we have different reference nucleotides; panic
                (ReferenceStatus::Original { base, original_position: _, counts: _ }, new_ref) if base != new_ref && *base != b'-' && *new_ref != b'-' => {
                    return Err(format!("Two mismatched reference nucleotides that are not gaps: {} and {}, pos {} and {}", *base as char, *new_ref as char, existing_index, incoming_ref_index));
                }
                // easy -- two reference aligned bases -- make sure they agree and then add to the counts
                (ReferenceStatus::Original { base, original_position: _, counts }, new_ref) if base == new_ref && *base != b'-' && *new_ref != b'-' => {
                    //println!("step3 {}", *incoming_read_base as char);
                    counts.update(*incoming_read_base, Some(*incoming_read_qual));
                    incoming_ref_index += 1;
                    existing_index += 1;
                    if *incoming_read_base != b'-' {
                        incoming_read_qual_index += 1;
                    }
                }
                (x, y) => {
                    return Err(format!("Unmanaged alignment merging issue for {} and {}", x, y));
                }
            }
        };
        Ok(())
    }

    pub fn to_consensus(&mut self, gap_call_threshold: &f64) -> AlignmentResult {
        assert!(self.read_names.len() > 0);
        let mut resulting_alignmented_read = Vec::new();
        let mut resulting_alignmented_ref = Vec::new();
        let mut resulting_alignmented_qual = Vec::new();
        let mut cigar_tokens = Vec::new();

        self.reference.iter().for_each(|rb| {
            match rb {
                ReferenceStatus::Original { base, original_position: _, counts } => {
                    //println!("original {} {:?}", *base as char, counts);
                    let base_and_qual = counts.consensus_base(gap_call_threshold);
                    resulting_alignmented_ref.push(*base);
                    resulting_alignmented_read.push(base_and_qual.0);

                    match base_and_qual.0 {
                        b'-' => {cigar_tokens.push(AlignmentTag::Del(1))}
                        _ => {
                            resulting_alignmented_qual.push(base_and_qual.1.unwrap() + 33);
                            cigar_tokens.push(AlignmentTag::MatchMismatch(1));
                        }
                    }

                }
                ReferenceStatus::Insertion { base, counts } if counts.proportion(base,&self.read_names.len()) >= *gap_call_threshold => {
                    //println!("added insert {} {}", *base as char, counts.proportion(base, &self.read_names.len()));
                    let base_qual = counts.consensus_base(gap_call_threshold);
                    resulting_alignmented_ref.push(b'-');
                    resulting_alignmented_read.push(base_qual.0);

                    match base_qual.0 {
                        b'-' => {panic!("Can't insert a deletion")}
                        _ => {
                            cigar_tokens.push(AlignmentTag::Ins(1));
                            resulting_alignmented_qual.push(base_qual.1.unwrap() + 33);
                        }
                    }
                }
                ReferenceStatus::Insertion { base: _, counts: _ } => {
                    //println!("dropped insert {} {}", *base as char, counts.proportion(base,&self.read_names.len()));
                    // do nothing, we're not going to include this gap in the reference as it's not supported by enough reads
                }
            }
        });
        //if resulting_alignmented_read.len() != resulting_alignmented_qual.len() {
        //    println!("ref {} read seq {} quals {}-----",self.reference_name.clone(),u8s(&resulting_alignmented_read.clone()), u8s(&resulting_alignmented_qual.clone()));
        //}
        //assert_eq!(resulting_alignmented_read.len(),resulting_alignmented_qual.len());
        //

        AlignmentResult {
            reference_name: self.reference_name.clone(),
            read_name: self.read_names.get(0).unwrap_or(&"UnnamedRead".to_string()).clone(),
            reference_aligned: resulting_alignmented_ref,
            read_aligned: resulting_alignmented_read,
            read_quals: Some(resulting_alignmented_qual),
            cigar_string: simplify_cigar_string(&cigar_tokens),
            path: vec![],
            score: 0.0,
            reference_start: 0,
            read_start: 0,
            bounding_box: None,
        }
    }
}
// AAAGAACAGGCCTTGTTGAAAGAGACTAGATCGATAGCGATCCCGGCCGATAGGGGATGTACCTGTCGTCTTAGCTAAGATGACAGGACATGTCCAGGAAGTGCTCGTGTACTTCCTGGCCCATGTACTCTGCGTTGATACCACTGCTT------------------------------------------------------------------
// ###########################################################################################################"#########################################
#[cfg(test)]
mod tests {
    use utils::read_utils::u8s;
    use super::*;

    fn create_alignment_result(read_bases: &str, ref_bases: &str) -> AlignmentResult {
        AlignmentResult {
            reference_name: "testref".to_string(),
            read_name: "testread".to_string(),
            reference_aligned: ref_bases.as_bytes().to_vec(),
            read_aligned: read_bases.as_bytes().to_vec(),
            read_quals: Some(read_bases.as_bytes().iter().filter(|x| **x != b'-').map(|_x| DEFAULT_QUAL_FOR_UNKNOWN_QUAL).collect::<Vec<u8>>()), // not right
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
        let mut candidate = AlignmentCandidate::new(ref_bases.as_bytes(),"ref_name".as_bytes());
        candidate.add_alignment(&create_alignment_result(read_bases, ref_bases)).expect("Failed to add candidate to alignment");
        let conc = candidate.to_consensus(&0.75);
        assert_eq!(u8s(&conc.reference_aligned),u8s(&ref_bases.as_bytes().to_vec()));
        assert_eq!(u8s(&conc.read_aligned),u8s(&read_bases.as_bytes().to_vec()));

        let ref_bases = "ACGT-ACGT";
        let read_bases = "ACGTAACGT";
        candidate.add_alignment(&create_alignment_result(read_bases, ref_bases)).expect("Failed to add candidate to alignment");
        let conc = candidate.to_consensus(&0.75);
        assert_eq!(u8s(&conc.reference_aligned),u8s(&"ACGTACGT".as_bytes().to_vec())); // we don't have enough evidence for the insertion
        assert_eq!(u8s(&conc.read_aligned),u8s(&"ACGTACGT".as_bytes().to_vec())); // we don't have enough evidence for the insertion

        let ref_bases = "ACGTACGT";
        let read_bases = "ACGTACGT";
        candidate.add_alignment(&create_alignment_result(read_bases, ref_bases)).expect("Failed to add candidate to alignment");
        let conc = candidate.to_consensus(&0.75);
        assert_eq!(u8s(&conc.reference_aligned),u8s(&"ACGTACGT".as_bytes().to_vec()));
        assert_eq!(u8s(&conc.read_aligned),u8s(&"ACGTACGT".as_bytes().to_vec()));


        let ref_bases = "ACGTACGT";
        let read_bases = "--------";
        candidate.add_alignment(&create_alignment_result(read_bases, ref_bases)).expect("Failed to add candidate to alignment");
        let conc = candidate.to_consensus(&0.75);
        assert_eq!(u8s(&conc.reference_aligned),u8s(&"ACGTACGT".as_bytes().to_vec()));
        assert_eq!(u8s(&conc.read_aligned),u8s(&"ACGTACGT".as_bytes().to_vec()));

        for _ in 0..20 {
            //     t ref_bases =   "ACGTACGT";
            let ref_bases = "ACGT----ACGT";
            let read_bases = "ACGTAGGAACGT";

            candidate.add_alignment(&create_alignment_result(read_bases, ref_bases)).expect("Unable to unwrap candidate");
        }
        let conc = candidate.to_consensus(&0.75);
        assert_eq!(u8s(&conc.reference_aligned),u8s(&"ACGT----ACGT".as_bytes().to_vec()));
        assert_eq!(u8s(&conc.read_aligned),u8s(&"ACGTAGGAACGT".as_bytes().to_vec()));
    }

}