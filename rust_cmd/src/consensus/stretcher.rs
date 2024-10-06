use std::fmt;
use itertools::enumerate;
use crate::alignment::alignment_matrix::AlignmentResult;
use crate::alignment::fasta_bit_encoding::{encoding_to_u8, FASTA_A, FASTA_N, FastaBase};

// take a series of reads aligned to a common reference and create a consensus using their alignments
struct StretchConsensus {
    max_alignment_distance: f64,
}

struct AlignmentPile {
    //Vec<Alignment>
}

struct NucCounts {
    pub a: usize,
    pub c: usize,
    pub g: usize,
    pub t: usize,
    pub n: usize,
    pub gap: usize,
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
            .finish()
    }
}
impl NucCounts {
    pub fn new() -> NucCounts {
        NucCounts{
            a: 0,
            c: 0,
            g: 0,
            t: 0,
            n: 0,
            gap: 0,
        }
    }

    pub fn new_from(base: u8) -> NucCounts {
        let mut new_nc = NucCounts::new();
        new_nc.update(base);
        new_nc
    }

    pub fn update(&mut self, base: u8) {
        match base {
            b'a' | b'A' => {
                self.a += 1;
            },
            b'c' | b'C' => {
                self.c += 1;
            },
            b'g' | b'G' => {
                self.g += 1;
            },
            b't' | b'T' => {
                self.t += 1;
            },
            b'-' => {
                self.gap += 1;
            },
            _ => {
                self.n += 1;
            },
        }
    }

}

enum ReferenceStatus {
    Original{base: u8, original_position: usize, counts: NucCounts},
    Insertion{base: u8, counts: NucCounts},
}

struct AlignmentCandidate {
    reference: Vec<ReferenceStatus>,
    aligned_reads: Vec<NucCounts>,
    alignment_count: usize,
}

impl fmt::Debug for AlignmentCandidate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Point")
            .field("ref", &String::from_utf8(self.reference.clone()).unwrap())
            .field("alignments", &self.aligned_reads)
            .finish()
    }
}
impl AlignmentCandidate {
    pub fn new(alignment: &AlignmentResult) -> AlignmentCandidate {
        AlignmentCandidate{
            reference: FastaBase::string(&alignment.reference_aligned).into_bytes().map,
            aligned_reads: alignment.read_aligned.iter().map(|x| NucCounts::new_from(encoding_to_u8(x))).collect(),
            alignment_count: 1,
        }
    }

    pub fn add_alignment(&mut self, alignment: &AlignmentResult) {
        let mut existing_index = 0;
        let mut incoming_index = 0;

        // alignment.reference_aligned.iter().zip(alignment.read_aligned.iter()).collect::<Vec<(FastaBase,FastaBase)>>();


        while existing_index < alignment.reference_aligned.len() && incoming_index < alignment.read_aligned.len() {
            let incoming_ref_base = encoding_to_u8(&alignment.reference_aligned[incoming_index]);
            let incoming_read_base = encoding_to_u8(&alignment.read_aligned[incoming_index]);

            let existing_ref_base = &self.reference[existing_index];
            println!("Two mismatched reference nucleotides that are not gaps: {} and {}, pos {} and {}",*existing_ref_base as char,incoming_ref_base as char,existing_index,incoming_index);
            match (existing_ref_base, incoming_ref_base) {
                // this should not happen, where we're at the same position but we have different reference nucleotides; panic
                (existing_ref, new_ref)
                if *existing_ref != new_ref && *existing_ref != b'-' && new_ref != b'-' => {
                    panic!("Two mismatched reference nucleotides that are not gaps: {} and {}, pos {} and {}",*existing_ref as char,new_ref as char,existing_index,incoming_index);
                }
                // if the references agree it's easy -- we just add a count to the base storage corresponding to the read's nucleotide
                (existing_ref, new_ref)
                if *existing_ref == new_ref => {
                    self.aligned_reads.get_mut(existing_index).unwrap().update(incoming_read_base);
                    incoming_index += 1;
                    existing_index += 1; // advance the existing position,

                }
                // a gap in the existing reference and a base in the incoming reference -- simply skip over this site, adding no supporting bases
                (existing_ref, new_ref)
                if *existing_ref == b'-' => {
                    existing_index += 1; // advance the existing position,
                }
                // a gap in the incoming reference that's not in the existing reference -- make a new gap in the existing reference, and add a single dash of support (+1 to dash count)
                (existing_ref, new_ref)
                if incoming_ref_base == b'-' => {
                    // TODO: add a gap to the existing reference, and add a count for the incoming base
                    self.aligned_reads.insert(existing_index,NucCounts::new_from(incoming_read_base));
                    self.reference.insert(existing_index, b'-');
                    incoming_index += 1;
                    existing_index += 2; // advance the existing position over two places (the current base and the added base)

                }
                (x, y) => {
                    panic!("Unmanaged alignment merging issue for {} and {}",x,y);
                }
            }
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_alignment_result(read_bases: &str, ref_bases: &str) -> AlignmentResult {
        AlignmentResult{
            reference_name: "testref".to_string(),
            read_name: "testread".to_string(),
            reference_aligned: FastaBase::from_u8_slice(ref_bases.as_bytes()),
            read_aligned: FastaBase::from_u8_slice(read_bases.as_bytes()),
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
        let ref_bases =  "ACGTACGT";
        let read_bases = "ACG--CGT";
        let mut candidate = AlignmentCandidate::new(&create_alignment_result(read_bases,ref_bases));

        let ref_bases =  "ACGT-CGT";
        let read_bases = "ACGTACGT";
        //  ref_bases =        "ACGTACGT";

        candidate.add_alignment(&create_alignment_result(read_bases,ref_bases));
        println!("Alignment {:?}",candidate);

        // ref at this point     ACGT-ACGT
        let ref_bases =   "ACG--CGT";
        let read_bases =  "ACGTACGT";
        //  ref_bases =         "ACGTACGT";

        candidate.add_alignment(&create_alignment_result(read_bases,ref_bases));
        println!("Alignment {:?}",candidate);
    }
}