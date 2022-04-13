use std::collections::{HashMap, BTreeMap};
use std::fmt::Write; // needed by the `write!` macro

use bio::alignment::{Alignment, AlignmentMode, AlignmentOperation, TextSlice};
use bio::alignment::pairwise::{Aligner, MIN_SCORE, Scoring};

// Two sets of known characters: our standard DNA alphabet and a second version with known gaps.
// These are used to mask known values when looking for extractable UMI/ID/barcode sequences
lazy_static! {
    static ref KNOWNBASES: HashMap<u8, char> = {
        let mut hashedvalues = HashMap::new();
        hashedvalues.insert(b'a', 'A');
        hashedvalues.insert(b'A', 'A');
        hashedvalues.insert(b'c', 'C');
        hashedvalues.insert(b'C', 'C');
        hashedvalues.insert(b'g', 'G');
        hashedvalues.insert(b'G', 'G');
        hashedvalues.insert(b't', 'T');
        hashedvalues.insert(b'T', 'T');
        hashedvalues
    };

    static ref KNOWNBASESPLUSINSERT: HashMap<u8, u8> = {
        let mut hashedvalues = HashMap::new();
        hashedvalues.insert(b'a', b'A');
        hashedvalues.insert(b'A', b'A');
        hashedvalues.insert(b'c', b'C');
        hashedvalues.insert(b'C', b'C');
        hashedvalues.insert(b'g', b'G');
        hashedvalues.insert(b'G', b'G');
        hashedvalues.insert(b't', b'T');
        hashedvalues.insert(b'T', b'T');
        hashedvalues.insert(b'-', b'-');
        hashedvalues
    };
}

lazy_static! {
    static ref REVERSECOMP: HashMap<u8, u8> = {
            let mut hashedvalues: HashMap<u8,u8> = HashMap::new();
            hashedvalues.insert(b'a', b'T');
            hashedvalues.insert(b'A', b'T');
            hashedvalues.insert(b'c', b'G');
            hashedvalues.insert(b'C', b'G');
            hashedvalues.insert(b'g', b'C');
            hashedvalues.insert(b'G', b'C');
            hashedvalues.insert(b't', b'A');
            hashedvalues.insert(b'T', b'A');
            hashedvalues
        };
}

// reverse complement a string of DNA nucleotides, with non canonical (ACGT) bases converted to N
// when reverse complemented
pub fn reverse_complement_string(sequence: &String) -> String {
    let mut rev_comp: String = String::with_capacity(sequence.len());
    // iterate through the input &str
    let mut bytes = sequence.as_bytes().to_owned();
    bytes.reverse();
    for c in bytes {
        // test the input
        match c {
            x if REVERSECOMP.contains_key(&c) => write!(rev_comp, "{:02x}", REVERSECOMP[&x]).unwrap(),
            _ => rev_comp.push('N')
        }
    }
    rev_comp
}


// reverse complement a string of DNA nucleotides, with non canonical (ACGT) bases converted to N
// when reverse complemented
pub fn reverse_complement_u8(sequence: &[u8]) -> String {
    let mut rev_comp: String = String::with_capacity(sequence.len());
    // iterate through the input &str

    let mut chars = sequence.to_owned();
    chars.reverse();

    for c in chars {
        // test the input
        match c {
            x if REVERSECOMP.contains_key(&c) => rev_comp.push(REVERSECOMP[&x] as char),
            _ => rev_comp.push('N')
        }
    }
    rev_comp
}

// return the aligned strings for a forward read and the reference
pub fn align_unknown_orientation_read(read: &String, reference: &String) -> (Alignment, String, String) {

    let fwd = align_forward_read(read,reference);
    let rev = align_forward_read(&reverse_complement_string(read), reference);

    if rev.0.score > fwd.0.score {
        rev
    } else {
        fwd
    }
}


// return the aligned strings for a forward read and the reference
pub fn align_unknown_orientation_read_u8_ref(read: &[u8], reference: &[u8]) -> (Alignment, String, String) {

    let fwd = align_forward_read_u8(read,reference);
    let rev = align_forward_read_u8(&reverse_complement_u8(read).as_bytes(), reference);

    if rev.0.score > fwd.0.score {
        rev
    } else {
        fwd
    }
}

// reverse complement a read before aligning it to the reference
pub fn align_reverse_read(read2: &String, reference: &String) -> (Alignment, String, String) {
    align_forward_read(&reverse_complement_string(read2), reference)
}

// return the aligned strings for a forward read and the reference
pub fn align_forward_read(read1: &String, reference: &String) -> (Alignment, String, String) {
    let ref_bytes = reference.as_bytes();
    let read_bytes = read1.as_bytes();

    align_forward_read_u8(read_bytes,ref_bytes)
}

pub fn align_forward_read_u8(read1: &[u8], reference: &[u8]) -> (Alignment, String, String) {

        let scoring = Scoring::new(-20, -1, &custom_umi_score) // Gap open, gap extend and our custom match score function
        .xclip(MIN_SCORE) // Clipping penalty for x set to 'negative infinity', hence global in x
        .yclip(0); // Clipping penalty for y set to 0, hence local in y


    let mut aligner = Aligner::with_scoring(scoring);

    let alignment = aligner.custom(reference, read1); // The custom aligner invocation
    let alignment_string = alignment_strings(&alignment, reference, read1);

    (alignment, alignment_string.0.replace(" ", "-"), alignment_string.1.replace(" ", "-"))
}

// return a list of
pub fn extract_tagged_sequences(aligned_read: &String, aligned_ref: &String) -> BTreeMap<String, String> {
    let mut special_values: BTreeMap<u8, Vec<u8>> = BTreeMap::new();
    let empty = &Vec::new(); // ugh this is dumb
    for (reference_base, read_base) in std::iter::zip(aligned_ref.as_bytes(), aligned_read.as_bytes()) {
        if !KNOWNBASESPLUSINSERT.contains_key(&reference_base) {
            let mut current_code = special_values.get(&reference_base).unwrap_or(empty).clone();
            current_code.push(read_base.clone());

            special_values.insert(*reference_base, current_code.clone());
        } else if reference_base.is_ascii_uppercase() {
            let mut current_code = special_values.get(&('r' as u8)).unwrap_or(empty).clone();
            current_code.push(read_base.clone());

            special_values.insert('r' as u8, current_code.clone());
        }
    }
    special_values.iter().map(|(key, value)| {
        (String::from_utf8(vec![*key]).unwrap(), String::from_utf8(value.clone()).unwrap())
    }).collect()
}


// A mangled version of the Rust::bio alignment code, to just return the strings from the two alignments
pub fn alignment_strings(alignment: &Alignment, x: TextSlice, y: TextSlice) -> (String, String) {
    let mut x_pretty = String::new();
    let mut y_pretty = String::new();
    let mut inb_pretty = String::new();

    if !alignment.operations.is_empty() {
        let mut x_i: usize;
        let mut y_i: usize;

        // If the alignment mode is one of the standard ones, the prefix clipping is
        // implicit so we need to process it here
        match alignment.mode {
            AlignmentMode::Custom => {
                x_i = 0;
                y_i = 0;
            }
            _ => {
                x_i = alignment.xstart;
                y_i = alignment.ystart;
                for k in x.iter().take(alignment.xstart) {
                    x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[*k])));
                    inb_pretty.push(' ');
                    y_pretty.push(' ')
                }
                for k in y.iter().take(alignment.ystart) {
                    y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[*k])));
                    inb_pretty.push(' ');
                    x_pretty.push(' ')
                }
            }
        }

        // Process the alignment.
        for i in 0..alignment.operations.len() {
            match alignment.operations[i] {
                AlignmentOperation::Match => {
                    x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[x[x_i]])));
                    x_i += 1;

                    inb_pretty.push('|');

                    y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[y[y_i]])));
                    y_i += 1;
                }
                AlignmentOperation::Subst => {
                    x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[x[x_i]])));
                    x_i += 1;

                    inb_pretty.push('\\');

                    y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[y[y_i]])));
                    y_i += 1;
                }
                AlignmentOperation::Del => {
                    x_pretty.push('-');

                    inb_pretty.push('x');

                    y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[y[y_i]])));
                    y_i += 1;
                }
                AlignmentOperation::Ins => {
                    x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[x[x_i]])));
                    x_i += 1;

                    inb_pretty.push('+');

                    y_pretty.push('-');
                }
                AlignmentOperation::Xclip(len) => {
                    for k in x.iter().take(len) {
                        x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[*k])));
                        x_i += 1;

                        inb_pretty.push(' ');

                        y_pretty.push(' ')
                    }
                }
                AlignmentOperation::Yclip(len) => {
                    for k in y.iter().take(len) {
                        y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[*k])));
                        y_i += 1;

                        inb_pretty.push(' ');

                        x_pretty.push(' ')
                    }
                }
            }
        }

        // If the alignment mode is one of the standard ones, the suffix clipping is
        // implicit so we need to process it here
        match alignment.mode {
            AlignmentMode::Custom => {}
            _ => {
                for k in x.iter().take(alignment.xlen).skip(x_i) {
                    x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[*k])));
                    inb_pretty.push(' ');
                    y_pretty.push(' ')
                }
                for k in y.iter().take(alignment.ylen).skip(y_i) {
                    y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[*k])));
                    inb_pretty.push(' ');
                    x_pretty.push(' ')
                }
            }
        }
    }

    assert_eq!(x_pretty.len(), inb_pretty.len());
    assert_eq!(y_pretty.len(), inb_pretty.len());

    (x_pretty, y_pretty)
}

// a custom scoring function for matching nucleotide bases to each other or
// zero-cost matches to other characters. This allows us to encode
pub fn custom_umi_score(a: u8, b: u8) -> i32 {
    match (a, b) {
        (a, b) if KNOWNBASES.contains_key(&a) && KNOWNBASES.contains_key(&b) && KNOWNBASES[&a] == KNOWNBASES[&b] => { 10 }
        (a, b) if KNOWNBASES.contains_key(&a) && KNOWNBASES.contains_key(&b) => { -8 }
        _ => { 0 }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alignment_basic_test() {
        let reference = String::from("AATGATACGGCGACCACCGAGATCTACAC##########ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNNN##########CTGTAGGTAGTTTGTC");
        let test_read = String::from("AATGATACGGCGACCAGATCTACACACCCCTTTGCACACTCTTTCCCTACACGACGCTCTTCCGATCTAAAAAAAATTTTTTTTTTCTGTAGGTAGTTTGTC");
        let test_read_aligned = String::from("AATGATACGGCGACC----AGATCTACACACCCCTTTGCACACTCTTTCCCTACACGACGCTCTTCCGATCTAAAAAAAATTTTTTTTTTCTGTAGGTAGTTTGTC");

        let aligned_string = align_forward_read(&test_read, &reference);
        assert_eq!(aligned_string.1, reference);
        assert_eq!(aligned_string.2, test_read_aligned);
    }
    //

    #[test]
    fn tagged_sequence_test() {
        let reference = String::from("AATGATACGGCGACCACCGAGATCTACAC##########ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNNN##########CTGTAGGTAGTTTGTC");
        let test_read = String::from("AATGATACGGCGACCAGATCTACACACCCCTTTGCACACTCTTTCCCTACACGACGCTCTTCCGATCTAAAAAAAATTTTTTTTTTCTGTAGGTAGTTTGTC");
        let test_read_aligned = String::from("AATGATACGGCGACC----AGATCTACACACCCCTTTGCACACTCTTTCCCTACACGACGCTCTTCCGATCTAAAAAAAATTTTTTTTTTCTGTAGGTAGTTTGTC");

        let keyvalues = extract_tagged_sequences(&test_read, &reference);
        for (key, value) in keyvalues {
            println!("{} / {}", key, value);
        }
    }

    #[test]
    fn tagged_sequence_test_space() {
        for n in 1..100 {
            let reference = String::from("AAATACTTGTACTTCGTTCAGTTACGTATTGCTAAGCAGTGGTAT*********GAGTACC------TTA--CAGTTCGATCTA");
            let test_read = String::from("                               CT-AGCAG----ATCACCGTAAGGACTACCAGACGTTTAGCC           ");

            let keyvalues = extract_tagged_sequences(&test_read, &reference);
            for (key, value) in keyvalues {
                assert_eq!(key, "*");
                assert_eq!(value, "CACCGTAAG");
            }
        }
    }

}