use num_traits::Pow;
use rand::{seq::IteratorRandom, rng};
use crate::read_strategies::read_set::ReadSetContainer;

#[allow(dead_code)]
pub fn phred_to_prob(phred: &u8) -> f64 {
    let phred_f64 = ((*phred as usize) - 33) as f64;
    (10.0).pow((-1.0 * phred_f64)/10.0)
}

#[allow(dead_code)]
pub fn prob_to_phred(qual: f64) -> u8 {
    (((-10.0) * qual.log10()) + 33.0) as u8
}

pub fn u8s(u8s: &Vec<u8>) -> String {
    String::from_utf8(u8s.clone()).unwrap()
}

#[allow(dead_code)]
// TODO: BUG - In the `false` (disagree) branch, the formula `1.0 - ((1.0 - prob2) * (1.0 * prob1))`
// simplifies to `1.0 - prob1 + prob1*prob2`. The `1.0 * prob1` should likely be `(1.0 - prob1)`
// to correctly compute `1.0 - (1.0 - prob1) * (1.0 - prob2)`, which is the probability that at
// least one of the two independent error events occurred. The current formula uses `1.0 * prob1`
// which is just `prob1` and produces incorrect disagreement quality scores.
pub fn combine_phred_scores(phred_one: &u8, phred_two: &u8, agree: bool) -> u8 {
    let prob1 = phred_to_prob(phred_one);
    let prob2 = phred_to_prob(phred_two);

    match agree {
        true => {
            prob_to_phred(prob1 * prob2)
        }
        false => {
            prob_to_phred(1.0 - ((1.0 - prob2) * ( 1.0 * prob1)))
        }
    }
}

pub fn strip_gaps(bases: &Vec<u8>) -> Vec<u8> {
    bases.iter().filter(|x| **x != b'-').map(|x|*x).collect()
}

pub fn pad_right(v: &Vec<u8>, target_len: usize, pad_byte: u8) -> Vec<u8> {
    let mut vv = v.clone();
    vv.resize(target_len, pad_byte);
    vv
}

pub fn reverse_complement(dna: &[u8]) -> Vec<u8> {
    dna.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            b'R' => b'Y', // purine <-> pyrimidine
            b'Y' => b'R',
            b'S' => b'S',
            b'W' => b'W',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'D' => b'H',
            b'H' => b'D',
            b'V' => b'B',
            b'N' => b'N',
            other => other, // unrecognized base: leave unchanged
        })
        .collect()
}
// TODO: BUG - `choose_multiple` samples WITHOUT replacement from the `bases` vector of length 4.
// This means if `length` > 4, it will only return 4 elements (the entire vector), silently
// truncating the output. For `length` <= 4, it returns a subset without repeats, which means
// it cannot generate sequences like "AAAA". This should use `choose` with replacement in a loop,
// or use `(0..length).map(|_| *bases.choose(&mut rng).unwrap()).collect()`.
pub fn random_sequence(length: usize) -> String {
    let bases = vec![b'A', b'C', b'G', b'T'];
    let mut rng = rng();

    String::from_utf8(bases.iter().choose_multiple(&mut rng, length).iter().map(|c| **c).collect::<Vec<u8>>()).unwrap()
}

pub fn all_combinations(n: usize) -> Vec<String> {
    let characters = vec!["A", "C", "G", "T"];

    (2..n).fold(
        characters.iter().map(|c| characters.iter().map(move |&d| d.to_owned() + *c)).flatten().collect(),
        |acc, _| acc.into_iter().map(|c| characters.iter().map(move |&d| d.to_owned() + &*c)).flatten().collect(),
    )
}

pub fn create_fake_quality_scores(length: usize) -> Vec<u8> {
    vec![b'H'; length]
}

#[allow(dead_code)]
pub fn fake_reads(full_length: usize, permutation_leader_size: usize) -> Vec<ReadSetContainer> {
    //sort_string: &Vec<u8>, read: &ReadSetContainer
    let mut fake_reads = Vec::new();
    let all_perm = all_combinations(permutation_leader_size);
    println!("All perm size {} for size {}", all_perm.len(), permutation_leader_size);
    for sequence_leader in all_perm {
        let mut read_seq = sequence_leader.clone();
        read_seq.push_str(random_sequence(full_length - permutation_leader_size).as_str());
        let qual = create_fake_quality_scores(full_length);
        let record = bio::io::fastq::Record::with_attrs("fakeRead", None, read_seq.as_bytes(), qual.as_slice());
        let fake_rsc = ReadSetContainer::new_from_read1(record);
        fake_reads.push(fake_rsc);
    }
    fake_reads
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phred_to_qual_test() {
        assert_eq!(phred_to_prob(&b'I'), 0.0001);
        assert_eq!(phred_to_prob(&b'H'), 0.00012589254117941674);
        assert_eq!(phred_to_prob(&b'+'), 0.1);
        assert_eq!(phred_to_prob(&b'5'), 0.01);
    }

    #[test]
    fn qual_to_phred_test() {
        assert_eq!(prob_to_phred(0.0001), b'I');
        assert_eq!(prob_to_phred(0.00012589254117941674), b'H');
        assert_eq!(prob_to_phred(0.1), b'+');
        assert_eq!(prob_to_phred(0.01), b'5');
    }

    #[test]
    fn combine_qual_test() {
        assert_eq!(combine_phred_scores(&b'H',&b'+', false), b'!');
        assert_eq!(combine_phred_scores(&b'H',&b'+', true), b'R');
    }

    #[test]
    fn test_reverse_complement_standard() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"TTTT"), b"AAAA");
        assert_eq!(reverse_complement(b"CCCC"), b"GGGG");
        assert_eq!(reverse_complement(b"GGGG"), b"CCCC");
    }

    #[test]
    fn test_reverse_complement_palindrome() {
        assert_eq!(reverse_complement(b"AATT"), b"AATT");
        assert_eq!(reverse_complement(b"GCGC"), b"GCGC");
    }

    #[test]
    fn test_reverse_complement_single_base() {
        assert_eq!(reverse_complement(b"A"), b"T");
        assert_eq!(reverse_complement(b"T"), b"A");
        assert_eq!(reverse_complement(b"G"), b"C");
        assert_eq!(reverse_complement(b"C"), b"G");
        assert_eq!(reverse_complement(b"N"), b"N");
    }

    #[test]
    fn test_reverse_complement_empty() {
        assert_eq!(reverse_complement(b""), b"");
    }

    #[test]
    fn test_reverse_complement_degenerate_bases() {
        // R (A/G) -> complement Y (C/T), reversed
        assert_eq!(reverse_complement(b"R"), b"Y");
        assert_eq!(reverse_complement(b"Y"), b"R");
        assert_eq!(reverse_complement(b"S"), b"S"); // S = C/G, complement is also S
        assert_eq!(reverse_complement(b"W"), b"W"); // W = A/T, complement is also W
        assert_eq!(reverse_complement(b"K"), b"M");
        assert_eq!(reverse_complement(b"M"), b"K");
        assert_eq!(reverse_complement(b"B"), b"V");
        assert_eq!(reverse_complement(b"V"), b"B");
        assert_eq!(reverse_complement(b"D"), b"H");
        assert_eq!(reverse_complement(b"H"), b"D");
    }

    #[test]
    fn test_reverse_complement_lowercase() {
        assert_eq!(reverse_complement(b"acgt"), b"ACGT");
    }

    #[test]
    fn test_reverse_complement_double_application_is_identity() {
        let seq = b"ACGTRYSWKMBDHVN";
        let rc = reverse_complement(seq);
        let rc_rc = reverse_complement(&rc);
        // Double RC should give back the original (uppercased)
        let uppercased: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
        assert_eq!(rc_rc, uppercased);
    }

    #[test]
    fn test_strip_gaps() {
        assert_eq!(strip_gaps(&vec![b'A', b'-', b'C', b'-', b'G']), vec![b'A', b'C', b'G']);
        assert_eq!(strip_gaps(&vec![b'A', b'C', b'G']), vec![b'A', b'C', b'G']);
        assert_eq!(strip_gaps(&vec![b'-', b'-', b'-']), Vec::<u8>::new());
        assert_eq!(strip_gaps(&vec![]), Vec::<u8>::new());
    }

    #[test]
    fn test_pad_right() {
        assert_eq!(pad_right(&vec![b'A', b'C'], 5, b'-'), vec![b'A', b'C', b'-', b'-', b'-']);
        assert_eq!(pad_right(&vec![b'A', b'C'], 2, b'-'), vec![b'A', b'C']);
        assert_eq!(pad_right(&vec![], 3, b'N'), vec![b'N', b'N', b'N']);
    }

    #[test]
    fn test_pad_right_shorter_target() {
        // When target_len < current length, Vec::resize truncates
        assert_eq!(pad_right(&vec![b'A', b'C', b'G'], 1, b'-'), vec![b'A']);
    }

    #[test]
    fn test_u8s() {
        assert_eq!(u8s(&vec![b'A', b'C', b'G', b'T']), "ACGT");
        assert_eq!(u8s(&vec![]), "");
    }

    #[test]
    fn test_create_fake_quality_scores() {
        let quals = create_fake_quality_scores(5);
        assert_eq!(quals.len(), 5);
        assert!(quals.iter().all(|q| *q == b'H'));
    }

    #[test]
    fn test_create_fake_quality_scores_zero() {
        let quals = create_fake_quality_scores(0);
        assert_eq!(quals.len(), 0);
    }

    #[test]
    fn test_all_combinations_length_2() {
        let combos = all_combinations(2);
        // For n=2, fold range is 2..2 which is empty, so result is the initial
        // characters.iter().map(...).flatten() which gives all 2-mers = 16
        assert_eq!(combos.len(), 16);
        assert!(combos.contains(&"AA".to_string()));
        assert!(combos.contains(&"TT".to_string()));
        assert!(combos.contains(&"AC".to_string()));
    }

    #[test]
    fn test_all_combinations_length_3() {
        let combos = all_combinations(3);
        // For n=3, should give 4^3 = 64 combinations
        assert_eq!(combos.len(), 64);
    }

    #[test]
    fn test_phred_roundtrip() {
        // Converting to prob and back should yield the same phred
        for phred in [b'!', b'+', b'5', b'I'] {
            let prob = phred_to_prob(&phred);
            let back = prob_to_phred(prob);
            assert_eq!(back, phred, "Roundtrip failed for phred {}", phred);
        }
    }

    #[test]
    fn test_phred_to_prob_boundaries() {
        // Phred 33 (ASCII '!') = quality 0 => error prob = 1.0
        assert_eq!(phred_to_prob(&b'!'), 1.0);
    }
}