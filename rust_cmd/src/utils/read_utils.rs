use num_traits::Pow;
use rand::{seq::IteratorRandom, thread_rng};
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
pub fn random_sequence(length: usize) -> String {
    let bases = vec![b'A', b'C', b'G', b'T'];
    let mut rng = thread_rng();

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

}