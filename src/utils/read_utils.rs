use rand::{seq::IteratorRandom, thread_rng};
use crate::read_strategies::read_set::ReadSetContainer;

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
