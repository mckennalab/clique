use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use bio::io::fasta::*;
use itertools::Itertools;
use needletail::Sequence;

use suffix::SuffixTable;
use crate::alignment::fasta_bit_encoding::{FastaBase};
use crate::read_strategies::read_set::ReadSetContainer;


pub fn reference_file_to_structs(reference_file: &String, kmer_size: usize) -> Vec<Reference> {

    let reader = Reader::from_file(reference_file).unwrap();
    let fasta_entries: Vec<Record> = reader.records().map(|f| f.unwrap()).collect();

    let mut references = Vec::new();

    for ref_entry in fasta_entries {
        let seeds = ReferenceManager::find_seeds(&ref_entry.seq().to_vec(), kmer_size);

        references.push(Reference {
            sequence: FastaBase::from_vec_u8_default_ns(&ref_entry.seq().to_vec()),
            sequence_u8: ref_entry.seq().to_vec(),
            name: str::as_bytes(ref_entry.id()).to_vec(),
            suffix_table: seeds,
        });
    }
    references
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SuffixTableLookup<'s, 't> {
    pub suffix_table: SuffixTable<'s, 't>,
    pub seed_size: usize,
}


#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Reference<'s, 't> {
    pub sequence: Vec<FastaBase>,
    pub sequence_u8: Vec<u8>,
    pub name: Vec<u8>,
    pub suffix_table: SuffixTableLookup<'s, 't>,
}

impl Hash for Reference<'_, '_> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.sequence_u8.hash(state);
        self.name.hash(state);
    }
}

#[allow(dead_code)]
pub struct ReferenceManager<'a, 's, 't> {
    pub references: HashMap<usize, Reference<'a, 'a>>,
    pub reference_name_to_ref: HashMap<Vec<u8>, usize>,
    pub unique_kmers: UniqueKmerLookup<'s, 't>,
    pub kmer_size: usize,
    pub kmer_skip: usize,
    pub longest_ref: usize,
}

#[allow(dead_code)]
pub struct UniqueKmerLookup<'s, 't> {
    kmer_length: usize,
    kmer_to_reference: HashMap<Vec<u8>,Reference<'s, 't>>,
    reference_to_kmer: HashMap<Reference<'s, 't>,Vec<Vec<u8>>>,
    all_have_unique_mappings: bool,

}

impl <'a, 's, 't>ReferenceManager<'a, 's, 't> {

    pub fn from(fasta: &String, kmer_size: usize, kmer_spacing: usize) -> ReferenceManager {

        let references = reference_file_to_structs(&fasta, kmer_size);
        let longest_ref = references.iter().map(|r| r.sequence.len()).max().unwrap_or(0);
        let unique_kmers = ReferenceManager::unique_kmers(&references, &kmer_size, &kmer_spacing);
        let references = references.into_iter().enumerate().collect::<HashMap<usize,Reference>>();
        let reference_name_to_ref = references.iter().map(|(i,r)| (r.name.clone(),*i)).collect();
        ReferenceManager{ references, reference_name_to_ref, unique_kmers, kmer_size: kmer_size.clone(), kmer_skip: kmer_spacing.clone(), longest_ref }
    }

    /// Find the suffix array 'seeds' given a reference sequence
    ///
    /// # Arguments
    ///
    /// * `name` - a u8 Vec representing the reference sequence
    /// * `seed_size` - not used in the suffix array creation, but tracked for each analysis
    pub fn find_seeds(reference: &Vec<u8>, seed_size: usize) -> SuffixTableLookup<'s, 't> {
        return SuffixTableLookup { suffix_table: SuffixTable::new(String::from_utf8(reference.clone()).unwrap()), seed_size };
    }

    pub fn sequence_to_kmers(reference: &Vec<u8>, kmer_size: &usize, kmer_spacing: &usize) -> Vec<(Vec<u8>, usize)> {
        let mut kmers: Vec<(Vec<u8>,usize)> = reference.to_ascii_uppercase().windows(*kmer_size).step_by(*kmer_spacing).dedup_with_count().map(|(c,w)| (w.to_vec(),c)).collect::<Vec<(Vec<u8>,usize)>>();
        kmers.extend(reference.reverse_complement().to_ascii_uppercase().windows(*kmer_size).step_by(*kmer_spacing).dedup_with_count().map(|(c,w)| (w.to_vec(),c)).collect::<Vec<(Vec<u8>,usize)>>());
        kmers
    }

    /// find a list of unique kmers per reference
    pub fn unique_kmers(references: &Vec<Reference<'s, 't>>, kmer_size: &usize, kmer_spacing: &usize) -> UniqueKmerLookup<'s, 't> {
        let mut kmer_counts = HashMap::new();

        for reference in references {
            let kmers = ReferenceManager::sequence_to_kmers(&reference.sequence_u8, kmer_size, kmer_spacing);
            for kmer in kmers {
                kmer_counts.insert(kmer.0.clone(),if kmer_counts.contains_key(&kmer.0) {kmer_counts.get(&kmer.0).unwrap()} else {&0} + kmer.1);
            }
        }

        let mut reference_to_unique: HashMap<Reference,Vec<Vec<u8>>> = HashMap::new();
        let mut unique_kmer_to_reference = HashMap::new();

        let mut all_unique = true;
        for reference in references {
            let kmers = ReferenceManager::sequence_to_kmers(&reference.sequence_u8, kmer_size, kmer_spacing);
            let unique_kmers = kmers.iter().filter(|(k,_c)| kmer_counts.contains_key(k) && *kmer_counts.get(k).unwrap() == 1).collect_vec();

            if unique_kmers.is_empty() {
                warn!("Unique kmer count for reference {} is empty!", String::from_utf8(reference.name.clone()).unwrap());
                all_unique = false;
            }
            for (kmer,_count) in &unique_kmers {
                unique_kmer_to_reference.insert(kmer.clone(),reference.clone());
            }
            reference_to_unique.insert(reference.clone(), unique_kmers.into_iter().map(|(v,_c)|v.clone()).collect());

        }

        UniqueKmerLookup{ kmer_length: 0, kmer_to_reference: unique_kmer_to_reference, reference_to_kmer: reference_to_unique, all_have_unique_mappings: all_unique }
    }

    #[allow(dead_code)]
    pub fn match_references(&self, read: ReadSetContainer) -> Vec<&Reference> {
        let read_kmers = ReferenceManager::sequence_to_kmers(&read.read_one.seq().to_vec(), &self.kmer_size, &self.kmer_skip );

        // now collect reference that have unique kmers matching this sequence
        let mut votes = HashMap::new();
        for (kmer,_count) in read_kmers {
            if self.unique_kmers.kmer_to_reference.contains_key(&kmer) {
                let reference= self.unique_kmers.kmer_to_reference.get(&*kmer).unwrap();
                *votes.entry(reference).or_insert(0) += 1;

            }
        }
        votes.keys().into_iter().map(|k|*k).collect()
    }

}


#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use super::*;
    use crate::read_strategies::read_set::ReadIterator;

    #[test]
    fn test_kmer_creation_from_large_library() {
        // basic does it work test
        let order_fastas = String::from("test_data/18guide1_pcr_sequence.fasta");
        let rm = ReferenceManager::from(&order_fastas,15, 5 );

        for (reference,kmers) in rm.unique_kmers.reference_to_kmer {
            println!("Reference name {} and count {}",String::from_utf8(reference.name).unwrap(),kmers.len());
        }
        assert_eq!(rm.references.len(),180);
    }

    #[test]
    fn test_kmer_creation_from_two_libs() {
        // basic does it work test
        let order_fastas = String::from("test_data/two_references_just_one.fa");
        let rm = ReferenceManager::from(&order_fastas,15, 5 );
        assert_eq!(rm.references.len(),1);

        for (reference,kmers) in rm.unique_kmers.reference_to_kmer {
            println!("kmers {} {}",String::from_utf8(reference.name.clone()).unwrap(),kmers.iter().map(|c| String::from_utf8(c.clone()).unwrap()).join(","));
            assert!(kmers.contains(&"GGGCGAGATCAAGCA".as_bytes().to_vec()));
        }
        let order_fastas = String::from("test_data/two_references.fa");
        let rm = ReferenceManager::from(&order_fastas,15, 5 );
        assert_eq!(rm.references.len(),2);

        for (_reference_index,kmers) in rm.unique_kmers.reference_to_kmer.iter().enumerate() {
            println!("kmers {} {}",String::from_utf8(kmers.0.name.clone()).unwrap(),kmers.1.iter().map(|c| String::from_utf8(c.clone()).unwrap()).join(","));
            assert!(!kmers.1.contains(&"TCACCTATTAGCGGCTAA".as_bytes().to_vec()));
            if kmers.0.name == "cas_tag".as_bytes().to_vec() {
                assert!(kmers.1.contains(&"TTTTTTTTTTTTTTC".as_bytes().to_vec()));
            } else {
                assert!(kmers.1.contains(&"AAAAAAAAAAAATTC".as_bytes().to_vec()));
            }
        }
    }

    #[test]
    fn test_alignment_to_large_library() {
        let order_fastas = String::from("test_data/18guide1_pcr_sequence.fasta");
        let rm = ReferenceManager::from(&order_fastas,8, 5 );

        let read_iterator = ReadIterator::new(PathBuf::from("test_data/PAM_TWIST_1_018_S20_merged_001.fastq.gz"),
                                              None,
                                              None,
                                              None);

        println!("Running reads");
        let mut total = 0;
        let mut found = 0;
        let mut conflicted = 0;
        let mut missing = 0;
        for read in read_iterator {
            //println!("Length of hits {} ", rm.best_reference(read).len());
            let cnt = rm.match_references(read.clone()).len();
            total += 1;
            if cnt > 0 {
                found += 1;
                if cnt > 1 {
                    conflicted += 1;
                }
            } else {
                missing += 1;
                if missing < 200 {
                    println!("Not found {}",String::from_utf8(read.read_one.seq().to_vec()).unwrap());
                }
            }
        }
        println!("Found {} from total {}, conflicted {} missing {}",found,total,conflicted,missing);


    }
}