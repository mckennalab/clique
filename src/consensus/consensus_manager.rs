use crate::umis::sequenceclustering::BestHits;
use std::collections::HashMap;

pub struct ConsensusManager {
    observed_identifiers: HashMap<Vec<u8>,BestHits>,
    observed_identifiers_counts: HashMap<Vec<u8>,u64>,
    total_reads: u64,
}

impl ConsensusManager {
    pub fn new() -> ConsensusManager {
        let ids : HashMap<Vec<u8>,BestHits> = HashMap::new();
        let counts : HashMap<Vec<u8>,u64> = HashMap::new();
        ConsensusManager{ observed_identifiers: ids, observed_identifiers_counts: counts, total_reads: 0 }
    }

    pub fn add_hit(&mut self, barcode: &Vec<u8>, best_hit: BestHits) {
        self.total_reads = self.total_reads + 1;
        let new_total = *self.observed_identifiers_counts.entry(barcode.clone()).or_insert(0) + 1 as u64;
        self.observed_identifiers_counts.insert(barcode.clone(),new_total);
        self.observed_identifiers.insert(barcode.clone(),best_hit);
    }

    pub fn unified_consensus_list(&self) -> ConsensusManager {
        let mut total = 0;
        let mut unresolved = 0;
        let mut resolved = 0;
        let mut clean_mapping : HashMap<Vec<u8>,BestHits> = HashMap::new();
        let mut clean_ids : HashMap<Vec<u8>,u64> = HashMap::new();

        self.observed_identifiers.iter().for_each(|(k,v)| {
            let count = *self.observed_identifiers_counts.get(&k.clone()).unwrap();
            total += count;
            if v.hits.len() > 1 {
                let matching: Vec<&Vec<u8>> = v.hits.iter().filter(|hit| self.observed_identifiers_counts.contains_key(*hit)).collect::<Vec<&Vec<u8>>>();
                if matching.len() != 1 {
                    unresolved += count;
                } else {
                    resolved += count;
                    let best_hits = BestHits{ hits: vec![matching[0].clone()], distance: v.distance };
                    clean_mapping.insert(k.clone(),best_hits);
                    clean_ids.insert(k.clone(),count);
                }
            } else {
                clean_mapping.insert(k.clone(),v.clone());
                clean_ids.insert(k.clone(),count);
            }

        });
        println!("total = {}, unresolved = {}, resolved = {}",total, unresolved, resolved);
        ConsensusManager{
            observed_identifiers: clean_mapping,
            observed_identifiers_counts: clean_ids,
            total_reads: total,
        }
    }
}