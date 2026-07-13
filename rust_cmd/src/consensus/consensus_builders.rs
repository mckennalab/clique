//! Turning a group of reads that share all their UMI tags into a single
//! per-molecule consensus read, and writing it out with its annotation tags
//! (read count, alignment rate, score, and the `ce` called-edit string).

use crate::alignment::alignment_matrix::AlignmentTag;
use crate::alignment::scoring_functions::AffineScoring;
use crate::alignment_manager::{simplify_cigar_string, OutputAlignmentWriter};
use crate::read_strategies::read_disk_sorter::SortingReadSetContainer;
use crate::read_strategies::sequence_layout::SequenceLayout;
use crate::events::call_read_events;
use crate::reference::fasta_reference::ReferenceManager;
use counter::Counter;
use rust_htslib::bam::record::CigarString;
use shardio::{Range, ShardReader};
use std::cmp::Ordering;
use std::collections::{HashMap, VecDeque};
use std::convert::TryFrom;
use std::sync::{Arc, Mutex};
use num_traits::{Pow, ToPrimitive};
use ::{FASTA_N, FASTA_UNSET};

#[allow(dead_code)]
const PHRED_OFFSET: u8 = 32;

pub enum ReadOutputApproach {
    Collapse,
    Correct,
}

#[allow(dead_code)]
pub enum MergeStrategy {
    StrictConsensus,
    Hybrid,
    Stretcher,
}

// this function isn't really in the right place, but it matches the consensus output function
pub fn write_corrected_reads(
    reader: &ShardReader<SortingReadSetContainer>,
    writer: &mut dyn OutputAlignmentWriter,
    levels: usize,
    reference_manager: &ReferenceManager,
    read_structure: &SequenceLayout,
) {
    let mut processed_reads = 0;

    reader.iter_range(&Range::all()).unwrap().for_each(|x| {
        let x = x.unwrap();
        assert_eq!(x.ordered_sorting_keys.len(), levels);

        let mut read_pile = VecDeque::new();
        read_pile.push_front(x);

        let new_read = create_consensus_sam_read(
            reference_manager,
            read_structure,
            &0,
            &read_pile,
            &AffineScoring::default_dna(),
            &MergeStrategy::Stretcher,
        );

        match new_read {
            None => {}
            Some(x) => {
                writer.write_read(&x.read, &x.added_tags).expect("Unable to write a read to the arc writer (LOC1)");
                processed_reads += 1;
                if processed_reads % 1000000 == 0 {
                    info!("Processed {} reads into output file", processed_reads);
                }
            }
        }
    });

    writer.close().unwrap();

}


pub fn write_consensus_reads(
    reader: &ShardReader<SortingReadSetContainer>,
    writer: &mut dyn OutputAlignmentWriter,
    levels: usize,
    reference_manager: &ReferenceManager,
    read_structure: &SequenceLayout,
    maximum_reads_before_downsampling: &usize,
    merge_strategy: &MergeStrategy,
    processing_threads: &usize,
) {
    let mut last_read: Option<SortingReadSetContainer> = None;
    let mut buffered_reads = VecDeque::new();

    let processed_reads = Arc::new(Mutex::new(0));

    let score = AffineScoring::default_dna();

    let arc_output = Arc::new(Mutex::new(writer));

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(*processing_threads)
        .build()
        .unwrap();

    pool.scope(|s| {
        let _handled_reads = 0;
        reader.iter_range(&Range::all()).unwrap().for_each(|x| {
            let x = x.unwrap();
            assert_eq!(x.ordered_sorting_keys.len(), levels);
            if !(last_read.is_some() && x.cmp(last_read.as_ref().unwrap()) == Ordering::Equal)
                && !buffered_reads.is_empty()
            {
                let my_buffered_reads = buffered_reads.clone();
                buffered_reads = VecDeque::new();
                s.spawn(|_y| {
                    let my_buffered_reads = my_buffered_reads;

                    let new_read = create_consensus_sam_read(
                        reference_manager,
                        read_structure,
                        maximum_reads_before_downsampling,
                        &my_buffered_reads,
                        &AffineScoring::default_dna(),
                        merge_strategy,
                    );

                    match new_read {
                        None => {}
                        Some(x) => {
                            let arc_writer = arc_output.clone();
                            let mut arc_writer = arc_writer.lock().expect("Unable to access multi-threaded writer");
                            arc_writer.write_read(&x.read, &x.added_tags).expect("Unable to write a read to the arc writer (LOC1)");

                            let processed_reads = processed_reads.clone();
                            let mut processed_reads = processed_reads.lock().expect("Unable to lock processed read count");
                            let current_proc_read = *processed_reads;
                            *processed_reads += my_buffered_reads.len();
                            if (*processed_reads as f64 / 100000.0).floor() - (current_proc_read as f64 / 100000.0).floor() >= 1.0 {
                                info!("Processed {} reads into their consensus", processed_reads);
                            }
                        }
                    }
                });
            }

            buffered_reads.push_back(x.clone());
            last_read = Some(x);
        });
    });
    if !buffered_reads.is_empty() {
        let new_read = create_consensus_sam_read(
            reference_manager,
            read_structure,
            maximum_reads_before_downsampling,
            &buffered_reads,
            &score,
            merge_strategy,
        );

        match new_read {
            None => {}
            Some(x) => {
                let arc_writer = arc_output.clone();
                let mut arc_writer = arc_writer.lock().expect("Unable to access multi-threaded writer");
                arc_writer.write_read(&x.read, &x.added_tags).expect("Unable to write a read to the arc writer (LOC2)");
            }
        }
    }

    let arc_writer = arc_output.clone();
    let mut arc_writer = arc_writer.lock().expect("Unable to access multi-threaded writer");
    arc_writer.close().unwrap();
}

pub struct SamReadyOutput {
    pub read: SortingReadSetContainer,
    pub added_tags: HashMap<[u8; 2], String>,
}

/// BAM tag under which called CRISPR edit events are recorded (McKenna
/// per-target event string, targets joined by `_`). Uses a `c`-prefix to stay
/// clear of the `e<symbol>` / `o<symbol>` extracted-UMI tag namespace.
pub const EVENT_TAG: [u8; 2] = [b'c', b'e'];

/// Select the reads actually consumed by consensus generation. A zero limit
/// disables downsampling, as used by the corrected-read output path.
fn select_reads_for_consensus<'a>(
    buffered_reads: &'a VecDeque<SortingReadSetContainer>,
    maximum_reads_before_downsampling: usize,
) -> Vec<&'a SortingReadSetContainer> {
    let limit = if maximum_reads_before_downsampling == 0 {
        buffered_reads.len()
    } else {
        maximum_reads_before_downsampling.min(buffered_reads.len())
    };

    buffered_reads.iter().take(limit).collect()
}

/// Call CRISPR edit events for the given aligned pair against the reference's
/// declared targets and, if any targets exist, record them under [`EVENT_TAG`].
fn insert_event_tag(
    added_tags: &mut HashMap<[u8; 2], String>,
    read_structure: &SequenceLayout,
    reference_name: &str,
    reference_aligned: &[u8],
    read_aligned: &[u8],
) {
    if let Some(reference_record) = read_structure.references.get(reference_name) {
        let events = call_read_events(reference_aligned, read_aligned, reference_record);
        if !events.is_empty() {
            added_tags.insert(EVENT_TAG, events);
        }
    }
}


#[allow(deprecated)]
pub fn create_consensus_sam_read(
    reference_manager: &ReferenceManager,
    read_structure: &SequenceLayout,
    maximum_reads_before_downsampling: &usize,
    buffered_reads: &VecDeque<SortingReadSetContainer>,
    _my_aff_score: &AffineScoring,
    merge_strategy: &MergeStrategy,
) -> Option<SamReadyOutput> {
    assert!(
        !buffered_reads.is_empty(),
        "Cannot create a consensus from an empty read collection"
    );
    let consensus_reads = select_reads_for_consensus(
        buffered_reads,
        *maximum_reads_before_downsampling,
    );

    let mut added_tags = HashMap::new();
    added_tags.insert([b'r', b'c'], buffered_reads.len().to_string());
    added_tags.insert([b'd', b'c'], consensus_reads.len().to_string());

    if consensus_reads.len() > 1 {
        let consensus_reference = Counter::<Vec<u8>, usize>::init(
            consensus_reads
                .iter()
                .map(|x| x.aligned_read.reference_name.clone().as_bytes().to_vec())
                .collect::<Vec<Vec<u8>>>(),
        ).most_common_ordered();

        let top_ref = &consensus_reference.get(0).unwrap().clone().0;
        let reference_pointer = reference_manager
            .references
            .get(
                reference_manager
                    .reference_name_to_ref
                    .get(top_ref)
                    .unwrap(),
            )
            .unwrap();

        let consensus_alignment = match merge_strategy {
            MergeStrategy::StrictConsensus => {
                unimplemented!("Removing SPOA");
            }
            MergeStrategy::Hybrid => {
                unimplemented!("Removing SPOA");
            }
            MergeStrategy::Stretcher => {
                let mut candidate = crate::consensus::stretcher::AlignmentCandidate::new(reference_pointer.sequence.as_slice(), reference_pointer.name.as_slice());

                let valid: usize = consensus_reads.iter().map(|x| {
                    match candidate.add_alignment(&x.aligned_read) {
                        Ok(_) => { 0 }
                        Err(_) => { 1 }
                    }
                }).sum();

                // Skip the whole molecule if ANY read failed to merge: `add_alignment` mutates the
                // candidate on its error path (pushing the read name and partial counts), so a
                // tolerated failed read would pollute the consensus. Discarding is safer.
                if valid >= 1 {
                    None
                } else {
                    Some(candidate.to_consensus(&0.75))
                }
            }
        };

        match consensus_alignment {
            None => {
                // No usable consensus (a read failed to merge): skip this molecule. Both callers
                // treat a `None` return as "skip", so this must not panic and abort the run.
                return None;
            }
            Some(con) => {
                let read_names = consensus_reads
                    .iter()
                    .map(|x| x.aligned_read.read_name.clone())
                    .collect::<Vec<String>>();

                added_tags.insert([b'a', b'r'], read_names.join(","));
                added_tags.insert(
                    [b'r', b'm'],
                    get_reference_alignment_rate(
                        &con.reference_aligned,
                        &con.read_aligned,
                    )
                        .to_string(),
                );


                added_tags.insert([b'a', b's'], con.score.to_string());

                let consensus_ref_name = String::from_utf8_lossy(top_ref).into_owned();
                insert_event_tag(
                    &mut added_tags,
                    read_structure,
                    &consensus_ref_name,
                    &con.reference_aligned,
                    &con.read_aligned,
                );

                let new_sorting_read = consensus_reads[0].with_new_alignment(con);

                Some(SamReadyOutput { read: new_sorting_read, added_tags })
            }
        }
    } else {
        let single_read = (*consensus_reads[0]).clone();
        added_tags.insert([b'a', b'r'], single_read.aligned_read.read_name.clone());
        added_tags.insert(
            [b'r', b'm'],
            get_reference_alignment_rate(
                &single_read.aligned_read.reference_aligned,
                &single_read.aligned_read.read_aligned,
            )
                .to_string(),
        );


        added_tags.insert([b'a', b's'], single_read.aligned_read.score.to_string());

        insert_event_tag(
            &mut added_tags,
            read_structure,
            &single_read.aligned_read.reference_name,
            &single_read.aligned_read.reference_aligned,
            &single_read.aligned_read.read_aligned,
        );

        Some(SamReadyOutput { read: single_read, added_tags })
    }
}

pub fn get_reference_alignment_rate(reference: &[u8], read: &[u8]) -> f64 {
    let mut matches = 0;
    let mut mismatches = 0;

    //println!("reference: {} read: {}", FastaBase::string(&reference),FastaBase::string(&read));
    for index in 0..reference.len() {
        if *reference.get(index).unwrap() > 64 as u8 &&
            reference.get(index).unwrap() != &FASTA_N &&
            *read.get(index).unwrap() > 64 as u8
        {
            if reference.get(index).unwrap() == read.get(index).unwrap() {
                matches += 1;
            } else {
                mismatches += 1;
            }
        }
    }
    //println!("matches: {} mismatches: {}", matches, mismatches);
    (matches as f64) / ((matches + mismatches) as f64)
}

#[allow(dead_code)]
pub fn reference_read_to_cigar_string(
    reference_seq: &Vec<u8>,
    read_seq: &[u8],
) -> CigarString {
    let mut result: Vec<AlignmentTag> = Vec::new();

    for index in 0..reference_seq.len() {
        if reference_seq.get(index).unwrap().eq(&FASTA_UNSET) {
            result.push(AlignmentTag::Ins(1))
        } else if read_seq.get(index).unwrap().eq(&FASTA_UNSET) {
            result.push(AlignmentTag::Del(1))
        } else {
            result.push(AlignmentTag::MatchMismatch(1))
        }
    }

    let simplied_cigar = simplify_cigar_string(&result);
    CigarString::try_from(
        simplied_cigar
            .iter()
            .map(|m| format!("{}", m))
            .collect::<Vec<String>>()
            .join("")
            .as_bytes(),
    )
        .expect("Unable to parse cigar string.")
}


#[allow(dead_code)]
pub fn calculate_conc_qual_score(alignments: &Vec<Vec<u8>>, quality_scores: &Vec<Vec<u8>>) -> (Vec<u8>, Vec<u8>) {
    assert_eq!(alignments.len() - 1, quality_scores.len());

    let mut conc = Vec::new();
    let mut final_quals = Vec::new();

    let mut sequence_indexes = vec![0_usize; alignments.len()];

    let ln = alignments.get(0).unwrap().len();

    let reference = alignments.get(0).unwrap();

    (0..ln).for_each(|index| {
        let mut bases = Vec::new();
        let mut quals = Vec::new();
        let _ = &alignments[1..alignments.len()].iter().enumerate().for_each(|(sequence_index, x)| {
            assert_eq!(ln, x.len());

            let base = x[index];
            let qual = match base {
                b'-' => 20,
                _ => quality_scores[sequence_index][sequence_indexes[sequence_index]]
            };

            sequence_indexes[sequence_index] += match base {
                b'-' => { 0 }
                _ => { 1 }
            };
            bases.push(base);
            quals.push(qual);
        });
        let qual_scores = combine_qual_scores(vec![bases.as_slice()].as_slice(), vec![quals.as_slice()].as_slice(), &reference[index], &0.99);
        let index_of_max: usize = qual_scores
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.total_cmp(b))
            .map(|(index, _)| index).unwrap();

        if index_of_max < 5 {
            let prob = prob_to_phred(&qual_scores[index_of_max]);

            final_quals.push(prob);

            match index_of_max {
                0 => conc.push(b'A'),
                1 => conc.push(b'C'),
                2 => conc.push(b'G'),
                3 => conc.push(b'T'),
                4 => conc.push(b'-'),
                _ => panic!("Unknown index"),
            };
        }
    });
    (conc, final_quals)
}

pub fn phred_to_error_prob(phred: &u8) -> f64 {
    (10.0_f64).pow((phred.to_f64().unwrap()) / (-10.0)) // 32.9999 to avoid Inf powers
}

pub fn prob_to_phred(prob: &f64) -> u8 {
    // the upper bound is a bit arbitrary, but 33 + 93 = 126, the highest possible value reported by Illumina
    if f64::is_nan(*prob) {
        return 0;
    }

    assert!(prob >= &0.0_f64 && prob <= &1.0_f64, "{}", format!("Unable to format prob {}", prob));
    if prob < &0.00000001_f64 {
        return 0_u8;
    }

    // TODO: we dont deal with phred + 64 format data
    let ret = (-10.0) * (1.00000000001 - prob).log10(); // again to prevent zero getting in, we subtract from 1 + epsilon
    let ret = ret.round().to_u8().unwrap();
    let ret = if ret > 40 { // cap PHRED at 40; Noodles doesn't like higher TODO: fix this
        40_u8
    } else {
        ret as u8
    };
    //assert!(ret >= 0 && ret <= 40_u8);
    ret
}

pub(crate) fn combine_qual_scores(bases: &[&[u8]], scores: &[&[u8]], reference_base: &u8, reference_prob: &f64) -> [f64; 5] {
    let mut allele_props = [((1.0_f64 - reference_prob) / 4.0).log2(); 5];
    match *reference_base {
        b'A' | b'a' => { allele_props[0] = (*reference_prob).log2(); }
        b'C' | b'c' => { allele_props[1] = (*reference_prob).log2(); }
        b'G' | b'g' => { allele_props[2] = (*reference_prob).log2(); }
        b'T' | b't' => { allele_props[3] = (*reference_prob).log2(); }
        b'N' | b'n' => { allele_props[4] = (*reference_prob).log2(); }
        _ => {
            debug!("unaccounted for quality score");
        }
    };

    assert_eq!(bases.len(), scores.len());

    bases.iter().zip(scores.iter()).for_each(|(base_set, qual_set)| {
        assert_eq!(base_set.len(), qual_set.len());
        for i in 0..base_set.len() {
            let base = base_set[i];
            let qs = qual_set[i];
            //base_set.iter().zip(qual_set).for_each(|(base, qs)|{
            let base_id = match base {
                b'A' | b'a' => { 0 }
                b'C' | b'c' => { 1 }
                b'G' | b'g' => { 2 }
                b'T' | b't' => { 3 }
                b'N' | b'n' => { 4 }
                _ => {
                    debug!("unaccounted for quality score");
                    5
                }
            };
            if base_id < 5 {
                (0..5).for_each(|i| {
                    if i == base_id {
                        allele_props[i] = allele_props[i] + (1.0 - phred_to_error_prob(&qs)).log2();
                        //println!("MT id {} new prop {} {} {}", i, allele_props[i],&qs,phred_to_error_prob(&qs));
                    } else {
                        allele_props[i] = allele_props[i] + (phred_to_error_prob(&qs) / 3.0_f64).log2();
                        //println!("NM id {} new prop {}", i, allele_props[i]);
                    }
                });
            }
            //println!();
        };
    });
    //println!("combined : {:?} -- {:?}",allele_props,allele_props.iter().map(|x| 2.0_f64.pow(x)));
    //println!("combined 3: {:?}",calculate_qual_scores(&mut allele_props));
    calculate_qual_scores(&mut allele_props)
}

pub fn calculate_qual_scores(allele_props: &mut [f64; 5]) -> [f64; 5] {
    let total: f64 = allele_props.iter().map(|x| 2.0_f64.pow(x)).sum();
    [2.0_f64.pow(allele_props[0]) / total,
        2.0_f64.pow(allele_props[1]) / total,
        2.0_f64.pow(allele_props[2]) / total,
        2.0_f64.pow(allele_props[3]) / total,
        2.0_f64.pow(allele_props[4]) / total]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::alignment_matrix::AlignmentResult;
    use crate::read_strategies::read_disk_sorter::SortingReadSetContainer;
    use crate::read_strategies::sequence_layout::ReferenceRecord;
    use rust_htslib::bam::record::Cigar;
    use std::collections::{BTreeMap, VecDeque};

    #[cfg(feature = "spoa")]
    use FASTA_A;
    #[cfg(feature = "spoa")]
    use read_strategies::read_disk_sorter::CorrectedKey;
    #[cfg(feature = "spoa")]
    use utils::read_utils::u8s;

    fn downsampling_layout() -> SequenceLayout {
        let mut references = BTreeMap::new();
        references.insert(
            "reference".to_string(),
            ReferenceRecord {
                sequence: "AAAA".to_string(),
                umi_configurations: BTreeMap::new(),
                targets: vec![],
                target_types: vec![],
                target_locations: Some(vec![]),
            },
        );

        SequenceLayout {
            aligner: None,
            merge: None,
            reads: vec![],
            known_strand: true,
            references,
        }
    }

    fn downsampling_read(read_name: &str, bases: &[u8]) -> SortingReadSetContainer {
        SortingReadSetContainer::empty_tags(AlignmentResult {
            reference_name: "reference".to_string(),
            read_name: read_name.to_string(),
            reference_aligned: b"AAAA".to_vec(),
            read_aligned: bases.to_vec(),
            read_quals: Some(vec![40; bases.len()]),
            cigar_string: vec![AlignmentTag::MatchMismatch(bases.len())],
            path: vec![],
            score: 4.0,
            reference_start: 0,
            read_start: 0,
            bounding_box: None,
        })
    }

    #[test]
    fn test_consensus_downsampling_caps_consumed_reads() {
        let layout = downsampling_layout();
        let reference_manager = ReferenceManager::from_yaml_input(&layout, 2, 1);
        let mut reads = VecDeque::new();
        reads.push_back(downsampling_read("included_1", b"AAAA"));
        reads.push_back(downsampling_read("included_2", b"AAAA"));
        for index in 0..6 {
            reads.push_back(downsampling_read(
                &format!("excluded_{}", index),
                b"TTTT",
            ));
        }

        let result = create_consensus_sam_read(
            &reference_manager,
            &layout,
            &2,
            &reads,
            &AffineScoring::default_dna(),
            &MergeStrategy::Stretcher,
        )
        .unwrap();

        assert_eq!(result.read.aligned_read.read_aligned, b"AAAA");
        assert_eq!(result.added_tags.get(&[b'r', b'c']).unwrap(), "8");
        assert_eq!(result.added_tags.get(&[b'd', b'c']).unwrap(), "2");
        assert_eq!(
            result.added_tags.get(&[b'a', b'r']).unwrap(),
            "included_1,included_2"
        );
    }

    #[test]
    fn test_zero_downsampling_limit_uses_all_reads() {
        let reads = VecDeque::from(vec![
            downsampling_read("read_1", b"AAAA"),
            downsampling_read("read_2", b"AAAA"),
        ]);

        assert_eq!(select_reads_for_consensus(&reads, 0).len(), 2);
    }

    #[cfg(feature = "spoa")]
    fn create_test_alignment_result(
        reference_name: String,
        read_name: String,
        reference_aligned: Vec<u8>,
        read_aligned: Vec<u8>,
        read_quals: Option<Vec<u8>>,
        score: f64,
    ) -> AlignmentResult {
        AlignmentResult {
            reference_name,
            read_name,
            reference_aligned,
            read_aligned,
            read_quals,
            cigar_string: vec![],
            path: vec![],
            score,
            reference_start: 0,
            read_start: 0,
            bounding_box: None,
        }
    }

    #[cfg(feature = "spoa")]
    fn create_test_sorting_read(alignment: AlignmentResult) -> SortingReadSetContainer {
        SortingReadSetContainer {
            ordered_sorting_keys: vec![('*', CorrectedKey::new('*', vec![FASTA_A, FASTA_A], vec![FASTA_A, FASTA_A]))],
            ordered_unsorted_keys: VecDeque::new(),
            aligned_read: alignment,
        }
    }

    #[test]
    #[cfg(feature = "spoa")] // Pre-existing: references removed SPOA function
    fn test_create_poa_consensus_basic() {
        let mut reads = VecDeque::new();

        let alignment1 = create_test_alignment_result(
            "test_ref".to_string(),
            "read1".to_string(),
            b"ACGTACGTACGTACGT".to_vec(),
            b"ACGTACGTACGTACGT".to_vec(),
            Some(vec![b'I' - PHRED_OFFSET; 16]),
            100.0,
        );

        let alignment2 = create_test_alignment_result(
            "test_ref".to_string(),
            "read2".to_string(),
            b"ACGTACGTACGTACGT".to_vec(),
            b"ACGTACGTACGTACGT".to_vec(),
            Some(vec![b'I' - PHRED_OFFSET; 16]),
            100.0,
        );

        reads.push_back(create_test_sorting_read(alignment1));
        reads.push_back(create_test_sorting_read(alignment2));

        let max_reads = 10;
        let result = create_poa_consensus(&reads, &max_reads);

        assert!(!result.0.is_empty());
        assert!(!result.1.is_empty());
    }

    #[test]
    #[cfg(feature = "spoa")] // Pre-existing: references removed SPOA function
    fn test_create_poa_consensus_with_gaps() {
        let mut reads = VecDeque::new();

        let alignment1 = create_test_alignment_result(
            "test_ref".to_string(),
            "read1".to_string(),
            b"ACGTACGTACGTACGT".to_vec(),
            b"ACGTACGTACGTACGT".to_vec(),
            Some(vec![b'I' - PHRED_OFFSET; 16]),
            100.0,
        );

        let alignment2 = create_test_alignment_result(
            "test_ref".to_string(),
            "read2".to_string(),
            b"ACGTACGTACGTACGT".to_vec(),
            b"ACGT----ACGTACGT".to_vec(),
            Some(vec![b'I' - PHRED_OFFSET; 12]),
            90.0,
        );

        reads.push_back(create_test_sorting_read(alignment1));
        reads.push_back(create_test_sorting_read(alignment2));

        let max_reads = 10;
        let result = create_poa_consensus(&reads, &max_reads);

        assert!(!result.0.is_empty());
        assert!(!result.1.is_empty());
    }

    #[test]
    #[cfg(feature = "spoa")] // Pre-existing: references removed SPOA function
    fn test_create_poa_consensus_with_mismatches() {
        let mut reads = VecDeque::new();

        let alignment1 = create_test_alignment_result(
            "test_ref".to_string(),
            "read1".to_string(),
            b"ACGTACGTACGTACGT".to_vec(),
            b"ACGTACGTACGTACGT".to_vec(),
            Some(vec![b'I' - PHRED_OFFSET; 16]),
            100.0,
        );

        let alignment2 = create_test_alignment_result(
            "test_ref".to_string(),
            "read2".to_string(),
            b"ACGTACGTACGTACGT".to_vec(),
            b"TCGTACGTACGTACGA".to_vec(), // Mismatches at positions 0 and 15
            Some(vec![b'I' - PHRED_OFFSET; 16]),
            90.0,
        );

        let alignment3 = create_test_alignment_result(
            "test_ref".to_string(),
            "read3".to_string(),
            b"ACGTACGTACGTACGT".to_vec(),
            b"ACGTACGTACGTACGT".to_vec(),
            Some(vec![b'I' - PHRED_OFFSET; 16]),
            100.0,
        );

        reads.push_back(create_test_sorting_read(alignment1));
        reads.push_back(create_test_sorting_read(alignment2));
        reads.push_back(create_test_sorting_read(alignment3));

        let max_reads = 10;
        let result = create_poa_consensus(&reads, &max_reads);

        assert!(!result.0.is_empty());
        assert!(!result.1.is_empty());
        // Consensus should favor the majority sequence (2 vs 1)
    }

    #[test]
    #[cfg(feature = "spoa")] // Pre-existing: references removed SPOA function
    fn test_create_poa_consensus_with_downsampling() {
        let mut reads = VecDeque::new();

        // Create more reads than the downsample limit
        for i in 0..15 {
            let alignment = create_test_alignment_result(
                "test_ref".to_string(),
                format!("read{}", i),
                b"ACGTACGTACGTACGT".to_vec(),
                b"ACGTACGTACGTACGT".to_vec(),
                Some(vec![b'I' - PHRED_OFFSET; 16]),
                100.0,
            );
            reads.push_back(create_test_sorting_read(alignment));
        }

        let max_reads = 5; // Downsample to 5 reads
        let result = create_poa_consensus(&reads, &max_reads);

        assert!(!result.0.is_empty());
        assert!(!result.1.is_empty());
        // The algorithm should only use the first 5 reads due to downsampling
    }

    #[test]
    #[cfg(feature = "spoa")] // Pre-existing: references removed SPOA function
    fn test_create_poa_consensus_single_read() {
        let mut reads = VecDeque::new();

        let alignment = create_test_alignment_result(
            "test_ref".to_string(),
            "read1".to_string(),
            b"ACGTACGTACGTACGT".to_vec(),
            b"ACGTACGTACGTACGT".to_vec(),
            Some(vec![b'I' - PHRED_OFFSET; 16]),
            100.0,
        );

        reads.push_back(create_test_sorting_read(alignment));

        let max_reads = 10;
        let result = create_poa_consensus(&reads, &max_reads);

        assert!(!result.0.is_empty());
        assert!(!result.1.is_empty());
        // Should handle single read case
    }

    #[test]
    #[cfg(feature = "spoa")] // Pre-existing: references removed SPOA function
    fn test_create_poa_consensus_complex_gaps() {
        let mut reads = VecDeque::new();

        let alignment1 = create_test_alignment_result(
            "test_ref".to_string(),
            "read1".to_string(),
            b"ACGTACGTACGTACGT".to_vec(),
            b"ACGTACGTACGTACGT".to_vec(),
            Some(vec![b'I' - PHRED_OFFSET; 16]),
            100.0,
        );

        let alignment2 = create_test_alignment_result(
            "test_ref".to_string(),
            "read2".to_string(),
            b"ACGTACGTACGTACGT".to_vec(),
            b"AC--AC--ACGTACGT".to_vec(), // Multiple gaps
            Some(vec![b'I' - PHRED_OFFSET; 12]),
            90.0,
        );

        let alignment3 = create_test_alignment_result(
            "test_ref".to_string(),
            "read3".to_string(),
            b"ACGTACGTACGTACGT".to_vec(),
            b"ACGTACGTAC--ACGT".to_vec(), // Different gap pattern
            Some(vec![b'I' - PHRED_OFFSET; 14]),
            95.0,
        );

        reads.push_back(create_test_sorting_read(alignment1));
        reads.push_back(create_test_sorting_read(alignment2));
        reads.push_back(create_test_sorting_read(alignment3));

        let max_reads = 10;
        let result = create_poa_consensus(&reads, &max_reads);

        assert!(!result.0.is_empty());
        assert!(!result.1.is_empty());
        // Should handle complex gap patterns
    }

    #[test]
    #[cfg(feature = "spoa")] // Pre-existing: references removed SPOA function
    fn test_poa_consensus_direct() {
        let quals = vec![
            vec![b'I' - PHRED_OFFSET; 8],
            vec![b'I' - PHRED_OFFSET; 8]];

        let reference = "ACGTACGT\0".as_bytes().to_vec();
        let read1 = "ACGTACGT\0".as_bytes().to_vec();
        let read2 = "ACGTACGT\0".as_bytes().to_vec();
        let vec_of_reads = vec![reference, read1, read2];
        let result = poa_consensus(&vec_of_reads, &quals);
        println!("{:?}", u8s(&result.0));
        assert_eq!(result.0, "ACGTACGT".as_bytes().to_vec());
        assert_eq!(result.1.len(), 8);
        assert!(result.1.iter().all(|&q| q > 0)); // All quality scores should be positive
    }

    #[test]
    #[cfg(feature = "spoa")] // Pre-existing: references removed SPOA function
    fn test_poa_consensus_with_disagreement() {
        let quals = vec![
            vec![b'I' - PHRED_OFFSET; 8],
            vec![b'I' - PHRED_OFFSET; 8]];

        let reference = "ACGTACGT\0".as_bytes().to_vec();
        let read1 = "TCGTACGA\0".as_bytes().to_vec(); // Different first and last base
        let read2 = "ACGTACGT\0".as_bytes().to_vec();
        let vec_of_reads = vec![reference, read1, read2];
        let result = poa_consensus(&vec_of_reads, &quals);

        assert!(!result.0.is_empty());
        assert_eq!(result.1.len(), result.0.len());
        // Should handle disagreements and pick consensus
    }

    #[test]
    fn test_get_reference_alignment_rate() {
        let reference = b"ACGTACGT";
        let read_perfect = b"ACGTACGT";
        let rate_perfect = get_reference_alignment_rate(reference, read_perfect);
        assert_eq!(rate_perfect, 1.0);

        let read_half_match = b"ACGTTTTT";
        let rate_half = get_reference_alignment_rate(reference, read_half_match);
        assert_eq!(rate_half, 0.625); // 5 matches out of 8: positions 0,1,2,3,7

        let read_some_match = b"AAAAAAAA"; // Matches at positions 0,4: A vs A
        let rate_some = get_reference_alignment_rate(reference, read_some_match);
        assert_eq!(rate_some, 0.25); // 2 matches out of 8

        // Test with gaps and Ns - gaps are skipped in calculation
        let reference_with_gap = b"ACG-TACGT";
        let read_with_gap = b"ACG-TACGT";
        let rate_gap = get_reference_alignment_rate(reference_with_gap, read_with_gap);
        assert_eq!(rate_gap, 1.0);

        // Test mixed case
        let reference_mixed = b"ACGTACGT";
        let read_mixed = b"ACGTTTCG"; // 4 matches out of 8: positions 0,1,2,3
        let rate_mixed = get_reference_alignment_rate(reference_mixed, read_mixed);
        assert_eq!(rate_mixed, 0.5);
    }

    #[test]
    #[cfg(feature = "spoa")] // Pre-existing: references removed SPOA function
    fn test_consensus_with_downsampling() {
        let mut reads = VecDeque::new();

        // Create more reads than the downsample limit
        for i in 0..15 {
            let alignment = create_test_alignment_result(
                "test_ref".to_string(),
                format!("read{}", i),
                b"ACGTACGTACGTACGT".to_vec(),
                b"ACGTACGTACGTACGT".to_vec(),
                Some(vec![b'I' - PHRED_OFFSET; 16]),
                100.0,
            );
            reads.push_back(create_test_sorting_read(alignment));
        }

        let max_reads = 5; // Downsample to 5 reads
        let result = create_poa_consensus(&reads, &max_reads);

        assert!(!result.0.is_empty());
        assert!(!result.1.is_empty());
        // The algorithm should only use the first 5 reads due to downsampling
    }

    #[test]
    #[cfg(feature = "spoa")] // Pre-existing: references removed SPOA function
    fn test_consensus_with_mismatches() {
        let mut reads = VecDeque::new();

        let alignment1 = create_test_alignment_result(
            "test_ref".to_string(),
            "read1".to_string(),
            b"ACGTACGTACGTACGT".to_vec(),
            b"ACGTACGTACGTACGT".to_vec(),
            Some(vec![b'I' - PHRED_OFFSET; 16]),
            100.0,
        );

        let alignment2 = create_test_alignment_result(
            "test_ref".to_string(),
            "read2".to_string(),
            b"ACGTACGTACGTACGT".to_vec(),
            b"TCGTACGTACGTACGA".to_vec(), // Mismatches at positions 0 and 15
            Some(vec![b'I' - PHRED_OFFSET; 16]),
            90.0,
        );

        let alignment3 = create_test_alignment_result(
            "test_ref".to_string(),
            "read3".to_string(),
            b"ACGTACGTACGTACGT".to_vec(),
            b"ACGTACGTACGTACGT".to_vec(),
            Some(vec![b'I' - PHRED_OFFSET; 16]),
            100.0,
        );

        reads.push_back(create_test_sorting_read(alignment1));
        reads.push_back(create_test_sorting_read(alignment2));
        reads.push_back(create_test_sorting_read(alignment3));

        let max_reads = 10;
        let result = create_poa_consensus(&reads, &max_reads);

        assert!(!result.0.is_empty());
        assert!(!result.1.is_empty());
        // Consensus should favor the majority sequence (2 vs 1)
    }

    #[test]
    fn test_cigar_string() {
        let reference =
            &"CGTACGCTAGACATTGTGCCGCATCGATTGTAGTGACAATAGGAAA-------TATACAAG".as_bytes().to_vec();
        let read = &"CGT-----AGACATTGTGCCGCATCGATTGTAGTGACAATAGGAAATGACGGCTATACAAG".as_bytes().to_vec();
        let cigar = reference_read_to_cigar_string(&reference, &read);
        let true_cigar = CigarString(vec![
            Cigar::Match(3),
            Cigar::Del(5),
            Cigar::Match(38),
            Cigar::Ins(7),
            Cigar::Match(8),
        ]);
        assert_eq!(cigar, true_cigar);
    }

    #[test]
    #[cfg(feature = "spoa")] // Pre-existing: references removed SPOA function
    fn test_consensus_string() {
        let quals = vec![
            vec![b'I' - PHRED_OFFSET; 8],
            vec![b'I' - PHRED_OFFSET; 8],
            vec![b'I' - PHRED_OFFSET; 7]];

        let reference = "ACGTACGT\0".as_bytes().to_vec();
        let read1 = "ACGTACGT\0".as_bytes().to_vec();
        let read2 = "ACGTACGT\0".as_bytes().to_vec();
        let read3 = "ACGTAC-T\0".as_bytes().to_vec();
        let vec_of_reads = vec![reference, read1, read2, read3];
        let result = poa_consensus(&vec_of_reads, &quals);
        assert_eq!(result.0, "ACGTACGT".as_bytes().to_vec());

        let quals = vec![
            vec![b'I' - PHRED_OFFSET; 8],
            vec![b'I' - PHRED_OFFSET; 7],
            vec![b'I' - PHRED_OFFSET; 7]];

        let reference = "ACGTACGT\0".as_bytes().to_vec();
        let read1 = "ACGTAC-T\0".as_bytes().to_vec();
        let read2 = "ACGTAC-T\0".as_bytes().to_vec();
        let read3 = "ACGTAC-T\0".as_bytes().to_vec();
        let vec_of_reads = vec![reference, read1, read2, read3];
        let result = poa_consensus(&vec_of_reads, &quals);
        assert_eq!(result.0, "ACGTAC-T".as_bytes().to_vec());

        let quals = vec![
            vec![b'I' - PHRED_OFFSET; 8],
            vec![b'I' - PHRED_OFFSET; 12],
            vec![b'I' - PHRED_OFFSET; 7]];

        let reference = "ACGTACGT\0".as_bytes().to_vec();
        let read1 = "ACGTACGT\0".as_bytes().to_vec();       //      ACGTACGT
        let read2 = "AAAAAACGTAC-T\0".as_bytes().to_vec();  // AAAAAACGTAC-T
        let read3 = "ACGTAC-T\0".as_bytes().to_vec();       //      ACGTAC-T
        let vec_of_reads = vec![reference, read1, read2, read3];
        let result = poa_consensus(&vec_of_reads, &quals);
        println!("result {}", u8s(&result.0));
        assert_eq!(result.0, "-----ACGTACGT".as_bytes().to_vec());

        let quals = vec![
            vec![b'5' - PHRED_OFFSET; 11],
            vec![b'5' - PHRED_OFFSET; 15],
            vec![b'5' - PHRED_OFFSET; 7]];

        let reference = "ACGTACGT\0".as_bytes().to_vec();
        let read1 = "ACGTACGTTTT\0".as_bytes().to_vec();      //      ACGTACGTTTT
        let read2 = "AAAAAACGTACTTTT\0".as_bytes().to_vec();  // AAAAAACGTAC-TTTT
        let read3 = "ACGTACT\0".as_bytes().to_vec();          //      ACGTAC-T
        let vec_of_reads = vec![reference, read1, read2, read3];
        let result = poa_consensus(&vec_of_reads, &quals);
        println!("result {}", u8s(&result.0));

        assert_eq!(result.0, "-----ACGTACGTTTT".as_bytes().to_vec());
        //assert_eq!(result.1, [50, 50, 50, 50, 50, 90, 90, 90, 90, 90, 90, 4, 3, 3, 3, 90]); // we max out at Q40
        assert_eq!(result.1, [40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 4, 3, 3, 3, 40]); // we max out at Q40

        // "     ACGTACGTTTT\0".as_bytes().to_vec();
        // "AAAAAACGTAC TTTT\0".as_bytes().to_vec();
        // "     ACGTAC T\0".as_bytes().to_vec();*/
    }

    #[test]
    fn test_phred_to_prob() {
        assert_eq!(0.0001, phred_to_error_prob(&(b'I' - 33)));
        assert_eq!(1.0, phred_to_error_prob(&(b'!' - 33)));
        assert_eq!(0.1, phred_to_error_prob(&(b'+' - 33)));
    }

    //fn within_delta(x: &f64,y: &f64, delta: f64) -> bool {

    //}
    #[test]
    fn test_combine_qual_scores() {
        let bases = Vec::from(&[b'A', b'A', b'A', b'A']);
        let quals = Vec::from(&[b'I' - PHRED_OFFSET, b'I' - PHRED_OFFSET, b'I' - PHRED_OFFSET, b'I' - PHRED_OFFSET]);

        // rounding to 1
        assert_eq!(1.0, combine_qual_scores(vec![bases.as_slice()].as_slice(), vec![quals.as_slice()].as_slice(), &b'A', &0.1_f64)[0]); // we're ~ 1.0, fully confident it's an 'A'

        let bases = Vec::from(&[b'A', b'C', b'G', b'T']);

        // recover the priors
        let qual_combined = combine_qual_scores(vec![bases.as_slice()].as_slice(), vec![quals.as_slice()].as_slice(), &b'A', &0.99_f64);
        //println!("qual combined {} {} {} {} {}", qual_combined[0], qual_combined[1], qual_combined[2], qual_combined[3], qual_combined[4]);

        assert!((0.9924811371413187 - qual_combined[0]).abs() < 0.0001);
    }

    #[test]
    fn test_phred_to_error_prob_high_quality() {
        // Phred 40 => error prob ~0.0001
        let prob = phred_to_error_prob(&40);
        assert!((prob - 0.0001).abs() < 0.00001);
    }

    #[test]
    fn test_phred_to_error_prob_low_quality() {
        // Phred 10 => error prob 0.1
        let prob = phred_to_error_prob(&10);
        assert!((prob - 0.1).abs() < 0.001);
    }

    #[test]
    fn test_phred_to_error_prob_zero() {
        // Phred 0 => error prob 1.0
        let prob = phred_to_error_prob(&0);
        assert!((prob - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_prob_to_phred_high_confidence() {
        // prob close to 1.0 => high phred (high quality): -10*log10(1-0.9999) = 40
        let phred = prob_to_phred(&0.9999);
        assert_eq!(phred, 40);
    }

    #[test]
    fn test_prob_to_phred_low_confidence() {
        // prob = 0.5 => phred ~3
        let phred = prob_to_phred(&0.5);
        assert_eq!(phred, 3);
    }

    #[test]
    fn test_prob_to_phred_nan() {
        let phred = prob_to_phred(&f64::NAN);
        assert_eq!(phred, 0);
    }

    #[test]
    fn test_prob_to_phred_very_small() {
        let phred = prob_to_phred(&0.000000001);
        assert_eq!(phred, 0);
    }

    #[test]
    fn test_prob_to_phred_capped_at_40() {
        // prob close to 0 should be capped at 40
        let phred = prob_to_phred(&0.001);
        assert!(phred <= 40);
    }

    #[test]
    fn test_calculate_qual_scores_uniform() {
        // All equal log probs => should produce uniform distribution
        let mut props = [0.0_f64; 5];
        let result = calculate_qual_scores(&mut props);
        assert!((result[0] - 0.2).abs() < 0.001);
        assert!((result[1] - 0.2).abs() < 0.001);
        assert!((result[2] - 0.2).abs() < 0.001);
        assert!((result[3] - 0.2).abs() < 0.001);
        assert!((result[4] - 0.2).abs() < 0.001);
    }

    #[test]
    fn test_calculate_qual_scores_sum_to_one() {
        let mut props = [-1.0, -2.0, -3.0, -4.0, -5.0];
        let result = calculate_qual_scores(&mut props);
        let sum: f64 = result.iter().sum();
        assert!((sum - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_calculate_qual_scores_dominant() {
        // One very high, rest very low
        let mut props = [0.0, -100.0, -100.0, -100.0, -100.0];
        let result = calculate_qual_scores(&mut props);
        assert!(result[0] > 0.99);
    }

    #[test]
    fn test_get_reference_alignment_rate_all_match() {
        let reference = b"ACGTACGT";
        let read = b"ACGTACGT";
        assert_eq!(get_reference_alignment_rate(reference, read), 1.0);
    }

    #[test]
    fn test_get_reference_alignment_rate_no_match() {
        let reference = b"AAAA";
        let read = b"TTTT";
        assert_eq!(get_reference_alignment_rate(reference, read), 0.0);
    }

    #[test]
    fn test_get_reference_alignment_rate_with_gaps() {
        // Gaps (b'-' = 45) are below ASCII 64, so they're skipped
        let reference = b"A-A";
        let read = b"A-A";
        assert_eq!(get_reference_alignment_rate(reference, read), 1.0);
    }

    #[test]
    fn test_combine_qual_scores_all_same_base() {
        let bases = vec![b'A', b'A', b'A'];
        let quals = vec![30, 30, 30];
        let result = combine_qual_scores(
            vec![bases.as_slice()].as_slice(),
            vec![quals.as_slice()].as_slice(),
            &b'A',
            &0.75,
        );
        // A should dominate
        assert!(result[0] > result[1]);
        assert!(result[0] > result[2]);
        assert!(result[0] > result[3]);
    }

    #[test]
    fn test_combine_qual_scores_all_different_bases() {
        let bases = vec![b'A', b'C', b'G', b'T'];
        let quals = vec![30, 30, 30, 30];
        let result = combine_qual_scores(
            vec![bases.as_slice()].as_slice(),
            vec![quals.as_slice()].as_slice(),
            &b'N',
            &0.25, // uniform prior
        );
        // With uniform prior and one observation each, should be roughly equal
        let sum: f64 = result[0..4].iter().sum();
        assert!((sum - 1.0).abs() < 0.01 || result[4] < 0.01);
    }
}
