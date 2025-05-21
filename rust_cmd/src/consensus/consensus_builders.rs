extern crate spoa;

use crate::alignment::alignment_matrix::{
    create_scoring_record_3d, Alignment, AlignmentTag, AlignmentType as LocalAlignmentType,
};
use crate::alignment::scoring_functions::AffineScoring;
use crate::alignment_manager::{align_two_strings, simplify_cigar_string, OutputAlignmentWriter};
use crate::read_strategies::read_disk_sorter::SortingReadSetContainer;
use crate::reference::fasta_reference::ReferenceManager;
use counter::Counter;
use ndarray::Ix3;
use rust_htslib::bam::record::CigarString;
use shardio::{Range, ShardReader};
use spoa::{AlignmentEngine, AlignmentType, Graph};
use std::cmp;
use std::cmp::Ordering;
use std::collections::{HashMap, VecDeque};
use std::convert::TryFrom;
use std::ffi::CString;
use std::sync::{Arc, Mutex};
use num_traits::{Pow, ToPrimitive};
use ::{FASTA_N, FASTA_UNSET};
use alignment::alignment_matrix::AlignmentResult;


pub enum MergeStrategy {
    STRICT_CONSENSUS,
    HYBRID,
    STRETCHER,
}

pub fn write_consensus_reads(
    reader: &ShardReader<SortingReadSetContainer>,
    writer: &mut dyn OutputAlignmentWriter,
    levels: usize,
    reference_manager: &ReferenceManager,
    maximum_reads_before_downsampling: &usize,
    merge_strategy: &MergeStrategy,
) {

    let mut last_read: Option<SortingReadSetContainer> = None;
    let mut buffered_reads = VecDeque::new();

    let processed_reads = Arc::new(Mutex::new(0));

    let score = AffineScoring::default_dna();

    let arc_output = Arc::new(Mutex::new(writer));

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(8)
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

                    // TODO fix to pooled approach like alignment
                    let mut alignment_mat: Alignment<Ix3> = create_scoring_record_3d(
                        (reference_manager.longest_ref + 1) * 2,
                        (reference_manager.longest_ref + 1) * 2,
                        LocalAlignmentType::Affine,
                        false,
                    );

                    let my_buffered_reads = my_buffered_reads;
                    let new_read = create_sam_read(
                        reference_manager,
                        maximum_reads_before_downsampling,
                        &my_buffered_reads,
                        &AffineScoring::default_dna(),
                        &mut alignment_mat,
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
                            if (*processed_reads as f64 / 100000.0).floor() - (current_proc_read as f64 / 1000.0).floor() >= 1.0 {
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
        let mut alignment_mat: Alignment<Ix3> = create_scoring_record_3d(
            (reference_manager.longest_ref + 1) * 2,
            (reference_manager.longest_ref + 1) * 2,
            LocalAlignmentType::Affine,
            false,
        );

        let new_read = create_sam_read(
            reference_manager,
            maximum_reads_before_downsampling,
            &buffered_reads,
            &score,
            &mut alignment_mat,
            merge_strategy
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
}

pub struct SamReadyOutput {
    pub read: SortingReadSetContainer,
    pub added_tags: HashMap<[u8; 2], String>,
}

fn create_sam_read(
    reference_manager: &ReferenceManager,
    maximum_reads_before_downsampling: &usize,
    buffered_reads: &VecDeque<SortingReadSetContainer>,
    my_aff_score: &AffineScoring,
    _alignment_mat: &mut Alignment<Ix3>,
    merge_strategy: &MergeStrategy,
) -> Option<SamReadyOutput> {
    let mut added_tags = HashMap::new();
    added_tags.insert([b'r', b'c'], buffered_reads.len().to_string());
    added_tags.insert(
        [b'd', b'c'],
        cmp::min(*maximum_reads_before_downsampling, buffered_reads.len()).to_string(),
    );

    if buffered_reads.len() > 1 {

        let consensus_reference = Counter::<Vec<u8>, usize>::init(
            buffered_reads
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

        let read_name = buffered_reads
            .iter()
            .next()
            .unwrap()
            .aligned_read
            .read_name
            .clone();

        let consensus_reads = match merge_strategy {
            MergeStrategy::STRICT_CONSENSUS => {
                let consensus_reads = create_poa_consensus(buffered_reads, maximum_reads_before_downsampling);

                let gapless_quals = consensus_reads.1.iter().enumerate().filter(|(i, _x)| consensus_reads.0.get(*i).unwrap() != &b'-').map(|(_i, x)| x.clone()).collect();

                let mut new_alignment = align_two_strings(
                    &reference_pointer.sequence,
                    &consensus_reads.0,
                    None, // TODO: fix with quality scores
                    my_aff_score,
                    false,
                    &String::from_utf8(reference_pointer.name.clone()).unwrap(),
                    &read_name,
                    None,
                );


                new_alignment.read_quals = Some(gapless_quals);
                Some(new_alignment)

            }
            MergeStrategy::HYBRID => {
                // try the stretcher first, if that fails out, try the POA
                let mut candidate = crate::consensus::stretcher::AlignmentCandidate::new(reference_pointer.sequence.as_slice(), reference_pointer.name.as_slice());

                let valid : usize = buffered_reads.iter().map(|x| {
                    match candidate.add_alignment(&x.aligned_read) {
                        Ok(_) => {0}
                        Err(_) => {1}
                    }
                }).sum();

                if valid > 1 {
                    let consensus_reads = create_poa_consensus(buffered_reads, maximum_reads_before_downsampling);

                    let gapless_quals = consensus_reads.1.iter().enumerate().filter(|(i, _x)| consensus_reads.0.get(*i).unwrap() != &b'-').map(|(_i, x)| x.clone()).collect();

                    let mut new_alignment = align_two_strings(
                        &reference_pointer.sequence,
                        &consensus_reads.0,
                        None, // TODO: fix with quality scores
                        my_aff_score,
                        false,
                        &String::from_utf8(reference_pointer.name.clone()).unwrap(),
                        &read_name,
                        None,
                    );


                    new_alignment.read_quals = Some(gapless_quals);
                    Some(new_alignment)
                } else {
                    Some(candidate.to_consensus(&0.75))
                }
            }
            MergeStrategy::STRETCHER => {
                let mut candidate = crate::consensus::stretcher::AlignmentCandidate::new(reference_pointer.sequence.as_slice(), reference_pointer.name.as_slice());

                let valid : usize = buffered_reads.iter().map(|x| {
                    match candidate.add_alignment(&x.aligned_read) {
                        Ok(_) => {0}
                        Err(_) => {1}
                    }
                }).sum();

                if valid > 1 {
                    None
                } else {
                    Some(candidate.to_consensus(&0.75))
                }

            }
        };

        match consensus_reads {
            None => {
                None
            }
            Some(con) => {
                //println!("New alignment: \n{}\n{}\n{:?}", FastaBase::string(&new_alignment.read_aligned),FastaBase::string(&new_alignment.reference_aligned),simplify_cigar_string(&new_alignment.cigar_string));
                let read_names = buffered_reads
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
                let new_sorting_read = buffered_reads
                    .get(0)
                    .unwrap()
                    .with_new_alignment(con);

                Some(SamReadyOutput { read: new_sorting_read, added_tags })
            }
        }

    } else {
        let single_read = buffered_reads.get(0).unwrap().clone();
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

pub fn create_poa_consensus(
    sequences: &VecDeque<SortingReadSetContainer>,
    downsample_to: &usize,
) -> (Vec<u8>, Vec<u8>) {
    let mut base_sequences = Vec::new();
    let mut quals_sequences = Vec::new();

    sequences
        .iter().take(*downsample_to)
        .for_each(|n| {
            let mut y = n.aligned_read.read_aligned.iter().filter(|x| **x != b'-').map(|x| *x).collect::<Vec<u8>>();
            y.push(b'\0');
            base_sequences.push(y);
            quals_sequences.push(n.aligned_read.read_quals.as_ref().unwrap().clone());
            // TODO we were panic'ing here --- why? oh right, we should move to stretcher
        });

    poa_consensus(&base_sequences, &quals_sequences)
}


fn poa_consensus(base_sequences: &Vec<Vec<u8>>, qual_sequences: &Vec<Vec<u8>>) -> (Vec<u8>, Vec<u8>) {
    let mut eng = AlignmentEngine::new(AlignmentType::kNW, 5, -4, -3, -1, -3, -1);
    let mut graph = Graph::new();

    for seq in base_sequences {
        assert!(!seq.is_empty());

        let cseq = CString::from_vec_with_nul(seq.clone()).unwrap_or_else(|_| {
            panic!(
                "CString::new failed from {}",
                String::from_utf8(seq.clone()).unwrap()
            )
        });
        let cqual = CString::new(vec![34u8; seq.len() - 1]).unwrap();

        let aln = eng.align(cseq.as_c_str(), &graph);
        graph.add_alignment(&aln, cseq.as_c_str(), &cqual);
    }

    let alignment = graph.multiple_sequence_alignment(false).into_iter().map(|x| x.to_str().unwrap().to_owned().into_bytes()).collect::<Vec<Vec<u8>>>();

    calculate_conc_qual_score(&alignment, qual_sequences)

    //graph.consensus().to_str().unwrap().to_owned().into_bytes()
}

pub fn calculate_conc_qual_score(alignments: &Vec<Vec<u8>>, quality_scores: &Vec<Vec<u8>>) -> (Vec<u8>, Vec<u8>) {
    // create a consensus as we go
    let mut conc = Vec::new();
    let mut final_quals = Vec::new();

    let mut sequence_indexes = vec![0_usize; alignments.len()];
    assert_eq!(alignments.len(), quality_scores.len());
    let ln = alignments.get(0).unwrap().len();

    (0..ln).for_each(|index| {
        let mut bases = Vec::new();
        let mut quals = Vec::new();
        alignments.iter().enumerate().for_each(|(sequence_index, x)| {
            assert_eq!(ln, x.len());

            let base = x.get(index).unwrap();
            let qual = quality_scores.get(sequence_index).unwrap();
            let qual = qual.get(*sequence_indexes.get(sequence_index).unwrap());
            let qual = qual.unwrap();

            sequence_indexes[sequence_index] = sequence_index + match base {
                b'-' => { 0 }
                _ => { 1 }
            };
            bases.push(*base);
            quals.push(*qual);
        });

        let qual_scores = combine_qual_scores(vec![bases.as_slice()].as_slice(), vec![quals.as_slice()].as_slice(), &0.75, &true);
        let index_of_max: usize = qual_scores
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.total_cmp(b))
            .map(|(index, _)| index).unwrap();

        if index_of_max < 4 {
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

pub fn phred_to_error_prob(phred: &u8, floor_zero_at_33: &bool) -> f64 {
    let phred =
        match (phred < &33_u8) && *floor_zero_at_33 {
            true => 33_u8,
            false => *phred,
        };


    assert!(phred >= 33 && phred <= 98, "{}", format!("Unable to format phred {}", phred)); // the upper bound is a bit arbitrary, but 33 + 93 = 126, the highest possible value reported by Illumina
    // TODO: we dont deal with phred + 64 format data -- at some point there will be legacy data that comes through; we should at least document this
    (10.0_f64).pow((phred.to_f64().unwrap() - 32.99999999999999999) / (-10.0)) // 32.9999 to avoid Inf powers
}

pub fn prob_to_phred(prob: &f64) -> u8 {
    // the upper bound is a bit arbitrary, but 33 + 93 = 126, the highest possible value reported by Illumina
    assert!(prob >= &0.0_f64 && prob <= &1.0_f64, "{}", format!("Unable to format prob {}", prob));
    if prob < &0.00000001_f64 {
        return 33_u8;
    }

    // TODO: we dont deal with phred + 64 format data
    let ret = 33.0 + ((-10.0) * (1.00000000001 - prob).log10()); // again to prevent zero getting in, we subtract from 1 + epsilon
    //println!("prob {} ret {}",prob, ret);
    assert!(ret >= 0.0_f64 && ret <= 256.0_f64, "{}", format!("Unable to format phred {}", ret));

    let ret = ret.round().to_u8().unwrap();
    let ret = if ret > 40 { // cap PHRED at 40; Noodles doesn't like higher TODO: fix this
        40_u8
    } else {
        ret as u8
    };
    assert!(ret >= 33_u8 && ret <= 40_u8);
    ret
}

pub(crate) fn combine_qual_scores(bases: &[&[u8]], scores: &[&[u8]], error_prior: &f64, phred_floor_at_33: &bool) -> [f64; 5] {
    // setup the priors
    let mut allele_props = [
        (error_prior / 3.0).log2(),
        (error_prior / 3.0).log2(),
        (error_prior / 3.0).log2(),
        (error_prior / 3.0).log2(),
        (error_prior / 3.0).log2()];

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
                b'-' => { 4 }
                _ => {
                    info!("unaccounted for quality score");
                    5
                }
            };
            if base_id < 5 {
                (0..5).for_each(|i| {
                    if i == base_id {
                        allele_props[i] = allele_props[i] + (1.0 - phred_to_error_prob(&qs, phred_floor_at_33)).log2();
                        //println!("MT id {} new prop {}", i, allele_props[i]);
                    } else {
                        allele_props[i] = allele_props[i] + (phred_to_error_prob(&qs, phred_floor_at_33) / 3.0_f64).log2();
                        //println!("NM id {} new prop {}", i, allele_props[i]);
                    }
                });
            }
        };
    });
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
    use rust_htslib::bam::record::Cigar;

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
        // CGTACGCTAGACATTGTGCCGCATCGATTGTAGTGACAATAGGAAATGACGGCTATACAAG-----------------------------------------------------------------------------------------------------------------------------------------------------------------GTAGGAGCCGTTACCAGGATGA------AGGTTATTAGGGGATCCGCTTTAAGGCCGGTCCTAGCAACAAGCTAACGGTGCAGGATCTTGGGTTTCTGTCTCTTATTCACATCTGACGCTGCCGACGACGAGGTATAAGTGTAGATATCGGTGGTCGCCGTATCATT
        // AAACCCAAGATCCTGCACCGTTAGCTTGCGTACGCTAGACATTGTGCCGCATCGATTGTAGTGACAATAGGAAATGACGGCTATACAAGGTAGGAGCCGTTACCAGGATGAAGGTTATTAGGGGATCCGCTTTAAGGCCGGTCCTAGCAACAAGCTAACGGTGCAGGATCTTGGGTTTCTGTCTCTTATTCACATCTGACGCTGCCGACGACGAGGTATAAGTGTAGATATCGGTGGTCGCCGTATCATTACAAAAGTGGGTGGGGGGGGGGGGGGGGC
        // CGTACGCTAGACATTGTGCCGCATC22222222222222TAGGAAATGACGGCTATACAAGGCATCGCGGTGTCTCGTCAATACACCTTACGGAGGCATTGGATGATAATGTCGCAAGGAGGTCTCAAGATTCTGTACCACACGTCGGCACGCGATTGAACCAATGGACAGAGGACAGGATACGTAGGATCACCAACTAGGTCATTAGGTGGAAGGTGATACGTAGGAGCCGTTACCAGGATGAACGATGAGGTTATTAGGGGATCCGCTTTAAGGCCGGTCCTAGCAANNNNNNNNNNNNNNNNNNNNNNNNNNNNCTGTCTCTTATACACATCTGACGCTGCCGACGANNNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT
    }

    #[test]
    fn test_consensus_string() {
        let quals = vec![vec![b'I'; 8], vec![b'I'; 8], vec![b'I'; 7]];

        let read1 = "ACGTACGT\0".as_bytes().to_vec();
        let read2 = "ACGTACGT\0".as_bytes().to_vec();
        let read3 = "ACGTAC-T\0".as_bytes().to_vec();
        let vec_of_reads = vec![read1, read2, read3];
        let result = poa_consensus(&vec_of_reads, &quals);
        assert_eq!(result.0, "ACGTACGT".as_bytes().to_vec());

        let quals = vec![vec![b'I'; 8], vec![b'I'; 7], vec![b'I'; 7]];

        let read1 = "ACGTACGT\0".as_bytes().to_vec();
        let read2 = "ACGTAC-T\0".as_bytes().to_vec();
        let read3 = "ACGTAC-T\0".as_bytes().to_vec();
        let vec_of_reads = vec![read1, read2, read3];
        let result = poa_consensus(&vec_of_reads, &quals);
        assert_eq!(result.0, "ACGTACT".as_bytes().to_vec());

        let quals = vec![vec![b'I'; 8], vec![b'I'; 12], vec![b'I'; 7]];

        let read1 = "ACGTACGT\0".as_bytes().to_vec();       //      ACGTACGT
        let read2 = "AAAAAACGTAC-T\0".as_bytes().to_vec();  // AAAAAACGTAC-T
        let read3 = "ACGTAC-T\0".as_bytes().to_vec();       //      ACGTAC-T
        let vec_of_reads = vec![read1, read2, read3];
        let result = poa_consensus(&vec_of_reads, &quals);
        assert_eq!(result.0, "ACGTACT".as_bytes().to_vec());

        let quals = vec![vec![b'5'; 11], vec![b'5'; 15], vec![b'5'; 7]];

        let read1 = "ACGTACGTTTT\0".as_bytes().to_vec();      //      ACGTACGTTTT
        let read2 = "AAAAAACGTACTTTT\0".as_bytes().to_vec();  // AAAAAACGTAC-TTTT
        let read3 = "ACGTACT\0".as_bytes().to_vec();          //      ACGTAC-T
        let vec_of_reads = vec![read1, read2, read3];
        let result = poa_consensus(&vec_of_reads, &quals);
        assert_eq!(result.0, "ACGTACTTTT".as_bytes().to_vec());
        assert_eq!(result.1, [40, 40, 40, 40, 40, 40, 40, 40, 40, 40]); // we max out at Q40

        // "     ACGTACGTTTT\0".as_bytes().to_vec();
        // "AAAAAACGTAC TTTT\0".as_bytes().to_vec();
        // "     ACGTAC T\0".as_bytes().to_vec();
    }

    #[test]
    fn test_phred_to_prob() {
        assert_eq!(0.0001, phred_to_error_prob(&b'I', &true));
        assert_eq!(1.0, phred_to_error_prob(&b'!', &true));
        assert_eq!(0.1, phred_to_error_prob(&b'+', &true));
    }

    //fn within_delta(x: &f64,y: &f64, delta: f64) -> bool {

    //}
    #[test]
    fn test_combine_qual_scores() {
        let bases = Vec::from(&[b'A', b'A', b'A', b'A']);
        let quals = Vec::from(&[b'I', b'I', b'I', b'I']);

        // rounding to 1
        assert_eq!(1.0, combine_qual_scores(vec![bases.as_slice()].as_slice(), vec![quals.as_slice()].as_slice(), &0.1_f64, &true)[0]); // we're ~ 1.0, fully confident it's an 'A'

        let bases = Vec::from(&[b'A', b'C', b'G', b'T']);

        // recover the priors
        let qual_combined = combine_qual_scores(vec![bases.as_slice()].as_slice(), vec![quals.as_slice()].as_slice(), &0.1_f64, &true);
        println!("qual combined {} {} {} {} {}", qual_combined[0], qual_combined[1], qual_combined[2], qual_combined[3], qual_combined[4]);

        assert!((0.25 - qual_combined[1]).abs() < 0.0001);
    }
}
