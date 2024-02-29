extern crate spoa;

use crate::alignment::alignment_matrix::{
    create_scoring_record_3d, Alignment, AlignmentTag, AlignmentType as LocalAlignmentType,
};
use crate::alignment::fasta_bit_encoding::{FastaBase, FASTA_N, FASTA_UNSET};
use crate::alignment::scoring_functions::AffineScoring;
use crate::alignment_manager::{align_two_strings, simplify_cigar_string, OutputAlignmentWriter};
use crate::read_strategies::read_disk_sorter::SortingReadSetContainer;
use crate::reference::fasta_reference::ReferenceManager;
use counter::Counter;
use ndarray::Ix3;
use rand::prelude::*;
use rust_htslib::bam::record::CigarString;
use shardio::{Range, ShardReader};
use spoa::{AlignmentEngine, AlignmentType, Graph};
use std::cmp;
use std::cmp::Ordering;
use std::collections::{HashMap, VecDeque};
use std::ffi::CString;
use std::sync::{Arc, Mutex};
use indicatif::style::ProgressTracker;

pub fn write_consensus_reads(
    reader: &ShardReader<SortingReadSetContainer>,
    writer: &mut dyn OutputAlignmentWriter,
    levels: usize,
    reference_manager: &ReferenceManager,
    maximum_reads_before_downsampling: &usize,
) {
    let mut last_read: Option<SortingReadSetContainer> = None;
    let mut buffered_reads = VecDeque::new();

    let mut written_buffers = 0;
    let mut processed_reads = 0;

    let score = AffineScoring::default_DNA();


    let arc_output = Arc::new(Mutex::new(writer));

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(8)
        .build()
        .unwrap();

    pool.scope(|s| {
        let mut handled_reads = 0;
        reader.iter_range(&Range::all()).unwrap().for_each(|x| {
            let x = x.unwrap();
            assert_eq!(x.ordered_sorting_keys.len(), levels);
            if !(last_read.is_some() && x.cmp(last_read.as_ref().unwrap()) == Ordering::Equal)
                && !buffered_reads.is_empty()
            {
                let my_buffered_reads = buffered_reads.clone();
                buffered_reads = VecDeque::new();
                s.spawn( |y| {
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
                        &AffineScoring::default_DNA(),
                        &mut alignment_mat,
                    );

                    let arc_writer = arc_output.clone();
                    let mut arc_writer = arc_writer.lock().expect("Unable to access multi-threaded writer");
                    arc_writer.write_read(&new_read.read, &new_read.added_tags);
                });
                written_buffers += 1;
            }
            processed_reads += 1;
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
        );
        let arc_writer = arc_output.clone();
        let mut arc_writer = arc_writer.lock().expect("Unable to access multi-threaded writer");
        arc_writer.write_read(&new_read.read, &new_read.added_tags);
        written_buffers += 1;
    }
    info!(
        "Processed {} reads into {} consensus reads",
        processed_reads, written_buffers
    );
}

pub struct SamReadyOutput {
    read: SortingReadSetContainer,
    added_tags: HashMap<(u8, u8), String>,
}

fn create_sam_read(
    reference_manager: &ReferenceManager,
    maximum_reads_before_downsampling: &usize,
    buffered_reads: &VecDeque<SortingReadSetContainer>,
    my_aff_score: &AffineScoring,
    mut alignment_mat: &mut Alignment<Ix3>,
) -> SamReadyOutput {
    let mut added_tags = HashMap::new();
    added_tags.insert((b'r', b'c'), buffered_reads.len().to_string());
    added_tags.insert(
        (b'd', b'c'),
        cmp::min(*maximum_reads_before_downsampling, buffered_reads.len()).to_string(),
    );

    if buffered_reads.len() > 1 {
        let consensus_reads =
            create_poa_consensus(buffered_reads, maximum_reads_before_downsampling);
        let consensus_reference = Counter::<Vec<u8>, usize>::init(
            buffered_reads
                .iter()
                .map(|x| x.aligned_read.reference_name.clone().as_bytes().to_vec())
                .collect::<Vec<Vec<u8>>>(),
        )
            .most_common_ordered();

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

        let new_alignment = align_two_strings(
            &reference_pointer.sequence,
            &FastaBase::from_vec_u8(&consensus_reads),
            my_aff_score,
            false,
            &String::from_utf8(reference_pointer.name.clone()).unwrap(),
            &read_name,
            None,
        );
        //println!("New alignment: \n{}\n{}\n{:?}", FastaBase::string(&new_alignment.read_aligned),FastaBase::string(&new_alignment.reference_aligned),simplify_cigar_string(&new_alignment.cigar_string));
        let read_names = buffered_reads
            .iter()
            .map(|x| x.aligned_read.read_name.clone())
            .collect::<Vec<String>>();

        added_tags.insert((b'a', b'r'), read_names.join(","));
        added_tags.insert(
            (b'r', b'm'),
            get_reference_alignment_rate(
                &new_alignment.reference_aligned,
                &new_alignment.read_aligned,
            )
                .to_string(),
        );

        added_tags.insert((b'a', b's'), new_alignment.score.to_string());
        let new_sorting_read = buffered_reads
            .get(0)
            .unwrap()
            .with_new_alignment(new_alignment);

        SamReadyOutput { read: new_sorting_read, added_tags }
    } else {
        let single_read = buffered_reads.get(0).unwrap().clone();
        added_tags.insert((b'a', b'r'), single_read.aligned_read.read_name.clone());
        added_tags.insert(
            (b'r', b'm'),
            get_reference_alignment_rate(
                &single_read.aligned_read.reference_aligned,
                &single_read.aligned_read.read_aligned,
            )
                .to_string(),
        );
        added_tags.insert((b'a', b's'), single_read.aligned_read.score.to_string());
        SamReadyOutput { read: single_read, added_tags }
    }
}

pub fn get_reference_alignment_rate(reference: &[FastaBase], read: &[FastaBase]) -> f64 {
    let mut matches = 0;
    let mut mismatches = 0;

    //println!("reference: {} read: {}", FastaBase::string(&reference),FastaBase::string(&read));
    for index in 0..reference.len() {
        if reference.get(index).unwrap() != &FASTA_UNSET
            && !reference.get(index).unwrap().strict_identity(&FASTA_N)
            && read.get(index).unwrap() != &FASTA_UNSET
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

pub fn reference_read_to_cigar_string(
    reference_seq: &Vec<FastaBase>,
    read_seq: &[FastaBase],
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
) -> Vec<u8> {
    let mut base_sequences = sequences
        .iter()
        .map(|n| {
            let mut y = FastaBase::vec_u8(&FastaBase::strip_gaps(&n.aligned_read.read_aligned));
            y.push(b'\0');
            y
        })
        .collect::<Vec<Vec<u8>>>();

    // downsample if needed -- it's not the best idea to do downsampling here after all the work above,
    // but sometimes the borrow checker is a pain when your brain is small
    if base_sequences.len() > *downsample_to {
        let mut rng = thread_rng();
        base_sequences = base_sequences
            .into_iter()
            .choose_multiple(&mut rng, *downsample_to);
    }
    poa_consensus(base_sequences)
}

fn poa_consensus(base_sequences: Vec<Vec<u8>>) -> Vec<u8> {
    let mut eng = AlignmentEngine::new(AlignmentType::kNW, 5, -4, -3, -1, -3, -1);
    let mut graph = Graph::new();

    for seq in &base_sequences {
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

    graph.consensus().to_str().unwrap().to_owned().into_bytes()
}

// reference_read_to_cigar_string

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::Cigar;

    #[test]
    fn test_cigar_string() {
        let reference = FastaBase::from_string(
            &"CGTACGCTAGACATTGTGCCGCATCGATTGTAGTGACAATAGGAAA-------TATACAAG".to_string(),
        );
        let read = FastaBase::from_string(
            &"CGT-----AGACATTGTGCCGCATCGATTGTAGTGACAATAGGAAATGACGGCTATACAAG".to_string(),
        );
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
        let read1 = "ACGTACGT\0".as_bytes().to_vec();
        let read2 = "ACGTACGT\0".as_bytes().to_vec();
        let read3 = "ACGTAC-T\0".as_bytes().to_vec();
        let vec_of_reads = vec![read1, read2, read3];
        let result = poa_consensus(vec_of_reads);
        assert_eq!(result, "ACGTACGT".as_bytes().to_vec());

        let read1 = "ACGTACGT\0".as_bytes().to_vec();
        let read2 = "ACGTAC-T\0".as_bytes().to_vec();
        let read3 = "ACGTAC-T\0".as_bytes().to_vec();
        let vec_of_reads = vec![read1, read2, read3];
        let result = poa_consensus(vec_of_reads);
        assert_eq!(result, "ACGTAC-T".as_bytes().to_vec());

        let read1 = "ACGTACGT\0".as_bytes().to_vec();
        let read2 = "AAAAAACGTAC-T\0".as_bytes().to_vec();
        let read3 = "ACGTAC-T\0".as_bytes().to_vec();
        let vec_of_reads = vec![read1, read2, read3];
        let result = poa_consensus(vec_of_reads);
        assert_eq!(result, "AAAAAACGTAC-T".as_bytes().to_vec());

        let read1 = "ACGTACGTTTT\0".as_bytes().to_vec();
        let read2 = "AAAAAACGTAC-TTTT\0".as_bytes().to_vec();
        let read3 = "ACGTAC-T\0".as_bytes().to_vec();
        let vec_of_reads = vec![read1, read2, read3];
        let result = poa_consensus(vec_of_reads);
        assert_eq!(result, "AAAAAACGTAC-TTTT".as_bytes().to_vec());
    }
}
