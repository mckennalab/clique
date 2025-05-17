use crate::consensus::consensus_builders::write_consensus_reads;
use crate::extractor::{extract_tag_sequences, extract_tagged_sequences, recover_soft_clipped_align_sequences, SoftClipResolution, stretch_sequence_to_alignment};
use crate::read_strategies::read_disk_sorter::SortingReadSetContainer;
use crate::read_strategies::sequence_layout::{ReferenceRecord, SequenceLayout, UMIConfiguration, UMISortType};
use crate::reference::fasta_reference::ReferenceManager;
use crate::umis::known_list::KnownList;

use crate::InstanceLivedTempDir;
use indicatif::ProgressBar;

use ::{FASTA_UNSET, noodles_bam as bam};
use noodles_bam::{bai, Record};
use noodles_sam::Header;

use shardio::{Range, ShardReader, ShardWriter};
use std::cmp::{min, Ordering};
use std::collections::{HashMap};
use std::path::PathBuf;
use itertools::Itertools;

use noodles_sam::alignment::record::QualityScores;
use petgraph::visit::Walker;
use FASTA_N;
use utils::read_utils::u8s;
use crate::alignment::alignment_matrix::{AlignmentResult, AlignmentTag};
use crate::alignment_manager::BamFileAlignmentWriter;
use crate::umis::degenerate_tags::DegenerateBuffer;

pub fn collapse(
    final_output: &String,
    temp_directory: &mut InstanceLivedTempDir,
    read_structure: &SequenceLayout,
    bam_file: &String,
) {
    // load up the reference files
    let rm = ReferenceManager::from_yaml_input(read_structure, 8, 4);

    // validate that each reference has the specified capture groups
    let validated_references = rm
        .references
        .iter()
        .map(|rf| {
            let reference_config = read_structure
                .references
                .get(String::from_utf8(rf.1.name.clone()).unwrap().as_str())
                .unwrap();
            SequenceLayout::validate_reference_sequence(
                &rf.1.sequence,
                &reference_config.umi_configurations,
            )
        })
        .all(|x| x == true);

    assert!(validated_references, "The reference sequences do not match the capture groups specified in the read structure file.");

    let mut known_level_lookups = get_known_level_lookups(read_structure);

    let mut writer = BamFileAlignmentWriter::new(&PathBuf::from(final_output), &rm);

    let _levels = 0;
    let mut read_count = 0;

    // for each reference, we fetch aligned reads, pull the sorting tags, and output the collapsed reads to a BAM file
    rm.references.iter().for_each(|(_id, reference)| {
        let ref_name = String::from_utf8(reference.name.clone()).unwrap();
        info!("processing reads from input BAM file: {}", bam_file);

        let sorted_reads_option=
            sort_reads_from_bam_file(bam_file, &ref_name, &rm, read_structure, temp_directory);
        read_count = sorted_reads_option.read_stats.passing_reads();

        let mut levels = 0;

        match sorted_reads_option.bam {
            None => {
                warn!("No valid reads found for reference {}", ref_name);
            }
            Some(mut sorted_reads) => {

                read_structure
                    .get_sorted_umi_configurations(&ref_name)
                    .iter()
                    .for_each(|tag| {

                        match tag.sort_type {
                            UMISortType::KnownTag => {
                                let ret = sort_known_level(
                                    temp_directory,
                                    &sorted_reads,
                                    &tag,
                                    &read_count,
                                    &mut known_level_lookups,
                                );
                                sorted_reads = ret.1;
                                read_count = ret.0;
                            }

                            UMISortType::DegenerateTag => {
                                let ret = sort_degenerate_level(
                                    temp_directory,
                                    &sorted_reads,
                                    &tag,
                                    &levels,
                                    &read_count,
                                );
                                sorted_reads = ret.1;
                                read_count = ret.0;
                            }
                        }

                        levels += 1;
                    });

                info!("writing consensus reads for reference {}", ref_name);
                // collapse the final reads down to a single sequence and write everything to the disk
                write_consensus_reads(&sorted_reads, &mut writer, levels, &rm, &40);
            }
        }
    });
}

#[allow(dead_code)]
struct JointLookup {
    bam_index: usize,
    reference_manager_index: usize,
}

#[allow(dead_code)]
struct NamedRef {
    name: String,
    sequence: String,
}

#[allow(dead_code)]
struct ReferenceLookupTable {
    bam_reference_id_to_name: HashMap<usize, String>,
    bam_reference_name_to_id: HashMap<String, usize>,
    fasta_reference_id_to_name: HashMap<usize, String>,
    fasta_reference_name_to_id: HashMap<String, usize>,
    name_to_joint_index: HashMap<String, JointLookup>,
    unified_name_to_seq: HashMap<String, String>,
}

#[allow(dead_code)]
impl ReferenceLookupTable {
    pub fn new(reference_manager: &ReferenceManager, bam_header: &Header) -> ReferenceLookupTable {
        let mut bam_reference_id_to_name = HashMap::new();
        let mut bam_reference_name_to_id = HashMap::new();
        let mut fasta_reference_id_to_name = HashMap::new();
        let mut fasta_reference_name_to_id = HashMap::new();
        let mut name_to_joint_index = HashMap::new();
        let mut unified_name_to_seq = HashMap::new();

        // TODO: the header uses an in-order map for storing reference sequences, think through this
        bam_header
            .reference_sequences()
            .iter()
            .enumerate()
            .for_each(|(index, (bstr_name, _ref_map))| {
                let string_name = bstr_name.to_string();

                // we need to have this reference sequence stored in our database as well
                if reference_manager
                    .reference_name_to_ref
                    .contains_key(string_name.as_bytes())
                {
                    let our_index = reference_manager
                        .reference_name_to_ref
                        .get(string_name.as_bytes())
                        .unwrap();
                    assert!(!bam_reference_id_to_name.contains_key(&index));
                    assert!(!fasta_reference_id_to_name.contains_key(our_index));

                    bam_reference_id_to_name.insert(index, string_name.clone());
                    bam_reference_name_to_id.insert(string_name.clone(), index);
                    fasta_reference_id_to_name.insert(*our_index, string_name.clone());
                    fasta_reference_name_to_id.insert(string_name.clone(), *our_index);

                    name_to_joint_index.insert(
                        string_name.clone(),
                        JointLookup {
                            bam_index: index,
                            reference_manager_index: *our_index,
                        },
                    );

                    unified_name_to_seq.insert(
                        string_name,
                        String::from_utf8(
                            reference_manager
                                .references
                                .get(our_index)
                                .unwrap()
                                .clone()
                                .sequence,
                        )
                        .unwrap()
                        .clone(),
                    );
                }
            });

        reference_manager
            .references
            .iter()
            .for_each(|(id, reference)| {
                if !fasta_reference_id_to_name.contains_key(id) {
                    warn!(
                        "We dont have an entry in the BAM file for reference {}",
                        String::from_utf8(reference.name.clone()).unwrap()
                    );
                }
            });

        ReferenceLookupTable {
            bam_reference_id_to_name,
            bam_reference_name_to_id,
            fasta_reference_id_to_name,
            fasta_reference_name_to_id,
            name_to_joint_index,
            unified_name_to_seq,
        }
    }
}

trait AlignmentFilter {
    fn keep(&self, read: &SortingReadSetContainer) -> bool;
}

pub struct AlignmentCheck {
    min_aligned_bases: usize,
    min_aligned_identical_proportion: f64,
}

impl AlignmentFilter for AlignmentCheck {
    fn keep(&self, read: &SortingReadSetContainer) -> bool {
        let mut alignment_count =0;
        let mut alignable_bases = 0;

        read.aligned_read.read_aligned.iter().zip(read.aligned_read.reference_aligned.iter()).for_each(|(x,y)| {
            if *x > 59 && x != &FASTA_N {
                alignable_bases += 1;
                if x == y {
                    alignment_count += 1;
                }
            }
        });

        let ret = (alignment_count as f64 / alignable_bases as f64 >= self.min_aligned_identical_proportion) && (alignable_bases >= self.min_aligned_bases);
        if !ret {
            println!("aligning {} {}\n{}\n{}",alignment_count,alignable_bases,u8s(&read.aligned_read.read_aligned), u8s(&read.aligned_read.reference_aligned));
        }
        ret

    }
}



/// We want to be extra confident in the alignments around our 'tags'.
/// This filters out reads where we have mismatches and gaps around the
/// degenerate sequences we recover
pub struct FlankingDegenerateBaseFilter {
    min_flanking_indentity: f64,
    flanking_window_size: usize,
}

impl AlignmentFilter for FlankingDegenerateBaseFilter {
    fn keep(&self, read: &SortingReadSetContainer) -> bool {
        // create a sliding window set to the flanking window size. When we hit a degenerate sequence
        // check that our window meets the criteria
        let mut pushed_binary_comp = Vec::new();
        let mut ret = true;
        let mut count_down_check = usize::MAX;
        read.aligned_read.read_aligned.iter().zip(read.aligned_read.reference_aligned.iter()).for_each(|(read_base, reference_base)| {
            if count_down_check == 0 {
                count_down_check = usize::MAX;
                let lookback_length = min(pushed_binary_comp.len(),self.flanking_window_size);
                let sum: u32 = pushed_binary_comp[pushed_binary_comp.len() - lookback_length..pushed_binary_comp.len()].iter().sum();
                let matching_prop = sum as f64 / lookback_length as f64;
                pushed_binary_comp.clear();
                if matching_prop < self.min_flanking_indentity {
                    ret = false;
                }
            }
            else if *reference_base > 58 && reference_base != &FASTA_N {
                count_down_check -= 1;
                if read_base == reference_base {
                    pushed_binary_comp.push(1)
                } else {
                    pushed_binary_comp.push(0)
                }
            }
                // lookback case for start of Ns
            else if reference_base == &FASTA_N && pushed_binary_comp.len() > 0 {
                let lookback_length = min(pushed_binary_comp.len(),self.flanking_window_size);
                let sum: u32 = pushed_binary_comp[pushed_binary_comp.len() - lookback_length..pushed_binary_comp.len()].iter().sum();
                let matching_prop = sum as f64 / lookback_length as f64;
                pushed_binary_comp.clear();
                if matching_prop < self.min_flanking_indentity {
                    ret = false;
                }
            }
            else if reference_base == &FASTA_N && pushed_binary_comp.len() == 0 {
                count_down_check = self.flanking_window_size;
            }
        });
        ret
    }
}

#[derive(Default,Copy,Clone,Debug)]
struct BamReadFiltering {
    total_reads: usize,
    unmapped_flag_reads: usize,
    secondary_flag_reads: usize,
    failed_alignment_filters: usize,
    failed_alignment_creation: usize,
    duplicate_reads: usize,
    invalid_tags: usize, 
}

impl BamReadFiltering {
    pub fn passing_reads(&self) -> usize {
        self.total_reads - self.unmapped_flag_reads - self.secondary_flag_reads - self.failed_alignment_filters - self.duplicate_reads - self.invalid_tags
    }

    pub fn results(&self, filters_counts: &HashMap<String,u64>) {
        let filter_summary = filters_counts.iter().map(|x| format!("Name: {} failed {}",x.0.clone(),x.1)).join(", ");
        info!(
            "Total reads processed: {}, Unmapped: {}, Secondary: {}, [Failed: {}, Failed alignment filters: {}, Duplicate: {}, Invalid_tags: {}, Passing: {} filter summary {}",
            self.total_reads,
            self.unmapped_flag_reads,
            self.secondary_flag_reads,
            self.failed_alignment_creation,
            self.failed_alignment_filters,
            self.duplicate_reads,
            self.invalid_tags,
            self.passing_reads(),
            filter_summary,
        );
    }
}

pub struct SortedReadsFromBam {
    pub bam: Option<ShardReader<SortingReadSetContainer>>,
    pub read_stats: BamReadFiltering,
}

pub fn sort_reads_from_bam_file(
    bam_file: &String,
    reference_name: &String,
    reference_manager: &ReferenceManager,
    read_structure: &SequenceLayout,
    temp_directory: &mut InstanceLivedTempDir,
) -> SortedReadsFromBam {

    let aligned_temp = temp_directory.temp_file("bam.reads.sorted.sharded");

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let mut bai_file = bam_file.clone();
    bai_file.push_str(".bai");

    let mut read_stats = BamReadFiltering::default();

    let filters : Vec<(String, &dyn AlignmentFilter)> = vec![("FlankingDegenerateBaseFilter".to_string(),&FlankingDegenerateBaseFilter{ min_flanking_indentity: 0.80, flanking_window_size: 10}),
                                                             ("AlignmentCheck".to_string(),&AlignmentCheck{ min_aligned_bases: 45, min_aligned_identical_proportion: 0.8})];
    let mut filter_counts: HashMap<String,u64> = HashMap::default();
    filter_counts.insert("FlankingDegenerateBaseFilter".to_string(),0);
    filter_counts.insert("AlignmentCheck".to_string(),0);

    let index = bai::read(bai_file).expect("Unable to open BAM BAI file");
    let header = reader.read_header().unwrap();
    {
        let mut sharded_output: ShardWriter<SortingReadSetContainer> =
            ShardWriter::new(&aligned_temp, 32, 256, 1 << 16).unwrap();
        let mut sender = sharded_output.get_sender();

        let reference_sequence_id = reference_manager
            .reference_name_to_ref
            .get(reference_name.as_bytes())
            .unwrap();

        let reference_sequence = reference_manager
            .references
            .get(reference_sequence_id)
            .unwrap()
            .sequence
            .clone();

        let reference_config = read_structure.references.get(reference_name).unwrap();

        let region = &reference_name.parse().expect("Unable to parse chromosome");

        let records = reader
            .query(&header, &index, &region)
            .map(Box::new)
            .expect("Unable to parse out region information");

        for result in records {
            read_stats.total_reads += 1;
            if read_stats.total_reads % 1000000 == 0 {
                read_stats.results(&filter_counts);
            }

            let record = result.unwrap();

            if !record.flags().is_secondary() && !record.flags().is_unmapped() {
                let read = create_sorted_read_container(reference_name, &reference_manager, &mut read_stats, &reference_sequence_id, &reference_sequence, reference_config, &record);


                match read {
                    Some(x) => {
                        let survives_filtering = filters.iter().map(|t| {
                            let x = t.1.keep(&x);
                            if !x {
                                filter_counts.insert(t.0.clone(),filter_counts.get(&t.0).unwrap_or(&0) + 1);
                            }
                            x
                        }).filter(|b| !*b).count() == 0;
                        if survives_filtering {
                            sender.send(x).unwrap();
                        } else {
                            read_stats.failed_alignment_filters += 1;
                        }
                    },
                    None => {
                        read_stats.failed_alignment_creation += 1;
                    }
                }
            } else {
                if record.flags().is_secondary() {
                    read_stats.secondary_flag_reads += 1;
                }
                if record.flags().is_unmapped() {
                    read_stats.unmapped_flag_reads += 1;
                }
            }
        }
        sender.finished().unwrap();
        sharded_output.finish().unwrap();

    }
    read_stats.results(&filter_counts);
    if read_stats.passing_reads() > 0 {
        SortedReadsFromBam{
            bam: Some(ShardReader::open(aligned_temp).unwrap()),
            read_stats,
        }
    } else {
        SortedReadsFromBam{
            bam: None,
            read_stats,
        }
    }
}

/// create a read container from the read and reference sequence,
/// returning Some(read) if successful, or None if we couldn't
/// extract the tag
fn create_sorted_read_container(reference_name: &String,
                                reference_manager: &&ReferenceManager,
                                _read_stats: &mut BamReadFiltering,
                                reference_sequence_id: &&usize,
                                reference_sequence: &Vec<u8>,
                                reference_config: &ReferenceRecord, record: &Record) -> Option<SortingReadSetContainer> {

    let seq: Vec<u8> = record.sequence().iter().collect();
    let start_pos = record.alignment_start().unwrap().unwrap().get();
    let cigar = record.cigar();
    let read_name: bam::record::Name<'_> = record.name().unwrap();
    let read_qual = record.quality_scores().iter().collect();
    let ref_slice = reference_sequence.as_slice();

    let aligned_read =
        recover_soft_clipped_align_sequences(&seq, start_pos, &cigar.iter().map(|x| x.unwrap()).collect(), &SoftClipResolution::Realign, ref_slice);

    let stretched_alignment = stretch_sequence_to_alignment(
        &aligned_read.aligned_ref,
        &reference_manager
            .references
            .get(reference_sequence_id)
            .unwrap()
            .sequence,
    );

    let extracted_tags =
        extract_tagged_sequences(&aligned_read.aligned_read, &stretched_alignment);

    let (valid_tags_extracted, read_tags_ordered) =
        extract_tag_sequences(reference_config, extracted_tags);

    if !valid_tags_extracted {
        Some(SortingReadSetContainer {
            ordered_sorting_keys: vec![], // for future use during sorting
            ordered_unsorted_keys: read_tags_ordered, // the current unsorted tag collection
            aligned_read: AlignmentResult {
                reference_name: reference_name.clone(),
                read_aligned: aligned_read.aligned_read,
                read_quals: Some(read_qual),
                cigar_string: cigar
                    .iter()
                    .map(|op| AlignmentTag::from(op.unwrap()))
                    .collect(),
                path: vec![],
                score: 0.0,
                reference_start: start_pos,
                read_start: 0,
                reference_aligned: aligned_read.aligned_ref,
                read_name: String::from_utf8(read_name.as_bytes().to_vec()).unwrap(),
                bounding_box: None,
            },
        })

    } else {
        None
    }
}

fn get_known_level_lookups(read_structure: &SequenceLayout) -> HashMap<String, KnownList> {
    let mut ret: HashMap<String, KnownList> = HashMap::new();

    read_structure
        .references
        .iter()
        .for_each(|(_name, reference)| {
            reference
                .umi_configurations
                .iter()
                .for_each(|(_name, config)| match &config.file {
                    None => {}
                    Some(x) => {
                        if !ret.contains_key(x.as_str()) {
                            let known_lookup = KnownList::new(config, &8);
                            ret.insert(x.clone(), known_lookup);
                        }
                    }
                })
        });
    ret
}

/// Sorts the reads by the degenerate tag
///
/// we group reads into a container where previous tags all match. We then determine the clique of
/// degenerate tags within the container and correct the sequences to the consensuses of cliques within
/// the group. For example, if we've sorted by a 10X cell ID tag, that we would collect cells that have the
/// same cell ID, and then correct the UMI sequences within the cell to the consensuses of the UMI sequences
///
/// # Arguments
///     * `temp_directory` - The temporary directory to use for sorting
///    * `reader` - The reader to sort
///   * `tag` - The tag to sort by
///
/// # Returns
///    * `ShardReader` - The sorted reader
///
/// # Errors
///     * `std::io::Error` - If there is an error writing to the temporary directory
///
/// # Panics
///    * If the tag is not a known tag
///
/// # Examples
///
pub fn sort_degenerate_level(
    temp_directory: &mut InstanceLivedTempDir,
    reader: &ShardReader<SortingReadSetContainer>,
    tag: &UMIConfiguration,
    iteration: &usize,
    read_count: &usize,
) -> (usize, ShardReader<SortingReadSetContainer>) {

    info!("Sorting degenerate level {}", tag.symbol);

    let mut all_read_count: usize = 0;
    let mut output_reads: usize = 0;

    let aligned_temp = temp_directory.temp_file(&*(tag.order.to_string() + ".sorted.sharded"));
    let mut sharded_output: ShardWriter<SortingReadSetContainer> =
        ShardWriter::new(&aligned_temp, 32, 256, 1 << 16).unwrap();

    let mut sender = sharded_output.get_sender();
    let mut bar: Option<ProgressBar> = match *read_count > 100000 {
        true => Some(ProgressBar::new(read_count.clone() as u64)),
        false => None,
    };

    let mut last_read: Option<SortingReadSetContainer> = None;

    let maximum_reads_per_bin = if tag.maximum_subsequences.is_some() {
        tag.maximum_subsequences.clone().unwrap()
    } else {
        10000 // TODO make this a constant somewhere
    };
    info!("Starting to sort degenerate level {}", tag.symbol);

    let mut current_sorting_bin: Option<DegenerateBuffer> = None;

    reader.iter_range(&Range::all()).unwrap().for_each(|current_read| {
        all_read_count += 1;
        //if all_read_count % 10000 == 0 { // TODO for now to figure out what's going on
            bar.as_mut().map(|b| b.set_position(all_read_count as u64));
        //}
        let mut current_read = current_read.unwrap();
        let next_last_read = current_read.clone();

        match current_sorting_bin.as_mut() {
            None => {
                let mut bin = DegenerateBuffer::new(
                    temp_directory.temp_file(format!("{}.fasta", tag.order).as_str()),
                    &maximum_reads_per_bin,
                    tag.clone(),
                );
                bin.push(current_read);
                current_sorting_bin = Some(bin);
            }

            Some(bin) => {
                let reads_equal = last_read.as_ref().unwrap().cmp(&mut current_read) == Ordering::Equal;

                match reads_equal {
                    true => {
                        // add the current read to the bin
                        bin.push(current_read);
                    }
                    false => {
                        // write the previous bin, and add the current read to the next bin
                        output_reads += bin.close_and_write_to_shard_writer(&mut sender);
                        bin.push(current_read);
                    }
                }
            }
        };

        last_read = Some(next_last_read);
    });

    match current_sorting_bin {
        None => {}
        Some(mut bin) => {
            output_reads += bin.close_and_write_to_shard_writer(&mut sender);

        }
    }

    bar.as_mut().map(|b| b.set_position(all_read_count as u64));

    info!("For degenerate tag {} (iteration {}) we processed {} reads, of which {} were passed to the next level", &tag.symbol, iteration, all_read_count, output_reads);
    sender.finished().unwrap();

    sharded_output.finish().unwrap();

    (output_reads, ShardReader::open(aligned_temp).unwrap())
}

pub fn sort_known_level(
    temp_directory: &mut InstanceLivedTempDir,
    reader: &ShardReader<SortingReadSetContainer>,
    tag: &UMIConfiguration,
    read_count: &usize,
    known_lookup_obj: &mut HashMap<String, KnownList>,
) -> (usize, ShardReader<SortingReadSetContainer>) {
    info!("Sorting known level {}", tag.symbol);

    info!(
        "Loading the known lookup table for tag {}, this can take some time",
        tag.symbol
    );
    let known_lookup = known_lookup_obj
        .get_mut(&tag.file.as_ref().unwrap().clone())
        .expect(
            format!(
                "Unable to find pre-cached lookup table {}",
                &tag.file.as_ref().clone().unwrap()
            )
            .as_str(),
        );
    let mut processed_reads = 0;
    let mut dropped_reads = 0;
    let mut collided_reads = 0;
    let mut sent_reads = 0;

    info!("Sorting {} reads", read_count);
    let mut bar: Option<ProgressBar> = match *read_count > 100000 {
        true => Some(ProgressBar::new(read_count.clone() as u64)),
        false => Some(ProgressBar::new(read_count.clone() as u64)),
    };

    // create a new output
    let aligned_temp = temp_directory.temp_file(&*(tag.order.to_string() + ".sorted.sharded"));

    // scoping required here to ensure proper handing of the senders
    {
        let mut sharded_output: ShardWriter<SortingReadSetContainer> =
            ShardWriter::new(&aligned_temp, 32, 256, 1 << 16).unwrap();
        let mut sender = sharded_output.get_sender();

        reader.iter_range(&Range::all()).unwrap().for_each(|x| {
            processed_reads += 1;
            //if processed_reads % 10000 == 0 {
                bar.as_mut().map(|b| b.set_position(processed_reads as u64));
            //}
            let mut sorting_read_set_container = x.unwrap();
            assert_eq!(
                sorting_read_set_container.ordered_sorting_keys.len(),
                tag.order
            );

            let next_key = sorting_read_set_container
                .ordered_unsorted_keys
                .pop_front()
                .unwrap();
            assert_eq!(next_key.0, tag.symbol);

            let corrected_hits =
            known_lookup.correct_to_known_list(&next_key.1,  &(tag.max_distance as u32));
            match (corrected_hits.hits.len(), corrected_hits.distance) {
                (x, _) if x < 1 => {
                    dropped_reads += 1;
                }
                (x, _) if x > 1 => {
                    collided_reads += 1;
                    dropped_reads += 1;
                }
                (_x, y) if y > (tag.max_distance as u32)=> {
                    dropped_reads += 1;
                }
                (_x, _y) => {
                    sent_reads += 1;
                    sorting_read_set_container
                        .ordered_sorting_keys
                        .push((next_key.0, corrected_hits.hits.get(0).unwrap().clone()));
                    sender.send(sorting_read_set_container).unwrap();
                }
            };
        });

        info!(
            "Dropped {} reads (of which {} were collided reads), {} total reads, sent reads: {}",
            dropped_reads, collided_reads, processed_reads, sent_reads
        );
        sender.finished().unwrap();
        sharded_output.finish().unwrap();
    }
    bar.as_mut().map(|b| b.set_position(processed_reads as u64));

    info!(
        "For known tag {} we processed {} reads",
        &tag.symbol, processed_reads
    );

    (
        processed_reads - dropped_reads,
        ShardReader::open(aligned_temp).unwrap(),
    )
}

#[cfg(test)]
mod tests {
    use std::collections::VecDeque;
    use ::{FASTA_A, FASTA_T};
    use super::*;
    use crate::alignment::alignment_matrix::AlignmentResult;
    use crate::utils::read_utils::fake_reads;

    pub fn consensus(input: &Vec<Vec<u8>>) -> Vec<u8> {
        let mut consensus = Vec::new();

        // for each position
        for i in 0..input[0].len() {
            let mut counter = HashMap::new();

            // for each input string
            input.iter().for_each(|vector| {
                assert_eq!(
                    vector.len(),
                    input[0].len(),
                    "string {} is not the same length as the first string {}",
                    String::from_utf8(vector.clone()).unwrap(),
                    String::from_utf8(input[0].clone()).unwrap()
                );

                *counter.entry(&vector[i]).or_insert(0) += 1;
            });

            let mut max = 0;
            let mut consensus_byte = b'N';

            //println!("consensus {:?}",counter);
            for (byte, count) in counter {
                // if we're the new maximum OR we're tied for the maximum and we're an N or a gap, then we'll take the new value
                if count > max
                    || (count == max && consensus_byte == b'N')
                    || (count == max && consensus_byte == b'-')
                {
                    max = count;
                    consensus_byte = *byte;
                }
            }

            consensus.push(consensus_byte);
        }

        consensus
    }

    #[test]
    fn test_alignment_check() {
        let alignment_check = AlignmentCheck{ min_aligned_bases: 10, min_aligned_identical_proportion: 0.8 };

        let fake_read_alignment = SortingReadSetContainer{
            ordered_sorting_keys: vec![],
            ordered_unsorted_keys: Default::default(),
            aligned_read: AlignmentResult {
                reference_name: "".to_string(),
                read_name: "".to_string(),
                reference_aligned: vec![FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A],
                read_aligned: vec![FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,],
                read_quals: None,
                cigar_string: vec![],
                path: vec![],
                score: 0.0,
                reference_start: 0,
                read_start: 0,
                bounding_box: None,
            },
        };

        assert!(alignment_check.keep(&fake_read_alignment));
    }

    // FlankingDegenerateBaseFilter

    #[test]
    fn test_flanking_degenerate_base_filter() {
        let alignment_check = FlankingDegenerateBaseFilter{ min_flanking_indentity: 0.9, flanking_window_size: 3 };

        let fake_read_alignment = SortingReadSetContainer{
            ordered_sorting_keys: vec![],
            ordered_unsorted_keys: Default::default(),
            aligned_read: AlignmentResult {
                reference_name: "".to_string(),
                read_name: "".to_string(),
                reference_aligned: vec![FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_N,FASTA_N,FASTA_N,FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A],
                read_aligned: vec![FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,],
                read_quals: None,
                cigar_string: vec![],
                path: vec![],
                score: 0.0,
                reference_start: 0,
                read_start: 0,
                bounding_box: None,
            },
        };

        assert!(alignment_check.keep(&fake_read_alignment));

        let fake_read_alignment = SortingReadSetContainer{
            ordered_sorting_keys: vec![],
            ordered_unsorted_keys: Default::default(),
            aligned_read: AlignmentResult {
                reference_name: "".to_string(),
                read_name: "".to_string(),
                reference_aligned: vec![FASTA_A,FASTA_A,FASTA_A,FASTA_A,  FASTA_N,FASTA_N,FASTA_N,  FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A],
                read_aligned:      vec![FASTA_A,FASTA_A,FASTA_A, FASTA_T,  FASTA_A,FASTA_A,FASTA_A,  FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,],
                read_quals: None,
                cigar_string: vec![],
                path: vec![],
                score: 0.0,
                reference_start: 0,
                read_start: 0,
                bounding_box: None,
            },
        };

        assert!(!alignment_check.keep(&fake_read_alignment));
        let fake_read_alignment = SortingReadSetContainer{
            ordered_sorting_keys: vec![],
            ordered_unsorted_keys: Default::default(),
            aligned_read: AlignmentResult {
                reference_name: "".to_string(),
                read_name: "".to_string(),
                reference_aligned: vec![FASTA_A,FASTA_A,FASTA_A,FASTA_A,  FASTA_N,FASTA_N,FASTA_N,  FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A],
                read_aligned:      vec![FASTA_A,FASTA_A,FASTA_A,FASTA_A,  FASTA_A,FASTA_A,FASTA_A,  FASTA_T,FASTA_A,FASTA_A,FASTA_A,FASTA_A,],
                read_quals: None,
                cigar_string: vec![],
                path: vec![],
                score: 0.0,
                reference_start: 0,
                read_start: 0,
                bounding_box: None,
            },
        };

        assert!(!alignment_check.keep(&fake_read_alignment));

        assert!(!alignment_check.keep(&fake_read_alignment));
        let fake_read_alignment = SortingReadSetContainer{
            ordered_sorting_keys: vec![],
            ordered_unsorted_keys: Default::default(),
            aligned_read: AlignmentResult {
                reference_name: "".to_string(),
                read_name: "".to_string(),
                reference_aligned: vec![FASTA_N,FASTA_N,FASTA_N,  FASTA_A,FASTA_A,FASTA_A,FASTA_A,  FASTA_N,FASTA_N,FASTA_N,  FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A],
                read_aligned:      vec![FASTA_A,FASTA_A,FASTA_A,  FASTA_A,FASTA_A,FASTA_A,FASTA_A,  FASTA_A,FASTA_A,FASTA_A,  FASTA_T,FASTA_A,FASTA_A,FASTA_A,FASTA_A,],
                read_quals: None,
                cigar_string: vec![],
                path: vec![],
                score: 0.0,
                reference_start: 0,
                read_start: 0,
                bounding_box: None,
            },
        };

        assert!(!alignment_check.keep(&fake_read_alignment));


        assert!(!alignment_check.keep(&fake_read_alignment));
        let fake_read_alignment = SortingReadSetContainer{
            ordered_sorting_keys: vec![],
            ordered_unsorted_keys: Default::default(),
            aligned_read: AlignmentResult {
                reference_name: "".to_string(),
                read_name: "".to_string(),
                reference_aligned: vec![FASTA_N,FASTA_N,FASTA_N,  FASTA_A,FASTA_A,FASTA_A,FASTA_A,  FASTA_N,FASTA_N,FASTA_N,  FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,  FASTA_N,FASTA_N,FASTA_N],
                read_aligned:      vec![FASTA_A,FASTA_A,FASTA_A,  FASTA_A,FASTA_A,FASTA_A,FASTA_A,  FASTA_A,FASTA_A,FASTA_A,  FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_T,  FASTA_A,FASTA_A,FASTA_A],
                read_quals: None,
                cigar_string: vec![],
                path: vec![],
                score: 0.0,
                reference_start: 0,
                read_start: 0,
                bounding_box: None,
            },
        };

        assert!(!alignment_check.keep(&fake_read_alignment));


        assert!(!alignment_check.keep(&fake_read_alignment));
        let fake_read_alignment = SortingReadSetContainer{
            ordered_sorting_keys: vec![],
            ordered_unsorted_keys: Default::default(),
            aligned_read: AlignmentResult {
                reference_name: "".to_string(),
                read_name: "".to_string(),
                reference_aligned: vec![FASTA_N,FASTA_N,FASTA_N,  FASTA_A,FASTA_A,FASTA_A,FASTA_A,  FASTA_N,FASTA_N,FASTA_N,  FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,  FASTA_N,FASTA_N,FASTA_N],
                read_aligned:      vec![FASTA_A,FASTA_A,FASTA_A,  FASTA_A,FASTA_A,FASTA_A,FASTA_A,  FASTA_A,FASTA_A,FASTA_A,  FASTA_A,FASTA_A,FASTA_A,FASTA_A,FASTA_A,  FASTA_A,FASTA_A,FASTA_A],
                read_quals: None,
                cigar_string: vec![],
                path: vec![],
                score: 0.0,
                reference_start: 0,
                read_start: 0,
                bounding_box: None,
            },
        };

        assert!(alignment_check.keep(&fake_read_alignment));


    }
    #[test]
    fn test_consensus() {
        let basic_seqs: Vec<Vec<u8>> = vec![
            String::from("ATCG").as_bytes().to_vec(),
            String::from("GCTA").as_bytes().to_vec(),
            String::from("ATCG").as_bytes().to_vec(),
        ];
        let cons = consensus(&basic_seqs);
        assert_eq!(cons, String::from("ATCG").as_bytes().to_vec());

        let basic_seqs: Vec<Vec<u8>> = vec![
            String::from("ATCG").as_bytes().to_vec(),
            String::from("ATC-").as_bytes().to_vec(),
        ];
        let cons = consensus(&basic_seqs);
        assert_eq!(cons, String::from("ATCG").as_bytes().to_vec());

        // reverse the above to check that order doesn't matter (it did at one point)
        let basic_seqs: Vec<Vec<u8>> = vec![
            String::from("ATC-").as_bytes().to_vec(),
            String::from("ATCG").as_bytes().to_vec(),
        ];
        let cons = consensus(&basic_seqs);
        assert_eq!(cons, String::from("ATCG").as_bytes().to_vec());

        // real world issue
        let basic_seqs: Vec<Vec<u8>> = vec![
            String::from("TGGTATGCTGG-").as_bytes().to_vec(),
            String::from("TGGTATGCTGGG").as_bytes().to_vec(),
        ];
        let cons = consensus(&basic_seqs);
        assert_eq!(cons, String::from("TGGTATGCTGGG").as_bytes().to_vec());

        // reverse the above to check that order doesn't matter (it did at one point)
        let basic_seqs: Vec<Vec<u8>> = vec![
            String::from("TGGTATGCTGG-").as_bytes().to_vec(),
            String::from("TGGTATGCTGGG").as_bytes().to_vec(),
        ];
        let cons = consensus(&basic_seqs);
        assert_eq!(cons, String::from("TGGTATGCTGGG").as_bytes().to_vec());
    }

    //TODO fix this #[test]
    fn test_consensus_real_world() {
        let reads = fake_reads(10, 1);

        let read_seq = reads
            .get(0)
            .unwrap()
            .read_one
            .seq()
            .iter()
            .map(|x| x.clone())
            .collect::<Vec<u8>>();

        let fake_read = AlignmentResult {
            reference_name: "".to_string(),
            read_name: "".to_string(),
            reference_aligned: read_seq.clone(),
            read_aligned: read_seq.clone(),
            read_quals: None,
            cigar_string: vec![],
            path: vec![],
            score: 0.0,
            reference_start: 0,
            read_start: 0,
            bounding_box: None,
        };

        let mut tbb = DegenerateBuffer::new(
            PathBuf::from("test_data/consensus_test.fastq"),
            &1000,
            UMIConfiguration {
                symbol: 'c',
                file: None,
                reverse_complement_sequences: None,
                sort_type: UMISortType::DegenerateTag,
                length: "TGGTATGCTGGG".len(),
                order: 0,
                pad: None,
                max_distance: 2, // TODO: the plus 1 here is confusing -- does STARCODE
                maximum_subsequences: None,
                max_gaps: Some(1),
                minimum_collapsing_difference: Some(10.0),
            },
        );

        // real example we hit
        let st1 = SortingReadSetContainer {
            ordered_sorting_keys: vec![
                ('a', "AAACCCATCAGCATTA".as_bytes().to_vec()),
                ('b', "TATTGACAACCT".as_bytes().to_vec()),
            ],
            ordered_unsorted_keys: VecDeque::from(vec![('c', "TGGTATGCTGG".as_bytes().to_vec())]),
            aligned_read: fake_read.clone(),
        };

        let st2 = SortingReadSetContainer {
            ordered_sorting_keys: vec![
                ('a', "AAACCCATCAGCATTA".as_bytes().to_vec()),
                ('b', "TATTGACAACCT".as_bytes().to_vec()),
            ],
            ordered_unsorted_keys: VecDeque::from(vec![('c', "TGGTATGCTGGG".as_bytes().to_vec())]),
            aligned_read: fake_read.clone(),
        };

        assert_eq!(st1.cmp(&st2) == Ordering::Equal, true);

        tbb.push(st1);
        for _ in 0..10 {
            tbb.push(st2.clone());
        }
        let correction = tbb.correct_list();

        correction.iter().for_each(|x| {
            println!("{} -> {}",String::from_utf8(x.0.clone()).unwrap(),String::from_utf8(x.1.clone()).unwrap());
            assert_eq!(
                String::from_utf8(x.1.clone()).unwrap(),
                String::from_utf8("TGGTATGCTGGG".as_bytes().to_vec()).unwrap()
            );
        });
    }
}
