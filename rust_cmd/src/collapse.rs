//! The `collapse` pipeline: for each reference, read its aligned records from a
//! BAM, filter and extract UMI/barcode tags, sort reads down through the UMI
//! hierarchy (correcting tags at each level), then emit one consensus (or
//! corrected) read per molecule.

use crate::consensus::consensus_builders::{write_consensus_reads, ConsensusWriteStats};
use crate::extractor::{
    extract_tag_sequences, extract_tagged_sequences, recover_soft_clipped_align_sequences,
    stretch_sequence_to_alignment, SoftClipResolution,
};
use crate::read_strategies::read_disk_sorter::SortingReadSetContainer;
use crate::read_strategies::sequence_layout::{ReferenceRecord, SequenceLayout, UMIConfiguration};
use crate::reference::fasta_reference::ReferenceManager;
use crate::run_summary::{RunSummary, SummaryTable};

use crate::InstanceLivedTempDir;
use indicatif::ProgressBar;

use noodles_bam as bam;
use noodles_bam::{bai, Record};
use noodles_sam::{alignment::record::Flags, Header};

use itertools::Itertools;
use shardio::{Range, ShardReader, ShardWriter};
use std::cmp::{min, Ordering};
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::time::Instant;

use crate::alignment::alignment_matrix::{AlignmentResult, AlignmentTag};
use crate::alignment_manager::{BamFileAlignmentWriter, OutputAlignmentWriter};
use crate::umis::correct_tags::SequenceCorrector;
use consensus::consensus_builders::{write_corrected_reads, MergeStrategy, ReadOutputApproach};
use read_strategies::sequence_layout::UMISortType;
use rust_star::Trie;
use umis::known_list::KnownList;
use utils::read_utils::{reverse_complement, u8s};
use FASTA_N;

const CONSENSUS_HISTOGRAM_WIDTH: usize = 40;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UmiLevelRunStats {
    reference: String,
    symbol: char,
    sort_type: String,
    input_reads: usize,
    output_reads: usize,
}

#[derive(Clone, Debug)]
pub struct CollapseReferenceStats {
    reference: String,
    filtering: BamReadFiltering,
    reads_after_umi_correction: usize,
    umi_levels: Vec<UmiLevelRunStats>,
    output: ConsensusWriteStats,
}

#[derive(Clone, Debug)]
pub struct CollapseRunStats {
    mode: String,
    references: Vec<CollapseReferenceStats>,
    elapsed_seconds: f64,
}

impl CollapseRunStats {
    pub fn summary(&self) -> RunSummary {
        let mut total_filtering = BamReadFiltering::default();
        let mut total_after_umi = 0usize;
        let mut total_output = ConsensusWriteStats::default();
        for reference in &self.references {
            total_filtering.accumulate(&reference.filtering);
            total_after_umi += reference.reads_after_umi_correction;
            total_output.accumulate(&reference.output);
        }

        let mut results = SummaryTable::new(
            "Collapse results",
            &[
                "Reference",
                "BAM records",
                "Passed filters",
                "After UMI",
                "UMI removed",
                "Groups",
                "Reads selected",
                "Downsampled",
                "Failed groups",
                "Output records",
            ],
        );
        push_collapse_result_row(
            &mut results,
            "All",
            &total_filtering,
            total_after_umi,
            &total_output,
        );
        let mut references = self.references.iter().collect::<Vec<_>>();
        references.sort_by(|left, right| left.reference.cmp(&right.reference));
        for reference in references {
            push_collapse_result_row(
                &mut results,
                &reference.reference,
                &reference.filtering,
                reference.reads_after_umi_correction,
                &reference.output,
            );
        }

        let mut filtering = SummaryTable::new(
            "Collapse filtering",
            &[
                "Reference",
                "Unmapped",
                "Secondary",
                "Supplementary",
                "Alignment filter",
                "Reconstruction/tag",
                "Duplicate",
                "Invalid tags",
            ],
        );
        push_filtering_row(&mut filtering, "All", &total_filtering);
        let mut references = self.references.iter().collect::<Vec<_>>();
        references.sort_by(|left, right| left.reference.cmp(&right.reference));
        for reference in references {
            push_filtering_row(&mut filtering, &reference.reference, &reference.filtering);
        }

        let mut umi = SummaryTable::new(
            "UMI correction levels",
            &["Reference / UMI", "Type", "Input reads", "Output reads", "Removed"],
        );
        let mut levels = self
            .references
            .iter()
            .flat_map(|reference| reference.umi_levels.iter())
            .collect::<Vec<_>>();
        levels.sort_by(|left, right| {
            left.reference
                .cmp(&right.reference)
                .then(left.symbol.cmp(&right.symbol))
        });
        for level in levels {
            umi.push_row([
                format!("{} / {}", level.reference, level.symbol),
                level.sort_type.clone(),
                level.input_reads.to_string(),
                level.output_reads.to_string(),
                level.input_reads.saturating_sub(level.output_reads).to_string(),
            ]);
        }

        let mut run = SummaryTable::new("Run details", &["Scope", "Mode", "Elapsed (s)"]);
        run.push_row([
            "All".to_string(),
            self.mode.clone(),
            format!("{:.2}", self.elapsed_seconds),
        ]);

        let mut summary = RunSummary::new("Collapse run summary");
        summary.add_table(results);
        summary.add_table(filtering);
        if !self.references.iter().all(|reference| reference.umi_levels.is_empty()) {
            summary.add_table(umi);
        }
        if self.mode == "collapse" && !total_output.reads_per_consensus.is_empty() {
            summary.add_table(consensus_size_histogram(&total_output.reads_per_consensus));
        }
        summary.add_table(run);
        summary
    }
}

fn consensus_size_histogram(read_counts: &BTreeMap<usize, usize>) -> SummaryTable {
    let mut binned_counts = BTreeMap::<(usize, usize), usize>::new();
    for (reads_per_consensus, consensus_count) in read_counts {
        let bounds = consensus_size_bin(*reads_per_consensus);
        *binned_counts.entry(bounds).or_insert(0) += consensus_count;
    }

    let total_consensuses = binned_counts.values().sum::<usize>();
    let maximum_bin_count = binned_counts.values().copied().max().unwrap_or(0);
    let mut histogram = SummaryTable::new(
        "Reads per consensus (rc tag)",
        &[
            "Reads per consensus",
            "Consensus records",
            "Consensus records (%)",
            "Histogram",
        ],
    );

    for ((lower, upper), consensus_count) in binned_counts {
        let label = if lower == upper {
            lower.to_string()
        } else {
            format!("{}-{}", lower, upper)
        };
        let percentage = if total_consensuses == 0 {
            0.0
        } else {
            100.0 * consensus_count as f64 / total_consensuses as f64
        };
        histogram.push_row([
            label,
            consensus_count.to_string(),
            format!("{:.2}", percentage),
            consensus_histogram_bar(consensus_count, maximum_bin_count),
        ]);
    }

    histogram
}

fn consensus_size_bin(read_count: usize) -> (usize, usize) {
    match read_count {
        0..=2 => (read_count, read_count),
        _ => {
            let upper = read_count.checked_next_power_of_two().unwrap_or(usize::MAX);
            (upper / 2 + 1, upper)
        }
    }
}

fn consensus_histogram_bar(count: usize, maximum: usize) -> String {
    if count == 0 || maximum == 0 {
        return String::new();
    }
    let bar_length = count
        .saturating_mul(CONSENSUS_HISTOGRAM_WIDTH)
        .saturating_add(maximum - 1)
        / maximum;
    "#".repeat(bar_length.max(1).min(CONSENSUS_HISTOGRAM_WIDTH))
}

fn push_collapse_result_row(
    table: &mut SummaryTable,
    reference: &str,
    filtering: &BamReadFiltering,
    reads_after_umi: usize,
    output: &ConsensusWriteStats,
) {
    table.push_row([
        reference.to_string(),
        filtering.total_reads.to_string(),
        filtering.passing_reads().to_string(),
        reads_after_umi.to_string(),
        filtering
            .passing_reads()
            .saturating_sub(reads_after_umi)
            .to_string(),
        output.groups_attempted.to_string(),
        output.reads_selected.to_string(),
        output.reads_downsampled.to_string(),
        output.failed_groups.to_string(),
        output.output_reads.to_string(),
    ]);
}

fn push_filtering_row(table: &mut SummaryTable, reference: &str, filtering: &BamReadFiltering) {
    table.push_row([
        reference.to_string(),
        filtering.unmapped_flag_reads.to_string(),
        filtering.secondary_flag_reads.to_string(),
        filtering.supplementary_flag_reads.to_string(),
        filtering.failed_alignment_filters.to_string(),
        filtering.failed_alignment_creation.to_string(),
        filtering.duplicate_reads.to_string(),
        filtering.invalid_tags.to_string(),
    ]);
}

/// Collapses aligned reads from a BAM file by processing UMI configurations and generating consensus sequences.
///
/// This function processes BAM file reads for each reference sequence, sorts them by UMI configurations,
/// and generates consensus reads that are written to an output BAM file. The collapse process involves
/// multiple levels of sorting based on the UMI structure defined in the sequence layout.
///
/// # Arguments
/// * `final_output` - Path to the output BAM file where collapsed reads will be written
/// * `temp_directory` - Temporary directory for intermediate file storage during processing
/// * `read_structure` - Configuration defining the sequence layout and UMI structure
/// * `bam_file` - Path to the input BAM file containing aligned reads
/// * `merge_strategy` - Strategy for merging reads during consensus generation
///
/// # Process
/// 1. Loads reference sequences from the sequence layout
/// 2. Creates lookup tables for known UMI sequences
/// 3. For each reference sequence:
///    - Sorts reads from BAM file by reference
///    - Applies multiple levels of UMI-based sorting
///    - Generates consensus reads using the specified merge strategy
/// 4. Writes collapsed consensus reads to the output BAM file
///
/// # Examples
/// ```
/// use clique::collapse::collapse;
/// use clique::consensus::consensus_builders::MergeStrategy;
/// 
/// collapse(
///     &"output.bam".to_string(),
///     &mut temp_dir,
///     &sequence_layout,
///     &"input.bam".to_string(),
///     &MergeStrategy::Stretcher,
///     &ReadOutputApproach::Collapse,
///     &AlignmentFilterConfig::default(),
///     &1,
///     &40,
/// );
/// ```
pub fn collapse(
    final_output: &String,
    temp_directory: &mut InstanceLivedTempDir,
    read_structure: &SequenceLayout,
    bam_file: &String,
    merge_strategy: &MergeStrategy,
    output_approach: &ReadOutputApproach,
    alignment_filter: &AlignmentFilterConfig,
    processing_threads: &usize,
    maximum_reads_before_downsampling: &usize,
) -> CollapseRunStats {
    assert!(
        *processing_threads > 0,
        "Collapse threads must be greater than zero"
    );
    let start = Instant::now();

    // load up the reference files
    let rm = ReferenceManager::from_yaml_input(read_structure, 8, 4);

    let mut known_level_lookups = get_known_level_lookups(read_structure);

    let mut writer = BamFileAlignmentWriter::new(&PathBuf::from(final_output), &rm);
    let mut reference_stats = Vec::new();

    // for each reference, we fetch aligned reads, pull the sorting tags, and output the collapsed reads to a BAM file
    for (_id, reference) in rm.references.iter() {
        let ref_name = String::from_utf8(reference.name.clone()).unwrap();
        info!("processing reads from input BAM file: {}", bam_file);

        let sorted_reads_option =
            sort_reads_from_bam_file(
                bam_file,
                &ref_name,
                &rm,
                read_structure,
                temp_directory,
                alignment_filter,
            );
        let filtering = sorted_reads_option.read_stats;
        let mut read_count = filtering.passing_reads();
        let mut umi_levels = Vec::new();
        let mut output_stats = ConsensusWriteStats::default();

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
                        let input_reads = read_count;
                        let ret = sort_level(
                            temp_directory,
                            &sorted_reads,
                            &tag,
                            &levels,
                            &read_count,
                            &mut known_level_lookups,
                        );
                        sorted_reads = ret.1;
                        read_count = ret.0;
                        umi_levels.push(UmiLevelRunStats {
                            reference: ref_name.clone(),
                            symbol: tag.symbol,
                            sort_type: format!("{:?}", tag.sort_type),
                            input_reads,
                            output_reads: read_count,
                        });

                        levels += 1;
                    });

                // collapse the final reads down to a single sequence and write everything to the disk
                
                match output_approach {
                    ReadOutputApproach::Collapse => {
                        info!("writing consensus reads for reference {}", ref_name);

                        output_stats = write_consensus_reads(
                            &sorted_reads,
                            &mut writer,
                            levels,
                            &rm,
                            read_structure,
                            maximum_reads_before_downsampling,
                            merge_strategy,
                            processing_threads,
                        );

                    }
                    ReadOutputApproach::Correct => {
                        info!("writing reads for reference {}", ref_name);

                        output_stats = write_corrected_reads(
                            &sorted_reads,
                            &mut writer,
                            levels,
                            &rm,
                            read_structure,
                        );
                    }
                }
                
            }
        }
        reference_stats.push(CollapseReferenceStats {
            reference: ref_name,
            filtering,
            reads_after_umi_correction: read_count,
            umi_levels,
            output: output_stats,
        });
    }

    writer.close().unwrap();

    CollapseRunStats {
        mode: match output_approach {
            ReadOutputApproach::Collapse => "collapse".to_string(),
            ReadOutputApproach::Correct => "correct-only".to_string(),
        },
        references: reference_stats,
        elapsed_seconds: start.elapsed().as_secs_f64(),
    }
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

        // TODO: the header uses an in-order map for storing reference sequences, think about this
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

#[derive(Clone, Copy, Debug)]
pub struct AlignmentFilterConfig {
    pub min_aligned_bases: usize,
    pub min_aligned_identical_proportion: f64,
}

impl AlignmentFilterConfig {
    pub fn new(min_aligned_bases: usize, min_aligned_identical_proportion: f64) -> Self {
        assert!(
            min_aligned_identical_proportion.is_finite()
                && (0.0..=1.0).contains(&min_aligned_identical_proportion),
            "Minimum aligned identity must be between 0.0 and 1.0"
        );
        Self {
            min_aligned_bases,
            min_aligned_identical_proportion,
        }
    }

    fn for_reference(&self, reference: &ReferenceRecord) -> AlignmentCheck {
        let alignable_reference_bases = reference
            .sequence
            .as_bytes()
            .iter()
            .filter(|base| **base > 59 && **base != FASTA_N)
            .count();

        AlignmentCheck {
            min_aligned_bases: self.min_aligned_bases.min(alignable_reference_bases),
            min_aligned_identical_proportion: self.min_aligned_identical_proportion,
        }
    }
}

impl Default for AlignmentFilterConfig {
    fn default() -> Self {
        Self::new(45, 0.8)
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
        let mut alignment_count = 0;
        let mut alignable_bases = 0;

        read.aligned_read
            .read_aligned
            .iter()
            .zip(read.aligned_read.reference_aligned.iter())
            .for_each(|(x, y)| {
                if *y > 59 && *x > 59 && y != &FASTA_N {
                    alignable_bases += 1;
                    if x == y {
                        alignment_count += 1;
                    }
                }
            });

        if alignable_bases == 0 {
            return false;
        }

        let ret = (alignment_count as f64 / alignable_bases as f64
            >= self.min_aligned_identical_proportion)
            && (alignable_bases >= self.min_aligned_bases);
        ret
    }
}

/// We want to be extra confident in the alignments around our 'tags'.
/// This filters out reads where we have mismatches and gaps around the
/// degenerate sequences we recover
#[allow(dead_code)]
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

        read.aligned_read
            .read_aligned
            .iter()
            .zip(read.aligned_read.reference_aligned.iter())
            .for_each(|(read_base, reference_base)| {
                // we're at the end of the countdown window - check the mating proportion
                if count_down_check == 0 {
                    count_down_check = usize::MAX;
                    let lookback_length = min(pushed_binary_comp.len(), self.flanking_window_size);
                    let sum: u32 = pushed_binary_comp
                        [pushed_binary_comp.len() - lookback_length..pushed_binary_comp.len()]
                        .iter()
                        .sum();
                    let matching_prop = sum as f64 / lookback_length as f64;
                    pushed_binary_comp.clear();
                    if matching_prop < self.min_flanking_indentity {
                        ret = false;
                    }
                }
                //
                else if *reference_base > 58 && reference_base != &FASTA_N {
                    count_down_check -= 1;
                    if read_base == reference_base {
                        pushed_binary_comp.push(1)
                    } else {
                        pushed_binary_comp.push(0)
                    }
                }
                // lookback case for start of Ns
                else if *reference_base < 59 && pushed_binary_comp.len() > 0 {
                    let lookback_length = min(pushed_binary_comp.len(), self.flanking_window_size);
                    let sum: u32 = pushed_binary_comp
                        [pushed_binary_comp.len() - lookback_length..pushed_binary_comp.len()]
                        .iter()
                        .sum();
                    let matching_prop = sum as f64 / lookback_length as f64;
                    pushed_binary_comp.clear();
                    if matching_prop < self.min_flanking_indentity {
                        ret = false;
                    }
                } else if reference_base == &FASTA_N && pushed_binary_comp.len() == 0 {
                    count_down_check = self.flanking_window_size;
                }
            });
        //println!("aligning {}\n{}\n{}",ret,u8s(&read.aligned_read.read_aligned), u8s(&read.aligned_read.reference_aligned));
        ret
    }
}

#[derive(Default, Copy, Clone, Debug)]
pub struct BamReadFiltering {
    total_reads: usize,
    unmapped_flag_reads: usize,
    secondary_flag_reads: usize,
    supplementary_flag_reads: usize,
    failed_alignment_filters: usize,
    failed_alignment_creation: usize,
    duplicate_reads: usize,
    invalid_tags: usize,
}

impl BamReadFiltering {
    pub fn passing_reads(&self) -> usize {
        self.total_reads
            - self.unmapped_flag_reads
            - self.secondary_flag_reads
            - self.supplementary_flag_reads
            - self.failed_alignment_filters
            - self.failed_alignment_creation
            - self.duplicate_reads
            - self.invalid_tags
    }

    fn accumulate(&mut self, other: &Self) {
        self.total_reads += other.total_reads;
        self.unmapped_flag_reads += other.unmapped_flag_reads;
        self.secondary_flag_reads += other.secondary_flag_reads;
        self.supplementary_flag_reads += other.supplementary_flag_reads;
        self.failed_alignment_filters += other.failed_alignment_filters;
        self.failed_alignment_creation += other.failed_alignment_creation;
        self.duplicate_reads += other.duplicate_reads;
        self.invalid_tags += other.invalid_tags;
    }

    fn count_flag_exclusion(&mut self, flags: Flags) -> bool {
        if flags.is_unmapped() {
            self.unmapped_flag_reads += 1;
        } else if flags.is_secondary() {
            self.secondary_flag_reads += 1;
        } else if flags.is_supplementary() {
            self.supplementary_flag_reads += 1;
        } else {
            return false;
        }

        true
    }

    pub fn results(&self, filters_counts: &HashMap<String, u64>) {
        let filter_summary = filters_counts
            .iter()
            .map(|x| format!("Name: {} failed {}", x.0.clone(), x.1))
            .join(", ");
        info!(
            "Total reads processed: {}, Unmapped: {}, Secondary: {}, Supplementary: {}, [Failed: {}, Failed alignment filters: {}, Duplicate: {}, Invalid_tags: {}, Passing: {} filter summary {}",
            self.total_reads,
            self.unmapped_flag_reads,
            self.secondary_flag_reads,
            self.supplementary_flag_reads,
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

/// Sorts and filters reads from a BAM file for a specific reference sequence.
///
/// This function reads aligned reads from a BAM file, filters them based on quality criteria,
/// extracts UMI tag information, and writes the valid reads to a temporary sharded file for
/// further processing. Indexed BAMs use a region query; unindexed BAMs are streamed and filtered
/// by reference ID so output from `align` can be passed directly to `collapse`. Only primary,
/// mapped reads that pass alignment filters are retained.
///
/// # Arguments
/// * `bam_file` - Path to the input BAM file
/// * `reference_name` - Name of the reference sequence to process
/// * `reference_manager` - Manager containing reference sequence information
/// * `read_structure` - Configuration defining the sequence layout and UMI structure
/// * `temp_directory` - Temporary directory for intermediate file storage
///
/// # Returns
/// * `SortedReadsFromBam` - Container with optional sharded reader and filtering statistics
///
/// # Filtering Criteria
/// * Excludes unmapped reads (unmapped flag set)
/// * Excludes secondary alignments
/// * Excludes supplementary alignments
/// * Applies alignment quality filters (minimum aligned bases and identity proportion)
/// * Validates UMI tag extraction
///
/// # Examples
/// ```
/// let sorted_reads = sort_reads_from_bam_file(
///     &"input.bam".to_string(),
///     &"chr1".to_string(),
///     &reference_manager,
///     &sequence_layout,
///     &mut temp_directory,
/// );
/// ```
pub fn sort_reads_from_bam_file(
    bam_file: &String,
    reference_name: &String,
    reference_manager: &ReferenceManager,
    read_structure: &SequenceLayout,
    temp_directory: &mut InstanceLivedTempDir,
    alignment_filter: &AlignmentFilterConfig,
) -> SortedReadsFromBam {

    let aligned_temp = temp_directory.temp_file("bam.reads.sorted.sharded");

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();

    let mut read_stats = BamReadFiltering::default();

    let reference_config = read_structure.references.get(reference_name).unwrap();
    let alignment_check = alignment_filter.for_reference(reference_config);
    let filters: Vec<(String, &dyn AlignmentFilter)> = vec![
        /*(
            "FlankingDegenerateBaseFilter".to_string(),
            &FlankingDegenerateBaseFilter {
                min_flanking_indentity: 0.50,
                flanking_window_size: 10,
            },
        ),*/
        (
            "AlignmentCheck".to_string(),
            &alignment_check,
        ),
    ];
    let mut filter_counts: HashMap<String, u64> = HashMap::default();
    //filter_counts.insert("FlankingDegenerateBaseFilter".to_string(), 0);
    filter_counts.insert("AlignmentCheck".to_string(), 0);

    let header = reader.read_header().unwrap();
    let bai_path = PathBuf::from(format!("{}.bai", bam_file));
    let index = if bai_path.exists() {
        Some(bai::fs::read(&bai_path).unwrap_or_else(|error| {
            panic!(
                "Unable to read BAM index {}: {}",
                bai_path.display(),
                error
            )
        }))
    } else {
        warn!(
            "No BAM index found at {}; scanning the input for reference '{}'",
            bai_path.display(),
            reference_name
        );
        None
    };
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

        let records: Box<dyn Iterator<Item = std::io::Result<Record>> + '_> =
            if let Some(index) = index.as_ref() {
                let region = reference_name.parse().expect("Unable to parse chromosome");
                Box::new(
                    reader
                        .query(&header, index, &region)
                        .expect("Unable to parse out region information"),
                )
            } else {
                let reference_lookup = ReferenceLookupTable::new(reference_manager, &header);
                let bam_reference_id = *reference_lookup
                    .bam_reference_name_to_id
                    .get(reference_name)
                    .unwrap_or_else(|| {
                        panic!(
                            "Reference '{}' is not present in BAM header",
                            reference_name
                        )
                    });

                Box::new(reader.records().filter_map(move |result| match result {
                    Ok(record) => match record.reference_sequence_id().transpose() {
                        Ok(Some(reference_id)) if reference_id == bam_reference_id => {
                            Some(Ok(record))
                        }
                        Ok(_) => None,
                        Err(error) => Some(Err(error)),
                    },
                    Err(error) => Some(Err(error)),
                }))
            };

        warn!("fetching reads for reference {} ", reference_name);
        let mut read_count = 0;
        let mut last_read_name: Option<Vec<u8>> = None;
        
        for result in records {
            read_stats.total_reads += 1;
            if read_stats.total_reads % 1000000 == 0 {
                read_stats.results(&filter_counts);
            }

            let record = match result {
                Ok(x) => {
                    last_read_name = Some(x.name().unwrap().to_vec());
                    x
                }
                Err(x) => {
                    println!("Read: {}", read_count);
                    println!("Read: {}", u8s(last_read_name.as_ref().unwrap_or(&"UKNOWN".as_bytes().to_vec())));
                    panic!("Unable to read record: {:?}", x);
                }
            };
            
            read_count += 1;
            
            if !read_stats.count_flag_exclusion(record.flags()) {
                let read = create_sorted_read_container(
                    reference_name,
                    &reference_manager,
                    &mut read_stats,
                    &reference_sequence_id,
                    &reference_sequence,
                    reference_config,
                    &record,
                );

                match read {
                    Some(x) => {
                        let survives_filtering = filters
                            .iter()
                            .map(|t| {
                                let x = t.1.keep(&x);
                                if !x {
                                    filter_counts.insert(
                                        t.0.clone(),
                                        filter_counts.get(&t.0).unwrap_or(&0) + 1,
                                    );
                                }
                                x
                            })
                            .filter(|b| !*b)
                            .count()
                            == 0;
                        if survives_filtering {
                            sender.send(x).unwrap();
                        } else {
                            read_stats.failed_alignment_filters += 1;
                        }
                    }
                    None => {
                        read_stats.failed_alignment_creation += 1;
                    }
                }
            }
        }
        sender.finished().unwrap();
        sharded_output.finish().unwrap();
    }
    read_stats.results(&filter_counts);
    if read_stats.passing_reads() > 0 {
        SortedReadsFromBam {
            bam: Some(ShardReader::open(aligned_temp).unwrap()),
            read_stats,
        }
    } else {
        SortedReadsFromBam {
            bam: None,
            read_stats,
        }
    }
}

/// Creates a sorted read container from a BAM record and reference sequence information.
///
/// This function processes a single BAM record, performs sequence alignment operations,
/// extracts UMI tag sequences, and creates a `SortingReadSetContainer` for downstream
/// processing. The function handles sequence alignment, tag extraction, and quality score
/// preservation.
///
/// # Arguments
/// * `reference_name` - Name of the reference sequence
/// * `reference_manager` - Manager containing reference sequence information
/// * `_read_stats` - Mutable reference to BAM read filtering statistics (unused)
/// * `reference_sequence_id` - ID of the reference sequence in the manager
/// * `reference_sequence` - The reference sequence as bytes
/// * `reference_config` - Configuration for this reference including UMI structure
/// * `record` - The BAM record to process
///
/// # Returns
/// * `Some(SortingReadSetContainer)` - if valid tags were extracted
/// * `None` - if tag extraction failed (invalid tags); the caller counts this as
///   `failed_alignment_creation`
///
/// # Process
/// 1. Extracts sequence, alignment position, CIGAR string, and quality scores from BAM record
/// 2. Recovers soft-clipped alignment sequences using realignment
/// 3. Stretches alignment to match reference sequence
/// 4. Extracts tagged sequences and validates tag extraction
/// 5. Creates alignment result with all necessary information
fn create_sorted_read_container(
    reference_name: &String,
    reference_manager: &&ReferenceManager,
    _read_stats: &mut BamReadFiltering,
    reference_sequence_id: &&usize,
    reference_sequence: &Vec<u8>,
    reference_config: &ReferenceRecord,
    record: &Record,
) -> Option<SortingReadSetContainer> {
    let seq: Vec<u8> = record.sequence().iter().collect();
    let one_based_start_pos = record.alignment_start().unwrap().unwrap().get();
    let zero_based_start_pos = one_based_start_pos - 1;
    let cigar = record.cigar();
    let read_name = record.name().unwrap();
    let read_quals = match record.quality_scores().iter().collect::<Vec<u8>>() {
        qualities if qualities.is_empty() => None,
        qualities => Some(qualities),
    };
    let ref_slice = reference_sequence.as_slice();

    let aligned_read = recover_soft_clipped_align_sequences(
        &seq,
        one_based_start_pos,
        &cigar.iter().map(|x| x.unwrap()).collect(),
        &SoftClipResolution::Realign,
        ref_slice,
    );

    let stretched_alignment = stretch_sequence_to_alignment(
        &aligned_read.aligned_ref,
        &reference_manager
            .references
            .get(reference_sequence_id)
            .unwrap()
            .sequence,
    );

    let extracted_tags = extract_tagged_sequences(&aligned_read.aligned_read, &stretched_alignment);

    let (invalid_tag, read_tags_ordered) =
        extract_tag_sequences(reference_config, extracted_tags);

    if !invalid_tag {
        Some(SortingReadSetContainer {
            ordered_sorting_keys: vec![], // for future use during sorting
            ordered_unsorted_keys: read_tags_ordered, // the current unsorted tag collection
            aligned_read: AlignmentResult {
                reference_name: reference_name.clone(),
                read_aligned: aligned_read.aligned_read,
                read_quals,
                cigar_string: cigar
                    .iter()
                    .map(|op| AlignmentTag::from(op.unwrap()))
                    .collect(),
                path: vec![],
                score: 0.0,
                reference_start: zero_based_start_pos,
                read_start: 0,
                reference_aligned: aligned_read.aligned_ref,
                read_name: String::from_utf8(read_name.to_vec()).unwrap(),
                bounding_box: None,
            },
        })
    } else {
        None
    }
}

/// Creates a vector of byte sequences from a text file.
///
/// This function reads a text file line by line, converts each line to a byte vector,
/// and optionally applies reverse complement transformation to DNA sequences.
/// Each line in the file is treated as a separate sequence.
///
/// # Arguments
/// * `filename` - Path to the input file containing sequences (one per line)
/// * `reverse_comp` - Whether to apply reverse complement transformation to sequences
///
/// # Returns
/// * `Vec<Vec<u8>>` - Vector of byte sequences, one for each line in the file
///
/// # Panics
/// * If the input file cannot be opened
///
/// # Examples
/// ```
/// // Read sequences without reverse complement
/// let sequences = create_input_set("sequences.txt", &false);
///
/// // Read sequences with reverse complement
/// let sequences = create_input_set("sequences.txt", &true);
/// ```
fn create_input_set(filename: &str, reverse_comp: &bool) -> Vec<Vec<u8>> {
    let raw_reader = BufReader::new(
        File::open(filename).expect(&format!("Unable to open input file {}", filename)),
    );
    let mut input_set = Vec::new();
    for line in raw_reader.lines() {
        let mut bytes = line.unwrap().into_bytes();
        if *reverse_comp {
            bytes = reverse_complement(&bytes);
        }
        input_set.push(bytes);
    }
    input_set
}

/// Extracts a list of known UMI sequences from a configuration file.
///
/// This function reads known UMI sequences from a file specified in the UMI configuration.
/// The sequences can optionally be reverse complemented based on the configuration settings.
/// This is typically used to create allowlists of valid UMI sequences for sequence correction.
///
/// # Arguments
/// * `umi_type` - UMI configuration containing file path and processing options
/// * `_starting_nmer_size` - Starting n-mer size parameter (currently unused)
///
/// # Returns
/// * `Vec<Vec<u8>>` - Vector of known UMI sequences as byte vectors
///
/// # Panics
/// * If no file is specified in the UMI configuration
/// * If the specified file cannot be opened
///
/// # Examples
/// ```
/// use clique::collapse::extract_known_list;
/// 
/// let known_sequences = extract_known_list(&umi_config, &8);
/// ```
pub fn extract_known_list(
    umi_type: &UMIConfiguration,
    _starting_nmer_size: &usize,
) -> Vec<Vec<u8>> {
    let filename = umi_type.file.clone().unwrap();
    let filename = filename.as_str();

    info!(
        "Reading known list from file {}; large files may take a long time",
        filename
    );

    let rev_comp = umi_type.reverse_complement_sequences.unwrap_or(false);

    let input_set = create_input_set(filename, &rev_comp);
    for (index, sequence) in input_set.iter().enumerate() {
        assert_eq!(
            sequence.len(),
            umi_type.length,
            "Allowlist '{}' line {} has length {}; expected {}",
            filename,
            index + 1,
            sequence.len(),
            umi_type.length
        );
    }
    input_set
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct KnownLookupKey {
    filename: String,
    length: usize,
    max_distance: usize,
    reverse_complement_sequences: bool,
    levenshtein_distance: bool,
}

impl KnownLookupKey {
    pub fn from_config(config: &UMIConfiguration) -> Self {
        Self {
            filename: config
                .file
                .clone()
                .expect("KnownTag UMI must specify an allowlist file"),
            length: config.length,
            max_distance: config.max_distance,
            reverse_complement_sequences: config.reverse_complement_sequences.unwrap_or(false),
            levenshtein_distance: config.uses_levenshtein_distance(),
        }
    }
}

pub struct LookupCollection {
    pub ret_trie: HashMap<KnownLookupKey, Trie>,
    pub ret_known_lookup: HashMap<KnownLookupKey, KnownList>,
}

/// Creates lookup collections for known UMI sequences from the sequence layout configuration.
///
/// This function processes all UMI configurations in the sequence layout and creates appropriate
/// lookup data structures (tries for Levenshtein distance matching or hash maps for exact/Hamming
/// distance matching) for efficient sequence correction and validation.
///
/// # Arguments
/// * `read_structure` - Sequence layout containing UMI configurations for all references
///
/// # Returns
/// * `LookupCollection` - Collection containing tries for Levenshtein matching and known lists for other matching
///
/// # Process
/// 1. Iterates through all reference sequences in the layout
/// 2. For each UMI configuration with a known sequence file:
///    - If Levenshtein distance is enabled: creates a Trie for fuzzy matching
///    - Otherwise: creates a KnownList for exact/Hamming distance matching
/// 3. Populates lookup structures with sequences from configuration files
///
/// # Examples
/// ```
/// let lookup_collection = get_known_level_lookups(&sequence_layout);
/// ```
fn get_known_level_lookups(read_structure: &SequenceLayout) -> LookupCollection {
    let mut ret_trie: HashMap<KnownLookupKey, Trie> = HashMap::new();
    let mut ret_known_lookup: HashMap<KnownLookupKey, KnownList> = HashMap::new();

    read_structure
        .references
        .iter()
        .for_each(|(reference_name, reference)| {
            reference
                .umi_configurations
                .iter()
                .for_each(|(umi_name, config)| match (&config.sort_type, &config.file) {
                    (UMISortType::KnownTag, None) => {
                        panic!(
                            "KnownTag UMI '{}' for reference '{}' must specify an allowlist file",
                            umi_name, reference_name
                        );
                    }
                    (UMISortType::KnownTag, Some(filename)) => {
                        let lookup_key = KnownLookupKey::from_config(config);
                        if config.uses_levenshtein_distance() {
                            let known_lookup = extract_known_list(config, &8);

                            let mut trie = Trie::new(config.length);

                            info!(
                                "creating known lookup tree for file {}",
                                filename
                            );
                            known_lookup.iter().for_each(|sequence| {
                                trie.insert(sequence, None, &config.max_distance);
                            });
                            debug!("creating kn");
                            ret_trie.insert(lookup_key, trie);
                        } else {
                            ret_known_lookup.insert(lookup_key, KnownList::new(config));
                        }
                    }
                    (UMISortType::DegenerateTag, _) => {}
                })
        });

    LookupCollection {
        ret_trie,
        ret_known_lookup,
    }
}

/// Sorts reads by UMI tags and performs sequence correction within groups.
///
/// This function processes reads in batches where all previous UMI tags match, then determines
/// cliques of similar sequences within each batch and corrects them to consensus sequences.
/// For example, when sorting by 10X cell barcodes, reads with the same cell barcode are grouped
/// together, and UMI sequences within each cell are corrected to consensus sequences.
///
/// # Arguments
/// * `temp_directory` - Temporary directory for intermediate file storage
/// * `reader` - Sharded reader containing sorted reads to process
/// * `tag` - UMI configuration specifying the tag to sort by and correction parameters
/// * `iteration` - Current iteration level in the sorting hierarchy
/// * `read_count` - Total number of reads being processed (for progress tracking)
/// * `known_sequence_lists` - Lookup collections for known sequence correction
///
/// # Returns
/// * `(usize, ShardReader<SortingReadSetContainer>)` - Tuple of:
///   - Number of reads written to output
///   - Sharded reader for the next processing level
///
/// # Process
/// 1. Groups reads by matching UMI tags from previous sorting levels
/// 2. For each group, applies sequence correction based on sort type:
///    - `KnownTag`: Corrects to known sequences using Levenshtein or Hamming distance
///    - `DegenerateTag`: Performs clustering and consensus calling on similar sequences
/// 3. Writes corrected reads to temporary sharded output
/// 4. Reports processing statistics
///
/// # Examples
/// ```
/// let (output_count, next_reader) = sort_level(
///     &mut temp_dir,
///     &input_reader,
///     &umi_config,
///     &iteration,
///     &read_count,
///     &mut lookup_collections,
/// );
/// ```
pub fn sort_level(
    temp_directory: &mut InstanceLivedTempDir,
    reader: &ShardReader<SortingReadSetContainer>,
    tag: &UMIConfiguration,
    iteration: &usize,
    read_count: &usize,
    known_sequence_lists: &mut LookupCollection,
) -> (usize, ShardReader<SortingReadSetContainer>) {
    info!("Sorting level {}", tag.symbol);

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
        1000000 // TODO make this a constant somewhere
    };
    info!("Starting to sort {:?} level {}", tag.sort_type, tag.symbol);

    let mut current_sorting_bin: Option<SequenceCorrector> = None;

    reader
        .iter_range(&Range::all())
        .unwrap()
        .for_each(|current_read| {
            all_read_count += 1;
            if all_read_count % 10000 == 0 {
                bar.as_mut().map(|b| b.set_position(all_read_count as u64));
            }
            let mut current_read = current_read.unwrap();
            let next_last_read = current_read.clone();

            match current_sorting_bin.as_mut() {
                None => {
                    let mut bin = SequenceCorrector::new(
                        temp_directory.temp_file(format!("{}.fasta", tag.order).as_str()),
                        &maximum_reads_per_bin,
                        tag.clone(),
                    );
                    bin.push(current_read);
                    current_sorting_bin = Some(bin);
                }

                Some(bin) => {
                    let reads_equal =
                        last_read.as_ref().unwrap().cmp(&mut current_read) == Ordering::Equal;

                    match reads_equal {
                        true => {
                            // add the current read to the bin
                            bin.push(current_read);
                        }
                        false => {
                            // write the previous bin, and add the current read to the next bin
                            match tag.sort_type {
                                UMISortType::KnownTag => {
                                    if tag.uses_levenshtein_distance() {
                                        output_reads += bin.close_trie_known_list(
                                            &mut sender,
                                            tag,
                                            known_sequence_lists,
                                        );
                                    } else {
                                        output_reads += bin.close_hamming_known_list(
                                            &mut sender,
                                            tag,
                                            known_sequence_lists,
                                        );
                                    }
                                }
                                UMISortType::DegenerateTag => {
                                    output_reads += bin.close_degenerate_list(&mut sender);
                                }
                            }
                            bin.push(current_read);
                        }
                    }
                }
            };

            last_read = Some(next_last_read);
        });

    // Mark the scan phase complete before the potentially long final correction/rewrite.
    if let Some(progress_bar) = bar.as_mut() {
        progress_bar.set_position(all_read_count as u64);
        progress_bar.finish();
    }
    info!(
        "Sorting complete for {:?} level {}; applying corrections and rewriting {} reads",
        tag.sort_type, tag.symbol, all_read_count
    );

    match current_sorting_bin {
        None => {}
        Some(mut bin) => match tag.sort_type {
            UMISortType::KnownTag => {
                if tag.uses_levenshtein_distance() {
                    output_reads +=
                        bin.close_trie_known_list(&mut sender, tag, known_sequence_lists);
                } else {
                    output_reads +=
                        bin.close_hamming_known_list(&mut sender, tag, known_sequence_lists);
                }
            }
            UMISortType::DegenerateTag => {
                output_reads += bin.close_degenerate_list(&mut sender);
            }
        },
    }

    info!("For tag {} ({:?}, iteration {}) we processed {} reads, of which {} were passed to the next level",
        &tag.symbol,
        &tag.sort_type,
        iteration,
        all_read_count,
        output_reads);

    sender.finished().unwrap();

    sharded_output.finish().unwrap();

    (output_reads, ShardReader::open(aligned_temp).unwrap())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::alignment_matrix::AlignmentResult;
    use crate::alignment_manager::OutputAlignmentWriter;
    use std::collections::BTreeMap;

    const FASTA_A: u8 = b'A';
    #[allow(dead_code)]
    const FASTA_T: u8 = b'T';

    #[test]
    fn test_collapse_summary_aggregates_filtering_umi_and_consensus_counts() {
        let stats = CollapseRunStats {
            mode: "collapse".to_string(),
            references: vec![CollapseReferenceStats {
                reference: "reference_a".to_string(),
                filtering: BamReadFiltering {
                    total_reads: 100,
                    unmapped_flag_reads: 5,
                    secondary_flag_reads: 2,
                    supplementary_flag_reads: 1,
                    failed_alignment_filters: 4,
                    failed_alignment_creation: 3,
                    duplicate_reads: 0,
                    invalid_tags: 0,
                },
                reads_after_umi_correction: 70,
                umi_levels: vec![UmiLevelRunStats {
                    reference: "reference_a".to_string(),
                    symbol: '0',
                    sort_type: "KnownTag".to_string(),
                    input_reads: 85,
                    output_reads: 70,
                }],
                output: ConsensusWriteStats {
                    input_reads: 70,
                    groups_attempted: 10,
                    reads_selected: 60,
                    reads_downsampled: 10,
                    output_reads: 9,
                    failed_groups: 1,
                    reads_per_consensus: BTreeMap::from([
                        (1, 4),
                        (2, 2),
                        (4, 2),
                        (8, 1),
                    ]),
                },
            }],
            elapsed_seconds: 2.5,
        };

        let rendered = stats.summary().render();
        assert!(rendered.contains("| All         | 100         | 85"));
        assert!(rendered.contains("| reference_a / 0 | KnownTag | 85"));
        assert!(rendered.contains("Reads per consensus (rc tag)"));
        assert!(rendered.contains("| 1                   | 4"));
        assert!(rendered.contains(&"#".repeat(CONSENSUS_HISTOGRAM_WIDTH)));
        assert!(rendered.contains("| All   | collapse | 2.50"));
    }

    #[test]
    fn test_consensus_size_histogram_uses_compact_power_of_two_bins() {
        assert_eq!(consensus_size_bin(1), (1, 1));
        assert_eq!(consensus_size_bin(2), (2, 2));
        assert_eq!(consensus_size_bin(3), (3, 4));
        assert_eq!(consensus_size_bin(4), (3, 4));
        assert_eq!(consensus_size_bin(5), (5, 8));
        assert_eq!(consensus_size_bin(8), (5, 8));
        assert_eq!(consensus_size_bin(9), (9, 16));

        let histogram = consensus_size_histogram(&BTreeMap::from([
            (1, 5),
            (2, 3),
            (3, 2),
            (4, 1),
            (8, 2),
        ]));
        let mut summary = RunSummary::new("Test summary");
        summary.add_table(histogram);
        let rendered = summary.render();
        assert!(rendered.contains("| 1                   | 5"));
        assert!(rendered.contains("| 2                   | 3"));
        assert!(rendered.contains("| 3-4                 | 3"));
        assert!(rendered.contains("| 5-8                 | 2"));
    }

    fn known_tag_layout(file: Option<&str>, levenshtein_distance: Option<bool>) -> SequenceLayout {
        let mut umi_configurations = BTreeMap::new();
        umi_configurations.insert(
            "cell_id".to_string(),
            UMIConfiguration {
                symbol: '0',
                file: file.map(str::to_string),
                reverse_complement_sequences: None,
                sort_type: UMISortType::KnownTag,
                length: 16,
                order: 0,
                pad: None,
                max_distance: 1,
                maximum_subsequences: None,
                max_gaps: None,
                minimum_collapsing_difference: None,
                levenshtein_distance,
            },
        );

        let mut references = BTreeMap::new();
        references.insert(
            "reference".to_string(),
            ReferenceRecord {
                sequence: "0000000000000000".to_string(),
                umi_configurations,
                targets: vec![],
                target_types: vec![],
                target_locations: Some(vec![]),
                prime_edits: BTreeMap::new(),
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

    fn unindexed_bam_layout() -> SequenceLayout {
        let references = [
            ("reference_a", "A".repeat(50)),
            ("reference_b", "C".repeat(50)),
        ]
        .iter()
        .map(|(name, sequence)| {
            (
                name.to_string(),
                ReferenceRecord {
                    sequence: sequence.clone(),
                    umi_configurations: BTreeMap::new(),
                    targets: vec![],
                    target_types: vec![],
                    target_locations: Some(vec![]),
                    prime_edits: BTreeMap::new(),
                },
            )
        })
        .collect();

        SequenceLayout {
            aligner: None,
            merge: None,
            reads: vec![],
            known_strand: true,
            references,
        }
    }

    fn aligned_read(reference_name: &str, read_name: &str, base: u8) -> SortingReadSetContainer {
        SortingReadSetContainer::empty_tags(AlignmentResult {
            reference_name: reference_name.to_string(),
            read_name: read_name.to_string(),
            reference_aligned: vec![base; 50],
            read_aligned: vec![base; 50],
            read_quals: Some(vec![40; 50]),
            cigar_string: vec![AlignmentTag::MatchMismatch(50)],
            path: vec![],
            score: 50.0,
            reference_start: 0,
            read_start: 0,
            bounding_box: None,
        })
    }

    #[test]
    fn test_sort_reads_accepts_unindexed_unsorted_bam() {
        let layout = unindexed_bam_layout();
        let reference_manager = ReferenceManager::from_yaml_input(&layout, 8, 4);
        let output_directory = tempfile::tempdir().unwrap();
        let bam_path = output_directory.path().join("aligned.bam");
        let bai_path = PathBuf::from(format!("{}.bai", bam_path.display()));
        std::fs::write(&bai_path, b"stale index").unwrap();

        {
            let mut writer = BamFileAlignmentWriter::new(&bam_path, &reference_manager);
            writer
                .write_read(
                    &aligned_read("reference_b", "read_b", b'C'),
                    &HashMap::new(),
                )
                .unwrap();
            writer
                .write_read(
                    &aligned_read("reference_a", "read_a", b'A'),
                    &HashMap::new(),
                )
                .unwrap();
            writer.close().unwrap();
        }

        let bam_file = bam_path.to_string_lossy().into_owned();
        assert!(!bai_path.exists());

        let mut processing_directory = InstanceLivedTempDir::new().unwrap();
        let result = sort_reads_from_bam_file(
            &bam_file,
            &"reference_a".to_string(),
            &reference_manager,
            &layout,
            &mut processing_directory,
            &AlignmentFilterConfig::default(),
        );

        assert_eq!(result.read_stats.total_reads, 1);
        assert_eq!(result.read_stats.passing_reads(), 1);
        assert!(result.bam.is_some());
    }

    #[test]
    fn test_sort_reads_preserves_nonzero_alignment_start() {
        let layout = unindexed_bam_layout();
        let reference_manager = ReferenceManager::from_yaml_input(&layout, 8, 4);
        let output_directory = tempfile::tempdir().unwrap();
        let bam_path = output_directory.path().join("aligned.bam");

        let mut input_read = aligned_read("reference_a", "read_a", b'A');
        input_read.aligned_read.reference_aligned = vec![b'A'; 45];
        input_read.aligned_read.read_aligned = vec![b'A'; 45];
        input_read.aligned_read.read_quals = Some(vec![40; 45]);
        input_read.aligned_read.cigar_string = vec![AlignmentTag::MatchMismatch(45)];
        input_read.aligned_read.reference_start = 5;

        {
            let mut writer = BamFileAlignmentWriter::new(&bam_path, &reference_manager);
            writer.write_read(&input_read, &HashMap::new()).unwrap();
            writer.close().unwrap();
        }

        let mut processing_directory = InstanceLivedTempDir::new().unwrap();
        let result = sort_reads_from_bam_file(
            &bam_path.to_string_lossy().into_owned(),
            &"reference_a".to_string(),
            &reference_manager,
            &layout,
            &mut processing_directory,
            &AlignmentFilterConfig::default(),
        );

        let sorted_reads = result.bam.unwrap();
        let sorted_read = sorted_reads.iter().unwrap().next().unwrap().unwrap();
        assert_eq!(sorted_read.aligned_read.reference_start, 5);
    }

    #[test]
    fn test_known_lookup_defaults_to_levenshtein() {
        let filename = "test_data/subset_barcode_list_500.txt";
        let layout = known_tag_layout(Some(filename), None);
        let config = layout.references["reference"].umi_configurations["cell_id"].clone();
        let lookups = get_known_level_lookups(&layout);

        assert!(lookups
            .ret_trie
            .contains_key(&KnownLookupKey::from_config(&config)));
        assert!(lookups.ret_known_lookup.is_empty());
    }

    #[test]
    fn test_known_lookup_uses_hamming_when_explicitly_disabled() {
        let filename = "test_data/subset_barcode_list_500.txt";
        let layout = known_tag_layout(Some(filename), Some(false));
        let config = layout.references["reference"].umi_configurations["cell_id"].clone();
        let lookups = get_known_level_lookups(&layout);

        assert!(lookups.ret_trie.is_empty());
        assert!(lookups
            .ret_known_lookup
            .contains_key(&KnownLookupKey::from_config(&config)));
    }

    #[test]
    fn test_known_lookup_keeps_distinct_configurations_for_same_file() {
        let filename = "test_data/subset_barcode_list_500.txt";
        let mut layout = known_tag_layout(Some(filename), None);
        let mut second_reference = layout.references["reference"].clone();
        second_reference
            .umi_configurations
            .get_mut("cell_id")
            .unwrap()
            .max_distance = 2;
        layout
            .references
            .insert("second_reference".to_string(), second_reference);

        let lookups = get_known_level_lookups(&layout);

        assert_eq!(lookups.ret_trie.len(), 2);
    }

    #[test]
    #[should_panic(expected = "must specify an allowlist file")]
    fn test_known_lookup_requires_file() {
        get_known_level_lookups(&known_tag_layout(None, None));
    }

    /// Generates a consensus sequence from multiple input sequences.
    ///
    /// This function takes a vector of byte sequences (representing DNA sequences) and generates
    /// a consensus sequence by determining the most frequent base at each position. Special rules
    /// apply for tie-breaking: 'N' (unknown) and '-' (gap) characters are deprioritized.
    ///
    /// # Arguments
    /// * `input` - Vector of byte sequences, all must be the same length
    ///
    /// # Returns
    /// * `Vec<u8>` - Consensus sequence as a byte vector
    ///
    /// # Panics
    /// * If input sequences are not all the same length
    ///
    /// # Consensus Rules
    /// 1. Most frequent base at each position wins
    /// 2. In case of ties, prefer non-'N' and non-'-' bases
    /// 3. If tied between 'N' and '-', prefer the new candidate
    ///
    /// # Examples
    /// ```
    /// let sequences = vec![
    ///     "ATCG".as_bytes().to_vec(),
    ///     "ATCG".as_bytes().to_vec(), 
    ///     "GCTA".as_bytes().to_vec(),
    /// ];
    /// let result = consensus(&sequences);
    /// assert_eq!(result, "ATCG".as_bytes());
    /// ```
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
    fn test_bam_read_filtering_passing_reads() {
        let stats = BamReadFiltering {
            total_reads: 100,
            unmapped_flag_reads: 10,
            secondary_flag_reads: 5,
            supplementary_flag_reads: 0,
            failed_alignment_filters: 3,
            failed_alignment_creation: 2,
            duplicate_reads: 1,
            invalid_tags: 4,
        };
        assert_eq!(stats.passing_reads(), 75);
    }

    #[test]
    fn test_bam_read_filtering_all_passing() {
        let stats = BamReadFiltering {
            total_reads: 50,
            unmapped_flag_reads: 0,
            secondary_flag_reads: 0,
            supplementary_flag_reads: 0,
            failed_alignment_filters: 0,
            failed_alignment_creation: 0,
            duplicate_reads: 0,
            invalid_tags: 0,
        };
        assert_eq!(stats.passing_reads(), 50);
    }

    #[test]
    fn test_bam_read_filtering_none_passing() {
        let stats = BamReadFiltering {
            total_reads: 10,
            unmapped_flag_reads: 4,
            secondary_flag_reads: 3,
            supplementary_flag_reads: 0,
            failed_alignment_filters: 1,
            failed_alignment_creation: 0,
            duplicate_reads: 1,
            invalid_tags: 1,
        };
        assert_eq!(stats.passing_reads(), 0);
    }

    #[test]
    fn test_bam_read_filtering_default() {
        let stats = BamReadFiltering::default();
        assert_eq!(stats.total_reads, 0);
        assert_eq!(stats.passing_reads(), 0);
    }

    #[test]
    fn test_bam_read_filtering_excludes_supplementary_once() {
        let mut stats = BamReadFiltering {
            total_reads: 2,
            ..Default::default()
        };

        assert!(stats.count_flag_exclusion(Flags::SUPPLEMENTARY));
        assert!(stats.count_flag_exclusion(
            Flags::SECONDARY | Flags::SUPPLEMENTARY,
        ));

        assert_eq!(stats.supplementary_flag_reads, 1);
        assert_eq!(stats.secondary_flag_reads, 1);
        assert_eq!(stats.passing_reads(), 0);
    }

    #[test]
    fn test_bam_read_filtering_keeps_primary_alignment() {
        let mut stats = BamReadFiltering::default();

        assert!(!stats.count_flag_exclusion(Flags::empty()));
        assert_eq!(stats.supplementary_flag_reads, 0);
    }

    #[test]
    fn test_consensus_all_same() {
        let seqs = vec![
            b"ACGT".to_vec(),
            b"ACGT".to_vec(),
            b"ACGT".to_vec(),
        ];
        assert_eq!(consensus(&seqs), b"ACGT".to_vec());
    }

    #[test]
    fn test_consensus_majority_wins() {
        let seqs = vec![
            b"A".to_vec(),
            b"A".to_vec(),
            b"T".to_vec(),
        ];
        assert_eq!(consensus(&seqs), b"A".to_vec());
    }

    #[test]
    fn test_consensus_gap_deprioritized() {
        // Tie between G and -, G should win
        let seqs = vec![
            b"G".to_vec(),
            b"-".to_vec(),
        ];
        assert_eq!(consensus(&seqs), b"G".to_vec());
    }

    #[test]
    fn test_consensus_n_deprioritized() {
        // Tie between A and N, A should win
        let seqs = vec![
            b"N".to_vec(),
            b"A".to_vec(),
        ];
        assert_eq!(consensus(&seqs), b"A".to_vec());
    }

    #[test]
    fn test_consensus_single_sequence() {
        let seqs = vec![b"ACGTACGT".to_vec()];
        assert_eq!(consensus(&seqs), b"ACGTACGT".to_vec());
    }

    #[test]
    fn test_alignment_check() {
        let alignment_check = AlignmentCheck {
            min_aligned_bases: 10,
            min_aligned_identical_proportion: 0.8,
        };

        let fake_read_alignment = SortingReadSetContainer {
            ordered_sorting_keys: vec![],
            ordered_unsorted_keys: Default::default(),
            aligned_read: AlignmentResult {
                reference_name: "".to_string(),
                read_name: "".to_string(),
                reference_aligned: vec![
                    FASTA_A, FASTA_A, FASTA_A, FASTA_A, FASTA_A, FASTA_A, FASTA_A, FASTA_A,
                    FASTA_A, FASTA_A, FASTA_A, FASTA_A,
                ],
                read_aligned: vec![
                    FASTA_A, FASTA_A, FASTA_A, FASTA_A, FASTA_A, FASTA_A, FASTA_A, FASTA_A,
                    FASTA_A, FASTA_A, FASTA_A, FASTA_A,
                ],
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
    fn test_default_alignment_filter_supports_short_references() {
        let reference = ReferenceRecord {
            sequence: "0000AAAAAAAAAAAA".to_string(),
            umi_configurations: BTreeMap::new(),
            targets: vec![],
            target_types: vec![],
            target_locations: Some(vec![]),
            prime_edits: BTreeMap::new(),
        };
        let alignment_check = AlignmentFilterConfig::default().for_reference(&reference);
        let fake_read_alignment = SortingReadSetContainer::empty_tags(AlignmentResult {
            reference_name: "short".to_string(),
            read_name: "read".to_string(),
            reference_aligned: vec![FASTA_A; 12],
            read_aligned: vec![FASTA_A; 12],
            read_quals: None,
            cigar_string: vec![],
            path: vec![],
            score: 0.0,
            reference_start: 0,
            read_start: 0,
            bounding_box: None,
        });

        assert_eq!(alignment_check.min_aligned_bases, 12);
        assert!(alignment_check.keep(&fake_read_alignment));
    }

    #[test]
    fn test_alignment_filter_thresholds_are_configurable() {
        let config = AlignmentFilterConfig::new(20, 0.95);

        assert_eq!(config.min_aligned_bases, 20);
        assert_eq!(config.min_aligned_identical_proportion, 0.95);
    }

    #[test]
    #[should_panic(expected = "between 0.0 and 1.0")]
    fn test_alignment_filter_rejects_invalid_identity() {
        AlignmentFilterConfig::new(20, 1.1);
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

}
