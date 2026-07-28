//! High-level alignment entry points: orient each read, pick the best-matching
//! reference among the candidates (a k-mer vote fast path with an exhaustive
//! fallback), run the aligner, extract tags, and hand results to the writer.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};
use std::sync::{Arc, Mutex};

use crate::alignment::alignment_matrix::{
    convex_alignment, create_scoring_record_3d, inversion_alignment, perform_3d_global_traceback, perform_affine_alignment,
    perform_affine_alignment_bandwidth, Alignment, AlignmentResult, AlignmentTag, AlignmentType,
};
use crate::alignment::scoring_functions::{AffineScoring, ConvexScoring, InversionScoring};
use crate::linked_alignment::{
    align_string_with_anchors, find_greedy_non_overlapping_segments, orient_by_longest_segment,
};
use crate::rayon::iter::ParallelBridge;
use crate::rayon::iter::ParallelIterator;
use crate::read_strategies::read_set::ReadIterator;
use crate::reference::fasta_reference::ReferenceManager;
use crate::reference::discriminating::DiscriminatingClassifier;
use crate::reference::idf::IdfIndex;
use crate::reference::poa::PoaGraph;
use crate::run_summary::{RunSummary, SummaryTable};
use ndarray::Ix3;
use std::time::Instant;
use bio::alignment::AlignmentOperation;
use crate::merger::{MergedReadSequence, UnifiedRead};

use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation::*;

use crate::read_strategies::sequence_layout::SequenceLayout;

use crate::alignment_manager::BamFileAlignmentWriter;
use crate::alignment_manager::OutputAlignmentWriter;
use crate::read_strategies::read_disk_sorter::SortingReadSetContainer;
use consensus::consensus_builders::get_reference_alignment_rate;
use extractor::extract_tagged_sequences;
use itertools::Itertools;
use utils::read_utils::{reverse_complement, u8s};
use ::{Aligner as RustAligner, FASTA_UNSET};


/// Compute the affine alignment score between `a` and `b` with the given substitution,
/// gap-open, and gap-extend penalties.
#[allow(dead_code)]
fn rust_bio_alignment(
    read: &[u8],
    reference: &[u8],
    mismatch: &i32,
    gap_open: &i32,
    gap_extend: &i32
) -> (Vec<AlignmentOperation>, i32) {
    // A match rewards +1; mismatch and gap penalties come from the (positive-magnitude)
    // parameters, applied as negative costs. `N` and digit-marked UMI positions match
    // any base, consistent with AffineScoring::match_mismatch.
    let mismatch_penalty = -mismatch.abs();
    let score = |a: u8, b: u8| {
        if a == b || a == b'N' || b == b'N' || a.is_ascii_digit() || b.is_ascii_digit() {
            1i32
        } else {
            mismatch_penalty
        }
    };
    let mut aligner = Aligner::with_capacity(read.len(), reference.len(), -gap_open.abs(), -gap_extend.abs(), &score);
    let alignment = aligner.global(reference,read);
    // x is global (target sequence) and y is local (reference sequence)
    (alignment.operations, alignment.score)
}

/// POA is the default router only for panels with at most this many references
/// (above this, building the graph is not worth it — the k-mer router is used).
const POA_DEFAULT_MAX_REFS: usize = 256;
/// A panel is "compact" (references share enough structure for POA to be the
/// right tool) when its graph has at most this many nodes per base of the longest
/// reference. Diverse panels blow past this and fall back to the k-mer router.
const POA_COMPACT_FACTOR: usize = 4;
const REFERENCE_HISTOGRAM_WIDTH: usize = 40;

fn reference_histogram_bar(count: usize, maximum: usize) -> String {
    if count == 0 || maximum == 0 {
        return String::new();
    }
    let bar_length = count
        .saturating_mul(REFERENCE_HISTOGRAM_WIDTH)
        .saturating_add(maximum - 1)
        / maximum;
    "#".repeat(bar_length.max(1).min(REFERENCE_HISTOGRAM_WIDTH))
}

#[derive(Default)]
struct AlignmentCounters {
    input_reads: AtomicUsize,
    aligned_reads: AtomicUsize,
    too_short: AtomicUsize,
    too_long: AtomicUsize,
    failed: AtomicUsize,
    ambiguous_reference_calls: AtomicUsize,
    aligned_by_reference: Mutex<BTreeMap<String, usize>>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AlignmentRunStats {
    pub input_reads: usize,
    pub aligned_reads: usize,
    pub too_short: usize,
    pub too_long: usize,
    pub failed: usize,
    pub ambiguous_reference_calls: usize,
    pub aligned_by_reference: BTreeMap<String, usize>,
    pub elapsed_seconds: f64,
}

impl AlignmentRunStats {
    pub fn summary(&self) -> RunSummary {
        let alignment_rate = if self.input_reads == 0 {
            0.0
        } else {
            100.0 * self.aligned_reads as f64 / self.input_reads as f64
        };
        let mut overall = SummaryTable::new(
            "Alignment totals",
            &[
                "Scope",
                "Input reads",
                "Aligned",
                "Aligned (%)",
                "Too short",
                "Too long",
                "Failed",
                "Ambiguous calls",
                "Elapsed (s)",
            ],
        );
        overall.push_row([
            "All".to_string(),
            self.input_reads.to_string(),
            self.aligned_reads.to_string(),
            format!("{:.2}", alignment_rate),
            self.too_short.to_string(),
            self.too_long.to_string(),
            self.failed.to_string(),
            self.ambiguous_reference_calls.to_string(),
            format!("{:.2}", self.elapsed_seconds),
        ]);

        let mut references = SummaryTable::new(
            "Aligned reads by reference",
            &["Reference", "Aligned reads", "Aligned reads (%)", "Histogram"],
        );
        let maximum_reference_count = self
            .aligned_by_reference
            .values()
            .copied()
            .max()
            .unwrap_or(0);
        for (reference, count) in &self.aligned_by_reference {
            let proportion = if self.aligned_reads == 0 {
                0.0
            } else {
                100.0 * *count as f64 / self.aligned_reads as f64
            };
            references.push_row([
                reference.clone(),
                count.to_string(),
                format!("{:.2}", proportion),
                reference_histogram_bar(*count, maximum_reference_count),
            ]);
        }

        let mut summary = RunSummary::new("Alignment run summary");
        summary.add_table(overall);
        summary.add_table(references);
        summary
    }
}

pub fn align_reads(
    read_structure: &SequenceLayout,
    rm: &ReferenceManager,
    output: &Path,
    max_reference_multiplier: &usize,
    min_read_length: &usize,
    read1: &String,
    read2: &String,
    index1: &String,
    index2: &String,
    threads: &usize,
    aligner: &RustAligner,
    use_discriminating: bool,
    discriminating_min_margin: usize,
    use_kmer_idf: bool,
    kmer_idf_min_margin: f64,
    use_poa: bool,
    poa_min_margin: usize,
    no_poa_default: bool,
) -> AlignmentRunStats {
    let read_iterator = ReadIterator::new(
        PathBuf::from(&read1),
        Some(PathBuf::from(&read2)),
        Some(PathBuf::from(&index1)),
        Some(PathBuf::from(&index2)),
    );

    let read_iterator = MergedReadSequence::new(read_iterator, read_structure);

    let writer = BamFileAlignmentWriter::new(&PathBuf::from(output), &rm);

    let output = Arc::new(Mutex::new(writer));

    // setup our thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(*threads)
        .build_global()
        .unwrap();

    // Inversion scoring tuned for `--aligner inversion`. The legacy defaults
    // (inversion_penalty -40 / min_inversion_length 20) never fired; these detect
    // inversions of >= 8 bp. LIMITATION: the aligner's inversion candidate search
    // degrades as the flanking match length grows, so detection is reliable only
    // for short-flank contexts and is NOT yet dependable on long (>=~20 bp flank)
    // amplicon reads — improving that is aligner-algorithm work beyond scoring.
    let my_score = InversionScoring {
        match_score: 10.0,
        mismatch_score: -11.0,
        gap_open: -15.0,
        gap_extend: -5.0,
        inversion_penalty: -10.0,
        min_inversion_length: 8,
    };

    let my_aff_score = AffineScoring {
        match_score: 10.0,
        mismatch_score: -9.0,
        special_character_score: 9.0,
        gap_open: -20.0,
        gap_extend: -2.0,
        final_gap_multiplier: 1.0,
    };

    // `--aligner inversion` routes the final alignment (single + multi reference)
    // through the inversion-aware aligner; detected inversions are emitted as `iv`
    // BAM tags. `--aligner convex` uses a convex (logarithmic) gap penalty instead
    // of affine, so one long indel is preferred over several short ones.
    let use_inversions = matches!(aligner, RustAligner::Inversion);
    let use_convex = matches!(aligner, RustAligner::Convex);
    let convex_score = ConvexScoring {
        match_score: 10.0,
        mismatch_score: -9.0,
        special_character_score: 9.0,
        gap_open: -10.0,
        gap_extend: -3.0,
    };
    if use_inversions {
        info!("Inversion-aware alignment enabled (--aligner inversion)");
    }
    if use_convex {
        info!("Convex-gap alignment enabled (--aligner convex); O(n*m*(n+m)) per read, slower than affine");
    }
    let start = Instant::now();
    let counters = Arc::new(AlignmentCounters::default());

    type SharedStore = Arc<Mutex<Option<Alignment<Ix3>>>>;

    lazy_static! {
        static ref STORE_CLONES: Mutex<Vec<SharedStore>> = Mutex::new(Vec::new());
    }
    thread_local!(static STORE: SharedStore = Arc::new(Mutex::new(None)));

    let max_read_size = (rm.longest_ref + 1) * max_reference_multiplier;
    info!(
        "Longest reference found: {}, max read size set at {}",
        rm.longest_ref, max_read_size
    );

    let alignment_mat: Alignment<Ix3> = create_scoring_record_3d(
        rm.longest_ref + 1,
        max_read_size,
        AlignmentType::Affine,
        false,
    );

    // Optional discriminating-position classifier for near-identical panels
    // (opt-in via --discriminating-classifier). Built once from the whole panel;
    // disabled if there is <2 references or the panel has no discriminating columns.
    let discriminating_classifier = if use_discriminating {
        match DiscriminatingClassifier::from_reference_manager(rm, discriminating_min_margin) {
            Ok(clf) if clf.n_positions() > 0 => {
                info!(
                    "Discriminating-position classifier enabled: anchor '{}', {} discriminating position(s), min margin {}",
                    clf.anchor_name, clf.n_positions(), discriminating_min_margin
                );
                Some(clf)
            }
            Ok(_) => {
                warn!("Discriminating classifier requested but the panel has no discriminating positions; using the default reference search.");
                None
            }
            Err(e) => {
                warn!("Could not build the discriminating classifier ({}); using the default reference search.", e);
                None
            }
        }
    } else {
        None
    };
    let discriminating_classifier = discriminating_classifier.as_ref();

    // Optional IDF-weighted k-mer router (opt-in via --kmer-idf). Built once;
    // disabled if there are <2 references or no discriminating k-mers.
    let idf_index = if use_kmer_idf {
        if rm.references.len() < 2 {
            warn!("IDF router requested but there are fewer than 2 references; using the default reference search.");
            None
        } else {
            let index = IdfIndex::from_reference_manager(rm, kmer_idf_min_margin);
            if index.has_discriminating_kmers() {
                info!(
                    "IDF-weighted k-mer router enabled: {} references, min margin ratio {}",
                    index.n_refs(), kmer_idf_min_margin
                );
                Some(index)
            } else {
                warn!("IDF router requested but no k-mer discriminates the panel; using the default reference search.");
                None
            }
        }
    } else {
        None
    };
    let idf_index = idf_index.as_ref();

    // POA-graph router. It is the DEFAULT for compact multi-reference panels
    // (near-identical / indel-recording, where references share structure), and
    // can be forced with --poa-classifier. It falls back to the k-mer router for
    // large or diverse panels (a graph would explode and the k-mer/IDF router is
    // the right tool there), when another router flag is set, or via
    // --no-poa-default.
    let explicit_router = use_poa || use_kmer_idf || use_discriminating;
    let want_poa = use_poa || (!explicit_router && !no_poa_default);
    let poa_graph = if want_poa && rm.references.len() >= 2 {
        if !use_poa && rm.references.len() > POA_DEFAULT_MAX_REFS {
            info!(
                "Panel has {} references (> {}); using the k-mer router instead of the POA default.",
                rm.references.len(), POA_DEFAULT_MAX_REFS
            );
            None
        } else {
            let compact_limit = POA_COMPACT_FACTOR * rm.longest_ref.max(1);
            match PoaGraph::from_reference_manager(rm, poa_min_margin) {
                // Explicit --poa-classifier bypasses the compactness gate.
                Ok(graph)
                    if graph.is_discriminable() && (use_poa || graph.node_count() <= compact_limit) =>
                {
                    info!(
                        "POA-graph reference router enabled ({}): {} nodes ({} backbone, {} branch) over {} references, min margin {}",
                        if use_poa { "requested" } else { "default for a compact panel" },
                        graph.node_count(), graph.backbone_node_count(), graph.branch_node_count(), graph.n_refs(), poa_min_margin
                    );
                    Some(graph)
                }
                Ok(graph) => {
                    let reason = if !graph.is_discriminable() {
                        "the panel has no branch nodes".to_string()
                    } else {
                        format!("the panel is not compact ({} nodes > {} = {}x the {}-bp longest reference)",
                            graph.node_count(), compact_limit, POA_COMPACT_FACTOR, rm.longest_ref)
                    };
                    if use_poa {
                        warn!("POA requested but {}; using the k-mer router.", reason);
                    } else {
                        info!("Not using the POA default ({}); using the k-mer router.", reason);
                    }
                    None
                }
                Err(e) => {
                    warn!("Could not build the POA graph ({}); using the k-mer router.", e);
                    None
                }
            }
        }
    } else {
        None
    };
    let poa_graph = poa_graph.as_ref();

    read_iterator.par_bridge().for_each(|mut xx: UnifiedRead| {
        counters.input_reads.fetch_add(1, AtomicOrdering::Relaxed);
        STORE.with(|arc_mtx| {
            let mut local_alignment = arc_mtx.lock().unwrap();
            if local_alignment.is_none() {
                *local_alignment = Some(alignment_mat.clone());
                STORE_CLONES.lock().unwrap().push(arc_mtx.clone());
            }

            let name = &String::from_utf8(xx.name().clone()).unwrap();

            let seq_len = xx.seq().len();
            // FASTQ qualities are ASCII Phred+33; strip the offset to raw Phred for the aligned
            // read (the BAM writer re-applies +33 on output).
            let qual = Some(xx.quals.as_ref().unwrap().iter().map(|b| b.saturating_sub(33)).collect::<Vec<u8>>());
            if seq_len < *min_read_length {
                counters.too_short.fetch_add(1, AtomicOrdering::Relaxed);
                debug!(
                    "Skipping read {} because its length {} is below the minimum {}",
                    name, seq_len, min_read_length
                );
            } else if seq_len < max_read_size {
                let aligned = align_to_reference_choices(
                    name,
                    xx.seq(),
                    qual,
                    rm,
                    &true,
                    read_structure,
                    local_alignment.as_mut().unwrap(),
                    &my_aff_score,
                    &my_score,
                    &use_inversions,
                    &convex_score,
                    &use_convex,
                    *max_reference_multiplier as f64,
                    *min_read_length,
                    &seq_len,
                    discriminating_classifier,
                    idf_index,
                    poa_graph,
                );

                match aligned {
                    None => {
                        counters.failed.fetch_add(1, AtomicOrdering::Relaxed);
                        debug!("Unable to create alignment for read {}", name);
                    }
                    Some(alignment_obj) => {
                        let results = alignment_obj.alignment;

                        let _orig_ref_seq = alignment_obj.ref_sequence;
                        let classification = alignment_obj.classification;
                        let idf_info = alignment_obj.idf_info;
                        let poa_info = alignment_obj.poa_info;
                        let ambiguous_reference_call = classification
                            .as_ref()
                            .map(|info| info.ambiguous)
                            .unwrap_or(false)
                            || idf_info
                                .as_ref()
                                .map(|info| info.ambiguous)
                                .unwrap_or(false)
                            || poa_info
                                .as_ref()
                                .map(|info| info.ambiguous)
                                .unwrap_or(false);
                        match results {
                            None => {
                                counters.failed.fetch_add(1, AtomicOrdering::Relaxed);
                                debug!("Unable to create alignment for read {}", name);
                            }
                            Some(mut aln) => {
                                assert_eq!(aln.reference_aligned.len(), aln.read_aligned.len());

                                // Item B(ii): flatten inversion markers out of the CIGAR into
                                // `<len>V+<pos>` events (emitted below as the `iv` tag) so the
                                // read serializes to a single, standard BAM record.
                                let inversion_events = aln.take_inversion_events();

                                let read = SortingReadSetContainer::empty_tags(aln);

                                let extracted_tags = extract_tagged_sequences(
                                    &read.aligned_read.read_aligned,
                                    &read.aligned_read.reference_aligned,
                                );

                                let mut added_tags: HashMap<[u8; 2], String> = HashMap::new();

                                let structure = read_structure
                                    .references
                                    .get(&read.aligned_read.reference_name)
                                    .unwrap();

                                extracted_tags.iter().for_each(|(x, y)| {
                                    structure.umi_configurations.iter().for_each(|xi| {
                                        if xi.1.symbol as u8 == *x {
                                            added_tags.insert([b'e', xi.1.symbol as u8], y.clone());
                                        }
                                    })
                                });

                                added_tags.insert([b'r', b'c'], 1.to_string());

                                // Discriminating-position classifier confidence, when it made the call:
                                // dm = top-2 margin, di = informative discriminating positions covered,
                                // da = 1 if the call was ambiguous (near-tie) else 0.
                                if let Some(info) = &classification {
                                    added_tags.insert([b'd', b'm'], info.margin.to_string());
                                    added_tags.insert([b'd', b'i'], info.informative_positions.to_string());
                                    added_tags.insert([b'd', b'a'], (info.ambiguous as u8).to_string());
                                }

                                // IDF-router confidence, when it made the call:
                                // ib = best summed IDF weight, im = top-2 margin ratio,
                                // ik = informative k-mers, ia = 1 if ambiguous else 0.
                                if let Some(info) = &idf_info {
                                    added_tags.insert([b'i', b'b'], format!("{:.3}", info.best_score));
                                    added_tags.insert([b'i', b'm'], format!("{:.3}", info.margin_ratio));
                                    added_tags.insert([b'i', b'k'], info.informative_kmers.to_string());
                                    added_tags.insert([b'i', b'a'], (info.ambiguous as u8).to_string());
                                }

                                // POA-classifier confidence, when it made the call:
                                // pb = branch-votes for the best ref, pm = top-2 margin,
                                // pi = discriminating columns matched, pa = 1 if ambiguous.
                                if let Some(info) = &poa_info {
                                    added_tags.insert([b'p', b'b'], info.best_score.to_string());
                                    added_tags.insert([b'p', b'm'], info.margin.to_string());
                                    added_tags.insert([b'p', b'i'], info.informative_columns.to_string());
                                    added_tags.insert([b'p', b'a'], (info.ambiguous as u8).to_string());
                                }

                                added_tags
                                    .insert([b'a', b'r'], read.aligned_read.read_name.clone());
                                added_tags.insert(
                                    [b'r', b'm'],
                                    get_reference_alignment_rate(
                                        &read.aligned_read.reference_aligned,
                                        &read.aligned_read.read_aligned,
                                    )
                                    .to_string(),
                                );
                                added_tags
                                    .insert([b'a', b's'], read.aligned_read.score.to_string());

                                let event_calls = crate::events::call_read_event_details(
                                    &read.aligned_read.reference_aligned,
                                    &read.aligned_read.read_aligned,
                                    structure,
                                );
                                if !event_calls.events.is_empty() {
                                    added_tags.insert(
                                        crate::consensus::consensus_builders::EVENT_TAG,
                                        event_calls.events,
                                    );
                                }
                                if let Some(prime_edits) = event_calls.prime_edits {
                                    added_tags.insert(
                                        crate::consensus::consensus_builders::PRIME_EDIT_TAG,
                                        prime_edits,
                                    );
                                }

                                // Item B(ii): called inversions as `iv:Z:<len>V+<pos>` (0-based,
                                // ungapped reference coordinates; multiple joined by `&`).
                                if !inversion_events.is_empty() {
                                    added_tags.insert([b'i', b'v'], inversion_events.join("&"));
                                }

                                let output = Arc::clone(&output);
                                let arc_writer = output.clone();
                                let mut arc_writer = arc_writer
                                    .lock()
                                    .expect("Unable to access multi-threaded writer");
                                arc_writer
                                    .write_read(&read, &added_tags)
                                    .expect("Unable to write a read to the arc writer (LOC1)");

                                let aligned_count = counters
                                    .aligned_reads
                                    .fetch_add(1, AtomicOrdering::Relaxed)
                                    + 1;
                                if ambiguous_reference_call {
                                    counters
                                        .ambiguous_reference_calls
                                        .fetch_add(1, AtomicOrdering::Relaxed);
                                }
                                *counters
                                    .aligned_by_reference
                                    .lock()
                                    .expect("Unable to lock per-reference alignment counts")
                                    .entry(read.aligned_read.reference_name.clone())
                                    .or_insert(0) += 1;
                                if aligned_count % 1000000 == 0 {
                                    info!(
                                        "Aligned {} reads in {:?}",
                                        aligned_count,
                                        start.elapsed()
                                    );
                                }
                            }
                        }
                    }
                }
            } else {
                counters.too_long.fetch_add(1, AtomicOrdering::Relaxed);
                warn!(
                    "Dropped read {} because its length {} meets or exceeds the maximum {}",
                    String::from_utf8(xx.name().clone()).unwrap(),
                    xx.seq().len(),
                    max_read_size
                );
            }
        });
    });

    let output = Arc::clone(&output);
    let arc_writer = output.clone();
    let mut arc_writer = arc_writer
        .lock()
        .expect("Unable to access multi-threaded writer");
    arc_writer.close().unwrap();

    let aligned_by_reference = counters
        .aligned_by_reference
        .lock()
        .expect("Unable to lock per-reference alignment counts")
        .clone();
    AlignmentRunStats {
        input_reads: counters.input_reads.load(AtomicOrdering::Relaxed),
        aligned_reads: counters.aligned_reads.load(AtomicOrdering::Relaxed),
        too_short: counters.too_short.load(AtomicOrdering::Relaxed),
        too_long: counters.too_long.load(AtomicOrdering::Relaxed),
        failed: counters.failed.load(AtomicOrdering::Relaxed),
        ambiguous_reference_calls: counters
            .ambiguous_reference_calls
            .load(AtomicOrdering::Relaxed),
        aligned_by_reference,
        elapsed_seconds: start.elapsed().as_secs_f64(),
    }
}

#[allow(dead_code)]
pub fn align_two_strings(
    read1_name: &String,
    read1_seq: &[u8],
    sequence_2_seq: &[u8],
    scoring_function: &AffineScoring,
    local: bool,
    ref_name: &String,
    reference_manager: Option<&ReferenceManager>,
) -> AlignmentResult {
    let mut alignment_mat = create_scoring_record_3d(
        read1_seq.len() + 1,
        sequence_2_seq.len() + 1,
        AlignmentType::Affine,
        local,
    );

    match (reference_manager, ref_name) {
        (Some(x), y) => {
            let ref_id = x.reference_name_to_ref.get(y.as_bytes()).unwrap();
            let shared_segments = &x.references.get(ref_id).unwrap().suffix_table;
            let ref_name =
                String::from_utf8(x.references.get(ref_id).unwrap().name.clone()).unwrap();

            let ref_seq = read1_seq;
            let read_seq = sequence_2_seq;

            let shared_segs =
                find_greedy_non_overlapping_segments(&ref_seq, &read_seq, shared_segments);

            align_string_with_anchors(
                read1_name,
                &ref_name,
                read1_seq,
                sequence_2_seq,
                &shared_segs,
                None,
                scoring_function,
                &mut alignment_mat,
            )
        }

        _ => {
            perform_affine_alignment(
                &mut alignment_mat,
                read1_seq,
                sequence_2_seq,
                scoring_function,
            );

            perform_3d_global_traceback(
                &mut alignment_mat,
                None,
                read1_seq,
                sequence_2_seq,
                read1_name,
                ref_name,
                None,
                None,
            )
        }
    }
}

/// Aligns two sequences based on a provided matrix, using either specialized or affine alignment methods.
///
/// This function aligns two DNA or RNA sequences (`read1_seq` and `read2_seq`) and utilizes an
/// affine scoring function. It can optionally use a reference manager and a reference name to
/// perform specialized alignment using shared segments. If these are not provided, it defaults
/// to affine alignment.
///
/// # Arguments
///
/// * `read1_seq` - A reference to the vector containing the first sequence (as `FastaBase` elements) to be aligned.
/// * `read2_seq` - A reference to the vector containing the second sequence to be aligned.
/// * `scoring_function` - A reference to an object implementing the `AffineScoringFunction` trait for scoring alignments.
/// * `ref_name` - An optional reference to a vector of bytes representing the name of the reference sequence.
/// * `reference_manager` - An optional reference to a `ReferenceManager` object managing reference sequences.
/// * `alignment_mat` - A mutable reference to an `Alignment` object, a 3D matrix (`Ix3`) used for storing alignment scores.
///
/// # Returns
///
/// An `AlignmentResult` that may contain details such as the alignment score, the aligned sequences, and other relevant information.
///
/// # Behavior
///
/// The function operates in two modes depending on the availability of `reference_manager` and `ref_name`:
///
/// 1. **With Reference Manager and Reference Name**:
///    - Retrieves a reference ID and shared segments from the reference manager using the provided `ref_name`.
///    - Converts both sequences to `Vec<u8>`.
///    - Finds greedy non-overlapping segments between the sequences using the shared segments.
///    - Performs alignment of `read2_seq` with `read1_seq` using these segments and the scoring function.
///
/// 2. **Without Reference Manager and Reference Name**:
///    - Performs an affine alignment between `read1_seq` and `read2_seq` using the provided scoring function.
///    - Executes a 3D global traceback on the alignment matrix to produce the final alignment result.
///
/// # Panics
///
/// The function may panic if invalid references or indices are encountered, particularly when unwrapping
/// `Option` types or retrieving references from the reference manager.
///
/// # Example
///
/// ```
/// let read1 = vec![/* ... FastaBase elements ... */];
/// let read2 = vec![/* ... FastaBase elements ... */];
/// let scoring_function = /* ... implementation of AffineScoringFunction ... */;
/// let ref_name = Some(&vec![/* reference name as Vec<u8> */]);
/// let reference_manager = Some(&/* ReferenceManager instance */);
/// let mut alignment_matrix = Alignment::new(/* ... dimensions ... */);
///
/// let alignment_result = align_two_strings_passed_matrix(
///     &read1,
///     &read2,
///     &scoring_function,
///     ref_name,
///     reference_manager,
///     &mut alignment_matrix
/// );
///
/// // Use `alignment_result` here
/// ```
pub fn align_two_strings_passed_matrix(
    read1_name: &String,
    read2_name: &String,
    read1_seq: &[u8],
    read2_seq: &[u8],
    qual_sequence: Option<Vec<u8>>,
    scoring_function: &AffineScoring,
    alignment_mat: &mut Alignment<Ix3>,
    max_indel: &usize,
) -> AlignmentResult {
    /*match (reference_manager, ref_name) {
    (Some(x), Some(y)) => {
        //let ref_id = x.reference_name_to_ref.get(y).unwrap();
        //let shared_segments = &x.references.get(ref_id).unwrap().suffix_table;

        //let ref_seq = FastaBase::to_vec_u8(read1_seq);
        //let read_seq = FastaBase::to_vec_u8(read2_seq);

        // TODO fix this
        //let shared_segs = find_greedy_non_overlapping_segments(
        //    &read_seq,
        //    &ref_seq,
        //    shared_segments);

        perform_affine_alignment(
            alignment_mat,
            read1_seq,
            read2_seq,
            scoring_function);

        perform_3d_global_traceback(
            alignment_mat,
            None,
            read1_seq,
            read2_seq,
            None)

        /*align_string_with_anchors(read2_seq,
                                  read1_seq,
                                  &shared_segs,
                                  None,
                                  scoring_function,
                                  alignment_mat)*/
    }

    _ => {*/
    perform_affine_alignment_bandwidth(
        alignment_mat,
        read1_seq,
        read2_seq,
        scoring_function,
        &max_indel,
    );

    perform_3d_global_traceback(
        alignment_mat,
        None,
        read1_seq,
        read2_seq,
        read1_name,
        read2_name,
        qual_sequence,
        None,
    )
    //}
    //}
}

/// Confidence attached to a reference call made by the discriminating-position
/// classifier (only present when that path selected the reference).
#[derive(Clone, Debug)]
pub struct DiscriminatingInfo {
    /// Top-2 score gap at the discriminating positions.
    pub margin: usize,
    /// Discriminating positions the read covered with a real base.
    pub informative_positions: usize,
    /// True when the margin was below the configured minimum (a near-tie call).
    pub ambiguous: bool,
}

/// Confidence attached to a reference call made by the IDF-weighted k-mer router.
#[derive(Clone, Debug)]
pub struct IdfInfo {
    /// Summed IDF weight for the chosen reference.
    pub best_score: f64,
    /// Top-2 relative margin `(best - second) / best`.
    pub margin_ratio: f64,
    /// Distinct read k-mers that hit a discriminating panel k-mer.
    pub informative_kmers: usize,
    /// True when the margin was below the configured floor (a near-tie call).
    pub ambiguous: bool,
}

/// Confidence attached to a reference call made by the POA-graph classifier.
#[derive(Clone, Debug)]
pub struct PoaInfo {
    /// Branch-votes for the chosen reference.
    pub best_score: usize,
    /// Top-2 branch-vote margin.
    pub margin: usize,
    /// Discriminating columns the read matched.
    pub informative_columns: usize,
    /// True when the margin was below the configured floor (a near-tie call).
    pub ambiguous: bool,
}

#[allow(dead_code)]
#[derive(Clone)]
pub struct AlignmentWithRef {
    alignment: Option<AlignmentResult>,
    ref_name: Vec<u8>,
    ref_sequence: Vec<u8>,
    /// Set when the reference was chosen by the discriminating-position classifier.
    classification: Option<DiscriminatingInfo>,
    /// Set when the reference was routed by the IDF-weighted k-mer index.
    idf_info: Option<IdfInfo>,
    /// Set when the reference was chosen by the POA-graph classifier.
    poa_info: Option<PoaInfo>,
}

fn alignment_score(candidate: &AlignmentWithRef) -> f64 {
    candidate
        .alignment
        .as_ref()
        .map(|alignment| alignment.score)
        .unwrap_or(f64::NEG_INFINITY)
}

fn select_best_alignment(
    forward: Option<AlignmentWithRef>,
    reverse_complemented: Option<AlignmentWithRef>,
) -> Option<AlignmentWithRef> {
    match (forward, reverse_complemented) {
        (Some(forward), Some(reverse_complemented)) => {
            if alignment_score(&reverse_complemented) > alignment_score(&forward) {
                Some(reverse_complemented)
            } else {
                Some(forward)
            }
        }
        (Some(forward), None) => Some(forward),
        (None, Some(reverse_complemented)) => Some(reverse_complemented),
        (None, None) => None,
    }
}

fn search_multiple_references(
    read_name: &String,
    read: &Vec<u8>,
    qual_sequence: Option<Vec<u8>>,
    rm: &ReferenceManager,
    fast_lookup: bool,
    alignment_mat: &mut Alignment<Ix3>,
    my_aff_score: &AffineScoring,
    classifier: Option<&DiscriminatingClassifier>,
    idf: Option<&IdfIndex>,
    poa: Option<&PoaGraph>,
) -> Option<AlignmentWithRef> {
    // POA-graph classifier (for near-identical panels): align the read to the
    // reference graph and vote, at each discriminating column, for the references
    // whose branch it took; then align to the chosen reference and attach the
    // confidence. Takes precedence over the other routers when enabled.
    if let Some(graph) = poa {
        let cls = graph.classify_read(read);
        let mut subset = HashSet::new();
        subset.insert(cls.best.clone());
        return exhaustive_alignment_search(
            read_name,
            read,
            qual_sequence,
            rm,
            alignment_mat,
            my_aff_score,
            Some(subset),
        )
        .map(|mut awr| {
            awr.poa_info = Some(PoaInfo {
                best_score: cls.best_score,
                margin: cls.margin,
                informative_columns: cls.informative_columns,
                ambiguous: cls.ambiguous,
            });
            awr
        });
    }

    // IDF-weighted k-mer router: pick the reference by summed IDF weight (the
    // shared backbone contributes 0). When the IDF call is ambiguous and a
    // discriminating classifier is available, defer the tie-break to it. Then
    // align the read to the chosen reference (or fall back to a full search when
    // IDF found no discriminating signal) and attach the confidence.
    if let Some(index) = idf {
        let sc = index.score(read);
        let mut refined_disc: Option<DiscriminatingInfo> = None;
        let chosen: Option<Vec<u8>> = match &sc.best {
            None => None, // no discriminating k-mers matched -> full fallback below
            Some(best) => {
                if sc.ambiguous {
                    if let Some(clf) = classifier {
                        let dc = clf.classify_read(read);
                        refined_disc = Some(DiscriminatingInfo {
                            margin: dc.margin,
                            informative_positions: dc.informative_positions,
                            ambiguous: dc.ambiguous,
                        });
                        Some(dc.best.into_bytes())
                    } else {
                        Some(best.clone())
                    }
                } else {
                    Some(best.clone())
                }
            }
        };
        let subset = chosen.map(|name| {
            let mut set = HashSet::new();
            set.insert(name);
            set
        });
        return exhaustive_alignment_search(
            read_name,
            read,
            qual_sequence,
            rm,
            alignment_mat,
            my_aff_score,
            subset,
        )
        .map(|mut awr| {
            awr.idf_info = Some(IdfInfo {
                best_score: sc.best_score,
                margin_ratio: sc.margin_ratio,
                informative_kmers: sc.informative_kmers,
                ambiguous: sc.ambiguous,
            });
            awr.classification = refined_disc.clone();
            awr
        });
    }

    // Discriminating-position path (for near-identical panels): let the classifier
    // pick the reference by comparing only the columns where the panel differs,
    // then align the read to that single reference and attach the confidence.
    if let Some(clf) = classifier {
        let cls = clf.classify_read(read);
        let mut subset = HashSet::new();
        subset.insert(cls.best.clone().into_bytes());
        return exhaustive_alignment_search(
            read_name,
            read,
            qual_sequence,
            rm,
            alignment_mat,
            my_aff_score,
            Some(subset),
        )
        .map(|mut awr| {
            awr.classification = Some(DiscriminatingInfo {
                margin: cls.margin,
                informative_positions: cls.informative_positions,
                ambiguous: cls.ambiguous,
            });
            awr
        });
    }
    if fast_lookup {
        quick_alignment_search(
            read_name,
            read,
            qual_sequence,
            rm,
            alignment_mat,
            my_aff_score,
            &0.90,
        )
    } else {
        exhaustive_alignment_search(
            read_name,
            read,
            qual_sequence,
            rm,
            alignment_mat,
            my_aff_score,
            None,
        )
    }
}

/// Aligns two DNA or RNA sequences using affine alignment or a specialized alignment with anchors.
///
/// This function aligns two sequences represented by `Vec<FastaBase>` using either affine alignment
/// or a specialized alignment with anchors, depending on the availability of a reference manager
/// and a reference name. It employs a scoring function to evaluate the alignments.
///
/// # Arguments
///
/// * `read1_seq` - A reference to the first sequence to be aligned, each element being of type `FastaBase`.
/// * `read2_seq` - A reference to the second sequence to be aligned.
/// * `scoring_function` - A reference to an object implementing the `AffineScoringFunction` trait,
///    used for scoring alignments.
/// * `ref_name` - An optional reference to a byte vector representing the name of the reference sequence.
/// * `reference_manager` - An optional reference to a `ReferenceManager` object, which manages reference sequences.
/// * `alignment_mat` - A mutable reference to an `Alignment` object (3D matrix) used for storing alignment scores.
///
/// # Returns
///
/// An `AlignmentResult` containing the details of the alignment, which may include the alignment score,
/// aligned sequences, and other relevant information.
///
/// # Behavior
///
/// If both `reference_manager` and `ref_name` are provided, the function:
/// - Retrieves a reference ID and shared segments from the reference manager.
/// - Converts input sequences to `Vec<u8>`.
/// - Finds greedy non-overlapping segments between the read and reference sequences.
/// - Aligns `read2_seq` with `read1_seq` using these segments, the scoring function, and the alignment matrix.
///
/// If either `reference_manager` or `ref_name` is `None`, the function defaults to:
/// - Performing an affine alignment between `read1_seq` and `read2_seq`.
/// - Conducting a 3D global traceback on the alignment matrix for the final alignment result.
///
/// # Panics
///
/// The function may panic if invalid references or indices are encountered, particularly when unwrapping
/// `Option` types or retrieving references from the reference manager. Proper exception handling should be
/// considered in the calling code.
///
/// # Example
///
/// ```
/// let read1 = vec![/* ... FastaBase elements ... */];
/// let read2 = vec![/* ... FastaBase elements ... */];
/// let scoring_function = /* ... implementation of AffineScoringFunction ... */;
/// let ref_name = Some(&vec![/* reference name as Vec<u8> */]);
/// let reference_manager = Some(&/* ReferenceManager instance */);
/// let mut alignment_matrix = Alignment::new(/* ... dimensions ... */);
///
/// let result = align_two_strings_passed_matrix(
///     &read1,
///     &read2,
///     &scoring_function,
///     ref_name,
///     reference_manager,
///     &mut alignment_matrix
/// );
///
/// // Use `result` here
/// ```
///
pub fn align_to_reference_choices(
    read_name: &String,
    read: &Vec<u8>,
    qual_sequence: Option<Vec<u8>>,
    rm: &ReferenceManager,
    fast_lookup: &bool,
    read_structure: &SequenceLayout,
    alignment_mat: &mut Alignment<Ix3>,
    my_aff_score: &AffineScoring,
    my_score: &InversionScoring,
    use_inversions: &bool,
    convex_score: &ConvexScoring,
    use_convex: &bool,
    _max_reference_multiplier: f64,
    min_read_length: usize,
    _max_indel: &usize,
    classifier: Option<&DiscriminatingClassifier>,
    idf: Option<&IdfIndex>,
    poa: Option<&PoaGraph>,
) -> Option<AlignmentWithRef> {
    if read.len() < min_read_length {
        debug!(
            "Skipping read {} because its length {} is below the minimum {}",
            read_name,
            read.len(),
            min_read_length
        );
        return None;
    }

    match rm.references.len() {
        0 => {
            // TODO: we should track this and provide a final summary
            warn!(
                "Unable to align read {} as it has no candidate references",
                u8s(read)
            );
            None
        }
        1 => {
            // exactly one reference: take the sole entry (its map key is not guaranteed to be 0).
            let ref_base = &rm.references.values().next().unwrap();
            let ref_name = String::from_utf8(ref_base.name.clone()).unwrap();
            let (forward_oriented_seq, oriented_quals) = if !read_structure.known_strand {
                let orientation_search =
                    orient_by_longest_segment(&read, &ref_base.sequence, &ref_base.suffix_table);
                let forward_seed_score = orientation_search
                    .1
                    .alignment_segments
                    .iter()
                    .map(|segment| segment.length)
                    .sum::<usize>();
                let reverse_seed_score = orientation_search
                    .2
                    .alignment_segments
                    .iter()
                    .map(|segment| segment.length)
                    .sum::<usize>();
                let orientation = match forward_seed_score.cmp(&reverse_seed_score) {
                    std::cmp::Ordering::Greater => true,
                    std::cmp::Ordering::Less => false,
                    std::cmp::Ordering::Equal => {
                        let reverse_complemented = reverse_complement(&read);
                        let forward_score = rust_bio_alignment(
                            &ref_base.sequence,
                            read,
                            &4,
                            &10,
                            &1,
                        )
                        .1;
                        let reverse_score = rust_bio_alignment(
                            &ref_base.sequence,
                            &reverse_complemented,
                            &4,
                            &10,
                            &1,
                        )
                        .1;
                        forward_score >= reverse_score
                    }
                };
                if orientation {
                    (read.clone(), qual_sequence)
                } else {
                    // reverse-complementing the read reverses base order, so the quals must be
                    // reversed too to stay paired with the aligned bases.
                    let reversed_quals = qual_sequence.map(|mut q| { q.reverse(); q });
                    (reverse_complement(&read), reversed_quals)
                }
            } else {
                (read.clone(), qual_sequence)
            };

            let result = if *use_inversions {
                // Item A: route the final single-reference alignment through the
                // inversion-aware aligner. Its CIGAR may carry InversionOpen/Close
                // markers; qualities are threaded through and reordered within
                // inverted blocks, and the markers are flattened into an `iv` event
                // (item B(ii)) by the writer.
                inversion_alignment(
                    &ref_base.sequence,
                    &forward_oriented_seq,
                    &ref_name,
                    read_name,
                    my_score,
                    my_aff_score,
                    false,
                    oriented_quals,
                )
            } else if *use_convex {
                // `--aligner convex`: global alignment with a convex (logarithmic)
                // gap penalty, so one long indel is preferred over several short ones.
                convex_alignment(
                    &ref_base.sequence,
                    &forward_oriented_seq,
                    &ref_name,
                    read_name,
                    convex_score,
                    oriented_quals,
                )
            } else {
                let alignment = rust_bio_alignment(&ref_base.sequence, &forward_oriented_seq, &4, &10, &1).0;
                //println!("{}", alignment);

                let alignment = cigar_to_alignment(
                    &ref_base.sequence,
                    &forward_oriented_seq,
                    &alignment,
                );

                AlignmentResult {
                    reference_name: ref_name.clone(),
                    read_name: read_name.clone(),
                    reference_aligned: alignment.0,
                    read_aligned: alignment.1,
                    read_quals: oriented_quals,
                    cigar_string: alignment.2,
                    path: vec!(),
                    score: 0.0,
                    reference_start: 0,
                    read_start: 0,
                    bounding_box: None,
                }
            };


            /*
            let aln = align_two_strings_passed_matrix(
                &ref_name,
                read_name,
                &ref_base.sequence,
                &forward_oriented_seq,
                qual_sequence,
                my_aff_score,
                alignment_mat,
                &100,
            );
*/
            Some(AlignmentWithRef {
                alignment: Some(result),
                ref_name: ref_base.name.clone(),
                ref_sequence: ref_base.sequence.clone(),
                classification: None,
                idf_info: None,
                poa_info: None,
            })
        }
        x if x > 1 => {
            let base = if read_structure.known_strand {
                search_multiple_references(
                    read_name,
                    read,
                    qual_sequence,
                    rm,
                    *fast_lookup,
                    alignment_mat,
                    my_aff_score,
                    classifier,
                    idf,
                    poa,
                )
            } else {
                let forward = search_multiple_references(
                    read_name,
                    read,
                    qual_sequence.clone(),
                    rm,
                    *fast_lookup,
                    alignment_mat,
                    my_aff_score,
                    classifier,
                    idf,
                    poa,
                );

                let reverse_complemented_read = reverse_complement(read);
                let reverse_complemented_qualities = qual_sequence.map(|mut qualities| {
                    qualities.reverse();
                    qualities
                });
                let reverse_complemented = search_multiple_references(
                    read_name,
                    &reverse_complemented_read,
                    reverse_complemented_qualities,
                    rm,
                    *fast_lookup,
                    alignment_mat,
                    my_aff_score,
                    classifier,
                    idf,
                    poa,
                );

                select_best_alignment(forward, reverse_complemented)
            };

            // Item A (multi-reference): once a router has chosen the reference, re-align
            // the read to that reference with the inversion-aware aligner and swap the
            // result in, preserving the router confidence tags. The read is taken in the
            // orientation the router settled on (read_aligned with gaps removed).
            if *use_inversions {
                base.map(|mut awr| {
                    if let Some(existing) = awr.alignment.take() {
                        let oriented_read: Vec<u8> = existing
                            .read_aligned
                            .iter()
                            .filter(|b| **b != FASTA_UNSET)
                            .cloned()
                            .collect();
                        let ref_name = String::from_utf8(awr.ref_name.clone()).unwrap();
                        let inv = inversion_alignment(
                            &awr.ref_sequence,
                            &oriented_read,
                            &ref_name,
                            read_name,
                            my_score,
                            my_aff_score,
                            false,
                            existing.read_quals.clone(),
                        );
                        awr.alignment = Some(inv);
                    }
                    awr
                })
            } else if *use_convex {
                // Multi-reference convex: re-align the router-chosen reference with the
                // convex-gap aligner (opt-in via --aligner convex).
                base.map(|mut awr| {
                    if let Some(existing) = awr.alignment.take() {
                        let oriented_read: Vec<u8> = existing
                            .read_aligned
                            .iter()
                            .filter(|b| **b != FASTA_UNSET)
                            .cloned()
                            .collect();
                        let ref_name = String::from_utf8(awr.ref_name.clone()).unwrap();
                        let cv = convex_alignment(
                            &awr.ref_sequence,
                            &oriented_read,
                            &ref_name,
                            read_name,
                            convex_score,
                            existing.read_quals.clone(),
                        );
                        awr.alignment = Some(cv);
                    }
                    awr
                })
            } else {
                base
            }
        }
        x => {
            panic!("we dont know what to do with a reference count of {}", x)
        }
    }
}

/// Performs a quick alignment search for a given read against a reference manager.
///
/// This function takes a read sequence and aligns it to the most likely reference sequence
/// in a reference manager using k-mers. It supports both affine scoring and inversion scoring,
/// and can optionally include inversions in the alignment process.
///
/// # Arguments
///
/// * `read` - A reference to the sequence to be aligned, represented as `Vec<FastaBase>`.
/// * `rm` - A reference to the `ReferenceManager` containing reference sequences and k-mer information.
/// * `read_structure` - A reference to the `SequenceLayoutDesign` that describes the layout of the read sequence.
/// * `alignment_mat` - A mutable reference to an `Alignment` object (3D matrix) used for storing alignment scores.
/// * `my_aff_score` - A reference to an `AffineScoring` object for scoring alignments.
/// * `my_score` - A reference to an `InversionScoring` object for scoring inversions, if used.
/// * `use_inversions` - A boolean reference indicating whether inversions should be considered in the alignment.
/// * `max_reference_multiplier` - A floating-point value representing a multiplier to determine the maximum reference size.
/// * `min_read_length` - The minimum length of the read sequence for performing the alignment.
///
/// # Returns
///
/// Returns `Option<(Option<AlignmentResult>, Vec<u8>, Vec<u8>)>`. The tuple inside the `Option` contains:
/// - An `Option<AlignmentResult>` which is `None` if no alignment is found, or contains the alignment result.
/// - A `Vec<u8>` representing the sequence of the matched reference.
/// - A `Vec<u8>` representing the name of the matched reference.
///
/// # Behavior
///
/// The function converts the read sequence into a `Vec<u8>` and generates k-mers. It then finds
/// the reference sequence in the reference manager with the highest count of matching k-mers.
/// If such a reference is found, it performs an alignment using `align_two_strings_passed_matrix`
/// and returns the alignment result along with the sequence and name of the matched reference.
/// If no matching reference is found, it returns `None`.
///
/// # Example
///
/// ```
/// let read = vec![/* ... FastaBase elements ... */];
/// let reference_manager = /* ... ReferenceManager instance ... */;
/// let read_structure = /* ... SequenceLayoutDesign instance ... */;
/// let mut alignment_matrix = Alignment::new(/* ... dimensions ... */);
/// let affine_scoring = /* ... AffineScoring instance ... */;
/// let inversion_scoring = /* ... InversionScoring instance ... */;
/// let use_inversions = &true;
/// let max_reference_multiplier = 1.5;
/// let min_read_length = 100;
///
/// let result = quick_alignment_search(
///     &read,
///     &reference_manager,
///     &read_structure,
///     &mut alignment_matrix,
///     &affine_scoring,
///     &inversion_scoring,
///     use_inversions,
///     max_reference_multiplier,
///     min_read_length
/// );
///
/// // Process `result` here
/// ```
const MIN_DISTINCT_UNIQUE_KMER_HITS: usize = 3;

fn has_confident_kmer_support(
    winner_hits: usize,
    informative_hits: usize,
    match_threshold: f64,
) -> bool {
    informative_hits > 0
        && winner_hits >= MIN_DISTINCT_UNIQUE_KMER_HITS
        && winner_hits as f64 / informative_hits as f64 > match_threshold
}

fn quick_alignment_search(
    read_name: &String,
    read: &Vec<u8>,
    qual_sequence: Option<Vec<u8>>,
    rm: &ReferenceManager,
    alignment_mat: &mut Alignment<Ix3>,
    my_aff_score: &AffineScoring,
    match_threshold: &f64,
) -> Option<AlignmentWithRef> {
    let read_u8 = read;
    let read_kmers = ReferenceManager::sequence_to_kmers(&read_u8, &rm.kmer_size, &rm.kmer_skip);

    let mut seen_kmers = HashSet::new();
    let reference_votes = read_kmers
        .iter()
        .filter(|(kmer, _count)| seen_kmers.insert(kmer.as_slice()))
        .filter_map(|(kmer, _count)| rm.unique_kmers.kmer_to_reference.get(kmer))
        .counts();

    let informative_hits = reference_votes.values().sum::<usize>();
    let max_ref = reference_votes.iter().max_by_key(|(_reference, hits)| *hits);

    match max_ref {
        None => {
            info!("No reference found; moving to exhaustive_alignment_search");
            exhaustive_alignment_search(
                read_name,
                read,
                qual_sequence,
                rm,
                alignment_mat,
                my_aff_score,
                None,
            )
        }
        Some((reference, winner_hits)) => {
            if has_confident_kmer_support(*winner_hits, informative_hits, *match_threshold) {
                let ref_name = String::from_utf8(reference.name.clone()).unwrap();
                Some(AlignmentWithRef {
                    alignment: Some(align_two_strings_passed_matrix(
                        &ref_name,
                        read_name,
                        &reference.sequence,
                        read,
                        qual_sequence,
                        my_aff_score,
                        alignment_mat,
                        &read.len(),
                    )),
                    ref_name: reference.name.clone(),
                    ref_sequence: reference.sequence.clone(),
                    classification: None,
                    idf_info: None,
                    poa_info: None,
                })
            } else {
                exhaustive_alignment_search(
                    read_name,
                    read,
                    qual_sequence,
                    rm,
                    alignment_mat,
                    my_aff_score,
                    None,
                )
            }
        }
    }
}

fn exhaustive_alignment_search(
    read_name: &String,
    read: &Vec<u8>,
    qual_sequence: Option<Vec<u8>>,
    rm: &ReferenceManager,
    alignment_mat: &mut Alignment<Ix3>,
    my_aff_score: &AffineScoring,
    reference_subset: Option<HashSet<Vec<u8>>>,
) -> Option<AlignmentWithRef> {
    let references = &rm.references;

    let ranked_alignments = references
        .iter()
        .map(|reference| {
            if reference_subset.is_none()
                || reference_subset
                    .as_ref()
                    .unwrap()
                    .contains(&reference.1.name)
            {
                let qual = qual_sequence.clone();
                let lt = align_two_strings_passed_matrix(
                    &String::from_utf8(reference.1.name.clone()).unwrap(),
                    read_name,
                    &reference.1.sequence,
                    read,
                    qual,
                    my_aff_score,
                    alignment_mat,
                    &read.len(),
                );

                Some((lt, reference.1.sequence.clone(), reference.1.name.clone()))
            } else {
                None
            }
        })
        .filter(|x| x.is_some())
        .map(|c| c.unwrap());

    let ranked_alignments = ranked_alignments.into_iter().enumerate().max_by(|al, al2| {
        let score1 = al.1 .0.score; // / al.1.0.reference_aligned.len() as f64;
        let score2 = al2.1 .0.score; // / al2.1.0.reference_aligned.len() as f64;
        score1.partial_cmp(&score2).unwrap()
    });

    match ranked_alignments.iter().next() {
        None => None,
        Some((_x, y)) => {
            //Some((Some(y.0.clone()), , y.2.clone()))
            //println!("---- {} {}",String::from_utf8(y.1.clone()).unwrap(),String::from_utf8(y.2.clone()).unwrap());
            Some(AlignmentWithRef {
                alignment: Some(y.0.clone()),
                ref_name: y.2.clone(),
                ref_sequence: y.1.clone(),
                classification: None,
                idf_info: None,
                poa_info: None,
            })
        }
    }
}

#[allow(dead_code)]
fn cigar_to_alignment(reference: &Vec<u8>, read: &Vec<u8>, cigar: &Vec<AlignmentOperation>) -> (Vec<u8>, Vec<u8>, Vec<AlignmentTag>) {
    let mut alignment_string1 = Vec::new();
    let mut alignment_string2 = Vec::new();
    let mut cigar_vec = Vec::new();

    let mut seq1_index = 0;
    let mut seq2_index = 0;

    for c in cigar {
        match c {
            Match => {
                alignment_string1.push(reference.get(seq1_index).unwrap().clone());
                alignment_string2.push(read.get(seq2_index).unwrap().clone());
                cigar_vec.push(AlignmentTag::MatchMismatch(1));
                seq1_index += 1;
                seq2_index += 1;
            }
            Del => {
                alignment_string1.push(reference.get(seq1_index).unwrap().clone());
                alignment_string2.push(FASTA_UNSET);
                cigar_vec.push(AlignmentTag::Del(1));
                seq1_index += 1;
            }
           Ins => {
                alignment_string1.push(FASTA_UNSET);
                alignment_string2.push(read.get(seq2_index).unwrap().clone());
                cigar_vec.push(AlignmentTag::Ins(1));
                seq2_index += 1;
            },
            Subst => {
                alignment_string1.push(reference.get(seq1_index).unwrap().clone());
                alignment_string2.push(read.get(seq2_index).unwrap().clone());
                cigar_vec.push(AlignmentTag::MatchMismatch(1));
                seq1_index += 1;
                seq2_index += 1;
            }
            _ => {
                panic!("Unknown cigar operation {:?}",c)
            }
        }
    }
    (alignment_string1, alignment_string2, simplify_cigar_string(&cigar_vec))
}

pub fn simplify_cigar_string(cigar_tokens: &Vec<AlignmentTag>) -> Vec<AlignmentTag> {
    let mut new_cigar = Vec::new();

    let mut last_token: Option<AlignmentTag> = None; // zero length, so combining won't affect the final cigar string

    cigar_tokens
        .iter()
        .for_each(|token| match (&last_token, token) {
            (None, _) => last_token = Some(token.clone()),
            (Some(AlignmentTag::InversionOpen), AlignmentTag::InversionOpen) => {
                panic!("Cannot have two inversion open tags in a row");
            }
            (Some(AlignmentTag::InversionClose), AlignmentTag::InversionClose) => {
                panic!("Cannot have two inversion closed tags in a row");
            }
            (
                Some(AlignmentTag::MatchMismatch(last_count)),
                AlignmentTag::MatchMismatch(this_count),
            ) => {
                last_token = Some(AlignmentTag::MatchMismatch(last_count + this_count));
            }
            (Some(AlignmentTag::Del(last_count)), AlignmentTag::Del(this_count)) => {
                last_token = Some(AlignmentTag::Del(last_count + this_count));
            }
            (Some(AlignmentTag::Ins(last_count)), AlignmentTag::Ins(this_count)) => {
                last_token = Some(AlignmentTag::Ins(last_count + this_count));
            }
            (Some(x), y) => {
                new_cigar.push(x.clone());
                last_token = Some(y.clone());
            }
        });

    if let Some(x) = last_token {
        new_cigar.push(x);
    }
    new_cigar
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use crate::alignment::alignment_matrix::{
        create_scoring_record_3d, AlignmentResult, AlignmentTag, AlignmentType,
    };
    use crate::alignment::scoring_functions::{AffineScoring, ConvexScoring, InversionScoring};

    fn test_convex_scoring() -> ConvexScoring {
        ConvexScoring {
            match_score: 10.0,
            mismatch_score: -9.0,
            special_character_score: 9.0,
            gap_open: -10.0,
            gap_extend: -3.0,
        }
    }
    use crate::alignment_functions::{
        align_to_reference_choices, cigar_to_alignment, exhaustive_alignment_search,
        has_confident_kmer_support, quick_alignment_search, simplify_cigar_string,
        reference_histogram_bar, AlignmentRunStats, REFERENCE_HISTOGRAM_WIDTH,
    };
    use crate::read_strategies::sequence_layout::{
        AlignedReadOrientation, ReadPosition, ReferenceRecord, SequenceLayout,
    };
    use crate::reference::fasta_reference::ReferenceManager;
    use crate::extractor::extract_tagged_sequences;
    use crate::utils::read_utils::reverse_complement;
    use bio::alignment::AlignmentOperation;

    #[test]
    fn test_alignment_run_summary_reconciles_outcomes() {
        let stats = AlignmentRunStats {
            input_reads: 10,
            aligned_reads: 6,
            too_short: 1,
            too_long: 2,
            failed: 1,
            ambiguous_reference_calls: 2,
            aligned_by_reference: BTreeMap::from([
                ("reference_a".to_string(), 4),
                ("reference_b".to_string(), 2),
            ]),
            elapsed_seconds: 1.25,
        };

        assert_eq!(
            stats.input_reads,
            stats.aligned_reads + stats.too_short + stats.too_long + stats.failed
        );
        let rendered = stats.summary().render();
        assert!(rendered.contains("| All   | 10          | 6       | 60.00"));
        assert!(rendered.contains("| reference_a | 4             | 66.67"));
        assert!(rendered.contains(&"#".repeat(REFERENCE_HISTOGRAM_WIDTH)));
    }

    #[test]
    fn test_reference_histogram_scales_to_largest_reference() {
        assert_eq!(reference_histogram_bar(10, 10), "#".repeat(40));
        assert_eq!(reference_histogram_bar(5, 10), "#".repeat(20));
        assert_eq!(reference_histogram_bar(1, 100), "#");
        assert_eq!(reference_histogram_bar(0, 100), "");
    }

    fn multi_reference_layout(known_strand: bool) -> SequenceLayout {
        let mut references = BTreeMap::new();
        references.insert(
            "forward_reference".to_string(),
            ReferenceRecord {
                sequence: "TTTTAAAACCCCGGGG".to_string(),
                umi_configurations: BTreeMap::new(),
                targets: vec![],
                target_types: vec![],
                target_locations: None,
                prime_edits: BTreeMap::new(),
            },
        );
        references.insert(
            "reverse_reference".to_string(),
            ReferenceRecord {
                sequence: "ACGTCAGTGGATCCAA".to_string(),
                umi_configurations: BTreeMap::new(),
                targets: vec![],
                target_types: vec![],
                target_locations: None,
                prime_edits: BTreeMap::new(),
            },
        );

        SequenceLayout {
            aligner: None,
            merge: None,
            reads: vec![ReadPosition::Read1 {
                orientation: AlignedReadOrientation::Forward,
            }],
            known_strand,
            references,
        }
    }

    fn align_single_unknown_strand(
        reference: &str,
        read: Vec<u8>,
        qualities: Vec<u8>,
    ) -> AlignmentResult {
        let mut references = BTreeMap::new();
        references.insert(
            "reference".to_string(),
            ReferenceRecord {
                sequence: reference.to_string(),
                umi_configurations: BTreeMap::new(),
                targets: vec![],
                target_types: vec![],
                target_locations: None,
                prime_edits: BTreeMap::new(),
            },
        );
        let layout = SequenceLayout {
            aligner: None,
            merge: None,
            reads: vec![ReadPosition::Read1 {
                orientation: AlignedReadOrientation::Forward,
            }],
            known_strand: false,
            references,
        };
        let reference_manager = ReferenceManager::from_yaml_input(&layout, 8, 4);
        let mut alignment_matrix = create_scoring_record_3d(
            reference_manager.longest_ref + 1,
            read.len() + 1,
            AlignmentType::Affine,
            false,
        );
        let inversion_score = InversionScoring {
            match_score: 9.0,
            mismatch_score: -21.0,
            gap_open: -25.0,
            gap_extend: -1.0,
            inversion_penalty: -40.0,
            min_inversion_length: 20,
        };

        align_to_reference_choices(
            &"read".to_string(),
            &read,
            Some(qualities),
            &reference_manager,
            &false,
            &layout,
            &mut alignment_matrix,
            &AffineScoring::default_dna(),
            &inversion_score,
            &false,
            &test_convex_scoring(),
            &false,
            2.0,
            0,
            &read.len(),
            None,
            None,
            None,
        )
        .unwrap()
        .alignment
        .unwrap()
    }

    /// Item A prototype helper: single forward reference, inversion-aware dispatch on.
    fn align_single_forward_with_inversions(reference: &str, read: Vec<u8>) -> AlignmentResult {
        let mut references = BTreeMap::new();
        references.insert(
            "reference".to_string(),
            ReferenceRecord {
                sequence: reference.to_string(),
                umi_configurations: BTreeMap::new(),
                targets: vec![],
                target_types: vec![],
                target_locations: None,
                prime_edits: BTreeMap::new(),
            },
        );
        let layout = SequenceLayout {
            aligner: None,
            merge: None,
            reads: vec![ReadPosition::Read1 {
                orientation: AlignedReadOrientation::Forward,
            }],
            known_strand: true,
            references,
        };
        let reference_manager = ReferenceManager::from_yaml_input(&layout, 8, 4);
        let mut alignment_matrix = create_scoring_record_3d(
            reference_manager.longest_ref + 1,
            read.len() + 1,
            AlignmentType::Affine,
            false,
        );
        // Lenient inversion scoring (matches the known-good aligner unit tests) so a
        // short synthetic inversion is reliably detected; production tuning lives in
        // align_reads.
        let inversion_score = InversionScoring {
            match_score: 10.0,
            mismatch_score: -11.0,
            gap_open: -15.0,
            gap_extend: -5.0,
            inversion_penalty: -2.0,
            min_inversion_length: 4,
        };
        align_to_reference_choices(
            &"read".to_string(),
            &read,
            None,
            &reference_manager,
            &false,
            &layout,
            &mut alignment_matrix,
            &AffineScoring::default_dna(),
            &inversion_score,
            &true, // use_inversions -> route to the inversion-aware aligner
            &test_convex_scoring(),
            &false,
            2.0,
            0,
            &read.len(),
            None,
            None,
            None,
        )
        .unwrap()
        .alignment
        .unwrap()
    }

    #[test]
    fn test_inversion_aligner_dispatch_detects_inverted_segment() {
        // Item A prototype: with use_inversions = true (the `--aligner inversion`
        // flag), the single-reference path routes to the inversion-aware aligner,
        // which should detect a reverse-complemented internal block and mark it with
        // InversionOpen/InversionClose. The plain affine path never emits those tags.
        // Known-good inversion fixture from the aligner unit tests: the read's
        // internal block is reverse-complemented relative to the reference.
        let reference = "CCAATCTACTACTGCTTGCA";
        let read = b"CCGTAGATTTACTGCTTGCA".to_vec();

        let alignment = align_single_forward_with_inversions(reference, read);

        let has_inversion = alignment.cigar_string.iter().any(|tag| {
            matches!(tag, AlignmentTag::InversionOpen | AlignmentTag::InversionClose)
        });
        assert!(
            has_inversion,
            "inversion-aware dispatch should mark the inverted block; cigar = {:?}",
            alignment.cigar_string
        );
    }

    #[test]
    fn test_single_reference_zero_seed_tie_uses_forward_alignment_score() {
        let reference = b"AACGTA".to_vec();
        let qualities = vec![10, 11, 12, 13, 14, 15];
        let alignment = align_single_unknown_strand(
            std::str::from_utf8(&reference).unwrap(),
            reference.clone(),
            qualities.clone(),
        );

        assert_eq!(alignment.read_aligned, reference);
        assert_eq!(alignment.read_quals, Some(qualities));
    }

    #[test]
    fn test_single_reference_zero_seed_tie_uses_reverse_alignment_score() {
        let reference = b"AACGTA".to_vec();
        let read = reverse_complement(&reference);
        let qualities = vec![10, 11, 12, 13, 14, 15];
        let mut expected_qualities = qualities.clone();
        expected_qualities.reverse();
        let alignment = align_single_unknown_strand(
            std::str::from_utf8(&reference).unwrap(),
            read,
            qualities,
        );

        assert_eq!(alignment.read_aligned, reference);
        assert_eq!(alignment.read_quals, Some(expected_qualities));
    }

    #[test]
    fn test_single_reference_full_score_tie_keeps_forward_quality_order() {
        let qualities = vec![10, 11, 12, 13];
        let alignment = align_single_unknown_strand(
            "ACGT",
            b"ACGT".to_vec(),
            qualities.clone(),
        );

        assert_eq!(alignment.read_quals, Some(qualities));
    }

    #[test]
    fn test_single_reference_symbolic_umi_extracts_bases_without_gaps() {
        let reference = "ACGT0000TTAG";
        let read = b"ACGTGCAATTAG".to_vec();
        let alignment = align_single_unknown_strand(
            reference,
            read.clone(),
            vec![30; read.len()],
        );

        assert_eq!(alignment.reference_aligned, reference.as_bytes());
        assert_eq!(alignment.read_aligned, read);
        assert_eq!(
            extract_tagged_sequences(&alignment.read_aligned, &alignment.reference_aligned)
                .get(&b'0')
                .map(String::as_str),
            Some("GCAA")
        );
    }

    #[test]
    fn test_alignment_rejects_reads_below_minimum_length() {
        let layout = multi_reference_layout(true);
        let reference_manager = ReferenceManager::from_yaml_input(&layout, 8, 4);
        let read = b"ACGT".to_vec();
        let mut alignment_matrix = create_scoring_record_3d(
            reference_manager.longest_ref + 1,
            read.len() + 1,
            AlignmentType::Affine,
            false,
        );

        let result = align_to_reference_choices(
            &"short_read".to_string(),
            &read,
            None,
            &reference_manager,
            &false,
            &layout,
            &mut alignment_matrix,
            &AffineScoring::default_dna(),
            &InversionScoring::default(),
            &false,
            &test_convex_scoring(),
            &false,
            2.0,
            read.len() + 1,
            &read.len(),
            None,
            None,
            None,
        );

        assert!(result.is_none());
    }

    #[test]
    fn test_multi_reference_unknown_strand_checks_reverse_complement() {
        let layout = multi_reference_layout(false);
        let reference = layout
            .references
            .get("reverse_reference")
            .unwrap()
            .sequence
            .as_bytes()
            .to_vec();
        let read = reverse_complement(&reference);
        let qualities = (0..read.len() as u8).collect::<Vec<u8>>();
        let mut expected_qualities = qualities.clone();
        expected_qualities.reverse();
        let reference_manager = ReferenceManager::from_yaml_input(&layout, 8, 4);
        let affine_score = AffineScoring::default_dna();
        let inversion_score = InversionScoring {
            match_score: 9.0,
            mismatch_score: -21.0,
            gap_open: -25.0,
            gap_extend: -1.0,
            inversion_penalty: -40.0,
            min_inversion_length: 20,
        };

        for fast_lookup in [true, false] {
            let mut alignment_matrix = create_scoring_record_3d(
                reference_manager.longest_ref + 1,
                read.len() + 1,
                AlignmentType::Affine,
                false,
            );
            let result = align_to_reference_choices(
                &"reverse_read".to_string(),
                &read,
                Some(qualities.clone()),
                &reference_manager,
                &fast_lookup,
                &layout,
                &mut alignment_matrix,
                &affine_score,
                &inversion_score,
                &false,
                &test_convex_scoring(),
                &false,
                2.0,
                0,
                &read.len(),
                None,
                None,
                None,
            )
            .expect("reverse-complemented read should align");
            let alignment = result.alignment.expect("alignment should be present");

            assert_eq!(result.ref_name, b"reverse_reference");
            assert_eq!(alignment.read_aligned, reference);
            assert_eq!(alignment.read_quals, Some(expected_qualities.clone()));
        }
    }

    #[test]
    fn test_fast_lookup_falls_back_on_single_unique_kmer_hit() {
        let layout = multi_reference_layout(true);
        let read = layout
            .references
            .get("reverse_reference")
            .unwrap()
            .sequence
            .as_bytes()
            .to_vec();
        let mut reference_manager = ReferenceManager::from_yaml_input(&layout, 8, 4);
        let false_reference = reference_manager
            .references
            .values()
            .find(|reference| reference.name == b"forward_reference")
            .unwrap()
            .clone();

        reference_manager.unique_kmers.kmer_to_reference.clear();
        reference_manager
            .unique_kmers
            .kmer_to_reference
            .insert(read[..8].to_vec(), false_reference);

        let mut alignment_matrix = create_scoring_record_3d(
            reference_manager.longest_ref + 1,
            read.len() + 1,
            AlignmentType::Affine,
            false,
        );
        let result = quick_alignment_search(
            &"sparse_evidence_read".to_string(),
            &read,
            None,
            &reference_manager,
            &mut alignment_matrix,
            &AffineScoring::default_dna(),
            &0.90,
        )
        .expect("exhaustive fallback should find a reference");

        assert_eq!(result.ref_name, b"reverse_reference");
    }

    #[test]
    fn test_fast_lookup_requires_three_distinct_dominant_hits() {
        assert!(!has_confident_kmer_support(1, 1, 0.90));
        assert!(!has_confident_kmer_support(2, 2, 0.90));
        assert!(has_confident_kmer_support(3, 3, 0.90));
        assert!(!has_confident_kmer_support(9, 10, 0.90));
    }

    #[test]
    fn test_find_best_reference() {
        let ref_location = &"test_data/test_best_alignment.fasta".to_string();
        let rm = ReferenceManager::from_fa_file(&ref_location, 8, 8);

        let read_one = "atggactatcatatgcttaccgtaacttgaaagtatttcgatttcttggctttatatatcttgtggaaaggacgaaacaccgGGTAGCAAACGTTTGGACGTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTTTCCTGCAGGAAACCCCGGGgaat".to_string().to_ascii_uppercase().into_bytes();

        let _read_structure = SequenceLayout {
            aligner: None,
            merge: None,
            reads: vec![ReadPosition::Read1 {
                orientation: AlignedReadOrientation::Forward,
            }],
            known_strand: true,
            references: BTreeMap::new(),
        };

        let mut read_mat = create_scoring_record_3d(
            read_one.len() + 100,
            read_one.len() + 100,
            AlignmentType::Affine,
            false,
        );

        let _my_score = InversionScoring {
            match_score: 9.0,
            mismatch_score: -21.0,
            gap_open: -25.0,
            gap_extend: -1.0,
            inversion_penalty: -40.0,
            min_inversion_length: 20,
        };

        let my_aff_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -9.0,
            special_character_score: 9.0,
            gap_open: -20.0,
            gap_extend: -1.0,
            final_gap_multiplier: 1.0,
        };

        let best_ref = exhaustive_alignment_search(
            &"testread".to_string(),
            &read_one,
            None,
            &&rm,
            &mut read_mat,
            &my_aff_score,
            None,
        );
        assert_eq!(
            String::from_utf8(best_ref.unwrap().ref_name).unwrap(),
            String::from_utf8(
                "1_AAACCCCGGG_GGTAGCAAACGTTTGGACGTG"
                    .to_string()
                    .into_bytes()
            )
            .unwrap()
        );

        let read_one = "atggactatcatatgcttaccgtaacttgaaagtatttcgatttcttggctttatatatcttgtggaaaggacgaaacaccgGGTGCCCTTACTCTCACCTGATTACTTAATCCGTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTTTCCTGCAGGAACGCCCTACgaattcgggcccattggtatggc".to_string().to_ascii_uppercase().into_bytes();
        let best_ref = exhaustive_alignment_search(
            &"testread".to_string(),
            &read_one,
            None,
            &&rm,
            &mut read_mat,
            &my_aff_score,
            None,
        );

        assert_eq!(
            String::from_utf8(best_ref.unwrap().ref_name).unwrap(),
            String::from_utf8(
                "2_AACGCCCTAC_GGTGCCCTTACTCTCACCTGATTACTTAATCCGTG"
                    .to_string()
                    .into_bytes()
            )
            .unwrap()
        );
    }

    #[test]
    fn test_find_best_reference2() {
        let ref_location = &"test_data/test_ref_alignment.fasta".to_string();
        let rm = ReferenceManager::from_fa_file(&ref_location, 8, 8);

        let read_one = "ATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCGGTAAATTTGAGGCTCCGGCATGCAGGAGGCCGTGGGGTTAGAGCTAGAAATAGCAAGTTAACCTAAGGCTAGTCCGTTATCAACTTG".to_string().to_ascii_uppercase().into_bytes();

        let _read_structure = SequenceLayout {
            aligner: None,
            merge: None,
            reads: vec![ReadPosition::Read1 {
                orientation: AlignedReadOrientation::Forward,
            }],
            known_strand: true,
            references: BTreeMap::new(),
        };

        let mut read_mat = create_scoring_record_3d(
            read_one.len() + 100,
            read_one.len() + 100,
            AlignmentType::Affine,
            false,
        );

        let _my_score = InversionScoring {
            match_score: 9.0,
            mismatch_score: -21.0,
            gap_open: -25.0,
            gap_extend: -1.0,
            inversion_penalty: -40.0,
            min_inversion_length: 20,
        };

        let my_aff_score = AffineScoring {
            match_score: 10.0,
            mismatch_score: -9.0,
            special_character_score: 9.0,
            gap_open: -20.0,
            gap_extend: -1.0,
            final_gap_multiplier: 1.0,
        };

        let best_ref = exhaustive_alignment_search(
            &"testread".to_string(),
            &read_one,
            None,
            &&rm,
            &mut read_mat,
            &my_aff_score,
            None,
        );
        assert_eq!(
            String::from_utf8(best_ref.unwrap().ref_name).unwrap(),
            String::from_utf8(
                "ref_48_GGTAAATTTGAGGCTCCGGCATGCAGGAGGCCGTG"
                    .to_string()
                    .into_bytes()
            )
            .unwrap()
        );
    }

    #[test]
    fn simplify_cigar_test() {
        let input_cigar = vec![
            AlignmentTag::MatchMismatch(1),
            AlignmentTag::MatchMismatch(1),
            AlignmentTag::MatchMismatch(1),
        ];
        let merged_cigar = vec![AlignmentTag::MatchMismatch(3)];
        let resulting_cigar = simplify_cigar_string(&input_cigar);
        assert_eq!(resulting_cigar, merged_cigar);

        let input_cigar = vec![
            AlignmentTag::MatchMismatch(1),
            AlignmentTag::Ins(1),
            AlignmentTag::MatchMismatch(1),
            AlignmentTag::MatchMismatch(1),
        ];
        let merged_cigar = vec![
            AlignmentTag::MatchMismatch(1),
            AlignmentTag::Ins(1),
            AlignmentTag::MatchMismatch(2),
        ];
        let resulting_cigar = simplify_cigar_string(&input_cigar);
        assert_eq!(resulting_cigar, merged_cigar);
    }

    #[test]
    fn test_simplify_cigar_empty() {
        let input: Vec<AlignmentTag> = vec![];
        let result = simplify_cigar_string(&input);
        assert!(result.is_empty());
    }

    #[test]
    fn test_simplify_cigar_single_element() {
        let input = vec![AlignmentTag::Del(5)];
        let result = simplify_cigar_string(&input);
        assert_eq!(result, vec![AlignmentTag::Del(5)]);
    }

    #[test]
    fn test_simplify_cigar_no_merge_needed() {
        let input = vec![
            AlignmentTag::MatchMismatch(3),
            AlignmentTag::Del(2),
            AlignmentTag::Ins(1),
            AlignmentTag::MatchMismatch(4),
        ];
        let result = simplify_cigar_string(&input);
        assert_eq!(result, input);
    }

    #[test]
    fn test_simplify_cigar_all_same_type() {
        let input = vec![
            AlignmentTag::Del(1),
            AlignmentTag::Del(2),
            AlignmentTag::Del(3),
        ];
        let result = simplify_cigar_string(&input);
        assert_eq!(result, vec![AlignmentTag::Del(6)]);
    }

    #[test]
    fn test_simplify_cigar_insertions() {
        let input = vec![
            AlignmentTag::Ins(1),
            AlignmentTag::Ins(1),
            AlignmentTag::Ins(1),
        ];
        let result = simplify_cigar_string(&input);
        assert_eq!(result, vec![AlignmentTag::Ins(3)]);
    }

    #[test]
    fn test_cigar_to_alignment_perfect_match() {
        let reference = b"ACGT".to_vec();
        let read = b"ACGT".to_vec();
        let cigar = vec![
            AlignmentOperation::Match,
            AlignmentOperation::Match,
            AlignmentOperation::Match,
            AlignmentOperation::Match,
        ];
        let (ref_aln, read_aln, tags) = cigar_to_alignment(&reference, &read, &cigar);
        assert_eq!(ref_aln, b"ACGT".to_vec());
        assert_eq!(read_aln, b"ACGT".to_vec());
        assert_eq!(tags, vec![AlignmentTag::MatchMismatch(4)]);
    }

    #[test]
    fn test_cigar_to_alignment_with_deletion() {
        let reference = b"ACGT".to_vec();
        let read = b"AT".to_vec();
        let cigar = vec![
            AlignmentOperation::Match,
            AlignmentOperation::Del,
            AlignmentOperation::Del,
            AlignmentOperation::Match,
        ];
        let (ref_aln, read_aln, tags) = cigar_to_alignment(&reference, &read, &cigar);
        assert_eq!(ref_aln, reference);
        assert_eq!(read_aln[0], b'A');
        assert_eq!(read_aln[3], b'T');
        assert_eq!(tags.len(), 3); // M, D, M
    }

    #[test]
    fn test_cigar_to_alignment_with_insertion() {
        let reference = b"AT".to_vec();
        let read = b"ACGT".to_vec();
        let cigar = vec![
            AlignmentOperation::Match,
            AlignmentOperation::Ins,
            AlignmentOperation::Ins,
            AlignmentOperation::Match,
        ];
        let (ref_aln, read_aln, tags) = cigar_to_alignment(&reference, &read, &cigar);
        assert_eq!(read_aln, read);
        assert_eq!(ref_aln[0], b'A');
        assert_eq!(ref_aln[3], b'T');
        assert_eq!(tags.len(), 3); // M, I, M
    }

    #[test]
    fn test_cigar_to_alignment_with_substitution() {
        let reference = b"ACGT".to_vec();
        let read = b"ATGT".to_vec();
        let cigar = vec![
            AlignmentOperation::Match,
            AlignmentOperation::Subst,
            AlignmentOperation::Match,
            AlignmentOperation::Match,
        ];
        let (ref_aln, read_aln, tags) = cigar_to_alignment(&reference, &read, &cigar);
        assert_eq!(ref_aln, reference);
        assert_eq!(read_aln, read);
        // Subst is treated as MatchMismatch, so all merge to M(4)
        assert_eq!(tags, vec![AlignmentTag::MatchMismatch(4)]);
    }

    #[test]
    fn test_convex_alignment_single_long_gap() {
        // A contiguous deletion should be called as one gap; quals pass through in
        // read order (convex never reverses a block).
        let reference = b"AGCTTGCATGCCTGCAGGTCGACTCTAGAGTCGACCTGCA".to_vec(); // 40 bp
        let mut read = reference[0..15].to_vec();
        read.extend_from_slice(&reference[27..]); // delete reference[15..27] (12 bp)
        let quals: Vec<u8> = (0..read.len() as u8).collect();

        let res = crate::alignment::alignment_matrix::convex_alignment(
            &reference, &read, &"r".to_string(), &"q".to_string(),
            &test_convex_scoring(), Some(quals.clone()),
        );
        let dels: Vec<usize> = res
            .cigar_string
            .iter()
            .filter_map(|t| if let AlignmentTag::Del(k) = t { Some(*k) } else { None })
            .collect();
        assert_eq!(dels, vec![12], "expected one 12bp deletion; cigar={:?}", res.cigar_string);
        assert_eq!(res.read_quals, Some(quals), "quals pass through in read order");
        // ungapped aligned strings reconstruct the inputs
        let ref_ung: Vec<u8> = res.reference_aligned.iter().filter(|b| **b != b'-').cloned().collect();
        let read_ung: Vec<u8> = res.read_aligned.iter().filter(|b| **b != b'-').cloned().collect();
        assert_eq!(ref_ung, reference);
        assert_eq!(read_ung, read);
    }

    fn align_single_forward_with_convex(reference: &str, read: Vec<u8>) -> AlignmentResult {
        let mut references = BTreeMap::new();
        references.insert(
            "reference".to_string(),
            ReferenceRecord {
                sequence: reference.to_string(),
                umi_configurations: BTreeMap::new(),
                targets: vec![],
                target_types: vec![],
                target_locations: None,
                prime_edits: BTreeMap::new(),
            },
        );
        let layout = SequenceLayout {
            aligner: None,
            merge: None,
            reads: vec![ReadPosition::Read1 { orientation: AlignedReadOrientation::Forward }],
            known_strand: true,
            references,
        };
        let reference_manager = ReferenceManager::from_yaml_input(&layout, 8, 4);
        let mut alignment_matrix = create_scoring_record_3d(
            reference_manager.longest_ref + 1,
            read.len() + 1,
            AlignmentType::Affine,
            false,
        );
        align_to_reference_choices(
            &"read".to_string(),
            &read,
            None,
            &reference_manager,
            &false,
            &layout,
            &mut alignment_matrix,
            &AffineScoring::default_dna(),
            &InversionScoring::default(),
            &false, // use_inversions
            &test_convex_scoring(),
            &true, // use_convex -> route to the convex-gap aligner
            2.0,
            0,
            &read.len(),
            None,
            None,
            None,
        )
        .unwrap()
        .alignment
        .unwrap()
    }

    #[test]
    fn test_convex_dispatch_aligns_deletion() {
        // `--aligner convex` (use_convex = true) routes the single-reference path to
        // the convex aligner and produces a valid alignment (one long gap).
        let reference = "AGCTTGCATGCCTGCAGGTCGACTCTAGAGTCGACCTGCA";
        let refb = reference.as_bytes();
        let mut read = refb[0..15].to_vec();
        read.extend_from_slice(&refb[27..]); // 12bp deletion
        let aln = align_single_forward_with_convex(reference, read);
        let dels: Vec<usize> = aln
            .cigar_string
            .iter()
            .filter_map(|t| if let AlignmentTag::Del(k) = t { Some(*k) } else { None })
            .collect();
        assert_eq!(dels, vec![12], "convex dispatch should call one 12bp deletion; cigar={:?}", aln.cigar_string);
    }

    #[test]
    fn test_inversion_detected_with_long_flanks() {
        // Regression for the long-flank detection fix: the candidate search uses a
        // near-gapless scoring so inversions are detected even with realistic
        // amplicon-length flanks (previously failed beyond ~10 bp of flank).
        let sc = InversionScoring{match_score:10.0,mismatch_score:-11.0,gap_open:-15.0,gap_extend:-5.0,inversion_penalty:-10.0,min_inversion_length:8};
        let aff = AffineScoring::default_dna();
        let bank = b"GCCTCCACGGCCACTAGTATTATGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACGTATTAGTCATCGCTATTACCATGTACTCA";
        for &(flank, midlen) in [(20usize, 20usize), (20, 25), (50, 25)].iter() {
            let mid = &bank[flank..flank + midlen];
            let f1 = &bank[0..flank];
            let f2 = &bank[flank + midlen..(flank + midlen + 20).min(bank.len())];
            let reference: Vec<u8> = [f1, mid, f2].concat();
            let read: Vec<u8> = [f1, reverse_complement(&mid.to_vec()).as_slice(), f2].concat();
            let res = crate::alignment::alignment_matrix::inversion_alignment(
                &reference, &read, &"r".to_string(), &"q".to_string(), &sc, &aff, false, None,
            );
            let detected = res.cigar_string.iter().any(|t| {
                matches!(t, AlignmentTag::InversionOpen | AlignmentTag::InversionClose)
            });
            assert!(
                detected,
                "inversion not detected at flank={} midlen={}; cigar={:?}",
                flank, midlen, res.cigar_string
            );
        }
    }

    #[test]
    fn test_inversion_threads_and_reverses_quals() {
        // Quality scores are threaded through the inversion aligner and, within an
        // inverted block, reversed so they stay paired with the (reverse-complemented)
        // aligned bases. Distinct per-base quals let us detect any loss or misorder.
        let bank = b"GCCTCCACGGCCACTAGTATTATGCCCAGTACATGACCTTATGGGACTTT";
        let (flank, midlen) = (5usize, 20usize);
        let mid = &bank[flank..flank + midlen];
        let f1 = &bank[0..flank];
        let f2 = &bank[flank + midlen..flank + midlen + 15];
        let reference: Vec<u8> = [f1, mid, f2].concat();
        let read: Vec<u8> = [f1, reverse_complement(&mid.to_vec()).as_slice(), f2].concat();
        let sc = InversionScoring{match_score:10.0,mismatch_score:-11.0,gap_open:-15.0,gap_extend:-5.0,inversion_penalty:-10.0,min_inversion_length:8};
        let quals: Vec<u8> = (0..read.len() as u8).collect();

        let res = crate::alignment::alignment_matrix::inversion_alignment(
            &reference, &read, &"r".to_string(), &"q".to_string(),
            &sc, &AffineScoring::default_dna(), false, Some(quals.clone()),
        );
        assert!(
            res.cigar_string.iter().any(|t| matches!(t, AlignmentTag::InversionOpen | AlignmentTag::InversionClose)),
            "fixture should detect an inversion; cigar={:?}", res.cigar_string
        );
        let rq = res.read_quals.expect("quals should be threaded through");
        // gapless clean inversion: one qual per read base, no loss
        assert_eq!(rq.len(), quals.len(), "qual length must match the read");
        let mut a = rq.clone(); a.sort();
        let mut b = quals.clone(); b.sort();
        assert_eq!(a, b, "quals must be a permutation of the input (none lost/duplicated)");
        assert_ne!(rq, quals, "the inverted block's quals must be reordered");
    }

    #[test]
    fn test_inversion_flattened_to_iv_event_and_serializes() {
        // Item B(ii): a detected inversion is flattened to an `iv` event and the
        // record serializes to a single BAM line (no to_op panic). Uses a
        // short-flank fixture the aligner reliably detects.
        let bank = b"GCCTCCACGGCCACTAGTATTATGCCCAGTACATGACCTTATGGGACTTT";
        let (flank, midlen) = (5usize, 20usize);
        let mid = &bank[flank..flank + midlen];
        let f1 = &bank[0..flank];
        let f2 = &bank[flank + midlen..flank + midlen + 15];
        let reference: Vec<u8> = [f1, mid, f2].concat();
        let read: Vec<u8> = [f1, reverse_complement(&mid.to_vec()).as_slice(), f2].concat();

        let mut alignment = align_single_forward_with_inversions(
            std::str::from_utf8(&reference).unwrap(),
            read,
        );
        // detection produced inversion markers
        assert!(
            alignment.cigar_string.iter().any(|t| matches!(
                t,
                AlignmentTag::InversionOpen | AlignmentTag::InversionClose
            )),
            "expected inversion markers before flattening; cigar = {:?}",
            alignment.cigar_string
        );
        // flatten -> `<len>V+<pos>` event(s), marker-free coalesced cigar
        let events = alignment.take_inversion_events();
        assert!(!events.is_empty(), "expected an inversion event");
        assert!(
            events[0].contains("V+"),
            "event should be <len>V+<pos>, got {:?}",
            events
        );
        assert!(
            !alignment.cigar_string.iter().any(|t| matches!(
                t,
                AlignmentTag::InversionOpen | AlignmentTag::InversionClose
            )),
            "flattened cigar must be marker-free; cigar = {:?}",
            alignment.cigar_string
        );
        // serializes to a single BAM record without panicking on the inversion tags
        let _record =
            alignment.to_sam_record(&0, &std::collections::HashMap::new(), None);
    }

}
