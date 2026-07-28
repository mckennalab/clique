//! End-to-end integration tests for the `clique` binary.
//!
//! These are black-box tests: they build small fixture inputs (a read-structure
//! YAML plus a gzipped FASTQ) in a temp directory, invoke the actual compiled
//! `clique` binary as a subprocess, and assert on its outputs — the tab-separated
//! run summary (always) and, when `samtools` is on PATH, the BAM records.
//!
//! Run with `cargo test --test integration`. The `samtools`-gated assertions
//! skip with a printed note if `samtools` is unavailable, so the summary-based
//! checks still exercise the full pipeline everywhere.

extern crate flate2;
extern crate tempfile;

use std::collections::HashMap;
use std::io::Write;
use std::process::Command;

use flate2::write::GzEncoder;
use flate2::Compression;
use tempfile::TempDir;

// --- fixture helpers -------------------------------------------------------

fn clique_bin() -> &'static str {
    env!("CARGO_BIN_EXE_clique")
}

fn write_text(dir: &TempDir, name: &str, contents: &str) -> String {
    let path = dir.path().join(name);
    std::fs::write(&path, contents).expect("write fixture");
    path.to_string_lossy().into_owned()
}

/// Build `n` reads named `<prefix><i>` all carrying `seq`.
fn reads_named(prefix: &str, n: usize, seq: &str) -> Vec<(String, String)> {
    (0..n).map(|i| (format!("{}{}", prefix, i), seq.to_string())).collect()
}

/// Write a gzipped FASTQ. Each read is `(name, sequence)`; qualities are all Q40.
fn write_fastq(dir: &TempDir, name: &str, reads: &[(String, String)]) -> String {
    let path = dir.path().join(name);
    let file = std::fs::File::create(&path).expect("create fastq");
    let mut gz = GzEncoder::new(file, Compression::default());
    for (rname, seq) in reads {
        let qual: String = "I".repeat(seq.len());
        write!(gz, "@{}\n{}\n+\n{}\n", rname, seq, qual).expect("write fastq record");
    }
    gz.finish().expect("finish gzip");
    path.to_string_lossy().into_owned()
}

/// A single-reference read structure (no UMIs, no targets).
fn single_ref_yaml(sequence: &str) -> String {
    format!(
        "---\nmerge: \"ConcatenateBothForward\"\nknown_strand: true\nreads:\n  - !Read1\n    orientation: Forward\nreferences:\n  amplicon:\n    sequence: \"{}\"\n    targets: []\n    target_types: []\n    umi_configurations:\n",
        sequence
    )
}

fn run_align(yaml: &str, r1: &str, out_bam: &str, summary: &str, extra: &[&str]) -> std::process::Output {
    let mut cmd = Command::new(clique_bin());
    cmd.arg("align")
        .arg("--read-structure").arg(yaml)
        .arg("--read1").arg(r1)
        .arg("--output-bam-file").arg(out_bam)
        .arg("--summary-output").arg(summary);
    for a in extra {
        cmd.arg(a);
    }
    cmd.output().expect("run clique align")
}

/// Parse the tidy TSV run summary into a map keyed by (section, scope, metric).
fn parse_summary(path: &str) -> HashMap<(String, String, String), String> {
    let text = std::fs::read_to_string(path).expect("read summary");
    let mut map = HashMap::new();
    for (i, line) in text.lines().enumerate() {
        if i == 0 {
            assert_eq!(line, "section\tscope\tmetric\tvalue", "unexpected summary header");
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() == 4 {
            map.insert(
                (cols[0].to_string(), cols[1].to_string(), cols[2].to_string()),
                cols[3].to_string(),
            );
        }
    }
    map
}

fn summary_val(m: &HashMap<(String, String, String), String>, section: &str, scope: &str, metric: &str) -> String {
    m.get(&(section.to_string(), scope.to_string(), metric.to_string()))
        .unwrap_or_else(|| panic!("summary missing [{} / {} / {}]; have keys: {:?}", section, scope, metric, m.keys().collect::<Vec<_>>()))
        .clone()
}

/// `samtools view` -> rows of columns, or None if samtools is unavailable.
fn samtools_view(bam: &str) -> Option<Vec<Vec<String>>> {
    let out = Command::new("samtools").arg("view").arg(bam).output().ok()?;
    if !out.status.success() {
        return None;
    }
    let text = String::from_utf8_lossy(&out.stdout);
    Some(
        text.lines()
            .filter(|l| !l.is_empty())
            .map(|l| l.split('\t').map(|s| s.to_string()).collect())
            .collect(),
    )
}

/// Find an optional SAM tag value (e.g. key "ce" -> the part after "ce:Z:").
fn sam_tag(row: &[String], key: &str) -> Option<String> {
    for col in row.iter().skip(11) {
        if col.len() > 5 && &col[0..2] == key && &col[2..3] == ":" {
            return Some(col[5..].to_string());
        }
    }
    None
}

fn assert_success(out: &std::process::Output, what: &str) {
    assert!(
        out.status.success(),
        "{} failed (status {:?}). stderr:\n{}",
        what,
        out.status.code(),
        String::from_utf8_lossy(&out.stderr)
    );
}

// A clean 72 bp amplicon used across several tests.
const AMPLICON: &str =
    "GCCTCCACGGCCACTAGTATTATGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACG";

fn rc(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T', 'T' => 'A', 'G' => 'C', 'C' => 'G', other => other,
        })
        .collect()
}

// --- tests -----------------------------------------------------------------

#[test]
fn align_single_reference_all_map() {
    let dir = TempDir::new().unwrap();
    let yaml = write_text(&dir, "ref.yaml", &single_ref_yaml(AMPLICON));
    let reads = reads_named("read", 5, AMPLICON);
    let r1 = write_fastq(&dir, "r1.fastq.gz", &reads);
    let bam = dir.path().join("out.bam").to_string_lossy().into_owned();
    let summary = dir.path().join("summary.tsv").to_string_lossy().into_owned();

    let out = run_align(&yaml, &r1, &bam, &summary, &[]);
    assert_success(&out, "align");

    let s = parse_summary(&summary);
    assert_eq!(summary_val(&s, "Alignment totals", "All", "Input reads"), "5");
    assert_eq!(summary_val(&s, "Alignment totals", "All", "Aligned"), "5");
    assert_eq!(summary_val(&s, "Alignment totals", "All", "Aligned (%)"), "100.00");
    assert_eq!(summary_val(&s, "Aligned reads by reference", "amplicon", "Aligned reads"), "5");

    if let Some(rows) = samtools_view(&bam) {
        assert_eq!(rows.len(), 5, "expected 5 BAM records");
        assert!(rows.iter().all(|r| r[2] == "amplicon"), "all reads should map to 'amplicon'");
    } else {
        eprintln!("(samtools not found: skipping BAM-record checks)");
    }
}

#[test]
fn align_multi_reference_routing() {
    // two near-identical alleles differing at three positions -> POA auto-default
    let left = AMPLICON;
    let mut r: Vec<char> = left.chars().collect();
    for &p in &[21usize, 40, 60] {
        r[p] = if r[p] == 'A' { 'C' } else { 'A' }; // guarantee a real difference
    }
    let right: String = r.into_iter().collect();
    assert_ne!(left, right.as_str());

    let dir = TempDir::new().unwrap();
    let yaml_text = format!(
        "---\nmerge: \"ConcatenateBothForward\"\nknown_strand: true\nreads:\n  - !Read1\n    orientation: Forward\nreferences:\n  left:\n    sequence: \"{}\"\n    targets: []\n    target_types: []\n    umi_configurations:\n  right:\n    sequence: \"{}\"\n    targets: []\n    target_types: []\n    umi_configurations:\n",
        left, right
    );
    let yaml = write_text(&dir, "panel.yaml", &yaml_text);

    let mut reads = reads_named("L", 3, left);
    reads.extend(reads_named("R", 2, &right));
    let r1 = write_fastq(&dir, "r1.fastq.gz", &reads);
    let bam = dir.path().join("out.bam").to_string_lossy().into_owned();
    let summary = dir.path().join("summary.tsv").to_string_lossy().into_owned();

    let out = run_align(&yaml, &r1, &bam, &summary, &[]);
    assert_success(&out, "align multi-ref");

    let s = parse_summary(&summary);
    assert_eq!(summary_val(&s, "Alignment totals", "All", "Aligned"), "5");
    assert_eq!(summary_val(&s, "Aligned reads by reference", "left", "Aligned reads"), "3");
    assert_eq!(summary_val(&s, "Aligned reads by reference", "right", "Aligned reads"), "2");

    if let Some(rows) = samtools_view(&bam) {
        for row in rows {
            let expected = if row[0].starts_with('L') { "left" } else { "right" };
            assert_eq!(row[2], expected, "read {} routed to {} (expected {})", row[0], row[2], expected);
        }
    } else {
        eprintln!("(samtools not found: skipping routing RNAME checks)");
    }
}

#[test]
fn align_calls_base_edit_event() {
    // reference = 26 bp flank + 20 bp Cas9ABE target + 26 bp flank
    let flank_l = "GCCTCCACGGCCACTAGTATTATGCC";
    let target = "CCAAATAGCTAAGATGACAGG"; // 21 bp; an 'A' sits in the ABE window
    let flank_r = "TAATTCGAATTTAAATCGGATCCGCG";
    let reference = format!("{}{}{}", flank_l, target, flank_r);
    let edit_pos = flank_l.len() + 4; // an 'A' inside the target's editing window

    assert_eq!(reference.as_bytes()[edit_pos], b'A', "sanity: edit site is an A");
    let mut edited: Vec<u8> = reference.clone().into_bytes();
    edited[edit_pos] = b'G'; // ABE A->G
    let edited = String::from_utf8(edited).unwrap();

    let dir = TempDir::new().unwrap();
    let yaml_text = format!(
        "---\nmerge: \"ConcatenateBothForward\"\nknown_strand: true\nreads:\n  - !Read1\n    orientation: Forward\nreferences:\n  amplicon:\n    sequence: \"{}\"\n    targets: [\"{}\"]\n    target_types: [\"Cas9ABE\"]\n    umi_configurations:\n",
        reference, target
    );
    let yaml = write_text(&dir, "abe.yaml", &yaml_text);
    let reads = vec![("unedited".to_string(), reference.clone()), ("edited".to_string(), edited)];
    let r1 = write_fastq(&dir, "r1.fastq.gz", &reads);
    let bam = dir.path().join("out.bam").to_string_lossy().into_owned();
    let summary = dir.path().join("summary.tsv").to_string_lossy().into_owned();

    let out = run_align(&yaml, &r1, &bam, &summary, &[]);
    assert_success(&out, "align ABE");
    assert_eq!(summary_val(&parse_summary(&summary), "Alignment totals", "All", "Aligned"), "2");

    if let Some(rows) = samtools_view(&bam) {
        for row in &rows {
            let ce = sam_tag(row, "ce").unwrap_or_else(|| panic!("read {} has no ce tag", row[0]));
            if row[0] == "edited" {
                assert!(ce.contains("S+"), "edited read should carry a substitution event, got ce={}", ce);
                assert!(ce.contains("+G"), "the base edit should be to G, got ce={}", ce);
            } else {
                assert_eq!(ce, "NONE", "unedited read should have ce=NONE, got {}", ce);
            }
        }
    } else {
        eprintln!("(samtools not found: skipping ce-tag checks)");
    }
}

#[test]
fn align_inversion_emits_iv_tag() {
    // flanks + a >=8 bp invertible block; total >= min_read_length (50 bp)
    let f1 = &AMPLICON[0..10];
    let mid = &AMPLICON[10..35]; // 25 bp block that inverts
    let f2 = &AMPLICON[35..60]; // 25 bp -> 60 bp read
    let reference = format!("{}{}{}", f1, mid, f2);
    let read = format!("{}{}{}", f1, rc(mid), f2);

    let dir = TempDir::new().unwrap();
    let yaml = write_text(&dir, "ref.yaml", &single_ref_yaml(&reference));
    let r1 = write_fastq(&dir, "r1.fastq.gz", &[("inv".to_string(), read)]);
    let bam = dir.path().join("out.bam").to_string_lossy().into_owned();
    let summary = dir.path().join("summary.tsv").to_string_lossy().into_owned();

    let out = run_align(&yaml, &r1, &bam, &summary, &["--aligner", "inversion"]);
    assert_success(&out, "align inversion");
    assert_eq!(summary_val(&parse_summary(&summary), "Alignment totals", "All", "Aligned"), "1");

    if let Some(rows) = samtools_view(&bam) {
        let iv = sam_tag(&rows[0], "iv").expect("inverting read should carry an iv tag");
        assert!(iv.contains("V+"), "iv tag should be <len>V+<pos>, got {}", iv);
    } else {
        eprintln!("(samtools not found: skipping iv-tag check)");
    }
}

#[test]
fn align_convex_deletion_single_gap() {
    // read = reference with a contiguous 15 bp deletion; convex should call one gap
    let reference = AMPLICON;
    let read = format!("{}{}", &reference[..30], &reference[45..]);

    let dir = TempDir::new().unwrap();
    let yaml = write_text(&dir, "ref.yaml", &single_ref_yaml(reference));
    let r1 = write_fastq(&dir, "r1.fastq.gz", &[("del".to_string(), read)]);
    let bam = dir.path().join("out.bam").to_string_lossy().into_owned();
    let summary = dir.path().join("summary.tsv").to_string_lossy().into_owned();

    let out = run_align(&yaml, &r1, &bam, &summary, &["--aligner", "convex"]);
    assert_success(&out, "align convex");
    assert_eq!(summary_val(&parse_summary(&summary), "Alignment totals", "All", "Aligned"), "1");

    if let Some(rows) = samtools_view(&bam) {
        let cigar = &rows[0][5];
        assert!(cigar.contains("15D"), "convex should call one 15 bp deletion, got CIGAR {}", cigar);
        assert_eq!(cigar.matches('D').count(), 1, "expected exactly one deletion op, got CIGAR {}", cigar);
    } else {
        eprintln!("(samtools not found: skipping convex CIGAR check)");
    }
}

#[test]
fn align_then_collapse_degenerate_umi() {
    // reference with a 10 bp degenerate-UMI slot marked by '0'; reads carry two
    // distinct UMIs so collapse groups them into two consensus molecules.
    let flank_l = "GCCTCCACGGCCACTAGTATTATGCCCAGT";
    let umi_slot = "0000000000";
    let flank_r = "ACATGACCTTATGGGACTTTCCTACTTGGC";
    let reference = format!("{}{}{}", flank_l, umi_slot, flank_r);

    let dir = TempDir::new().unwrap();
    let yaml_text = format!(
        "---\nmerge: \"ConcatenateBothForward\"\nknown_strand: true\nreads:\n  - !Read1\n    orientation: Forward\nreferences:\n  amplicon:\n    sequence: \"{}\"\n    targets: []\n    target_types: []\n    umi_configurations:\n      molecule:\n        symbol: '0'\n        sort_type: \"DegenerateTag\"\n        length: 10\n        order: 0\n        max_distance: 2\n        maximum_subsequences: 10000\n        minimum_collapsing_difference: 3.0\n",
        reference
    );
    let yaml = write_text(&dir, "umi.yaml", &yaml_text);

    let umi_a = "AAAAAAAAAA";
    let umi_b = "TTTTTTTTTT";
    let read_a = format!("{}{}{}", flank_l, umi_a, flank_r);
    let read_b = format!("{}{}{}", flank_l, umi_b, flank_r);
    let mut reads = reads_named("A", 4, &read_a);
    reads.extend(reads_named("B", 3, &read_b));
    let r1 = write_fastq(&dir, "r1.fastq.gz", &reads);

    let aligned_bam = dir.path().join("aligned.bam").to_string_lossy().into_owned();
    let align_summary = dir.path().join("align_summary.tsv").to_string_lossy().into_owned();
    let out = run_align(&yaml, &r1, &aligned_bam, &align_summary, &[]);
    assert_success(&out, "align (umi)");
    assert_eq!(summary_val(&parse_summary(&align_summary), "Alignment totals", "All", "Aligned"), "7");

    // now collapse
    let collapsed_bam = dir.path().join("collapsed.bam").to_string_lossy().into_owned();
    let collapse_summary = dir.path().join("collapse_summary.tsv").to_string_lossy().into_owned();
    let out = Command::new(clique_bin())
        .arg("collapse")
        .arg("--read-structure").arg(&yaml)
        .arg("--input-bam-file").arg(&aligned_bam)
        .arg("--output-bam-file").arg(&collapsed_bam)
        .arg("--summary-output").arg(&collapse_summary)
        .output()
        .expect("run clique collapse");
    assert_success(&out, "collapse");

    let s = parse_summary(&collapse_summary);
    assert_eq!(summary_val(&s, "Collapse results", "All", "BAM records"), "7", "collapse should read 7 input records");
    assert_eq!(summary_val(&s, "Collapse results", "All", "Output records"), "2", "two distinct UMIs -> two consensus molecules");

    if let Some(rows) = samtools_view(&collapsed_bam) {
        assert_eq!(rows.len(), 2, "expected 2 consensus BAM records");
    } else {
        eprintln!("(samtools not found: skipping consensus-record count)");
    }
}
