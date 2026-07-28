# Clique

Clique is the core command-line tool for processing structured amplicon and
lineage-tracing sequencing data. It reads Illumina or long-read FASTQs,
assembles the configured read components, assigns each read to an amplicon
reference, writes an annotated alignment BAM, and optionally collapses
barcode/UMI families into consensus records in a second BAM.

The sequence layout is supplied as YAML. It describes:

- which physical reads are present (`Read1`, `Read2`, and optionally `Index1`);
- how those reads are oriented and assembled;
- the expected amplicon reference or reference panel;
- fixed-position cell barcodes, integration markers, and molecule UMIs; and
- CRISPR targets and their editing chemistries.

Clique is intended for targeted reads with a known internal layout. It is not a
general whole-genome or transcriptome aligner.

## Documentation Map

- [Workflow](#workflow): what `align` and `collapse` do.
- [Quick Start](#quick-start): a complete align/sort/collapse command sequence.
- [Command Reference](#command-reference): current CLI options and router
  behavior.
- [Run Summaries](#run-summaries): terminal tables and tidy TSV metrics.
- [BAM Output](#bam-output): record conventions and auxiliary tags.
- [Read-Structure YAML](#read-structure-yaml): layout and chemistry schema.
- [Realistic Example Configuration](#realistic-example-configuration):
  validated R1/I1/R2 lineage example.

## Workflow

The primary pipeline is:

```text
FASTQ R1/R2/I1
       |
       v
clique align
       |
       +-- assemble configured read components
       +-- choose the best amplicon reference
       +-- align and extract capture sequences
       +-- call configured editing events
       |
       v
per-read annotated BAM
       |
       +-- sort/index with samtools
       |
       v
clique collapse
       |
       +-- filter alignments
       +-- correct known and degenerate capture tags in YAML order
       +-- group reads by corrected tags
       +-- downsample large groups when requested
       +-- construct a consensus per molecule
       |
       v
per-molecule annotated consensus BAM
```

Clique also provides `genbank-to-yaml` for creating a starting read-structure
YAML from annotated GenBank records.

## Build and Requirements

Clique currently requires a nightly Rust toolchain because the crate uses
nightly language features. A native compiler toolchain and `libclang` may also
be required by transitive dependencies.

Build the optimized executable from this directory:

```bash
cargo build --release
./target/release/clique --help
```

The resulting binary is:

```text
target/release/clique
```

`samtools` is not required by the Clique executable itself, but it is
recommended for sorting, indexing, inspecting, and validating BAM files.

Set `RUST_LOG` to control runtime logging:

```bash
RUST_LOG=debug ./target/release/clique align --help
```

Clique defaults to `RUST_LOG=info` when the variable is unset.

## Quick Start

### 1. Align FASTQs

```bash
CLIQUE=./target/release/clique

"$CLIQUE" align \
  --read-structure read_structure.yaml \
  --read1 sample_R1.fastq.gz \
  --read2 sample_R2.fastq.gz \
  --index1 sample_I1.fastq.gz \
  --output-bam-file sample.aligned.bam \
  --summary-output sample.align.summary.tsv \
  --threads 8
```

Only `--read1` is required. Omit `--read2`, `--index1`, and `--index2` when
they are not part of the YAML layout.

### 2. Sort and Index

Clique does not guarantee coordinate-sorted output and does not create an
index:

```bash
samtools sort -@ 8 -o sample.aligned.sorted.bam sample.aligned.bam
samtools index sample.aligned.sorted.bam
```

An indexed input lets `collapse` fetch each configured reference directly.

### 3. Correct Captures and Build Consensus Records

```bash
"$CLIQUE" collapse \
  --read-structure read_structure.yaml \
  --input-bam-file sample.aligned.sorted.bam \
  --output-bam-file sample.consensus.bam \
  --summary-output sample.collapse.summary.tsv \
  --threads 8 \
  --min-aligned-bases 90 \
  --min-aligned-identity 0.90 \
  --maximum-reads-before-downsampling 40
```

Sort and index `sample.consensus.bam` separately if a downstream tool requires
coordinate order.

### 4. Inspect the Result

```bash
samtools view -H sample.consensus.bam
samtools view sample.consensus.bam | cut -f1-6,12- | head
```

## Command Reference

Use the executable for the definitive option list:

```bash
clique align --help
clique collapse --help
clique genbank-to-yaml --help
```

### `clique align`

`align` consumes synchronized FASTQ records and writes one annotated BAM record
for every read that passes the length filters and can be assigned and aligned
to a configured reference.

Alignment consists of three separate operations:

1. **Read assembly:** construct one logical sequence from the physical reads in
   the YAML.
2. **Reference routing:** select a candidate reference when the YAML contains a
   panel.
3. **Final alignment:** align the assembled sequence to the selected reference,
   extract captures, and annotate editing events.

#### Input and Output Options

| Option | Required | Default | Description |
| --- | --- | --- | --- |
| `--read-structure <YAML>` | Yes | - | YAML containing read assembly, references, captures, and targets. |
| `--read1 <FASTQ>` | Yes | - | Read 1 FASTQ. |
| `--read2 <FASTQ>` | No | `NONE` | Read 2 FASTQ when declared by the layout. |
| `--index1 <FASTQ>` | No | `NONE` | Index 1 FASTQ when declared by the layout. |
| `--index2 <FASTQ>` | No | `NONE` | Accepted by the CLI, but layouts containing `Index2` are not currently assembled by the merger. |
| `--output-bam-file <BAM>` | Yes | - | Per-read annotated BAM output. |
| `--summary-output <TSV>` | No | - | Long-form run summary in addition to the terminal tables. |
| `--threads <N>` | No | `1` | Rayon worker threads. Must be greater than zero. |

All provided FASTQs must contain the same number of records in the same order.
Companion records must have matching canonical read names. Standard `/1` and
`/2` suffixes are normalized when names are compared.

#### Length and Alignment Options

| Option | Default | Behavior |
| --- | --- | --- |
| `--min-read-length <N>` | `50` | Drop an assembled read shorter than `N` bases before reference alignment. |
| `--max-reference-multiplier <N>` | `2` | Drop an assembled read whose length is greater than or equal to `(longest reference + 1) * N`. |
| `--aligner <MODE>` | `wfa` | Select the final alignment mode. Values are `wfa`, `degenerate`, and `inversion`. |

Current aligner behavior:

| Mode | Current implementation |
| --- | --- |
| `wfa` | Default standard affine/global alignment path. |
| `degenerate` | Currently follows the same non-inversion alignment path as `wfa`; it is retained for CLI compatibility. |
| `inversion` | Enables experimental inversion-aware final alignment. Called inversions are flattened into an `iv` tag so the result remains a standard single BAM record. Validate this mode against controls before production use. |

Digits used as capture placeholders and `N` bases are treated as compatible
with any read base during alignment. They do not require the read to contain
the literal digit or `N`.

When `known_strand: false`, Clique evaluates the assembled sequence in forward
and reverse-complement orientation and retains the better alignment. When
`known_strand: true`, it uses the orientation produced by read assembly.

#### Reference Routing Options

No router is needed for a single-reference YAML. For multiple references,
Clique normally enables a POA graph automatically when the panel is compact:

- the panel contains at most 256 references;
- the POA graph is discriminable; and
- the graph contains at most four times as many nodes as the longest reference
  has bases.

If those conditions are not met, Clique falls back to its k-mer reference
search. The explicit router options are:

| Option | Default | Description |
| --- | --- | --- |
| `--poa-classifier` | Off | Force POA graph routing. It votes over discriminating graph columns and takes precedence over the other routers when a graph can be built. |
| `--poa-min-margin <N>` | `1` | Minimum top-two POA vote margin for a non-ambiguous call. |
| `--kmer-idf` | Off | Use IDF-weighted k-mers, assigning no weight to shared backbone k-mers. |
| `--kmer-idf-min-margin <F>` | `0.05` | Minimum `(best - second) / best` score ratio for a non-ambiguous IDF call. |
| `--discriminating-classifier` | Off | Compare only aligned positions at which panel references differ. |
| `--discriminating-min-margin <N>` | `1` | Minimum top-two discriminating-position margin for a non-ambiguous call. |
| `--no-poa-default` | Off | Disable automatic POA routing and use the legacy k-mer search unless another router is requested. |

`--kmer-idf` and `--discriminating-classifier` can be combined. In that case,
the discriminating-position classifier refines an ambiguous IDF call.
Ambiguous calls are still written; they are counted in the alignment summary
and marked in router-specific BAM tags.

Router confidence tags are described under [BAM Output](#bam-output).

### `clique collapse`

`collapse` reads the aligned BAM, filters records, corrects the capture tags in
their configured hierarchy, and writes either one consensus per corrected
molecule group or one corrected record per input record.

| Option | Required | Default | Description |
| --- | --- | --- | --- |
| `--read-structure <YAML>` | Yes | - | The same layout used for alignment. |
| `--input-bam-file <BAM>` | Yes | - | Aligned input BAM. `<BAM>.bai` is used when available. |
| `--output-bam-file <BAM>` | Yes | - | Corrected or consensus BAM output. |
| `--summary-output <TSV>` | No | - | Long-form run summary in addition to terminal tables. |
| `--threads <N>` | No | `1` | Worker threads. Must be greater than zero. |
| `--temp-dir <DIR>` | No | System temp | Parent directory for sharded sort and correction files. |
| `--min-aligned-bases <N>` | No | `45` | Minimum aligned, non-capture bases. The threshold is capped to the usable bases in each reference. |
| `--min-aligned-identity <F>` | No | `0.8` | Minimum identity among aligned, called, non-capture bases. Must be between 0 and 1. |
| `--maximum-reads-before-downsampling <N>` | No | `40` | Maximum reads selected for one consensus; `0` disables downsampling. |
| `--correct-only` | No | Off | Correct tags but do not combine records into consensuses. |

Correction follows increasing YAML `order`, not the order of configuration
keys or symbols. Each corrected tag becomes part of the grouping key for the
next level. This supports hierarchies such as:

```text
cell barcode -> integration marker -> molecule UMI
```

### `clique genbank-to-yaml`

`genbank-to-yaml` converts selected features from one annotated GenBank record
into a starting read-structure YAML:

```bash
clique genbank-to-yaml \
  --genbank annotated_amplicon.gb \
  --output read_structure.yaml \
  --tag lineage_target \
  --reference-name experiment_1
```

| Option | Default | Description |
| --- | --- | --- |
| `--genbank <FILE>` | Required | Annotated GenBank input. |
| `--output <FILE>` | Required | YAML output path. |
| `--tag <TEXT>` | `lineage_target` | Case-insensitive text that must occur in a selected feature name. |
| `--reference-name <NAME>` | GenBank LOCUS | Override the YAML reference name. |

Selected barcode/UMI features may use `clique_symbol`,
`clique_sort_type`, `clique_max_distance`, and `clique_file` qualifiers.
CRISPR target features may use `clique_type`. Prime-edit targets additionally
require `clique_edit_offset`, `clique_ref`, and `clique_alt`; optional
qualifiers include `clique_strand`, `clique_call_flank`, `clique_rtt`, and
`clique_scaffold`.

The generated layout defaults to one forward `Read1` with
`ConcatenateBothForward`. Adjust the `reads`, `merge`, and orientation fields
for the actual sequencing design.

## Run Summaries

Both `align` and `collapse` always print their summary tables to stderr.
`--summary-output <FILE>` additionally writes every non-empty table cell as a
tidy TSV:

```text
section    scope    metric    value
```

`section` is the terminal table title, `scope` is the first column of a row,
and `metric`/`value` are produced from the remaining columns. Percentages and
histogram strings are included in the TSV as displayed values.

### Alignment Summary

The `Alignment totals` table contains:

| Metric | Meaning |
| --- | --- |
| `Input reads` | Synchronized read sets consumed from R1 and its companions. |
| `Aligned` | BAM records successfully written. |
| `Aligned (%)` | `100 * Aligned / Input reads`. |
| `Too short` | Assembled reads below `--min-read-length`. |
| `Too long` | Assembled reads at or above the multiplier-derived maximum. |
| `Failed` | Reads for which reference assignment or final alignment did not produce a record. |
| `Ambiguous calls` | Written reads whose active reference router reported an ambiguous top call. |
| `Elapsed (s)` | Wall-clock alignment time. |

These outcomes reconcile as:

```text
Input reads = Aligned + Too short + Too long + Failed
```

The `Aligned reads by reference` table reports, for every reference:

- aligned read count;
- percentage of all written alignments; and
- a `#` histogram normalized to the largest reference count, with 40
  characters at full scale.

### Collapse Summary

The `Collapse results` table has an aggregate `All` row followed by one row per
reference:

| Metric | Meaning |
| --- | --- |
| `BAM records` | Records examined for that reference. |
| `Passed filters` | Records remaining after BAM flags, alignment quality, reconstruction, duplicate, and tag checks. |
| `After UMI` | Reads retained after all configured correction levels. |
| `UMI removed` | `Passed filters - After UMI`. |
| `Groups` | Corrected molecule groups submitted to the output stage. |
| `Reads selected` | Reads selected for consensus attempts, including groups whose consensus later fails. |
| `Downsampled` | Reads omitted by `--maximum-reads-before-downsampling`. |
| `Failed groups` | Molecule groups that did not produce an output record. |
| `Output records` | Consensus records, or corrected records in `--correct-only` mode. |

The `Collapse filtering` table separates:

- unmapped, secondary, and supplementary BAM flags;
- failures of `--min-aligned-bases` or `--min-aligned-identity`;
- alignment reconstruction or capture-extraction failures;
- duplicate records; and
- invalid capture tags.

When captures are configured, `UMI correction levels` reports the type, input
reads, output reads, and removed reads for each `reference / symbol` level.

In consensus mode, `Reads per consensus (rc tag)` groups successful output
records into bins `1`, `2`, `3-4`, `5-8`, `9-16`, and so on. Its `#` histogram
is normalized to the largest bin, with 40 characters at full scale. The final
`Run details` table reports mode (`collapse` or `correct`) and elapsed time.

## BAM Output

Both commands write records against the amplicon references declared in the
YAML. `RNAME` is the YAML reference name and `POS` is the 1-based position in
that amplicon, not necessarily a genomic chromosome coordinate.

Clique currently writes its custom auxiliary tags as SAM/BAM `Z` strings.
Numeric-looking values such as `rc:Z:10` and `rm:Z:0.98` must be parsed by
downstream tools when numeric values are required.

### Alignment BAM

Each successful `align` record retains the assembled read alignment and
contains:

| Tag | Meaning |
| --- | --- |
| `e<symbol>` | Capture sequence extracted from the input read. |
| `rc` | `1` for an uncollapsed alignment. |
| `ar` | Original input QNAME. |
| `rm` | Fraction of called, non-`N`, non-capture reference bases matching the read. |
| `as` | Clique internal alignment score; not the standard SAM `AS` contract. |
| `ce` | Chemistry-filtered CRISPR event calls when targets are configured. |
| `pe` | Prime-edit classifications when any target is `PrimeEdit`. |
| `iv` | Experimental inversion calls from `--aligner inversion`, encoded as `<length>V+<position>` and joined by `&`. |

Reference routers add confidence tags only when that router made the call:

| Tags | Meaning |
| --- | --- |
| `dm`, `di`, `da` | Discriminating-position margin, informative positions, and ambiguity flag (`1`/`0`). |
| `ib`, `im`, `ik`, `ia` | IDF best score, margin ratio, informative k-mers, and ambiguity flag. |
| `pb`, `pm`, `pi`, `pa` | POA best branch-vote score, margin, informative columns, and ambiguity flag. |

### Consensus BAM

Each normal `collapse` record represents one reference plus one complete set
of corrected capture values:

| Tag | Meaning |
| --- | --- |
| `e<symbol>` | Corrected capture value used in the molecule key. |
| `o<symbol>` | Original value from the representative selected read before correction. |
| `rc` | All retained reads assigned to the corrected group before downsampling. |
| `dc` | Reads actually selected to construct the consensus. |
| `ar` | Comma-separated QNAMEs selected for consensus construction. |
| `rm` | Reference match fraction for the output alignment. |
| `as`, `rs` | Clique internal alignment scores, written as strings. |
| `ce` | Event calls made from the consensus alignment. |
| `pe` | Prime-edit classifications made from the consensus alignment. |

For every successful consensus:

```text
number of comma-separated names in ar = dc <= rc
```

For example, `rc:Z:10`, `dc:Z:4`, and four names in `ar` mean ten reads formed
the molecule group and four were selected for consensus. The other six were
downsampled. The output QNAME is the first selected read, not a generated
molecule identifier.

`ar` lists selected reads only. It does not include group members omitted by
downsampling. In `--correct-only` mode, Clique writes one corrected record per
retained input record, so `rc`, `dc`, and the number of names in `ar` are all
one.

### Capture Tag Naming

Capture tags are derived from the YAML `symbol`, not the configuration name:

```yaml
umi_configurations:
  cell_barcode:
    symbol: '0'
    sort_type: KnownTag
    length: 16
    order: 0
    max_distance: 1
  molecule_umi:
    symbol: '1'
    sort_type: DegenerateTag
    length: 12
    order: 1
    max_distance: 1
```

This produces `e0`/`o0` for the cell barcode and `e1`/`o1` for the UMI.
`KnownTag` values can be corrected to an allowlist entry; `DegenerateTag`
values can be corrected to a clustered representative.

### CRISPR Event Encoding

The `ce` tag uses 0-based ungapped reference coordinates:

| Event | Encoding |
| --- | --- |
| Deletion | `<length>D+<reference-position>` |
| Insertion | `<length>I+<reference-position>+<inserted-bases>` |
| Substitution | `1S+<reference-position>+<alternate-base>` |
| No qualifying event | `NONE` |

Multiple events at one target are joined with `&`. Calls from multiple targets
are joined with `_` in YAML target order. Event filtering depends on
`target_types`; `ce` is not an unfiltered list of all alignment differences.

The `pe` tag follows target order and contains `WT`, `PRECISE`, `PARTIAL`,
`PRECISE_PLUS_BYPRODUCT`, `SCAFFOLD_INCORPORATION`, `INDEL`, `OTHER`, or
`NO_CALL` for a prime-edit target. Non-prime targets in a mixed target list are
represented by `NA`.

### SAM Conventions

- Alignment QNAME is the input read name; consensus QNAME is the first selected
  source read.
- `FLAG` is currently `0`. Clique normalizes reads into reference orientation,
  so reverse complementation is not represented by flag `0x10`.
- `MAPQ` is unset (`255`) and must not be interpreted as mapping confidence.
- Output is not guaranteed to be coordinate sorted.
- Clique does not automatically create a BAM index.

To recover selected FASTQ sequences from a consensus `ar` tag, see
[`rewind/extract_collapsed_reads.py`](rewind/extract_collapsed_reads.py).

## Read-Structure YAML

### Top-Level Fields

| Field | Required | Description |
| --- | --- | --- |
| `known_strand` | Yes | If `false`, test forward and reverse-complement orientation after read assembly. |
| `merge` | Multi-read layouts | `Align`, `Concatenate`, or `ConcatenateBothForward`. |
| `reads` | Yes | Physical reads and/or spacers used to construct one logical sequence. |
| `references` | Yes | Map from BAM reference name to reference configuration. |
| `aligner` | No | Parsed for compatibility; the `align --aligner` CLI option controls the current run. |

### Read Assembly

Read entries are tagged YAML objects:

```yaml
reads:
  - !Read1
    orientation: Forward
  - !Index1
    orientation: Forward
  - !Read2
    orientation: ReverseComplement
  - !Spacer
    spacer_sequence: "ACGT"
```

Sequence-bearing entries accept `Forward`, `Reverse`, `ReverseComplement`, or
`Unknown`. `Unknown` cannot be used by concatenation; use `known_strand: false`
when the assembled molecule can occur in either overall orientation.

Current merge behavior:

| Strategy | Behavior |
| --- | --- |
| `Align` | Overlap-merge an R1/R2 pair. The current implementation supports exactly R1 plus R2 for this strategy. |
| `Concatenate` | Append configured components in YAML order after applying each orientation. |
| `ConcatenateBothForward` | Currently follows the same ordered concatenation path; retained as a descriptive/legacy strategy name. |

Currently supported physical-read patterns are:

- R1 only;
- R1 + R2; and
- R1 + R2 + I1.

`Index2` is represented in the schema and accepted by the CLI, but the current
merger does not assemble a layout containing I2. `Spacer` is usable within a
supported concatenated multi-read pattern.

### Reference Records

Each entry under `references` defines:

| Field | Required | Description |
| --- | --- | --- |
| `sequence` | Yes | Expected assembled reference. Capture positions are replaced by their digit symbols. |
| `umi_configurations` | Yes | Map of capture definitions; use `{}` when there are none. |
| `targets` | Yes | Target sequences in output-call order; use `[]` when there are none. |
| `target_types` | Yes | One chemistry per target. Must match `targets` length. |
| `target_locations` | No | Explicit 0-based target starts. Recommended and required to select a specific repeated occurrence. |
| `prime_edits` | Prime editing only | Map from zero-based target index to expected prime-edit specification. |

The reference `sequence` contains ordinary DNA bases plus digit placeholders.
For example, sixteen `0` characters identify the 16 bases extracted for the
capture whose symbol is `0`. Each symbol must appear exactly `length` times.

### Capture and UMI Configuration

| Field | Required | Description |
| --- | --- | --- |
| `symbol` | Yes | Unique ASCII digit `0`-`9`; also determines BAM tags such as `e0`. |
| `sort_type` | Yes | `KnownTag` for allowlist correction or `DegenerateTag` for abundance-based clustering. |
| `length` | Yes | Expected capture length and number of symbol occurrences in the reference. |
| `order` | Yes | Correction hierarchy, sequential from zero. |
| `max_distance` | Yes | Maximum edit distance accepted during correction and extraction-length validation. |
| `file` | `KnownTag` only | Allowlist with one sequence per line and no header. Every sequence must match `length`. |
| `reverse_complement_sequences` | No | Reverse-complement allowlist sequences when loading them. Default `false`. |
| `levenshtein_distance` | No | Use Levenshtein correction when `true` (default); use Hamming correction when `false`. |
| `max_gaps` | No | Maximum alignment gaps accepted while extracting this capture. |
| `maximum_subsequences` | No | In-memory read threshold for one parent group before correction spills to disk. Default `1,000,000`; this is not a read-removal limit. |
| `minimum_collapsing_difference` | No | Minimum abundance ratio used when collapsing degenerate-tag neighbors. Default `5.0`. |
| `pad` | No | Parsed as `Left` or `Right`, but currently has no correction behavior. |

### Target Types

Editing windows are inclusive offsets relative to the target start:

| Type | Calls |
| --- | --- |
| `Static` | No event calls. |
| `Cas9WT`, `Cas9Homing` | Insertions/deletions overlapping offsets 14-19. |
| `Cas12AWT` | Insertions/deletions overlapping offsets 14-23. |
| `Cas9ABE`, `Cas12ABE`, `Cas9ABEPalindrome` | A-to-G or complementary T-to-C substitutions at offsets 2-19. |
| `Cas9CBE`, `Cas12CBE` | C-to-T or complementary G-to-A substitutions at offsets 2-19. |
| `Cas9ABECBE`, `Cas12ABECBE` | Either accepted ABE or CBE substitution class at offsets 2-19. |
| `PrimeEdit` | Raw differences in `ce` plus expected-haplotype classification in `pe`. |

Cas12 base-editor types currently use the same 2-19 window as the corresponding
Cas9 type. Calibrate target-level interpretation against appropriate controls.

### Prime Editing

Every `PrimeEdit` target requires a `prime_edits` entry with the same
zero-based target index:

```yaml
references:
  prime_edit_reference:
    sequence: "AAAACCCCGGGGTTTT"
    umi_configurations: {}
    targets: ["AAAACCCCGGGGTTTT"]
    target_types: ["PrimeEdit"]
    target_locations: [0]
    prime_edits:
      0:
        edit_offset: 6
        reference: "CC"
        alternate: "TT"
        strand: Forward
        call_flank: 3
        rtt_sequence: "TTGG"          # optional
        scaffold_sequence: "AACCGG"  # optional
```

Coordinates and alleles use forward-reference orientation. Use an empty
`reference` for a programmed insertion and an empty `alternate` for a
programmed deletion. `strand` controls partial-incorporation direction.

## Realistic Example Configuration

This example models a targeted 10x-style lineage assay:

| Input | Content |
| --- | --- |
| R1 | 16-bp cell barcode followed by a 12-bp molecule UMI |
| I1 | Fixed 8-bp library index `ATCACGTA` |
| R2 | 8-bp integration marker followed by a two-target lineage cassette |

The assembled reference is therefore:

```text
[cell: 0 x 16][UMI: 1 x 12][I1 sequence][integration: 2 x 8][lineage cassette]
```

```yaml
---
known_strand: true
merge: "ConcatenateBothForward"
reads:
  - !Read1
    orientation: Forward
  - !Index1
    orientation: Forward
  - !Read2
    orientation: Forward

references:
  tenx_lineage:
    sequence: "0000000000000000111111111111ATCACGTA22222222GCTAGCTACGATCGTACGTAGACCTGATCGTACGATCGTACGGTTGCAACGTAGGCTAACGTACTAGGTCAGTACGCTAGTCAGGAAACGTTGACTGATCGTAGCTACGTA"

    umi_configurations:
      cell_barcode:
        symbol: '0'
        file: "cell_barcodes.txt"
        sort_type: "KnownTag"
        length: 16
        order: 0
        max_distance: 1

      integration_marker:
        symbol: '2'
        file: "integration_markers.txt"
        sort_type: "KnownTag"
        length: 8
        order: 1
        max_distance: 1

      molecule_umi:
        symbol: '1'
        sort_type: "DegenerateTag"
        length: 12
        order: 2
        max_distance: 1
        maximum_subsequences: 1000

    targets:
      - "GACCTGATCGTACGATCGTACGG"
      - "CTAGGTCAGTACGCTAGTCAGGA"
    target_types:
      - "Cas9WT"
      - "Cas9WT"
    target_locations:
      - 64
      - 107
```

This produces the correction hierarchy:

```text
cell barcode (e0) -> integration marker (e2) -> molecule UMI (e1)
```

The physical sequence order does not have to match correction order; `order`
controls the hierarchy.

This exact design is generated and tested by
[`test_harness/tenx_lineage`](test_harness/tenx_lineage):

```bash
./test_harness/tenx_lineage/run_harness.sh /tmp/clique-tenx-lineage
```

The harness creates synchronized FASTQs, allowlists, truth tables, the YAML,
aligned and collapsed BAMs, and validates capture correction, consensus depth,
source-read provenance, and insertion/deletion calls.

## Configuration Validation

Clique validates the layout when loading it:

- capture `order` values must be sequential from zero;
- capture symbols must be unique ASCII digits;
- each reference must contain exactly `length` copies of every symbol;
- `KnownTag` configurations must provide an allowlist for collapse;
- allowlist sequences must all equal the configured capture length;
- target and target-type lists must have the same length;
- explicit target locations must match the target list length and reference
  bases;
- inferred targets must occur in the reference;
- targets and capture spans must not overlap in GenBank conversion; and
- every `PrimeEdit` target must have a valid matching specification whose
  reference allele matches the configured reference.

Validate a new experimental layout on a small read subset before processing a
full run.
