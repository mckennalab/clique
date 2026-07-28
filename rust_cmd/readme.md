# Clique command line aligner

This directory contains the Clique alignment and collapsing tool. The goal of this tool is to take amplication sequencing data, from either Illumina or long-read technologies like Nanopore, and 
produce consensus sequences that are aligned to the genome. Often these amplicons have internal structure, such as unique molecular identifiers (UMIs), sequences that mark individual integration
locations within the genome ('static IDs'), or cell identifiers that come from many single-cell sequencing experiments. You provide this layout as a YAML file (detailed below) which Clique uses
to collapse down reads to a consensus sequence, accounting for errors or other issues 

## Synthetic 10x-Style Test Harness

[`test_harness/tenx_lineage`](test_harness/tenx_lineage) provides a deterministic
end-to-end example with a 16-bp cell barcode and 12-bp UMI in R1, an I1 library
index, and an integration marker plus two-target lineage cassette in R2. It
generates synchronized FASTQs and truth tables, runs alignment and collapse,
and validates capture correction, molecule depths, source-read provenance, and
simulated insertion/deletion calls:

```bash
./test_harness/tenx_lineage/run_harness.sh /tmp/clique-tenx-lineage
```

## Run Summaries

The `align` and `collapse` commands print a table to stderr when processing
finishes. Add `--summary-output <FILE>` to also write the results as a tidy TSV
with `section`, `scope`, `metric`, and `value` columns:

```bash
clique align \
  --read-structure layout.yaml \
  --read1 reads.fastq.gz \
  --output-bam-file aligned.bam \
  --summary-output aligned.summary.tsv

clique collapse \
  --read-structure layout.yaml \
  --input-bam-file aligned.bam \
  --output-bam-file collapsed.bam \
  --summary-output collapsed.summary.tsv
```

The alignment summary reports input read sets, written alignments, reads below
or above the configured length limits, failed alignments/reference assignments,
ambiguous router calls, elapsed time, and aligned reads per reference.
The per-reference table includes a `#` histogram normalized to the largest
reference count (40 characters at full scale).

The collapse summary reports records examined for each configured reference,
filtering outcomes, reads retained after each UMI correction level, molecule
groups attempted, reads selected for consensus, reads omitted by downsampling,
failed groups, and output records. `Reads selected` includes consensus attempts
that later fail; `Downsampled` is the exact number omitted by
`--maximum-reads-before-downsampling`. Successful consensus records are also
grouped into power-of-two `rc` read-count ranges (`1`, `2`, `3-4`, `5-8`, and
so on), with a `#` histogram normalized to the largest range (40 characters at
full scale).

## BAM Output Format

Both `align` and `collapse` write BAM records against the amplicon references
declared in the sequence-layout YAML. These are reference-relative amplicon
alignments, not genomic-coordinate alignments: `RNAME` is the YAML reference
name and `POS` is the 1-based position within that reference.

Clique currently writes all of its custom auxiliary tags as SAM/BAM `Z`
(string) values. This includes numeric-looking values such as `rc:Z:10` and
`rm:Z:0.98`; downstream programs should parse those strings when numeric values
are required. The tag names are lowercase and should not be confused with
similarly named standard uppercase SAM tags.

### Collapse Output

Each normal collapse record represents one corrected barcode/UMI group. The
group is defined by its reference and the ordered set of corrected capture tags
(`e<symbol>`). The most important collapse annotations are:

| Tag | Meaning |
|-----|---------|
| `e<symbol>` | Extracted tag after barcode/UMI correction. This is part of the key used to form the collapsed group. |
| `o<symbol>` | Original extracted value, before correction, from the representative read carried into the output record. It is not a list of every original value in the group. |
| `rc` | Read count: all retained reads assigned to the corrected group before consensus downsampling. |
| `dc` | Downsampled consensus depth: the number of reads actually selected to build the consensus. |
| `ar` | Comma-separated QNAMEs of the reads selected to build the consensus, in selection order. Reads omitted by downsampling are not included. |
| `rm` | Fraction of called, non-`N`, non-capture reference bases that match the output alignment. |
| `as`, `rs` | Clique's internal alignment score, written as a string. These are primarily diagnostic and are not standard SAM `AS` values. |
| `ce` | Chemistry-aware CRISPR event calls, when the reference has configured targets. |
| `pe` | Prime-edit classification, present when the reference has at least one `PrimeEdit` target. |

For every successfully written collapse record:

```text
number of comma-separated names in ar = dc <= rc
```

For example, `rc:Z:10`, `dc:Z:4`, and four names in `ar` means that ten reads
formed the molecule group, four were selected for consensus, and six were
omitted because of `--maximum-reads-before-downsampling`.

The output QNAME is the name of the first selected read. It is a representative
identifier, not a newly generated molecule ID. Likewise, `o<symbol>` comes from
that representative selected record. Use `ar` when complete selected-read
provenance is needed. In `--correct-only` mode, Clique writes one corrected
record per retained input record, so `rc`, `dc`, and the number of names in
`ar` are all one.

The collapse `SEQ`, quality values, and CIGAR describe the generated consensus
in reference orientation. Insertion columns are retained only when supported by
the consensus threshold, so an unsupported raw-read overhang or insertion can
be absent from the output sequence even though the originating read is listed
in `ar`.

### Capture Tags

Capture-tag names are derived from the single-digit `symbol` in each YAML
`umi_configurations` entry, not from the human-readable configuration key. For
example:

```yaml
umi_configurations:
  cell_barcode:
    symbol: '0'
    sort_type: KnownTag
    length: 16
    order: 0
  molecule_umi:
    symbol: '1'
    sort_type: DegenerateTag
    length: 12
    order: 1
```

produces tags such as `e0`, `o0`, `e1`, and `o1`. In an `align` BAM,
`e<symbol>` is the value extracted directly from that read and no
`o<symbol>` tag is written. In a `collapse` BAM, `e<symbol>` is the corrected
group value and `o<symbol>` preserves the representative read's value before
correction. A `KnownTag` value may therefore change to its allowlisted barcode,
and a `DegenerateTag` value may change to its clustered representative.

### CRISPR Event Annotations

The `ce` tag uses 0-based, ungapped reference coordinates and encodes events as:

| Event | Encoding |
|-------|----------|
| Deletion | `<length>D+<reference-position>` |
| Insertion | `<length>I+<reference-position>+<inserted-bases>` |
| Substitution | `1S+<reference-position>+<alternate-base>` |
| No called event | `NONE` |

Multiple events at one target are joined with `&`. Calls for multiple targets
are joined with `_` in the same order as `targets` in the YAML. Event filtering
uses the configured target chemistry, so `ce` is not an unfiltered list of
every alignment difference.

The `pe` tag follows the same underscore-separated target order and contains
one of `WT`, `PRECISE`, `PARTIAL`, `PRECISE_PLUS_BYPRODUCT`,
`SCAFFOLD_INCORPORATION`, `INDEL`, `OTHER`, or `NO_CALL` for a prime-edit
target. A non-prime target in a mixed target list is represented by `NA`.

### Align-Only Routing Annotations

An `align` BAM has `rc:Z:1`, the originating read name in `ar`, and an extracted
`e<symbol>` value for every configured capture. Depending on which optional
reference classifier made the routing decision, it can also contain:

| Tags | Meaning |
|------|---------|
| `dm`, `di`, `da` | Discriminating-position top-two margin, informative positions covered, and ambiguity flag (`1` or `0`). |
| `ib`, `im`, `ik`, `ia` | IDF best score, top-two margin ratio, informative k-mers, and ambiguity flag. |
| `pb`, `pm`, `pi`, `pa` | POA best branch-vote score, top-two margin, informative columns, and ambiguity flag. |

These routing diagnostics describe reference assignment and are not propagated
as collapse-group annotations.

### SAM Field Conventions and Inspection

- `QNAME` is the input read name for alignments and the first selected read name
  for consensuses.
- `FLAG` is currently `0`; reads are normalized into reference orientation, so
  the reverse-complement operation is not represented with SAM flag `0x10`.
- `MAPQ` is unset (`255`) and should not be interpreted as mapping confidence.
- `RNAME`, `POS`, and `CIGAR` refer to a configured amplicon reference rather
  than a chromosome.
- Output order is not guaranteed to be coordinate sorted, and Clique does not
  create a BAM index automatically.

The header and selected record fields/tags can be inspected with `samtools`:

```bash
samtools view -H collapsed.bam
samtools view collapsed.bam | cut -f1-6,12- | head
```

To recover the FASTQ sequences named in a collapsed record's `ar` tag, use the
included helper:

```bash
./rewind/extract_collapsed_reads.py \
  --read-name <consensus-qname> \
  --collapsed-bam collapsed.bam \
  --input-fastq reads.fastq.gz \
  --output contributing_sequences.txt
```

Sort before creating an index for coordinate-based tools:

```bash
samtools sort -o collapsed.sorted.bam collapsed.bam
samtools index collapsed.sorted.bam
```

## Sequence Layout YAML Configuration

This document describes the YAML configuration format used to define sequence layouts for read processing.

## Overview

The sequence layout configuration is used to specify how to extract information from sequencing reads, including UMIs (Unique Molecular Identifiers) and target sequences. This configuration is specific to each sequencing platform and type (e.g., 10X, sci, etc.).

## Configuration Structure

### Top-Level Fields

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `aligner` | String | Optional | Specifies which alignment tool to use |
| `merge` | Enum | Optional | Strategy to merge reads: `Align`, `Concatenate`, or `ConcatenateBothForward` |
| `reads` | Array | Required | Defines the read positions required for this configuration |
| `known_strand` | Boolean | Required | Indicates whether the strand orientation is known |
| `references` | Map | Required | Contains reference records (key is reference name) |

### Read Positions

The `reads` field contains an array of read positions, which can be one of:

- `!Read1`: First read
- `!Read2`: Second read
- `!Index1`: First index read
- `!Index2`: Second index read
- `!Spacer`: A spacer sequence

Read position entries can include:
- `chain_align`: Optional boolean indicating whether to chain alignment
- `orientation`: Orientation of the read (`Forward`, `Reverse`, `ReverseComplement`, or `Unknown`)
- `spacer_sequence`: For Spacer types, the actual sequence to use

### Reference Records

Each reference record contains:

| Field | Type | Description |
|-------|------|-------------|
| `sequence` | String | The reference sequence |
| `umi_configurations` | Map | UMI configurations (key is UMI name) |
| `targets` | Array | List of target sequence strings |
| `target_types` | Array | List of target types (must match length of targets) |
| `target_locations` | Array | Optional zero-based target starts; required to disambiguate a specific repeated occurrence |
| `prime_edits` | Map | Prime-edit specifications keyed by zero-based target index |

### UMI Configurations

Each UMI configuration contains:

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `symbol` | Char | Required | Unique ASCII digit (`0`-`9`) used in the reference and BAM tag |
| `file` | String | Optional | Path to file containing known sequences (one per line, no header) |
| `reverse_complement_sequences` | Boolean | Optional | Whether to reverse complement sequences from the file |
| `sort_type` | Enum | Required | Either `KnownTag` or `DegenerateTag` |
| `length` | Integer | Required | Length of the UMI sequence |
| `order` | Integer | Required | Order of UMIs (must be sequential starting at 0) |
| `pad` | Enum | Optional | Padding direction: `Left` or `Right` |
| `max_distance` | Integer | Required | Maximum edit distance for matching |
| `maximum_subsequences` | Integer | Optional | Maximum number of subsequences to consider |
| `max_gaps` | Integer | Optional | Maximum number of gaps allowed |

### Target Types

Supported target types:
- `Static`
- `Cas9WT`
- `Cas12AWT`
- `Cas9ABE`
- `Cas9CBE`
- `Cas9ABECBE`
- `Cas12ABE`
- `Cas12CBE`
- `Cas12ABECBE`
- `Cas9Homing`
- `Cas9ABEPalindrome`
- `PrimeEdit`

### Prime Editing

Prime-edit alleles use forward-reference coordinates and sequences. Each
`PrimeEdit` target must have a matching entry in `prime_edits`:

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
        strand: "Forward"
        call_flank: 3
        rtt_sequence: "TTGG"          # optional
        scaffold_sequence: "AACCGG"  # optional
```

Use an empty `reference` for a programmed insertion and an empty `alternate`
for a programmed deletion. `strand` controls the direction used to recognize
partial incorporation. The existing `ce` BAM tag retains raw events, while
`pe` contains one of `WT`, `PRECISE`, `PARTIAL`,
`PRECISE_PLUS_BYPRODUCT`, `SCAFFOLD_INCORPORATION`, `INDEL`, `OTHER`, or
`NO_CALL`. Non-prime targets are represented as `NA` when a read has mixed
target types.

## Example Configuration

```yaml
---
known_strand: true
merge: "Concatenate"
reads:
  - !Read1
    orientation: Forward
  - !Read2
    orientation: Forward
references:
  shorter_reference:
    sequence: "0000000000000000ATCG111111111111222222222222"
    targets: ["ATCG"]
    target_types: ["Cas9WT"]
    umi_configurations:
      cell_id:
        symbol: '0'
        file: "cell_barcodes.txt"
        sort_type: "KnownTag"
        length: 16
        order: 0
        max_distance: 2
      cell_umi:
        symbol: '1'
        sort_type: "DegenerateTag"
        length: 12
        order: 1
        max_distance: 2
      static_id:
        symbol: '2'
        sort_type: "DegenerateTag"
        length: 12
        order: 2
        max_distance: 2
```

## Validation Rules

- UMI configurations must have sequential order numbers starting at 0
- UMI symbols must be unique ASCII digits (`0`-`9`), allowing at most 10 UMIs per reference
- Target sequences and target type lists must be the same length
- Target locations, when supplied, must match the target list length and reference sequence
- Target sequences must be found within the reference sequence
- Reference sequence must contain all symbols used in UMI configurations
- Each `PrimeEdit` target must have exactly one valid `prime_edits` entry whose reference allele matches the configured reference

For a complete example, see the `test_data/test_layout.yaml` file. 
