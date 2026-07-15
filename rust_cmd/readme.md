# Clique command line aligner

This directory contains the Clique alignment and collapsing tool. The goal of this tool is to take amplication sequencing data, from either Illumina or long-read technologies like Nanopore, and 
produce consensus sequences that are aligned to the genome. Often these amplicons have internal structure, such as unique molecular identifiers (UMIs), sequences that mark individual integration
locations within the genome ('static IDs'), or cell identifiers that come from many single-cell sequencing experiments. You provide this layout as a YAML file (detailed below) which Clique uses
to collapse down reads to a consensus sequence, accounting for errors or other issues 

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
`--maximum-reads-before-downsampling`.

# Sequence Layout YAML Configuration

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
