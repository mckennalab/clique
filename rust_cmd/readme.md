# Clique command line aligner

This directory contains the Clique alignment and collapsing tool. The goal of this tool is to take amplication sequencing data, from either Illumina or long-read technologies like Nanopore, and 
produce consensus sequences that are aligned to the genome. Often these amplicons have internal structure, such as unique molecular identifiers (UMIs), sequences that mark individual integration
locations within the genome ('static IDs'), or cell identifiers that come from many single-cell sequencing experiments. You provide this layout as a YAML file (detailed below) which Clique uses
to collapse down reads to a consensus sequence, accounting for errors or other issues 

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

### UMI Configurations

Each UMI configuration contains:

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `symbol` | Char | Required | Symbol used in the reference to mark this UMI |
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
    sequence: "ATCG"
    targets: ["ATCG"]
    target_types: ["Cas9WT"]
    umi_configurations:
      cell_id:
        symbol: '*'
        sort_type: "KnownTag"
        length: 16
        order: 0
        max_distance: 2
      cell_umi:
        symbol: '&'
        sort_type: "DegenerateTag"
        length: 12
        order: 1
        max_distance: 2
      static_id:
        symbol: '$'
        sort_type: "DegenerateTag"
        length: 12
        order: 2
        max_distance: 2
```

## Validation Rules

- UMI configurations must have sequential order numbers starting at 0
- Target sequences and target type lists must be the same length
- Target sequences must be found within the reference sequence
- Reference sequence must contain all symbols used in UMI configurations

For a complete example, see the `test_data/test_layout.yaml` file. 
