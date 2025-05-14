![Clique](https://github.com/mckennalab/clique/raw/main/clique_logo.png)

Clique aims to be a high-performance, unified package developed for the efficient processing and analysis of lineage data from Illumina and long-read technologies like Nanopore. It effectively collapses diverse unique molecular identifiers (UMIs) and static identifiers (SIs) and generates consensus sequences. The package aligns these consensus sequences to amplicon sequences, matching them to the best candidate sequence through a variety of alignment strategies, such as affine gap, convex alignments, and Hidden Markov Models (HMMs).

Clique supports a broad range of single-cell sequencing technologies, and it provides comprehensive statistics on read recovery, UMI complexity, alignment quality, and single-cell alignment. The package also enables the creation of lineage trees and statistics, empowering researchers with a deeper understanding of the genetic relationships between individual cells.


# Quick‑start: Running *clique* on a New Amplicon Experiment

*clique* is a **Rust + Python** toolkit for high‑accuracy, amplicon‑based variant calling. It takes raw long‑read data (ONT / PacBio), collapses identical molecules, filters background contaminants, then emits per‑molecule consensus calls and a final VCF—all in a single command chain.

Below is the minimum “happy‑path” workflow. You only need the Rust toolchain and Python ≥ 3.9.

---

## 1. Prepare a Custom Reference YAML

Define every constant and variable segment of your amplicon in a YAML file.  
Variable regions are marked with successive integers (`{1}`, `{2}`, …) as described in the wiki.

```yaml
name:   MY_AMPLICON
pieces:
  - seq: "GGGCGT"           # left flank – fixed
  - len: 12                 # variable segment {1}
  - seq: "AACCGGTT"         # internal anchor – fixed
  - len: 8                  # variable segment {2}
  - seq: "AGCT"             # right flank – fixed
```

Full spec: <https://github.com/mckennalab/clique/wiki/YAML-input-file>.  
Save the file as **`my_amplicon.yaml`**.

---

## 2. *(Optional)* Build a Global Background Reference

If the run is noisy (e.g. whole‑cell prep), add off‑target genomes or transcriptomes so *clique* can subtract reads that align elsewhere.

```bash
# Bundle your amplicon with hg38 + plasmid FASTA
python scripts/make_bg_reference.py        --target my_amplicon.fa        --genome hg38.fa        --output bg_reference.fa
```

Skip this step for clean, single‑amplicon runs.

---

## 3. Align Reads

You can point *clique* at any SAM/BAM/CRAM produced by minimap2, but the repo ships two ultra‑fast Rust aligners (`clique-align` for noisy ONT; `clique-align-hifi` for PacBio‑HiFi):

```bash
# Rust aligner (fastest for ONT reads)
clique-align  -r my_amplicon.fa   -i reads.fastq   -o aligned.sam

# or classic aligner
minimap2 -ax map-ont my_amplicon.fa reads.fastq > aligned.sam
```

---

## 4. Collapse Identical Molecules & Score Run Quality

```bash
clique-collapse      --input aligned.sam     --yaml  my_amplicon.yaml     --out   collapsed.bam     --report collapse_metrics.json
```

`clique-collapse` groups reads by tag, removes PCR/sequencing duplicates, and writes QC metrics (coverage per variable site, duplication rate, read‑length histograms, etc.).

---

## 5. Call Variants / Barcodes

```bash
clique-call      --bam   collapsed.bam     --yaml  my_amplicon.yaml     --output variants.vcf
```

The caller walks each molecule, builds a consensus for every `{n}` variable region, and outputs a standard **VCF** plus optional per‑site CSV summaries.

---

## Next Steps

* Inspect `collapse_metrics.json`—low per‑site coverage often signals primer failure or barcode dropout.  
* Tune caller thresholds with `--min-depth` and `--min-freq` for ultra‑rare‑variant detection.  
* For multi‑amplicon panels, feed a **panel YAML** (array of amplicons) to the same commands; *clique* auto‑routes reads.

For deeper dives—error‑modeling flags, GPU options, and interactive notebooks—see the project wiki and `examples/` directory.

Happy collapsing!
