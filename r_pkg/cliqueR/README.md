# cliqueR

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

`cliqueR` is the R analysis and orchestration package for
[clique](../../rust_cmd), a Rust command-line toolkit for amplicon sequencing,
barcode and UMI correction, consensus generation, and CRISPR lineage tracing.

The division of work is intentional:

- The Rust CLI aligns reads, extracts and corrects captures, collapses molecule
  groups, builds consensus reads, and calls editing events.
- `cliqueR` creates basic configurations and reproducible run scripts, invokes
  the CLI, imports its BAM annotations, constructs lineage character matrices,
  reconstructs trees, and visualizes results.

This package is under active development. The
[current feature status](#current-feature-status) distinguishes working
interfaces from exported placeholders.

## Contents

- [Current feature status](#current-feature-status)
- [Installation and requirements](#installation-and-requirements)
- [Quick start](#quick-start)
- [Configuration](#configuration)
- [Running clique](#running-clique)
- [Reading BAM output](#reading-bam-output)
- [Character matrices](#character-matrices)
- [Tree reconstruction](#tree-reconstruction)
- [Tree visualization](#tree-visualization)
- [Examples and further documentation](#examples-and-further-documentation)
- [Function index](#function-index)
- [Current limitations](#current-limitations)

## Current feature status

| Area | Available now | Status |
| --- | --- | --- |
| CLI discovery | `clique_binary()`, `set_clique_binary()` | Implemented |
| Direct CLI execution | `clique_run()` | Implemented |
| Basic read-structure generation | `fasta_targets_to_yaml()` | Implemented |
| Reproducible align/collapse scripts | `generate_clique_script()` | Implemented |
| Clique BAM import | `read_longread_lineage()` | Implemented |
| Edit character matrices | `build_indel_matrix()`, `indel_distance()` | Implemented |
| Multi-method tree orchestration | `build_trees()` | Implemented |
| Tree reconstruction | NJ, BIONJ, parsimony, IQ-TREE 2, VINE, PHYLIP MIX, Cassiopeia | Implemented; dependencies vary |
| Tree visualization | `plot_lineage_tree()`, `clique_categorical_palette()` | Implemented |
| Typed `align()` and `collapse()` wrappers | Exported API placeholders | Not implemented |
| 10x/Cell Ranger import | `read_10x_lineage()` placeholder | Not implemented |
| Tree comparison API | RF, triplet, quartet, and ancestor metrics | Not implemented |

For alignment and collapse, use `generate_clique_script()` or `clique_run()`.
Do not use the placeholder `align()` and `collapse()` functions in this
development release.

## Installation and requirements

### Install the R package

From a clone of this repository:

```sh
R CMD INSTALL r_pkg/cliqueR
```

Or install the package subdirectory from GitHub:

```r
remotes::install_github("aaronmck/clique", subdir = "r_pkg/cliqueR")
```

Load the package:

```r
library(cliqueR)
```

### Build the Rust CLI

The R package does not build or install the `clique` executable. From the
repository root:

```sh
cargo build --release --manifest-path rust_cmd/Cargo.toml
export CLIQUE_BIN="$PWD/rust_cmd/target/release/clique"
```

`cliqueR` resolves the executable in this order:

1. A path set for the current R session with `set_clique_binary()`.
2. The `CLIQUE_BIN` environment variable.
3. A `clique` executable on `PATH`.
4. Common Cargo build locations relative to the working directory.

Inspect or override the selected executable:

```r
clique_binary()
set_clique_binary("/absolute/path/to/clique")
```

### Other dependencies

`samtools` is required by generated scripts when sorting/indexing is enabled
and by `read_longread_lineage()`.

Most analysis dependencies are optional R packages and are checked only when a
feature needs them:

| Feature | Dependency |
| --- | --- |
| NJ/BIONJ | `ape` |
| In-process parsimony and comparison in examples | `ape`, `phangorn` |
| Tree plotting | `ggplot2`, `ggtree`; annotation rings also use `ggnewscale` |
| Cassiopeia | `reticulate` plus Python `cassiopeia-lineage` |
| IQ-TREE reconstruction | `iqtree2` or `iqtree` executable |
| VINE reconstruction | `vine` executable |
| PHYLIP parsimony | `mix` executable |

External tree wrappers search `PATH` and accept an explicit executable path.
They return installation-oriented errors when a backend is unavailable.

## Quick start

This example covers the normal R-side workflow: create a basic read structure,
generate an align/collapse script, import the consensus BAM, create a character
matrix, and reconstruct trees.

### 1. Create a basic read structure

```r
targets <- data.frame(
  target = c(
    "GACGGCTATACAAGGCATCGCGG",
    "CTCGTCAATACACCTTACGGAGG"
  ),
  type = c("Cas9WT", "Cas9ABE")
)

matches <- fasta_targets_to_yaml(
  fasta_file = "amplicons.fa",
  targets = targets,
  output_file = "read_structure.yaml",
  merge = "ConcatenateBothForward",
  reads = list(
    list(type = "Read1", orientation = "Forward"),
    list(type = "Read2", orientation = "ReverseComplement")
  )
)
```

### 2. Generate the processing script

```r
generate_clique_script(
  read_structure = "read_structure.yaml",
  reference = "amplicons.fa",
  read1 = "sample_R1.fastq.gz",
  read2 = "sample_R2.fastq.gz",
  output_script = "run_sample.sh",
  output_dir = "results",
  sample_name = "sample",
  threads = 8,
  collapse = TRUE
)
```

Review and run it:

```sh
bash run_sample.sh
```

### 3. Import the collapsed BAM

```r
reads <- read_longread_lineage(
  "results/sample.consensus.bam",
  tags = c("e0", "e1", "ce", "pe", "rc", "dc", "ar")
)

reads$rc <- as.integer(reads$rc)
reads$dc <- as.integer(reads$dc)
```

### 4. Build lineage characters

```r
characters <- build_indel_matrix(
  reads,
  cell_col = "e0",
  event_col = "ce",
  missing_threshold = 0.5
)

characters$matrix
characters$state_map
```

### 5. Build several trees

```r
trees <- build_trees(
  characters$matrix,
  methods = c("nj", "parsimony"),
  output_dir = "results/trees"
)

trees$nj
trees$parsimony
```

`build_trees()` isolates failures by method. With its default
`on_error = "warn"`, an unavailable optional backend is reported and returned
as `NULL` while other requested methods continue.

## Configuration

### Generate YAML from FASTA and target sequences

`fasta_targets_to_yaml()` reads a multi-record FASTA, locates targets in each
reference, checks the reverse complement when requested, and writes a basic
read-structure YAML accepted by `clique align` and `clique collapse`.

The helper supports these target annotations:

- `Static`
- `Cas9WT`, `Cas9ABE`, `Cas9CBE`, `Cas9ABECBE`
- `Cas12AWT`, `Cas12ABE`, `Cas12CBE`, `Cas12ABECBE`
- `Cas9Homing`
- `Cas9ABEPalindrome`

The returned object records the target, type, and matched strand for each
reference:

```r
matches <- fasta_targets_to_yaml(
  fasta_file = "amplicons.fa",
  targets = data.frame(
    target = c("ACGT...", "TGCA..."),
    type = c("Cas12AWT", "Cas12ABE")
  ),
  output_file = "read_structure.yaml",
  keep_references_without_targets = FALSE
)

matches[["amplicon_1"]]
```

This helper intentionally generates a starting configuration. Edit the YAML
directly when an experiment needs:

- `umi_configurations` for cell barcodes, lineage barcodes, or UMIs;
- known-tag allowlists and capture correction parameters;
- `target_locations` for repeated target sequences;
- `prime_edits` definitions;
- complex per-reference differences; or
- layout fields that cannot be inferred from FASTA sequence alone.

Prime-edit configurations are supported by the Rust CLI but are not generated
by `fasta_targets_to_yaml()` because they require more than a target sequence.
See the [read-structure YAML documentation](../../rust_cmd/README.md#read-structure-yaml)
for the complete schema.

## Running clique

### Generate a reproducible pipeline

`generate_clique_script()` is the recommended high-level execution interface.
It writes a reviewable Bash script containing absolute input paths and the
selected parameters.

The generated workflow can:

- accept read 1 plus optional read 2, index 1, and index 2 FASTQs;
- use the `WFA`, `Degenerate`, or `Inversion` aligner;
- set threads, minimum read length, and maximum reference-length multiplier;
- cross-check FASTA and YAML reference names;
- sort and index the aligned BAM with `samtools`;
- run capture correction and consensus collapse; and
- run immediately with `execute = TRUE` or be retained for later execution.

For `sample_name = "sample"`, the possible outputs are:

| Output | Meaning |
| --- | --- |
| `sample.aligned.bam` | Direct `clique align` output |
| `sample.sorted.bam` | Coordinate-sorted aligned BAM |
| `sample.sorted.bam.bai` | BAM index |
| `sample.consensus.bam` | Collapsed per-molecule consensus records |

The `reference` FASTA passed to `generate_clique_script()` is used for a
reference-name consistency check and provenance. Clique aligns to the
references embedded in the read-structure YAML.

Set `correct_only = TRUE` to correct captures without producing molecule
consensuses. Collapse requires sorted input, so the generator automatically
enables sorting and indexing when `collapse = TRUE`.

### Invoke CLI subcommands directly

`clique_run()` exposes CLI functionality not represented by the script
generator:

```r
result <- clique_run("collapse", list(
  input_bam_file = "results/sample.sorted.bam",
  output_bam_file = "results/sample.consensus.bam",
  read_structure = "read_structure.yaml",
  threads = 8L,
  maximum_reads_before_downsampling = 100L,
  correct_only = FALSE
))

result$status
result$command
```

Argument-list conventions:

- Names use R-friendly snake case and are converted to kebab-case flags.
- `TRUE` produces a bare flag.
- `FALSE`, `NULL`, and `NA` are omitted.
- Other scalar values are passed as the flag value.
- Nonzero CLI exit status raises an R error with the end of stderr.

The invisible return value contains `status`, `stdout`, `stderr`, and the
resolved command vector.

## Reading BAM output

`read_longread_lineage()` streams an aligned or collapsed BAM through
`samtools view`. It returns one data-frame row per BAM record with:

- `read_name` from SAM `QNAME`;
- `reference` from SAM `RNAME`; and
- one character column per requested auxiliary tag.

```r
reads <- read_longread_lineage(
  "sample.consensus.bam",
  tags = c("e0", "o0", "ce", "pe", "rc", "dc", "ar", "rm")
)
```

Missing tags are `NA`, and all values remain character strings regardless of
their SAM type. Convert numeric tags explicitly before aggregation.

Passing `tags = NULL` discovers tags from the first record and adds `ce`.
Explicit tags are safer when annotations may be present only on later records.
Use `region = "amplicon_name"` to select one reference from an indexed BAM.

### Common clique tags

| Tag | Meaning |
| --- | --- |
| `e0` ... `e9` | Extracted capture on align output; corrected capture on collapse output |
| `o0` ... `o9` | Original representative-read capture before collapse correction |
| `ce` | Chemistry-filtered CRISPR event calls in configured target order |
| `pe` | Prime-edit classification in configured target order |
| `rc` | Reads assigned to the corrected molecule group |
| `dc` | Reads selected for consensus after downsampling |
| `ar` | Comma-separated names of reads selected for the consensus |
| `rm` | Matching-reference-base fraction |

For a normal collapsed record, `length(strsplit(ar, ",")[[1]]) == dc <= rc`.
The `ar` tag lists consensus-selected reads, not every member of a downsampled
group.

See the [BAM output vignette](vignettes/bam-output.Rmd) for capture symbols,
event encodings, collapse provenance, reference-routing diagnostics, and SAM
conventions.

## Character matrices

`build_indel_matrix()` converts the `_`-separated target fields in `ce` into a
cell-by-site integer matrix suitable for lineage solvers.

Its state model is:

- `0`: unedited (`NONE`);
- positive integers: distinct edited alleles assigned independently per site;
- `-1`: missing data by default.

If several records share a `cell_col`, the function selects the majority event
string independently at each target. An `&`-joined compound event at one target
remains one allele. Cells exceeding `missing_threshold` are dropped.

The function returns:

```r
list(
  matrix = ...,    # integer matrix: cells x targets
  state_map = ...  # per-target edit string -> integer state mappings
)
```

Pass explicit `sites` when biological target names matter. Without them,
columns are named `site1`, `site2`, and so on. Avoid combining references with
different target layouts into one matrix; the function warns but cannot infer
whether site positions are comparable.

`indel_distance()` computes a normalized, optionally site-weighted Hamming
distance over sites observed in both cells. A pair with no shared observed site
receives distance `1`.

## Tree reconstruction

All tree builders consume the integer character matrix described above or, for
`tree_nj()`, the `dist` returned by `indel_distance()`. Default results are
`ape::phylo` objects with the original matrix row names restored as tip labels.

### Build multiple methods together

```r
trees <- build_trees(
  characters$matrix,
  methods = c("nj", "parsimony", "iqtree", "vine"),
  args = list(
    nj = list(method = "bionj"),
    iqtree = list(model = "MK+FQ", threads = 8, seed = 1),
    vine = list(threads = 8, nj_only = TRUE)
  ),
  output_dir = "results/trees",
  on_error = "warn"
)
```

One `tree_<method>.nwk` file is written per successful method when
`output_dir` is set.

### Choose a reconstruction method

| Function | Method | Use and behavior | Dependency |
| --- | --- | --- | --- |
| `tree_nj()` | NJ or BIONJ | Fast distance-based exploratory tree | `ape` |
| `tree_parsimony()` | Unordered maximum parsimony | In-process ratchet over categorical multistate characters | `ape`, `phangorn` |
| `tree_mix()` | Camin-Sokal or Wagner parsimony | One-hot expands each positive `(site, allele)` state; can return tied trees | PHYLIP `mix` |
| `tree_iqtree()` | Maximum likelihood | NEXUS `MORPH` alignment with Mk-family model options | IQ-TREE 2 |
| `tree_vine()` | Variational Bayesian CRISPR model | Ultrametric mean tree; can retain posterior samples | VINE |
| `tree_cassiopeia()` | Lineage-specific reconstruction | Greedy, ILP, hybrid, UPGMA, or neighbor joining | Python Cassiopeia |

Important backend details:

- IQ-TREE recodes categorical states independently per site and supports at
  most 32 observed states at one site.
- VINE uses clique's `0`/positive/`-1` CRISPR state convention. Set
  `return_posterior = TRUE` to attach posterior trees as
  `attr(tree, "clique_posterior")`.
- MIX models irreversible gains with `method = "camin-sokal"` or reversible
  gains/losses with `method = "wagner"`.
- Cassiopeia ILP and the default hybrid bottom solver require a working Gurobi
  installation and license.

External wrappers use temporary workspaces by default. Supply `output_prefix`
to retain converted inputs, logs, trees, and other native output. Returned
trees record available provenance in:

- `attr(tree, "clique_backend")`;
- `attr(tree, "clique_command")`; and
- `attr(tree, "clique_files")`.

## Tree visualization

`plot_lineage_tree()` uses one interface for trees from every backend:

```r
cell_metadata <- data.frame(
  cell_barcode = rownames(characters$matrix),
  treatment = sample(c("control", "treated"), nrow(characters$matrix), TRUE),
  edit_count = rowSums(characters$matrix > 0)
)

plot_lineage_tree(
  trees$parsimony,
  annotations = cell_metadata,
  mapping_column = "cell_barcode",
  columns = c("treatment", "edit_count"),
  layout = "circular",
  tip_labels = FALSE,
  title = "Sample lineage tree"
)
```

Categorical and continuous annotation columns are supported. Available layouts
are `circular`, `fan`, and `rectangular`. Use
`clique_categorical_palette()` to inspect or extend the default categorical
colors.

## Examples and further documentation

### Interactive tree-building example

From the package source directory:

```r
devtools::load_all()
source("inst/examples/tree-building-introduction-interactive.R")

result <- run_tree_building_example()
result$comparison
result$pairwise_rf
result$trees$nj
plot_tree_building_result(result)
```

From an installed package:

```r
example_script <- system.file(
  "examples",
  "tree-building-introduction-interactive.R",
  package = "cliqueR"
)
source(example_script)
result <- run_tree_building_example()
```

The companion `inst/examples/tree-building-introduction.R` is a command-line
version with `--tree`, `--matrix`, `--methods`, and `--output` options.

### Package documentation

- [Getting started vignette](vignettes/getting-started.Rmd): configuration
  through lineage reconstruction.
- [BAM output vignette](vignettes/bam-output.Rmd): BAM fields, capture tags,
  collapse provenance, and event encodings.
- [`inst/examples/cliqueR-workflow.Rmd`](inst/examples/cliqueR-workflow.Rmd):
  longer workflow template.
- `help(package = "cliqueR")`: installed help index.
- `?function_name`: argument and return documentation for any exported
  function.

## Function index

| Function | Purpose |
| --- | --- |
| `clique_binary()` | Locate the Rust executable |
| `set_clique_binary()` | Override executable discovery for the R session |
| `clique_run()` | Invoke any clique subcommand with named arguments |
| `fasta_targets_to_yaml()` | Generate a basic read-structure YAML |
| `generate_clique_script()` | Write an align/sort/index/collapse Bash script |
| `read_longread_lineage()` | Import selected clique BAM tags |
| `build_indel_matrix()` | Encode event calls as cell-by-target states |
| `indel_distance()` | Compute missing-aware lineage distances |
| `build_trees()` | Run and write several reconstruction methods |
| `tree_nj()` | Build an NJ or BIONJ tree |
| `tree_parsimony()` | Build an in-process unordered parsimony tree |
| `tree_mix()` | Build PHYLIP MIX parsimony trees |
| `tree_iqtree()` | Build an IQ-TREE maximum-likelihood tree |
| `tree_vine()` | Build a VINE variational Bayesian tree |
| `tree_cassiopeia()` | Build a Cassiopeia lineage tree |
| `plot_lineage_tree()` | Plot a tree with tip-associated annotations |
| `clique_categorical_palette()` | Return the default categorical palette |

## Current limitations

- `align()` and `collapse()` are exported placeholders and currently error.
- `read_10x_lineage()` is reserved and currently errors.
- `rf_distance()`, `triplet_correctness()`, `quartet_distance()`, and
  `ancestor_recall()` are reserved and currently error. The introductory
  tree-building examples compute comparison summaries directly with
  `phangorn`.
- `fasta_targets_to_yaml()` does not infer capture/UMI configurations, repeated
  target locations, or prime-edit specifications.
- `read_longread_lineage()` imports record names, references, and selected
  auxiliary tags; it is not a general-purpose BAM alignment parser.
- The API and output contracts may change while the package remains
  experimental.

## Development

From the repository root:

```r
devtools::document("r_pkg/cliqueR")
devtools::test("r_pkg/cliqueR")
devtools::check("r_pkg/cliqueR", args = "--no-manual")
```
