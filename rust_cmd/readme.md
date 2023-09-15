# Clique command line aligner

This directory contains the Clique alignment and collapsing tool. The goal of this tool is to take amplication sequencing data, from either Illumina or long-read technologies like Nanopore, and 
produce consensus sequences that are aligned to the genome. Often these amplicons have internal structure, such as unique molecular identifiers (UMIs), sequences that mark individual integration
locations within the genome ('static IDs'), or cell identifiers that come from many single-cell sequencing experiments. You provide this layout as a YAML file (detailed below) which Clique uses
to collapse down reads to a consensus sequence, accounting for errors or other issues 
