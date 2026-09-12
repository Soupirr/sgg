# sgg - Soupirr's Global Genotyper

CLI to identify pathogen genotypes and predict pathogenicity from nucleotide sequences.

## Usage

```
sgg list-entries
sgg analyze seq.fasta --entry "Newcastle Disease Virus F gene" --method pairwise
sgg add-entry my-virus refs.fasta --cleavage-start 333 --motifs motifs.csv
```

## Requirements

`mafft`, `fasttree`, and `iqtree` (system packages) are required for the `tree` command.
