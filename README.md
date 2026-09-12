# sgg - Soupirr's Global Genotyper

Identify pathogen genotypes and predict pathogenicity from nucleotide sequences, from the terminal.

CLI port of the original [Streamlit app](https://github.com/Soupirr/soupirr-global-genotyper).

## Installation

```
yay -S sgg
```

`mafft`, `fasttree`, and `iqtree` are pulled in automatically and are only needed for the `tree`
command - everything else works without them.

## You can run `sgg` from any directory

`sgg` doesn't care where your shell's current directory is. Reference data isn't read relative to
your cwd - it's resolved from two fixed locations:

- **Bundled entries** (Newcastle Disease Virus, Avian Influenza, Bluetongue, Epizootic Hemorrhagic
  Disease Virus) ship inside the package itself, under `/usr/lib/pythonX.Y/site-packages/sgg/data/`.
- **Entries you add** with `sgg add-entry` are written under `$XDG_DATA_HOME/sgg`
  (usually `~/.local/share/sgg/sequences/`), so they survive package upgrades/reinstalls.

The only paths that matter are the ones you pass explicitly - your input FASTA file, and, if you
use them, `--out-dir` / `--csv-out` / `--report` for where results get written. Those are resolved
relative to wherever you happen to run the command from, same as any other CLI tool.

See what's available with:

```
sgg list-entries
```

```
Avian Influenza Alpha Virus HA, NA genes  [multi-gene]  (/usr/lib/python3.14/site-packages/sgg/data/sequences/Avian Influenza Alpha Virus HA, NA genes)
Bluetongue Virus VP2 gene  [mono-gene]  (/usr/lib/python3.14/site-packages/sgg/data/sequences/Bluetongue Virus VP2 gene)
Epizootic Hemmorhagic Disease Virus VP2 gene  [mono-gene]  (/usr/lib/python3.14/site-packages/sgg/data/sequences/Epizootic Hemmorhagic Disease Virus VP2 gene)
Newcastle Disease Virus F gene  [mono-gene]  (/usr/lib/python3.14/site-packages/sgg/data/sequences/Newcastle Disease Virus F gene)
```

A **mono-gene** entry has one reference dataset. A **multi-gene** entry (like Avian Influenza's
HA/NA genes) has one dataset per gene, and commands that need a specific gene take a `--gene` flag.

## Commands

### `sgg analyze` - identify a genotype

```
sgg analyze seq.fasta --entry "Newcastle Disease Virus F gene" --method pairwise
```

Reads from stdin if you omit the file: `cat seq.fasta | sgg analyze --entry "..."`. A FASTA file
with multiple sequences analyzes each one in turn.

For multi-gene entries, skip the positional file and pass one `--gene NAME=file.fasta` per gene
instead (all gene files need the same number of sequences, matched by order):

```
sgg analyze --entry "Avian Influenza Alpha Virus HA, NA genes" \
  --gene HA=ha.fasta --gene NA=na.fasta
```

Useful flags: `--method {hamming,pairwise}` (pairwise is slower but handles indels), `--top N`
matches to show, `--matrix` to also compare the input sequences against each other,
`--json` for machine-readable output, `--report out.html` for a self-contained HTML report.

### `sgg tree` - build a phylogenetic tree

```
sgg tree seq.fasta --entry "Newcastle Disease Virus F gene" --plot --out-dir ./trees
```

`--type per-query` (default) builds one tree per input sequence against its closest references;
`--type combined` pools all queries into a single tree. `--tree-method {fasttree,iqtree}` picks the
tree-building tool; `--plot` also writes a self-contained HTML visualization next to each `.nwk`
file. For multi-gene entries, pick a gene with `--gene NAME`.

### `sgg stats` - database statistics

```
sgg stats --entry "Newcastle Disease Virus F gene"
```

Prints sequence counts, genotype/host/year/country breakdowns, and a coverage-health check
(genotypes with too few reference sequences to be reliable). If the entry has a pathogenicity
config but no pathogenicity data generated yet, add `--generate-pathogenicity` to compute and
cache it (works even on the read-only bundled entries - the cache falls back to
`$XDG_DATA_HOME/sgg/pathogenicity_cache/`).

### `sgg validate` - cross-validate accuracy

```
sgg validate --entry "Newcastle Disease Virus F gene" --runs 10 --method hamming
```

Repeatedly holds out sequences from the entry's own reference dataset and tests the identifier
against them, to get a mean accuracy instead of one lucky/unlucky draw. `--csv-out DIR` saves the
full predictions, confusion matrix, and per-run accuracy as CSVs.

### `sgg add-entry` - add your own reference dataset

```
sgg add-entry my-virus refs.fasta --cleavage-start 333 --motifs motifs.csv
```

Only single-gene entries are supported via the CLI for now. `--motifs` is a CSV with columns
`motif,label,type` (`type` is `virulent` or `avirulent`) for pathogenicity prediction. `--migrate`
runs NCBI Virus header migration on the input first if it isn't already in
`virus|accession|genotype|host|country|region|year` format.

Every command supports `--help` for the full flag list, and most support `--json` for scripting.

This CLI was ported from the [Streamlit version](https://github.com/Soupirr/soupirr-global-genotyper) with Claude Code.
