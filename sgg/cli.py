"""sgg - Soupirr's Global Genotyper CLI."""

import argparse
import sys

from sgg import __version__
from sgg.commands import analyze as analyze_cmd
from sgg.commands import entries as entries_cmd
from sgg.commands import stats as stats_cmd
from sgg.commands import tree as tree_cmd
from sgg.commands import validate as validate_cmd


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sgg",
        description="Identify pathogen genotypes and predict pathogenicity from nucleotide sequences.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p_list = sub.add_parser("list-entries", help="List available reference entries")
    p_list.add_argument("--json", action="store_true", help="Output JSON")
    p_list.set_defaults(func=entries_cmd.list_entries)

    p_analyze = sub.add_parser("analyze", help="Analyze sequence(s) against a reference entry")
    p_analyze.add_argument("fasta", nargs="?", default="-", help="FASTA file (default: stdin). Ignored for multi-gene entries; use --gene instead.")
    p_analyze.add_argument("--entry", required=True, help="Reference entry name (see `sgg list-entries`)")
    p_analyze.add_argument(
        "--gene",
        action="append",
        metavar="NAME=file.fasta",
        help="For multi-gene entries: one per gene, repeatable",
    )
    p_analyze.add_argument("--method", choices=["hamming", "pairwise"], default="pairwise")
    p_analyze.add_argument("--top", type=int, default=3, help="Number of top genotype matches to show (default: 3)")
    p_analyze.add_argument("--matrix", action="store_true", help="Compute a pairwise similarity matrix across the input sequences")
    p_analyze.add_argument("--matrix-out", metavar="FILE.csv", help="Write the similarity matrix to a CSV file")
    p_analyze.add_argument("--verbose", action="store_true", help="Show all three reading frames instead of just the best match")
    p_analyze.add_argument("--json", action="store_true", help="Output machine-readable JSON instead of text")
    p_analyze.add_argument("--report", metavar="FILE.html", help="Write an HTML report to FILE")
    p_analyze.set_defaults(func=analyze_cmd.run)

    p_tree = sub.add_parser("tree", help="Build phylogenetic tree(s) for input sequence(s) against a reference entry")
    p_tree.add_argument("fasta", nargs="?", default="-", help="FASTA file with query sequence(s) (default: stdin)")
    p_tree.add_argument("--entry", required=True, help="Reference entry name (see `sgg list-entries`)")
    p_tree.add_argument("--gene", metavar="NAME", help="For multi-gene entries: which gene's reference set to use")
    p_tree.add_argument("--type", choices=["per-query", "combined"], default="per-query", help="per-query: one tree per input sequence. combined: all queries in a single tree (default: per-query)")
    p_tree.add_argument("--neighbours", type=int, default=20, metavar="N", help="Closest references to include, per-query mode (default: 20)")
    p_tree.add_argument("--per-query-neighbours", type=int, default=5, metavar="N", help="Closest references to pool per query, combined mode (default: 5)")
    p_tree.add_argument("--mode", choices=["cladogram", "phylogram"], default="cladogram", help="Cladogram: uniform branch lengths. Phylogram: real evolutionary distances.")
    p_tree.add_argument("--tree-method", choices=["fasttree", "iqtree"], default="fasttree", help="FastTree: fast, approximate. IQ-TREE2: full ML + ModelFinder + bootstrap, slower (default: fasttree)")
    p_tree.add_argument("--out-dir", default=".", metavar="DIR", help="Directory to write .nwk (and --plot .html) files to (default: current directory)")
    p_tree.add_argument("--plot", action="store_true", help="Also render a self-contained HTML plot next to each .nwk file")
    p_tree.set_defaults(func=tree_cmd.run)

    p_stats = sub.add_parser("stats", help="Show database statistics for a reference entry")
    p_stats.add_argument("--entry", required=True, help="Reference entry name (see `sgg list-entries`)")
    p_stats.add_argument("--gene", metavar="NAME", help="For multi-gene entries: which gene to show stats for (default: first gene)")
    p_stats.add_argument("--top-countries", type=int, default=15, metavar="N")
    p_stats.add_argument("--generate-pathogenicity", action="store_true", help="Compute and cache pathogenicity data for every reference sequence if not already generated")
    p_stats.add_argument("--json", action="store_true", help="Output machine-readable JSON instead of text")
    p_stats.set_defaults(func=stats_cmd.run)

    p_validate = sub.add_parser("validate", help="Cross-validate genotype identification accuracy against an entry's own reference dataset")
    p_validate.add_argument("--entry", required=True, help="Reference entry name (see `sgg list-entries`)")
    p_validate.add_argument("--gene", metavar="NAME", help="For multi-gene entries: which gene to validate (default: first gene)")
    p_validate.add_argument("--holdout-size", type=int, metavar="N", help="Sequences held out per run (default: 1/3 of the dataset, capped at 100)")
    p_validate.add_argument("--runs", type=int, default=5, metavar="N", help="Number of holdout runs (default: 5)")
    p_validate.add_argument("--method", choices=["hamming", "pairwise"], default="pairwise", help="Pairwise is much slower but more accurate; Hamming assumes pre-aligned sequences (default: pairwise)")
    p_validate.add_argument("--top-confusions", type=int, default=20, metavar="N")
    p_validate.add_argument("--csv-out", metavar="DIR", help="Write full predictions, confusion matrix, and per-run accuracy CSVs to DIR")
    p_validate.add_argument("--seed", type=int, help="Random seed for reproducible holdout draws")
    p_validate.add_argument("--json", action="store_true", help="Output machine-readable JSON instead of text")
    p_validate.set_defaults(func=validate_cmd.run)

    p_add = sub.add_parser("add-entry", help="Add a new single-gene reference dataset")
    p_add.add_argument("name", help="Entry name")
    p_add.add_argument("fasta", nargs="+", help="One or more reference FASTA files")
    p_add.add_argument("--cleavage-start", type=int, default=0, metavar="N", help="0-indexed nucleotide position where the virulence motif area starts")
    p_add.add_argument("--motifs", metavar="FILE.csv", help="CSV with columns: motif,label,type (type is 'virulent' or 'avirulent')")
    p_add.add_argument("--migrate", action="store_true", help="Run NCBI Virus header migration on the input FASTA(s) before storing")
    p_add.set_defaults(func=entries_cmd.add_entry)

    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args) or 0
    except KeyboardInterrupt:
        print("\nInterrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main())
