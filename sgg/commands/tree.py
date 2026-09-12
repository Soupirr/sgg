"""`sgg tree` - build phylogenetic tree(s) for input sequence(s) against a reference entry."""

import os
import re
import sys
import tempfile

from sgg.analyzer import (
    FASTAParser,
    align_sequences_mafft,
    build_tree_fasttree,
    build_tree_iqtree2,
    clean_sequence,
    find_closest_neighbours,
    load_all_references,
    load_entry_config,
    tree_to_newick,
    write_temp_fasta,
)
from sgg.paths import resolve_entry
from sgg.treeplot import build_tree_figure


def _read_fasta_text(path: str) -> str:
    if path == "-":
        return sys.stdin.read()
    with open(path, "r", encoding="utf-8-sig", errors="ignore") as f:
        return f.read()


def _safe_filename(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)[:80] or "tree"


def _build_tree(aln_file: str, method: str):
    return build_tree_iqtree2(aln_file) if method == "iqtree" else build_tree_fasttree(aln_file)


def _emit(tree, title: str, filename_base: str, args, multi_query: bool) -> None:
    nwk_path = os.path.join(args.out_dir, f"{filename_base}.nwk")
    with open(nwk_path, "w") as f:
        f.write(tree_to_newick(tree))
    print(f"{title}: {nwk_path}")

    if args.plot:
        fig = build_tree_figure(tree, title, mode=args.mode, multi_query=multi_query)
        html_path = os.path.join(args.out_dir, f"{filename_base}.html")
        fig.write_html(html_path, include_plotlyjs=True)
        print(f"  plot: {html_path}")


def run(args) -> int:
    entry_path = resolve_entry(args.entry)
    entry_config = load_entry_config(str(entry_path))
    is_multi = bool(entry_config.get("multi", False))
    references, _files, _total, errors = load_all_references(str(entry_path))
    for err in errors:
        print(f"warning: {err}", file=sys.stderr)

    if is_multi:
        if not args.gene:
            raise SystemExit(
                f"Entry '{args.entry}' is multi-gene ({', '.join(sorted(references))}). "
                f"Pass --gene NAME to select which gene's reference set to build a tree against."
            )
        if args.gene not in references:
            raise SystemExit(
                f"Unknown gene '{args.gene}' for entry '{args.entry}'. "
                f"Available: {', '.join(sorted(references))}"
            )
        active_refs = references[args.gene]
    else:
        active_refs = references

    if not active_refs:
        raise SystemExit(f"No reference sequences could be loaded for entry '{args.entry}'")

    input_text = _read_fasta_text(args.fasta)
    if not input_text.strip():
        raise SystemExit("No input sequence provided")
    all_sequences = FASTAParser.parse_text(input_text)
    if not all_sequences:
        raise SystemExit("Could not parse any sequence from input")

    os.makedirs(args.out_dir, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="sgg-tree-") as tmp_dir:
        if args.type == "combined":
            pool = {}
            for qseq in all_sequences.values():
                for h, s, _ in find_closest_neighbours(qseq, active_refs, n=args.per_query_neighbours):
                    pool[h] = s

            tmp_fasta = os.path.join(tmp_dir, "combined_input.fasta")
            tmp_align = os.path.join(tmp_dir, "combined_aligned.fasta")
            with open(tmp_fasta, "w") as f:
                f.writelines(
                    f">QUERY_{qh}\n{clean_sequence(qs)}\n" for qh, qs in all_sequences.items()
                )
                for h, s in pool.items():
                    f.write(f">{h}\n{clean_sequence(s)}\n")

            align_sequences_mafft(tmp_fasta, tmp_align)
            tree = _build_tree(tmp_align, args.tree_method)
            if tree is None:
                raise SystemExit("Could not build the combined tree.")
            _emit(tree, "Combined tree - All queries", "combined", args, multi_query=True)

        else:
            any_built = False
            for header, sequence in all_sequences.items():
                neighbours = find_closest_neighbours(sequence, active_refs, n=args.neighbours)
                tmp_fasta = os.path.join(tmp_dir, "tmp_input.fasta")
                tmp_aligned = os.path.join(tmp_dir, "tmp_aligned.fasta")
                write_temp_fasta(header, sequence, neighbours, tmp_fasta)
                align_sequences_mafft(tmp_fasta, tmp_aligned)
                tree = _build_tree(tmp_aligned, args.tree_method)
                if tree is None:
                    print(f"warning: could not build tree for {header}", file=sys.stderr)
                    continue
                _emit(tree, header, _safe_filename(header), args, multi_query=False)
                any_built = True

            if not any_built:
                raise SystemExit("No tree could be built for any input sequence.")

    return 0
