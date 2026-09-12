"""`sgg analyze` - identify genotype(s) and pathogenicity for input sequence(s)."""

import json
import re
import sys
import time

from sgg.analyzer import (
    FASTAParser,
    SequenceSimilarity,
    analyze_sequence,
    load_all_references,
    load_entry_config,
    unpack_top_match,
)
from sgg.paths import resolve_entry


def _read_fasta_text(path: str) -> str:
    if path == "-":
        return sys.stdin.read()
    with open(path, "r", encoding="utf-8-sig", errors="ignore") as f:
        return f.read()


def _parse_gene_args(gene_args: list[str]) -> dict[str, str]:
    genes = {}
    for spec in gene_args:
        if "=" not in spec:
            raise SystemExit(f"Invalid --gene value '{spec}', expected NAME=path/to.fasta")
        name, path = spec.split("=", 1)
        genes[name.strip()] = path.strip()
    return genes


def _combined_genotype(gene_results: dict, gene_configs: dict) -> str:
    parts = []
    for gene, result in gene_results.items():
        config = gene_configs.get(gene, {})
        pattern = config.get("genotype_pattern", "")
        if not result["genotype_matches"]:
            parts.append("?")
            continue
        top_geno = result["genotype_matches"][0][0]
        if pattern:
            m = re.search(pattern, top_geno)
            parts.append(m.group(0) if m else top_geno)
        else:
            parts.append(top_geno)
    return "".join(parts)


def _best_cleavage(result: dict):
    frames = [
        (result["cleavage_main"], "Main"),
        (result["cleavage_plus_one"], "+1"),
        (result["cleavage_minus_one"], "-1"),
    ]
    for cleavage, label in frames:
        if cleavage["pathogenicity"] not in ("Undetermined", "Not configured"):
            return cleavage, label
    return frames[0]


def _similarity_matrix(sequences: dict[str, str]) -> list[list[float]]:
    seqs = list(sequences.values())
    return [[SequenceSimilarity.pairwise_similarity(s1, s2) for s2 in seqs] for s1 in seqs]


def _print_specimen_text(header: str, result_or_genes, is_multi: bool, gene_configs=None, verbose=False):
    print(f"\n=== {header} ===")
    if is_multi:
        combined = _combined_genotype(result_or_genes, gene_configs or {})
        print(f"Combined genotype: {combined}")
        for gene, result in result_or_genes.items():
            print(f"\n  [{gene}]")
            _print_result_body(result, indent="  ", verbose=verbose)
    else:
        _print_result_body(result_or_genes, indent="", verbose=verbose)


def _print_result_body(result: dict, indent: str, verbose: bool):
    if result.get("error"):
        print(f"{indent}Error: {result['error']}")
        return
    print(f"{indent}Length: {result['sequence_length']} bp")
    matches = result["genotype_matches"]
    if matches:
        top = unpack_top_match(matches[0])
        print(
            f"{indent}Best match: {top['genotype']} "
            f"(avg {top['avg_similarity']}%, best {top['best_score']}% vs {top['best_header']})"
        )
        for genotype, avg, count, best_header, best_score in matches[1:]:
            print(f"{indent}  also: {genotype} (avg {avg}%, n={count})")
    else:
        print(f"{indent}No genotype matches found.")

    if not result.get("pathogenicity_configured"):
        print(f"{indent}Pathogenicity: not configured for this entry")
        return

    if verbose:
        for label, cleavage in (
            ("Main", result["cleavage_main"]),
            ("+1", result["cleavage_plus_one"]),
            ("-1", result["cleavage_minus_one"]),
        ):
            print(f"{indent}Reading frame {label}: {cleavage['pathogenicity']} ({cleavage['motif_type']})")
    else:
        cleavage, label = _best_cleavage(result)
        print(
            f"{indent}Pathogenicity: {cleavage['pathogenicity']} "
            f"(frame {label}, motif {cleavage['motif_type']}, category {cleavage['motif_category']})"
        )


def _result_to_row(header: str, result: dict) -> dict:
    matches = result["genotype_matches"]
    top = unpack_top_match(matches[0]) if matches else None
    cleavage, label = _best_cleavage(result)
    return {
        "Header": header,
        "Length": result["sequence_length"],
        "Best Genotype": top["genotype"] if top else "N/A",
        "Ref Seq Count": top["sample_count"] if top else "N/A",
        "Avg Match Score (%)": top["avg_similarity"] if top else "N/A",
        "Best Match Score (%)": top["best_score"] if top else "N/A",
        "Pathogenicity": cleavage["pathogenicity"],
        "Cleavage Motif": cleavage["motif_type"],
        "Motif Category": cleavage["motif_category"],
        "Reading Frame": label,
    }


def run(args) -> int:
    entry_path = resolve_entry(args.entry)
    entry_config = load_entry_config(str(entry_path))
    is_multi = bool(entry_config.get("multi", False))

    start_time = time.time()

    if is_multi:
        gene_paths = _parse_gene_args(args.gene or [])
        gene_configs = entry_config.get("genes", {})
        known_genes = set(gene_configs)
        missing = known_genes - set(gene_paths)
        if missing:
            raise SystemExit(
                f"Entry '{args.entry}' is multi-gene ({', '.join(sorted(known_genes))}). "
                f"Provide --gene NAME=file.fasta for each: missing {', '.join(sorted(missing))}"
            )

        references, _files, _total, errors = load_all_references(str(entry_path))
        for err in errors:
            print(f"warning: {err}", file=sys.stderr)

        gene_sequences = {
            gene: FASTAParser.parse_text(_read_fasta_text(path))
            for gene, path in gene_paths.items()
        }
        counts = {gene: len(seqs) for gene, seqs in gene_sequences.items()}
        if len(set(counts.values())) > 1:
            details = ", ".join(f"{g}: {c}" for g, c in counts.items())
            raise SystemExit(f"Each gene input must have the same number of sequences. Got: {details}")

        n_specimens = next(iter(counts.values())) if counts else 0
        specimens = []
        for idx in range(n_specimens):
            gene_results = {}
            for gene in gene_paths:
                headers = list(gene_sequences[gene].keys())
                seqs = list(gene_sequences[gene].values())
                header, seq = headers[idx], seqs[idx]
                gene_results[gene] = analyze_sequence(
                    input_fasta=f">{header}\n{seq}",
                    reference_sequences=references[gene],
                    top_matches=args.top,
                    similarity_method=args.method,
                    pathogenicity_config=gene_configs.get(gene, {}),
                )
            specimen_header = next(iter(gene_results.values()))["input_header"]
            specimens.append((specimen_header, gene_results))

        elapsed = time.time() - start_time

        if args.json:
            payload = [
                {"specimen": header, "combined_genotype": _combined_genotype(genes, gene_configs), "genes": genes}
                for header, genes in specimens
            ]
            print(json.dumps(payload, indent=2))
        else:
            for header, gene_results in specimens:
                _print_specimen_text(header, gene_results, is_multi=True, gene_configs=gene_configs, verbose=args.verbose)
            print(f"\n{len(specimens)} specimen(s) analyzed in {elapsed:.2f}s")

        return 0

    # ---- mono-gene entry ----
    references, _files, total_count, errors = load_all_references(str(entry_path))
    for err in errors:
        print(f"warning: {err}", file=sys.stderr)
    if not references:
        raise SystemExit(f"No reference sequences could be loaded for entry '{args.entry}'")

    input_text = _read_fasta_text(args.fasta)
    if not input_text.strip():
        raise SystemExit("No input sequence provided")
    if not input_text.strip().startswith(">"):
        raise SystemExit("Input must be in FASTA format (start with '>')")

    all_sequences = FASTAParser.parse_text(input_text)
    all_results = []
    for header, sequence in all_sequences.items():
        result = analyze_sequence(
            input_fasta=f">{header}\n{sequence}",
            reference_sequences=references,
            top_matches=args.top,
            similarity_method=args.method,
            pathogenicity_config=entry_config,
        )
        all_results.append((header, result))

    elapsed = time.time() - start_time

    matrix = None
    if args.matrix and len(all_sequences) > 1:
        matrix = _similarity_matrix(all_sequences)

    if args.json:
        payload = {
            "entry": args.entry,
            "method": args.method,
            "elapsed_seconds": round(elapsed, 3),
            "results": [{"header": h, **r} for h, r in all_results],
        }
        if matrix is not None:
            payload["similarity_matrix"] = {
                "headers": list(all_sequences.keys()),
                "matrix": matrix,
            }
        print(json.dumps(payload, indent=2, default=str))
    else:
        for header, result in all_results:
            _print_specimen_text(header, result, is_multi=False, verbose=args.verbose)
        print(f"\n{len(all_results)} sequence(s) analyzed in {elapsed:.2f}s")
        if matrix is not None:
            headers = list(all_sequences.keys())
            print("\nSimilarity matrix (pairwise %):")
            for h, row in zip(headers, matrix):
                print(f"  {h[:24]:<24} " + " ".join(f"{v:6.1f}" for v in row))

    if args.matrix_out and matrix is not None:
        import csv

        headers = list(all_sequences.keys())
        with open(args.matrix_out, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([""] + headers)
            for h, row in zip(headers, matrix):
                writer.writerow([h] + row)
        print(f"Similarity matrix written to {args.matrix_out}")

    if args.report:
        _write_report(args, all_results, matrix, all_sequences, elapsed)

    return 0


def _write_report(args, all_results, matrix, all_sequences, elapsed):
    import pandas as pd

    from sgg.report import build_report_html

    export_df = pd.DataFrame([_result_to_row(h, r) for h, r in all_results])

    matrix_fig = None
    if matrix is not None:
        import plotly.graph_objects as go

        headers = list(all_sequences.keys())
        matrix_fig = go.Figure(
            data=go.Heatmap(
                z=matrix,
                x=[h[:10] for h in headers],
                y=[h[:10] for h in headers],
                colorscale="RdYlGn_r",
                zmid=85,
                text=[[f"{val:.2f}%" for val in row] for row in matrix],
                texttemplate="%{text}",
                colorbar={"title": "Similarity %"},
            )
        )

    html = build_report_html(
        all_results=all_results,
        export_df=export_df,
        method=args.method,
        elapsed_time=elapsed,
        entry_name=args.entry,
        matrix_fig=matrix_fig,
    )
    with open(args.report, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"Report written to {args.report}")
