"""`sgg validate` - stratified holdout cross-validation of genotype identification accuracy.

Repeatedly holds out sequences from an entry's own reference dataset and tests the
genotyper against them, to get a mean accuracy instead of a single lucky/unlucky draw.
"""

import json
import os
import random
import re
import sys

import pandas as pd

from sgg.analyzer import GenotypeIdentifier, load_all_references, load_entry_config
from sgg.paths import resolve_entry


def _stratified_holdout(
    references: dict, holdout_size: int, min_remaining_pct: float = 0.20, well_covered_threshold: int = 15
) -> list[str]:
    """Random stratified draw:
    - at least 1 sequence per well-covered genotype (>= well_covered_threshold) is guaranteed
    - at least min_remaining_pct of each genotype's sequences stay in the matching pool
    """
    by_genotype: dict[str, list[str]] = {}
    for h in references:
        parts = h.split("|")
        geno = parts[2] if len(parts) >= 3 else "Unknown"
        by_genotype.setdefault(geno, []).append(h)

    guaranteed = []
    remaining_pool = {geno: list(headers) for geno, headers in by_genotype.items()}

    for geno, headers in by_genotype.items():
        max_holdout = max(0, len(headers) - max(1, int(len(headers) * min_remaining_pct)))
        if len(headers) >= well_covered_threshold and max_holdout > 0:
            chosen = random.choice(remaining_pool[geno])
            guaranteed.append(chosen)
            remaining_pool[geno].remove(chosen)

    if len(guaranteed) >= holdout_size:
        random.shuffle(guaranteed)
        return guaranteed[:holdout_size]

    extra_slots = holdout_size - len(guaranteed)
    candidates = []
    for geno, headers in remaining_pool.items():
        already_taken = len(by_genotype[geno]) - len(headers)
        max_holdout = (
            max(0, len(by_genotype[geno]) - max(1, int(len(by_genotype[geno]) * min_remaining_pct)))
            - already_taken
        )
        if max_holdout > 0:
            candidates.extend(random.sample(headers, min(max_holdout, len(headers))))

    random.shuffle(candidates)
    holdout_headers = guaranteed + candidates[:extra_slots]
    random.shuffle(holdout_headers)
    return holdout_headers


def _apply_pattern(genotype: str, pattern: str) -> str:
    if not pattern:
        return genotype
    m = re.search(pattern, genotype)
    return m.group(0) if m else genotype


def run_validation(
    references: dict, holdout_size: int, n: int, method: str, genotype_pattern: str = "", verbose: bool = False
) -> pd.DataFrame:
    rows = []
    for run_idx in range(n):
        if verbose:
            print(f"Run {run_idx + 1}/{n}...", file=sys.stderr)
        holdout_headers = _stratified_holdout(references, holdout_size)
        holdout_set = set(holdout_headers)
        pool = {h: s for h, s in references.items() if h not in holdout_set}
        identifier = GenotypeIdentifier(pool)

        for header in holdout_headers:
            true_genotype = _apply_pattern(header.split("|")[2], genotype_pattern)
            sequence = references[header]
            matches = identifier.identify(sequence, method=method, top_n=1)
            if not matches:
                continue
            predicted_genotype = _apply_pattern(matches[0][0], genotype_pattern)
            rows.append(
                {
                    "Run": run_idx + 1,
                    "Header": header,
                    "True Genotype": true_genotype,
                    "Predicted Genotype": predicted_genotype,
                    "Match": predicted_genotype == true_genotype,
                    "Score": matches[0][1],
                }
            )
    return pd.DataFrame(rows)


def _summarize(df: pd.DataFrame, top_n: int):
    run_acc = df.groupby("Run")["Match"].mean() * 100
    mean_acc = run_acc.mean()
    std_acc = run_acc.std() if len(run_acc) > 1 else 0.0

    df_miss = df[~df["Match"]]
    if df_miss.empty:
        top_confusion = pd.DataFrame(columns=["True Genotype", "Predicted Genotype", "Count"])
    else:
        top_confusion = (
            df_miss.groupby(["True Genotype", "Predicted Genotype"])
            .size()
            .reset_index(name="Count")
            .sort_values("Count", ascending=False)
            .head(top_n)
        )

    genos = sorted(set(df["True Genotype"]) | set(df["Predicted Genotype"]))
    confusion_matrix = pd.DataFrame(0, index=genos, columns=genos)
    for _, row in df.iterrows():
        confusion_matrix.loc[row["True Genotype"], row["Predicted Genotype"]] += 1

    return run_acc, mean_acc, std_acc, top_confusion, confusion_matrix


def run(args) -> int:
    entry_path = resolve_entry(args.entry)
    entry_config = load_entry_config(str(entry_path))
    is_multi = bool(entry_config.get("multi", False))
    db_references, _files, _total, errors = load_all_references(str(entry_path))
    for err in errors:
        print(f"warning: {err}", file=sys.stderr)

    if is_multi:
        gene_names = sorted(db_references)
        selected_gene = args.gene or gene_names[0]
        if selected_gene not in gene_names:
            raise SystemExit(f"Unknown gene '{selected_gene}'. Available: {', '.join(gene_names)}")
        references = db_references[selected_gene]
        gene_cfg = entry_config.get("genes", {}).get(selected_gene, {})
        genotype_pattern = gene_cfg.get("genotype_pattern", "")
    else:
        selected_gene = None
        references = db_references
        genotype_pattern = ""

    gene_total = len(references)
    if gene_total < 10:
        raise SystemExit(
            f"Only {gene_total} reference sequences in this entry - too few for a meaningful holdout test."
        )

    holdout_size = args.holdout_size or max(1, min(100, gene_total // 3))
    max_holdout = max(1, gene_total // 2)
    if holdout_size > max_holdout:
        raise SystemExit(f"--holdout-size {holdout_size} is too large for {gene_total} sequences (max {max_holdout}).")

    if args.seed is not None:
        random.seed(args.seed)

    print(
        f"Validating '{args.entry}'" + (f" / {selected_gene}" if selected_gene else "")
        + f" - {gene_total} references, holdout={holdout_size}, runs={args.runs}, method={args.method}",
        file=sys.stderr,
    )
    df = run_validation(references, holdout_size, args.runs, args.method, genotype_pattern, verbose=True)

    if df.empty:
        raise SystemExit("No scorable predictions were produced.")

    run_acc, mean_acc, std_acc, top_confusion, confusion_matrix = _summarize(df, args.top_confusions)

    if args.csv_out:
        os.makedirs(args.csv_out, exist_ok=True)
        df.to_csv(os.path.join(args.csv_out, "validation_predictions.csv"), index=False)
        confusion_matrix.to_csv(os.path.join(args.csv_out, "validation_confusion_matrix.csv"))
        pd.DataFrame({"Run": run_acc.index, "Accuracy (%)": run_acc.values}).to_csv(
            os.path.join(args.csv_out, "validation_accuracy_per_run.csv"), index=False
        )
        print(f"Wrote predictions, confusion matrix, and per-run accuracy to {args.csv_out}", file=sys.stderr)

    if args.json:
        payload = {
            "entry": args.entry,
            "gene": selected_gene,
            "reference_count": gene_total,
            "holdout_size": holdout_size,
            "runs": args.runs,
            "method": args.method,
            "mean_accuracy": round(float(mean_acc), 2),
            "std_accuracy": round(float(std_acc), 2),
            "total_predictions": len(df),
            "accuracy_per_run": [{"run": int(r), "accuracy": round(float(a), 2)} for r, a in run_acc.items()],
            "top_confusions": top_confusion.to_dict(orient="records"),
        }
        print(json.dumps(payload, indent=2))
        return 0

    print(f"\nMean accuracy: {mean_acc:.1f}% (+/- {std_acc:.1f})  |  {len(run_acc)} runs, {len(df)} predictions total")
    print("\nAccuracy per run:")
    for r, a in run_acc.items():
        print(f"  Run {r}: {a:.1f}%")

    if top_confusion.empty:
        print("\nNo misclassification across any run.")
    else:
        print(f"\nTop {len(top_confusion)} confusions:")
        print(top_confusion.to_string(index=False))

    if not args.csv_out:
        print("\nTip: pass --csv-out DIR to save the full confusion matrix and raw predictions.")

    return 0
