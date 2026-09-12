"""`sgg stats` - database statistics for a reference entry (genotype/host/year/country/pathogenicity)."""

import json
import os
import re
import sys

import pandas as pd

from sgg.analyzer import CleavageSiteAnalyzer, load_all_references, load_entry_config
from sgg.geography import build_location_dataframe
from sgg.paths import USER_DATA_DIR, resolve_entry

COVERAGE_GOOD = 15
COVERAGE_LOW = 5


def _print_table(df: pd.DataFrame, title: str | None = None) -> None:
    if title:
        print(f"\n{title}")
    if df.empty:
        print("  (no data)")
        return
    print(df.to_string(index=False))


def _pathogenicity_csv_path(entry_path, entry_name: str):
    """Where a generated pathogenicity CSV lives: next to the entry, or a user cache fallback
    when the entry directory isn't writable (e.g. a bundled, root-owned AUR-installed entry)."""
    primary = os.path.join(entry_path, f"{entry_name}_pathogenicity.csv")
    fallback = USER_DATA_DIR / "pathogenicity_cache" / f"{entry_name}_pathogenicity.csv"
    return primary, fallback


def _generate_pathogenicity(flat_refs: dict, gene_cfg: dict, entry_path, entry_name: str) -> str:
    rows = []
    for header, sequence in flat_refs.items():
        parts = header.split("|")
        genotype = parts[2] if len(parts) > 2 else "Unknown"
        cleavage, _, _ = CleavageSiteAnalyzer.analyze(
            sequence,
            cleavage_start=gene_cfg["cleavage_start"],
            motifs_by_type=gene_cfg.get("motifs_by_type"),
        )
        rows.append(
            {
                "Header": header,
                "Best Genotype": genotype,
                "Pathogenicity": cleavage["pathogenicity"],
                "Cleavage Motif": cleavage["motif_type"] or "N/A",
                "Motif Category": cleavage["motif_category"] or "N/A",
            }
        )
    df = pd.DataFrame(rows)

    primary, fallback = _pathogenicity_csv_path(entry_path, entry_name)
    try:
        df.to_csv(primary, index=False)
        return primary
    except OSError:
        fallback.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(fallback, index=False)
        return str(fallback)


def _load_pathogenicity(entry_path, entry_name: str):
    primary, fallback = _pathogenicity_csv_path(entry_path, entry_name)
    if os.path.exists(primary):
        return pd.read_csv(primary)
    if fallback.exists():
        return pd.read_csv(fallback)
    return None


def run(args) -> int:
    entry_path = resolve_entry(args.entry)
    entry_name = os.path.basename(str(entry_path).rstrip("/"))
    entry_config = load_entry_config(str(entry_path))
    is_multi = bool(entry_config.get("multi", False))

    db_references, _files, _total, errors = load_all_references(str(entry_path))
    for err in errors:
        print(f"warning: {err}", file=sys.stderr)

    combined = None
    if is_multi:
        gene_names = sorted(db_references)
        selected_gene = args.gene or gene_names[0]
        if selected_gene not in gene_names:
            raise SystemExit(f"Unknown gene '{selected_gene}'. Available: {', '.join(gene_names)}")
        flat_refs = db_references[selected_gene]
        gene_path = os.path.join(entry_path, selected_gene)
        gene_cfg = entry_config.get("genes", {}).get(selected_gene, {})

        total_sequences = sum(len(refs) for refs in db_references.values())
        total_bp = sum(len(seq) for refs in db_references.values() for seq in refs.values())
        unique_genotypes = len(
            {h.split("|")[2] for refs in db_references.values() for h in refs if len(h.split("|")) >= 3}
        )
        combined = {"total_sequences": total_sequences, "total_bp": total_bp, "unique_genotypes": unique_genotypes}
    else:
        selected_gene = None
        flat_refs = db_references
        gene_path = str(entry_path)
        gene_cfg = entry_config or {}

    if not flat_refs:
        raise SystemExit(f"No reference sequences could be loaded for entry '{args.entry}'")

    db_total_count = len(flat_refs)
    total_bp = sum(len(seq) for seq in flat_refs.values())
    unique_genotypes = len({h.split("|")[2] for h in flat_refs if len(h.split("|")) >= 3})

    geno_pattern = gene_cfg.get("genotype_pattern", "") if is_multi else ""

    genotype_counts: dict[str, int] = {}
    host_counts: dict[str, int] = {}
    year_counts: dict[str, int] = {}

    for header in flat_refs:
        parts = header.split("|")
        if len(parts) < 7:
            continue

        raw_geno = parts[2]
        if geno_pattern:
            m = re.search(geno_pattern, raw_geno)
            geno = m.group(0) if m else raw_geno
        else:
            geno = raw_geno
        genotype_counts[geno] = genotype_counts.get(geno, 0) + 1

        host = parts[3].replace("_", " ")
        host = host if host not in ("UNKNOWN", "?", "") else "Unspecified"
        host_counts[host] = host_counts.get(host, 0) + 1

        year = parts[6]
        if year.isdigit() and len(year) == 4:
            year_counts[year] = year_counts.get(year, 0) + 1

    df_genotypes = pd.DataFrame(
        sorted(genotype_counts.items(), key=lambda x: x[1], reverse=True),
        columns=["Genotype", "Sequences"],
    )

    df_years = pd.DataFrame(
        sorted((int(y), c) for y, c in year_counts.items()),
        columns=["Year", "Sequences"],
    )

    host_items = sorted(host_counts.items(), key=lambda x: x[1], reverse=True)
    host_items = [(k, v) for k, v in host_items if k == "Unspecified"] + [
        (k, v) for k, v in host_items if k != "Unspecified"
    ]
    df_hosts = pd.DataFrame(host_items, columns=["Host", "Sequences"])

    df_loc = build_location_dataframe(gene_path)
    df_loc = df_loc[df_loc["label"] != "Unknown"].copy()
    df_loc["country"] = df_loc["label"].apply(lambda x: x.split(",")[-1].strip())
    country_counts = df_loc["country"].value_counts()
    df_countries = pd.DataFrame(
        sorted(country_counts.items(), key=lambda x: x[1], reverse=True)[: args.top_countries],
        columns=["Country", "Sequences"],
    )

    years_with_data = [int(y) for y in year_counts if year_counts[y] > 0]
    year_range = f"{min(years_with_data)}-{max(years_with_data)}" if years_with_data else "N/A"

    df_health = df_genotypes.copy()
    df_health["Status"] = df_health["Sequences"].apply(
        lambda n: "Good" if n >= COVERAGE_GOOD else ("Low" if n >= COVERAGE_LOW else "Critical")
    )
    df_health["% of Total"] = (df_health["Sequences"] / db_total_count * 100).round(1)

    # ---- pathogenicity ----
    df_patho = _load_pathogenicity(entry_path, entry_name)
    patho_generated_path = None
    if df_patho is None and "cleavage_start" in gene_cfg and args.generate_pathogenicity:
        patho_generated_path = _generate_pathogenicity(flat_refs, gene_cfg, entry_path, entry_name)
        df_patho = _load_pathogenicity(entry_path, entry_name)

    patho_summary = None
    if df_patho is not None:
        patho_counts = (
            df_patho.groupby(["Best Genotype", "Pathogenicity"]).size().reset_index(name="Count")
        )
        motif_counts = (
            df_patho.groupby(["Cleavage Motif", "Motif Category"]).size().reset_index(name="Count")
        )
        motif_counts = motif_counts[motif_counts["Cleavage Motif"] != "N/A"].sort_values(
            "Count", ascending=False
        )
        overall = df_patho["Pathogenicity"].value_counts()
        df_overall = pd.DataFrame(
            {"Pathogenicity": overall.index, "Count": overall.values, "% of Total": (overall.values / overall.values.sum() * 100).round(1)}
        )
        patho_summary = {"by_genotype": patho_counts, "motifs": motif_counts, "overall": df_overall}

    if args.json:
        payload = {
            "entry": args.entry,
            "gene": selected_gene,
            "combined": combined,
            "total_sequences": db_total_count,
            "total_bp": total_bp,
            "average_length": round(total_bp / db_total_count, 1),
            "unique_genotypes": unique_genotypes,
            "unique_sub_genotypes": len(genotype_counts),
            "year_range": year_range,
            "genotype_counts": df_genotypes.to_dict(orient="records"),
            "year_counts": df_years.to_dict(orient="records"),
            "host_counts": df_hosts.to_dict(orient="records"),
            "top_countries": df_countries.to_dict(orient="records"),
            "coverage_health": df_health.to_dict(orient="records"),
            "pathogenicity": (
                {
                    "by_genotype": patho_summary["by_genotype"].to_dict(orient="records"),
                    "motifs": patho_summary["motifs"].to_dict(orient="records"),
                    "overall": patho_summary["overall"].to_dict(orient="records"),
                }
                if patho_summary
                else None
            ),
        }
        print(json.dumps(payload, indent=2, default=str))
        return 0

    print(f"=== {args.entry}" + (f" / {selected_gene}" if selected_gene else "") + " ===")
    if combined:
        print(
            f"Combined (all genes): {combined['total_sequences']:,} sequences, "
            f"{combined['total_bp']:,} bp, {combined['unique_genotypes']:,} unique genotypes"
        )
    print(f"Sequences: {db_total_count:,}  |  Total bp: {total_bp:,}  |  Avg length: {total_bp / db_total_count:.1f} bp")
    print(f"Unique genotypes: {unique_genotypes}  |  Unique sub-genotypes: {len(genotype_counts)}  |  Year range: {year_range}")

    _print_table(df_genotypes, f"Sequences per {'sub-genotype' if is_multi else 'genotype'}")
    _print_table(df_years, "Sequences by year")
    _print_table(df_hosts, "Host distribution")
    _print_table(df_countries, f"Top {args.top_countries} countries")

    good = int((df_health["Status"] == "Good").sum())
    low = int((df_health["Status"] == "Low").sum())
    critical = int((df_health["Status"] == "Critical").sum())
    print(
        f"\nCoverage health (>= {COVERAGE_GOOD} good, >= {COVERAGE_LOW} low, else critical): "
        f"{good} good, {low} low, {critical} critical"
    )
    if critical:
        _print_table(df_health[df_health["Status"] == "Critical"], "Genotypes needing urgent attention")

    if patho_summary is not None:
        if patho_generated_path:
            print(f"\nGenerated pathogenicity data: {patho_generated_path}")
        _print_table(patho_summary["overall"], "Pathogenicity distribution")
        _print_table(patho_summary["by_genotype"], "Pathogenicity by genotype")
        _print_table(patho_summary["motifs"], "Motif distribution")
    elif "cleavage_start" in gene_cfg:
        print("\nNo pathogenicity data generated yet. Re-run with --generate-pathogenicity.")
    else:
        print("\nNo pathogenicity configuration for this entry.")

    return 0
