"""`sgg list-entries` and `sgg add-entry`."""

import json

from sgg.analyzer import load_entry_config
from sgg.entries import add_mono_entry
from sgg.paths import list_entries as _list_entries


def list_entries(args) -> int:
    entries = _list_entries()
    if not entries:
        print("No entries found.")
        return 0

    if args.json:
        payload = []
        for name, path in sorted(entries.items()):
            config = load_entry_config(str(path))
            payload.append({"name": name, "path": str(path), "multi": bool(config.get("multi"))})
        print(json.dumps(payload, indent=2))
        return 0

    for name, path in sorted(entries.items()):
        config = load_entry_config(str(path))
        kind = "multi-gene" if config.get("multi") else "mono-gene"
        print(f"{name}  [{kind}]  ({path})")
    return 0


def add_entry(args) -> int:
    report = add_mono_entry(
        name=args.name,
        fasta_paths=args.fasta,
        cleavage_start=args.cleavage_start,
        motif_path=args.motifs,
        migrate=args.migrate,
    )
    print(f"Entry '{args.name}' created at {report['entry_path']}")
    for filename, mig_stats in report["files"]:
        if mig_stats is None:
            print(f"  - {filename}")
        else:
            print(
                f"  - {filename}: {mig_stats['converted']}/{mig_stats['input']} sequences kept "
                f"({mig_stats['duplicates']} duplicates, {mig_stats['no_genotype']} no genotype, "
                f"{mig_stats['malformed']} malformed)"
            )
    return 0
