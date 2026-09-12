"""Reference entry management: listing and adding entries under the user data dir."""

import json
import os
import shutil

from sgg.migration import migrate_fasta_text
from sgg.paths import USER_SEQ_DIR, list_entries, resolve_entry

__all__ = ["list_entries", "resolve_entry", "add_mono_entry"]


def add_mono_entry(
    name: str,
    fasta_paths: list[str],
    cleavage_start: int = 0,
    motif_path: str | None = None,
    migrate: bool = False,
) -> dict:
    """Add a single-gene reference entry under the user data dir.

    Returns a report dict: {"entry_path": str, "files": [(filename, mig_stats|None), ...]}.
    """
    name = name.strip()
    if not name:
        raise SystemExit("Entry name must not be empty")
    if not fasta_paths:
        raise SystemExit("At least one FASTA file is required")

    entry_path = USER_SEQ_DIR / name
    entry_path.mkdir(parents=True, exist_ok=True)

    file_reports = []
    for fasta_path in fasta_paths:
        with open(fasta_path, "r", encoding="utf-8-sig", errors="ignore") as f:
            raw_text = f.read()

        dest_name = os.path.basename(fasta_path)
        dest_path = entry_path / dest_name

        if migrate:
            migrated_text, mig_stats = migrate_fasta_text(raw_text)
            dest_path.write_text(migrated_text, encoding="utf-8")
            file_reports.append((dest_name, mig_stats))
        else:
            dest_path.write_text(raw_text, encoding="utf-8")
            file_reports.append((dest_name, None))

    if cleavage_start:
        config_path = entry_path / "_config.json"
        config_path.write_text(json.dumps({"cleavage_start": int(cleavage_start)}, indent=2))

    if motif_path:
        motif_dest = entry_path / f"{name}_motifs.csv"
        shutil.copy(motif_path, motif_dest)

    return {"entry_path": str(entry_path), "files": file_reports}
