"""Data locations.

Two tiers of reference data:
- bundled: shipped with the package (read-only, e.g. under /usr/lib/pythonX/site-packages/sgg/data
  once installed via pacman/AUR).
- user: entries added locally via `sgg add-entry`, stored under $XDG_DATA_HOME/sgg (writable,
  survives package upgrades/reinstalls).

User entries take precedence over bundled ones when names collide.
"""

import os
from pathlib import Path

PACKAGE_DIR = Path(__file__).resolve().parent
BUNDLED_DATA_DIR = PACKAGE_DIR / "data"
BUNDLED_SEQ_DIR = BUNDLED_DATA_DIR / "sequences"
BUNDLED_HOSTS_DIR = BUNDLED_DATA_DIR / "hosts"
BUNDLED_LOCATION_DIR = BUNDLED_DATA_DIR / "locations"


def _xdg_data_home() -> Path:
    return Path(os.environ.get("XDG_DATA_HOME") or (Path.home() / ".local" / "share"))


USER_DATA_DIR = _xdg_data_home() / "sgg"
USER_SEQ_DIR = USER_DATA_DIR / "sequences"


def list_entries() -> dict[str, Path]:
    """Map entry name -> directory, merging bundled and user entries (user wins on clash)."""
    entries: dict[str, Path] = {}
    for base in (BUNDLED_SEQ_DIR, USER_SEQ_DIR):
        if base.is_dir():
            for p in sorted(base.iterdir()):
                if p.is_dir():
                    entries[p.name] = p
    return entries


def resolve_entry(name: str) -> Path:
    entries = list_entries()
    if name not in entries:
        available = ", ".join(sorted(entries)) or "(none found)"
        raise SystemExit(f"Unknown entry '{name}'. Available entries: {available}")
    return entries[name]
