"""Resolve FASTA header location fields to country/coordinates (used by `sgg stats`)."""

import os

import pandas as pd

from sgg.paths import BUNDLED_LOCATION_DIR

_FASTA_EXTS = {".fasta", ".fas", ".fa", ".txt"}


def _load_locations_database():
    df_china = pd.read_csv(os.path.join(BUNDLED_LOCATION_DIR, "china_provinces.csv"))
    df_us = pd.read_csv(os.path.join(BUNDLED_LOCATION_DIR, "US_capitals.csv"))
    df_world = pd.read_csv(os.path.join(BUNDLED_LOCATION_DIR, "countries.csv"))
    df_russia = pd.read_csv(os.path.join(BUNDLED_LOCATION_DIR, "russia_regions.csv"))
    return df_china, df_us, df_world, df_russia


_df_china, _df_us, _df_world, _df_russia = _load_locations_database()

_df_continent = pd.read_csv(os.path.join(BUNDLED_LOCATION_DIR, "continent_map.csv"))
CONTINENT_MAP = {
    k.replace("_", " "): v for k, v in zip(_df_continent["country"], _df_continent["continent"])
}

CENTROIDS = {
    "United States": (37.09, -95.71),
    "China": (35.86, 104.19),
    "Russia": (61.52, 105.31),
}


def parse_header_location(header: str):
    """Return (lat, lon, label) for a `virus|acc|genotype|host|country|region|year` header."""
    parts = header.split("|")
    if len(parts) < 7:
        return None, None, "Unknown"

    country = parts[4].replace("_", " ")
    region = parts[5].replace("_", " ") if parts[5] != "?" else None

    if country in ("UNKNOWN", "?", ""):
        return None, None, "Unknown"

    if country == "United States":
        if region:
            for _, row in _df_us.iterrows():
                if row["state"].lower() == region.lower():
                    return row["lat"], row["lon"], f"{region}, United States"
        return (*CENTROIDS["United States"], "United States")

    if country == "China":
        if region:
            for _, row in _df_china.iterrows():
                if row["province"].lower() == region.lower():
                    return row["lat"], row["lon"], f"{region}, China"
        return (*CENTROIDS["China"], "China")

    if country == "Russia":
        if region:
            for _, row in _df_russia.iterrows():
                if row["region_ru"].replace("_", " ").lower() == region.lower():
                    return row["lat"], row["lon"], f"{region}, Russia"
        return (*CENTROIDS["Russia"], "Russia")

    for _, row in _df_world.iterrows():
        if row["country"].replace("_", " ").lower() == country.lower():
            return row["lat"], row["lon"], country

    return None, None, "Unknown"


def build_location_dataframe(gene_path: str) -> pd.DataFrame:
    """One row per reference sequence: header, genotype, lat, lon, country label."""
    rows = []
    for file in os.listdir(gene_path):
        if os.path.splitext(file)[1] in _FASTA_EXTS:
            with open(os.path.join(gene_path, file), "r") as f:
                for line in f:
                    if line.startswith(">"):
                        header = line.strip().lstrip(">")
                        parts = header.split("|")
                        genotype = parts[2] if len(parts) > 2 else "Unknown"
                        lat, lon, label = parse_header_location(header)
                        rows.append(
                            {"header": header, "genotype": genotype, "lat": lat, "lon": lon, "label": label}
                        )
    return pd.DataFrame(rows, columns=["header", "genotype", "lat", "lon", "label"])
