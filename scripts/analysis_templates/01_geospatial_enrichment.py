#!/usr/bin/env python3
"""
01_geospatial_enrichment.py

Documentation template for geospatial enrichment of MAG/sample metadata.

Purpose
-------
This script illustrates how sample coordinates were enriched with ocean, sea,
marine region, and EEZ information using polygon shapefiles.

Repository role
---------------
This script is provided as a documentation template and may require adaptation
for local execution. It is not intended to function as a fully reproducible
pipeline component within this repository.

Expected inputs
---------------
1. A CSV table containing at least:
   - sample_latitude
   - sample_longitude

2. A GOaS shapefile for ocean names
3. An EEZ/IHO shapefile for marine region and sea names

Expected output
---------------
A CSV table with appended geospatial annotations.

Notes
-----
- Coordinates are assumed to be in WGS84 (EPSG:4326).
- Rows with missing coordinates are retained but will not receive annotations.
"""

from pathlib import Path
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

BASE_DIR = Path(__file__).resolve().parents[2]

INPUT_FILE = BASE_DIR / "data" / "metadata" / "mag_coordinates_raw.csv"
OUTPUT_FILE = BASE_DIR / "data" / "processed" / "mag_coordinates_enriched.csv"

GOAS_DIR = BASE_DIR / "data" / "shapefiles" / "goas"
EEZ_IHO_DIR = BASE_DIR / "data" / "shapefiles" / "eez_iho"


# ---------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------

def find_first_shapefile(folder: Path) -> Path:
    """Return the first .shp file found in a folder."""
    shp_files = sorted(folder.glob("*.shp"))
    if not shp_files:
        raise FileNotFoundError(f"No .shp file found in: {folder}")
    return shp_files[0]


def validate_input_columns(df: pd.DataFrame, required_cols: list[str]) -> None:
    """Raise an informative error if required columns are missing."""
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(
            f"Input file is missing required columns: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )


# ---------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------

def main() -> None:
    print("Starting geospatial enrichment template...")

    if not INPUT_FILE.exists():
        raise FileNotFoundError(
            f"Input file not found: {INPUT_FILE}\n"
            "Place the coordinate table in data/metadata/ before running."
        )

    GOAS_SHP = find_first_shapefile(GOAS_DIR)
    EEZ_IHO_SHP = find_first_shapefile(EEZ_IHO_DIR)

    print(f"Using input table: {INPUT_FILE}")
    print(f"Using GOaS shapefile: {GOAS_SHP}")
    print(f"Using EEZ/IHO shapefile: {EEZ_IHO_SHP}")

    # Read metadata
    df = pd.read_csv(INPUT_FILE)
    validate_input_columns(df, ["sample_latitude", "sample_longitude"])

    # Keep a copy of original table
    df_out = df.copy()

    # Select rows with usable coordinates
    df_valid = df.dropna(subset=["sample_latitude", "sample_longitude"]).copy()

    print(f"Total rows: {len(df)}")
    print(f"Rows with valid coordinates: {len(df_valid)}")
    print(f"Rows with missing coordinates: {len(df) - len(df_valid)}")

    if df_valid.empty:
        raise ValueError("No valid coordinates found in input table.")

    # Convert to GeoDataFrame
    gdf_points = gpd.GeoDataFrame(
        df_valid,
        geometry=[
            Point(xy)
            for xy in zip(df_valid["sample_longitude"], df_valid["sample_latitude"])
        ],
        crs="EPSG:4326",
    )

    # Load shapefiles
    oceans = gpd.read_file(GOAS_SHP)
    seas = gpd.read_file(EEZ_IHO_SHP)

    # ---------------------------
    # Ocean name enrichment
    # ---------------------------
    ocean_join = gpd.sjoin(gdf_points, oceans, how="left", predicate="within")

    possible_ocean_name_cols = ["name", "NAME", "Ocean", "OCEAN"]
    ocean_col = next((c for c in possible_ocean_name_cols if c in ocean_join.columns), None)

    if ocean_col is None:
        print("Warning: No ocean-name column found in GOaS shapefile.")
        df_valid["Ocean_Name"] = pd.NA
    else:
        df_valid["Ocean_Name"] = ocean_join[ocean_col].values

    # ---------------------------
    # Sea / marine region enrichment
    # ---------------------------
    desired_cols = [
        "MarRegion",
        "IHO_Sea",
        "EEZ",
        "SOVEREIGN1",
        "SOVEREIGN2",
        "SOVEREIGN3",
        "TERRITORY1",
        "TERRITORY2",
        "TERRITORY3",
        "AREA_KM2",
    ]
    available_cols = [col for col in desired_cols if col in seas.columns]

    sea_join = gpd.sjoin(
        gdf_points,
        seas[available_cols + ["geometry"]],
        how="left",
        predicate="within",
    )

    for col in available_cols:
        df_valid[col] = sea_join[col].values

    # ---------------------------
    # Merge annotations back
    # ---------------------------
    new_cols = ["Ocean_Name"] + available_cols
    for col in new_cols:
        if col not in df_out.columns:
            df_out[col] = pd.NA

    df_out.loc[df_valid.index, new_cols] = df_valid[new_cols]

    # Ensure output directory exists
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

    # Save
    df_out.to_csv(OUTPUT_FILE, index=False)
    print(f"Geospatially enriched table written to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
