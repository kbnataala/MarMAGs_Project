################################################################################
# Script: 01_geospatial_enrichment.py
# Purpose:
#   Enrich MAG coordinate metadata with Ocean/Sea/EEZ names using shapefiles.
#
# Input:
#   data/metadata/mag_coordinates_raw.csv
#     Required columns: Genome_ID, sample_latitude, sample_longitude
#
# Shapefiles (stored locally; typically NOT committed to GitHub):
#   data/shapefiles/goas/<GOaS shapefile>.shp
#   data/shapefiles/eez_iho/<EEZ/IHO shapefile>.shp
#
# Output:
#   data/processed/mag_coordinates_enriched.csv
#
# Notes:
#   - Rows with missing coordinates are retained but will have NA enrichments.
#   - Uses EPSG:4326 (WGS84) coordinates.
################################################################################

import os
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point


# ------------------------------------------------------------------------------
# CONFIG (relative paths only)
# ------------------------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DATA_DIR = os.path.join(BASE_DIR, "data")
SHAPE_DIR = os.path.join(DATA_DIR, "shapefiles")
META_DIR = os.path.join(DATA_DIR, "metadata")
PROC_DIR = os.path.join(DATA_DIR, "processed")

INPUT_FILE = os.path.join(META_DIR, "mag_coordinates_raw.csv")
OUTPUT_FILE = os.path.join(PROC_DIR, "mag_coordinates_enriched.csv")

# Update these filenames to match the actual .shp names inside your folders
OCEAN_SHP = os.path.join(SHAPE_DIR, "goas", "goas_v01.shp")
SEA_SHP = os.path.join(SHAPE_DIR, "eez_iho", "Intersect_EEZ_IHO_v5.shp")

os.makedirs(PROC_DIR, exist_ok=True)


def _pick_first_shp(folder: str) -> str:
    """If the expected shapefile name is unknown, pick the first .shp in folder."""
    if not os.path.isdir(folder):
        raise FileNotFoundError(f"Shapefile folder not found: {folder}")
    candidates = [f for f in os.listdir(folder) if f.lower().endswith(".shp")]
    if not candidates:
        raise FileNotFoundError(f"No .shp files found in: {folder}")
    candidates.sort()
    return os.path.join(folder, candidates[0])


# If the hardcoded filenames above don't exist, auto-detect the .shp
if not os.path.exists(OCEAN_SHP):
    OCEAN_SHP = _pick_first_shp(os.path.join(SHAPE_DIR, "goas"))

if not os.path.exists(SEA_SHP):
    SEA_SHP = _pick_first_shp(os.path.join(SHAPE_DIR, "eez_iho"))


# ------------------------------------------------------------------------------
# LOAD INPUT
# ------------------------------------------------------------------------------
if not os.path.exists(INPUT_FILE):
    raise FileNotFoundError(
        f"Missing input file: {INPUT_FILE}\n"
        "Place your raw file at: data/metadata/mag_coordinates_raw.csv"
    )

df = pd.read_csv(INPUT_FILE)

required = {"Genome_ID", "sample_latitude", "sample_longitude"}
missing = required.difference(df.columns)
if missing:
    raise ValueError(f"Input file is missing required columns: {missing}")

# Keep all rows, but only spatial-join those with valid coordinates
df_valid = df.dropna(subset=["sample_latitude", "sample_longitude"]).copy()

print(f"Total rows: {len(df)}")
print(f"Rows with valid coordinates: {len(df_valid)}")
print(f"Rows with missing coordinates: {len(df) - len(df_valid)}")

# ------------------------------------------------------------------------------
# BUILD POINTS GDF
# ------------------------------------------------------------------------------
gdf_points = gpd.GeoDataFrame(
    df_valid[["Genome_ID", "sample_latitude", "sample_longitude"]].copy(),
    geometry=[
        Point(xy) for xy in zip(df_valid["sample_longitude"], df_valid["sample_latitude"])
    ],
    crs="EPSG:4326",
)

# ------------------------------------------------------------------------------
# LOAD SHAPEFILES
# ------------------------------------------------------------------------------
print(f"Using GOaS shapefile: {OCEAN_SHP}")
print(f"Using EEZ/IHO shapefile: {SEA_SHP}")

oceans = gpd.read_file(OCEAN_SHP)
seas = gpd.read_file(SEA_SHP)

# ------------------------------------------------------------------------------
# OCEAN JOIN
# ------------------------------------------------------------------------------
oceans_join = gpd.sjoin(gdf_points, oceans, how="left", predicate="within")

# Try common name columns (shapefiles vary)
ocean_name_col_candidates = ["name", "NAME", "Ocean", "OCEAN"]
ocean_name_col = next((c for c in ocean_name_col_candidates if c in oceans_join.columns), None)

df_valid["Ocean_Name"] = oceans_join[ocean_name_col].values if ocean_name_col else pd.NA
if ocean_name_col is None:
    print("Warning: No ocean name column found in GOaS shapefile. Ocean_Name will be NA.")

# ------------------------------------------------------------------------------
# SEA / EEZ JOIN
# ------------------------------------------------------------------------------
sea_cols_wanted = ["MarRegion", "IHO_Sea", "EEZ", "SOVEREIGN1", "TERRITORY1", "AREA_KM2"]
sea_cols = [c for c in sea_cols_wanted if c in seas.columns]

seas_join = gpd.sjoin(
    gdf_points,
    seas[sea_cols + ["geometry"]],
    how="left",
    predicate="within",
)

for c in sea_cols:
    df_valid[c] = seas_join[c].values

# ------------------------------------------------------------------------------
# MERGE BACK (preserve original rows, fill enrichments where valid coords exist)
# ------------------------------------------------------------------------------
df_out = df.copy()

# Add new columns with NA defaults
new_cols = ["Ocean_Name"] + sea_cols
for c in new_cols:
    if c not in df_out.columns:
        df_out[c] = pd.NA

# Fill enrichments for rows that had valid coords
df_out.loc[df_valid.index, new_cols] = df_valid[new_cols]

# ------------------------------------------------------------------------------
# SAVE
# ------------------------------------------------------------------------------
df_out.to_csv(OUTPUT_FILE, index=False)
print(f"Saved enriched metadata to: {OUTPUT_FILE}")
