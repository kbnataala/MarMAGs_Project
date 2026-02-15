#######################################################
### Retrieving names of Oceans

import geopandas as gpd
from shapely.geometry import Point

# Load the shapefile
shapefile_path = "/gpfs1/schlecker/home/kabiruna/coord/GOaS_v1_20211214/goas_v01.shp"
oceans_seas = gpd.read_file(shapefile_path)

# Print columns in the shapefile for reference
print("Columns in the shapefile:")
print(oceans_seas.columns)

# Example coordinates to match
coordinates = [
    {"lat": 40.0, "lon": -70.0},  # Atlantic Ocean, near the U.S. east coast
    {"lat": 38.84917, "lon": -75.1076},  # Tasman Sea, near Australia
    {"lat": -62.231932, "lon": -58.655087},    # North Sea, near Norway
]

# Convert coordinates to GeoDataFrame
points = gpd.GeoDataFrame(
    coordinates,
    geometry=[Point(coord["lon"], coord["lat"]) for coord in coordinates],
    crs="EPSG:4326"  # WGS84
)

# Perform a spatial join to find which polygon contains each point
results = gpd.sjoin(points, oceans_seas, how="left", predicate="within")

# Print results of the spatial join
print("\nResults of the spatial join:")
print(results)

# Print the names of the seas or oceans for each coordinate
for index, row in results.iterrows():
    if not pd.isna(row["name"]):  # Check if a match is found
        print(f"Coordinates ({row['lat']}, {row['lon']}): {row['name']}")
    else:
        print(f"Coordinates ({row['lat']}, {row['lon']}): No match found")


#########################################################
## Retrieving names of Seas

import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

# Load the input CSV file
input_file = "data.csv"  # Replace with the path to your CSV file
data = pd.read_csv(input_file)

# Remove rows with NA coordinates and extract unique coordinates
data_cleaned = data.dropna(subset=["sample_latitude", "sample_longitude"])
unique_coords = data_cleaned[["sample_latitude", "sample_longitude"]].drop_duplicates()

# Convert unique coordinates to a GeoDataFrame
unique_coords_gdf = gpd.GeoDataFrame(
    unique_coords,
    geometry=[
        Point(xy) for xy in zip(unique_coords["sample_longitude"], unique_coords["sample_latitude"])
    ],
    crs="EPSG:4326"  # WGS84
)

# Load the shapefile
shapefile_path = "/gpfs1/schlecker/home/kabiruna/coord/Intersect_EEZ_IHO_v5_20241010/Intersect_EEZ_IHO_v5_20241010/Intersect_EEZ_IHO_v5_20241010.shp"  # Replace with the path to your shapefile
oceans_seas = gpd.read_file(shapefile_path)

# Debugging Step 1: Check the columns in the shapefile
print("Columns in the shapefile:", oceans_seas.columns)

# Select columns of interest from the shapefile for inclusion in the output
columns_of_interest = [
    "MarRegion", "SOVEREIGN1", "SOVEREIGN2", "SOVEREIGN3",
    "TERRITORY1", "TERRITORY2", "TERRITORY3",
    "EEZ", "IHO_Sea", "AREA_KM2"
]

# Ensure these columns exist in the shapefile
valid_columns = [col for col in columns_of_interest if col in oceans_seas.columns]
if not valid_columns:
    print("ERROR: None of the expected columns were found in the shapefile.")
    exit(1)

print("Using the following columns from the shapefile:", valid_columns)

# Perform a spatial join to find which ocean/sea contains each unique coordinate
results = gpd.sjoin(unique_coords_gdf, oceans_seas[valid_columns + ["geometry"]], how="left", predicate="within")

# Debugging Step 2: Display the first few rows of the spatial join results
print("Results of spatial join (first 5 rows):")
print(results.head())

# Add all selected columns to the unique coordinates DataFrame
for col in valid_columns:
    unique_coords[col] = results[col].values

# Debugging Step 3: Check for any missing matches
missing_matches = unique_coords[unique_coords[valid_columns[0]].isna()]
if not missing_matches.empty:
    print("Coordinates with missing matches:")
    print(missing_matches)

# Merge the additional data back into the original dataset
merged_data = pd.merge(
    data,
    unique_coords,
    on=["sample_latitude", "sample_longitude"],
    how="left"
)

# Save the result to a new CSV file
output_file = "output_with_full_names.csv"
merged_data.to_csv(output_file, index=False)

print(f"Output saved to {output_file}")
