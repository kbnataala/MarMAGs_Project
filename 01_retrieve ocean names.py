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
