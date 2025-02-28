# layoutgenerator
This repository contains a Python application that optimizes the placement of solar panels within defined geographical areas. The tool analyzes KML/KMZ files to determine viable installation zones and calculates the maximum number of solar panels that can be arranged within these boundaries.

Main Functions

get_user_inputs(): Collects panel specifications like type, dimensions, power, and spacing from the user.

calculate_table_dimensions(): Determines table dimensions based on panel type (Fixed Tilt or Single Axis Tracker).

get_kml_polygons(): Extracts polygon data from KML/KMZ files that match specific search terms.

transform_to_utm(): Projects geographic coordinates to UTM for accurate spatial calculations.

process_kml(): Identifies enabled areas, subtracts disabled areas, and calculates the effective installation area.

create_table_geometry(): Creates rectangular geometries representing solar panel tables.

calculate_layout(): Algorithmically places table geometries within the effective area (already in UTM coordinates) and tallies the number of tables that fit.

export_to_kmz(): Exports the calculated panel layout back to KMZ format for visualization.

The process flow should be understood as:

Load KML/KMZ file containing site boundaries
Define usable and unusable areas through polygon selection
Project coordinates to UTM for accurate measurements
Calculate table dimensions based on panel specifications
Mathematically position tables within the effective area according to spacing requirements
Calculate total power output based on number of tables and modules
Convert these computational geometries back to geographic coordinates
Export the geometric arrangement to KMZ for visualization

The calculations are all performed on geometric objects (Shapely polygons) in memory, not directly on the KMZ file. The KMZ is only used as input at the beginning and output at the end.
