from xml.etree import ElementTree
import zipfile
import os
import numpy as np
from shapely.geometry import shape, Polygon, MultiPolygon, box
from shapely.ops import unary_union, transform
import pyproj
from functools import partial
from shapely.affinity import rotate, translate

class SolarLayoutCalculator:
    def __init__(self):
        # Areas
        self.enabled_polygons = None
        self.disabled_polygons = None
        self.effective_area = None
        
        # Panel parameters
        self.panel_type = None
        self.module_width = None
        self.module_length = None
        self.modules_per_table = None
        self.module_power = None
        self.gap = None
        self.pitch = None
        
        # Table dimensions
        self.table_width = None
        self.table_length = None
        
        # Results
        self.total_tables = 0
        self.total_modules = 0
        self.total_power = 0
        self.table_geometries = []
        
        # Projections
        self.proj_wgs84 = pyproj.CRS('EPSG:4326')
        self.proj_utm = None
        self.transformer_to_wgs84 = None
    
    def get_user_inputs(self):
        """Gets all necessary parameters from the user"""
        while True:
            self.panel_type = input("Panel type (1: Fixed Tilt, 2: Single Axis Tracker): ").strip()
            if self.panel_type in ['1', '2']:
                break
            print("Please enter 1 or 2")
        
        self.module_width = float(input("Module width (meters): "))
        self.module_length = float(input("Module length (meters): "))
        self.modules_per_table = int(input("Modules per table: "))
        self.module_power = float(input("Module power (W): "))
        self.gap = float(input("Gap - Space between tables in the same row (meters): "))
        self.pitch = float(input("Pitch - Space between rows (meters): "))
        
        self.calculate_table_dimensions()
    
    def calculate_table_dimensions(self):
        """Calculates table dimensions according to the panel type"""
        if self.panel_type == '1':  # Fixed Tilt
            self.module_power = self.module_power * 2 # 2 modules in portrait
            self.table_width = self.module_length * 2 # 2 modules in portrait
            self.table_length = self.module_width * self.modules_per_table + self.gap * (self.modules_per_table - 1)
        else:  # Single Axis Tracker
            self.table_width = self.module_length
            self.table_length = self.module_width * self.modules_per_table + self.gap * (self.modules_per_table - 1)
    
    def get_kml_polygons(self, file_path, search_term):
        """Extracts polygons from KMZ/KML that match the search term"""
        polygons = []
        found_names = set()
        
        # Get KML namespace
        ns = {'kml': 'http://www.opengis.net/kml/2.2'}
        
        print(f"Looking for elements containing '{search_term}':")
        
        # Read KMZ or KML file
        if file_path.lower().endswith('.kmz'):
            with zipfile.ZipFile(file_path, 'r') as kmz:
                print(f"  Found: {os.path.basename(file_path)}")
                kml_file = next(f for f in kmz.namelist() if f.lower().endswith('.kml'))
                kml_content = kmz.read(kml_file)
        else:
            with open(file_path, 'rb') as f:
                kml_content = f.read()
        
        # Parse XML content
        root = ElementTree.fromstring(kml_content)
        
        # Search for all Placemarks
        for placemark in root.findall('.//kml:Placemark', ns):
            name = placemark.find('kml:name', ns)
            if name is not None:
                print(f"    Found: {name.text}")
                if search_term.lower() in name.text.lower():
                    found_names.add(name.text)
                    # Search for coordinates in the Polygon
                    coords_elem = placemark.find('.//kml:coordinates', ns)
                    if coords_elem is not None:
                        # Convert coordinates to polygon format
                        coords_text = coords_elem.text.strip()
                        coords_list = []
                        for coord in coords_text.split():
                            lon, lat, _ = map(float, coord.split(','))
                            coords_list.append((lon, lat))
                        
                        # Create Shapely polygon
                        if coords_list:
                            poly = Polygon(coords_list)
                            polygons.append(poly)
        
        print(f"Summary for '{search_term}':")
        print(f"Total polygons found: {len(polygons)}")
        if found_names:
            print("Names found:")
            for name in found_names:
                print(f"  - {name}")
        
        return polygons
    
    def transform_to_utm(self, geometry):
        """Transforms geometry from WGS84 to UTM"""
        # Get centroid to determine UTM zone
        centroid = geometry.centroid
        utm_zone = int((centroid.x + 180) / 6) + 1
        hemisphere = 'north' if centroid.y >= 0 else 'south'
        
        print(f"Using UTM zone {utm_zone} {hemisphere}")
        
        # Create specific UTM projection for the location
        self.proj_utm = pyproj.CRS(f'+proj=utm +zone={utm_zone} +{hemisphere} +ellps=WGS84')
        
        # Create transformers
        project = pyproj.Transformer.from_crs(self.proj_wgs84, self.proj_utm, always_xy=True).transform
        self.transformer_to_wgs84 = pyproj.Transformer.from_crs(self.proj_utm, self.proj_wgs84, always_xy=True).transform
        
        # Transform geometry
        return transform(project, geometry)
    
    def process_kml(self, file_path):
        """Processes the KML/KMZ and gets the effective area"""
        print("Processing file...")
        
        # Get search terms
        enabled_terms = input("Enter names of enabled polygons (separated by comma): ").split(',')
        enabled_terms = [term.strip() for term in enabled_terms]
        
        disabled_terms = input("Enter names of disabled polygons (separated by comma): ").split(',')
        disabled_terms = [term.strip() for term in disabled_terms]
        
        # Process enabled polygons
        enabled_polygons = []
        for term in enabled_terms:
            polygons = self.get_kml_polygons(file_path, term)
            if polygons:
                enabled_polygons.extend(polygons)
                print(f"Found {len(polygons)} polygons for term '{term}'")
            else:
                print(f"Warning: No polygons found with term '{term}'")
        
        if not enabled_polygons:
            raise ValueError("No enabled polygons found")
        
        print("Processing enabled area...")
        # Merge enabled polygons and transform to UTM
        self.enabled_polygons = unary_union(enabled_polygons)
        self.enabled_polygons = self.transform_to_utm(self.enabled_polygons)
        print("Enabled polygons processed and transformed to UTM")
        
        # Process disabled polygons
        if disabled_terms and disabled_terms[0]:
            print("Processing disabled areas...")
            disabled_polygons = []
            for term in disabled_terms:
                polygons = self.get_kml_polygons(file_path, term)
                if polygons:
                    print(f"Processing {len(polygons)} polygons of type '{term}'")
                    for i, poly in enumerate(polygons, 1):
                        # Transform to UTM
                        poly_utm = self.transform_to_utm(poly)
                        # Subtract from enabled area
                        self.enabled_polygons = self.enabled_polygons.difference(poly_utm)
                        print(f"  Subtracted polygon {i} successfully")
                        disabled_polygons.append(poly_utm)
            
            if disabled_polygons:
                self.disabled_polygons = unary_union(disabled_polygons)
                print("\nProcessing summary:")
                print(f"Total disabled polygons: {len(disabled_polygons)}")
                print(f"Successful subtractions: {len(disabled_polygons)}")
        
        # Assign effective area
        self.effective_area = self.enabled_polygons
        
        # If tracker, rotate the effective area 90 degrees
        if self.panel_type == '2':  # Single Axis Tracker
            print("Rotating effective area for N-S orientation...")
            centroid = self.effective_area.centroid
            self.effective_area = rotate(self.effective_area, -90, origin=centroid)
        
        # Print effective area information
        bounds = self.effective_area.bounds
        print(f"\nEffective area:")
        print(f"Dimensions (meters):")
        print(f"Width: {bounds[2] - bounds[0]:.2f}m")
        print(f"Height: {bounds[3] - bounds[1]:.2f}m")
        print(f"Total area: {self.effective_area.area:.2f}m²")
        
        print("Processing completed")
    
    def create_table_geometry(self):
        """Creates the geometry of a table"""
        return box(0, 0, self.table_length, self.table_width)
    
    def calculate_layout(self):
        """Calculates the arrangement of tables and total power"""
        if not self.effective_area:
            raise ValueError("You must process the KML first")
        
        print("\nCalculating panel arrangement...")
        print(f"Table dimensions: {self.table_length}m x {self.table_width}m")
        
        # Get the bounding box of the effective area
        minx, miny, maxx, maxy = self.effective_area.bounds
        
        # Create base table (without rotation, use the same logic for both types)
        base_table = self.create_table_geometry()
        
        # Clear previous geometries list
        self.table_geometries = []
        
        # Initialize counters
        self.total_tables = 0
        row_count = 0
        current_y = miny
        
        # Use the same placement logic for both types
        while current_y + self.pitch <= maxy:
            current_x = minx
            tables_in_row = 0
            
            while current_x + self.table_length <= maxx:
                table = translate(base_table, current_x, current_y)
                
                if self.effective_area.contains(table):
                    self.total_tables += 1
                    tables_in_row += 1
                    self.table_geometries.append(table)
                
                current_x += self.table_length + self.gap
            
            if tables_in_row > 0:
                row_count += 1
            
            current_y += self.pitch
        
        print(f"Rows of tables: {row_count}")
        
        # Calculate totals
        self.total_modules = self.total_tables * self.modules_per_table
        self.total_power = self.total_modules * self.module_power
        
        print("Calculation completed")
        
        return {
            'total_tables': self.total_tables,
            'total_modules': self.total_modules,
            'total_power_w': self.total_power,
            'total_power_mw': self.total_power / 1_000_000
        }
    
    def export_to_kmz(self, output_filename):
        """Exports the table arrangement to a KMZ file"""
        if not self.table_geometries:
            raise ValueError("No table geometries to export")
        
        # Create a new KML document
        kml_doc = ElementTree.Element('kml', xmlns="http://www.opengis.net/kml/2.2")
        document = ElementTree.SubElement(kml_doc, 'Document')
        
        # Define style for tables
        style = ElementTree.SubElement(document, 'Style', id="tableStyle")
        line_style = ElementTree.SubElement(style, 'LineStyle')
        ElementTree.SubElement(line_style, 'color').text = 'ff0000ff'  # Red
        ElementTree.SubElement(line_style, 'width').text = '2'
        poly_style = ElementTree.SubElement(style, 'PolyStyle')
        ElementTree.SubElement(poly_style, 'color').text = '7f0000ff'  # Semi-transparent red
        
        # If tracker, we need to rotate the geometries back before exporting
        if self.panel_type == '2':
            centroid = self.effective_area.centroid
            rotation_angle = 90
        else:
            rotation_angle = 0
        
        # Create a Placemark for each table
        for i, table in enumerate(self.table_geometries, 1):
            # If tracker, rotate the table back
            if self.panel_type == '2':
                table = rotate(table, rotation_angle, origin=centroid)
            
            # Transform coordinates back to WGS84
            table_wgs84 = transform(self.transformer_to_wgs84, table)
            
            # Create Placemark
            placemark = ElementTree.SubElement(document, 'Placemark')
            ElementTree.SubElement(placemark, 'name').text = f'Table_{i}'
            ElementTree.SubElement(placemark, 'styleUrl').text = '#tableStyle'
            
            # Create polygon
            polygon = ElementTree.SubElement(placemark, 'Polygon')
            ElementTree.SubElement(polygon, 'extrude').text = '1'
            ElementTree.SubElement(polygon, 'altitudeMode').text = 'relativeToGround'
            
            # Create outer ring
            outer_boundary = ElementTree.SubElement(polygon, 'outerBoundaryIs')
            linear_ring = ElementTree.SubElement(outer_boundary, 'LinearRing')
            
            # Get polygon coordinates
            coords = []
            for x, y in table_wgs84.exterior.coords:
                coords.append(f"{x},{y},0")
            
            # Add coordinates to KML
            ElementTree.SubElement(linear_ring, 'coordinates').text = ' '.join(coords)
        
        # Create temporary KML file
        temp_kml = 'temp.kml'
        tree = ElementTree.ElementTree(kml_doc)
        tree.write(temp_kml, encoding='utf-8', xml_declaration=True)
        
        # Create KMZ file
        with zipfile.ZipFile(output_filename, 'w') as kmz:
            kmz.write(temp_kml, arcname=os.path.basename(temp_kml))
        
        # Remove temporary KML file
        os.remove(temp_kml)
        
        print(f"KMZ file saved as: {output_filename}")

def main():
    """Main program function"""
    try:
        # Initialize calculator
        calculator = SolarLayoutCalculator()
        
        # Get file path
        file_path = input("Enter the path to the KML/KMZ file: ")
        
        # Get user inputs
        calculator.get_user_inputs()
        
        # Process file
        calculator.process_kml(file_path)
        
        # Calculate layout and power
        results = calculator.calculate_layout()
        
        # Show results
        print("\nCalculation results:")
        print(f"Total number of tables: {results['total_tables']}")
        print(f"Total number of modules: {results['total_modules']}")
        print(f"Total power: {results['total_power_mw']:.2f} MW")
        
        # Export to KMZ
        output_filename = input("Enter the output KMZ filename (e.g., layout.kmz): ")
        print("Exporting to KMZ...")
        calculator.export_to_kmz(output_filename)
        
    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    main()