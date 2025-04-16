import pandas as pd
import networkx as nx
import math
import folium
from itertools import pairwise
import colorsys
import webbrowser
import os
import glob
from typing import Dict, List, Tuple, Optional

def parse_linestring(geom_str: str) -> List[Tuple[float, float]]:
    """Parse LINESTRING coordinates from WKT format"""
    if pd.isna(geom_str) or not isinstance(geom_str, str):
        return []

    try:
        coords_str = geom_str.split('(', 1)[1].rsplit(')', 1)[0]
        coord_pairs = [pair.strip() for pair in coords_str.split(',')]
        return [(float(lat), float(lon)) for lon, lat in [pair.split() for pair in coord_pairs]]
    except (IndexError, ValueError):
        return []

def load_osm_graph(nodes_file: str, edges_file: str) -> nx.MultiDiGraph:
    """Create a NetworkX graph from local OSM CSV files with proper type handling"""
    G = nx.MultiDiGraph()

    # Load nodes with proper type conversion
    print(f"Loading nodes from {nodes_file}...")
    nodes_df = pd.read_csv(nodes_file)
    
    # Convert osmid to string and coordinates to float
    nodes_df['osmid'] = nodes_df['osmid'].astype(str)
    nodes_df['x'] = pd.to_numeric(nodes_df['x'], errors='coerce')
    nodes_df['y'] = pd.to_numeric(nodes_df['y'], errors='coerce')
    
    # Drop any rows with invalid coordinates
    nodes_df = nodes_df.dropna(subset=['x', 'y'])
    
    for _, row in nodes_df.iterrows():
        G.add_node(str(row['osmid']), y=float(row['y']), x=float(row['x']))
    print(f"Loaded {len(nodes_df)} nodes")

    # Load edges with proper type conversion
    print(f"Loading edges from {edges_file}...")
    edges_df = pd.read_csv(edges_file, low_memory=False)
    
    # Convert edge endpoints to strings and length to float
    edges_df['u'] = edges_df['u'].astype(str)
    edges_df['v'] = edges_df['v'].astype(str)
    edges_df['length'] = pd.to_numeric(edges_df['length'], errors='coerce')
    
    # Drop any rows with invalid lengths
    edges_df = edges_df.dropna(subset=['length'])
    
    for _, row in edges_df.iterrows():
        G.add_edge(
            str(row['u']),
            str(row['v']),
            key=row['key'],
            osmid=str(row['osmid']),
            highway=row['highway'],
            name=row['name'],
            length=float(row['length']),  # Ensure length is float
            oneway=row['oneway'],
            geometry=row['geometry']
        )

        if not row['oneway']:
            G.add_edge(
                str(row['v']),
                str(row['u']),
                key=row['key'],
                osmid=str(row['osmid']),
                highway=row['highway'],
                name=row['name'],
                length=float(row['length']),  # Ensure length is float
                oneway=False,
                geometry=row['geometry']
            )
    print(f"Loaded {len(edges_df)} edges ({G.number_of_edges()} directed edges)")

    return G

def generate_warehouse_colors(num_warehouses: int) -> List[str]:
    """Generate visually distinct colors for different warehouses"""
    colors = [
        '#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00',
        '#FFFF33', '#A65628', '#F781BF', '#999999'
    ]
    
    if num_warehouses > len(colors):
        for i in range(len(colors), num_warehouses):
            h = i / num_warehouses
            s, v = 0.8, 0.9
            r, g, b = colorsys.hsv_to_rgb(h, s, v)
            colors.append(f'#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}')
    
    return colors[:num_warehouses]

def calculate_total_route_distance(G: nx.MultiDiGraph, route: List[str]) -> float:
    """Calculate the total distance of a route"""
    total_distance = 0.0
    
    for i in range(len(route) - 1):
        try:
            path = nx.shortest_path(G, str(route[i]), str(route[i+1]), weight='length')
            segment_distance = 0.0
            
            for j in range(len(path) - 1):
                edge_data = G.get_edge_data(path[j], path[j+1])
                if edge_data:
                    for key, data in edge_data.items():
                        if 'length' in data:
                            segment_distance += float(data['length'])
                            break
                            
            total_distance += segment_distance
        except (nx.NetworkXNoPath, KeyError) as e:
            print(f"Warning: Could not calculate distance between {route[i]} and {route[i+1]}: {e}")
    
    return total_distance

def extract_area_labels(G: nx.MultiDiGraph) -> List[Tuple[float, float, str]]:
    """Extract area labels from edge data (road names)"""
    area_labels = {}
    
    # Collect all road/area names from the graph
    for u, v, data in G.edges(data=True):
        if 'name' in data and data['name'] and isinstance(data['name'], str):
            name = data['name'].strip()
            if name and len(name) > 3:  # Filter out very short names
                # Get coordinates for the middle of the edge
                u_node, v_node = G.nodes[u], G.nodes[v]
                if 'x' in u_node and 'y' in u_node and 'x' in v_node and 'y' in v_node:
                    lat = (float(u_node['y']) + float(v_node['y'])) / 2
                    lon = (float(u_node['x']) + float(v_node['x'])) / 2
                    
                    # Use name as key to avoid duplicates in same area
                    key = f"{name}_{int(lat*100)}_{int(lon*100)}"
                    if key not in area_labels:
                        area_labels[key] = (lat, lon, name)
    
    # Convert to list and limit to a reasonable number to avoid clutter
    all_labels = list(area_labels.values())
    
    # Choose labels that are well-distributed
    # First sort by y-coordinate
    all_labels.sort(key=lambda x: x[0])
    
    # Take every Nth label to distribute them evenly
    n = max(1, len(all_labels) // 30)  # Limit to ~30 labels
    return all_labels[::n]

def plot_multiple_routes(
    G: nx.MultiDiGraph,
    warehouse_routes: Dict[str, List[str]],
    warehouse_names: Dict[str, str],
    output_file: str = 'multi_warehouse_routes.html',
    default_map_style: str = 'OpenStreetMap'
) -> folium.Map:
    """Plot multiple warehouse routes with OpenStreetMap as default base laye"""
    try:
        import folium.plugins
        HAS_PLUGINS = True
    except ImportError:
        HAS_PLUGINS = False

    # Verify all warehouse nodes exist
    missing_warehouses = [wid for wid in warehouse_routes if wid not in G]
    if missing_warehouses:
        print(f"Warning: Missing warehouse nodes: {missing_warehouses}")

    warehouse_base_colors = generate_warehouse_colors(len(warehouse_routes))
    
    # Calculate map center
    center_points = []
    for warehouse_id, route in warehouse_routes.items():
        if route and warehouse_id in G.nodes:
            node = G.nodes[warehouse_id]
            center_points.append((float(node['y']), float(node['x'])))
    
    map_center = (28.6139, 77.2090)  # Default Delhi center
    if center_points:
        avg_lat = sum(p[0] for p in center_points) / len(center_points)
        avg_lon = sum(p[1] for p in center_points) / len(center_points)
        map_center = (float(avg_lat), float(avg_lon))

    # Create base map with title
    m = folium.Map(location=map_center, zoom_start=12, tiles='OpenStreetMap')
    
    # Add multiple basemap options with proper attribution
    folium.TileLayer('OpenStreetMap', attr='© OpenStreetMap contributors').add_to(m)
    folium.TileLayer('Stamen Terrain', attr='Map tiles by Stamen Design').add_to(m)
    folium.TileLayer('CartoDB dark_matter', attr='© CartoDB').add_to(m)
    folium.TileLayer('CartoDB positron', attr='© CartoDB').add_to(m)
    
    # Add ESRI satellite layer
    folium.TileLayer(
        tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
        attr='Tiles © Esri',
        name='Satellite View'
    ).add_to(m)

    # Add title
    title_html = """
    <div style="position: fixed; 
                top: 10px; left: 50px; width: 90%; 
                background-color: rgba(255, 255, 255, 0.8);
                border-radius: 10px; box-shadow: 0 0 15px rgba(0, 0, 0, 0.2);
                text-align: center; z-index: 9999; padding: 10px;">
        <h3 style="margin: 0; color: #1A5276; font-family: Arial, sans-serif;">
            <b>Route Optimization for Delivery Services</b>
        </h3>
    </div>
    """
    m.get_root().html.add_child(folium.Element(title_html))
    
    # Store total distances for each warehouse
    warehouse_distances = {}
    
    # Process each warehouse route
    for w_idx, (warehouse_id, node_sequence) in enumerate(warehouse_routes.items()):
        warehouse_name = warehouse_names.get(warehouse_id, f"Warehouse {w_idx+1}")
        route_group = folium.FeatureGroup(name=f"Route: {warehouse_name}").add_to(m)
        base_color = warehouse_base_colors[w_idx]
        
        route_segments = []
        segment_distances = []
        segment_nodes = []    # Store start and end nodes for each segment
        
        # Calculate route segments and distances
        for i, (u, v) in enumerate(pairwise(node_sequence)):
            segment_coords = []
            segment_distance = 0.0
            
            try:
                if str(u) not in G or str(v) not in G:
                    raise KeyError(f"Nodes not found: {u} or {v}")

                path = nx.shortest_path(G, str(u), str(v), weight='length')
                
                for u_edge, v_edge in pairwise(path):
                    edge_data = G.get_edge_data(u_edge, v_edge)
                    if edge_data:
                        for key, data in edge_data.items():
                            if 'length' in data:
                                segment_distance += float(data['length'])
                            
                            if 'geometry' in data:
                                segment_coords.extend(parse_linestring(data['geometry']))
                                break
                        else:
                            u_node = G.nodes[u_edge]
                            v_node = G.nodes[v_edge]
                            segment_coords.extend([
                                (float(u_node['y']), float(u_node['x'])),
                                (float(v_node['y']), float(v_node['x']))
                            ])

            except (nx.NetworkXNoPath, KeyError) as e:
                print(f"Warning: {e}. Drawing straight line between {u} and {v}")
                try:
                    u_node = G.nodes[str(u)]
                    v_node = G.nodes[str(v)]
                    segment_distance = math.sqrt(
                        (float(u_node['x']) - float(v_node['x']))**2 + 
                        (float(u_node['y']) - float(v_node['y']))**2
                    ) * 111000
                    
                    segment_coords.extend([
                        (float(u_node['y']), float(u_node['x'])),
                        (float(v_node['y']), float(v_node['x']))
                    ])
                except KeyError:
                    print(f"Could not draw segment between nodes {u} and {v}")

            if segment_coords:
                route_segments.append(segment_coords)
                segment_distances.append(segment_distance)
                
                # Store node IDs for this segment
                if i == 0:  # First segment starts at warehouse
                    node_type_u = "Warehouse"
                else:
                    node_type_u = f"Delivery {i}"
                
                node_type_v = f"Delivery {i+1}"
                segment_nodes.append((node_type_u, node_type_v))

        # Calculate total distance
        total_distance = sum(segment_distances)
        warehouse_distances[warehouse_name] = total_distance
        
        # Add segments to map with simplified tooltips
        for idx, (segment, segment_distance, node_pair) in enumerate(
            zip(route_segments, segment_distances, segment_nodes)
        ):
            # Create simplified tooltip with HTML
# Replace the tooltip_html in the segment plotting section with this:
            tooltip_html = f"""
            <div style="font-family: Arial, sans-serif; 
                        padding: 8px; 
                        min-width: 200px;
                        max-width: 300px;
                        white-space: nowrap;">
                <div style="display: flex; align-items: center; margin-bottom: 5px;">
                    <div style="width: 12px; height: 12px; background-color: {base_color}; 
                                margin-right: 8px;"></div>
                    <h4 style="margin: 0; color: {base_color}; font-size: 14px;">
                        {warehouse_name}
                    </h4>
                </div>
                <div style="border-top: 1px solid #eee; padding-top: 5px;">
                    <div style="display: flex; justify-content: space-between;">
                        <span style="font-weight: bold;">From:</span>
                        <span>{node_pair[0]}</span>
                    </div>
                    <div style="display: flex; justify-content: space-between;">
                        <span style="font-weight: bold;">To:</span>
                        <span>{node_pair[1]}</span>
                    </div>
                    <div style="display: flex; justify-content: space-between;
                                margin-top: 5px;">
                        <span style="font-weight: bold;">Distance:</span>
                        <span>{segment_distance/1000:.2f} km ({segment_distance:.0f} m)</span>
                    </div>
                </div>
            </div>
            """
            
            # Add segment to map with the simplified tooltip
            if HAS_PLUGINS:
                folium.plugins.AntPath(
                    locations=segment,
                    color=base_color,
                    weight=5,
                    opacity=0.8,
                    tooltip=folium.Tooltip(tooltip_html),
                    delay=800
                ).add_to(route_group)
            else:
                folium.PolyLine(
                    locations=segment,
                    color=base_color,
                    weight=4,
                    opacity=0.7,
                    tooltip=folium.Tooltip(tooltip_html)
                ).add_to(route_group)

        # Add markers with hover tooltips
        for i, node_id in enumerate(node_sequence):
            try:
                node = G.nodes[str(node_id)]
                if i == 0:  # Warehouse
                    folium.Marker(
                        location=[float(node['y']), float(node['x'])],
                        popup=f"{warehouse_name}<br>Total Distance: {total_distance/1000:.2f} km",
                        tooltip=f"<strong>{warehouse_name}</strong>",
                        icon=folium.Icon(color='red', icon='home')
                    ).add_to(m)
                else:  # Delivery point
                    folium.CircleMarker(
                        location=[float(node['y']), float(node['x'])],
                        radius=6,
                        popup=f"{warehouse_name}: Delivery {i}",
                        tooltip=f"<strong>Delivery {i}</strong>",
                        color=base_color,
                        fill=True
                    ).add_to(m)
            except KeyError:
                print(f"Warning: Node {node_id} not found in graph")

    # Add legend with total distances (positioned bottom left)
    legend_html = f"""
    <div style="position: fixed; 
                bottom: 50px; left: 50px; 
                border: 2px solid grey; z-index: 9999; 
                background-color: white; padding: 10px;
                border-radius: 5px; box-shadow: 0 0 15px rgba(0, 0, 0, 0.2);
                max-width: 250px;">
        <h4 style="margin-top: 0; margin-bottom: 5px;">Route Distances</h4>
        <table style="width: 100%;">
            {"".join(f"""
            <tr>
                <td style="padding: 2px 10px 2px 0; text-align: left;">{name}</td>
                <td style="padding: 2px 0; text-align: right;">{dist/1000:.2f} km</td>
            </tr>
            """ for name, dist in warehouse_distances.items())}
        </table>
    </div>
    """
    m.get_root().html.add_child(folium.Element(legend_html))
    
    # Add instructions popup for the hover feature
    instructions_html = """
    <div style="position: fixed; 
                top: 120px; right: 10px;
                background-color: rgba(255, 255, 255, 0.9);
                border-radius: 5px; box-shadow: 0 0 10px rgba(0, 0, 0, 0.2);
                padding: 10px; z-index: 1000; max-width: 250px;
                font-family: Arial, sans-serif; font-size: 12px;">
        <p style="margin: 0; font-weight: bold;">👆 Hover over a route segment to see distance information</p>
    </div>
    """
    m.get_root().html.add_child(folium.Element(instructions_html))
    
    # Add layer control (positioned lower)
    folium.LayerControl(
        position='bottomright',
        collapsed=True,
        autoZIndex=True
    ).add_to(m)
    
    # Add fullscreen option
    if HAS_PLUGINS:
        folium.plugins.Fullscreen(
            position='topright',
            title='Expand view',
            title_cancel='Exit fullscreen',
            force_separate_button=True
        ).add_to(m)
    
    # Add measure tool
    if HAS_PLUGINS:
        folium.plugins.MeasureControl(
            position='topright',
            primary_length_unit='kilometers',
            active_color='red',
            completed_color='green'
        ).add_to(m)
    
    m.save(output_file)
    print(f"Map saved to {output_file}")
    webbrowser.open(f'file://{os.path.abspath(output_file)}')
    return m

def process_route_files(route_file_pattern: str) -> Tuple[Dict[str, List[str]], Dict[str, str]]:
    """Process route files with proper node ID handling"""
    warehouse_routes = {}
    warehouse_names = {}
    
    route_files = glob.glob(route_file_pattern)
    print(f"Found {len(route_files)} route files")
    
    for i, file_path in enumerate(route_files):
        try:
            warehouse_name = f"Warehouse {i+1}"
            filename = os.path.basename(file_path)
            
            if "warehouse" in filename.lower():
                try:
                    w_num = int(''.join(filter(str.isdigit, filename)))
                    warehouse_name = f"Warehouse {w_num}"
                except ValueError:
                    pass
            
            with open(file_path, 'r') as f:
                node_ids = [str(line.strip()) for line in f if line.strip()]
            
            if node_ids:
                warehouse_id = str(node_ids[0])
                warehouse_routes[warehouse_id] = node_ids
                warehouse_names[warehouse_id] = warehouse_name
                print(f"Loaded route for {warehouse_name} with {len(node_ids)} nodes")
                
        except Exception as e:
            print(f"Error processing route file {file_path}: {e}")
    
    return warehouse_routes, warehouse_names

if __name__ == "__main__":
    # Configuration
    NODES_FILE = '../../data/nodes.csv'
    EDGES_FILE = '../../data/edges.csv'
    OUTPUT_FILE = 'multi_warehouse_routes.html'
    DEFAULT_MAP_STYLE = 'OpenStreetMap'
    ROUTE_FILE_PATTERN = '../cpp/route_warehouse_*.txt'

    try:
        print("Starting route visualization...")
        
        # Load routes
        warehouse_routes, warehouse_names = process_route_files(ROUTE_FILE_PATTERN)
        if not warehouse_routes:
            print("Error: No valid warehouse routes found")
            exit(1)
        
        # Load road network with proper numeric conversion
        road_network = load_osm_graph(NODES_FILE, EDGES_FILE)
        
        # Generate and display map
        route_map = plot_multiple_routes(
            road_network,
            warehouse_routes,
            warehouse_names,
            OUTPUT_FILE,
            default_map_style=DEFAULT_MAP_STYLE
        )

        print("\nVisualization completed successfully!")
        print("Additional features:")
        print("- Multiple map styles available via the layer control")
        print("- Toggle area labels and major roads using the buttons")
        print("- Interactive distance information on hover")
        print("- Full screen option in the top right corner")
        
    except Exception as e:
        print(f"\nError: {str(e)}")
        print("Please check:")
        print("- CSV files exist and are properly formatted")
        print("- Numeric columns (x, y, length) contain valid numbers")
        print("- Node IDs are consistent between files")