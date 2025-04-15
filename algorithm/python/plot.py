import pandas as pd
import networkx as nx
import folium
from itertools import pairwise
import colorsys
import webbrowser
import os
import glob
from typing import Dict, List, Tuple

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

def plot_multiple_routes(
    G: nx.MultiDiGraph,
    warehouse_routes: Dict[str, List[str]],
    warehouse_names: Dict[str, str],
    output_file: str = 'multi_warehouse_routes.html',
    map_style: str = 'OpenStreetMap'
) -> folium.Map:
    """Plot multiple warehouse routes with proper type checking"""
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
            center_points.append((float(node['y']), float(node['x'])))  # Ensure float
    
    map_center = (28.6139, 77.2090)  # Default Delhi center
    if center_points:
        avg_lat = sum(p[0] for p in center_points) / len(center_points)
        avg_lon = sum(p[1] for p in center_points) / len(center_points)
        map_center = (float(avg_lat), float(avg_lon))  # Ensure float

    # Create base map
    m = folium.Map(location=map_center, zoom_start=12, tiles=map_style)

    # Process each warehouse route
    for w_idx, (warehouse_id, node_sequence) in enumerate(warehouse_routes.items()):
        warehouse_name = warehouse_names.get(warehouse_id, f"Warehouse {w_idx+1}")
        route_group = folium.FeatureGroup(name=f"Route: {warehouse_name}").add_to(m)
        base_color = warehouse_base_colors[w_idx]
        
        route_segments = []
        
        # Calculate route segments
        for i, (u, v) in enumerate(pairwise(node_sequence)):
            segment_coords = []
            try:
                if str(u) not in G or str(v) not in G:
                    raise KeyError(f"Nodes not found: {u} or {v}")

                path = nx.shortest_path(G, str(u), str(v), weight='length')
                
                # Get coordinates for each segment
                for u_edge, v_edge in pairwise(path):
                    edge_data = G.get_edge_data(u_edge, v_edge)
                    if edge_data:
                        for key, data in edge_data.items():
                            if 'geometry' in data:
                                segment_coords.extend(parse_linestring(data['geometry']))
                                break
                        else:
                            # Fallback to straight line
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
                    segment_coords.extend([
                        (float(u_node['y']), float(u_node['x'])),
                        (float(v_node['y']), float(v_node['x']))
                    ])
                except KeyError:
                    print(f"Could not draw segment between nodes {u} and {v}")

            if segment_coords:
                route_segments.append(segment_coords)

        # Add segments to map
        if HAS_PLUGINS:
            for segment in route_segments:
                folium.plugins.AntPath(
                    locations=segment,
                    color=base_color,
                    weight=5,
                    opacity=0.8,
                    tooltip=warehouse_name,
                    delay=800
                ).add_to(route_group)
        else:
            all_coords = [coord for segment in route_segments for coord in segment]
            folium.PolyLine(
                locations=all_coords,
                color=base_color,
                weight=4,
                opacity=0.7,
                tooltip=warehouse_name
            ).add_to(route_group)

        # Add markers
        for i, node_id in enumerate(node_sequence):
            try:
                node = G.nodes[str(node_id)]
                if i == 0:  # Warehouse
                    folium.Marker(
                        location=[float(node['y']), float(node['x'])],
                        popup=warehouse_name,
                        icon=folium.Icon(color='red', icon='home')
                    ).add_to(m)
                else:  # Delivery point
                    folium.CircleMarker(
                        location=[float(node['y']), float(node['x'])],
                        radius=6,
                        popup=f"Delivery {i}",
                        color=base_color,
                        fill=True
                    ).add_to(m)
            except KeyError:
                print(f"Warning: Node {node_id} not found in graph")

    folium.LayerControl().add_to(m)
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
    MAP_STYLE = 'OpenStreetMap'
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
            map_style=MAP_STYLE
        )

        print("\nVisualization completed successfully!")
        
    except Exception as e:
        print(f"\nError: {str(e)}")
        print("Please check:")
        print("- CSV files exist and are properly formatted")
        print("- Numeric columns (x, y, length) contain valid numbers")
        print("- Node IDs are consistent between files")