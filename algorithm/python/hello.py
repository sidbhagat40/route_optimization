import pandas as pd
import networkx as nx
import folium
from itertools import pairwise
import ast
from typing import List, Tuple
import colorsys
import webbrowser
import os
import time

def parse_linestring(geom_str: str) -> List[Tuple[float, float]]:
    """
    Parse LINESTRING coordinates from WKT format
    Returns list of (lat, lon) tuples
    """
    if pd.isna(geom_str) or not isinstance(geom_str, str):
        return []

    try:
        coords_str = geom_str.split('(', 1)[1].rsplit(')', 1)[0]
        coord_pairs = [pair.strip() for pair in coords_str.split(',')]
        return [(float(lat), float(lon)) for lon, lat in [pair.split() for pair in coord_pairs]]
    except (IndexError, ValueError):
        return []

def load_osm_graph(nodes_file: str, edges_file: str) -> nx.MultiDiGraph:
    """
    Create a NetworkX graph from local OSM CSV files
    Args:
        nodes_file: Path to nodes.csv
        edges_file: Path to edges.csv
    Returns:
        NetworkX MultiDiGraph representing the road network
    """
    G = nx.MultiDiGraph()

    # Load nodes
    print(f"Loading nodes from {nodes_file}...")
    nodes_df = pd.read_csv(nodes_file)
    for _, row in nodes_df.iterrows():
        G.add_node(row['osmid'], y=row['y'], x=row['x'])
    print(f"Loaded {len(nodes_df)} nodes")

    # Load edges
    print(f"Loading edges from {edges_file}...")
    edges_df = pd.read_csv(edges_file, low_memory=False)  # Added low_memory=False to prevent DtypeWarning
    for _, row in edges_df.iterrows():
        # Add forward edge
        G.add_edge(
            row['u'],
            row['v'],
            key=row['key'],
            osmid=row['osmid'],
            highway=row['highway'],
            name=row['name'],
            length=row['length'],
            oneway=row['oneway'],
            geometry=row['geometry']
        )

        # Add reverse edge if not oneway
        if not row['oneway']:
            G.add_edge(
                row['v'],
                row['u'],
                key=row['key'],
                osmid=row['osmid'],
                highway=row['highway'],
                name=row['name'],
                length=row['length'],
                oneway=False,
                geometry=row['geometry']
            )
    print(f"Loaded {len(edges_df)} edges ({G.number_of_edges()} directed edges)")

    return G

def generate_gradient_colors(segments_count: int, start_color: tuple = (0, 0.7, 1), end_color: tuple = (0, 1, 0.7)) -> List[str]:
    """
    Generate a smooth gradient of colors for route segments
    Args:
        segments_count: Number of segments to generate colors for
        start_color: Starting color in HSV format (h, s, v)
        end_color: Ending color in HSV format (h, s, v)
    Returns:
        List of hex color codes
    """
    colors = []

    # Generate gradient in HSV space for smoother transitions
    for i in range(segments_count):
        # Calculate position in gradient (0 to 1)
        t = i / max(1, segments_count - 1)

        # Interpolate between start and end color
        h = start_color[0] + t * (end_color[0] - start_color[0])
        s = start_color[1] + t * (end_color[1] - start_color[1])
        v = start_color[2] + t * (end_color[2] - start_color[2])

        # Convert HSV to RGB
        r, g, b = colorsys.hsv_to_rgb(h, s, v)

        # Convert RGB to hex
        hex_color = f'#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}'
        colors.append(hex_color)

    return colors

def plot_route_with_geometry(
    G: nx.MultiDiGraph,
    node_sequence: List[int],
    output_file: str = 'osm_route.html',
    map_style: str = 'OpenStreetMap'
) -> folium.Map:
    """
    Plot route through specified nodes using actual road geometries
    with animated flow to show direction of travel
    Args:
        G: Road network graph
        node_sequence: List of OSM node IDs to visit in order
        output_file: Path to save HTML map
        map_style: Base map style to use ('OpenStreetMap', 'Stamen Terrain', etc.)
    Returns:
        folium.Map object with the route plotted
    """
    # Try to import folium plugins
    try:
        import folium.plugins
        HAS_PLUGINS = True
    except ImportError:
        HAS_PLUGINS = False
        print("folium.plugins not available. Some map features will be disabled.")
        print("To enable all features, install using: pip install folium")

    # Store route segments separately for animation
    route_segments = []
    segment_names = []
    missing_nodes = set()
    missing_edges = set()

    # Calculate route segments between each pair of nodes
    for i, (u, v) in enumerate(pairwise(node_sequence)):
        segment_coords = []
        try:
            # Verify nodes exist
            if u not in G or v not in G:
                raise KeyError(f"Nodes not found: {u} or {v}")

            # Find shortest path edges
            path = nx.shortest_path(G, u, v, weight='length')
            edges = list(zip(path[:-1], path[1:]))

            # Try to get edge name for the segment
            try:
                edge_name = G[path[0]][path[1]][0].get('name', f'Segment {i+1}')
                if pd.isna(edge_name) or edge_name == '':
                    edge_name = f'Segment {i+1}'
            except:
                edge_name = f'Segment {i+1}'

            segment_names.append(edge_name)

            # Get coordinates from edge geometries
            for u_edge, v_edge in edges:
                edge_data = G.get_edge_data(u_edge, v_edge)
                if not edge_data:
                    missing_edges.add((u_edge, v_edge))
                    continue

                for key, data in edge_data.items():
                    if 'geometry' in data:
                        segment_coords.extend(parse_linestring(data['geometry']))
                        break
                else:
                    # Fallback to straight line if no geometry
                    u_data = G.nodes[u_edge]
                    v_data = G.nodes[v_edge]
                    segment_coords.extend([(u_data['y'], u_data['x']),
                                        (v_data['y'], v_data['x'])])

        except (nx.NetworkXNoPath, KeyError) as e:
            print(f"Warning: {e}. Drawing straight line between {u} and {v}")
            try:
                u_data = G.nodes[u]
                v_data = G.nodes[v]
                segment_coords.extend([(u_data['y'], u_data['x']),
                                    (v_data['y'], v_data['x'])])
            except:
                print(f"Could not draw segment between nodes {u} and {v}")

        # Add the segment if it has coordinates
        if segment_coords:
            route_segments.append(segment_coords)

    # Create map centered on first segment's starting point
    if not route_segments or not route_segments[0]:
        raise ValueError("No valid route coordinates found")

    # Create the base map with the primary style
    m = folium.Map(
        location=route_segments[0][0],
        zoom_start=14,
        tiles=map_style,
        control_scale=True
    )

    # Add alternative tile layers with a layer control
    tile_options = {
        'OpenStreetMap': 'OpenStreetMap',
        'Cartodb Positron': 'CartoDB Positron',
        'Cartodb Dark Matter': 'CartoDB Dark Matter',
        'Stamen Terrain': 'Stamen Terrain',
        'Stamen Toner': 'Stamen Toner',
        'Stamen Watercolor': 'Stamen Watercolor',
        'ESRI World Street': 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Street_Map/MapServer/tile/{z}/{y}/{x}',
        'ESRI World Imagery': 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}'
    }

    # Add all tile layers except the base one
    for name, url in tile_options.items():
        if name != map_style:
            folium.TileLayer(url, name=name, attr="Map tiles").add_to(m)

    # Create a feature group for routes
    route_group = folium.FeatureGroup(name="Route Segments").add_to(m)

    # Generate beautiful gradient colors for segments
    colors = generate_gradient_colors(len(route_segments))

    # Get all coordinates for full route (for non-animated fallback)
    all_coords = [coord for segment in route_segments for coord in segment]

    # Add the animated segments (if plugins available)
    if HAS_PLUGINS:
        # Add animated segments with AntPath
        for i, (segment, name) in enumerate(zip(route_segments, segment_names)):
            # Use gradient color
            color = colors[i]

            # Add the animated path with improved styling
            folium.plugins.AntPath(
                locations=segment,
                color=color,
                weight=6,
                opacity=0.9,
                tooltip=f"Segment {i+1}: {name}",
                delay=800,
                dash_array=[10, 15],
                popup=f"Segment {i+1}: {name}",
                pulse_color='#FFFFFF'
            ).add_to(route_group)

    else:
        # Fallback: Add non-animated polyline if plugins not available
        folium.PolyLine(
            locations=all_coords,
            color='#0078FF',
            weight=5,
            opacity=0.8,
            tooltip="Route"
        ).add_to(route_group)

    # Add markers for each node with custom icons and colors
    for i, node_id in enumerate(node_sequence):
        try:
            node = G.nodes[node_id]

            if i == 0:  # Start point
                folium.Marker(
                    location=[node['y'], node['x']],
                    popup=f"Start: Node {node_id}",
                    tooltip="Start",
                    icon=folium.DivIcon(
                        html=f"""
                        <div style="background-color:#00c853; width:20px; height:20px;
                             border-radius:50%; display:flex; align-items:center; justify-content:center;
                             border:3px solid white; box-shadow:0 0 10px rgba(0,0,0,0.3);">
                        </div>
                        """
                    )
                ).add_to(m)
            elif i == len(node_sequence)-1:  # End point
                folium.Marker(
                    location=[node['y'], node['x']],
                    popup=f"End: Node {node_id}",
                    tooltip="End",
                    icon=folium.DivIcon(
                        html=f"""
                        <div style="background-color:#d50000; width:20px; height:20px;
                             border-radius:50%; display:flex; align-items:center; justify-content:center;
                             border:3px solid white; box-shadow:0 0 10px rgba(0,0,0,0.3);">
                        </div>
                        """
                    )
                ).add_to(m)
            else:  # Waypoints
                waypoint_color = colors[int((i / (len(node_sequence)-2)) * (len(colors)-1))]
                folium.Marker(
                    location=[node['y'], node['x']],
                    popup=f"Waypoint {i}: Node {node_id}",
                    tooltip=f"Waypoint {i}",
                    icon=folium.DivIcon(
                        html=f"""
                        <div style="background-color:{waypoint_color}; width:16px; height:16px;
                             border-radius:50%; display:flex; align-items:center; justify-content:center;
                             border:2px solid white; box-shadow:0 0 8px rgba(0,0,0,0.3);">
                        </div>
                        """
                    )
                ).add_to(m)
        except KeyError:
            missing_nodes.add(node_id)

    # Add warnings to map if needed
    if missing_nodes:
        folium.Marker(
            location=route_segments[0][0],
            icon=folium.DivIcon(
                html=f"""<div style="background-color: rgba(255, 0, 0, 0.7); color: white; padding: 5px; border-radius: 5px; font-weight: bold">
                    Warning: Missing nodes - {', '.join(map(str, missing_nodes))}
                </div>"""
            )
        ).add_to(m)

    if missing_edges:
        folium.Marker(
            location=route_segments[0][0],
            icon=folium.DivIcon(
                html=f"""<div style="background-color: rgba(255, 165, 0, 0.7); color: white; padding: 5px; border-radius: 5px; font-weight: bold; margin-top: 30px;">
                    Warning: Missing edges - {len(missing_edges)} pairs
                </div>"""
            )
        ).add_to(m)

    # Calculate total distance
    total_distance = sum(G[u][v][0]['length'] for u, v in pairwise(node_sequence)
                         if u in G and v in G and v in G[u])

    # Add extra features if plugins available
    if HAS_PLUGINS:
        folium.plugins.MiniMap(toggle_display=True).add_to(m)
        folium.plugins.Fullscreen().add_to(m)
        folium.plugins.MeasureControl(position='bottomleft', primary_length_unit='kilometers').add_to(m)

    # Add layer control
    folium.LayerControl(position='topright').add_to(m)

    # Save to HTML (this will overwrite existing file)
    m.save(output_file)
    print(f"Map saved to {output_file}")

    # Open in browser
    file_path = f'file://{os.path.abspath(output_file)}'
    webbrowser.open(file_path, new=2)  # new=2 opens in new tab if possible

    return m

if __name__ == "__main__":
    # Configuration
    NODES_FILE = 'nodes.csv'
    EDGES_FILE = 'edges.csv'
    OUTPUT_FILE = 'osm_route.html'  # Consistent filename for overwriting

    # Map style options
    MAP_STYLE = 'OpenStreetMap'

    # Your OSM node path
    NODE_PATH = [
        1827697182, 1826909611, 4233302722, 9884874655, 4229518792,
        4065531426, 9880357040, 12110038882, 2575648598, 920829184,
    ]

    try:
        # Load the road network
        road_network = load_osm_graph(NODES_FILE, EDGES_FILE)

        # Generate and display map
        route_map = plot_route_with_geometry(
            road_network,
            NODE_PATH,
            OUTPUT_FILE,
            map_style=MAP_STYLE
        )

        print("\nRoute plotting completed successfully!")
        print("Map automatically opened in your default browser")
        print(f"Note: Each run overwrites {OUTPUT_FILE} with the latest version")

    except Exception as e:
        print(f"\nError: {str(e)}")
        print("Please check:")
        print("- Your CSV files exist in the correct location")
        print("- The node IDs in NODE_PATH exist in your nodes.csv")
        print("- The edges connect the nodes in your path")
        print("- Required packages are installed: pandas, networkx, folium")