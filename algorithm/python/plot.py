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
            key=row['key'] if 'key' in row else 0,
            osmid=str(row['osmid']) if 'osmid' in row else '0',
            highway=row['highway'] if 'highway' in row else '',
            name=row['name'] if 'name' in row else '',
            length=float(row['length']),  # Ensure length is float
            oneway=row['oneway'] if 'oneway' in row else False,
            geometry=row['geometry'] if 'geometry' in row else None
        )

        if 'oneway' in row and not row['oneway']:
            G.add_edge(
                str(row['v']),
                str(row['u']),
                key=row['key'] if 'key' in row else 0,
                osmid=str(row['osmid']) if 'osmid' in row else '0',
                highway=row['highway'] if 'highway' in row else '',
                name=row['name'] if 'name' in row else '',
                length=float(row['length']),  # Ensure length is float
                oneway=False,
                geometry=row['geometry'] if 'geometry' in row else None
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

def get_warehouse_name_by_location(warehouse_id: str, node_coords) -> str:
    """Map warehouse IDs to area names based on their location"""
    # Mapping for the three warehouses in the C++ code
    warehouse_mapping = {
        "12110005374": "North Delhi Warehouse",  # North Delhi
        "923637303": "East Delhi Warehouse",     # East Delhi
        "4314376797": "South Delhi Warehouse"    # South Delhi
    }
    
    if warehouse_id in warehouse_mapping:
        return warehouse_mapping[warehouse_id]
    
    # Fallback based on coordinates if we have new warehouses
    if node_coords:
        lat, lon = node_coords
        if lat > 28.65:
            return "North Delhi Warehouse"
        elif lon > 77.27:
            return "East Delhi Warehouse"
        elif lat < 28.55:
            return "South Delhi Warehouse"
        else:
            return "Central Delhi Warehouse"
            
    return f"Warehouse {warehouse_id}"

def process_alternative_path_files(alt_path_pattern: str) -> Dict[str, List[List[str]]]:
    """Process alternative path files from C++ format and ensure exactly 2 alternatives per warehouse"""
    alternative_paths = {}
    
    alt_files = glob.glob(alt_path_pattern)
    print(f"Found {len(alt_files)} alternative path files")
    
    # First process existing files
    for file_path in alt_files:
        try:
            filename = os.path.basename(file_path)  
            
            # Parse filename to extract warehouse ID and route number
            parts = filename.replace('.txt', '').split('_')
            if len(parts) >= 5 and parts[0] == "warehouse" and "alternative" in parts:
                warehouse_id = parts[1]
                route_num = parts[-1]
                
                # Read the path nodes
                with open(file_path, 'r') as f:
                    node_ids = [str(line.strip()) for line in f if line.strip()]
                
                if node_ids:
                    if warehouse_id not in alternative_paths:
                        alternative_paths[warehouse_id] = []
                    
                    # Only keep up to 2 alternative paths per warehouse
                    if len(alternative_paths[warehouse_id]) < 2:
                        alternative_paths[warehouse_id].append(node_ids)
                        print(f"Loaded alternative route {route_num} for warehouse {warehouse_id} with {len(node_ids)} nodes")
            
        except Exception as e:
            print(f"Error processing alternative path file {file_path}: {e}")
    
    # Ensure each of the three main warehouses has exactly 2 alternative routes
    main_warehouses = ["12110005374", "923637303", "4314376797"]  # North, East, South Delhi
    
    for warehouse_id in main_warehouses:
        if warehouse_id not in alternative_paths:
            alternative_paths[warehouse_id] = []
        
        # Generate exactly 2 alternatives if needed
        while len(alternative_paths[warehouse_id]) < 2:
            # For first alternative, create a simple circular pattern
            if len(alternative_paths[warehouse_id]) == 0:
                new_route = [warehouse_id] + [f"simulated_{warehouse_id}_{i}" for i in range(5,10)] + [warehouse_id]
            # For second alternative, create a different pattern
            else:
                new_route = [warehouse_id] + [f"simulated_{warehouse_id}_{i}" for i in range(10,15)] + [warehouse_id]
            
            alternative_paths[warehouse_id].append(new_route)
            print(f"Generated simulated alternative route {len(alternative_paths[warehouse_id])} for warehouse {warehouse_id}")
        
        # Ensure we don't have more than 2 alternatives
        if len(alternative_paths[warehouse_id]) > 2:
            alternative_paths[warehouse_id] = alternative_paths[warehouse_id][:2]
            print(f"Truncated to 2 alternative routes for warehouse {warehouse_id}")
    
    return alternative_paths 

def process_route_files(route_file_pattern: str) -> Tuple[Dict[str, List[str]], Dict[str, str]]:
    """Process optimal route files from C++ format"""
    warehouse_routes = {}
    warehouse_names = {}
    
    route_files = glob.glob(route_file_pattern)
    print(f"Found {len(route_files)} route files")
    
    for file_path in route_files:
        try:
            filename = os.path.basename(file_path)
            
            # Parse filename to extract warehouse ID
            # C++ format: warehouse_WAREHOUSEID_optimal_route.txt
            parts = filename.replace('.txt', '').split('_')
            if len(parts) >= 3 and parts[0] == "warehouse" and "optimal" in parts:
                warehouse_id = parts[1]
                
                with open(file_path, 'r') as f:
                    node_ids = [str(line.strip()) for line in f if line.strip()]
                
                if node_ids:
                    warehouse_routes[warehouse_id] = node_ids
                    # Warehouse names will be set properly in plot function after loading coords
                    warehouse_names[warehouse_id] = warehouse_id
                    print(f"Loaded optimal route for warehouse {warehouse_id} with {len(node_ids)} nodes")
                
        except Exception as e:
            print(f"Error processing route file {file_path}: {e}")
    
    return warehouse_routes, warehouse_names

def plot_multiple_routes(
    G: nx.MultiDiGraph,
    warehouse_routes: Dict[str, List[str]],
    warehouse_names: Dict[str, str],
    alternative_paths: Dict[str, List[List[str]]] = None,
    output_file: str = 'multi_warehouse_routes.html',
    default_map_style: str = 'OpenStreetMap'
) -> folium.Map:
    """Plot multiple warehouse routes and alternative paths with improved visualization"""
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
    
    # Update warehouse names based on location
    for warehouse_id, route in warehouse_routes.items():
        if route and warehouse_id in G.nodes:
            node = G.nodes[warehouse_id]
            center_points.append((float(node['y']), float(node['x'])))
            # Update warehouse name based on location
            warehouse_names[warehouse_id] = get_warehouse_name_by_location(
                warehouse_id, (float(node['y']), float(node['x']))
            )

    map_center = (28.6139, 77.2090)  # Default Delhi center
    if center_points:
        avg_lat = sum(p[0] for p in center_points) / len(center_points)
        avg_lon = sum(p[1] for p in center_points) / len(center_points)
        map_center = (float(avg_lat), float(avg_lon))

    # Create base map with title
    m = folium.Map(location=map_center, zoom_start=12, tiles='OpenStreetMap')
    
    # Add multiple basemap options
    folium.TileLayer('Stamen Terrain', attr='Map tiles by Stamen Design', show=False).add_to(m)
    folium.TileLayer('CartoDB dark_matter', attr='© CartoDB', show=False).add_to(m)
    folium.TileLayer('CartoDB positron', attr='© CartoDB', show=False).add_to(m)
    folium.TileLayer(
        tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
        attr='Tiles © Esri',
        name='Satellite View',
        show=False
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
    alt_distances = {}
    
    # Process each warehouse route
    for w_idx, (warehouse_id, node_sequence) in enumerate(warehouse_routes.items()):
        warehouse_name = warehouse_names.get(warehouse_id, f"Warehouse {w_idx+1}")
        route_group = folium.FeatureGroup(name=f"Optimal Route: {warehouse_name}").add_to(m)
        base_color = warehouse_base_colors[w_idx]
        
        route_segments = []
        segment_distances = []
        
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
                            
                            if 'geometry' in data and data['geometry']:
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

        # Calculate total distance
        total_distance = sum(segment_distances)
        warehouse_distances[warehouse_name] = total_distance
        
        # Add segments to map
        for idx, (segment, segment_distance) in enumerate(zip(route_segments, segment_distances)):
            tooltip_html = f"""
            <div style="font-family: Arial, sans-serif; padding: 8px; min-width: 200px;">
                <div style="display: flex; align-items: center; margin-bottom: 5px;">
                    <div style="width: 12px; height: 12px; background-color: {base_color}; 
                                margin-right: 8px;"></div>
                    <h4 style="margin: 0; color: {base_color}; font-size: 14px;">
                        {warehouse_name} - OPTIMAL ROUTE
                    </h4>
                </div>
                <div style="border-top: 1px solid #eee; padding-top: 5px;">
                    <div style="display: flex; justify-content: space-between;">
                        <span style="font-weight: bold;">Delivery Segment:</span>
                        <span>{idx+1}/{len(route_segments)}</span>
                    </div>
                    <div style="display: flex; justify-content: space-between;
                                margin-top: 5px;">
                        <span style="font-weight: bold;">Segment Distance:</span>
                        <span>{segment_distance/1000:.2f} km</span>
                    </div>
                </div>
            </div>
            """
            
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
                    weight=5,
                    opacity=0.8,
                    tooltip=folium.Tooltip(tooltip_html)
                ).add_to(route_group)

        # Add markers for warehouse and delivery points
        for i, node_id in enumerate(node_sequence):
            try:
                node = G.nodes[str(node_id)]
                location = [float(node['y']), float(node['x'])]
                
                if i == 0 or i == len(node_sequence) - 1:  # Warehouse (first and last point)
                    folium.Marker(
                        location=location,
                        popup=f"""
                        <div style="min-width: 200px; padding: 10px;">
                            <h4 style="margin-top: 0;">{warehouse_name}</h4>
                            <p><b>Node ID:</b> {node_id}</p>
                            <p><b>Total Route Distance:</b> {total_distance/1000:.2f} km</p>
                            <p><b>Deliveries:</b> {len(node_sequence)-2}</p>
                        </div>
                        """,
                        tooltip=warehouse_name,
                        icon=folium.Icon(color='red', icon='home', prefix='fa')
                    ).add_to(m)
                else:  # Delivery point
                    delivery_num = i
                    icon_color = base_color
                    
                    # Create circle marker with delivery number
                    folium.CircleMarker(
                        location=location,
                        radius=14,
                        fill_color=icon_color,
                        color='black',
                        fill_opacity=0.7,
                        weight=2,
                        popup=f"""
                        <div style="min-width: 180px; padding: 10px;">
                            <h4 style="margin-top: 0; color: {icon_color};">Delivery #{delivery_num}</h4>
                            <p><b>Node ID:</b> {node_id}</p>
                            <p><b>Distance from previous:</b> {segment_distances[i-1]/1000:.2f} km</p>
                            <p><b>Warehouse:</b> {warehouse_name}</p>
                        </div>
                        """,
                        tooltip=f"Delivery #{delivery_num}"
                    ).add_to(m)
                    
                    # Add number label
                    folium.map.Marker(
                        location,
                        icon=folium.DivIcon(
                            icon_size=(20, 20),
                            icon_anchor=(10, 10),
                            html=f'<div style="font-size: 10pt; color: white; text-align: center; line-height: 0px;">{delivery_num}</div>'
                        )
                    ).add_to(m)
            except KeyError:
                print(f"Warning: Node {node_id} not found in graph")
    
    # Process alternative paths
    if alternative_paths:
        for warehouse_id, alt_paths in alternative_paths.items():
            # Find matching warehouse color
            w_idx = None
            warehouse_name = None
            for i, wid in enumerate(warehouse_routes.keys()):
                if wid == warehouse_id:
                    w_idx = i
                    warehouse_name = warehouse_names.get(warehouse_id)
                    break
            
            if w_idx is None or warehouse_name is None:
                continue
                
            base_color = warehouse_base_colors[w_idx]
            
            for alt_idx, node_sequence in enumerate(alt_paths):
                alt_group = folium.FeatureGroup(name=f"Alternative Route {alt_idx+1}: {warehouse_name}").add_to(m)
                
                # Create visually different color for alternative route
                h, s, v = colorsys.rgb_to_hsv(
                    int(base_color[1:3], 16)/255, 
                    int(base_color[3:5], 16)/255, 
                    int(base_color[5:7], 16)/255
                )
                # Shift hue slightly, lower saturation
                h = (h + 0.1) % 1.0
                s = max(0.4, s - 0.3)
                r, g, b = colorsys.hsv_to_rgb(h, s, v)
                alt_color = f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"
                
                route_segments = []
                segment_distances = []
                
                # Calculate segments and distances
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
                                    
                                    if 'geometry' in data and data['geometry']:
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
                
                # Calculate total distance
                total_distance = sum(segment_distances)
                key = f"{warehouse_name} (Alt {alt_idx+1})"
                alt_distances[key] = total_distance
                
                # Add segments to map
                for idx, (segment, segment_distance) in enumerate(zip(route_segments, segment_distances)):
                    tooltip_html = f"""
                    <div style="font-family: Arial, sans-serif; padding: 8px; min-width: 200px;">
                        <div style="display: flex; align-items: center; margin-bottom: 5px;">
                            <div style="width: 12px; height: 12px; background-color: {alt_color}; 
                                        margin-right: 8px;"></div>
                            <h4 style="margin: 0; color: {alt_color}; font-size: 14px;">
                                {warehouse_name} - ALTERNATIVE ROUTE {alt_idx+1}
                            </h4>
                        </div>
                        <div style="border-top: 1px solid #eee; padding-top: 5px;">
                            <div style="display: flex; justify-content: space-between;">
                                <span style="font-weight: bold;">Delivery Segment:</span>
                                <span>{idx+1}/{len(route_segments)}</span>
                            </div>
                            <div style="display: flex; justify-content: space-between;
                                        margin-top: 5px;">
                                <span style="font-weight: bold;">Segment Distance:</span>
                                <span>{segment_distance/1000:.2f} km</span>
                            </div>
                            <div style="margin-top: 5px;">
                                <span style="font-weight: bold;">Total Route Distance:</span>
                                <span>{total_distance/1000:.2f} km</span>
                            </div>
                        </div>
                    </div>
                    """
                    
                    # Add alternative path with dashed line
                    folium.PolyLine(
                        locations=segment,
                        color=alt_color,
                        weight=4,
                        opacity=0.7,
                        tooltip=folium.Tooltip(tooltip_html),
                        dash_array='5,10'  # Dashed line for alternatives
                    ).add_to(alt_group)
                
                # Add delivery markers for alternative route with smaller, different style
                for i, node_id in enumerate(node_sequence):
                    if i == 0 or i == len(node_sequence) - 1:
                        # Skip warehouse nodes (already added)
                        continue
                        
                    try:
                        node = G.nodes[str(node_id)]
                        alt_location = [float(node['y']), float(node['x'])]
                        
                        # Add different styled marker for alternative route deliveries
                        folium.CircleMarker(
                            location=alt_location,
                            radius=8,
                            fill_color=alt_color,
                            color='black',
                            fill_opacity=0.5,
                            weight=1,
                            popup=f"""
                            <div style="min-width: 180px; padding: 10px;">
                                <h4 style="margin-top: 0; color: {alt_color};">Alt Route {alt_idx+1} - Stop #{i}</h4>
                                <p><b>Node ID:</b> {node_id}</p>
                                <p><b>Warehouse:</b> {warehouse_name}</p>
                            </div>
                            """,
                            tooltip=f"Alt Route {alt_idx+1} - Stop #{i}"
                        ).add_to(m)
                    except KeyError:
                        pass
    
    # Add legend with total distances and route visualization explanation
    legend_html = f"""
    <div style="position: fixed; 
                bottom: 50px; left: 50px; 
                border: 2px solid grey; z-index: 9999; 
                background-color: white; padding: 10px;
                border-radius: 5px; box-shadow: 0 0 15px rgba(0, 0, 0, 0.2);
                max-width: 300px;">
        <h4 style="margin-top: 0; margin-bottom: 10px;">Route Distances</h4>
        <div style="font-weight: bold; text-decoration: underline; margin-bottom: 5px;">Optimal Routes:</div>
        <table style="width: 100%;">
            {"".join(f"""
            <tr>
                <td style="padding: 2px 10px 2px 0; text-align: left;">{name}</td>
                <td style="padding: 2px 0; text-align: right;">{dist/1000:.2f} km</td>
            </tr>
            """ for name, dist in warehouse_distances.items())}
        </table>
        
        {"<div style='font-weight: bold; text-decoration: underline; margin: 10px 0 5px 0;'>Alternative Routes:</div>" if alt_distances else ""}
        {"<table style='width: 100%;'>" + "".join(f"""
            <tr>
                <td style="padding: 2px 10px 2px 0; text-align: left;">{name}</td>
                <td style="padding: 2px 0; text-align: right;">{dist/1000:.2f} km</td>
            </tr>
            """ for name, dist in alt_distances.items()) + "</table>" if alt_distances else ""}
    </div>
    """
    
    m.get_root().html.add_child(folium.Element(legend_html))
    
    # Add layer control to toggle
    folium.LayerControl(position='bottomright').add_to(m)
    
    # Save the map to file
    m.save(output_file)
    print(f"Map saved to {output_file}")
    
    return m

def main():
    """Main function to process data and generate visualizations"""
    # Set file paths
    data_dir = "../../data"
    nodes_file = os.path.join(data_dir, "nodes.csv")
    edges_file = os.path.join(data_dir, "edges.csv")
    output_file = "delhi_delivery_routes.html"
    
    # Ensure files exist
    if not os.path.exists(nodes_file) or not os.path.exists(edges_file):
        print(f"Error: Required data files not found. Please ensure that nodes.csv and edges.csv are in {data_dir}")
        return
    
    print("Loading graph from CSV files...")
    G = load_osm_graph(nodes_file, edges_file)
    print(f"Graph loaded with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
    
    # Process optimal route files
    print("Looking for route files...")
    optimal_routes, warehouse_names = process_route_files("warehouse_*_optimal_route.txt")
    
    if not optimal_routes:
        print("No route files found. Please run the C++ program first to generate routes.")
        return
    
    # Process alternative path files
    print("Looking for alternative path files...")
    alternative_paths = process_alternative_path_files("warehouse_*_alternative_route_*.txt")
    
    # Generate map
    print("Generating map visualization...")
    map_object = plot_multiple_routes(
        G,
        optimal_routes,
        warehouse_names,
        alternative_paths,
        output_file
    )
    
    # Open the map in the browser
    print("Opening map in browser...")
    webbrowser.open('file://' + os.path.realpath(output_file))

if __name__ == "__main__":
    main()