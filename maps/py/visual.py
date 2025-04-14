import pandas as pd
import folium

# Load nodes from CSV
nodes_df = pd.read_csv("nodes.csv")

# Your OSM node path
osm_node_path = [
    1827697182, 1826909611, 4233302722, 9884874655, 4229518792,
    4065531426, 9880357040, 12110038882, 2575648598, 920829184, 1827697182
]

# Create a dictionary for quick lookup
id_to_coords = {
    row['id']: (row['lat'], row['lon']) for _, row in nodes_df.iterrows()
}

# Build path coordinates from node IDs
path_coords = []
for node_id in osm_node_path:
    if node_id in id_to_coords:
        path_coords.append(id_to_coords[node_id])
    else:
        print(f"Node ID {node_id} not found in CSV")

# Create map centered on the first node
m = folium.Map(location=path_coords[0], zoom_start=13)

# Add markers
for i, (lat, lon) in enumerate(path_coords):
    folium.Marker([lat, lon], popup=f"Node {i}").add_to(m)

# Add path polyline
folium.PolyLine(path_coords, color="blue", weight=3, opacity=0.8).add_to(m)

# Save to HTML
m.save("route_map_local.html")
