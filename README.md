# Delivery Route Optimization System

## Overview

This project is a C++ and Python-based solution for optimizing delivery routes from multiple warehouses to a set of destinations using graph algorithms and interactive visualization. It enables efficient path computation and clear spatial insight using OpenStreetMap (OSM) data.

C++ handles route optimization logic using Dijkstra’s algorithm and writes results to structured files.

Python reads the computed paths, builds the city graph, and visualizes the optimal and alternative delivery routes on an interactive web map using Folium.

## Features

**Multi-Warehouse Support :** Calculates separate optimal routes for each warehouse.

**Interactive Map Visualization :** Renders delivery paths and markers with hoverable tooltips, distance data, and custom color coding.

**Alternative Path Support :** Displays backup delivery routes with visual distinction.

**Custom OSM Input :** Graph created using OSM node and edge data in CSV format.

## Dependencies

Python

pandas

networkx

folium

colorsys

**Install all dependencies using:**

```bash
pip install pandas networkx folium
```
**C++**

Standard C++ compiler supporting C++11 or above

## Usage

**Step 1 :**  Generate Routes (C++)
- Go to directory

```bash
route_optimization\algorithm\cpp>
```
- Run the command
```bash
g++ main.cpp && .\a.exe
```
This generates route files in the current directory or specified outputs/ folder and generates delhi_delivery_routes.html and opens it in your default web browser.


##  Map Output

The final visualization includes:

- Total route distance and delivery counts per warehouse.

- Color-coded markers for delivery stops.

- Toggleable layers for each route.

