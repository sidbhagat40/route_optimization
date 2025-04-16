#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <queue>
#include <cmath>
#include <limits>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <random>
#include <unordered_set>
#include <ctime>
#include <map>

using namespace std;

// Node structure
struct Node {
    double x, y;
    string id;
};

// Edge structure
struct Edge {
    string source, target;
    double length;
    bool oneway;
};

// Warehouse structure
struct Warehouse {
    string id;
    string name;
    vector<string> assignedDeliveries;
};

class RoadNetwork {
private:
    unordered_map<string, Node> nodes;
    unordered_map<string, vector<pair<string, double>>> adjList;
    
public:
    void addNode(const string& id, double x, double y) {
        nodes[id] = {x, y, id};
    }

    void addEdge(const string& u, const string& v, double length, bool oneway) {
        adjList[u].emplace_back(v, length);
        if (!oneway) {
            adjList[v].emplace_back(u, length);
        }
    }

    bool hasNode(const string& id) const {
        return nodes.find(id) != nodes.end();
    }

    const Node& getNode(const string& id) const {
        return nodes.at(id);
    }

    // Calculate straight-line distance between two nodes
    double calculateHaversineDistance(const string& node1, const string& node2) const {
        if (nodes.find(node1) == nodes.end() || nodes.find(node2) == nodes.end()) {
            return numeric_limits<double>::infinity();
        }

        // Get coordinates
        double lat1 = nodes.at(node1).y;
        double lon1 = nodes.at(node1).x;
        double lat2 = nodes.at(node2).y;
        double lon2 = nodes.at(node2).x;

        // Convert to radians
        lat1 *= M_PI / 180.0;
        lon1 *= M_PI / 180.0;
        lat2 *= M_PI / 180.0;
        lon2 *= M_PI / 180.0;

        // Haversine formula
        double dlon = lon2 - lon1;
        double dlat = lat2 - lat1;
        double a = pow(sin(dlat/2), 2) + cos(lat1) * cos(lat2) * pow(sin(dlon/2), 2);
        double c = 2 * asin(sqrt(a));
        double r = 6371000; // Earth radius in meters
        
        return r * c;
    }

    unordered_map<string, double> dijkstra(const string& start) {
        unordered_map<string, double> distances;
        for (const auto& node : nodes) {
            distances[node.first] = numeric_limits<double>::infinity();
        }
        distances[start] = 0.0;

        priority_queue<pair<double, string>, 
                     vector<pair<double, string>>,
                     greater<pair<double, string>>> pq;
        pq.push({0.0, start});

        while (!pq.empty()) {
            auto current = pq.top();
            pq.pop();
            double currentDist = current.first;
            string currentNode = current.second;

            if (currentDist > distances[currentNode]) continue;

            if (adjList.find(currentNode) == adjList.end()) continue;

            for (const auto& neighbor : adjList.at(currentNode)) {
                string nextNode = neighbor.first;
                double edgeLength = neighbor.second;
                double newDist = currentDist + edgeLength;

                if (newDist < distances[nextNode]) {
                    distances[nextNode] = newDist;
                    pq.push({newDist, nextNode});
                }
            }
        }

        return distances;
    }

    vector<vector<double>> createDistanceMatrix(const vector<string>& deliveryNodes) {
        size_t n = deliveryNodes.size();
        vector<vector<double>> matrix(n, vector<double>(n, -1));

        for (size_t i = 0; i < n; ++i) {
            cout << "Calculating routes from node " << deliveryNodes[i] << " (" << i+1 << "/" << n << ")..." << endl;
            auto distances = dijkstra(deliveryNodes[i]);
            
            for (size_t j = 0; j < n; ++j) {
                if (i == j) {
                    matrix[i][j] = 0;
                } else {
                    if (distances.find(deliveryNodes[j]) != distances.end()) {
                        matrix[i][j] = distances[deliveryNodes[j]];
                    } else {
                        matrix[i][j] = numeric_limits<double>::infinity();
                    }
                }
            }
        }

        return matrix;
    }

    vector<string> selectRandomDeliveryNodes(int count, const vector<string>& excludeNodes = {}) {
        vector<string> allNodes;
        unordered_set<string> excludeSet(excludeNodes.begin(), excludeNodes.end());
        
        for (const auto& node : nodes) {
            if (excludeSet.find(node.first) == excludeSet.end()) {
                allNodes.push_back(node.first);
            }
        }
        
        if (allNodes.empty()) {
            cerr << "Error: No nodes available in the graph!" << endl;
            return {};
        }
        
        count = min(count, (int)allNodes.size());
        
        // Reservoir sampling algorithm for random selection
        vector<string> result(count);
        random_device rd;
        mt19937 gen(rd());
        
        // Fill the reservoir
        for (int i = 0; i < count; i++) {
            result[i] = allNodes[i];
        }
        
        // Replace elements with gradually decreasing probability
        for (int i = count; i < allNodes.size(); i++) {
            uniform_int_distribution<> dis(0, i);
            int j = dis(gen);
            if (j < count) {
                result[j] = allNodes[i];
            }
        }
        
        return result;
    }

    void printGraphStats() const {
        cout << "Graph Statistics:" << endl;
        cout << "Number of nodes: " << nodes.size() << endl;
        cout << "Number of edges: " << adjList.size() << endl;
        
        if (!nodes.empty()) {
            cout << "Sample node: " << nodes.begin()->first 
                 << " (" << nodes.begin()->second.x 
                 << ", " << nodes.begin()->second.y << ")" << endl;
        }
        
        if (!adjList.empty()) {
            cout << "Sample edge from: " << adjList.begin()->first 
                 << " to " << adjList.begin()->second[0].first 
                 << " length " << adjList.begin()->second[0].second << endl;
        }
    }

    // Assign deliveries to nearest warehouse based on Dijkstra distances
    map<string, vector<string>> assignDeliveriesToWarehouses(
        const vector<string>& warehouseIds, 
        const vector<string>& deliveryNodes
    ) {
        map<string, vector<string>> assignments;
        
        cout << "Assigning deliveries to warehouses..." << endl;
        
        // Initialize assignments
        for (const auto& warehouseId : warehouseIds) {
            assignments[warehouseId] = {};
        }
        
        // For each delivery node, find the closest warehouse
        for (const auto& deliveryNode : deliveryNodes) {
            string closestWarehouse;
            double minDistance = numeric_limits<double>::infinity();
            
            // Calculate distances to all warehouses using Dijkstra
            // First try network distance
            for (const auto& warehouseId : warehouseIds) {
                auto distances = dijkstra(warehouseId);
                
                // Check if there's a path to this delivery node
                if (distances.find(deliveryNode) != distances.end() && 
                    distances[deliveryNode] < minDistance) {
                    minDistance = distances[deliveryNode];
                    closestWarehouse = warehouseId;
                }
            }
            
            // If no path found with Dijkstra, use straight-line distance
            if (minDistance == numeric_limits<double>::infinity()) {
                cout << "Warning: No path found to delivery node " << deliveryNode 
                     << " using Dijkstra. Using straight-line distance instead." << endl;
                
                for (const auto& warehouseId : warehouseIds) {
                    double dist = calculateHaversineDistance(warehouseId, deliveryNode);
                    if (dist < minDistance) {
                        minDistance = dist;
                        closestWarehouse = warehouseId;
                    }
                }
            }
            
            // Assign to the closest warehouse
            if (!closestWarehouse.empty()) {
                assignments[closestWarehouse].push_back(deliveryNode);
                cout << "Assigned delivery " << deliveryNode << " to warehouse " 
                     << closestWarehouse << " (distance: " << minDistance << "m)" << endl;
            } else {
                cerr << "Warning: Could not assign delivery " << deliveryNode 
                     << " to any warehouse!" << endl;
            }
        }
        
        return assignments;
    }
};

void loadNodes(RoadNetwork& graph, const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening nodes file: " << filename << endl;
        return;
    }

    string line;
    getline(file, line); // Skip header
    int count = 0;

    while (getline(file, line)) {
        stringstream ss(line);
        string osmid, y, x;
        
        getline(ss, osmid, ',');
        getline(ss, y, ',');
        getline(ss, x, ',');
        
        try {
            graph.addNode(osmid, stod(x), stod(y));
            count++;
        } catch (...) {
            cerr << "Warning: Failed to parse node line: " << line << endl;
            continue;
        }
    }
    cout << "Successfully loaded " << count << " nodes" << endl;
}

void loadEdges(RoadNetwork& graph, const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening edges file: " << filename << endl;
        return;
    }

    string line;
    getline(file, line); // Skip header
    int count = 0;

    while (getline(file, line)) {
        stringstream ss(line);
        vector<string> tokens;
        string token;
        
        bool inQuotes = false;
        string field;
        for (char c : line) {
            if (c == '"') {
                inQuotes = !inQuotes;
            } else if (c == ',' && !inQuotes) {
                tokens.push_back(field);
                field.clear();
            } else {
                field += c;
            }
        }
        tokens.push_back(field);
        
        if (tokens.size() < 9) {
            cerr << "Warning: Insufficient columns in edge line: " << line << endl;
            continue;
        }
        
        try {
            string u = tokens[0];
            string v = tokens[1];
            double length = stod(tokens[8]);
            bool oneway = (tokens[6] == "True" || tokens[6] == "TRUE");
            
            graph.addEdge(u, v, length, oneway);
            count++;
        } catch (...) {
            cerr << "Warning: Failed to parse edge line: " << line << endl;
            continue;
        }
    }
    cout << "Successfully loaded " << count << " edges" << endl;
}

vector<int> nearestNeighborTSP(const vector<vector<double>>& distMatrix) {
    int n = distMatrix.size();
    
    if (n <= 1) {
        return {0};
    }
    
    vector<int> tour = {0};
    vector<bool> visited(n, false);
    visited[0] = true;
    
    for (int i = 1; i < n; i++) {
        int closest = -1;
        double minDist = numeric_limits<double>::max();
        
        for (int j = 0; j < n; j++) {
            if (!visited[j] && distMatrix[tour.back()][j] < minDist && 
                distMatrix[tour.back()][j] != numeric_limits<double>::infinity()) {
                minDist = distMatrix[tour.back()][j];
                closest = j;
            }
        }
        
        if (closest == -1) {
            // If we can't find a reachable node, try to find any unvisited node
            for (int j = 0; j < n; j++) {
                if (!visited[j]) {
                    closest = j;
                    break;
                }
            }
        }
        
        if (closest == -1) break; // All nodes visited
        
        tour.push_back(closest);
        visited[closest] = true;
    }
    
    // Try to complete the cycle back to the start
    if (tour.size() > 1 && distMatrix[tour.back()][0] != numeric_limits<double>::infinity()) {
        tour.push_back(0);
    }
    
    return tour;
}

void twoOptOptimize(vector<int>& tour, const vector<vector<double>>& distMatrix) {
    // Skip optimization for very small tours
    if (tour.size() <= 3) return;
    
    bool improved = true;
    int iterCount = 0;
    const int MAX_ITER = 1000; // Limit iterations to prevent infinite loops
    
    while (improved && iterCount < MAX_ITER) {
        improved = false;
        iterCount++;
        
        for (int i = 1; i < tour.size()-2; i++) {
            for (int j = i+1; j < tour.size()-1; j++) {
                // Skip if any distances are infinity
                if (distMatrix[tour[i-1]][tour[j]] == numeric_limits<double>::infinity() ||
                    distMatrix[tour[i]][tour[j+1]] == numeric_limits<double>::infinity() ||
                    distMatrix[tour[i-1]][tour[i]] == numeric_limits<double>::infinity() ||
                    distMatrix[tour[j]][tour[j+1]] == numeric_limits<double>::infinity()) {
                    continue;
                }
                
                double delta = distMatrix[tour[i-1]][tour[j]] 
                            + distMatrix[tour[i]][tour[j+1]]
                            - distMatrix[tour[i-1]][tour[i]] 
                            - distMatrix[tour[j]][tour[j+1]];
                            
                if (delta < -0.001) {
                    reverse(tour.begin()+i, tour.begin()+j+1);
                    improved = true;
                }
            }
        }
    }
    
    if (iterCount >= MAX_ITER) {
        cout << "Warning: Two-opt optimization reached maximum iterations" << endl;
    }
}

double calculateTourDistance(const vector<int>& tour, 
                           const vector<vector<double>>& distMatrix) {
    double total = 0;
    for (int i = 1; i < tour.size(); i++) {
        if (distMatrix[tour[i-1]][tour[i]] != numeric_limits<double>::infinity()) {
            total += distMatrix[tour[i-1]][tour[i]];
        } else {
            // If there's no path between consecutive points
            cerr << "Warning: No path between tour points " 
                 << tour[i-1] << " and " << tour[i] << endl;
        }
    }
    return total;
}

vector<int> solveTSP(const vector<vector<double>>& distMatrix) {
    if (distMatrix.empty()) {
        return {};
    }
    
    vector<int> tour = nearestNeighborTSP(distMatrix);
    
    if (tour.size() <= 1) {
        return tour;
    }
    
    double initialDist = calculateTourDistance(tour, distMatrix);
    
    twoOptOptimize(tour, distMatrix);
    double optimizedDist = calculateTourDistance(tour, distMatrix);
    
    cout << "Initial tour distance: " << initialDist << " meters\n";
    cout << "Optimized tour distance: " << optimizedDist << " meters\n";
    cout << "Improvement: " << (initialDist - optimizedDist) << " meters (" 
         << ((initialDist - optimizedDist)/initialDist*100) << "%)\n";
    
    return tour;
}

int main() {
    // Seed the random number generator
    srand(time(0));
    
    RoadNetwork graph;
    
    cout << "Loading nodes..." << endl;
    loadNodes(graph, "../../data/nodes.csv");
    
    cout << "Loading edges..." << endl;
    loadEdges(graph, "../../data/edges.csv");
    
    graph.printGraphStats();
    
    // Define warehouse nodes for Delhi
    // These are example OSM IDs - replace with actual IDs from your Delhi dataset
    vector<string> warehouseIds = {
        "12110005374",  // Warehouse 1 (North Delhi)
        "923637303",  // Warehouse 2 (East Delhi)
        "4314376797"   // Warehouse 3 (South Delhi)
    };
    
    // Verify warehouse nodes exist
    vector<string> validWarehouseIds;
    for (const auto& id : warehouseIds) {
        if (graph.hasNode(id)) {
            validWarehouseIds.push_back(id);
            cout << "Warehouse " << id << " is valid." << endl;
        } else {
            cerr << "Warning: Warehouse " << id << " not found in the graph!" << endl;
        }
    }
    
    if (validWarehouseIds.empty()) {
        cerr << "Error: No valid warehouses found. Exiting." << endl;
        return 1;
    }
    
    // Ask user for the number of delivery points
    int numDeliveries;
    cout << "Enter the number of deliveries to make: ";
    cin >> numDeliveries;
    
    if (numDeliveries <= 0) {
        cerr << "Error: Number of deliveries must be positive. Exiting." << endl;
        return 1;
    }
    
    cout << "Selecting " << numDeliveries << " random delivery locations..." << endl;
    vector<string> deliveryNodes = graph.selectRandomDeliveryNodes(numDeliveries, validWarehouseIds);
    
    if (deliveryNodes.empty()) {
        cerr << "Error: No delivery nodes selected. Exiting." << endl;
        return 1;
    }
    
    cout << "\nSelected " << deliveryNodes.size() << " delivery nodes." << endl;
    
    // Assign deliveries to warehouses
    auto assignments = graph.assignDeliveriesToWarehouses(validWarehouseIds, deliveryNodes);
    
    // Process each warehouse's deliveries
    for (size_t w = 0; w < validWarehouseIds.size(); ++w) {
        string warehouseId = validWarehouseIds[w];
        vector<string> warehouseDeliveries = assignments[warehouseId];
        
        cout << "\n--- Processing Warehouse " << warehouseId << " ---" << endl;
        cout << "Assigned " << warehouseDeliveries.size() << " deliveries" << endl;
        
        if (warehouseDeliveries.empty()) {
            cout << "No deliveries assigned to this warehouse. Skipping." << endl;
            continue;
        }
        
        // Create a list with warehouse first, then deliveries
        vector<string> routeNodes = {warehouseId};
        routeNodes.insert(routeNodes.end(), warehouseDeliveries.begin(), warehouseDeliveries.end());
        
        cout << "\nCalculating distance matrix for warehouse " << warehouseId << "..." << endl;
        auto distanceMatrix = graph.createDistanceMatrix(routeNodes);
        
        cout << "\nSolving TSP for warehouse " << warehouseId << "..." << endl;
        vector<int> optimalTour = solveTSP(distanceMatrix);
        
        cout << "\nOptimal Delivery Route for Warehouse " << warehouseId << ":\n";
        for (size_t i = 0; i < optimalTour.size(); i++) {
            int nodeIndex = optimalTour[i];
            cout << i << ": Node " << routeNodes[nodeIndex];
            if (i > 0) {
                cout << " | Distance from previous: " 
                     << distanceMatrix[optimalTour[i-1]][nodeIndex] << "m";
            }
            cout << endl;
        }
        
        cout << "\nTotal route distance: " 
             << calculateTourDistance(optimalTour, distanceMatrix) 
             << " meters" << endl;

        // Create output filename for this warehouse's route
        string outputFile = "route_warehouse_" + to_string(w + 1) + ".txt";
        ofstream outFile(outputFile);
        if (!outFile.is_open()) {
            cerr << "Error: Could not open " << outputFile << " for writing!" << endl;
            continue;
        }

        // Save the route with proper OSM IDs
        for (int idx : optimalTour) {
            outFile << routeNodes[idx] << endl;
        }
        outFile.close();

        cout << "\nRoute successfully saved to " << outputFile << endl;
    }

        cout << "\nRunning visualization script..." << endl;
        string command = "python ../python/plot.py";
        int result = system(command.c_str());
        if (result != 0) {
            cerr << "Error: Failed to run visualization script" << endl;
        }
    
    return 0;
}