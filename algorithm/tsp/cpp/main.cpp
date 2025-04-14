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

class RoadNetwork {
private:
    unordered_map<string, Node> nodes;
    unordered_map<string, vector<pair<string, double>>> adjList;
    
public:
    void addNode(const string& id, double x, double y);
    void addEdge(const string& u, const string& v, double length, bool oneway);
    vector<vector<double>> createDistanceMatrix(const vector<string>& deliveryNodes);
    unordered_map<string, double> dijkstra(const string& start);
    vector<string> selectRandomDeliveryNodes(int count);
    void printGraphStats() const;
};

void RoadNetwork::addNode(const string& id, double x, double y) {
    nodes[id] = {x, y, id};
}

void RoadNetwork::addEdge(const string& u, const string& v, double length, bool oneway) {
    adjList[u].emplace_back(v, length);
    if (!oneway) {
        adjList[v].emplace_back(u, length);
    }
}

unordered_map<string, double> RoadNetwork::dijkstra(const string& start) {
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

        // Check if node exists in adjacency list
        if (adjList.find(currentNode) == adjList.end()) continue;

        for (const auto& neighbor : adjList.at(currentNode)){
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

vector<vector<double>> RoadNetwork::createDistanceMatrix(const vector<string>& deliveryNodes) {
    size_t n = deliveryNodes.size();
    vector<vector<double>> matrix(n, vector<double>(n, -1)); // Initialize with -1

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

vector<string> RoadNetwork::selectRandomDeliveryNodes(int count) {
    vector<string> allNodes;
    for (const auto& node : nodes) {
        allNodes.push_back(node.first);
    }
    
    if (allNodes.empty()) {
        cerr << "Error: No nodes available in the graph!" << endl;
        return {};
    }
    
    count = min(count, (int)allNodes.size());
    
    // Better random shuffling
    random_device rd;
    mt19937 g(rd());
    shuffle(allNodes.begin(), allNodes.end(), g);
    
    return vector<string>(allNodes.begin(), allNodes.begin() + count);
}

void RoadNetwork::printGraphStats() const {
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
        
        // Comma-separated format
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
        
        // Split by commas but handle quoted fields
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
        tokens.push_back(field); // Add last field
        
        if (tokens.size() < 9) {
            cerr << "Warning: Insufficient columns in edge line: " << line << endl;
            continue;
        }
        
        try {
            string u = tokens[0];
            string v = tokens[1];
            double length = stod(tokens[8]); // length is 9th column (0-indexed 8)
            bool oneway = (tokens[6] == "True" || tokens[6] == "TRUE"); // case-insensitive check
            
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
    vector<int> tour = {0}; // Start at first location
    vector<bool> visited(n, false);
    visited[0] = true;
    
    for (int i = 1; i < n; i++) {
        int closest = -1;
        double minDist = numeric_limits<double>::max();
        
        // Find nearest unvisited neighbor
        for (int j = 0; j < n; j++) {
            if (!visited[j] && distMatrix[tour.back()][j] < minDist) {
                minDist = distMatrix[tour.back()][j];
                closest = j;
            }
        }
        
        if (closest == -1) break; // No unvisited nodes left
        tour.push_back(closest);
        visited[closest] = true;
    }
    
    // Return to starting point
    tour.push_back(0);
    return tour;
}

// 2. 2-Opt Local Search Optimization
void twoOptOptimize(vector<int>& tour, const vector<vector<double>>& distMatrix) {
    bool improved = true;
    while (improved) {
        improved = false;
        for (int i = 1; i < tour.size()-2; i++) {
            for (int j = i+1; j < tour.size()-1; j++) {
                double delta = distMatrix[tour[i-1]][tour[j]] 
                            + distMatrix[tour[i]][tour[j+1]]
                            - distMatrix[tour[i-1]][tour[i]] 
                            - distMatrix[tour[j]][tour[j+1]];
                if (delta < -0.001) { // Found improvement
                    reverse(tour.begin()+i, tour.begin()+j+1);
                    improved = true;
                }
            }
        }
    }
}

// 3. Calculate Total Tour Distance
double calculateTourDistance(const vector<int>& tour, 
                           const vector<vector<double>>& distMatrix) {
    double total = 0;
    for (int i = 1; i < tour.size(); i++) {
        total += distMatrix[tour[i-1]][tour[i]];
    }
    return total;
}

// 4. Main TSP Solver Function
vector<int> solveTSP(const vector<vector<double>>& distMatrix) {
    // Step 1: Get initial solution
    vector<int> tour = nearestNeighborTSP(distMatrix);
    double initialDist = calculateTourDistance(tour, distMatrix);
    
    // Step 2: Optimize with 2-opt
    twoOptOptimize(tour, distMatrix);
    double optimizedDist = calculateTourDistance(tour, distMatrix);
    
    cout << "Initial tour distance: " << initialDist << " meters\n";
    cout << "Optimized tour distance: " << optimizedDist << " meters\n";
    cout << "Improvement: " << (initialDist - optimizedDist) << " meters (" 
         << ((initialDist - optimizedDist)/initialDist*100) << "%)\n";
    
    return tour;
}

int main() {
    RoadNetwork graph;
    
    // Load data from files with verbose output
    cout << "Loading nodes..." << endl;
    loadNodes(graph, "nodes.csv");
    
    cout << "Loading edges..." << endl;
    loadEdges(graph, "edges.csv");
    
    // Print graph statistics for verification
    graph.printGraphStats();
    
    // Select delivery nodes
    cout << "Selecting delivery locations..." << endl;
    vector<string> deliveryNodes = graph.selectRandomDeliveryNodes(10);
    
    if (deliveryNodes.empty()) {
        cerr << "Error: No delivery nodes selected. Exiting." << endl;
        return 1;
    }
    
    cout << "\nSelected delivery nodes:" << endl;
    for (size_t i = 0; i < deliveryNodes.size(); ++i) {
        cout << i << ": " << deliveryNodes[i] << endl;
    }
    
    // Create distance matrix
    cout << "\nCalculating distance matrix..." << endl;
    auto distanceMatrix = graph.createDistanceMatrix(deliveryNodes);
    
    // Print distance matrix with better formatting
    cout << "\nDistance Matrix (in meters):" << endl;
    cout << fixed << setprecision(2);
    
    // Column headers
    cout << setw(10) << " ";
    for (size_t j = 0; j < deliveryNodes.size(); ++j) {
        cout << setw(10) << j;
    }
    cout << endl;
    
    // Matrix rows
    for (size_t i = 0; i < distanceMatrix.size(); ++i) {
        cout << setw(4) << i << " -> ";
        for (size_t j = 0; j < distanceMatrix[i].size(); ++j) {
            if (distanceMatrix[i][j] == -1) {
                cout << setw(10) << "N/A";
            } else if (isinf(distanceMatrix[i][j])) {
                cout << setw(10) << "INF";
            } else {
                cout << setw(10) << distanceMatrix[i][j];
            }
        }
        cout << endl;
    }

    vector<int> optimalTour = solveTSP(distanceMatrix);
    
    // Print optimal tour
    cout << "\nOptimal Delivery Route:\n";
    for (size_t i = 0; i < optimalTour.size(); i++) {
        int nodeIndex = optimalTour[i];
        cout << i << ": Node " << deliveryNodes[nodeIndex] 
             << " (ID: " << nodeIndex << ")";
        if (i > 0) {
            cout << " | Distance from previous: " 
                 << distanceMatrix[optimalTour[i-1]][nodeIndex] << "m";
        }
        cout << endl;
    }
    
    // Print total distance
    cout << "\nTotal route distance: " 
         << calculateTourDistance(optimalTour, distanceMatrix) 
         << " meters" << endl;

    
    
    return 0;
}