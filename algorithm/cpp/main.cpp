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
#include <chrono>

using namespace std;

// Forward declarations
double calculateTourDistance(const vector<int>& tour, const vector<vector<double>>& distMatrix);
vector<int> solveTSP(const vector<vector<double>>& distMatrix);

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

struct WarehouseDistInfo {
    string warehouseId;
    double distance;
    vector<string> path;
    bool isReachable;
};

struct PathOption {
    string warehouseId;
    vector<string> path;
    double distance;
};

class RoadNetwork {
private:
    unordered_map<string, Node> nodes;
    unordered_map<string, vector<pair<string, double>>> adjList;
    unordered_map<string, unordered_map<string, double>> distanceCache;
    unordered_map<string, unordered_map<string, string>> pathCache;

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

    double calculateHaversineDistance(const string& node1, const string& node2) const {
        if (nodes.find(node1) == nodes.end() || nodes.find(node2) == nodes.end()) {
            return numeric_limits<double>::infinity();
        }

        const auto& n1 = nodes.at(node1);
        const auto& n2 = nodes.at(node2);
        
        double lat1 = n1.y * M_PI / 180.0;
        double lon1 = n1.x * M_PI / 180.0;
        double lat2 = n2.y * M_PI / 180.0;
        double lon2 = n2.x * M_PI / 180.0;

        double dlon = lon2 - lon1;
        double dlat = lat2 - lat1;
        double a = pow(sin(dlat/2), 2) + cos(lat1) * cos(lat2) * pow(sin(dlon/2), 2);
        return 6371000 * 2 * asin(sqrt(a));
    }

    pair<unordered_map<string, double>, unordered_map<string, string>> dijkstraWithPath(const string& start) {
        if (distanceCache.find(start) != distanceCache.end()) {
            return make_pair(distanceCache[start], pathCache[start]);
        }

        unordered_map<string, double> distances;
        unordered_map<string, string> previousNode;
        
        for (const auto& node : nodes) {
            distances[node.first] = numeric_limits<double>::infinity();
        }
        distances[start] = 0.0;

        priority_queue<pair<double, string>, vector<pair<double, string>>, greater<pair<double, string>>> pq;
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
                    previousNode[nextNode] = currentNode;
                    pq.push({newDist, nextNode});
                }
            }
        }

        distanceCache[start] = distances;
        pathCache[start] = previousNode;
        return make_pair(distances, previousNode);
    }

    vector<string> getShortestPath(const string& start, const string& end) {
        if (pathCache.find(start) == pathCache.end()) {
            dijkstraWithPath(start);
        }
        
        vector<string> path;
        string current = end;
        
        while (current != start) {
            path.push_back(current);
            if (pathCache[start].find(current) == pathCache[start].end()) {
                return {};
            }
            current = pathCache[start].at(current);
        }
        path.push_back(start);
        reverse(path.begin(), path.end());
        return path;
    }

    vector<vector<double>> createDistanceMatrix(const vector<string>& deliveryNodes) {
        size_t n = deliveryNodes.size();
        vector<vector<double>> matrix(n, vector<double>(n, -1));

        #pragma omp parallel for
        for (size_t i = 0; i < n; ++i) {
            auto result = dijkstraWithPath(deliveryNodes[i]);
            auto& distances = result.first;
            
            for (size_t j = 0; j < n; ++j) {
                if (i == j) {
                    matrix[i][j] = 0;
                } else if (distances.find(deliveryNodes[j]) != distances.end()) {
                    matrix[i][j] = distances[deliveryNodes[j]];
                } else {
                    matrix[i][j] = numeric_limits<double>::infinity();
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
        
        random_device rd;
        unsigned int seed = rd() ^ static_cast<unsigned int>(chrono::system_clock::now().time_since_epoch().count());
        mt19937 gen(seed);
        
        shuffle(allNodes.begin(), allNodes.end(), gen);
        return vector<string>(allNodes.begin(), allNodes.begin() + count);
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
    }

    map<string, vector<string>> assignDeliveriesToWarehousesOptimized(
        const vector<string>& warehouseIds, 
        const vector<string>& deliveryNodes
    ) {
        map<string, vector<string>> assignments;
        
        cout << "\nOptimized delivery assignment..." << endl;
        
        for (const auto& warehouseId : warehouseIds) {
            assignments[warehouseId] = {};
        }
        
        for (const auto& deliveryNode : deliveryNodes) {
            WarehouseDistInfo bestWarehouse;
            bestWarehouse.isReachable = false;
            bestWarehouse.distance = numeric_limits<double>::infinity();
            
            for (const auto& warehouseId : warehouseIds) {
                auto result = dijkstraWithPath(warehouseId);
                auto distances = result.first;
                auto previous = result.second;
                
                if (distances.find(deliveryNode) != distances.end() && 
                    distances[deliveryNode] < numeric_limits<double>::infinity()) {
                    
                    double distance = distances[deliveryNode];
                    
                    if (distance < bestWarehouse.distance) {
                        bestWarehouse.warehouseId = warehouseId;
                        bestWarehouse.distance = distance;
                        bestWarehouse.path = getShortestPath(warehouseId, deliveryNode);
                        bestWarehouse.isReachable = true;
                    }
                }
            }
            
            if (bestWarehouse.isReachable) {
                assignments[bestWarehouse.warehouseId].push_back(deliveryNode);
            } else {
                string closestWarehouse;
                double minDistance = numeric_limits<double>::infinity();
                
                for (const auto& warehouseId : warehouseIds) {
                    double dist = calculateHaversineDistance(warehouseId, deliveryNode);
                    if (dist < minDistance) {
                        minDistance = dist;
                        closestWarehouse = warehouseId;
                    }
                }
                
                if (!closestWarehouse.empty()) {
                    assignments[closestWarehouse].push_back(deliveryNode);
                }
            }
        }
        
        cout << "\nAssignment Summary:" << endl;
        for (const auto& assignment : assignments) {
            cout << "Warehouse " << assignment.first << ": " << assignment.second.size() << " deliveries" << endl;
        }
        
        return assignments;
    }

    vector<PathOption> generateAlternativeRoutes(const string& warehouseId, 
        const vector<string>& deliveryNodes,
        const vector<int>& optimalTour,
        const vector<vector<double>>& distMatrix,
        int numAlternatives = 3) 
    {
        vector<PathOption> alternatives;

        PathOption optimalPath;
        optimalPath.warehouseId = warehouseId;
        optimalPath.distance = calculateTourDistance(optimalTour, distMatrix);
        for (int idx : optimalTour) {
            optimalPath.path.push_back(deliveryNodes[idx]);
        }
        alternatives.push_back(optimalPath);

        if (deliveryNodes.size() <= 3) {
            return alternatives;
        }

        random_device rd;
        mt19937 gen(rd());

        // Strategy 1: Random shuffle
        int attempts = 0;
        const int maxAttempts = 100;
        while (alternatives.size() < numAlternatives + 1 && attempts < maxAttempts) {
            attempts++;
            vector<int> permutedTour = optimalTour;
            
            if (permutedTour.size() > 3) {
                shuffle(permutedTour.begin() + 1, permutedTour.end() - 1, gen);
            }

            double distance = calculateTourDistance(permutedTour, distMatrix);
            if (distance > optimalPath.distance * 1.15) {
                PathOption altPath;
                altPath.warehouseId = warehouseId;
                altPath.distance = distance;
                for (int idx : permutedTour) {
                    altPath.path.push_back(deliveryNodes[idx]);
                }
                alternatives.push_back(altPath);
            }
        }

        // Strategy 2: Reverse tour
        if (alternatives.size() < numAlternatives + 1 && optimalTour.size() > 3) {
            vector<int> reversedTour = {optimalTour[0]};
            for (int i = optimalTour.size() - 2; i > 0; i--) {
                reversedTour.push_back(optimalTour[i]);
            }
            reversedTour.push_back(optimalTour.back());

            double distance = calculateTourDistance(reversedTour, distMatrix);
            if (distance > optimalPath.distance * 1.1) {
                PathOption altPath;
                altPath.warehouseId = warehouseId;
                altPath.distance = distance;
                for (int idx : reversedTour) {
                    altPath.path.push_back(deliveryNodes[idx]);
                }
                alternatives.push_back(altPath);
            }
        }

        return alternatives;
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
    if (n <= 1) return {0};
    
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
            for (int j = 0; j < n; j++) {
                if (!visited[j]) {
                    closest = j;
                    break;
                }
            }
        }
        
        if (closest == -1) break;
        tour.push_back(closest);
        visited[closest] = true;
    }
    
    if (tour.size() > 1 && distMatrix[tour.back()][0] != numeric_limits<double>::infinity()) {
        tour.push_back(0);
    }
    
    return tour;
}

void twoOptOptimize(vector<int>& tour, const vector<vector<double>>& distMatrix) {
    if (tour.size() <= 3) return;
    
    bool improved = true;
    int iterCount = 0;
    const int MAX_ITER = 1000;
    
    while (improved && iterCount < MAX_ITER) {
        improved = false;
        iterCount++;
        
        for (int i = 1; i < tour.size()-2; i++) {
            for (int j = i+1; j < tour.size()-1; j++) {
                double current1 = distMatrix[tour[i-1]][tour[i]];
                double current2 = distMatrix[tour[j]][tour[j+1]];
                double proposed1 = distMatrix[tour[i-1]][tour[j]];
                double proposed2 = distMatrix[tour[i]][tour[j+1]];
                
                if (proposed1 == numeric_limits<double>::infinity() ||
                    proposed2 == numeric_limits<double>::infinity()) {
                    continue;
                }
                
                double delta = (proposed1 + proposed2) - (current1 + current2);
                if (delta < -0.001) {
                    reverse(tour.begin()+i, tour.begin()+j+1);
                    improved = true;
                }
            }
        }
    }
}

double calculateTourDistance(const vector<int>& tour, const vector<vector<double>>& distMatrix) {
    double total = 0;
    for (int i = 1; i < tour.size(); i++) {
        if (distMatrix[tour[i-1]][tour[i]] != numeric_limits<double>::infinity()) {
            total += distMatrix[tour[i-1]][tour[i]];
        }
    }
    return total;
}

vector<int> solveTSP(const vector<vector<double>>& distMatrix) {
    if (distMatrix.empty()) return {};
    
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

void saveRouteToFile(const string& filename, const vector<string>& route) {
    ofstream outFile(filename);
    if (!outFile.is_open()) {
        cerr << "Error: Could not open " << filename << " for writing!" << endl;
        return;
    }

    for (const string& nodeId : route) {
        outFile << nodeId << endl;
    }
    outFile.close();
    cout << "Route saved to " << filename << endl;
}

int main() {
    RoadNetwork graph;
    
    cout << "Loading nodes..." << endl;
    loadNodes(graph, "../../data/nodes.csv");
    
    cout << "Loading edges..." << endl;
    loadEdges(graph, "../../data/edges.csv");
    
    graph.printGraphStats();
    
    vector<string> warehouseIds = {
        "12110005374",  // North Delhi
        "923637303",    // East Delhi
        "4314376797"    // South Delhi
    };
    
    vector<string> validWarehouses;
    for (const auto& id : warehouseIds) {
        if (graph.hasNode(id)) {
            validWarehouses.push_back(id);
            cout << "Warehouse " << id << " is valid." << endl;
        }
    }
    
    if (validWarehouses.empty()) {
        cerr << "Error: No valid warehouses!" << endl;
        return 1;
    }
    
    int numDeliveries;
    cout << "Enter number of deliveries: ";
    cin >> numDeliveries;
    
    auto deliveries = graph.selectRandomDeliveryNodes(numDeliveries, validWarehouses);
    auto assignments = graph.assignDeliveriesToWarehousesOptimized(validWarehouses, deliveries);
    
    for (const auto& warehouseId : validWarehouses) {
        if (assignments[warehouseId].empty()) continue;
        
        vector<string> routeNodes = {warehouseId};
        routeNodes.insert(routeNodes.end(), 
                         assignments[warehouseId].begin(), 
                         assignments[warehouseId].end());
        
        cout << "\nCalculating distance matrix for warehouse " << warehouseId << "..." << endl;
        auto distanceMatrix = graph.createDistanceMatrix(routeNodes);
        
        cout << "\nSolving TSP for warehouse " << warehouseId << "..." << endl;
        auto optimalTour = solveTSP(distanceMatrix);
        
        cout << "\nOptimal Delivery Route for Warehouse " << warehouseId << ":\n";
        for (size_t i = 0; i < optimalTour.size(); i++) {
            cout << i << ": " << routeNodes[optimalTour[i]] << endl;
        }
        
        double optimalDistance = calculateTourDistance(optimalTour, distanceMatrix);
        cout << "\nTotal route distance: " << optimalDistance << " meters" << endl;

        // Save optimal route
        string optimalFile = "warehouse_" + warehouseId + "_optimal_route.txt";
        vector<string> optimalRoute;
        for (int idx : optimalTour) {
            optimalRoute.push_back(routeNodes[idx]);
        }
        saveRouteToFile(optimalFile, optimalRoute);
        
        // Generate and save alternative routes
        cout << "\nGenerating alternative routes..." << endl;
        auto alternatives = graph.generateAlternativeRoutes(warehouseId, routeNodes, optimalTour, distanceMatrix, 3);
        
        for (size_t i = 1; i < alternatives.size(); i++) {
            string altFile = "warehouse_" + warehouseId + "_alternative_route_" + to_string(i) + ".txt";
            cout << "Alternative " << i << " distance: " << alternatives[i].distance << " meters (" 
                 << ((alternatives[i].distance - optimalDistance)/optimalDistance*100) << "% longer)" << endl;
            saveRouteToFile(altFile, alternatives[i].path);
        }
    }
    
    cout << "\nRunning visualization script..." << endl;
    string command = "python ../python/plot.py";
    int result = system(command.c_str());
    if (result != 0) {
        cerr << "Error: Failed to run visualization script" << endl;
    }
    
    return 0;
}