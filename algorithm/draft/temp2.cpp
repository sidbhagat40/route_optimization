#include <bits/stdc++.h>
#include <chrono>
#include <thread>
#include <queue>
#include <limits>

// Define a structure to represent a geo-coordinate (latitude, longitude)
struct GeoCoordinate {
    double latitude;
    double longitude;

    // Default constructor
    GeoCoordinate() : latitude(0.0), longitude(0.0) {}

    // Parameterized constructor
    GeoCoordinate(double lat, double lon) : latitude(lat), longitude(lon) {}

    // Overload equality operator for unordered_map
    bool operator==(const GeoCoordinate& other) const {
        return latitude == other.latitude && longitude == other.longitude;
    }

    // Add inequality operator
    bool operator!=(const GeoCoordinate& other) const {
        return !(*this == other);
    }

    // Overload less than operator for map
    bool operator<(const GeoCoordinate& other) const {
        if (latitude == other.latitude) return longitude < other.longitude;
        return latitude < other.latitude;
    }
};

// Hash function for GeoCoordinate (required for unordered_map)
struct GeoCoordinateHash {
    std::size_t operator()(const GeoCoordinate& coord) const {
        return std::hash<double>()(coord.latitude) ^ std::hash<double>()(coord.longitude);
    }
};

// Function to calculate the Haversine distance between two geo-coordinates (in kilometers)
double haversineDistance(const GeoCoordinate& a, const GeoCoordinate& b) {
    const double R = 6371.0; // Earth radius in kilometers
    double dLat = (b.latitude - a.latitude) * M_PI / 180.0;
    double dLon = (b.longitude - a.longitude) * M_PI / 180.0;
    double lat1 = a.latitude * M_PI / 180.0;
    double lat2 = b.latitude * M_PI / 180.0;

    double x = sin(dLat / 2) * sin(dLat / 2) +
               sin(dLon / 2) * sin(dLon / 2) * cos(lat1) * cos(lat2);
    double y = 2 * atan2(sqrt(x), sqrt(1 - x));
    return R * y;
}

// Define a structure to represent a delivery point with an optional time window
struct Delivery {
    GeoCoordinate location;
    std::pair<int, int> time_window; // Time window in hours (e.g., 15-17 for 3 PM - 5 PM)
    bool has_time_window; // Flag to indicate if the delivery has a time window
    int priority; // Priority level (1-5, 1 being highest)
    int delivery_id; // Unique identifier for the delivery
    int estimated_service_time; // Estimated time to complete delivery (in minutes)

    Delivery(GeoCoordinate loc, std::pair<int, int> window, bool has_window, int prio = 3, int id = -1, int service_time = 15)
        : location(loc), time_window(window), has_time_window(has_window), 
          priority(prio), delivery_id(id), estimated_service_time(service_time) {}
};

// Define a structure to represent a vehicle
struct Vehicle {
    int id;
    GeoCoordinate current_location;
    std::vector<GeoCoordinate> route;
    double total_distance;
    std::vector<std::pair<int, int>> schedule;
    int current_time; // Current time in hours (e.g., 14 for 2 PM)
    int capacity; // Vehicle capacity (e.g., number of packages it can carry)
    int current_load; // Current number of packages loaded
    std::vector<int> delivery_ids; // IDs of deliveries assigned to this vehicle

    Vehicle(int id, GeoCoordinate depot, int cap = 20) 
        : id(id), current_location(depot), total_distance(0), current_time(9), // Start at 9 AM
          capacity(cap), current_load(0) {}

    void addToRoute(const GeoCoordinate& next_location, double distance, int travel_time, int delivery_id = -1) {
        route.push_back(next_location);
        total_distance += distance;
        current_location = next_location;
        current_time += travel_time; // Update current time based on travel time
        if (delivery_id != -1) {
            delivery_ids.push_back(delivery_id);
            current_load++;
        }
    }

    void addToSchedule(std::pair<int, int> time_window) {
        schedule.push_back(time_window);
    }

    bool canAcceptDelivery() const {
        return current_load < capacity;
    }
};

// Function to calculate travel time based on distance and speed
int calculateTravelTime(double distance, double speed = 50.0) {
    // Speed is in km/h, distance is in km
    return static_cast<int>((distance / speed) * 60); // Convert hours to minutes
}

// Graph representation for route planning
class RoadNetwork {
private:
    std::map<GeoCoordinate, std::vector<std::pair<GeoCoordinate, double>>> adj_list;

public:
    void addEdge(const GeoCoordinate& u, const GeoCoordinate& v, double distance) {
        adj_list[u].emplace_back(v, distance);
        adj_list[v].emplace_back(u, distance); // Assuming bidirectional roads
    }

    // Dijkstra's algorithm for shortest path
    std::pair<std::vector<GeoCoordinate>, double> dijkstra(const GeoCoordinate& start, const GeoCoordinate& end) const {
        std::priority_queue<std::pair<double, GeoCoordinate>, 
                          std::vector<std::pair<double, GeoCoordinate>>, 
                          std::greater<std::pair<double, GeoCoordinate>>> pq;
        
        std::map<GeoCoordinate, double> dist;
        std::map<GeoCoordinate, GeoCoordinate> prev;
        
        for (const auto& node : adj_list) {
            dist[node.first] = std::numeric_limits<double>::infinity();
        }
        
        dist[start] = 0.0;
        pq.push({0.0, start});
        
        while (!pq.empty()) {
            auto current = pq.top();
            pq.pop();
            
            double current_dist = current.first;
            GeoCoordinate u = current.second;
            
            if (u == end) break;
            if (current_dist > dist[u]) continue;
            
            for (const auto& neighbor : adj_list.at(u)) {
                GeoCoordinate v = neighbor.first;
                double weight = neighbor.second;
                
                if (dist[v] > dist[u] + weight) {
                    dist[v] = dist[u] + weight;
                    prev[v] = u;
                    pq.push({dist[v], v});
                }
            }
        }
        
        // Reconstruct path
        std::vector<GeoCoordinate> path;
        if (dist[end] == std::numeric_limits<double>::infinity()) {
            return {path, -1.0}; // No path found
        }
        
        for (GeoCoordinate at = end; at != start; at = prev[at]) {
            path.push_back(at);
        }
        path.push_back(start);
        std::reverse(path.begin(), path.end());
        
        return {path, dist[end]};
    }

    // A* algorithm with Haversine distance as heuristic
    std::pair<std::vector<GeoCoordinate>, double> aStar(const GeoCoordinate& start, const GeoCoordinate& end) const {
        auto heuristic = [](const GeoCoordinate& a, const GeoCoordinate& b) {
            return haversineDistance(a, b);
        };
        
        std::priority_queue<std::pair<double, GeoCoordinate>, 
                          std::vector<std::pair<double, GeoCoordinate>>, 
                          std::greater<std::pair<double, GeoCoordinate>>> open_set;
        
        std::map<GeoCoordinate, double> g_score;
        std::map<GeoCoordinate, double> f_score;
        std::map<GeoCoordinate, GeoCoordinate> came_from;
        
        for (const auto& node : adj_list) {
            g_score[node.first] = std::numeric_limits<double>::infinity();
            f_score[node.first] = std::numeric_limits<double>::infinity();
        }
        
        g_score[start] = 0.0;
        f_score[start] = heuristic(start, end);
        open_set.push({f_score[start], start});
        
        while (!open_set.empty()) {
            auto current = open_set.top();
            open_set.pop();
            
            GeoCoordinate u = current.second;
            
            if (u == end) {
                // Reconstruct path
                std::vector<GeoCoordinate> path;
                for (GeoCoordinate at = end; at != start; at = came_from[at]) {
                    path.push_back(at);
                }
                path.push_back(start);
                std::reverse(path.begin(), path.end());
                return {path, g_score[end]};
            }
            
            for (const auto& neighbor : adj_list.at(u)) {
                GeoCoordinate v = neighbor.first;
                double weight = neighbor.second;
                
                double tentative_g_score = g_score[u] + weight;
                
                if (tentative_g_score < g_score[v]) {
                    came_from[v] = u;
                    g_score[v] = tentative_g_score;
                    f_score[v] = g_score[v] + heuristic(v, end);
                    open_set.push({f_score[v], v});
                }
            }
        }
        
        return {std::vector<GeoCoordinate>(), -1.0}; // No path found
    }

    // Floyd-Warshall algorithm for all-pairs shortest paths
    std::map<std::pair<GeoCoordinate, GeoCoordinate>, double> floydWarshall() const {
        std::map<std::pair<GeoCoordinate, GeoCoordinate>, double> dist;
        
        // Initialize distances
        for (const auto& u : adj_list) {
            for (const auto& v : adj_list) {
                if (u.first == v.first) {
                    dist[{u.first, v.first}] = 0.0;
                } else {
                    dist[{u.first, v.first}] = std::numeric_limits<double>::infinity();
                }
            }
        }
        
        // Set initial edge weights
        for (const auto& u : adj_list) {
            for (const auto& neighbor : u.second) {
                dist[{u.first, neighbor.first}] = neighbor.second;
            }
        }
        
        // Floyd-Warshall main algorithm
        for (const auto& k : adj_list) {
            for (const auto& i : adj_list) {
                for (const auto& j : adj_list) {
                    if (dist[{i.first, j.first}] > dist[{i.first, k.first}] + dist[{k.first, j.first}]) {
                        dist[{i.first, j.first}] = dist[{i.first, k.first}] + dist[{k.first, j.first}];
                    }
                }
            }
        }
        
        return dist;
    }
};

// Function to optimize delivery routes with multiple constraints
std::vector<Vehicle> optimizeDeliveryRoutes(const GeoCoordinate& depot,
                                          std::vector<Delivery>& deliveries,
                                          int num_vehicles,
                                          const RoadNetwork& road_network) {
    std::vector<Vehicle> vehicles;
    for (int i = 0; i < num_vehicles; ++i) {
        vehicles.push_back(Vehicle(i + 1, depot));
    }

    // Sort deliveries by priority (highest first) and then by time window
    std::sort(deliveries.begin(), deliveries.end(), [](const Delivery& a, const Delivery& b) {
        if (a.priority != b.priority) return a.priority < b.priority;
        if (a.has_time_window != b.has_time_window) return a.has_time_window;
        if (a.has_time_window) return a.time_window.first < b.time_window.first;
        return false;
    });

    // Assign deliveries to vehicles
    for (auto& delivery : deliveries) {
        bool assigned = false;
        
        // Try to assign to existing vehicles first
        for (auto& vehicle : vehicles) {
            if (!vehicle.canAcceptDelivery()) continue;
            
            // Get the best path using A* algorithm
            auto result = road_network.aStar(vehicle.current_location, delivery.location);
            const auto& path = result.first;
            double distance = result.second;
            if (path.empty()) continue; // No path found
            
            int travel_time = calculateTravelTime(distance);
            int arrival_time = vehicle.current_time + travel_time;
            
            // Check constraints
            bool meets_constraints = true;
            if (delivery.has_time_window) {
                if (arrival_time > delivery.time_window.second) {
                    meets_constraints = false;
                }
            }
            
            if (meets_constraints) {
                vehicle.addToRoute(delivery.location, distance, travel_time, delivery.delivery_id);
                if (delivery.has_time_window) {
                    vehicle.addToSchedule(delivery.time_window);
                }
                assigned = true;
                break;
            }
        }
        
        // If not assigned, try to create a new route (if we have vehicle capacity)
        if (!assigned) {
            for (auto& vehicle : vehicles) {
                if (vehicle.route.empty()) { // Unused vehicle
                    auto result = road_network.aStar(depot, delivery.location);
                    const auto& path = result.first;
                    double distance = result.second;
                    if (path.empty()) continue;
                    
                    int travel_time = calculateTravelTime(distance);
                    int arrival_time = vehicle.current_time + travel_time;
                    
                    bool meets_constraints = true;
                    if (delivery.has_time_window) {
                        if (arrival_time > delivery.time_window.second) {
                            meets_constraints = false;
                        }
                    }
                    
                    if (meets_constraints) {
                        vehicle.addToRoute(delivery.location, distance, travel_time, delivery.delivery_id);
                        if (delivery.has_time_window) {
                            vehicle.addToSchedule(delivery.time_window);
                        }
                        assigned = true;
                        break;
                    }
                }
            }
        }
        
        // If still not assigned, we have a problem (not enough vehicles or impossible constraints)
        if (!assigned) {
            std::cerr << "Warning: Delivery " << delivery.delivery_id << " could not be assigned to any vehicle.\n";
        }
    }

    return vehicles;
}

// Function to dynamically re-route vehicles based on new information
void dynamicReRouting(Vehicle& vehicle, const std::vector<Delivery>& new_deliveries, 
                     const GeoCoordinate& depot, const RoadNetwork& road_network) {
    // Create a list of remaining deliveries (including new ones)
    std::vector<Delivery> remaining_deliveries;
    
    // Add new deliveries
    for (const auto& delivery : new_deliveries) {
        remaining_deliveries.push_back(delivery);
    }
    
    // Add any unvisited deliveries from the current route
    // (This would require tracking which deliveries have been completed)
    
    // Optimize the remaining route
    // For simplicity, we'll just find the nearest neighbor
    GeoCoordinate current_location = vehicle.current_location;
    std::vector<GeoCoordinate> new_route;
    double new_total_distance = 0.0;
    int new_current_time = vehicle.current_time;
    
    while (!remaining_deliveries.empty()) {
        double min_distance = std::numeric_limits<double>::infinity();
        auto nearest_it = remaining_deliveries.begin();
        
        for (auto it = remaining_deliveries.begin(); it != remaining_deliveries.end(); ++it) {
            auto result = road_network.aStar(current_location, it->location);
            const auto& path = result.first;
            double distance = result.second;
            
            if (path.empty()) continue;
            
            if (distance < min_distance) {
                min_distance = distance;
                nearest_it = it;
            }
        }
        
        if (min_distance == std::numeric_limits<double>::infinity()) {
            break; // No reachable deliveries left
        }
        
        auto result = road_network.aStar(current_location, nearest_it->location);
        const auto& path = result.first;
        double distance = result.second;
        int travel_time = calculateTravelTime(distance);
        
        new_route.insert(new_route.end(), path.begin() + 1, path.end());
        new_total_distance += distance;
        new_current_time += travel_time;
        current_location = nearest_it->location;
        
        remaining_deliveries.erase(nearest_it);
    }
    
    // Return to depot if needed
    if (!new_route.empty()) {
        auto result = road_network.aStar(current_location, depot);
        const auto& path = result.first;
        double distance = result.second;
        if (!path.empty()) {
            new_route.insert(new_route.end(), path.begin() + 1, path.end());
            new_total_distance += distance;
        }
    }
    
    // Update vehicle's route
    vehicle.route = new_route;
    vehicle.total_distance = new_total_distance;
    // Note: Would also need to update schedule and other parameters
}

int main() {
    // Example usage
    
    // 1. Create a road network
    RoadNetwork road_network;
    
    // Define some locations (in a real app, these would be actual coordinates)
    GeoCoordinate depot(26.466471045524983, 73.11222367206733);
    GeoCoordinate loc1(26.47523876708196, 73.11710999950559);
    GeoCoordinate loc2(26.475575331906825, 73.12028960504121);
    GeoCoordinate loc3(26.478363194734854, 73.11141385027695);
    GeoCoordinate loc4(26.409415827510088, 73.05786881559627);
    GeoCoordinate loc5(26.474657296139867, 73.11134106607929);
    
    // Add edges to the road network (with distances in km)
    road_network.addEdge(depot, loc1, 1.5);
    road_network.addEdge(depot, loc3, 1.2);
    road_network.addEdge(loc1, loc2, 0.8);
    road_network.addEdge(loc1, loc5, 0.7);
    road_network.addEdge(loc2, loc3, 1.1);
    road_network.addEdge(loc3, loc4, 8.5);
    road_network.addEdge(loc3, loc5, 0.5);
    road_network.addEdge(loc5, depot, 1.0);
    
    // 2. Create deliveries with various constraints
    std::vector<Delivery> deliveries = {
        {loc1, {15, 17}, true, 1, 1}, // High priority, time window 3-5 PM
        {loc2, {0, 0}, false, 3, 2},   // Medium priority, no time window
        {loc3, {10, 12}, true, 2, 3}, // Medium-high priority, time window 10 AM-12 PM
        {loc4, {0, 0}, false, 5, 4},   // Low priority, no time window
        {loc5, {0, 0}, false, 4, 5}    // Low-medium priority, no time window
    };
    
    // 3. Optimize delivery routes
    int num_vehicles = 2;
    std::vector<Vehicle> vehicles = optimizeDeliveryRoutes(depot, deliveries, num_vehicles, road_network);
    
    // 4. Print the optimized routes and schedules for each vehicle
    for (const auto& vehicle : vehicles) {
        std::cout << "Vehicle " << vehicle.id << " Route:\n";
        for (const auto& coord : vehicle.route) {
            std::cout << "(" << coord.latitude << ", " << coord.longitude << ")\n";
        }
        std::cout << "Total Distance: " << vehicle.total_distance << " km\n";
        std::cout << "Schedule:\n";
        for (const auto& window : vehicle.schedule) {
            std::cout << window.first << ":00 - " << window.second << ":00\n";
        }
        std::cout << "Current Time: " << vehicle.current_time << ":00\n";
        std::cout << "Assigned Deliveries: ";
        for (auto id : vehicle.delivery_ids) {
            std::cout << id << " ";
        }
        std::cout << "\n\n";
    }
    
    // 5. Demonstrate dynamic re-routing
    if (!vehicles.empty()) {
        std::cout << "=== Dynamic Re-Routing Example ===\n";
        
        // Simulate a new delivery coming in
        GeoCoordinate new_loc(26.470123, 73.115678);
        Delivery new_delivery(new_loc, {16, 18}, true, 1, 6); // High priority, time window 4-6 PM
        
        // Add this location to the road network
        road_network.addEdge(loc1, new_loc, 0.6);
        road_network.addEdge(new_loc, loc5, 0.9);
        
        // Perform dynamic re-routing for the first vehicle
        std::vector<Delivery> new_deliveries = {new_delivery};
        dynamicReRouting(vehicles[0], new_deliveries, depot, road_network);
        
        // Print the updated route
        std::cout << "Vehicle 1 Updated Route:\n";
        for (const auto& coord : vehicles[0].route) {
            std::cout << "(" << coord.latitude << ", " << coord.longitude << ")\n";
        }
        std::cout << "Updated Total Distance: " << vehicles[0].total_distance << " km\n";
    }
    
    // 6. Demonstrate all-pairs shortest paths with Floyd-Warshall
    std::cout << "\n=== Floyd-Warshall All-Pairs Shortest Paths ===\n";
    auto all_pairs_dist = road_network.floydWarshall();
    for (const auto& entry : all_pairs_dist) {
    const auto& pair = entry.first;
    double dist = entry.second;

    std::cout << "From (" << pair.first.latitude << ", " << pair.first.longitude << ") to "
              << "(" << pair.second.latitude << ", " << pair.second.longitude << "): "
              << dist << " km\n";
}


    
    return 0;
}