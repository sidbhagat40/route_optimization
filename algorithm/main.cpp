#include <bits/stdc++.h>
#include <chrono>
#include <thread>
#include <queue>
#include <limits>

struct GeoCoordinate {
    double latitude;
    double longitude;

    GeoCoordinate() : latitude(0.0), longitude(0.0) {}

    GeoCoordinate(double lat, double lon) : latitude(lat), longitude(lon) {}

    bool operator==(const GeoCoordinate& other) const {
        return latitude == other.latitude && longitude == other.longitude;
    }

    bool operator!=(const GeoCoordinate& other) const {
        return !(*this == other);
    }

    bool operator<(const GeoCoordinate& other) const {
        if (latitude == other.latitude) return longitude < other.longitude;
        return latitude < other.latitude;
    }
};

struct GeoCoordinateHash {
    std::size_t operator()(const GeoCoordinate& coord) const {
        return std::hash<double>()(coord.latitude) ^ std::hash<double>()(coord.longitude);
    }
};

double haversineDistance(const GeoCoordinate& a, const GeoCoordinate& b) {
    const double R = 6371.0;
    double dLat = (b.latitude - a.latitude) * M_PI / 180.0;
    double dLon = (b.longitude - a.longitude) * M_PI / 180.0;
    double lat1 = a.latitude * M_PI / 180.0;
    double lat2 = b.latitude * M_PI / 180.0;

    double x = sin(dLat / 2) * sin(dLat / 2) +
               sin(dLon / 2) * sin(dLon / 2) * cos(lat1) * cos(lat2);
    double y = 2 * atan2(sqrt(x), sqrt(1 - x));
    return R * y;
}

struct Delivery {
    GeoCoordinate location;
    std::pair<int, int> time_window;
    bool has_time_window;
    int priority;
    int delivery_id;
    int estimated_service_time;

    Delivery(GeoCoordinate loc, std::pair<int, int> window, bool has_window, int prio = 3, int id = -1, int service_time = 15)
        : location(loc), time_window(window), has_time_window(has_window), 
          priority(prio), delivery_id(id), estimated_service_time(service_time) {}
};

struct Vehicle {
    int id;
    GeoCoordinate current_location;
    std::vector<GeoCoordinate> route;
    double total_distance;
    std::vector<std::pair<int, int>> schedule;
    int current_time;
    int capacity;
    int current_load;
    std::vector<int> delivery_ids;

    Vehicle(int id, GeoCoordinate depot, int cap = 20) 
        : id(id), current_location(depot), total_distance(0), current_time(9),
          capacity(cap), current_load(0) {}

    void addToRoute(const GeoCoordinate& next_location, double distance, int travel_time, int delivery_id = -1) {
        route.push_back(next_location);
        total_distance += distance;
        current_location = next_location;
        current_time += travel_time;
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

int calculateTravelTime(double distance, double speed = 50.0) {
    return static_cast<int>((distance / speed) * 60);
}

class RoadNetwork {
private:
    std::map<GeoCoordinate, std::vector<std::pair<GeoCoordinate, double>>> adj_list;

public:
    void addEdge(const GeoCoordinate& u, const GeoCoordinate& v, double distance) {
        adj_list[u].emplace_back(v, distance);
        adj_list[v].emplace_back(u, distance);
    }

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
        
        std::vector<GeoCoordinate> path;
        if (dist[end] == std::numeric_limits<double>::infinity()) {
            return {path, -1.0};
        }
        
        for (GeoCoordinate at = end; at != start; at = prev[at]) {
            path.push_back(at);
        }
        path.push_back(start);
        std::reverse(path.begin(), path.end());
        
        return {path, dist[end]};
    }

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
        
        return {std::vector<GeoCoordinate>(), -1.0};
    }

    std::map<std::pair<GeoCoordinate, GeoCoordinate>, double> floydWarshall() const {
        std::map<std::pair<GeoCoordinate, GeoCoordinate>, double> dist;
        
        for (const auto& u : adj_list) {
            for (const auto& v : adj_list) {
                if (u.first == v.first) {
                    dist[{u.first, v.first}] = 0.0;
                } else {
                    dist[{u.first, v.first}] = std::numeric_limits<double>::infinity();
                }
            }
        }
        
        for (const auto& u : adj_list) {
            for (const auto& neighbor : u.second) {
                dist[{u.first, neighbor.first}] = neighbor.second;
            }
        }
        
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

std::vector<Vehicle> optimizeDeliveryRoutes(const GeoCoordinate& depot,
                                          std::vector<Delivery>& deliveries,
                                          int num_vehicles,
                                          const RoadNetwork& road_network) {
    std::vector<Vehicle> vehicles;
    for (int i = 0; i < num_vehicles; ++i) {
        vehicles.push_back(Vehicle(i + 1, depot));
    }

    std::sort(deliveries.begin(), deliveries.end(), [](const Delivery& a, const Delivery& b) {
        if (a.priority != b.priority) return a.priority < b.priority;
        if (a.has_time_window != b.has_time_window) return a.has_time_window;
        if (a.has_time_window) return a.time_window.first < b.time_window.first;
        return false;
    });

    for (auto& delivery : deliveries) {
        bool assigned = false;
        
        for (auto& vehicle : vehicles) {
            if (!vehicle.canAcceptDelivery()) continue;
            
            auto result = road_network.aStar(vehicle.current_location, delivery.location);
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
        
        if (!assigned) {
            for (auto& vehicle : vehicles) {
                if (vehicle.route.empty()) {
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
        
        if (!assigned) {
            std::cerr << "Warning: Delivery " << delivery.delivery_id << " could not be assigned to any vehicle.\n";
        }
    }

    return vehicles;
}

void dynamicReRouting(Vehicle& vehicle, const std::vector<Delivery>& new_deliveries, 
                     const GeoCoordinate& depot, const RoadNetwork& road_network) {
    std::vector<Delivery> remaining_deliveries;
    
    for (const auto& delivery : new_deliveries) {
        remaining_deliveries.push_back(delivery);
    }
    
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
            break;
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
    
    if (!new_route.empty()) {
        auto result = road_network.aStar(current_location, depot);
        const auto& path = result.first;
        double distance = result.second;
        if (!path.empty()) {
            new_route.insert(new_route.end(), path.begin() + 1, path.end());
            new_total_distance += distance;
        }
    }
    
    vehicle.route = new_route;
    vehicle.total_distance = new_total_distance;
}

int main() {
    RoadNetwork road_network;
    
    GeoCoordinate depot(26.466471045524983, 73.11222367206733);
    GeoCoordinate loc1(26.47523876708196, 73.11710999950559);
    GeoCoordinate loc2(26.475575331906825, 73.12028960504121);
    GeoCoordinate loc3(26.478363194734854, 73.11141385027695);
    GeoCoordinate loc4(26.409415827510088, 73.05786881559627);
    GeoCoordinate loc5(26.474657296139867, 73.11134106607929);
    
    road_network.addEdge(depot, loc1, 1.5);
    road_network.addEdge(depot, loc3, 1.2);
    road_network.addEdge(loc1, loc2, 0.8);
    road_network.addEdge(loc1, loc5, 0.7);
    road_network.addEdge(loc2, loc3, 1.1);
    road_network.addEdge(loc3, loc4, 8.5);
    road_network.addEdge(loc3, loc5, 0.5);
    road_network.addEdge(loc5, depot, 1.0);
    
    std::vector<Delivery> deliveries = {
        {loc1, {15, 17}, true, 1, 1},
        {loc2, {0, 0}, false, 3, 2},
        {loc3, {10, 12}, true, 2, 3},
        {loc4, {0, 0}, false, 5, 4},
        {loc5, {0, 0}, false, 4, 5}
    };
    
    int num_vehicles = 2;
    std::vector<Vehicle> vehicles = optimizeDeliveryRoutes(depot, deliveries, num_vehicles, road_network);
    
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
    
    if (!vehicles.empty()) {
        std::cout << "=== Dynamic Re-Routing Example ===\n";
        
        GeoCoordinate new_loc(26.470123, 73.115678);
        Delivery new_delivery(new_loc, {16, 18}, true, 1, 6);
        
        road_network.addEdge(loc1, new_loc, 0.6);
        road_network.addEdge(new_loc, loc5, 0.9);
        
        std::vector<Delivery> new_deliveries = {new_delivery};
        dynamicReRouting(vehicles[0], new_deliveries, depot, road_network);
        
        std::cout << "Vehicle 1 Updated Route:\n";
        for (const auto& coord : vehicles[0].route) {
            std::cout << "(" << coord.latitude << ", " << coord.longitude << ")\n";
        }
        std::cout << "Updated Total Distance: " << vehicles[0].total_distance << " km\n";
    }
    
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