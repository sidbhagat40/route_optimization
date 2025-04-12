#include <bits/stdc++.h>

using namespace std;

// Constants
const double INF = numeric_limits<double>::infinity();

// Data structures
struct Location {
    int id;
    double latitude;
    double longitude;
    string name;
    
    Location(int id, double lat, double lon, string name) 
        : id(id), latitude(lat), longitude(lon), name(name) {}
};

struct Road {
    int from;
    int to;
    double distance; // in km
    double base_time; // in minutes (without traffic)
    double current_time; // in minutes (with traffic)
    
    Road(int from, int to, double dist, double time)
        : from(from), to(to), distance(dist), base_time(time), current_time(time) {}
};

struct DeliveryOrder {
    int id;
    int pickup_location;
    int delivery_location;
    double ready_time; // when the order is ready for pickup
    double due_time; // when the order must be delivered by
    int priority; // higher number = higher priority
    
    DeliveryOrder(int id, int pickup, int delivery, double ready, double due, int prio)
        : id(id), pickup_location(pickup), delivery_location(delivery), 
          ready_time(ready), due_time(due), priority(prio) {}
};

struct Vehicle {
    int id;
    int current_location;
    double current_time; // current time in minutes since start
    double capacity; // remaining capacity
    vector<int> assigned_orders;
    
    Vehicle(int id, int loc, double cap)
        : id(id), current_location(loc), current_time(0), capacity(cap) {}
};

class TrafficSimulator {
private:
    mutex mtx;
    default_random_engine generator;
    uniform_real_distribution<double> dist;
    
public:
    TrafficSimulator() : dist(0.5, 2.0) {
        generator.seed(chrono::system_clock::now().time_since_epoch().count());
    }
    
    void updateTraffic(vector<Road>& roads) {
        lock_guard<mutex> lock(mtx);
        for (auto& road : roads) {
            double factor = dist(generator); // random traffic factor
            road.current_time = road.base_time * factor;
        }
    }
};

class DeliveryOptimizer {
private:
    vector<Location> locations;
    vector<Road> roads;
    vector<DeliveryOrder> orders;
    vector<Vehicle> vehicles;
    unordered_map<int, vector<pair<int, double>>> graph; // adjacency list: location -> [(neighbor, time)]
    
    mutex data_mutex;
    TrafficSimulator traffic_simulator;
    bool running = false;
    thread optimization_thread;
    thread traffic_thread;
    
    // For real-time updates
    condition_variable cv;
    bool new_data_available = false;
    
    void buildGraph() {
        graph.clear();
        for (const auto& road : roads) {
            graph[road.from].emplace_back(road.to, road.current_time);
            graph[road.to].emplace_back(road.from, road.current_time); // assuming bidirectional roads
        }
    }
    
    // Dijkstra's algorithm for shortest path
    pair<vector<double>, unordered_map<int, int>> dijkstra(int start) {
        priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;
        vector<double> dist(locations.size(), INF);
        unordered_map<int, int> prev;
        
        dist[start] = 0;
        pq.emplace(0, start);
        
        while (!pq.empty()) {
            auto [current_dist, u] = pq.top();
            pq.pop();
            
            if (current_dist > dist[u]) continue;
            
            for (const auto& [v, time] : graph[u]) {
                if (dist[v] > dist[u] + time) {
                    dist[v] = dist[u] + time;
                    prev[v] = u;
                    pq.emplace(dist[v], v);
                }
            }
        }
        
        return {dist, prev};
    }
    
    // A* algorithm with heuristic (haversine distance)
    pair<vector<double>, unordered_map<int, int>> astar(int start, int goal) {
        auto heuristic = [&](int a, int b) {
            const auto& loc_a = locations[a];
            const auto& loc_b = locations[b];
            
            // Haversine distance (simplified for this example)
            double lat1 = loc_a.latitude * M_PI / 180.0;
            double lon1 = loc_a.longitude * M_PI / 180.0;
            double lat2 = loc_b.latitude * M_PI / 180.0;
            double lon2 = loc_b.longitude * M_PI / 180.0;
            
            double dlat = lat2 - lat1;
            double dlon = lon2 - lon1;
            
            double a = sin(dlat / 2) * sin(dlat / 2) + 
                       cos(lat1) * cos(lat2) * 
                       sin(dlon / 2) * sin(dlon / 2);
            double c = 2 * atan2(sqrt(a), sqrt(1 - a));
            double distance = 6371 * c; // Earth radius in km
            
            // Convert to time estimate (assuming average speed of 50 km/h)
            return distance * 1.2; // 1.2 = 60/50 (convert km to minutes at 50 km/h)
        };
        
        priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;
        vector<double> g_score(locations.size(), INF);
        vector<double> f_score(locations.size(), INF);
        unordered_map<int, int> came_from;
        
        g_score[start] = 0;
        f_score[start] = heuristic(start, goal);
        pq.emplace(f_score[start], start);
        
        while (!pq.empty()) {
            auto [current_f, u] = pq.top();
            pq.pop();
            
            if (u == goal) {
                break;
            }
            
            if (current_f > f_score[u]) continue;
            
            for (const auto& [v, time] : graph[u]) {
                double tentative_g = g_score[u] + time;
                
                if (tentative_g < g_score[v]) {
                    came_from[v] = u;
                    g_score[v] = tentative_g;
                    f_score[v] = g_score[v] + heuristic(v, goal);
                    pq.emplace(f_score[v], v);
                }
            }
        }
        
        return {g_score, came_from};
    }
    
    vector<int> reconstructPath(int start, int goal, const unordered_map<int, int>& came_from) {
        vector<int> path;
        if (came_from.find(goal) == came_from.end()) {
            return path; // no path found
        }
        
        int current = goal;
        while (current != start) {
            path.push_back(current);
            current = came_from.at(current);
        }
        path.push_back(start);
        reverse(path.begin(), path.end());
        return path;
    }
    
    // Simple greedy assignment for demonstration
    void assignOrdersToVehicles() {
        lock_guard<mutex> lock(data_mutex);
        
        // Reset vehicle assignments (in a real system, we'd be more careful)
        for (auto& vehicle : vehicles) {
            vehicle.assigned_orders.clear();
        }
        
        // Sort orders by priority then due time
        vector<DeliveryOrder> sorted_orders = orders;
        sort(sorted_orders.begin(), sorted_orders.end(), 
            [](const DeliveryOrder& a, const DeliveryOrder& b) {
                if (a.priority != b.priority) return a.priority > b.priority;
                return a.due_time < b.due_time;
            });
        
        // Assign to vehicles (simple round-robin for this example)
        size_t vehicle_idx = 0;
        for (const auto& order : sorted_orders) {
            if (vehicle_idx >= vehicles.size()) vehicle_idx = 0;
            vehicles[vehicle_idx].assigned_orders.push_back(order.id);
            vehicle_idx++;
        }
    }
    
    // Optimize routes for all vehicles
    void optimizeRoutes() {
        lock_guard<mutex> lock(data_mutex);
        
        for (auto& vehicle : vehicles) {
            if (vehicle.assigned_orders.empty()) continue;
            
            // For simplicity, we'll just find a path that picks up and delivers all assigned orders
            // In a real system, we'd use a more sophisticated algorithm like VRP with time windows
            
            vector<int> stops;
            stops.push_back(vehicle.current_location);
            
            // Collect all pickup and delivery locations
            unordered_set<int> pickup_locations;
            unordered_set<int> delivery_locations;
            for (int order_id : vehicle.assigned_orders) {
                auto it = find_if(orders.begin(), orders.end(), 
                                 [order_id](const DeliveryOrder& o) { return o.id == order_id; });
                if (it != orders.end()) {
                    pickup_locations.insert(it->pickup_location);
                    delivery_locations.insert(it->delivery_location);
                }
            }
            
            // Add pickups first, then deliveries (simplified approach)
            for (int loc : pickup_locations) {
                stops.push_back(loc);
            }
            for (int loc : delivery_locations) {
                stops.push_back(loc);
            }
            
            // Find a path through all stops (TSP approximation)
            vector<int> route;
            int current = vehicle.current_location;
            unordered_set<int> visited;
            visited.insert(current);
            route.push_back(current);
            
            while (route.size() < stops.size()) {
                // Find nearest unvisited stop
                double min_dist = INF;
                int next_stop = -1;
                auto [distances, _] = dijkstra(current);
                
                for (int stop : stops) {
                    if (visited.count(stop) == 0 && distances[stop] < min_dist) {
                        min_dist = distances[stop];
                        next_stop = stop;
                    }
                }
                
                if (next_stop == -1) break; // no more reachable stops
                
                // Get the path to the next stop
                auto [_, came_from] = astar(current, next_stop);
                vector<int> path_segment = reconstructPath(current, next_stop, came_from);
                
                // Add the path segment to the route (excluding the current location)
                for (size_t i = 1; i < path_segment.size(); i++) {
                    route.push_back(path_segment[i]);
                }
                
                current = next_stop;
                visited.insert(current);
            }
            
            // Update vehicle's planned route
            cout << "Vehicle " << vehicle.id << " optimized route: ";
            for (int loc : route) {
                cout << locations[loc].name << " -> ";
            }
            cout << "END" << endl;
            
            // In a real system, we'd store the route and estimated times
        }
    }
    
    void optimizationLoop() {
        while (running) {
            unique_lock<mutex> lock(data_mutex);
            cv.wait(lock, [this] { return new_data_available || !running; });
            
            if (!running) break;
            
            // Process new data
            buildGraph();
            assignOrdersToVehicles();
            optimizeRoutes();
            
            new_data_available = false;
            lock.unlock();
            
            // Sleep for a bit to simulate periodic optimization
            this_thread::sleep_for(chrono::seconds(5));
        }
    }
    
    void trafficUpdateLoop() {
        while (running) {
            {
                lock_guard<mutex> lock(data_mutex);
                traffic_simulator.updateTraffic(roads);
                new_data_available = true;
            }
            cv.notify_one();
            
            this_thread::sleep_for(chrono::seconds(10));
        }
    }
    
public:
    DeliveryOptimizer() {
        // Initialize with some sample data
        locations.emplace_back(0, 40.7128, -74.0060, "New York");
        locations.emplace_back(1, 34.0522, -118.2437, "Los Angeles");
        locations.emplace_back(2, 41.8781, -87.6298, "Chicago");
        locations.emplace_back(3, 29.7604, -95.3698, "Houston");
        locations.emplace_back(4, 39.9526, -75.1652, "Philadelphia");
        
        // Add roads between locations (simplified)
        roads.emplace_back(0, 1, 3940, 4740); // NY to LA
        roads.emplace_back(0, 2, 1143, 1372); // NY to Chicago
        roads.emplace_back(0, 4, 137, 164);   // NY to Philly
        roads.emplace_back(1, 3, 2192, 2630); // LA to Houston
        roads.emplace_back(2, 3, 1500, 1800); // Chicago to Houston
        roads.emplace_back(2, 4, 1020, 1224); // Chicago to Philly
        roads.emplace_back(3, 4, 2090, 2508); // Houston to Philly
        
        // Add some vehicles
        vehicles.emplace_back(0, 0, 100); // Vehicle 0 in NY, capacity 100
        vehicles.emplace_back(1, 1, 100); // Vehicle 1 in LA, capacity 100
        
        buildGraph();
    }
    
    ~DeliveryOptimizer() {
        stop();
    }
    
    void start() {
        if (running) return;
        
        running = true;
        optimization_thread = thread(&DeliveryOptimizer::optimizationLoop, this);
        traffic_thread = thread(&DeliveryOptimizer::trafficUpdateLoop, this);
    }
    
    void stop() {
        if (!running) return;
        
        running = false;
        cv.notify_all();
        
        if (optimization_thread.joinable()) optimization_thread.join();
        if (traffic_thread.joinable()) traffic_thread.join();
    }
    
    void addOrder(const DeliveryOrder& order) {
        lock_guard<mutex> lock(data_mutex);
        orders.push_back(order);
        new_data_available = true;
        cv.notify_one();
    }
    
    void printStatus() {
        lock_guard<mutex> lock(data_mutex);
        
        cout << "\n=== System Status ===" << endl;
        cout << "Locations: " << locations.size() << endl;
        cout << "Roads: " << roads.size() << endl;
        cout << "Active Orders: " << orders.size() << endl;
        cout << "Vehicles: " << vehicles.size() << endl;
        
        for (const auto& vehicle : vehicles) {
            cout << "\nVehicle " << vehicle.id << " at " << locations[vehicle.current_location].name;
            cout << " with " << vehicle.assigned_orders.size() << " orders" << endl;
        }
    }
    
    // For demonstration: simulate some traffic updates
    void simulateTrafficEvent(int road_index, double delay_factor) {
        lock_guard<mutex> lock(data_mutex);
        if (road_index >= 0 && road_index < roads.size()) {
            roads[road_index].current_time = roads[road_index].base_time * delay_factor;
            new_data_available = true;
            cv.notify_one();
        }
    }
};

int main() {
    DeliveryOptimizer optimizer;
    
    // Add some sample orders
    optimizer.addOrder(DeliveryOrder(1, 0, 2, 0, 1440, 1)); // NY to Chicago, due in 24h, priority 1
    optimizer.addOrder(DeliveryOrder(2, 1, 3, 0, 720, 2));  // LA to Houston, due in 12h, priority 2
    optimizer.addOrder(DeliveryOrder(3, 4, 0, 0, 360, 3));  // Philly to NY, due in 6h, priority 3
    
    // Start the optimization system
    optimizer.start();
    
    // Let the system run for a while
    this_thread::sleep_for(chrono::seconds(2));
    optimizer.printStatus();
    
    // Simulate a traffic jam
    cout << "\nSimulating traffic jam on NY-Chicago road..." << endl;
    optimizer.simulateTrafficEvent(1, 3.0); // Road index 1 (NY-Chicago) now 3x slower
    
    this_thread::sleep_for(chrono::seconds(5));
    optimizer.printStatus();
    
    // Add another order
    cout << "\nAdding new high-priority order..." << endl;
    optimizer.addOrder(DeliveryOrder(4, 2, 4, 0, 180, 5)); // Chicago to Philly, due in 3h, priority 5
    
    this_thread::sleep_for(chrono::seconds(5));
    optimizer.printStatus();
    
    // Stop the system
    optimizer.stop();
    
    return 0;
}