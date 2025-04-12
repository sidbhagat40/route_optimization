#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <limits>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sstream>

using namespace std;

const double INF = numeric_limits<double>::max();
const double EARTH_RADIUS = 6371.0; // in kilometers

// ===================== LOCATION STRUCT ======================
struct Location {
    int id;
    string name;
    double latitude;
    double longitude;

    // Default constructor
    Location() : id(-1), name(""), latitude(0.0), longitude(0.0) {}

    // Parameterized constructor
    Location(int id, const string& name, double lat, double lon)
        : id(id), name(name), latitude(lat), longitude(lon) {}
};

// ===================== HAVERSINE DISTANCE ======================
double haversineDistance(const Location& a, const Location& b) {
    double lat1 = a.latitude * M_PI / 180.0;
    double lon1 = a.longitude * M_PI / 180.0;
    double lat2 = b.latitude * M_PI / 180.0;
    double lon2 = b.longitude * M_PI / 180.0;

    double dlat = lat2 - lat1;
    double dlon = lon2 - lon1;

    double h = sin(dlat / 2) * sin(dlat / 2) +
               cos(lat1) * cos(lat2) * sin(dlon / 2) * sin(dlon / 2);

    double c = 2 * atan2(sqrt(h), sqrt(1 - h));
    return EARTH_RADIUS * c;
}

// ===================== TIME FORMAT ======================
string formatTime(int minutes) {
    int h = minutes / 60;
    int m = minutes % 60;
    ostringstream oss;
    oss << setw(2) << setfill('0') << h << ":"
        << setw(2) << setfill('0') << m;
    return oss.str();
}

// ===================== GRAPH CLASS ======================
class Graph {
public:
    unordered_map<int, Location> locations;
    unordered_map<int, vector<pair<int, double>>> adjList;

    void addLocation(Location loc) {
        locations[loc.id] = loc;
    }

    void addEdge(int from, int to) {
        double dist = haversineDistance(locations[from], locations[to]);
        adjList[from].emplace_back(to, dist);
        adjList[to].emplace_back(from, dist); // Undirected graph
    }

    void dijkstra(int src, unordered_map<int, double>& dist) {
        priority_queue<pair<double, int>, vector<pair<double, int>>, greater<>> pq;
        dist.clear();

        for (const auto& pair : locations) {
            dist[pair.first] = INF;
        }

        dist[src] = 0;
        pq.emplace(0.0, src);

        while (!pq.empty()) {
            pair<double, int> current = pq.top(); pq.pop();
            double d = current.first;
            int u = current.second;

            if (d > dist[u]) continue;

            for (const auto& neighbor : adjList[u]) {
                int v = neighbor.first;
                double weight = neighbor.second;
                if (dist[u] + weight < dist[v]) {
                    dist[v] = dist[u] + weight;
                    pq.emplace(dist[v], v);
                }
            }
        }
    }
};

// ===================== DELIVERY STRUCT ======================
struct Delivery {
    int locationID;
    int priority;        // Lower = higher priority
    int earliestTime;    // In minutes
    int latestTime;      // In minutes
};

// ===================== VEHICLE STRUCT ======================
struct Vehicle {
    int id;
    int capacity;
    vector<pair<int, int>> assignedLocations; // locationID and ETA
};

// ===================== DISPATCHER CLASS ======================
// ===================== DISPATCHER CLASS ======================
class Dispatcher {
    Graph& graph;
    int depotID;

public:
    Dispatcher(Graph& g, int depotID) : graph(g), depotID(depotID) {}

    void assignDeliveries(vector<Vehicle>& vehicles, vector<Delivery>& deliveries) {
        // Sort deliveries by priority (lower priority number means higher priority)
        sort(deliveries.begin(), deliveries.end(), [](const Delivery& a, const Delivery& b) {
            return a.priority < b.priority;
        });

        // Assign deliveries to vehicles based on available capacity
        int v = 0;
        for (Delivery& d : deliveries) {
            vehicles[v % vehicles.size()].assignedLocations.push_back({d.locationID, -1});
            v++;  // Move to the next vehicle
        }
    }

    void planRoutes(vector<Vehicle>& vehicles, const vector<Delivery>& deliveries) {
        unordered_map<int, double> dist;
        unordered_map<int, Delivery> deliveryMap;
        for (const auto& d : deliveries) deliveryMap[d.locationID] = d;

        const int deliveryTimePerStop = 15;  // Time per delivery in minutes
        const int startTime = 540;           // 9:00 AM in minutes
        const int endTime = 1020;            // 5:00 PM in minutes

        // Plan routes for each vehicle
        for (Vehicle& v : vehicles) {
            vector<pair<int, int>> stops = v.assignedLocations;
            vector<pair<int, int>> route;
            int curr = depotID;
            int currentTime = startTime;

            while (!stops.empty()) {
                graph.dijkstra(curr, dist);

                int nextIdx = -1;
                double minDist = INF;

                // Find the next best delivery that can be completed within time constraints
                for (int i = 0; i < stops.size(); ++i) {
                    int candidate = stops[i].first;
                    double arrivalTime = currentTime + dist[candidate];

                    const Delivery& d = deliveryMap[candidate];
                    if (arrivalTime <= d.latestTime && arrivalTime >= d.earliestTime &&
                        arrivalTime + deliveryTimePerStop <= endTime) {
                        if (dist[candidate] < minDist) {
                            minDist = dist[candidate];
                            nextIdx = i;
                        }
                    }
                }

                if (nextIdx == -1) {
                    cout << " Vehicle " << v.id << " cannot service all deliveries due to time constraints.\n";
                    break;
                }

                int nextLoc = stops[nextIdx].first;
                int eta = static_cast<int>(currentTime + dist[nextLoc]);
                currentTime = eta + deliveryTimePerStop;
                route.push_back({nextLoc, eta});
                curr = nextLoc;
                stops.erase(stops.begin() + nextIdx); // Remove the assigned stop
            }

            v.assignedLocations = route;
        }
    }

    void displayRoutes(const vector<Vehicle>& vehicles, const unordered_map<int, Location>& locs) {
        for (const Vehicle& v : vehicles) {
            cout << "\n Vehicle " << v.id << " Route:\n";
            cout << "  Start at: " << locs.at(depotID).name << " (09:00)\n";
            for (const auto& assignment : v.assignedLocations) {
                int locID = assignment.first;
                int eta = assignment.second;
                cout << "  -> Deliver to: " << locs.at(locID).name
                     << " (ETA: " << formatTime(eta) << ")\n";
            }
            cout << "  -> Return to: " << locs.at(depotID).name << "\n";
        }
    }
};

// ===================== MAIN FUNCTION ======================
int main() {
    Graph city;

    // Add Locations (ID, Name, Latitude, Longitude)
    city.addLocation({0, "Depot", 26.466284587413128, 73.11177376078871});         // SF
    city.addLocation({1, "Location A", 26.23837109674883, 73.00619526985768});    // Oakland
    city.addLocation({2, "Location B", 26.219393744183503, 73.04135354904878});    // Daly City
    city.addLocation({3, "Location C",26.47520225789069, 73.11708515743788});    // Palo Alto
    city.addLocation({4, "Location D", 26.476457578760822, 73.11990278504588});    // San Jose
    city.addLocation({5, "Location E", 26.474167116435424, 73.11186391712427});    // Tracy (fake)

    // Define Edges (connections between locations)
    city.addEdge(0, 1);
    city.addEdge(0, 2);
    city.addEdge(1, 3);
    city.addEdge(2, 3);
    city.addEdge(3, 4);
    city.addEdge(4, 5);
    city.addEdge(2, 5);

    // Deliveries: {locationID, priority, earliestTime, latestTime}
    vector<Delivery> deliveries = {
        {1, 1, 540, 1440},   // 09:00 - 10:30
        {2, 1, 540, 1440},   // 10:00 - 12:30
        {3, 1, 540, 1440},   // 09:30 - 11:30
        {4, 1, 540, 1440},   // 11:30 - 15:00
        {5, 1, 540, 1440}    // 10:30 - 14:00
    };

    // Define Vehicles: {id, capacity}
    vector<Vehicle> vehicles = {
        {1, 3}, {2, 3}  // Vehicle 1 and Vehicle 2
    };

    Dispatcher dispatcher(city, 0);
    dispatcher.assignDeliveries(vehicles, deliveries);
    dispatcher.planRoutes(vehicles, deliveries);
    dispatcher.displayRoutes(vehicles, city.locations);

    return 0;
}

