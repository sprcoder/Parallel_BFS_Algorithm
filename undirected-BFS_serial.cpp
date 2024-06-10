#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <unordered_map>
#include <limits.h>
#include <algorithm>
#include <stdint.h>
#include <queue>
#include <chrono>

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cout << "Input file not provided. Please provide an input file during execution" << endl;
        cout << "Use: " << argv[0] << " <input_file>" << endl;
        return 1;
    }

    char *filename = argv[1];
    ifstream file(filename);
    if (!file.is_open())
    {
        cout << "Unable to read " << filename << " provided as input" << endl;
        return 1;
    }

    // n - number of vertices
    // m - number of edges
    long n, m;
    file >> n >> m;

    // Create a hash table to store the adjacency list
    std::unordered_map<int, std::vector<int>> graph;

    // offsets info for each vertex
    std::vector<int> offsets(n);
    for (long i = 0; i < n; i++)
    {
        file >> offsets[i];
    }
    std::vector<int>  pointToVertex(m);
    // Read the edges and map the values accordingly
    for (long i = 0; i < m; i++)
    {
        file >> pointToVertex[i];
    }

    for (long i = 0; i < n; i++){
        for (long j = offsets[i]; j < (i + 1 < n ? offsets[i + 1] : m); j++){
            graph[i].push_back(pointToVertex[j]);
            // For undirected graph, you need to add edge in both directions
            graph[pointToVertex[j]].push_back(i);
        }
    }
    // Print the graph
    // cout << "Graph:" << endl;
    // for (const auto &entry : graph) {
    //     cout << "Vertex " << entry.first << ": ";
    //     for (int neighbor : entry.second) {
    //         cout << neighbor << " ";
    //     }
    //     cout << endl;
    // }
    auto starttime = high_resolution_clock::now();

    // Calculate eccentricity - Naive BFS
    std::vector<int> eccentricity(n, 0);

    // Get eccentricity values for every vertex
    for (long i = 0; i < n; i++) {
        std::vector<int> distances(n, 0);
        std::queue<int> q;

        // Get distance of vertices for every vertex i
        distances[i] = 0;
        q.push(i);

        while (!q.empty()) {
            int currVertex = q.front();
            q.pop();

            for (int neighbor : graph[currVertex]) {
                if (distances[neighbor] == 0) {
                    distances[neighbor] = distances[currVertex] + 1;
                    q.push(neighbor);
                }
            }
        }

        int maxDistance = 0;
        for (int d : distances) {
            maxDistance = max(maxDistance, d);
        }
        eccentricity[i] = maxDistance;
    }

    // // Weighted graph eccentricity calculation for each vertex
    // std::vector<int> eccentricity(n, std::numeric_limits<int>::min());
    // for (long i = 0; i < n; i++) {
    //     // Do Dijkstra's algorithm from vertex i to find distances to all other vertices
    //     std::vector<int> distance(n, std::numeric_limits<int>::max());
    //     std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pq;
    //     distance[i] = 0;
    //     pq.push({0, i});
    //     while (!pq.empty()) {
    //         int u = pq.top().second;
    //         pq.pop();
    //         for (const auto& edge : graph[u]) {
    //             int v = edge.destination;
    //             int weight = edge.weight;
    //             if (distance[v] > distance[u] + weight) {
    //                 distance[v] = distance[u] + weight;
    //                 pq.push({distance[v], v});
    //             }
    //         }
    //     }
    //     int maxDistance = *std::max_element(distance.begin(), distance.end());
    //     eccentricity[i] = maxDistance;
    // }
    auto endtime = high_resolution_clock::now();
    auto executiontime = duration_cast<duration<double>>(endtime - starttime);
 
    cout << "Time taken by function: "
         << executiontime.count() << " seconds" << endl;
}