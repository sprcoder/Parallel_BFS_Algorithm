#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <set>
#include <unordered_map>
#include <limits.h>
#include <algorithm>
#include <stdint.h>
#include <omp.h>
#include <random>

using namespace std;

const int inf = INT_MAX; // Infinity value

// Constants
const int w = 64;  // Word size
int k = 10; // Sample size
long skipped = 0;
long orOperations = 0;
long swaps = 0;

// Global variables
unordered_map<int, vector<pair<int, int>>> graph;  // Adjacency list representation of the graph
vector<uint64_t> Visited;     // Visited array
vector<uint64_t> NextVisited; // NextVisited array
vector<int> Ecc;                        // Eccentricity array
vector<int> AccumulatedWeight;          // Accumulated weight array
int round_num = 0;

// Atomic OR operation
void AtomicOR(uint64_t *location, uint64_t value)
{
#pragma omp atomic
   *location |= value;
}

// Compare and swap
bool CAS(int *location, int old_value, int new_value)
{
   bool success = false;
   {
      
      if (*location == old_value)
      {
         // #pragma omp atomic
         swaps++;
         *location = new_value;
         success = true;
      }
   }
   return success;
}

// Update atomic function
bool updateAtomic(int s, int d, int length, vector<uint64_t>& VisitedArray, vector<uint64_t>& NextVisitedArray, vector<int>& ecc, int weight)
{
   bool changed = false; 
   for(int i = 0; i < length; i++) {
      uint64_t toWrite = VisitedArray[d * length + (i / w)] | VisitedArray[s * length + (i / w)];
      orOperations++;
      if (VisitedArray[d * length + (i / w)] != toWrite)
      {
         AtomicOR(&NextVisitedArray[d * length + (i / w)], toWrite);
         int oldEcc = ecc[d];
         if(ecc[d] < weight) {
            // cout << "Ecc is different from the one to be updated, entering CAS\n";
            #pragma omp atomic
            changed |= CAS(&Ecc[d], oldEcc, AccumulatedWeight[d] + weight);
         }
      } else {
         skipped++;
      }
   }
   return changed;
}

// Copy function
void COPY(int i)
{
   Visited[i] = NextVisited[i];
}

// Initialization function
void INIT(int i, int j)
{
   unsigned long long bit = 1ULL << (j % w);
   Visited[i] |= bit;
   // #pragma omp critical
   // cout << "Given j: "<< j << " value: "<< Visited[i]<<endl;
   NextVisited[i] |= bit;
   Ecc[i] = 0;
   AccumulatedWeight[i] = 0;
}

void printFrontier(vector<int> Frontier){
   for (int i = 0; i < Frontier.size(); i++){
      cout << Frontier[i] << "\t";
   }
   cout << "\n";
}

// Compute eccentricity
void COMPUTE_ECC(vector<int> &Frontier)
{
   // Initialize frontier
   #pragma omp parallel for
   for (int i = 0; i < Frontier.size(); i++)
   {
      INIT(Frontier[i], i);
   }

   // Calculate length
   long length = (k + w - 1) / w;

   while (!Frontier.empty())
   {
      round_num++;
      // cout << "Round " << round_num << "\n";
      vector<int> NextFrontier;
      // omp_set_num_threads(Frontier.size());
      // Process each vertex in the current frontier
      #pragma omp parallel for
      for (int i = 0; i < Frontier.size(); i++)
      {
         int v = Frontier[i];
         // cout << "Current Node " << v << "\n";
         #pragma omp parallel for
         for (int j = 0; j < graph[v].size(); j++)
         {
            int neighbor = graph[v][j].first; // Neighbor vertex
            int weight = graph[v][j].second; // Edge weight

            if (updateAtomic(v, neighbor, length, Visited, NextVisited, Ecc, AccumulatedWeight[v] + weight))
            {
               // cout << "Neighbor is getting added\n";
               #pragma omp critical
               {
                  NextFrontier.push_back(neighbor);
                  AccumulatedWeight[neighbor] = AccumulatedWeight[v] + weight;
               }
            }
         }
      }

      Frontier = NextFrontier;
      sort(Frontier.begin(), Frontier.end());

      // Remove duplicates
      Frontier.erase(unique(Frontier.begin(), Frontier.end()), Frontier.end());
      // random_device rd;
      // mt19937 g(rd());
      // shuffle(Frontier.begin(), Frontier.end(), g);

      // Select first 5 elements
      // vector<int> randomValues(Frontier.begin(), Frontier.begin() + min<size_t>(k, Frontier.size()));
      // Frontier = randomValues;
      // cout << "Next frontier size: " << Frontier.size() << "\n";
      
      #pragma omp parallel for
      for (int i = 0; i < Frontier.size(); i++)
      {
         // cout << Frontier[i] << " Updating \n";
         COPY(Frontier[i]);
         // cout << "Nextnodes " << Frontier[i] << "\n";
      }
   }
}

// Function to find k vertices with maximum eccentricity
vector<int> findMaxEccentricityVertices(vector<int>& Ecc, int k) {
   vector<int> result;
   vector<pair<int, int>> indices;
   for (int i = 0; i < Ecc.size(); ++i) {
      indices.push_back({Ecc[i], i});
   }
   sort(indices.begin(), indices.end(), greater<pair<int, int>>());
   for (int i = 0; i < k; ++i) {
      // cout << indices[i].first << "\t";
      result.push_back(indices[i].second);
   }
   // cout << "\n";
   return result;
}

int main(int argc, char *argv[])
{
   if (argc != 3)
   {
      cout << "Input file/No.of Process not provided. Please provide an input file during execution" << endl;
      cout << "Use: " << argv[0] << "<No.of Process> <input_file>" << endl;
      return 1;
   }

   k = atoi(argv[1]);
   char *filename = argv[2];
   ifstream file(filename);
   if (!file.is_open())
   {
      cout << "Unable to read " << filename << " provided as input" << endl;
      return 1;
   }
   // omp_set_num_threads(4);
   // n - number of vertices
   // m - number of edges
   long n, m;
   file >> n >> m;

   // offsets info for each vertex
   vector<int> offsets(n);
   for (long i = 0; i < n; i++)
   {
      file >> offsets[i];
   }
   vector<int>  pointToVertex(m);
   // Read the edges and map the values accordingly
   for (long i = 0; i < m; i++)
   {
      file >> pointToVertex[i];
   }

   for (long i = 0; i < n; i++){
      for (long j = offsets[i]; j < (i + 1 < n ? offsets[i + 1] : m); j++){
         int weight;
         file >> weight;
         graph[i].push_back(make_pair(pointToVertex[j], weight));
         // For undirected graph, you need to add edge in both directions
         graph[pointToVertex[j]].push_back(make_pair(i, weight));
      }
   }

   // Print the adjacency list
   for (long i = 0; i < n; i++)
   {
      // cout << "Vertex " << i << " :";
      for (auto edge : graph[i])
      {
         // cout << " " << edge.first << "(" << edge.second << ")";
      }
      // cout << endl;
   }
   // std::cout << "Adjacency list created... \n Find random k frontier vertices" << std::endl;

   omp_set_num_threads(k);

   // Initialize Visited, NextVisited, Ecc, and AccumulatedWeight arrays
   // Calculate length
   long length = (k + w - 1) / w;
   Visited.resize(n * length, 0);
   NextVisited.resize(n * length, 0);
   Ecc.resize(graph.size(), 0);
   AccumulatedWeight.resize(graph.size(), 0);

   // Compute eccentricity
   double start_time = omp_get_wtime();
   // Randomly sample k vertices
   vector<int> Frontier;
   // cout << "Calculating random vertex values\n";
   srand(100);
   for (int i = 0; i < k; i++)
   {
      Frontier.push_back(rand() % graph.size());
   }
   COMPUTE_ECC(Frontier);
   double end_time = omp_get_wtime();
   cout << "Phase 1 Time taken: "<< end_time - start_time <<" seconds\n";
   cout << "skipped nodes as they are already visited: " << skipped << "\n";
   // cout << "Vertex\tEccentricity\tGraph Size " << graph.size();
   vector<int> Ecc1 = Ecc;

   Visited.resize(n * length, 0);
   NextVisited.resize(n * length, 0);
   Ecc.resize(graph.size(), 0);
   cout << Ecc1[0] << "Safe to create Ecc1";
   AccumulatedWeight.resize(graph.size(), 0);
   // Reinitialize Frontier with k vertices with maximum eccentricity
   double start_time_phase2 = omp_get_wtime();
   Frontier = findMaxEccentricityVertices(Ecc, k);
   // printFrontier(Frontier);
   // Finalize eccentricity calculation
   
   COMPUTE_ECC(Frontier);
   double end_time_phase2 = omp_get_wtime();
   cout << "Phase 2 Time taken: "<< end_time_phase2 - start_time_phase2 <<" seconds\n";
   cout << "Total Time taken: "<< end_time_phase2 - start_time <<" seconds\n";
   cout << "OrOperations: " << orOperations << endl;
   cout << "swaps: " << swaps << endl;
   cout << "skipped: " << skipped << endl;

   // Output eccentricity
   ofstream outputFile("output.txt");
   if (!outputFile.is_open()) {
        std::cerr << "Error opening the file." << std::endl;
        return 1; // Return error code
    }
   outputFile << "Vertex\tEccentricity\tGraph Size " << graph.size() << "\n";
   for (int i = 0; i < graph.size(); i++)
   {
      outputFile << i << "\t" << max(Ecc[i], Ecc1[i]) << endl;
   }
}
