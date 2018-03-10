/*********************************************************
 * Trevor Bergstrom
 * CS325 Final Project
 * Nearest Neighbor Algorithm for TSP
 * 3/9/18
 *********************************************************/

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdio.h>
#include <time.h>
#include <limits.h>

std::vector<std::vector<int>> distance_builder(std::string fileName) {
    std::ifstream inFile(fileName);
    std::vector<int> city;
    std::vector<int> xCoords;
    std::vector<int> yCoords;
    int xIn, yIn, cityIn;
    
    while(inFile >> cityIn >> xIn >> yIn ){
        city.push_back(cityIn);
        xCoords.push_back(xIn);
        yCoords.push_back(yIn);
    }
    int cityCount = (int) city.size();
    
    std::vector<std::vector<int>> Distances(cityCount,std::vector<int>(cityCount, 0));
    double xDiff, yDiff;
    for(int ii=0; ii < cityCount; ii++){
        for(int jj=0; jj < cityCount; jj++){
            xDiff = pow(xCoords[ii] - xCoords[jj], 2);
            yDiff = pow(yCoords[ii] - yCoords[jj], 2);
            Distances[ii][jj] = round(sqrt(xDiff+yDiff));
            
        }
    }
    
    return Distances;
}

void swapV(std::vector<int> &v, int x, int y) {
    int temp = v[x];
    v[x] = v[y];
    v[y] = temp;
}

void quicksort2Vectors(std::vector<int> &distances, std::vector<int> &neighbors, int low, int high) {
    int i, j, mid, p;
    
    i = low;
    j = high;
    mid = low + (high - low) / 2;
    p = distances[mid];
    
    while (i < high || j > low) {
        while(distances[i] < p) {
            i++;
        }
        while(distances[j] > p) {
            j--;
        }
        
        if(i <= j) {
            swapV(distances, i, j);
            swapV(neighbors, i, j);
            i++;
            j--;
        } else {
            if(i < high) {
                quicksort2Vectors(distances, neighbors, i, high);
            }
            if(j > low) {
                quicksort2Vectors(distances, neighbors, low, j);
            }
            return;
        }
    }
}

void quicksortVector(std::vector<int> &v, int low, int high) {
    int i, j, mid, p;
    
    i = low;
    j = high;
    mid = low + (high - low) / 2;
    p = v[mid];
    
    while (i < high || j > low) {
        while(v[i] < p) {
            i++;
        }
        while(v[j] > p) {
            j--;
        }
        
        if(i <= j) {
            swapV(v, i, j);
            i++;
            j--;
        } else {
            if(i < high) {
                quicksortVector(v, i, high);
            }
            if(j > low) {
                quicksortVector(v, low, j);
            }
            return;
        }
    }
}

std::vector<std::vector<int>> build_closest_cities(std::vector<std::vector<int>> distances) {
    
    int numCities = (int) distances.size();
    std::vector<std::vector<int>> closest_cities(numCities,std::vector<int>(numCities, 0));
    
    for(int i = 0; i < numCities; i++) {
        for(int j = 0; j < numCities; j++) {
            closest_cities[i][j] = j;
        }
    }
    
    for(int i = 0; i < numCities; i++) {
        quicksort2Vectors(distances[i], closest_cities[i], 0, (int)distances[i].size() - 1);
    }
    return closest_cities;
}

std::vector<int> greedy_path(std::vector<std::vector<int>> &neighbors) {
    std::vector<int> path;
    int numCities = (int) neighbors.size();
    std::vector<int> is_visited(numCities, 0);
    int i = 0;
    int numVisited = 0;
    
    while(path.size() < numCities - 1) {
        //cout << "for vertex "<< i << ": \n";
        numVisited++;
        path.push_back(i); // add the vertex to the greedy path
        is_visited[i] = 1; //check its visited
        
        int j = 0;
        int next = neighbors[i][j];
        while(is_visited[next] == 1) { //find the next unvisited city
            j++;
            //cout << j << " ";
            next = neighbors[i][j];
        }
        i = next;
        
    }
    path.push_back(i);
    path.push_back(0);
    return path;
}

int find_shortest_path(std::vector<int> distances, std::vector<int> is_visited) {
    int shortest = std::numeric_limits<int>::max();
    int position = -1;
    
    for(int i = 0; i < distances.size(); i++) {
        if(distances[i] < shortest && is_visited[i] != 1) {
            shortest = distances[i];
            position = i;
        }
    }
    return position;
}

std::vector<int> greedy_path2(std::vector<std::vector<int>> distances) {
    std::vector<int> path;
    int numCities = (int) distances.size();
    std::vector<int> is_visited(numCities, 0);
    int i = 0;
    
    while(path.size() < numCities - 1) {
        path.push_back(i);
        is_visited[i] = 1;
        
        int j = find_shortest_path(distances[i], is_visited);
        
        
        i = j;
    }
    path.push_back(i);
    path.push_back(0);
    return path;
}


int length_of_path(std::vector<int> path, std::vector<std::vector<int>> distances) {
    int total_distance = 0;
    for(int i = 0; i < path.size() - 1; i++) {
        int start = path[i];
        int end = path[i + 1];
        
        total_distance += distances[start][end];
    }
    return total_distance;
}

std::vector<int> tourSwap(std::vector<int> tour, int i, int k) {
    int size = (int) tour.size();
    std::vector<int> tourCopy;
    
    for(int j = 0; j <= i - 1; j++) {
        tourCopy.push_back(tour[j]);
    }
    
    int d = 0;
    
    for(int j = i; j <= k; j++) {
        tourCopy.push_back(tour[k - d]);
        d++;
    }
    for(int j = k + 1; j < size; j++) {
        tourCopy.push_back(tour[j]);
    }
    
    return tourCopy;
}

std::vector<int> tourOpt(std::vector<int> tour, std::vector<std::vector<int>> distances, int lengthOfPath) {
    int size = (int) tour.size();
    
    int improve = 0;
    int bestDistance = lengthOfPath;
    
    while(improve < 5) {
        
        for(int i = 0; i < size - 1; i++) {
            for(int k = i + 1; k < size; k++) {
                std::vector<int> newTour = tourSwap(tour, i, k);
                int newDist = length_of_path(newTour, distances);
                
                if(newDist < bestDistance) {
                    improve = 0;
                    tour = newTour;
                    bestDistance = newDist;
                }
            }
        }
        improve++;
    }
    return tour;
}

int main(int argc, const char * argv[]) {
    
    std::string file = "ex1.txt";
    clock_t start, stop;
    start = clock();
    std::vector<std::vector<int>> distances = distance_builder(file);
    std::vector<std::vector<int>> closest = build_closest_cities(distances);
    std::vector<int> myPath = greedy_path(closest);
    
    
    int length = length_of_path(myPath, distances);
    stop = clock();
    float diff ((float)stop-(float)start);
    std::cout << "Final path length for Greedy1: " << length << " found in : " << diff / CLOCKS_PER_SEC << " seconds\n";
    
    start = clock();
    
    std::vector<int> myPath2 = greedy_path2(distances);
    
    
    int length2 = length_of_path(myPath, distances);
    stop = clock();
    float diff2 ((float)stop-(float)start);
    std::cout << "Final path length for Greedy2: " << length2 << " found in : " << diff2 / CLOCKS_PER_SEC << " seconds\n";
    
     start = clock();
    std::vector<int> optimized_tour = tourOpt(myPath2, distances, length2);
    int length3 = length_of_path(optimized_tour, distances);
    stop = clock();
    float diff3 ((float)stop-(float)start);
    std::cout << "Final path length for Greedy2 with 2-opt optimization: " << length3 << " found in : " << diff3 / CLOCKS_PER_SEC << " seconds\n";

    return 0;
}
