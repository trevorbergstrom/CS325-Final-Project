#include <stdio.h>
#include <iostream>
#include <vector>
#include <stack>

int getDegree(std::vector<int> adj) {
    int sum = 0;
    for(int i = 0; i < adj.size(); i++) {
        sum += adj[i];
    }
    return sum;
}

/* THIS FUNCTION DOES NOT WORK WITH INPUT FILES */
void eulerTour(std::vector<std::vector<int>> &m, std::vector<int> &tour) {
    std::vector<int> numNeighbors; //vector that contains the number of neighbors of each node
    bool eulerTourPossible = true;
    
    for(int i = 0; i < m.size(); i++) { //fill num neighbors while checking the degree of each node
        int degree = getDegree(m[i]);
        numNeighbors.push_back(degree);
        if((degree % 2) != 0) {
            eulerTourPossible = false;
        }
    }
    
    std::stack<int> eulerStack; //stack to hold nodes passed through
    int i = 0;
    int j = 0;
    
    while(!eulerStack.empty() || numNeighbors[i] != 0) { //loop while the stack is not empty and the current node has neighbors
        
        if(numNeighbors[i] == 0) { //if the current node has no neighbors
            tour.push_back(i); //add this node to the tour
            i = eulerStack.top(); //pop a node off the stack as the current node
            eulerStack.pop();
        } else { //if the current node has neighbors
            bool found = false;
            
            while(found != true) { //loop through its adjacency vector
                if(j >= m[i].size()) {
                    j = 0;
                }
                
                if(m[i][j] == 2) { //if the node has a loop with itself
                    eulerStack.push(i); //push this node on the stack
                    m[i][j] = 0; //remove edge in adj matrix
                    m[j][i] = 0;
                    numNeighbors[i] -= 2; //reduce the number of nieghbors
                    found = true;
                } else if(m[i][j] == 1) { //if the node has a edge between another node
                    if(numNeighbors[j] > 1 || numNeighbors[i] == 1) {
                        eulerStack.push(i); //push current node on the stach
                        numNeighbors[i] -= 1;
                        numNeighbors[j] -= 1; //reduce the number of neighbors
                        m[i][j] = 0;//remove edge in adj matrix
                        m[j][i] = 0;
                        i = j; //set the current node as the neighbor selected
                        found = true;
                    } else {
                        j++;
                    }
                } else {
                    j++;
                }
            }
        }
    }
    tour.push_back(0);
}

int findNumEdges(std::vector<std::vector<int>> &m) {
    int size = (int) m.size();
    int sum = 0;
    
    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            sum += m[i][j];
        }
    }
    return sum;
}

/* THIS FUNCTION WORKS WITH INPUT FILES and main from christofides.cpp */
void eulerTour2(std::vector<std::vector<int>> &m, std::vector<int> &tour) {
    int numNodes = (int) m[0].size();
    int i = 0;
    int j = 0;
    bool done = false;
    int numEdgesTotal = (findNumEdges(m) / 2);
    if((numEdgesTotal % 2) != 0) {
        numEdgesTotal -= 1;
    }
    int numEdgesVisited = 0;
    
    while(done == false) {
        
        while(m[i][j] == 0) {
            if(numEdgesVisited >= numEdgesTotal) {
                done = true;
                break;
            }
            if(j == numNodes) {
                j = 0;
            } else {
                j++;
            }
        }
        
        m[i][j] = 0;
        m[j][i] = 0;
        numEdgesVisited++;
        tour.push_back(i);
        i = j;
    }
    tour.push_back(0);
}

/*DRIVER TO TEST*/
int main(int argc, const char * argv[]) {
    
    std::vector<std::vector<int>> m;
    
    /*Adjacency matrix for graph with a euler tour*/
    std::vector<int> a = {0,1,1,0,0,0};
    std::vector<int> b = {1,0,1,1,1,0};
    std::vector<int> c = {1,1,0,1,1,0};
    std::vector<int> d = {0,1,1,0,1,1};
    std::vector<int> e = {0,1,1,1,0,1};
    std::vector<int> f = {0,0,0,1,1,0};
    
    m = {a,b,c,d,e,f};
    
    std::vector<int> tour;
    
    eulerTour(m, tour);
    for(int i = 0; i < tour.size(); i++) {
        std::cout << tour[i] << ", ";
    }
    
    
    return 0;
}

