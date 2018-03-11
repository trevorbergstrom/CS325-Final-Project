/******************************************
 * ** Author: Tim Ip, Rich Reese, Trevor Bergstrom
 * ** Date: 3/4/2018 
 * ** Description CS 325 TSP Project 
 * This program will use christofides 
 * method to solve the TSP problem 
 * ***************************************/

#include <iostream> 
#include <limits> 
#include <fstream> 
#include <string> 
#include <vector> 
#include <deque> 
#include <sstream> 
#include <cmath> 
#include <ctime>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <list> 
#include <stack>

//Function that takes the Adjacency matrix of edge weights D as arg and
//a min spanning tree as arg. 
//Modifies the MST by pairing all odd degree vertices s.t. 
//there are no vertices with odd degree after
void perfectMatching(const std::vector<std::vector<int> >& D, std::vector<std::vector<int> >& M){
	//Create list of odd degree vertices from MST
	std::list<int> odds; 
	std::cout << "Odds list contains " << std::endl; 
	for(int ii=0; ii < D.size(); ii++){
		int edgeSum = 0; 
		for(int jj=0; jj < D.size(); jj++){
			edgeSum += M[ii][jj]; 	
		}
		if(edgeSum % 2 != 0){
			odds.push_front(ii); 
			std::cout << ii << ", " ;
		}
	}
	std::cout << "\n"; 
	//Match and remove all odd vertices from list	
	while( !odds.empty() ){
		int v1 = odds.front();
		odds.pop_front(); 
		int edgeLength = std::numeric_limits<int>::max(); 
		int v2;
		for(std::list<int>::iterator ii=odds.begin(); ii != odds.end(); ii++){
			if(D[v1][*ii] < edgeLength){
				edgeLength = D[v1][*ii]; 
				v2 = *ii; 
			}
		}
		M[v1][v2] += 1; 
		M[v2][v1] += 1; 
		odds.remove(v2); 
	}
}


//Function that calculates a minimum spanning tree
// Takes a 2d adjacency matrix D and writes to the 2d matrix M passed in
void MST(const std::vector<std::vector<int> >& D, std::vector<std::vector<int> >& M){
    std::vector<int> C(D.size(), std::numeric_limits<int>::max()); //Cheapest edge cost to vertex
    std::vector<int> E(D.size(), -1); //-1 no edge to indexed vertex, otherwise value is city that edge is to
    std::vector<int> Q(D.size(), 0); //Vertices not yet included in graph = 0 
    int totalVertices = D.size(); 
    Q[0] = 1; //Add random vertex to graph, I chose vertex 0
    totalVertices --; 
    while (totalVertices > 0){
        for(int ii=0; ii < D.size(); ii++){
            if(Q[ii] == 1){
                for(int jj=0; jj < D.size(); jj++){
                    if(D[ii][jj] < C[ii] && Q[jj]==0){
                        C[jj] = D[ii][jj]; 
                        E[jj] = ii; 
                    }
                }
            }
        }
        int minCost = std::numeric_limits<int>::max(); 
        int minCostVert = -1; 
        int minCostEdgeTo = -1; 
        for(int ii=0; ii < D.size(); ii++){
            if(C[ii] < minCost){
                minCost = C[ii]; 
                minCostVert = ii;
                minCostEdgeTo = E[ii]; 
            }
        }
        Q[minCostVert] = 1; 
        C[minCostVert] = std::numeric_limits<int>::max(); 
        M[minCostVert][minCostEdgeTo] = 1;
        M[minCostEdgeTo] [minCostVert] = 1;
        totalVertices --; 
    }
}

int getDegree(std::vector<int> adj) {
    int sum = 0;
    for(int i = 0; i < adj.size(); i++) {
        sum += adj[i];
    }
    return sum;
}

/* THIS FUNCTION DOES NOT WORK WITH INPUT FILES */
void eulerTour(std::vector<std::vector<int> > &m, std::vector<int> &tour) {
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

int findNumEdges(std::vector<std::vector<int> > &m) {
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
void eulerTour2(std::vector<std::vector<int> > &m, std::vector<int> &tour) {
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

/*********************************************************************
** Function: ham_path
** Description: Takes in a Euler circuit and finds a Hamiltonian path
** for the circuit
** Input: 	pass by reference the vector holding the Euler circuit
**			pass by reference the empty vector to hold the ham path
**			which is of length n+1, where n is the number of vertices
**			pass by reference an int to hold the length of the ham path
**			pass by reference the matrix of distances between vertices
**			pass by number the number of vertices in the graph
** Output:	the ham path vector will be populated
**			the tspLength int will be populated
*********************************************************************/

void ham_path(std::vector<int> &EC, std::vector<int> &HP, std::vector< std::vector<int> > &D, int numOfVertices, int* &hamLength)
{
	/*****make the Hamiltonian path*****/
	
	//make a vector to hold a true or false depending on if the vertex has been visited
	std::vector<bool> visited (numOfVertices,false);
	
	//set the current vertex to -1 since it will be updated to 0 in the while loop
	int ecIndex = -1;
	
	//make a counter to countdown until there all vertices have been visited
	int verticesLeft = numOfVertices;
	
	//continue until all the vertices have been visited
	while (verticesLeft != 0)
	{
		//look at the next vertex in the Euler circuit
		ecIndex++;
		//check to see if it has been visited, if not, then add it to the Hamiltonian Path
		if (visited[EC[ecIndex]] == false)
		{
			visited[EC[ecIndex]] = true;
			HP.push_back (EC[ecIndex]);
			verticesLeft--;
		}
	}
	
	//add the first vertex in the Euler circuit to close the Hamiltonian path
	HP.push_back (EC[0]);
	
	/*****calculate the length of the path*****/
	
	//initialize the length of the Hamiltonian path to 0
	hamLength = 0;
	
	//run through each edge in the Hamiltonian path and add the value to the overall Hamiltonian path length
	for (int i = 0; i < HP.size()-1; i++)
	{
	hamLength = hamLength + D[HP[i]][i+1];
	}
}

int main(int argc, char** argv)
{
    clock_t begin = std::clock(); 
    std::ifstream inFile(argv[1]);
    std::vector<int> city;  
    std::vector<int> xCoords;  
    std::vector<int> yCoords;  
    int xIn, yIn, cityIn; 

//Import data from file given as argument
    while(inFile >> cityIn >> xIn >> yIn ){
    	city.push_back(cityIn); 
	    xCoords.push_back(xIn); 
    	yCoords.push_back(yIn); 
    }
    int cityCount = city.size(); 

    //Compute matrix of Distances between cities
    //Initilized Eta, desirability of state transitions
    std::vector<std::vector<int> > Distances(
		cityCount, std::vector<int>(cityCount, 0 )); 
    double xDiff, yDiff;
    for(int ii=0; ii < cityCount; ii++){
        for(int jj=0; jj < cityCount; jj++){
            xDiff = pow(xCoords[ii] - xCoords[jj], 2); 
            yDiff = pow(yCoords[ii] - yCoords[jj], 2); 
            Distances[ii][jj] = round(sqrt(xDiff+yDiff)); 
        //std::cout << Distances[ii][jj] << ", "; 
        }
        //std::cout << '\n';  
    }
    //std::cout << '\n';  

    std::vector<std::vector<int> > M(
		cityCount, std::vector<int>(cityCount, 0 )); 

    //Call Fxn to get minimum spanning tree
    MST(Distances, M); 

    //Print MST to check
	std::cout << "Adjacency Matrix" << std::endl;
    for(int ii=0; ii < M.size(); ii++){
	for(int jj=0; jj < M.size(); jj++){
		std::cout << M[ii][jj] <<", ";
	}
    std::cout << "\n";  
    }
    std::cout << "\n";  

    perfectMatching(Distances, M); 
    //Print matched adjacency graph to check
	std::cout << "Adjacency graph" << std::endl;
    for(int ii=0; ii < M.size(); ii++){
	for(int jj=0; jj < M.size(); jj++){
		std::cout << M[ii][jj] <<", ";
	}
    std::cout << "\n";  
    }
	
	//Create the Euler tour
	std::vector<int> tour;
	std::cout << "starting Euler tour" << std::endl;
	eulerTour2(M, tour);
	
	//Print the Euler tour
	std::cout << "Euler tour" << std::endl;
    for(int i = 0; i < tour.size(); i++) {
        std::cout << tour[i] << ", ";
    }
	std::cout << std::endl;
	
	//Create the Hamiltonian path
	std::vector<int> hamPath;
	std::cout << "starting Hamiltonian path" << std::endl;
	int* hamLength = 0;
	ham_path(tour, hamPath, Distances, city.size(), hamLength);
	
	//print solution
	std::cout << "tsp length: " << hamLength << std::endl;
	for (int i = 0; i < hamPath.size(); i++)
	{
		std::cout << hamPath[i] << ", ";
	}
	
	delete hamLength;
}