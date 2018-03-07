/******************************************
 * ** Author: Tim Ip
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

//Function that takes the Adjacency matris of edge weights D as arg and
//a min spanning tree as arg. 
//Modifies the MST by pairing all odd degree vertices s.t. 
//there are no vertices with odd degree after
void perfectMatching(const std::vector<std::vector<int> >& D, std::vector<std::vector<int> >& M){
	//Create list of odd degree vertices from MST
	std::list<int> odds; 
	std::cout << "Odds list contains " ; 
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
    for(int ii=0; ii < M.size(); ii++){
	for(int jj=0; jj < M.size(); jj++){
		std::cout << M[ii][jj] <<", ";
	}
    std::cout << "\n";  
    }
    std::cout << "\n";  

    perfectMatching(Distances, M); 
    //Print matched adjacency graph to check
    for(int ii=0; ii < M.size(); ii++){
	for(int jj=0; jj < M.size(); jj++){
		std::cout << M[ii][jj] <<", ";
	}
    std::cout << "\n";  
    }
}







