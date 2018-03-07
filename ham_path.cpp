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

void ham_path(vector<int> &EC, vector<int> &HP, vector< vector<int> > &D, int numOfVertices, int* &hamLength)
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
		edIndex++;
		//check to see if it has been visited, if not, then add it to the Hamiltonian Path
		if (visited[EC[edIndex]] == false)
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
	for (int i = 0; i < HP.size-1; i++)
	{
		hamLength = hamLength + D[HP[i]][[i+1]];
	}
}

/*********************************************************************
** Code to write the data to file
*********************************************************************/

	std::ofstream outputFile;
	std::string newFileName;
	newFileName = inputFileName + ".tour";
	outputFile.open(newFileName.c_str());
	
	for (int i = 0; i < numSets; i++)
	{
		outputFile << hamLength;
		for (int i = 0; i < HP.size(); i++)
		{
			outputFile << HP[i] << "\n";
		}
	}
	outputFile.close();