//CSC 399 MPI Project Sequential KDD Tree
//William Fawcet, Justin Mizelle
//This program randomly generates a KDD tree, takes in a point from a txt file,
//then computes the nearest points and distance and outputs this to a txt file.

#include "stdafx.h"
#include "KDTree.h"
#include <iostream>
#include "mpi.h"
#include <fstream>
#include <stdlib.h>
#include <string>
#include <ctime>

using namespace std;

KDTree<int> Tree;

void generate_random_point(int* p, int sd)
{
	for(int k=0; k<sd; k++)
		p[k] = rand() % RMAX ;
}

void GenerateRandomTree(int nNodes)
{
	int p [SD]; // random point
	for(int n=0; n<nNodes; n++)
	{
		generate_random_point(p, SD);
		Tree.add(p); 
	}
}


int main(int argc, char* argv[])
{
	//start time measurment
	clock_t begin = clock();

	// Initialize MPI environment
	MPI_Init(&argc, &argv);

	// Get number of processes
	int numProcesses;
	MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);

	//Get rank of process
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Send the random number generator to get different results each time for each processor
	srand((int)time(NULL)*rank);
	GenerateRandomTree(POINTS_NUM);

	cout << "KD Tree Created by process " << rank << endl;
	
	// generate random point to query tree
	int pTarget[SD];
	float disKD;
	bool loop = true;

	while (loop) {

		//get the input file from the user, hard coding for speed testing
		/*string filename;
		cout << "Input the name of file with points(don't forget .txt): ";
		getline(cin, filename);*/
		ifstream points("blah.txt");

		int point;

		//put points in array
		if (points.is_open())
		{
			int i = 0;
			while (!points.eof())
			{
				points >> point;
				pTarget[i++] = point;
			}
			points.close();
			loop = false;
		}
		else cout << "Error: Could not open file.\n";
	}

	//output appropriate values to the console and file
	ofstream near_output("X_NN.txt");

	cout << "Points from file: ";
	for (int i = 0; i < SD; i++)
	{
		cout << pTarget[i] << " ";
		near_output << pTarget[i] << " ";
	}
	cout << endl;
	near_output << endl;

	KDNode<int>* nearest = Tree.find_nearest(pTarget);
	disKD = (float)sqrt(Tree.d_min);

	cout << "Nearest points: ";
	for (int i = 0; i < SD; i++)
	{
		cout << nearest->x[i] << " ";
		near_output << nearest->x[i] << " ";
	}
	
	near_output << endl;
	cout << endl;
	cout << "Distance: " << disKD << endl;
	near_output << disKD;

	clock_t end = clock();
	double ex_time = double(end - begin) / CLOCKS_PER_SEC;
	cout << "Exicution time in seconds: " << ex_time << endl;
	system("pause");
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}