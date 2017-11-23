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
#include <sstream>

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

	while (loop) 
	{

		ifstream points("blah.txt");

		string point;

		//put points in array
		if (points.is_open())
		{
			int i = 0;
			const char* fileName = "_NN";
			const char* fileType = ".txt";
			char name_buffer[512];
			while (!points.eof())
			{
				//get line as string then cast into int array
				getline(points, point);
				//cout << "Point " << i << ": " << point << endl;
				stringstream stream(point);
				for (int j = 0; j < SD; j++)
					stream >> pTarget[j];
				
				//generate correct file name
				sprintf(name_buffer, "%d%s%s", i, fileName, fileType);
				ofstream near_output(name_buffer);
				
				//put point in first line of file
				near_output << point << endl;

				//compute nearest and make it the second line in file
				KDNode<int>* nearest = Tree.find_nearest(pTarget);
				disKD = (float)sqrt(Tree.d_min);
				//cout << "Nearest point for point " << i << ": ";
				for (int i = 0; i < SD; i++)
				{
					//cout << nearest->x[i] << " ";
					near_output << nearest->x[i] << " ";
				}
				//cout << endl;
				near_output << endl;

				//put discance in third line of file
				//cout << "Distance: " << disKD << endl;
				near_output << disKD;

				i++;
			}
			points.close();
			loop = false;
		}
		else cout << "Error: Could not open file.\n";
	}

	//display exicution time
	clock_t end = clock();
	double ex_time = double(end - begin) / CLOCKS_PER_SEC;
	cout << "Exicution time in seconds: " << ex_time << endl;
	system("pause");
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}