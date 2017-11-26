//CSC 399 MPI Project Sequential KDD Tree
//William Fawcet, Justin Mizelle
//This program randomly generates a KDD tree, takes in points from a txt file,
//then computes the nearest points and distance and outputs this to a txt file for each point.

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

void parallel_execution(int &numProcesses, int &rank, int &numPoints) {
	//init and finalize happen outside this function
	//MPI 4 processer version goes here


	cout << "numProcess:" << numProcesses << " rank:" << rank << " numPoints:" << numPoints << endl;

	/*READ DATA FILES*/



	
	return;
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

	int numPoints = 0;
	//find out how many points we want to process (i.e. which data file(s) to use)
	if (rank == 0) {
		//check input
		while (numPoints != 100 && numPoints != 1000 && numPoints != 10000 && numPoints != 100000) {
			cout << "Enter number of points -- 100, 1000, 10000, or 100000: ";
			cin >> numPoints;
		}
	}
	if (numProcesses > 1) {
		//make sure all processors have the same value in numPoints as 0
		MPI_Bcast(&numPoints, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}

	// Parallel or Serial (or exit)
	if (numProcesses == 4) {
		parallel_execution(numProcesses, rank, numPoints);
	}
	else if (numProcesses != 1) {
		//We only support 1 or 4 processes
		//Exit with error message
		if(rank == 0)
			fprintf(stderr, "Please try again using either 1 or 4 processes.");
		MPI_Barrier(MPI_COMM_WORLD); //synch processes
		MPI_Finalize(); //exit MPI gracefully
		exit(1); //end app prematurely
	}
	else {
		// SERIAL VERSION

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

			//input file is hard coded to keep I/O from affecting exicution time
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
		} // end while loop
	}// END if/elseif/else

	//display exicution time
	clock_t end = clock();
	double ex_time = double(end - begin) / CLOCKS_PER_SEC;
	if (rank == 0) {
		cout << "Exicution time in seconds: " << ex_time << endl;
		//system("pause");
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}