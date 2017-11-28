//CSC 399 MPI Project Sequential Tree
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
#include "parallel_version.h"

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

	// Initialize MPI environment
	MPI_Init(&argc, &argv);

	// Get number of processes
	int numProcesses;
	MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);

	//Get rank of process
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int serial = 1;

	if (rank == 0) {
		cout << "Input 0 for concurent version and 1 for sequential: ";
		cin >> serial;
	}

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

	//start time measurment after getting user input
	clock_t begin = clock();

	// Parallel or Serial (or exit)
	if (serial == 0) {
		if (numProcesses != 4) {
		//We only support 1 or 4 processes
		//Exit with error message
		if (rank == 0)
		{
			fprintf(stderr, "Please try again using either 1 or 4 processes.\n");
		}
		MPI_Barrier(MPI_COMM_WORLD); //synch processes
		system("pause");
		exit(1); //end app prematurely
		}
		else
		{
			parallel_execution(numProcesses, rank, numPoints);
			MPI_Barrier(MPI_COMM_WORLD);
		}
		
	}
		
	else
	{
		// SERIAL VERSION

		// Send the random number generator to get different results each time for each processor
		srand((int)time(NULL)*rank);
		GenerateRandomTree(POINTS_NUM);

		cout << "KD Tree Created by process " << rank << endl;

		// generate random point to query tree
		int pTarget[SD];
		float disKD;
		string infile_name;

		if (numPoints == 100)
			infile_name = "input_files//eval_100_SEQ.txt";
		else if (numPoints == 1000)
			infile_name = "input_files//eval_1000_SEQ.txt";
		else if (numPoints == 10000)
			infile_name = "input_files//eval_10000_SEQ.txt";
		else if (numPoints == 100000)
			infile_name = "input_files//eval_100000_SEQ.txt";

		ifstream points(infile_name);
		ofstream near_output("sequential_output.txt");

		if (points.is_open())
		{
			while (!points.eof())
			{
				string point;
				//get line as string then cast into int array
				getline(points, point);
				//cout << "Point " << i << ": " << point << endl;
				stringstream stream(point);
				for (int j = 0; j < SD; j++)
					stream >> pTarget[j];

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
				near_output << disKD << endl;
				}
			points.close();
			}
			else cout << "Error: Could not open file.\n";
		} 

	//display exicution time
	MPI_Barrier(MPI_COMM_WORLD);
	clock_t end = clock();
	double ex_time = double(end - begin) / CLOCKS_PER_SEC;
	if (rank == 0) {
		cout << "Exicution time in seconds: " << ex_time << endl;
		system("pause");
	}
	MPI_Finalize();

	return 0;
}