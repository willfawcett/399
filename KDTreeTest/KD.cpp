// KD.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include "KDTree.h"
#include <iostream>
#include "mpi.h"
#include <fstream>
#include <stdlib.h>
#include <string>

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

	// Send the random number generator to get different results each time for each processor
	srand(time(NULL)*rank);
	GenerateRandomTree(POINTS_NUM);

	cout << "KD Tree Created by process " << rank << endl;
	
	// generate random point to query tree
	int pTarget[SD];
	float disKD;
	bool loop = true;

	while (loop)
	{
		string filename;
		cout << "Input the name of file with points(don't forget .txt): ";
		getline(cin, filename);
		ifstream points(filename);

		int point;

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
	
	cout << "Points from file: ";
	for (int i = 0; i < SD; i++)
	{
		cout << pTarget[i] << " ";
	}
	cout << endl;

	ofstream near_output("nearest.txt");

	KDNode<int>* nearest = Tree.find_nearest(pTarget);
	disKD = (float)sqrt(Tree.d_min);

	//make out put one file
	//first line is original point
	//secdond line is nearest point
	//third line is distance between points

	for (int i = 0; i < SD; i++)
	{
		cout << nearest->x[i] << " ";
		near_output << nearest->x[i] << endl;
	}
	
	cout << endl;
	ofstream dist_output("distance.txt");
	cout << "Distance: " << disKD << endl;
	dist_output << disKD;
	system("pause");
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}