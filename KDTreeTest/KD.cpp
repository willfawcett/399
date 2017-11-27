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

#define RANK_ZERO 1u << 0;
#define RANK_ONE 1u << 1;
#define RANK_TWO 1u << 2;
#define RANK_THREE 1u << 3;

#define POINTS_PER_MPI_MESSAGE = 10; //for best results, keep this evenly divisible by number of points being evaluated 

struct parallel_point
{
	int x[SD]; //point from input file
	int y[SD]; //nearest point found in tree(s)
	double distance; //distance between points
	unsigned int visited; //bitmask
} ;

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

	//the array of points we'll be working with
	parallel_point *points = new parallel_point[numPoints];

	cout << "numProcess:" << numProcesses << " rank:" << rank << " numPoints:" << numPoints << endl;
	MPI_Barrier(MPI_COMM_WORLD);

	/* GENERATE KD TREE */
	srand((int)time(NULL)*rank);
	GenerateRandomTree(POINTS_NUM);

	cout << "KD Tree Created by process " << rank << endl;

	/*READ DATA FILES*/
	string input_filename;
	input_filename = "..\\input_files\\eval_" + to_string(numPoints) + "_" + to_string(rank) + ".txt";
	ifstream input(input_filename);

	//cout << input_filename << endl;

	//load all points from input files for evaluation
	int lineNumber = 0;
	for (string line; getline(input, line); )
	{
		stringstream stream(line);
		//BOUNDS CHECKING SHOULD GO HERE
		parallel_point *point = &points[lineNumber];
		int j = 0;
		int n; while (stream >> n) { 
			//cout << "Found integer: " << n; 
			//BOUNDS CHECKING SHOULD GO HERE
			point->x[j++] = n;
		}
		lineNumber++;
		//cout << i++ << " : " << line << endl;
		//...for each line in input...
	}





	//Open Output File
	string output_file = "output_eval_" + to_string(numPoints) + "_" + to_string(rank) + ".txt";
	ofstream output(output_file);



	//Loop Thru ALL Points (debugging at this point)
	for (int i = 0; i < numPoints; i++) {
		//BOUNDS CHECKING SHOULD GO HERE
		parallel_point* point = &points[i];
		//compute nearest and get distance
		KDNode<int>* nearest = Tree.find_nearest(point->x);
		double disKD = (double)sqrt(Tree.d_min);
		point->distance = disKD;

		for (int j = 0; j < SD; j++) {
			//BOUNDS CHECKING SHOULD GO HERE
			 point->y[j] = nearest->x[j]; //copy nearest point from kd tree to our struct array
		}
		
		//OUTPUT TO FILE
		//Output X (point read in for evaluation)
		for (int j = 0; j < SD; j++) {
			//BOUNDS CHECKING SHOULD GO HERE
			output << point->x[j] << " ";
		}
		output << endl;
		//Output Y (point read in for evaluation)
		for (int j = 0; j < SD; j++) {
			//BOUNDS CHECKING SHOULD GO HERE
			output << point->y[j] << " ";
		}
		output << endl;
		//Output Distance
		output << point->distance << endl;
	}

	//close output file
	output.close();



	MPI_Barrier(MPI_COMM_WORLD);
	return;
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
	if (numProcesses == 4) {
		parallel_execution(numProcesses, rank, numPoints);
		MPI_Barrier(MPI_COMM_WORLD);
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
	MPI_Barrier(MPI_COMM_WORLD);
	clock_t end = clock();
	double ex_time = double(end - begin) / CLOCKS_PER_SEC;
	if (rank == 0) {
		cout << "Exicution time in seconds: " << ex_time << endl;
		//system("pause");
	}
	MPI_Finalize();

	return 0;
}