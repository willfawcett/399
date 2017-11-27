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

#define POINTS_PER_MPI_MESSAGE 10; //for best results, keep this evenly divisible by number of points being evaluated 

struct parallel_point
{
	int x[SD];			//point from input file
	int y[SD];			//nearest point found in tree(s)
	double distance;	//distance between points
	int origin;			//starting rank
	int index;			//index in origin ranks's array of points
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

	// the next processor.  This is who we will be sending points to
	int next_rank = (rank + 1) % numProcesses;

	//the previous processor.  This is who we will be asking for work from
	int prev_rank = (rank - 1) % numProcesses;

	int points_per_message = POINTS_PER_MPI_MESSAGE; //for those times when the constant isn't gladly accepted

	cout << "numProcess:" << numProcesses << " rank:" << rank << " numPoints:" << numPoints << endl;
	MPI_Barrier(MPI_COMM_WORLD);

	/* GENERATE KD TREE */
	srand((int)time(NULL)*rank);
	GenerateRandomTree(POINTS_NUM);

	cout << "KD Tree Created by process " << rank << endl;

	/*READ DATA FILES*/
	string input_filename;
	input_filename = "..\\input_files\\eval_" + to_string(numPoints) + "_" + to_string(rank) + ".txt";

	cout << "loading file " << input_filename << endl;
	ifstream input(input_filename);
	MPI_Barrier(MPI_COMM_WORLD);

	//cout << input_filename << endl;

	cout << "loading points from file " << endl;
	//load points array from files
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
		point->origin = rank;
		point->index = lineNumber;
		lineNumber++;
		//cout << i++ << " : " << line << endl;
		//...for each line in input...
	}


	cout << "going to create datatype" << endl;

	//Set up the MPI Datatype
	MPI_Datatype parallel_point_mpi_t;
	int array_of_blocklengths[5] = {SD, SD, 1, 1, 1};
	int array_of_displacements[5];
	parallel_point *point = &points[0];
	MPI_Aint a_addr, b_addr, c_addr, d_addr, e_addr;
	MPI_Get_address(&point->x[0], &a_addr);
	array_of_displacements[0] = 0;
	MPI_Get_address(&point->y[0], &b_addr);
	array_of_displacements[1] = b_addr - a_addr;
	MPI_Get_address(&point->distance, &c_addr);
	array_of_displacements[2] = c_addr - b_addr;
	MPI_Get_address(&point->origin, &d_addr);
	array_of_displacements[3] = d_addr - c_addr;
	MPI_Get_address(&point->index, &e_addr);
	array_of_displacements[4] = e_addr - d_addr;
	MPI_Datatype array_of_types[5] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_INT};
	MPI_Type_create_struct(5, array_of_blocklengths, array_of_displacements, array_of_types, &parallel_point_mpi_t);
	MPI_Type_commit(&parallel_point_mpi_t);


	cout << "created datatype" << endl;
	MPI_Barrier(MPI_COMM_WORLD);


	//Loop Thru batches of POINTS_PER_MPI_MESSAGE
	int numberOfBatches = numPoints / POINTS_PER_MPI_MESSAGE;
	for (int batch = 0; batch < numberOfBatches; batch++) {
		
		//loop through the points in the batch and calculate distance
		for (int i = 0; i < points_per_message; i++) {
			//BOUNDS CHECKING SHOULD GO HERE
			int pointIndex = batch * POINTS_PER_MPI_MESSAGE + i;
			parallel_point* point = &points[pointIndex];

			//compute nearest and get distance
			KDNode<int>* nearest = Tree.find_nearest(point->x);
			double disKD = (double)sqrt(Tree.d_min);
			point->distance = disKD;

			for (int j = 0; j < SD; j++) {
				//BOUNDS CHECKING SHOULD GO HERE
				point->y[j] = nearest->x[j]; //copy nearest point from kd tree to our struct array
			}
		}

		//send the last POINTS_PER_MPI_MESSAGE points to the next processor
		int startingPointIndex = batch * POINTS_PER_MPI_MESSAGE;
		parallel_point* startingPoint = &points[startingPointIndex];
		MPI_Request request;
		MPI_Status sendStatus, recvStatus;
		parallel_point *incoming_points = new parallel_point[points_per_message];
		int result = MPI_Isend( &startingPoint, points_per_message, parallel_point_mpi_t, next_rank, 0, MPI_COMM_WORLD, &request);
		/*
		cout << result << endl;
		if (result == MPI_SUCCESS)
			cout << "No error; MPI routine completed successfully." << endl;
		else if (result == MPI_ERR_COMM)
			cout << "Invalid communicator.A common error is to use a null communicator in a call(not even allowed in MPI_Comm_rank).." << endl;
		else if (result == MPI_ERR_COUNT)
			cout << "Invalid count argument.Count arguments must be non - negative; a count of zero is often valid.." << endl;
		else if (result == MPI_ERR_TYPE)
			cout << "Invalid datatype argument.May be an uncommitted MPI_Datatype(see MPI_Type_commit).." << endl;
		else if (result == MPI_ERR_TAG)
			cout << "Invalid tag argument.Tags must be non - negative; tags in a receive(MPI_Recv, MPI_Irecv, MPI_Sendrecv, etc.) may also be MPI_ANY_TAG.The largest tag value is available through the the attribute MPI_TAG_UB.." << endl;
		else if (result == MPI_ERR_RANK)
			cout << "Invalid source or destination rank.Ranks must be between zero and the size of the communicator minus one; ranks in a receive(MPI_Recv, MPI_Irecv, MPI_Sendrecv, etc.) may also be MPI_ANY_SOURCE.." << endl;
		else if (result == MPI_ERR_INTERN)
			cout << "This error is returned when some part of the MPICH implementation is unable to acquire memory.." << endl;
		*/
		int recvresult = MPI_Recv(&incoming_points, SD, parallel_point_mpi_t, prev_rank, 0, MPI_COMM_WORLD, &recvStatus);
		cout << "receive: " << recvresult << endl;
		//MPI_Wait(&request, &sendStatus); //don't move on until we know our sent points were received

		int points_origin = -1;
		while (points_origin != rank) {
			//process the incoming points
			for (int i = 0; i < points_per_message; i++) {
				//BOUNDS CHECKING SHOULD GO HERE
				int pointIndex = batch * POINTS_PER_MPI_MESSAGE + i;
				parallel_point* point = &incoming_points[pointIndex];
				if (i == 0) {
					points_origin = point->origin; //we will exit the while loop after processing this set of points
				}

				cout << rank << " received points from " << point->origin << endl;;

				//compute nearest and get distance
				KDNode<int>* nearest = Tree.find_nearest(point->x);
				double disKD = (double)sqrt(Tree.d_min);

				if (disKD < point->distance && point->origin != rank) {
					//update the point in the incoming points
					point->distance = disKD;
					cout << "Updating foreign point distance to " << point->distance;

					for (int j = 0; j < SD; j++) {
						//BOUNDS CHECKING SHOULD GO HERE
						point->y[j] = nearest->x[j]; //copy nearest point from kd tree to our struct array
					}
				}
				else if (disKD < point->distance && point->origin == rank) {
					//point has travelled to all 4 processors and a better distance was found elsewhere
					//copy point data back into local points
					int localPointIndex = point->index;
					parallel_point *localPoint = &points[localPointIndex];
					localPoint->distance = point->distance;
					for (int j = 0; j < SD; j++) {
						//BOUNDS CHECKING SHOULD GO HERE
						localPoint->y[j] = nearest->x[j]; //copy nearest point from kd tree to our struct array
					}
				}
			}// end for each incomingpoints
			if (points_origin != rank) {
				//pass the points along to the next processor
				//reusing a lot of variables created just before entering the while loop

				//COMMENTED OUT MPI BELOW BECAUSE OF AN ERROR

				//MPI_Isend(&incoming_points, points_per_message, parallel_point_mpi_t, next_rank, 0, MPI_COMM_WORLD, &request);
				//MPI_Recv(&incoming_points, SD, parallel_point_mpi_t, prev_rank, 0, MPI_COMM_WORLD, &recvStatus);
			//	MPI_Wait(&request, &status); //don't move on until we know the points we sent were received
			}
		}

	}

	/*
	//LOOPS THROUGH ALL POINTS
	for (int i = 0; i < numPoints); i++) {
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
	}
	*/



	//Open Output File
	string output_file = "output_eval_" + to_string(numPoints) + "_" + to_string(rank) + ".txt";
	ofstream output(output_file);
		
	//OUTPUT TO FILE
	for (int i = 0; i < numPoints; i++) {
		//Output X (point from input file)
		for (int j = 0; j < SD; j++) {
			//BOUNDS CHECKING SHOULD GO HERE
			output << point->x[j] << " ";
		}
		output << endl;
		//Output Y (nearest point from kdtree)
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

	//Free the MPI datatype
	MPI_Type_free(&parallel_point_mpi_t);

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