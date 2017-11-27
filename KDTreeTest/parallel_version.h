#pragma once
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

#define TESTING true

#define POINTS_PER_MPI_MESSAGE 10 //for best results, keep this evenly divisible by number of points being evaluated 

/*
struct parallel_point
{
int x[SD];			//point from input file
int y[SD];			//nearest point found in tree(s)
double distance;	//distance between points
int origin;			//starting rank
int index;			//index in origin ranks's array of points
} ;
*/

typedef struct parallel_point_s
{
	int x[SD];			//point from input file
	int y[SD];			//nearest point found in tree(s)
	double distance;	//distance between points
	int origin;			//starting rank
	int index;			//index in origin ranks's array of points
} parallel_point;

KDTree<int> Tree_parallel;

//duplicated and renamed to avoid conflict
void generate_random_point_parallel(int* p, int sd)
{
	for (int k = 0; k<sd; k++)
		p[k] = rand() % RMAX;
}

void GenerateRandomTree_parallel(int nNodes)
{
	int p[SD]; // random point
	for (int n = 0; n<nNodes; n++)
	{
		generate_random_point_parallel(p, SD);
		Tree_parallel.add(p);
	}
}

void print_point(parallel_point &point, int &rank) {
	cout << rank << ": " << " distance " << point.distance << " origin " << point.origin << " index " << point.index << endl;
	cout << rank << " x :";
	for (int i = 0; i < SD; i++)
		cout << " " << point.x[i];
	cout << endl;
	cout << rank << " y :";
	for (int i = 0; i < SD; i++)
		cout << " " << point.y[i];
	cout << endl;
}

void parallel_execution(int &numProcesses, int &rank, int &numPoints) {
	//init and finalize happen outside this function
	//MPI 4 processer version goes here


	//the array of points we'll be working with
	parallel_point *points = new parallel_point[numPoints];

	// the next processor.  This is who we will be sending points to
	int next_rank = (rank + 1) % numProcesses;

	//the previous processor.  This is who we will be asking for work from
	int prev_rank = (rank + numProcesses - 1) % numProcesses;

	cout << "rank: " << rank << " next_rank: " << next_rank << " prev_rank " << prev_rank << endl;
	//system("pause");
	//MPI_Barrier(MPI_COMM_WORLD);

	int points_per_message = POINTS_PER_MPI_MESSAGE; //for those times when the constant isn't gladly accepted

	cout << "numProcess:" << numProcesses << " rank:" << rank << " numPoints:" << numPoints << endl;
	MPI_Barrier(MPI_COMM_WORLD);

	/* GENERATE KD TREE */
	srand((int)time(NULL)*rank);
	GenerateRandomTree_parallel(POINTS_NUM);

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



	//Set up the MPI Datatype
	MPI_Datatype parallel_point_mpi_t;
	int array_of_blocklengths[5] = { SD, SD, 1, 1, 1 };
	int array_of_displacements[5];
	/*
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
	*/
	array_of_displacements[0] = offsetof(parallel_point, x);
	array_of_displacements[1] = offsetof(parallel_point, y);
	array_of_displacements[2] = offsetof(parallel_point, distance);
	array_of_displacements[3] = offsetof(parallel_point, origin);
	array_of_displacements[4] = offsetof(parallel_point, index);
	MPI_Datatype array_of_types[5] = { MPI_INT, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_INT };
	MPI_Type_create_struct(5, array_of_blocklengths, array_of_displacements, array_of_types, &parallel_point_mpi_t);
	MPI_Type_commit(&parallel_point_mpi_t);


	cout << "created datatype" << endl;
	MPI_Barrier(MPI_COMM_WORLD);



	//BROADCAST A SINGLE POINT AS A TEST
	if (TESTING) {
		parallel_point bcastPoint;

		if (rank == 0) {
			bcastPoint.distance = (double)100;
			for (int i = 0; i < SD; i++)
				bcastPoint.x[i] = i;
			for (int i = 0; i < SD; i++)
				bcastPoint.y[i] = i + SD;
			bcastPoint.origin = 21;
			bcastPoint.index = 22;


			//DEBUG OUTPUT
			cout << rank << ": " << " distance " << bcastPoint.distance << " origin " << bcastPoint.origin << " index " << bcastPoint.index << endl;
			cout << rank << " x :";
			for (int i = 0; i < SD; i++)
				cout << " " << bcastPoint.x[i];
			cout << endl;
			cout << rank << " y :";
			for (int i = 0; i < SD; i++)
				cout << " " << bcastPoint.y[i];
			cout << endl;

		}

		MPI_Bcast(&bcastPoint, 1, parallel_point_mpi_t, 0, MPI_COMM_WORLD);

		//DEBUG OUTPUT
		cout << rank << ": " << " distance " << bcastPoint.distance << " origin " << bcastPoint.origin << " index " << bcastPoint.index << endl;
		cout << rank << " x :";
		for (int i = 0; i < SD; i++)
			cout << " " << bcastPoint.x[i];
		cout << endl;
		cout << rank << " y :";
		for (int i = 0; i < SD; i++)
			cout << " " << bcastPoint.y[i];
		cout << endl;

		MPI_Barrier(MPI_COMM_WORLD);
		//MPI_Finalize();
		//exit(1);
	}




	//HOW IT WORKS
	//1 ) calculate a set of local points (number of points being - POINTS_PER_MPI_MESSAGE)
	//2 ) send that set of points to the next processor
	//3 ) wait for the previous processor to send us points
	//4 ) wait for confirmation that our points were received
	//5 ) calculate the received points
	//6.1 ) if these are not our points, calculate against our local KDTree and update the foreign points with nearest distances found
	//6.2 ) if these are our points, update the local struct with the nearest distances found
	//7.1 ) if they aren't our points, forward the set of points to next processor (basically 2, 3, 4 again)
	//7.2 ) do any cleanup and then return to 1


	//Loop Thru batches of POINTS_PER_MPI_MESSAGE
	int numberOfBatches = numPoints / POINTS_PER_MPI_MESSAGE;
	for (int batch = 0; batch < numberOfBatches; batch++) {

		//#1 loop through the points in the batch and calculate distance
		for (int i = 0; i < points_per_message; i++) {
			//BOUNDS CHECKING SHOULD GO HERE
			int pointIndex = batch * POINTS_PER_MPI_MESSAGE + i;
			parallel_point* point = &points[pointIndex];

			//compute nearest and get distance
			KDNode<int>* nearest = Tree_parallel.find_nearest(point->x);
			double disKD = (double)sqrt(Tree_parallel.d_min);
			point->distance = disKD;

			//copy nearest point from kd tree to our struct array
			for (int j = 0; j < SD; j++) {
				//BOUNDS CHECKING SHOULD GO HERE
				point->y[j] = nearest->x[j]; 
			}
		}


		//send the last POINTS_PER_MPI_MESSAGE points to the next processor
		int startingPointIndex = batch * POINTS_PER_MPI_MESSAGE;
		parallel_point* startingPoint = &points[startingPointIndex];
		MPI_Request request;
		MPI_Status sendStatus, recvStatus;
		parallel_point points_buffer_a[POINTS_PER_MPI_MESSAGE];
		parallel_point points_buffer_b[POINTS_PER_MPI_MESSAGE]; //buffer used later on for forwarding 
		bool processing_points_buffer_a = true;
		//#2 USING ISEND so all procs can send data without deadlock
		int result = MPI_Isend(startingPoint, points_per_message, parallel_point_mpi_t, next_rank, 0, MPI_COMM_WORLD, &request);
		//#3 USING BLOCKING RECEIVE because we need to keep things synchronized
		int recvresult = MPI_Recv(&points_buffer_a, points_per_message, parallel_point_mpi_t, prev_rank, 0, MPI_COMM_WORLD, &recvStatus);
		//#4 We'll be altering the buffers, so we need to make sure communication is complted
		MPI_Wait(&request, &sendStatus); //don't move on until we know our sent points were received

		int points_origin = -1;
		while (points_origin != rank) {
			//process the incoming points
			for (int i = 0; i < points_per_message; i++) {
				//BOUNDS CHECKING SHOULD GO HERE
				int pointIndex = batch * POINTS_PER_MPI_MESSAGE + i;

				parallel_point* point;
				if(processing_points_buffer_a)
					point = &points_buffer_a[pointIndex];
				else
					point = &points_buffer_b[pointIndex];

				//update points_origin with data from the first point -- they will all be of the same origin
				if (i == 0) {
					points_origin = point->origin; //we will exit the while loop after processing this set of points
				}

				cout << rank << " received points from " << point->origin << endl;;

				//compute nearest and get distance
				KDNode<int>* nearest = Tree_parallel.find_nearest(point->x);
				double disKD = (double)sqrt(Tree_parallel.d_min);

				if (disKD < point->distance && point->origin != rank) {
					//update the point in the incoming points
					point->distance = disKD;
					cout << rank << " Updating foreign point distance to " << point->distance;

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

					cout << rank << " Updating local point distance to " << point->distance;
					for (int j = 0; j < SD; j++) {
						//BOUNDS CHECKING SHOULD GO HERE
						localPoint->y[j] = nearest->x[j]; //copy nearest point from kd tree to our struct array
					}
				}
			}// end for each incomingpoints
			if (points_origin != rank) {
				// #7.1
				//pass the points along to the next processor
				//we can reuse the receiving buffer as the sending buffer, but we will need a new receiving buffer
				

				//ALTERNATE SEND/RECEIVE BUFFERS
				if (processing_points_buffer_a) {
					//reusing request/recv/status variables
					//we're going to send buffer_a and load into buffer b
					MPI_Isend(&points_buffer_a, POINTS_PER_MPI_MESSAGE, parallel_point_mpi_t, next_rank, 0, MPI_COMM_WORLD, &request);
					MPI_Recv(&points_buffer_b, SD, parallel_point_mpi_t, prev_rank, 0, MPI_COMM_WORLD, &recvStatus);
					MPI_Wait(&request, &sendStatus); //don't move on until we know the points we sent were received

					//let the next iteration know to process buffer b
					processing_points_buffer_a = false; //we've loaded the next set into points_buffer_b
				}
				else {
					//reusing request/recv/status variables - 
					//we're going to send buffer_b and load into buffer a
					MPI_Isend(&points_buffer_b, POINTS_PER_MPI_MESSAGE, parallel_point_mpi_t, next_rank, 0, MPI_COMM_WORLD, &request);
					MPI_Recv(&points_buffer_a, SD, parallel_point_mpi_t, prev_rank, 0, MPI_COMM_WORLD, &recvStatus);
					MPI_Wait(&request, &sendStatus); //don't move on until we know the points we sent were received

					//let the next iteration know to process buffer a
					processing_points_buffer_a = true; //we've loaded the next set into points_buffer_b
				}
			}
			else {
				//we're exiting the loop. need to do any cleanup?
				//#7.2
				cout << "We've received our own points" << endl;
				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Finalize();
				exit(0);
			}
		}

	}

	/*
	//LOOPS THROUGH ALL POINTS
	for (int i = 0; i < numPoints); i++) {
	//BOUNDS CHECKING SHOULD GO HERE
	parallel_point* point = &points[i];

	//compute nearest and get distance
	KDNode<int>* nearest = Tree_parallel.find_nearest(point->x);
	double disKD = (double)sqrt(Tree_parallel.d_min);
	point->distance = disKD;

	for (int j = 0; j < SD; j++) {
	//BOUNDS CHECKING SHOULD GO HERE
	point->y[j] = nearest->x[j]; //copy nearest point from kd tree to our struct array
	}
	}
	*/



	//Open Output File
	/*
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
	*/



	MPI_Barrier(MPI_COMM_WORLD);

	//Free the MPI datatype
	MPI_Type_free(&parallel_point_mpi_t);

	return;
}