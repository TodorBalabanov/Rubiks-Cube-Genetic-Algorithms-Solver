#include <cmath>
#include <vector>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <ostream>
#include <sstream>
#include <iostream>

#include <mpi.h>
#include <unistd.h>

#include "Constants.h"
#include "RubiksCube.h"
#include "GeneticAlgorithm.h"
#include "GeneticAlgorithmOptimizer.h"

namespace std {
template < typename T > std::string to_string( const T& n ) {
	std::ostringstream stm;
	stm << n ;
	return stm.str() ;
}
}

static int rank = -1;
static int size = 0;

/*
 * Receive buffer.
 */
static char buffer[RECEIVE_BUFFER_SIZE];

static RubiksCube solved;

static void master() {
	unsigned long counter = 0;

	if(rank != ROOT_NODE) {
		return;
	}

	/*
	 * Cube to be solved.
	 */
	RubiksCube shuffled;
	shuffled.shuffle(CUBE_SHUFFLING_STEPS);
	std::cout << "Sender : " << std::to_string(shuffled.compare(solved)) << std::endl;

	std::map<int,GeneticAlgorithm> populations;

	/*
	 * Send shffled cube to all other nodes.
	 */{
		const std::string &value = shuffled.toString();
		for(int r=0; r<size; r++) {
			/*
			 * Root node is not included.
			 */
			if(r == ROOT_NODE) {
				continue;
			}

			MPI_Send(value.c_str(), value.size(), MPI_BYTE, r, DEFAULT_TAG, MPI_COMM_WORLD);
		}
	}

	do {
		std::cout << "Round : " << (counter+1) << std::endl;

		/*
		 * Send GA population to all other nodes.
		 */
		for(int r=0; r<size; r++) {
			/*
			 * Root node is not included.
			 */
			if(r == ROOT_NODE) {
				continue;
			}

			GeneticAlgorithm ga;
			if(counter == 0) {
				GeneticAlgorithmOptimizer::addEmptyCommand(ga, solved, shuffled);
				GeneticAlgorithmOptimizer::addRandomCommands(ga, solved, shuffled, LOCAL_POPULATION_SIZE);
				populations[r] = ga;
			} else {
				/*
				 * Ring migration strategy.
				 */
				int next = (r+1) % size;
				while(next == ROOT_NODE) {
					next = (next+1) % size;
				}
				populations[r].replaceWorst(populations[next].getBestChromosome());
			}
			const std::string &value = populations[r].toString();
			MPI_Send(value.c_str(), value.size(), MPI_BYTE, r, DEFAULT_TAG, MPI_COMM_WORLD);
		}

		/*
		 * Collect results from all other nodes.
		 */
		for(int r=0; r<size; r++) {
			/*
			 * Root node is not included.
			 */
			if(r == ROOT_NODE) {
				continue;
			}

			GeneticAlgorithm ga;
			MPI_Recv(buffer, RECEIVE_BUFFER_SIZE, MPI_BYTE, r, DEFAULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			ga.fromString(buffer);
			populations[r] = ga;
			std::cout << "Worker " << r << " : " << ga.getBestChromosome().fitness << std::endl;
		}

		counter++;
	} while(counter < NUMBER_OF_BROADCASTS);
}

static void slave() {
	unsigned long counter = 0;

	if(rank == ROOT_NODE) {
		return;
	}

	RubiksCube shuffled;
	MPI_Recv(buffer, RECEIVE_BUFFER_SIZE, MPI_BYTE, ROOT_NODE, DEFAULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	shuffled.fromString(buffer);

	do {
		GeneticAlgorithm ga;
		MPI_Recv(buffer, RECEIVE_BUFFER_SIZE, MPI_BYTE, ROOT_NODE, DEFAULT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		ga.fromString(buffer);

		/*
		 * Calculate as regular node.
		 */
		GeneticAlgorithmOptimizer::optimize(ga, solved, shuffled, LOCAL_OPTIMIZATION_EPOCHES);

		std::string result = ga.toString();
		MPI_Send(result.c_str(), result.size(), MPI_BYTE, ROOT_NODE, DEFAULT_TAG, MPI_COMM_WORLD);

		counter++;
	} while(counter < NUMBER_OF_BROADCASTS);
}

int main(int argc, char **argv) {
	MPI_Init (&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	srand( time(NULL)^getpid() );

	/*
	 * Firs process will distribute the working tasks.
	 */
	master();
	slave();

	MPI_Finalize();

	return( EXIT_SUCCESS );
}
