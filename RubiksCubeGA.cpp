#include <cmath>
#include <vector>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <ostream>
#include <iostream>

#include <mpi.h>
#include <unistd.h>

enum RubiksColor {
	GREEN = 1,
	PURPLE = 2,
	YELLOW = 3,
	RED = 4,
	WHITE = 5,
	BLUE = 6
};

enum RotationDirection {
	COUNTERCLOCKWISE = -1,
	CLOCKWISE = +1
};

enum RubiksSide {
	TOP = 'T',
	LEFT = 'L',
	BACK = 'B',
	RIGHT = 'R',
	FRONT = 'F',
	DOWN = 'D',
};

class RubiksCube {
private:
	int top[3][3];
	int left[3][3];
	int right[3][3];
	int front[3][3];
	int back[3][3];
	int down[3][3];

	void spinSide(RubiksSide side) {
		static int buffer[ 3 ];

		if (side == TOP) {
			for (int i = 0; i < 3; i++) {
				buffer[i] = left[i][2];
			}
			for (int i = 0; i < 3; i++) {
				left[i][2] = front[0][i];
			}
			for (int i = 0; i < 3; i++) {
				front[0][i] = right[3 - i - 1][0];
			}
			for (int i = 0; i < 3; i++) {
				right[i][0] = back[2][i];
			}
			for (int i = 0; i < 3; i++) {
				back[2][3 - i - 1] = buffer[i];
			}
		} else if (side == LEFT) {
			for (int i = 0; i < 3; i++) {
				buffer[i] = down[i][2];
			}
			for (int i = 0; i < 3; i++) {
				down[3 - i - 1][2] = front[i][0];
			}
			for (int i = 0; i < 3; i++) {
				front[i][0] = top[i][0];
			}
			for (int i = 0; i < 3; i++) {
				top[i][0] = back[i][0];
			}
			for (int i = 0; i < 3; i++) {
				back[3 - i - 1][0] = buffer[i];
			}
		} else if (side == BACK) {
			for (int i = 0; i < 3; i++) {
				buffer[i] = down[0][i];
			}
			for (int i = 0; i < 3; i++) {
				down[0][i] = left[0][i];
			}
			for (int i = 0; i < 3; i++) {
				left[0][i] = top[0][i];
			}
			for (int i = 0; i < 3; i++) {
				top[0][i] = right[0][i];
			}
			for (int i = 0; i < 3; i++) {
				right[0][i] = buffer[i];
			}
		} else if (side == RIGHT) {
			for (int i = 0; i < 3; i++) {
				buffer[i] = down[i][0];
			}
			for (int i = 0; i < 3; i++) {
				down[i][0] = back[3 - i - 1][2];
			}
			for (int i = 0; i < 3; i++) {
				back[i][2] = top[i][2];
			}
			for (int i = 0; i < 3; i++) {
				top[i][2] = front[i][2];
			}
			for (int i = 0; i < 3; i++) {
				front[3 - i - 1][2] = buffer[i];
			}
		} else if (side == FRONT) {
			for (int i = 0; i < 3; i++) {
				buffer[i] = down[2][i];
			}
			for (int i = 0; i < 3; i++) {
				down[2][i] = right[2][i];
			}
			for (int i = 0; i < 3; i++) {
				right[2][i] = top[2][i];
			}
			for (int i = 0; i < 3; i++) {
				top[2][i] = left[2][i];
			}
			for (int i = 0; i < 3; i++)
				left[2][i] = buffer[i];
		} else if (side == DOWN) {
			for (int i = 0; i < 3; i++) {
				buffer[i] = front[2][i];
			}
			for (int i = 0; i < 3; i++) {
				front[2][i] = left[i][0];
			}
			for (int i = 0; i < 3; i++) {
				left[i][0] = back[0][3 - i - 1];
			}
			for (int i = 0; i < 3; i++) {
				back[0][i] = right[i][2];
			}
			for (int i = 0; i < 3; i++) {
				right[3 - i - 1][2] = buffer[i];
			}
		}
	}

	void spinClockwise(int side[3][3], int times, RubiksSide index) {
		static int buffer[3][3];
		static int newarray[3][3];

		if (times == 0) {
			return;
		}

		/*
		 * Transponse.
		 */
		for (int j = 0; j < 3; j++) {
			for (int i = 0; i < 3; i++) {
				newarray[j][i] = side[i][j];
			}
		}
		/*
		 * Rearrange.
		 */
		for (int i = 0; i < 3; i++) {
			static int cache = 0;
			cache = newarray[i][0];
			newarray[i][0] = newarray[i][2];
			newarray[i][2] = cache;
		}

		spinSide(index);
		memcpy(buffer, newarray, sizeof(int)*3*3);

		for (int t = 1; t < times; t++) {
			for (int j = 0; j < 3; j++) {
				for (int i = 0; i < 3; i++) {
					newarray[j][i] = buffer[i][j];
				}
			}
			for (int i = 0; i < 3; i++) {
				static int cache = 0;
				cache = newarray[i][0];
				newarray[i][0] = newarray[i][2];
				newarray[i][2] = cache;
			}

			spinSide(index);

			memcpy(buffer, newarray, sizeof(int)*3*3);
		}

		memcpy(side, buffer, sizeof(int)*3*3);
	}

	friend std::ostream& operator<< (std::ostream &out, const RubiksCube &cube);

public:
	RubiksCube() {
		reset();
	}

	void reset() {
		for(int i=0; i<3; i++) {
			for(int j=0; j<3; j++) {
				top[i][j] = GREEN;
				left[i][j] = PURPLE;
				right[i][j] = RED;
				front[i][j] = WHITE;
				back[i][j] = YELLOW;
				down[i][j] = BLUE;
			}
		}
	}

	long compare(const RubiksCube &cube) const {
		/**/
		long difference = 0;

		for(int i=0; i<3; i++) {
			for(int j=0; j<3; j++) {
				difference += abs(top[i][j]-cube.top[i][j]);
				difference += abs(left[i][j]-cube.left[i][j]);
				difference += abs(right[i][j]-cube.right[i][j]);
				difference += abs(front[i][j]-cube.front[i][j]);
				difference += abs(back[i][j]-cube.back[i][j]);
				difference += abs(down[i][j]-cube.down[i][j]);
			}
		}

		return(difference);
		/**/

		//TODO Find better distance measure (for example Hausdorff distance).
		/*
		long ha = 0;
		long hb = 0;
		long result = 0;

		for(int m=0; m<3; m++) {
			for(int n=0; n<3; n++) {
				int distances[] = {0, 0, 0, 0, 0, 0, 0, 0, 0 , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

				for(int i=0, d=0; i<3; i++) {
					for(int j=0; j<3; j++) {
						distances[d++] = abs(top[m][n]-cube.top[i][j]);
						distances[d++] = abs(left[m][n]-cube.left[i][j]);
						distances[d++] = abs(right[m][n]-cube.right[i][j]);
						distances[d++] = abs(front[m][n]-cube.front[i][j]);
						distances[d++] = abs(back[m][n]-cube.back[i][j]);
						distances[d++] = abs(down[m][n]-cube.down[i][j]);
					}
				}

				int min = distances[0];
				for(int d=0; d<54; d++) {
					if(distances[d] < min) {
						min = distances[d];
					}
				}

				if(min > ha) {
					ha = min;
				}
			}
		}

		for(int m=0; m<3; m++) {
			for(int n=0; n<3; n++) {
				int distances[] = {0, 0, 0, 0, 0, 0, 0, 0, 0 , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

				for(int i=0, d=0; i<3; i++) {
					for(int j=0; j<3; j++) {
						distances[d++] = abs(top[i][j]-cube.top[m][n]);
						distances[d++] = abs(left[i][j]-cube.left[m][n]);
						distances[d++] = abs(right[i][j]-cube.right[m][n]);
						distances[d++] = abs(front[i][j]-cube.front[m][n]);
						distances[d++] = abs(back[i][j]-cube.back[m][n]);
						distances[d++] = abs(down[i][j]-cube.down[m][n]);
					}
				}

				int min = distances[0];
				for(int d=0; d<54; d++) {
					if(distances[d] < min) {
						min = distances[d];
					}
				}

				if(min > hb) {
					hb = min;
				}
			}
		}

		result = std::max(ha, hb);

		return(result);
		/**/
	}

	void callSpin(RubiksSide side, RotationDirection direction, int numberOfTimes) {
		if (numberOfTimes < 0) {
			numberOfTimes = -numberOfTimes;
			if(direction == CLOCKWISE) {
				direction = COUNTERCLOCKWISE;
			} else if(direction == COUNTERCLOCKWISE) {
				direction = CLOCKWISE;
			}
		}

		numberOfTimes %= 4;

		if (direction == CLOCKWISE) {
			if (side == TOP) {
				spinClockwise(top, numberOfTimes, TOP);
			}
			if (side == LEFT) {
				spinClockwise(left, numberOfTimes, LEFT);
			}
			if (side == RIGHT) {
				spinClockwise(right, numberOfTimes, RIGHT);
			}
			if (side == FRONT) {
				spinClockwise(front, numberOfTimes, FRONT);
			}
			if (side == BACK) {
				spinClockwise(back, numberOfTimes, BACK);
			}
			if (side == DOWN) {
				spinClockwise(down, numberOfTimes, DOWN);
			}
		}
	}

	void execute(std::string commands) {
		for(int i=0; i<commands.length(); i++) {
			callSpin((RubiksSide)commands[i], CLOCKWISE, 1);
		}
	}

	std::string shuffle(int numberOfMoves=0) {
		std::string commands;

		for(int i=0; i<numberOfMoves; i++) {
			switch(rand()%6) {
			case 0:
				commands+=(char)TOP;
				break;
			case 1:
				commands+=(char)LEFT;
				break;
			case 2:
				commands+=(char)RIGHT;
				break;
			case 3:
				commands+=(char)FRONT;
				break;
			case 4:
				commands+=(char)BACK;
				break;
			case 5:
				commands+=(char)DOWN;
				break;
			}
		}

		execute(commands);

		return(commands);
	}
};

std::ostream& operator<< (std::ostream &out, const RubiksCube &cube) {
	for(int i=0; i<3; i++) {
		out << "      ";
		for(int j=0; j<3; j++) {
			out << cube.back[i][j] << " ";
		}
		out << std::endl;
	}

	for(int i=0; i<3; i++) {
		for(int j=0; j<3; j++) {
			out << cube.left[i][j] << " ";
		}
		for(int j=0; j<3; j++) {
			out << cube.top[i][j] << " ";
		}
		for(int j=0; j<3; j++) {
			out << cube.right[i][j] << " ";
		}
		for(int j=0; j<3; j++) {
			out << cube.down[i][j] << " ";
		}
		out << std::endl;
	}

	for(int i=0; i<3; i++) {
		out << "      ";
		for(int j=0; j<3; j++) {
			out << cube.front[i][j] << " ";
		}
		out << std::endl;
	}
}

class GeneticAlgorithm {
private:
	std::vector<std::string> population;
	std::vector<int> fitness;
	int resultIndex;
	int firstIndex;
	int secondIndex;
	int bestIndex;
	int worstIndex;

	void selectRandom() {
		do {
			resultIndex = rand() % population.size();
			firstIndex = rand() % population.size();
			secondIndex = rand() % population.size();
		} while(resultIndex==firstIndex || resultIndex==secondIndex || population[firstIndex].length()==0 || population[secondIndex].length()==0);
	}

	friend std::ostream& operator<< (std::ostream &out, const GeneticAlgorithm &ga);

public:
	static const int INVALID_FITNESS_VALUE = -1;
	static const bool KEEP_ELITE = true;

public:
	GeneticAlgorithm(int populationSize=0) {
		if(populationSize < 0) {
			populationSize = 0;
		}
		population.resize(populationSize);
		fitness.resize(populationSize);
		resultIndex = 0;
		firstIndex = 0;
		secondIndex = 0;
		bestIndex = 0;
		worstIndex = 0;
	}

	int getResultIndex() {
		return( resultIndex );
	}

	int getBestIndex() {
		return( bestIndex );
	}

	void setChromosome(std::string chromosome, int index=-1) {
		if(index < -1) {
			return;
		}

		if(index == -1) {
			population.push_back( chromosome );
			//fitness.push_back( GeneticAlgorithm::INVALID_FITNESS_VALUE );
			fitness.push_back( -1 );
			index = population.size() - 1;
		} else if(index < population.size()) {
			population[index] = chromosome;
		}

		//TODO Calculate fitness.
	}

	const std::string& getChromosome(int index) const {
		const static std::string EMPTY_STRING = "";

		if(population.size() <= index || index <= -1) {
			return( EMPTY_STRING );
		}

		return( population[index] );
	}

	void setFitness(int fitness, int index) {
		if(population.size() <= index || index <= -1) {
			return;
		}

		this->fitness[index] = fitness;
		if(fitness < this->fitness[bestIndex]) {
			bestIndex = index;
		}
		if(fitness > this->fitness[worstIndex]) {
			worstIndex = index;
		}
	}

	int getFitness(int index) {
		if(population.size() <= index || index <= -1) {
			return( INVALID_FITNESS_VALUE );
		}

		return( fitness[index] );
	}

	int size() {
		return( population.size() );
	}

	void selection() {
		static const int CROSSOVER_RESULT_INTO_BEST_PERCENT = 1;
		static const int CROSSOVER_RESULT_INTO_MIDDLE_PERCENT = 9;
		static const int CROSSOVER_RESULT_INTO_WORST_PERCENT = 90;

		static int percent = -1;
		percent = rand()
				  % (CROSSOVER_RESULT_INTO_WORST_PERCENT
					 + CROSSOVER_RESULT_INTO_MIDDLE_PERCENT
					 + CROSSOVER_RESULT_INTO_BEST_PERCENT);

		if (percent < CROSSOVER_RESULT_INTO_WORST_PERCENT) {
			do {
				selectRandom();
			} while (fitness[resultIndex] < fitness[firstIndex]
					 || fitness[resultIndex] < fitness[secondIndex]);
		} else if (percent
				   < (CROSSOVER_RESULT_INTO_WORST_PERCENT
					  + CROSSOVER_RESULT_INTO_MIDDLE_PERCENT)) {
			if (fitness[secondIndex] < fitness[firstIndex]) {
				int index = secondIndex;
				secondIndex = firstIndex;
				firstIndex = index;
			}
			do {
				selectRandom();
			} while (fitness[resultIndex] < fitness[firstIndex]
					 || fitness[resultIndex] > fitness[secondIndex]);
		} else if (percent
				   < (CROSSOVER_RESULT_INTO_WORST_PERCENT
					  + CROSSOVER_RESULT_INTO_MIDDLE_PERCENT
					  + CROSSOVER_RESULT_INTO_BEST_PERCENT)) {
			do {
				selectRandom();
			} while (fitness[resultIndex] > fitness[firstIndex]
					 || fitness[resultIndex] > fitness[secondIndex]);
		}

		if (resultIndex == bestIndex && KEEP_ELITE==true) {
			resultIndex = worstIndex;
		}
	}

	void crossover() {
		population[resultIndex] = population[firstIndex].substr(0, rand()%(population[firstIndex].length())+1);
		population[resultIndex] += population[secondIndex].substr(rand()%population[secondIndex].length(), population[secondIndex].length());
		fitness[resultIndex] = INVALID_FITNESS_VALUE;
	}

	void mutation() {
		int index = rand() % population[resultIndex].length();

		switch(rand()%6) {
		case 0:
			population[resultIndex][index]=(char)TOP;
			break;
		case 1:
			population[resultIndex][index]=(char)LEFT;
			break;
		case 2:
			population[resultIndex][index]=(char)RIGHT;
			break;
		case 3:
			population[resultIndex][index]=(char)FRONT;
			break;
		case 4:
			population[resultIndex][index]=(char)BACK;
			break;
		case 5:
			population[resultIndex][index]=(char)DOWN;
			break;
		}

		fitness[resultIndex] = INVALID_FITNESS_VALUE;
	}
};

std::ostream& operator<< (std::ostream &out, const GeneticAlgorithm &ga) {
	for(int p=0; p<ga.population.size(); p++) {
		out << ga.fitness[p];
		out << "\t";
		for(int i=0; i<ga.population[p].length(); i++) {
			out << ga.population[p][i];
		}
		out << std::endl;
	}
}

class GeneticAlgorithmOptimizer {
private:
	static int evaluate(const RubiksCube &solved, const RubiksCube &shuffled, const std::string &commands) {
		static RubiksCube used;

		used = shuffled;
		used.execute(commands);
		return( solved.compare(used) );
	}

	GeneticAlgorithmOptimizer() {
	}

public:
	static void optimize(GeneticAlgorithm &ga, RubiksCube &solved, RubiksCube &shuffled, int populationSize=0, long epoches=0) {
		for(int p=0; p<populationSize; p++) {
			RubiksCube mixed = solved;
			mixed.shuffle(1000);
			ga.setChromosome(mixed.shuffle(1000));
			ga.setFitness(evaluate(solved, shuffled, ga.getChromosome(p)), p);
		}

		for(long e=0L; e<epoches*ga.size(); e++) {
			ga.selection();
			ga.crossover();
			ga.mutation();
			int index = ga.getResultIndex();
			ga.setFitness(evaluate(solved, shuffled, ga.getChromosome(index)), index);
		}

		shuffled.execute(ga.getChromosome(ga.getBestIndex()));
	}
};

int main(int argc, char **argv) {
	int rank, size;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &size);

	srand( time(NULL)^getpid() );

	/*
	 * Firs process will distribute the working tasks.
	 */
	if(rank == 0) {
		RubiksCube solved, shuffled;
		shuffled.shuffle(10000);
		std::cout << "Sender difference: " << shuffled.compare(solved);
		std::cout << std::endl;
		for(int destination=1; destination<rank; destination++) {
			//MPI_Send(&buffer,count, MPI_PACKED, destination, 0, MPI_COMM_WORLD);
		}
	} else if(rank > 0) {
		GeneticAlgorithm ga;
		RubiksCube solved, shuffled;
		shuffled.shuffle(3);
		GeneticAlgorithmOptimizer::optimize(ga, solved, shuffled, 100, 10000);
		std::cout << "Worker " << rank << " : " << shuffled.compare(solved) << std::endl;
	}

	MPI_Finalize();

	return( EXIT_SUCCESS );
}

