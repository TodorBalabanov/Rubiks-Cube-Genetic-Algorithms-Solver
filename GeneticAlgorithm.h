#ifndef GENETICALGORITHM_H_INCLUDED
#define GENETICALGORITHM_H_INCLUDED

#include "Chromosome.h"

class GeneticAlgorithm {
private:
	std::vector<Chromosome> population;
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
		} while(resultIndex==firstIndex || resultIndex==secondIndex || (resultIndex == bestIndex && KEEP_ELITE==true) || population[firstIndex].command.length()==0 || population[secondIndex].command.length()==0);
	}

	friend std::ostream& operator<< (std::ostream &out, const GeneticAlgorithm &ga);

public:
	static const bool KEEP_ELITE = true;

public:
	GeneticAlgorithm(int populationSize=0) {
		if(populationSize < 0) {
			populationSize = 0;
		}
		population.resize(populationSize);
		resultIndex = 0;
		firstIndex = 0;
		secondIndex = 0;
		bestIndex = 0;
		worstIndex = 0;
	}

	GeneticAlgorithm(const GeneticAlgorithm &ga) {
		(*this) = ga;
	}

	int getResultIndex() {
		return( resultIndex );
	}

	int getBestIndex() {
		return( bestIndex );
	}

	void setChromosome(Chromosome chromosome, int index=-1) {
		if(index < -1) {
			return;
		}

		if(index == -1) {
			population.push_back( chromosome );
			index = population.size() - 1;
		} else if(index < population.size()) {
			population[index] = chromosome;
		}

		if(population[index].fitness < population[bestIndex].fitness) {
			bestIndex = index;
		}
		if(population[index].fitness > population[worstIndex].fitness) {
			worstIndex = index;
		}
	}

	const Chromosome& getChromosome(int index) const {
		const static std::string EMPTY_STRING = "";

		if(population.size() <= index || index <= -1) {
			//TODO Handle exception.
		}

		return( population[index] );
	}

	const Chromosome& getBestChromosome() const {
		return( population[bestIndex] );
	}

	const Chromosome& getRandomChromosome() const {
		return( population[rand()%population.size()] );
	}

	const Chromosome& getWorstChromosome() const {
		return( population[worstIndex] );
	}

	void replaceWorst(const Chromosome& chromosome) {
		population[worstIndex] = chromosome;

		bestIndex = 0;
		worstIndex = 0;
		for(int i=0; i<population.size(); i++) {
			if(population[i].fitness < population[bestIndex].fitness) {
				bestIndex = i;
			}
			if(population[i].fitness > population[worstIndex].fitness) {
				worstIndex = i;
			}
		}
	}

	void setFitness(double fitness, int index=-1) {
		if(index == -1) {
			index = population.size()-1;
		}

		if(population.size() <= index) {
			//TODO Handle exception.
			return;
		}

		population[index].fitness = fitness;
		if(fitness < population[bestIndex].fitness) {
			bestIndex = index;
		}
		if(fitness > population[worstIndex].fitness) {
			worstIndex = index;
		}
	}

	double getFitness(int index) {
		if(population.size() <= index || index <= -1) {
			//TODO Handle exception.
			return( INVALID_FITNESS_VALUE );
		}

		return( population[index].fitness );
	}

	double getBestFitness() const {
		return( population[bestIndex].fitness );
	}

	int size() {
		return( population.size() );
	}

	void subset(GeneticAlgorithm &ga, const int size) const {
		if(population.size() <= 0) {
			return;
		}

		for(int i=0; i<size; i++) {
			ga.setChromosome( population[ rand()%population.size() ] );
		}
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
			} while (population[resultIndex].fitness < population[firstIndex].fitness
					 || population[resultIndex].fitness < population[secondIndex].fitness);
		} else if (percent
				   < (CROSSOVER_RESULT_INTO_WORST_PERCENT
					  + CROSSOVER_RESULT_INTO_MIDDLE_PERCENT)) {
			do {
				selectRandom();
			} while (population[resultIndex].fitness < population[firstIndex].fitness
					 || population[resultIndex].fitness > population[secondIndex].fitness);
		} else if (percent
				   < (CROSSOVER_RESULT_INTO_WORST_PERCENT
					  + CROSSOVER_RESULT_INTO_MIDDLE_PERCENT
					  + CROSSOVER_RESULT_INTO_BEST_PERCENT)) {
			do {
				selectRandom();
			} while (population[resultIndex].fitness > population[firstIndex].fitness
					 || population[resultIndex].fitness > population[secondIndex].fitness);
		}
	}

	void crossover() {
		population[resultIndex].command = population[firstIndex].command.substr(0, rand()%(population[firstIndex].command.length())+1);
		population[resultIndex].command += population[secondIndex].command.substr(rand()%population[secondIndex].command.length(), population[secondIndex].command.length());
		population[resultIndex].fitness = INVALID_FITNESS_VALUE;
	}

	void mutation() {
		int index = rand() % population[resultIndex].command.length();

		switch(rand()%6) {
		case 0:
			population[resultIndex].command[index]=(char)TOP;
			break;
		case 1:
			population[resultIndex].command[index]=(char)LEFT;
			break;
		case 2:
			population[resultIndex].command[index]=(char)RIGHT;
			break;
		case 3:
			population[resultIndex].command[index]=(char)FRONT;
			break;
		case 4:
			population[resultIndex].command[index]=(char)BACK;
			break;
		case 5:
			population[resultIndex].command[index]=(char)DOWN;
			break;
		}

		population[resultIndex].fitness = INVALID_FITNESS_VALUE;
	}

	void reduction() {
		static const char nop[] = {NONE, '\0'};

		if(COMMANDS_REDUCTION == false) {
			return;
		}

		bool done = true;
		std::string &value = population[resultIndex].command;

		do {
			done = true;

			for(int i=1, j=0; i<value.size(); i++) {
				if(value[i] == value[i-1]) {
					j++;
				} else {
					j = 0;
				}

				if(j == 3) {
					done = false;
					value = value.substr(0,i-3) + nop + value.substr(i+1);

					j = 0;
					i -= 3;
				}
			}
		} while(done == false);
	}

	const std::string& toString() {
		static std::string result;
		result = "";

		/*
		 * Keep population size.
		 */
		result += std::to_string(population.size());
		result += " ";

		for(int i=0; i<population.size(); i++) {
			result += std::to_string(population[i].fitness);
			result += " ";
			if(population[i].command == "") {
				result += NONE;
			} else {
				result += population[i].command;
			}
			result += " ";
		}

		/*
		 * Trim spaces.
		 */
		result.erase(result.size()-1, 1);
		result += '\0';

		return result;
	}

	void fromString(const char text[]) {
		std::string buffer(text);
		std::istringstream in(buffer);

		population.clear();
		bestIndex = 0;
		worstIndex = 0;

		int size = 0;
		in >> size;

		double value;
		std::string commands;
		for(int i=0; i<size; i++) {
			in >> value;
			in >> commands;

			setChromosome(Chromosome(commands,value));

			if(population[bestIndex].fitness > population[i].fitness) {
				bestIndex = i;
			}
			if(population[worstIndex].fitness < population[i].fitness) {
				worstIndex = i;
			}
		}
	}

	void operator=(const GeneticAlgorithm &ga) {
		this->population.clear();

		this->population = ga.population;
		this->resultIndex = ga.resultIndex;
		this->firstIndex = ga.firstIndex;
		this->secondIndex = ga.secondIndex;
		this->bestIndex = ga.bestIndex;
		this->worstIndex = ga.worstIndex;
	}
};

std::ostream& operator<< (std::ostream &out, const GeneticAlgorithm &ga) {
	for(int p=0; p<ga.population.size(); p++) {
		out << ga.population[p].fitness;
		out << "\t";
		for(int i=0; i<ga.population[p].command.length(); i++) {
			out << ga.population[p].command[i];
		}
		out << std::endl;
	}

	return out;
}

#endif
