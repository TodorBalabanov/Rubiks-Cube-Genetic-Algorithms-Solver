#ifndef CHROMOSOME_H_INCLUDED
#define CHROMOSOME_H_INCLUDED

class Chromosome {
public:
	double fitness;
	std::string command;

	Chromosome(std::string command, double fitness) {
		this->command = command;
		this->fitness = fitness;
	}

	Chromosome(const Chromosome &chromosome) {
		(*this) = chromosome;
	}

	Chromosome() {
		this->command = "";
		this->fitness = INVALID_FITNESS_VALUE;
	}

	void operator=(const Chromosome &chromosome) {
		this->command = chromosome.command;
		this->fitness = chromosome.fitness;
	}
};

#endif
