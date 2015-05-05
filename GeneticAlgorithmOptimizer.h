#ifndef GENETICALGORITHMOPTIMIZER_H_INCLUDED
#define GENETICALGORITHMOPTIMIZER_H_INCLUDED

class GeneticAlgorithmOptimizer {
private:
	static double evaluate(const RubiksCube &solved, const RubiksCube &shuffled, const std::string &commands) {
		static RubiksCube used;

		used = shuffled;
		used.execute(commands);
		return( solved.compare(used) );
	}

	GeneticAlgorithmOptimizer() {
	}

public:
	static void addRandomCommands(GeneticAlgorithm &ga, const RubiksCube &solved, const RubiksCube &shuffled, int populationSize=0) {
		for(int p=0; p<populationSize; p++) {
			RubiksCube mixed;
			std::string commands = mixed.shuffle(CHROMOSOMES_INITIAL_SIZE);
			ga.setChromosome( Chromosome(commands,INVALID_FITNESS_VALUE) );
			ga.setFitness(evaluate(solved, shuffled, commands));
		}
	}

	static void addEmptyCommand(GeneticAlgorithm &ga, const RubiksCube &solved, const RubiksCube &shuffled) {
		ga.setChromosome(Chromosome("",INVALID_FITNESS_VALUE));
		ga.setFitness(evaluate(solved, shuffled, ""));
	}

	static void optimize(GeneticAlgorithm &ga, RubiksCube &solved, RubiksCube &shuffled, long epoches=0) {
		for(long e=0L; e<epoches*ga.size(); e++) {
			ga.selection();
			ga.crossover();
			ga.mutation();
			int index = ga.getResultIndex();
			ga.setFitness(evaluate(solved, shuffled, ga.getChromosome(index).command), index);
		}

		shuffled.execute(ga.getChromosome(ga.getBestIndex()).command);
	}
};

#endif
