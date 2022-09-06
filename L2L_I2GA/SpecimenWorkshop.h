

#ifndef H_SPECIMENWORKSHOP
#define H_SPECIMENWORKSHOP

#include "Specimen.h"
#include <vector>
// Class that performs main actions with the specimens and the population
class SpecimenWorkshop
{
public:

	// Number of specimens in the population
	static int PopulationSize;

	// Number of specimens in the selection 
	static int SelectionSize;

	// Power of mutations
	// 0 - no mutations
	// 100 - mutations up to whole gene value
	// 200 - mutations up to twice gene value
	static double MaxMutationPercent;

	// Likelyhood of mutation
	// 0 - no mutations, 100 - mutations for all cases
	static int MutationLikelyhoodPercent;

	// Maximal affinity that is considered 
	// for the solution found
	static double Epsilon;

	// Generate initial population
	static Specimen** GeneratePopulation(vector<MatrixXd> Basis, MatrixXd C, VectorXd D, double E, double EE);

	// Generate population by reproduction of selection
	static Specimen** GeneratePopulation(Specimen** selection, vector<MatrixXd> Basis, MatrixXd C, VectorXd D, double E, double EE);

	// Reproduce new specimen on base of two parents
	static Specimen* ReproduceNew(Specimen* a, Specimen* b, vector<MatrixXd> Basis, MatrixXd C, VectorXd D, double E, double EE);

	// Select the best specimens from the population
	static Specimen** SpecimenWorkshop::Select(Specimen** population);

	// Sort the population
	static void SpecimenWorkshop::Sort(Specimen** population);

	// Mutate the specimen
	static void SpecimenWorkshop::Mutate(Specimen* sp);
};

#endif
