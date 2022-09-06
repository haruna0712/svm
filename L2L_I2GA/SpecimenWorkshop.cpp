#include "SpecimenWorkshop.h"
#include <stdlib.h>
#include <vector>

// Number of specimens in the population
int SpecimenWorkshop::PopulationSize;

// Number of specimens in the selection 
int SpecimenWorkshop::SelectionSize;

// Power of mutations
// 0 - no mutations
// 100 - mutations up to whole gene value
// 200 - mutations up to twice gene value
double SpecimenWorkshop::MaxMutationPercent;

// Likelyhood of mutation
// 0 - no mutations, 100 - mutations for all cases
int SpecimenWorkshop::MutationLikelyhoodPercent;

// Maximal affinity that is considered 
// for the solution found
double SpecimenWorkshop::Epsilon;

// Generate initial population
Specimen** SpecimenWorkshop::GeneratePopulation(vector<MatrixXd> Basis, MatrixXd C, VectorXd D, double E, double EE)
{
	// Creates array representing the population
	Specimen** p = new Specimen*[PopulationSize];

	// Creates specimens
	// Mutation of all specimens in initial population
	// increases variance that increases chance to
	// get better instance.
	for (int i = 0; i < PopulationSize; i++)
	{
		p[i] = new Specimen();
		Mutate(p[i]);

		// Calculate Affinity for new specimens
		p[i]->CalculateAffinity(Basis, C, D, E, EE);
	}

	return p;
}

// Generate population by reproduction of selection
Specimen** SpecimenWorkshop::GeneratePopulation(Specimen** selection, vector<MatrixXd> Basis, MatrixXd C, VectorXd D, double E, double EE)
{
	// Creates array representing the population
	Specimen** p = new Specimen*[PopulationSize];

	// Copy instances from the selection to keep them
	// in new generation of population
	for (int i = 0; i < SelectionSize; i++)
	{
		p[i] = selection[i]->Clone();
	}

	// Creates new specimens by reproducing two parents
	// Parents are selected randomly from the selection.
	int chield_index = SelectionSize;
	int parent1_index;
	int parent2_index;

	while (chield_index < PopulationSize)
	{
		// Slect two parents randomly in way
		// they are different instances
		do
		{
			parent1_index = rand() % SelectionSize;
			parent2_index = rand() % SelectionSize;
		} while (parent1_index == parent2_index);

		// Creates new specimen
		p[chield_index] = ReproduceNew(selection[parent1_index], selection[parent2_index], Basis, C, D, E, EE);

		chield_index++;
	}

	return p;
}

// Reproduce new specimen on base of two parents
Specimen* SpecimenWorkshop::ReproduceNew(Specimen* a, Specimen* b, vector<MatrixXd> Basis, MatrixXd C, VectorXd D, double E, double EE)
{
	Specimen* s = new Specimen();

	// Iherit genes as the average oh the parents' genes
	s->Genes[0] = (a->Genes[0] + b->Genes[0]) / 2;
	s->Genes[1] = 0*(a->Genes[1] + b->Genes[1]) / 2;
	s->Genes[2] = (a->Genes[2] + b->Genes[2]) / 2;
	// Mutate if likelyhoo allows
	int ml = rand() % 101;
	if (ml <= MutationLikelyhoodPercent)
	{
		Mutate(s);
	}

	// Calculate Affinity for new specimen
	s->CalculateAffinity(Basis, C, D, E, EE);

	return s;
}

// Select the best specimens from the population
Specimen** SpecimenWorkshop::Select(Specimen** population)
{
	// Sort population by increasing the affinity
	// The best specimens are moving to start of the array
	Sort(population);

	// Create set of selected specimens
	Specimen** selected = new Specimen*[SelectionSize];

	// Copy best specimens into the selection
	for (int i = 0; i < SelectionSize; i++)
	{
		selected[i] = population[i]->Clone();
	}

	return selected;
}

// Sort the population
void SpecimenWorkshop::Sort(Specimen** population)
{
	Specimen* temp;

	for (int i = 0; i < PopulationSize; i++)
	{
		for (int j = 0; j < PopulationSize; j++)
		{
			if (population[i]->Affinity < population[j]->Affinity)
			{
				temp = population[i];
				population[i] = population[j];
				population[j] = temp;
			}
		}
	}
}

// Mutate the specimen
void SpecimenWorkshop::Mutate(Specimen* sp)
{
	// Calculate Mutation Factor that is random value between 0 and MaxMutationPercent
	// calculates as ratio between 0 and 1
	double MutationFactor0 = (MaxMutationPercent / 100.0) * (rand() % 5001 / 5000.0);

	/*
	// Set mutation to negative with 50% likelyhood
	if ((rand() % 10) < 5)
	{
		MutationFactor0 = -MutationFactor0;
	}
	*/
	double MutationFactor1 = (MaxMutationPercent / 100.0) * (rand() % 5001 / 5000.0);
	
	/*
	// Set mutation to negative with 50% likelyhood
	if ((rand() % 10) < 5)
	{
		MutationFactor1 = -MutationFactor1;
	}
	*/
	double MutationFactor2 = (MaxMutationPercent / 100.0) * (rand() % 5001 / 5000.0);

	/*
	// Set mutation to negative with 50% likelyhood
	if ((rand() % 10) < 5)
	{
		MutationFactor2 = -MutationFactor2;
	}
	*/

	// Calculate new gene
	sp->Genes[0] = sp->Genes[0] * MutationFactor0;
	sp->Genes[1] = sp->Genes[1] * (MutationFactor1)*0;
	sp->Genes[2] = sp->Genes[2] * (MutationFactor2);
}
