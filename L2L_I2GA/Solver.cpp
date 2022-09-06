#include "Solver.h"
#include "SpecimenWorkshop.h"
#include <iostream>
#include <vector>
#include "Specimen.h"
// Current Population
Specimen** Solver::Current_Population;

// Current Selection
Specimen** Solver::Current_Selection;

// Current Proximity
double Solver::Current_Proximity;

// Current Iteration
int Solver::Current_Iteration;

// Initialize the algorithm
void Solver::Initialize(vector<MatrixXd> Basis, MatrixXd C, VectorXd D, double E, double EE)
{
	// Set up options for the algorithm
	SpecimenWorkshop::PopulationSize = 1000;
	SpecimenWorkshop::SelectionSize = 70;
	SpecimenWorkshop::MaxMutationPercent = 500;
	SpecimenWorkshop::MutationLikelyhoodPercent = 20;
	SpecimenWorkshop::Epsilon = 0.01;

	// Generate initial population
	Current_Population = SpecimenWorkshop::GeneratePopulation(Basis, C, D, E, EE);

	// Set Current_Selection to zero for correct work delete[] operator
	Current_Selection = 0;

	Current_Iteration = 0;
}

// Run the algorithm
MatrixXd Solver::GAsolve(vector<MatrixXd> Basis, MatrixXd C, VectorXd D, double E, double EE)
{

	MatrixXd newstate_matrix = MatrixXd::Zero(2, 2);
	int Bsize = Basis.size();
	int count = 0;
	// Set Current_Proximity to the biggest value
	Current_Proximity = 1.7976931348623157E+308;

	// Loop while Current_Proximity is not less than the Epsilon
	while (1)
	{
		double Previous_Proximity = Current_Proximity;
		// Delete old selection
		if (Current_Selection != 0)
		{
			for (int i = 0; i < SpecimenWorkshop::SelectionSize; i++)
			{
				// Calculate Affinity for new specimens
				delete Current_Selection[i];
			}

			delete[] Current_Selection;
		}

		// Select the best specimens
		Current_Selection = SpecimenWorkshop::Select(Current_Population);

		// Calculate proximity for the top-selected (the best) specimen
		Current_Proximity = Current_Selection[0]->Affinity;
		Current_Selection[0]->Genes[0];

		// End the calculations if Current_Proximity is less than the Epsilon
		
		
		double diff = Previous_Proximity - Current_Proximity;
		if (diff>=0 && diff< SpecimenWorkshop::Epsilon)
		{
			count++;
		}
		
		if (count == 3) {
			newstate_matrix(0, 0) = Current_Selection[0]->Genes[0];
			newstate_matrix(0, 1) = Current_Selection[0]->Genes[1];
			newstate_matrix(1, 0) = Current_Selection[0]->Genes[1];
			newstate_matrix(1, 1) = Current_Selection[0]->Genes[2];
				return newstate_matrix;
			break;
		}
		
		// Delete old population
		for (int i = 0; i < SpecimenWorkshop::PopulationSize; i++)
		{
			// Calculate Affinity for new specimens
			delete Current_Population[i];
		}

		delete[] Current_Population;

		// Generate new population by reproducing specimens from the selection
		Current_Population = SpecimenWorkshop::GeneratePopulation(Current_Selection, Basis, C, D, E, EE);

		Current_Iteration++;
		
	}
}
