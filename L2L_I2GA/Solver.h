#ifndef H_SOLVER
#define H_SOLVER

#include "Specimen.h"
#include <vector>
#include "SVM.h"
// Class that implements the Genetic Algorithm
// at the heighest absraction level
class Solver
{
private:
	SVM sv;
public:
	// Current Population
	static Specimen** Current_Population;

	// Current Selection
	static Specimen** Current_Selection;

	// Current Proximity
	static double Current_Proximity;

	// Current Iteration
	static int Current_Iteration;

	static void Solver::Initialize(vector<MatrixXd> Basis, MatrixXd C, VectorXd D, double E, double EE);
	
	MatrixXd Solver::GAsolve(vector<MatrixXd> Basis, MatrixXd C, VectorXd D, double E, double EE);
};

#endif
