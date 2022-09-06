#ifndef H_SPECIMEN
#define H_SPECIMEN
#include "MatrixElement.h"
#include <vector>
#include "Eigen/Core"
#include "SVM.h"

using namespace Eigen;
using namespace std;
// Class that represents the Specimen
class Specimen
{
private:
	SVM sv;
public:
	// The genes
	double* Genes;

	// Affinity to the solution
	double Affinity;

	// Constructor that creates the Specimen
	// instance with initial genes
	Specimen();

	// Destructor. Frees memory for the Genes
	~Specimen();

	// Clone the Specimen for simple memory management
	Specimen* Clone();

	// Calculates affinity of this instance to the solution
	// Contains the equation formula
	void CalculateAffinity(vector<MatrixXd> Basis, MatrixXd C, VectorXd D, double E, double EE);
};

#endif
