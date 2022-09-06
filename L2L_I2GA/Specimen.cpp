#include "Specimen.h"
#include <math.h>
#include <cmath>
#include <time.h>
#include "Rand.h"
#include "SVM.h"
#include <vector>
#include "SVM.h"

// Constructor that creates the Specimen
// instance with initial genes

Specimen::Specimen():sv()
{
	// Assume one double value as the genes
	Genes = new double[3];

	// Set up initial value to the genes
	Genes[0] = 2;
	Genes[1] = 10;
	Genes[2] = 5;
}

// Destructor. Frees memory for the Genes
Specimen::~Specimen()
{
	delete[] Genes;
}

// Clone the Specimen for simple memory management
Specimen* Specimen::Clone()
{
	Specimen* s = new Specimen();
	s->Genes[0] = Genes[0];
	s->Genes[1] = Genes[1]=0;
	s->Genes[2] = Genes[2];
	s->Affinity = Affinity;

	return s;
}

// Calculates affinity of this instance to the solution
// Contains the equation formula


void Specimen::CalculateAffinity(vector<MatrixXd> Basis, MatrixXd C, VectorXd D, double E, double EE)
{
	
	MatrixXd Genesmatrix = MatrixXd::Zero(2, 2);
	Genesmatrix(0, 0) = Genes[0];
	Genesmatrix(0, 1) = Genes[1];
	Genesmatrix(1, 1) = Genes[2];
	Genesmatrix(1, 0) = Genes[1];

	int Bsize = Basis.size();
	Basis.push_back(Genesmatrix);
	
	double NewE;
	if (sv.CheckOverlap(Basis) == 1) {
		NewE = sv.NewEnergy(Basis, C, D, E, EE);
	}
	else { NewE=10000000; }
	Basis.resize(Bsize);
	
	Affinity = NewE;
	
}
