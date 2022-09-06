#pragma once
#define SVM_H
#include "MatrixElement.h"
#include "Eigen/Core"
#include <vector>
using namespace Eigen;
using namespace std;

class SVM
{
private:
	MatrixElement me;
	int kk0;
	double bmin, bmax;

public:
	SVM();
	int CheckOverlap(vector<MatrixXd> &Basis);
	MatrixXd FirstNewState();
	double NewEnergy(vector<MatrixXd> &Basis, MatrixXd C, VectorXd D, double E, double EE);
	MatrixXd NewState(vector<MatrixXd>& Basis, MatrixXd C, VectorXd D, double E, double EE);
	MatrixXd Dmatrix();
	MatrixXd NormMatrix(vector<MatrixXd> &Basis);
	MatrixXd HamiltonianMatrix(vector<MatrixXd> &Basis);
	MatrixXd Hmatrix, Nmatrix;
	void UpdateNorm(vector<MatrixXd> &Basis);
	void UpdateHamiltonian(vector<MatrixXd> &Basis);
};