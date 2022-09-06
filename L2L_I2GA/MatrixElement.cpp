#include "MatrixElement.h"
#include <cmath>
#include <iostream>
#include "Eigen/Core"
#include "Eigen/Dense"
#include <vector>
using namespace Eigen;
using namespace std;

//=============================================================================
MatrixElement::MatrixElement()
{
}
//=============================================================================
double MatrixElement::overlap(MatrixXd A1x, MatrixXd A2x)
{
	MatrixXd matrix00 = MatrixXd::Zero(2, 2);
	matrix00(0, 0) = 1.0;
	MatrixXd matrix11 = MatrixXd::Zero(2, 2);
	matrix11(1, 1) = 1.0;
	double det;
	double overlap = 0.0;
	det = (A1x + A2x).determinant();
	double Inv_sqrt_det = 1.0 / sqrt(det);
	if (Inv_sqrt_det == 0.0) { std::cout << "det=0" << endl; overlap = 100000; }
	else {
		overlap = Inv_sqrt_det;
		//overlap = Inv_sqrt_det* Inv_sqrt_det* Inv_sqrt_det;
	}
	return overlap;
}

//==============================================================================================
double MatrixElement::energy(MatrixXd A1x, MatrixXd A2x) {
	
	MatrixXd matrix00= MatrixXd::Zero(2, 2);
	matrix00(0, 0) = 1.0;
	MatrixXd matrix11 = MatrixXd::Zero(2, 2);
	matrix11(1, 1) = 1.0;

	// kinetic energy and harmonic potential energy of coodinate1
	double KinEnergy;
	double harmonicEnergy1;
	MatrixXd InvA = (A1x + A2x).inverse();
	double det = (A1x + A2x).determinant();
	MatrixXd identity_matrix = MatrixXd::Zero(2, 2);
	identity_matrix(0, 0) = 0.5;
	identity_matrix(1, 1) = 2.0;
	// Kinetic Energy of coodinate1
	if (det == 0.0)
		KinEnergy = 0.0;
	else {
		double Inv_sqrt_det = 1 / sqrt(det);
		MatrixXd T_kin1 = A2x * InvA * A1x*identity_matrix;
		KinEnergy= 0.5*T_kin1.trace() * overlap(A1x, A2x);
		//KinEnergy = 1.5 * T_kin1.trace() * overlap(A1x, A2x);
	}

	// harmonic Energy of
	double harmonicEnergy;

	if (det == 0.0)
		harmonicEnergy = 10000000.0;
	else {
		double Inv_sqrt_det = 1 / sqrt(det);
		MatrixXd T_harm1 =  InvA  * identity_matrix;
		harmonicEnergy = 0.5*T_harm1.trace() * overlap(A1x, A2x);
		//harmonicEnergy = 1.5 * T_harm1.trace() * overlap(A1x, A2x);
	}

	// 2 body energy (1 - 2)
	double r0 = 0.01;
	double asc = 2.0;
	// 2 body energy (1 - 3)
	double v13 = -2.0 / asc;
	VectorXd w13(2);
	w13(0) = 1.0;
	w13(1) = -1.0;
	double c13_inv = w13.transpose() * (A1x + A2x).inverse() * w13;
	//double h13 = sqrt(2 * 3.141592 * (r0 * r0 + c13_inv));
	double h13= sqrt(2 * 3.141592 * (r0 * r0 + 1/(A1x(1,1)+A2x(1,1))));
	double potEnergy_13 = v13/h13 * overlap(A1x, A2x);

	//double potEnergy_13 = v13 / (sqrt(h13)* sqrt(h13)* sqrt(h13)) * overlap(A1x, A2x);
	//return KinEnergy  +harmonicEnergy+ potEnergy_13;
	return KinEnergy + potEnergy_13;
}