#include "SVM.h"
#include "MatrixElement.h"
#include <iostream>
#include <vector>
#include <cmath>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <time.h>
#include"XorShift.h"
using namespace Eigen;
using namespace XorShifts;
GeneralizedSelfAdjointEigenSolver<MatrixXd> ges;
VectorXd D;
double E;
//=============================================================================
SVM::SVM() : me()
{
	bmin = 0;
	bmax = 20;
	kk0 = 10000;
	int ndb = 100;
	Hmatrix = MatrixXd::Zero(ndb, ndb);
	Nmatrix = MatrixXd::Zero(ndb, ndb);
}
//=============================================================================
int SVM::CheckOverlap(vector <MatrixXd> &Basis)
{
	int itr = Basis.size() - 1;
	double vnorm, vdotv;
	MatrixXd normmatrix = NormMatrix(Basis);
	vnorm = normmatrix(itr, itr);
	if (vnorm < 1e-8) return 0;
	if (itr > 0) {
		for (int i = 0; i < itr; i++) {
			vdotv = normmatrix(itr, i) / sqrt(normmatrix(i, i) * normmatrix(itr, itr));
			if (vdotv > 0.995) {
				return 0;
			}
		}
	}
	return 1;
}
//=============================================================================
double EigenValuesEquation(int itr, VectorXd D, VectorXd q, double aa, double xx)
{
	double vv = 1;
	double ww = 1;
	double yy = 0;
	double zz = 0;

	for (int n1 = 0; n1 < itr; n1++)
	{
		vv = vv * (D(n1) - xx);
		ww = 1;
		for (int n2 = 0; n2 < itr; n2++)
		{
			if (n2 != n1) ww = ww * (D(n2) - xx);
		}
		yy = yy + q(n1) * q(n1) * ww;

	}
	zz = (aa - xx) * vv - yy;

	return zz;
}
//=============================================================================


/*Basis includes new state*/
double SVM::NewEnergy(vector <MatrixXd> &Basis, MatrixXd C, VectorXd D, double E, double EE)
{

	int itr = Basis.size() - 1;
	vector<double> voverlap(itr);
	for (int k1 = 0; k1 < itr; k1++)
	{
		voverlap[k1] = me.overlap(Basis[itr], Basis[k1]);
	}
	vector<double> eoverlap(itr);
	for (int k1 = 0; k1 < itr; k1++)
	{
		eoverlap[k1] = me.energy(Basis[itr], Basis[k1]);
	}
	//cout << "itr is " << itr;
	VectorXd c = VectorXd::Zero(itr + 1);
	VectorXd Overlap = VectorXd::Zero(itr);
	VectorXd q = VectorXd::Zero(itr);
	double aa = 0;
	double NN = 0;

	//solving the equation for the new eigenvalue===============================================  
	int count = 0;
	double e1 = E;
	double e2 = E - abs(0.5 * (E - EE));
	double e3 = E;
	
	NN = 0;
	VectorXd Overlap2 = VectorXd::Zero(itr);
	q = VectorXd::Zero(itr);
	c = VectorXd::Zero(itr + 1);
	for (int k1 = 0; k1 < itr; k1++)
	{
		for (int k2 = 0; k2 < itr; k2++)
		{
			//cout <<"Nmatrix(k1,k2) is " << Nmatrix(itr, k2) << "\t"<<"me.overlap is" << me.overlap(Basis[itr], Basis[k2]) << endl;
			Overlap2(k1) = Overlap2(k1) + C(k2, k1) * me.overlap(Basis[itr],Basis[k2]);
		}
	}

	c(itr) = 1;
	for (int k1 = 0; k1 < itr; k1++)
	{
		for (int k2 = 0; k2 < itr; k2++)
		{
			c(k1) = c(k1) - Overlap2(k2) * C(k1, k2);
		}
	}
	for (int k1 = 0; k1 < itr+1; k1++)
	{
		for (int k2 = 0; k2 < itr+1; k2++)
		{

			NN = NN + c(k1) * c(k2) * me.overlap(Basis[k1],Basis[k2]);
		}
	}
	//cout << "NN calculated by overlap" << NN<<endl;


	NN = 0;
	Overlap2 = VectorXd::Zero(itr);
	q = VectorXd::Zero(itr);
	c = VectorXd::Zero(itr + 1);
	for (int k1 = 0; k1 < itr; k1++)
	{
		for (int k2 = 0; k2 < itr; k2++)
		{
			//cout <<"Nmatrix(k1,k2) is " << Nmatrix(itr, k2) << "\t"<<"me.overlap is" << me.overlap(Basis[itr], Basis[k2]) << endl;
			Overlap2(k1) = Overlap2(k1) + C(k2, k1) * voverlap[k2];
		}
	}

	c(itr) = 1;
	//cout << "c" << c << endl;
	for (int k1 = 0; k1 < itr; k1++)
	{
		for (int k2 = 0; k2 < itr; k2++)
		{
			c(k1) = c(k1) - Overlap2(k2) * C(k1, k2);
		}
	}
	for (int k1 = 0; k1 < itr; k1++)
	{
		for (int k2 = 0; k2 < itr; k2++)
		{

			NN = NN + c(k1) * c(k2) * Nmatrix(k1, k2);
		}
	}

	for (int k1 = 0; k1 < itr; k1++)
	{
		NN = NN + 2 * c(k1) * c(itr) * voverlap[k1];
	}
	NN = NN + me.overlap(Basis[itr], Basis[itr]);
	//cout << "NN calculated by Nmatrix" << NN<<endl;
	//voverlap.erase(voverlap.begin(), voverlap.end());
	
	for (int k1 = 0; k1 < itr + 1; k1++)
	{
		c(k1) = c(k1) / sqrt(NN);
	}
	//cout<<"c" << c << endl;
	//cout << "NN" << NN << endl;
	/*
	for (int k1 = 0; k1 < itr; k1++)
	{
		for (int k2 = 0; k2 < itr; k2++)
		{
			for (int k3 = 0; k3 < itr + 1; k3++)
			{
				q(k1) = q(k1) + Hmatrix(k2, k3) * c(k3) * C(k2, k1);
			}
		}
	}



	for (int k1 = 0; k1 < itr; k1++)
	{
		for (int k2 = 0; k2 < itr; k2++)
		{
			aa = aa + c(k1) * c(k2) * Hmatrix(k1, k2);
		}
	}
	for (int k1 = 0; k1 < itr; k1++)
	{
		aa = aa + 2 * c(k1) * c(itr) * eoverlap[k1];
	}
	aa = aa + me.energy(Basis[itr], Basis[itr]);
	cout << "aa calculated by Hmatrix" << aa << endl;
	*/
	aa = 0;
	q = VectorXd::Zero(itr);
	for (int k1 = 0; k1 < itr; k1++)
	{
		for (int k2 = 0; k2 < itr; k2++)
		{
			for (int k3 = 0; k3 < itr + 1; k3++)
			{
				q(k1) = q(k1) + me.energy(Basis[k2],Basis[k3]) * c(k3) * C(k2, k1);
				//cout << "c(k3)" << c(k3) << endl;
				//cout << "me.energy" << me.energy(Basis[k2], Basis[k3])<<endl;
			}
		}
	}
	//cout << q;


	for (int k1 = 0; k1 < itr+1; k1++)
	{
		for (int k2 = 0; k2 < itr+1; k2++)
		{
			aa = aa + c(k1) * c(k2) * me.energy(Basis[k1],Basis[k2]);
		}
	}
	
	//cout << "aa calculated by me.energy only" << aa << endl;
	//eoverlap.erase(eoverlap.begin(), eoverlap.end());
	
	double Ee1 = EigenValuesEquation(itr, D, q, aa, e1);
	while (count < 101)
	{

		if (Ee1 * EigenValuesEquation(itr, D, q, aa, e2) < 0)  break;
		else
		{
			e1 = e2;
			e2 = e2 - abs(0.5 * (E - EE));
			count++;
		}

	}
	if (count <= 100)
	{
		count = 0;
		while (abs((e1 - e2) / e2) > abs(1e-3 * (E - EE) / EE))
		{
			e3 = (e1 + e2) / 2;
			if (std::signbit(EigenValuesEquation(itr, D, q, aa, e3)) != std::signbit(EigenValuesEquation(itr, D, q, aa, e2)))  e1 = e3;
			else e2 = e3;
			count++;
			if (count > 100)  break;
		}

	}

	return e3;
	//==============================================================================
}


//=============================================================================
MatrixXd SVM::Dmatrix()
{
	MatrixXd d = MatrixXd::Zero(2, 2);
	auto r_const = XorShift(10);
	auto r = XorShift();

	d(0, 0) = bmax * r.randDouble();
	d(1, 1) = bmax * r.randDouble();
	d(0, 1) = 0;
	d(1, 0) = d(0, 1);
	/*condition that d matrix is non-negative quadratic form*/


	return d;
}


//=============================================================================
MatrixXd SVM::FirstNewState() {
	MatrixXd NewState = MatrixXd::Zero(2, 2);

	double e_overlap = 0;
	int count = 0;
	while (count < 10) {
		count++;
		NewState = Dmatrix();
		//NewState.print();
		e_overlap = me.overlap(NewState, NewState);
		/*
		cout<<e_overlap;
		cout << "success to calsulate e_overlap";
		*/
		if (e_overlap > 1e-8) break;

	}

	//if (count >= 8) { NewState.notdefined = true; return NewState; }

	double MinE = me.energy(NewState, NewState) / e_overlap;
	//cout << MinE;
	double NewE;
	MatrixXd State;
	State = NewState;
	count = 0;
	while (count < 10000) {
		count++;
		MatrixXd matrixx = Dmatrix();
		NewState = matrixx;
		e_overlap = me.overlap(NewState, NewState);
		if (e_overlap < 1e-8) continue;
		NewE = me.energy(NewState, NewState) / e_overlap;
		if (NewE < MinE) { MinE = NewE; State = NewState; }
	}

	cout << "first new state is  " << endl;
	IOFormat OctaveFmt(StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
	cout << "Ax=\n" << State.format(OctaveFmt) << "\n";

	return State;
}
//=============================================================================

//=============================================================================
MatrixXd SVM::NormMatrix(vector<MatrixXd> &Basis)
{
	int itr = Basis.size();
	MatrixXd Norm = MatrixXd::Zero(itr, itr);
	for (int n1 = 0; n1 < itr; n1++)
	{
		for (int n2 = n1; n2 < itr; n2++)
		{

			Norm(n1, n2) = me.overlap(Basis[n1], Basis[n2]);
			Norm(n2, n1) = Norm(n1, n2);
		}
	}
	return Norm;
}
MatrixXd SVM::NewState(vector<MatrixXd>& Basis, MatrixXd C, VectorXd D, double E, double EE)
{

	int Bsize = Basis.size();


	MatrixXd dx = Dmatrix();
	MatrixXd mindx = dx;
	Basis.push_back(dx);

	int count1;
	int xx = 0;
	double minE, NewE;

	count1 = 0;
	minE = E; NewE = E;
	while (count1 < kk0) {
		Basis[Bsize] = Dmatrix();
		

		if (CheckOverlap(Basis) == 1) {
			NewE = NewEnergy(Basis, C, D, E, EE);
			//cout <<"newstate" << NewE << endl;
			if (NewE < minE) {
				minE = NewE;
				mindx = Basis[Bsize];
			}
		}
		count1++;

	}

	if (isfinite(minE)) {
		std::cout << "SVM.cpp New state minE is" << minE << "mindx is" << mindx << endl;
	}
	Basis.resize(Bsize);
	return mindx;
}
//=============================================================================
MatrixXd SVM::HamiltonianMatrix(vector<MatrixXd> &Basis)
{
	int it = Basis.size();
	MatrixXd H = MatrixXd::Zero(it, it);
	for (int n1 = 0; n1 < it; n1++)
	{
		for (int n2 = n1; n2 < it; n2++)
		{
			H(n1, n2) = me.energy(Basis[n1], Basis[n2]);
			H(n2, n1) = H(n1, n2);
		}
	}
	return H;
}
//=============================================================================
void SVM::UpdateNorm(vector<MatrixXd> &Basis)
{
	//  cout << Basis.size() << "\n";
	int ncur = Basis.size() - 1;
	for (int n1 = 0; n1 <= ncur; n1++)
	{
		//      cout << Basis.size() << " " << n1 << "\n";
		Nmatrix(n1, ncur) = me.overlap(Basis[n1], Basis[ncur]);
		Nmatrix(ncur, n1) = Nmatrix(n1, ncur);
	}
}
//=============================================================================
void SVM::UpdateHamiltonian(vector<MatrixXd> &Basis)
{
	int ncur = Basis.size() - 1;
	for (int n1 = 0; n1 <= ncur; n1++)
	{
		Hmatrix(n1, ncur) = me.energy(Basis[n1], Basis[ncur]);
		Hmatrix(ncur, n1) = Hmatrix(n1, ncur);
	}
}
