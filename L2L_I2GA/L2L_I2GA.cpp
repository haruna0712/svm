
#include "SVM.h"
#include "MatrixElement.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

#include "Eigen/Core"
#include "Eigen/Eigen"
#include <iomanip>
#include <string>
#include <stdlib.h> 
#include <ctime> 
#include <cstdio>
#include <cstdlib>
#include <stdio.h>
#include <time.h>
#include "Solver.h"
using namespace Eigen;
using namespace std;

int main()
{
    
    int maxbasis = 20;
    SVM svm;
    Solver solver;
    MatrixXd NewState;
    GeneralizedSelfAdjointEigenSolver<MatrixXd> ges;
    vector<MatrixXd> Basis;
    

    NewState = svm.FirstNewState();
    cout << "first new state is  " << endl;
    cout << NewState << "\n";

    /* start SVM iterations */
    Basis.push_back(NewState);
    svm.UpdateNorm(Basis);
    svm.UpdateHamiltonian(Basis);

    MatrixXd Norm;
    MatrixXd H;
    MatrixXd C;
    VectorXd D;
    VectorXd D_sort;
    double E;
    double EE;
    vector<double> dE;

    cout << "\t Start SVM iters\n\n";
    int itr = 1;
    double n_accuracy = 1;
    clock_t start1 = clock();
    //-------------------------------------------------
    while (itr < maxbasis)
    {

        if (itr == 8) { clock_t end1 = clock(); 
        std::cout << "duration = " << (double)(end1 - start1) / CLOCKS_PER_SEC << "sec.\n";
        }
        

        Norm = svm.NormMatrix(Basis);
        H = svm.HamiltonianMatrix(Basis);
        ges.compute(H, Norm);
        C = ges.eigenvectors();
        D = ges.eigenvalues();
        sort(D.begin(), D.end());

        if (itr < 2) E = D[0];
        else
            E = D[0];
        D = ges.eigenvalues();
        if (itr == 1)  EE = E + abs(E / 2);
        printf("\t iter = %4d     E = %14.8f  \n", itr, E);
        
        Solver::Initialize(Basis, C, D, E, EE);
        NewState = solver.GAsolve(Basis, C, D, E, EE);
       
        //NewState = svm.NewState(Basis,C,D,E,EE);
        cout << "newstate is"<<endl << NewState<<endl;
        Basis.push_back(NewState);


        dE.push_back(abs((EE - E) / E));
        if (dE[itr - 1] > 0.0005) n_accuracy = 1;
        if (dE[itr - 1] < 0.0005) n_accuracy++;
        if (n_accuracy == 6) break;
        EE = E;

        svm.UpdateNorm(Basis);
        svm.UpdateHamiltonian(Basis);
        itr = itr + 1;
        

    }
    return 0;
}