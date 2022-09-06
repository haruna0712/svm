#pragma once
#define MatrixElement_H
#include <vector>
#include "Eigen/Core"
using namespace Eigen;
//=============================================================================

class MatrixElement
{

public:
    MatrixElement();
    double overlap(MatrixXd state1, MatrixXd state2);
    double energy(MatrixXd state1, MatrixXd state2);
};