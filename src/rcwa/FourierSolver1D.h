#pragma once
#include "defns.h"
#include <Eigen/Core>

class FourierSolver1D
{
private:
    scalar d_;
    Eigen::VectorXs x_;

public:
    FourierSolver1D(const Eigen::VectorXs& x, scalar d): 
        d_(d), x_(x) { }

    void solve(const Eigen::VectorXcs& f,
        int M,
        Eigen::VectorXcs& Ff);

    void dSolve(const Eigen::VectorXcs& f,
        int M,
        int layout,
        Eigen::VectorXcs& dFf);
};
