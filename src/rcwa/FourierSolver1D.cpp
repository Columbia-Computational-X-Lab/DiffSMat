#include "FourierSolver1D.h"

#include <iostream>
#include "MathFunction.h"

using namespace std;
using namespace Eigen;

void FourierSolver1D::solve(const Eigen::VectorXcs& f,
    int M,
    Eigen::VectorXcs& Ff)
{
    if (f.size() != x_.size() - 1)
    {
        std::cerr << "The function configuration does not match the coordinates!\n";
    }

    Ff.resize(M);
    Ff.setZero();

    for (int i = 0; i < M; ++i)
    {
        int m = i - (M-1) / 2;

        for (int j = 1; j < x_.size(); ++j)
        {
            scalar x1 = x_(j);
            scalar x0 = x_(j-1);
            if (m == 0)
            {
                Ff(i) += f(j-1) / d_ * (x1 - x0);
            }
            else
            {
                Ff(i) += f(j-1) / d_ * exp(-scalex(0, m*Pi*(x1 + x0)/d_))
                * sinc(m * (x1 - x0) / d_) * (x1 - x0);     
            }
        }
    }
}

void FourierSolver1D::dSolve(const Eigen::VectorXcs& f,
    int M,
    int layout,
    Eigen::VectorXcs& dFf)
{
    if (f.size() != x_.size() - 1)
    {
        cerr << "[ERROR]: The function configuration does not match the coordinates!\n";
    }

    if (layout < 0 || layout > x_.size() - 1)
    {
        cerr << "[ERROR]: invalid layout!\n";
    }

    dFf.resize(M);
    dFf.setZero();

    for (int i = 0; i < M; ++i)
    {
        int m = i - (M-1) / 2;

        if (m == 0)
        {
            if (layout > 0)
            {
                dFf(i) += f(layout-1) / d_;
                
            }
            if (layout < x_.size() - 1)
            {
                dFf(i) += -f(layout) / d_;
            }
        }
        else
        {
            if (layout > 0)
            {
                dFf(i) += f(layout-1) / d_ 
                                * exp(-scalex(0, 2*m*Pi*x_(layout)/d_));
            }

            if (layout < x_.size() - 1)
            {
                dFf(i) += -f(layout) / d_ 
                                * exp(-scalex(0, 2*m*Pi*x_(layout)/d_));                
            }

        }
    }   
}