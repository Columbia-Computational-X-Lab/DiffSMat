// This file is part of DiffSMat.
// 
// DiffSMat is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// DiffSMat is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with DiffSMat. If not, see <https://www.gnu.org/licenses/>.

#include <iostream>
#include "WaveEquationSolver.h"
#include "WaveEquationCoeff.h"
#include "MKLEigenSolver.h"

using namespace std;
using namespace Eigen;

struct GammaSqrtOp
{
    GammaSqrtOp() {} 

    inline const scalex operator() (const scalex& g) const
    {
        scalex gv = sqrt(g);
        if ( gv.imag() < 0 ) gv *= -1.;
        // now gv is at the upper half of the complex plane

        if (abs(gv.imag()) < abs(gv.real()) && gv.real() < 0)
            gv *= -1.;
        return gv;
    }
};

void WaveEquationSolver::solve(
    scalar lambda, scalar wx, scalar wy,
    Eigen::MatrixXcs& eigvecE,
    Eigen::MatrixXcs& eigvecH,
    Eigen::VectorXcs& keff) const
{
    WaveEquationCoeff coeff(nx_, ny_, Lx_, Ly_);
    MatrixXcs P, Q;
    coeff.solve(lambda, wx, wy, P, Q);
    MKLEigenSolver ces;
    ces.compute(P*Q);

#if 0
    VectorXcs gamma_ = ces.eigenvalues() * invk02;
    for (int i = 0; i < gamma_.size(); ++i)
    {
        scalex gamma_i = sqrt(gamma_(i));

        // energy gain is impossible for dielectric and metal
        // sometimes we could get a negative imag part of gamma,
        // this would cause NaN in the final results
        if (gamma_i.imag() < 0)
        {
            gamma_i = -gamma_i;
        }
        
        if (abs(gamma_i.imag()) < abs(gamma_i.real()))
        {
            if (gamma_i.real() < 0)
            {
                gamma_i = -gamma_i;
            }
        }
        gamma_(i) = gamma_i;
    }

    keff = gamma_;
#else
    keff = ces.eigenvalues().unaryExpr(GammaSqrtOp());    // eigenvalue square root
#endif

    eigvecE = ces.eigenvectors();
    eigvecH = Q * eigvecE * keff.cwiseInverse().asDiagonal();
}

void WaveEquationSolver::solve(
    scalar lambda, const Eigen::MatrixXcs& P, const Eigen::MatrixXcs& Q,
    Eigen::MatrixXcs& PQ,
    Eigen::MatrixXcs& eigvecE, Eigen::MatrixXcs& eigvecH, 
    Eigen::VectorXcs& keff) const
{
    PQ = P*Q;
    MKLEigenSolver ces;
    ces.compute(PQ);

    keff = ces.eigenvalues().unaryExpr(GammaSqrtOp());    // eigenvalue square root
    eigvecE = ces.eigenvectors();
    eigvecH = Q * eigvecE * keff.cwiseInverse().asDiagonal();
}
