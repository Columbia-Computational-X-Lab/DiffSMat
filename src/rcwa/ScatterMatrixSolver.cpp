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

#include "ScatterMatrixSolver.h"
#include "WaveEquationCoeff.h"
#include "WaveEquationSolver.h"

#include <iostream>
#include <Eigen/Dense>
#include "utils/timer.hpp"

using namespace std;
using namespace Eigen;

void ScatterMatrixSolver::solve(
    scalar lambda,
    scalar wx,
    scalar wy,
    scalar wz,
    Eigen::MatrixXcs& Tuu,
    Eigen::MatrixXcs& Rud,
    Eigen::MatrixXcs& Rdu,
    Eigen::MatrixXcs& Tdd)
{
    int nDim = 2 * nx_ * ny_;
    
    MatrixXcs W0, V0;
    VectorXcs gamma0;
    evaluateHomogeneousLayer(lambda,
        e_Vacuum,
        W0,
        V0,
        gamma0);

    WaveEquationSolver wes(nx_, ny_, Lx_, Ly_);
    MatrixXcs W, V;
    VectorXcs gamma;
    wes.solve(lambda, wx, wy, W, V, gamma);
    
    // -----------------------------------------------
    auto t1 = sploosh::now();
#if 0
    VectorXcs X(nDim);
    for (int j = 0; j < nDim; ++j)
    {
        X(j) = exp(scalex(0, 1) * gamma(j) * wz);
    }

    CompleteOrthogonalDecomposition<MatrixXcs> solverE(W);
    MatrixXcs F0 = solverE.solve(W0);
    CompleteOrthogonalDecomposition<MatrixXcs> solverH(V);
    MatrixXcs G0 = solverH.solve(V0);

    const MatrixXcs& A = F0 + G0;
    const MatrixXcs& B = F0 - G0;

    PartialPivLU<MatrixXcs> solverA(A);
    MatrixXcs invA = solverA.inverse();

    MatrixXcs invAX = invA * X.asDiagonal();
    MatrixXcs invAB = invA * B;

    PartialPivLU<MatrixXcs> solverLhs(A - X.asDiagonal() * B * invAX * B);
    MatrixXcs invLhs = solverLhs.inverse();

    Tuu = invLhs * X.asDiagonal() * (A - B * invAB);
    Rud = invLhs * (X.asDiagonal() * B * invAX * A - B);

    // -----------------------------------------------
#else
    // exp^{j*\Lambda*wz}
    
    const scalar invk0 = lambda / (2.*Pi);
    const scalar ss = invk0*wz;
    VectorXcs X = gamma.unaryExpr([&ss](const scalex& g) { return exp(scalex(0,1)*g*ss); });

    PartialPivLU<MatrixXcs> invSolver(nDim);

    invSolver.compute(W);
    MatrixXcs F0 = invSolver.solve(W0); // W^{-1}W_0

    invSolver.compute(V);
    MatrixXcs G0 = invSolver.solve(V0); // V^{-1}V_0

    const MatrixXcs& A = F0 + G0;
    const MatrixXcs& B = F0 - G0;
    const MatrixXcs Xm = X.asDiagonal();

    invSolver.compute(A); // invA

    MatrixXcs invAX = invSolver.solve( Xm );    // invA*X
    MatrixXcs invAB = invSolver.solve( B );     // invA*B

    invSolver.compute(A - X.asDiagonal() * B * invAX * B);

    Tuu = invSolver.solve(X.asDiagonal() * (A - B * invAB));
    Rud = invSolver.solve(X.asDiagonal() * B * invAX * A - B);
#endif
    auto t2 = sploosh::now();
    cout << "TIME SPEND: " << sploosh::duration_milli_d(t1, t2) << " ms" << endl;

    Tdd = Tuu;
    Rdu = Rud;


    // RTCM algorithm
    // // one interface

    // PartialPivLU<MatrixXcs> solverTau0(F0 + G0);
    // MatrixXcs tau0 = solverTau0.inverse();
    
    // Rud = - 2. * G0 * tau0;
    // Rud.diagonal().array() += scalar(1.0);
    // Tdd = 2. * tau0;
    // Tuu = F0 * tau0 * G0 + G0 * tau0 * F0;
    // Rdu = tau0 * (G0 - F0);

    // // another interface

    // MatrixXcs waveRud = X.asDiagonal() * Rud * X.asDiagonal();
    // MatrixXcs waveTdd = Tdd * X.asDiagonal();
    // MatrixXcs waveTuu = X.asDiagonal() * Tuu;

    // CompleteOrthogonalDecomposition<MatrixXcs> solverE1(W0);
    // MatrixXcs K1 = solverE1.solve(W);
    // CompleteOrthogonalDecomposition<MatrixXcs> solverH1(V0);
    // MatrixXcs K2 = solverH1.solve(V);

    // MatrixXcs F = K1 + K1 * waveRud;
    // MatrixXcs G = K2 - K2 * waveRud;
    // PartialPivLU<MatrixXcs> solverTau(F + G);
    // MatrixXcs tau = solverTau.inverse();
    
    // Rud = - 2. * G * tau;
    // Rud.diagonal().array() += scalar(1.0);
    // Tdd = 2. * waveTdd * tau;
    // Tuu = (F * tau * K2 + G * tau * K1) * waveTuu;
    // Rdu = Rdu + waveTdd * tau * (K2 - K1) * waveTuu;
}


void ScatterMatrixSolver::evaluateHomogeneousLayer(
    scalar lambda,
    scalex eps,
    Eigen::MatrixXcs& eigvecE,
    Eigen::MatrixXcs& eigvecH,
    Eigen::VectorXcs& keff)
{
    WaveEquationCoeff coeff(nx_, ny_, Lx_, Ly_);
    MatrixXcs Kx, Ky;
    coeff.evaluateKMatrices(Kx, Ky);
    
    scalar k0 = 2 * Pi / lambda;
    scalex k02 = scalex(k0 * k0);
    int blockDim = nx_ * ny_;

    MatrixXcs Q(2*blockDim, 2*blockDim);
    Q.topLeftCorner(blockDim, blockDim) = -Kx * Ky;
    Q.topRightCorner(blockDim, blockDim) = Kx * Kx;
    Q.topRightCorner(blockDim, blockDim).diagonal().array() -= k02 * eps;

    Q.bottomLeftCorner(blockDim, blockDim) = -Ky * Ky;
    Q.bottomLeftCorner(blockDim, blockDim).diagonal().array() += k02 * eps;

    Q.bottomRightCorner(blockDim, blockDim) = Ky * Kx;

    eigvecE.resize(2*blockDim, 2*blockDim);
    eigvecE.setIdentity();

    keff.resize(2*blockDim);
    for (int i = 0; i < blockDim; ++i)
    {
        scalex gamma_i = sqrt(eps * k02 - Kx(i, i) * Kx(i, i) - Ky(i, i) * Ky(i, i));

        if (gamma_i.imag() < 0)
        {
            gamma_i = -gamma_i;
        }

        
        if (abs(gamma_i.imag()) < abs(gamma_i.real()) && gamma_i.real() < 0)
        {
            gamma_i = -gamma_i;
        }

        keff(i) = keff(i + blockDim) = gamma_i;
    }

    //eigvecH.resize(2*blockDim, 2*blockDim);
    eigvecH = Q * keff.cwiseInverse().asDiagonal() * (1. / k0);
}
