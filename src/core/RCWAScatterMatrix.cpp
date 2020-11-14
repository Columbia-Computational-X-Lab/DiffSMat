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
#include "RCWAScatterMatrix.h"
#include "rcwa/WaveEquationCoeff.h"
#include "utils/timer.hpp"

namespace dtmm
{

void RCWAScatterMatrix::compute(scalar lambda, scalar wz,
        const RCWAScatterMatrix::MatrixType& P,
        const RCWAScatterMatrix::MatrixType& Q)
{
    using namespace Eigen;

    int nDim = 2 * nx_ * ny_;
    VectorXcs gamma0;
    evaluateHomogeneousLayer(lambda, e_Vacuum, W0_, V0_, gamma0);

    WaveEquationSolver wes(nx_, ny_, Lx_, Ly_);
    wes.solve(lambda, P, Q, PQ_, W_, V_, gamma_);

    // -----------------------------------------------
    //auto t1 = sploosh::now();

    const scalar invk0 = lambda / (2.*Pi);
    normalizedLz_ = invk0*wz;   // wz / k0
    const scalar ss = normalizedLz_;
    X_ = gamma_.unaryExpr([&ss](const scalex& g) { return exp(scalex(0,1)*g*ss); });

    invWSol_.compute(W_);
    MatrixXcs F0 = invWSol_.solve(W0_);  // W^{-1}W_0

    PartialPivLU<MatrixXcs> invSolver(nDim);
    invSolver.compute(V_);
    MatrixXcs G0 = invSolver.solve(V0_); // V^{-1}V_0

    // const MatrixXcs& A = F0 + G0;
    // const MatrixXcs& B = F0 - G0;
    A_ = F0 + G0;
    B_ = F0 - G0;

    // invSolver.compute(A); // invA
    invASol_.compute(A_);

    invAB_  = invASol_.solve( B_ );  // invA*B
    MatrixXcs XA = X_.asDiagonal() * A_;       // XA
    MatrixXcs XB = X_.asDiagonal() * B_;       // XB
    invAXA_ = invASol_.solve(XA);
    invAXB_ = invASol_.solve(XB);

    const MatrixXcs I = MatrixXcs::Identity(nDim, nDim);
    I_minus_invAXB2_ = I - invAXB_*invAXB_;
    inv_I_minus_invAXB2_Sol_.compute(I_minus_invAXB2_);

    invAXB_invAXA_invAB_ = invAXB_*invAXA_ - invAB_;
    Rud_ = inv_I_minus_invAXB2_Sol_.solve(invAXB_invAXA_invAB_);
    invAXA_invAXB_invAB_ = invAXA_ - invAXB_*invAB_;
    Tuu_ = inv_I_minus_invAXB2_Sol_.solve(invAXA_invAXB_invAB_);

    isComputed_ = true;
    //auto t2 = sploosh::now();
    //std::cout << "TIME SPEND: " << sploosh::duration_milli_d(t1, t2) << " ms" << std::endl;
}

void RCWAScatterMatrix::compute_dz(scalar lambda, scalar wz)
{
    using namespace Eigen;

    if (!isComputed_) {
        std::cerr << "The scattering matrix has not been computed!" << std::endl;
        exit(-1);
    }

    const scalar invk0 = lambda / (2.*Pi);
    const scalar ss = normalizedLz_;   
    VectorXcs dX = gamma_.unaryExpr(
        [&ss, &invk0](const scalex& g) { 
            return scalex(0.,1.) * g * invk0 * exp(scalex(0,1)*g*ss); 
        });

    MatrixXcs invAdXA = invASol_.solve(dX.asDiagonal() * A_);
    MatrixXcs invAdXB = invASol_.solve(dX.asDiagonal() * B_);

    MatrixXcs dI_minus_invAXB2 = - invAdXB*invAXB_ - invAXB_*invAdXB;
    MatrixXcs dinvAXB_invAXA_invAB = invAdXB*invAXA_ + invAXB_*invAdXA;
    MatrixXcs dinvAXA_invAXB_invAB = invAdXA - invAdXB * invAB_;
    // -(X)^-1 dX (X)^-1 Y + (X)^-1 dY

    dRudz_ = inv_I_minus_invAXB2_Sol_.solve(
            dinvAXB_invAXA_invAB - dI_minus_invAXB2 * Rud_) ;
    dTuuz_ = inv_I_minus_invAXB2_Sol_.solve(
            dinvAXA_invAXB_invAB - dI_minus_invAXB2 * Tuu_);
}

void RCWAScatterMatrix::evaluateHomogeneousLayer(
    scalar lambda, scalex eps,
    RCWAScatterMatrix::MatrixType& eigvecE,
    RCWAScatterMatrix::MatrixType& eigvecH,
    RCWAScatterMatrix::VectorType& keff)
{
    WaveEquationCoeff coeff(nx_, ny_, Lx_, Ly_);
    MatrixType Kx, Ky;
    coeff.evaluateKMatrices(Kx, Ky);
    
    scalar k0 = 2 * Pi / lambda;
    scalex k02 = scalex(k0 * k0);
    int blockDim = nx_ * ny_;

    MatrixType Q(2*blockDim, 2*blockDim);
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

} // end of namespace dtmm
