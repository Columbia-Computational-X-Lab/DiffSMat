#include "WaveEquationCoeff.h"
#include "FourierSolver2D.h"
#include "FourierSolverPolygon.h"
#include <iostream>
#include <Eigen/Dense>
#include "dispersion.h"

using namespace std;
using namespace Eigen;

void WaveEquationCoeff::solve(
    scalar lambda,
    scalar wx, // width of scatterer along x direction
    scalar wy, // width of scatterer along y direction
    Eigen::MatrixXcs& P,
    Eigen::MatrixXcs& Q,
    Eigen::MatrixXcs* dPx, // partial P / partial wx
    Eigen::MatrixXcs* dPy, // partial P / partial wy
    Eigen::MatrixXcs* dQx, // partial Q / partial wx
    Eigen::MatrixXcs* dQy) // partial Q / partial wy
{
    MatrixXcs Kx, Ky;
    evaluateKMatrices(Kx, Ky);

    VectorXs coordX, coordY;
    MatrixXcs eps;
    constructLayer_(lambda, coordX, coordY, eps, wx, wy);

    // Q is only dependent on fe11, fe22
    // P is only dependent on fe33

    FourierSolver2D e11Solver(coordX, coordY);
    MatrixXcs fe11;
    e11Solver.solve(eps, nx_, ny_, DiscontinuousX, fe11);

    FourierSolver2D e22Solver(coordX, coordY);
    MatrixXcs fe22;
    e22Solver.solve(eps, nx_, ny_, DiscontinuousY, fe22);

    FourierSolver2D e33Solver(coordX, coordY);
    MatrixXcs fe33;
    e33Solver.solve(eps, nx_, ny_, ContinuousXY, fe33);

    // dfe -> coord(1), coord(2) -> w
    MatrixXcs dfe11x, dfe11y;
    if (dQx != nullptr || dQy != nullptr)
    {
        MatrixXcs dfe11x_1, dfe11x_2, dfe11y_1, dfe11y_2;
        e11Solver.solve(eps, nx_, ny_, DiscontinuousX, fe11, 
            1, &dfe11x_1, 1, &dfe11y_1);
        e11Solver.solve(eps, nx_, ny_, DiscontinuousX, fe11,
            2, &dfe11x_2, 2, &dfe11y_2);
        dfe11x = 0.5 * (dfe11x_2 - dfe11x_1);
        dfe11y = 0.5 * (dfe11y_2 - dfe11y_1);
    }

    MatrixXcs dfe22x, dfe22y;
    if (dQx != nullptr || dQy != nullptr)
    {
        MatrixXcs dfe22x_1, dfe22x_2, dfe22y_1, dfe22y_2;
        e22Solver.solve(eps, nx_, ny_, DiscontinuousY, fe22, 
            1, &dfe22x_1, 1, &dfe22y_1);
        e22Solver.solve(eps, nx_, ny_, DiscontinuousY, fe22,
            2, &dfe22x_2, 2, &dfe22y_2);
        dfe22x = 0.5 * (dfe22x_2 - dfe22x_1);
        dfe22y = 0.5 * (dfe22y_2 - dfe22y_1);
    }

    MatrixXcs dfe33x, dfe33y;
    if (dPx != nullptr || dPy != nullptr)
    {
        MatrixXcs dfe33x_1, dfe33x_2, dfe33y_1, dfe33y_2;
        e33Solver.solve(eps, nx_, ny_, ContinuousXY, fe33, 
            1, &dfe33x_1, 1, &dfe33y_1);
        e33Solver.solve(eps, nx_, ny_, ContinuousXY, fe33,
            2, &dfe33x_2, 2, &dfe33y_2);
        dfe33x = 0.5 * (dfe33x_2 - dfe33x_1);
        dfe33y = 0.5 * (dfe33y_2 - dfe33y_1);        
    }

    BDCSVD<MatrixXcs> SVDSolver(fe33, ComputeThinU | ComputeThinV);
    const MatrixXcs& F1h = SVDSolver.solve(Kx);
    const MatrixXcs& F2h = SVDSolver.solve(Ky);
    
    scalar k0 = 2 * Pi / lambda;
    scalex k02 = scalex(k0 * k0);

    int blockDim = nx_ * ny_;

    P.resize(2*blockDim, 2*blockDim);
    P.topLeftCorner(blockDim, blockDim) = Kx * F2h;
    P.topRightCorner(blockDim, blockDim) = - Kx * F1h;
    P.topRightCorner(blockDim, blockDim).diagonal().array() += k02;
    P.bottomLeftCorner(blockDim, blockDim) = Ky * F2h;
    P.bottomLeftCorner(blockDim, blockDim).diagonal().array() -= k02;
    P.bottomRightCorner(blockDim, blockDim) = -Ky * F1h;

    if (dPx != nullptr)
    {
        dPx->resize(2*blockDim, 2*blockDim);

        const MatrixXcs& temp = SVDSolver.solve(dfe33x);

        dPx->topLeftCorner(blockDim, blockDim) = Kx * (-temp * F2h);
        dPx->topRightCorner(blockDim, blockDim) = - Kx * (-temp * F1h);
        dPx->bottomLeftCorner(blockDim, blockDim) = Ky * (-temp * F2h);
        dPx->bottomRightCorner(blockDim, blockDim) = -Ky * (-temp * F1h);        
    }

    if (dPy != nullptr)
    {
        dPy->resize(2*blockDim, 2*blockDim);

        const MatrixXcs& temp = SVDSolver.solve(dfe33y);
        dPy->topLeftCorner(blockDim, blockDim) = Kx * (-temp * F2h);
        dPy->topRightCorner(blockDim, blockDim) = - Kx * (-temp * F1h);
        dPy->bottomLeftCorner(blockDim, blockDim) = Ky * (-temp * F2h);
        dPy->bottomRightCorner(blockDim, blockDim) = -Ky * (-temp * F1h);
    }

    Q.resize(2*blockDim, 2*blockDim);
    Q.topLeftCorner(blockDim, blockDim) = -Kx * Ky;
    Q.topRightCorner(blockDim, blockDim) = Kx * Kx - k02 * fe22;
    Q.bottomLeftCorner(blockDim, blockDim) = -Ky * Ky + k02 * fe11;
    Q.bottomRightCorner(blockDim, blockDim) = Ky * Kx;

    if (dQx != nullptr)
    {
        dQx->resize(2*blockDim, 2*blockDim);
        dQx->setZero();

        dQx->topRightCorner(blockDim, blockDim) = -k02 * dfe22x;
        dQx->bottomLeftCorner(blockDim, blockDim) = k02 * dfe11x;
    }

    if (dQy != nullptr)
    {
        dQy->resize(2*blockDim, 2*blockDim);
        dQy->setZero();

        dQy->topRightCorner(blockDim, blockDim) = -k02 * dfe22y;
        dQy->bottomLeftCorner(blockDim, blockDim) = k02 * dfe11y;
    }
}


void WaveEquationCoeff::solve(
    scalar lambda,
    const Eigen::VectorXs& alphas,
    Eigen::MatrixXcs& P,
    Eigen::MatrixXcs& Q,        
    const std::vector<int>* layouts,
    std::vector<Eigen::MatrixXcs>* dPs,
    std::vector<Eigen::MatrixXcs>* dQs)
{
    MatrixXcs Kx, Ky;
    evaluateKMatrices(Kx, Ky);

    // here we use the simplest factorization without considering
    // about convergence
    FourierSolverPolygon epsSolver(alphas, Lx_, Ly_);
   
    MatrixXcs feps;
    std::vector<MatrixXcs> dfeps;

    scalex e_Si = permSi(lambda);
    epsSolver.solve(e_Si, e_Vacuum, nx_, ny_, feps, layouts, &dfeps);


    BDCSVD<MatrixXcs> SVDSolver(feps, ComputeThinU | ComputeThinV);
    const MatrixXcs& F1h = SVDSolver.solve(Kx);
    const MatrixXcs& F2h = SVDSolver.solve(Ky);
    
    scalar k0 = 2 * Pi / lambda;
    scalex k02 = scalex(k0 * k0);

    int blockDim = nx_ * ny_;

    P.resize(2*blockDim, 2*blockDim);
    P.topLeftCorner(blockDim, blockDim) = Kx * F2h;
    P.topRightCorner(blockDim, blockDim) = - Kx * F1h;
    P.topRightCorner(blockDim, blockDim).diagonal().array() += k02;
    P.bottomLeftCorner(blockDim, blockDim) = Ky * F2h;
    P.bottomLeftCorner(blockDim, blockDim).diagonal().array() -= k02;
    P.bottomRightCorner(blockDim, blockDim) = -Ky * F1h;

    if (layouts != nullptr)
    {
        int Nd = layouts->size();
        dPs->resize(Nd);

        for (int k = 0; k < Nd; ++k)
        {
            auto& dP = (*dPs)[k];
            dP.resize(2*blockDim, 2*blockDim);
            const MatrixXcs& temp = SVDSolver.solve(dfeps[k]);
            dP.topLeftCorner(blockDim, blockDim) = Kx * (-temp * F2h);
            dP.topRightCorner(blockDim, blockDim) = - Kx * (-temp * F1h);
            dP.bottomLeftCorner(blockDim, blockDim) = Ky * (-temp * F2h);
            dP.bottomRightCorner(blockDim, blockDim) = -Ky * (-temp * F1h);             
        }
    }

    Q.resize(2*blockDim, 2*blockDim);
    Q.topLeftCorner(blockDim, blockDim) = -Kx * Ky;
    Q.topRightCorner(blockDim, blockDim) = Kx * Kx - k02 * feps;
    Q.bottomLeftCorner(blockDim, blockDim) = -Ky * Ky + k02 * feps;
    Q.bottomRightCorner(blockDim, blockDim) = Ky * Kx;

    if (layouts != nullptr)
    {
        int Nd = layouts->size();
        dQs->resize(Nd);

        for (int k = 0; k < Nd; ++k)
        {
            auto& dQ = (*dQs)[k];
            dQ.resize(2*blockDim, 2*blockDim);
            dQ.setZero();
            
            dQ.topRightCorner(blockDim, blockDim) = -k02 * dfeps[k];
            dQ.bottomLeftCorner(blockDim, blockDim) = k02 * dfeps[k];                      
        }
    }
}

void WaveEquationCoeff::solve(
    scalar lambda,
    const struct LayerStructure& layer,
    int layoutX,
    int layoutY,
    Eigen::MatrixXcs& P,
    Eigen::MatrixXcs& Q,
    Eigen::MatrixXcs* dPx, // partial P / partial x_i
    Eigen::MatrixXcs* dPy, // partial P / partial y_i
    Eigen::MatrixXcs* dQx, // partial Q / partial x_i
    Eigen::MatrixXcs* dQy) // partial Q / partial y_i)
{
    MatrixXcs Kx, Ky;
    evaluateKMatrices(Kx, Ky);

    VectorXs coordX, coordY;
    MatrixXcs eps;
    
    coordX = layer.coordX;
    coordY = layer.coordY;
    eps = layer.eps;

    // Q is only dependent on fe11, fe22
    // P is only dependent on fe33

    FourierSolver2D e11Solver(coordX, coordY);
    MatrixXcs fe11;
    e11Solver.solve(eps, nx_, ny_, DiscontinuousX, fe11);

    FourierSolver2D e22Solver(coordX, coordY);
    MatrixXcs fe22;
    e22Solver.solve(eps, nx_, ny_, DiscontinuousY, fe22);

    FourierSolver2D e33Solver(coordX, coordY);
    MatrixXcs fe33;
    e33Solver.solve(eps, nx_, ny_, ContinuousXY, fe33);

    // dfe -> coord(1), coord(2) -> w
    MatrixXcs dfe11x, dfe11y;
    if (dQx != nullptr || dQy != nullptr)
    {
        e11Solver.solve(eps, nx_, ny_, DiscontinuousX, fe11, 
            layoutX, &dfe11x, layoutY, &dfe11y);
    }

    MatrixXcs dfe22x, dfe22y;
    if (dQx != nullptr || dQy != nullptr)
    {
        e22Solver.solve(eps, nx_, ny_, DiscontinuousY, fe22, 
            layoutX, &dfe22x, layoutY, &dfe22y);
    }

    MatrixXcs dfe33x, dfe33y;
    if (dPx != nullptr || dPy != nullptr)
    {
        e33Solver.solve(eps, nx_, ny_, ContinuousXY, fe33, 
            layoutX, &dfe33x, layoutY, &dfe33y);
    }

    BDCSVD<MatrixXcs> SVDSolver(fe33, ComputeThinU | ComputeThinV);
    const MatrixXcs& F1h = SVDSolver.solve(Kx);
    const MatrixXcs& F2h = SVDSolver.solve(Ky);
    
    scalar k0 = 2 * Pi / lambda;
    scalex k02 = scalex(k0 * k0);

    int blockDim = nx_ * ny_;

    P.resize(2*blockDim, 2*blockDim);
    P.topLeftCorner(blockDim, blockDim) = Kx * F2h;
    P.topRightCorner(blockDim, blockDim) = - Kx * F1h;
    P.topRightCorner(blockDim, blockDim).diagonal().array() += k02;
    P.bottomLeftCorner(blockDim, blockDim) = Ky * F2h;
    P.bottomLeftCorner(blockDim, blockDim).diagonal().array() -= k02;
    P.bottomRightCorner(blockDim, blockDim) = -Ky * F1h;

    if (dPx != nullptr)
    {
        dPx->resize(2*blockDim, 2*blockDim);

        const MatrixXcs& temp = SVDSolver.solve(dfe33x);

        dPx->topLeftCorner(blockDim, blockDim) = Kx * (-temp * F2h);
        dPx->topRightCorner(blockDim, blockDim) = - Kx * (-temp * F1h);
        dPx->bottomLeftCorner(blockDim, blockDim) = Ky * (-temp * F2h);
        dPx->bottomRightCorner(blockDim, blockDim) = -Ky * (-temp * F1h);       
    }

    if (dPy != nullptr)
    {
        dPy->resize(2*blockDim, 2*blockDim);

        const MatrixXcs& temp = SVDSolver.solve(dfe33y);
        dPy->topLeftCorner(blockDim, blockDim) = Kx * (-temp * F2h);
        dPy->topRightCorner(blockDim, blockDim) = - Kx * (-temp * F1h);
        dPy->bottomLeftCorner(blockDim, blockDim) = Ky * (-temp * F2h);
        dPy->bottomRightCorner(blockDim, blockDim) = -Ky * (-temp * F1h);
    }

    Q.resize(2*blockDim, 2*blockDim);
    Q.topLeftCorner(blockDim, blockDim) = -Kx * Ky;
    Q.topRightCorner(blockDim, blockDim) = Kx * Kx - k02 * fe22;
    Q.bottomLeftCorner(blockDim, blockDim) = -Ky * Ky + k02 * fe11;
    Q.bottomRightCorner(blockDim, blockDim) = Ky * Kx;

    if (dQx != nullptr)
    {
        dQx->resize(2*blockDim, 2*blockDim);
        dQx->setZero();

        dQx->topRightCorner(blockDim, blockDim) = -k02 * dfe22x;
        dQx->bottomLeftCorner(blockDim, blockDim) = k02 * dfe11x;
    }

    if (dQy != nullptr)
    {
        dQy->resize(2*blockDim, 2*blockDim);
        dQy->setZero();

        dQy->topRightCorner(blockDim, blockDim) = -k02 * dfe22y;
        dQy->bottomLeftCorner(blockDim, blockDim) = k02 * dfe11y;
    }       
}

void WaveEquationCoeff::evaluateKMatrices(
    Eigen::MatrixXcs& Kx,
    Eigen::MatrixXcs& Ky)
{
    DiagonalMatrixXcs alpha;
    DiagonalMatrixXcs beta;

    evaluateAlphaBeta_(alpha, beta);

    int nDim = nx_ * ny_;

    // assume no PML
    Kx.resize(nDim, nDim);
    Ky.resize(nDim, nDim);
    Kx.setZero();
    Ky.setZero();
    Kx.diagonal() = alpha.diagonal();
    Ky.diagonal() = beta.diagonal();
}


void WaveEquationCoeff::evaluateAlphaBeta_(
    Eigen::DiagonalMatrixXcs& alpha,
    Eigen::DiagonalMatrixXcs& beta)
{
    alpha.resize(nx_ * ny_);
    beta.resize(nx_ * ny_);

    int Nx = (nx_ - 1) / 2;
    int Ny = (ny_ - 1) / 2;

    scalar Kx = 2 * Pi / Lx_;
    scalar Ky = 2 * Pi / Ly_;

    for (int i = 0;  i < nx_; ++i)
    {
        for (int j = 0; j < ny_; ++j)
        {
            scalar m = i - Nx;
            scalar n = j - Ny;

            alpha.diagonal()(ny_ * i + j) = Kx * m;
            beta.diagonal()(ny_ * i + j) = Ky * n;
        }
    }

}

void WaveEquationCoeff::constructLayer_(
     scalar lambda,
    Eigen::VectorXs& coordX,
    Eigen::VectorXs& coordY,
    Eigen::MatrixXcs& eps,
    scalar wx, scalar wy)
{
    scalex e_Si = permSi(lambda);

    if (wx > Lx_ || wy > Ly_)
    {
        cerr << "[ERROR]: The size of scatterer is too large!\n";
        return;
    }

    coordX.resize(4);
    coordY.resize(4);

    coordX(0) = .0;
    coordX(1) = 0.5 * Lx_ - 0.5 * wx;
    coordX(2) = 0.5 * Lx_ + 0.5 * wx;
    coordX(3) = Lx_;

    coordY(0) = .0;
    coordY(1) = 0.5 * Ly_ - 0.5 * wy;
    coordY(2) = 0.5 * Ly_ + 0.5 * wy;
    coordY(3) = Ly_;

    eps.resize(3, 3);
    eps(0, 0) = eps(0, 1) = eps(0, 2) = 
    eps(1, 0) =             eps(1, 2) =
    eps(2, 0) = eps(2, 1) = eps(2, 2) = e_Vacuum;

    eps(1, 1) = e_Si;
}
