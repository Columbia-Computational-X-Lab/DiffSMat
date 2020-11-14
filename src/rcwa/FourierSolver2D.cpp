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

#include "FourierSolver2D.h"

#include "FourierSolver1D.h"
#include "ToeplitzMatrixXcs.h"
#include "BlockToeplitzMatrixXcs.h"
#include <iostream>


using namespace std;
using namespace Eigen;

void FourierSolver2D::solve(const Eigen::MatrixXcs& f,
    int nx, int ny, FourierSolverOption opt,
    Eigen::MatrixXcs& Ff, int layoutX,
    Eigen::MatrixXcs * dFfx, int layoutY,
    Eigen::MatrixXcs * dFfy)
{
    Ff.resize(nx*ny, nx*ny);
    Ff.setZero();

    if (isValidX(layoutX, dFfx))
    {
        dFfx->resize(nx*ny, nx*ny);
        dFfx->setZero();
    }

    if (isValidY(layoutY, dFfy))
    {
        dFfy->resize(nx*ny, nx*ny);
        dFfy->setZero();
    }

    if (opt == DiscontinuousX)
    {
        FourierSolver1D xSolver(coordX_, dx_);
        for (int i = 0; i < coordY_.size()-1; ++i)
        {
            VectorXcs invFi = f.row(i).cwiseInverse();
            VectorXcs fInvFi;
            xSolver.solve(invFi, 2*nx-1, fInvFi);
            ToeplitzMatrixXcs tFInvFi(fInvFi);
            MatrixXcs tFFi = tFInvFi.inverse();

            MatrixXcs dtFFi;
            if (isValidX(layoutX, dFfx))
            {
                VectorXcs dfInvFi;
                xSolver.dSolve(invFi, 2*nx-1, layoutX, dfInvFi);
                ToeplitzMatrixXcs dtFInvFi(dfInvFi);
                dtFFi = -tFFi * dtFInvFi.toDense() * tFFi;
            }

            VectorXs tempCoordY(2);
            tempCoordY(0) = coordY_(i);
            tempCoordY(1) = coordY_(i+1);
            FourierSolver1D ySolver(tempCoordY, dy_);
            VectorXcs fY1;
            VectorXcs identity(1);
            identity(0) = scalex(1.0, 0.0);
            ySolver.solve(identity, 2*ny-1, fY1);

            VectorXcs dfY1(2*ny-1);
            dfY1.setZero();

            if (isValidY(layoutY, dFfy))
            {
                if (layoutY == i)
                {
                    VectorXcs dfY1_0;
                    ySolver.dSolve(identity, 2*ny-1, 0, dfY1_0);
                    dfY1 += dfY1_0;
                }
                else if (layoutY == i+1)
                {
                    VectorXcs dfY1_1;
                    ySolver.dSolve(identity, 2*ny-1, 1, dfY1_1);
                    dfY1 += dfY1_1;
                }
            }

            for (int ii = 0; ii < nx; ++ii)
            {
                for (int jj = 0; jj < ny; ++jj)
                {
                    int m = ii - (nx - 1) / 2;
                    int n = jj - (ny - 1) / 2;

                    for (int kk = 0; kk < nx; ++kk)
                    {
                        for (int ll = 0; ll < ny; ++ll)
                        {
                            int j = kk - (nx - 1) / 2;
                            int l = ll - (ny - 1) / 2;
                            Ff(ii * ny + jj, kk * ny + ll) += 
                            fY1(n-l + ny - 1) * tFFi(m +(nx-1)/2, j+(nx-1)/2);

                            if (isValidX(layoutX, dFfx))
                            {
                                (*dFfx)(ii * ny + jj, kk * ny + ll) += 
                                fY1(n-l + ny - 1) * dtFFi(m +(nx-1)/2, j+(nx-1)/2);
                            }

                            if (isValidY(layoutY, dFfy))
                            {
                                (*dFfy)(ii * ny + jj, kk * ny + ll) += 
                                dfY1(n-l + ny - 1) * tFFi(m +(nx-1)/2, j+(nx-1)/2);
                            }
                        }
                    }
                }
            }
        }       
    }
    else if (opt == DiscontinuousY)
    {
        FourierSolver1D ySolver(coordY_, dy_);
        for (int i = 0; i < coordX_.size()-1; ++i)
        {
            VectorXcs invFi = f.col(i).cwiseInverse();
            VectorXcs fInvFi;
            ySolver.solve(invFi, 2*ny-1, fInvFi);
            ToeplitzMatrixXcs tFInvFi(fInvFi);
            MatrixXcs tFFi = tFInvFi.inverse();

            MatrixXcs dtFFi;
            if (isValidY(layoutY, dFfy))
            {
                VectorXcs dfInvFi;
                ySolver.dSolve(invFi, 2*ny-1, layoutY, dfInvFi);
                ToeplitzMatrixXcs dtFInvFi(dfInvFi);
                dtFFi = -tFFi * dtFInvFi.toDense() * tFFi;
            }

            VectorXs tempCoordX(2);
            tempCoordX(0) = coordX_(i);
            tempCoordX(1) = coordX_(i+1);
            FourierSolver1D xSolver(tempCoordX, dx_);
            VectorXcs fX1;
            VectorXcs identity(1);
            identity(0) = scalex(1.0, 0.0);
            xSolver.solve(identity, 2*nx-1, fX1);

            VectorXcs dfX1(2*nx-1);
            dfX1.setZero();

            if (isValidX(layoutX, dFfx))
            {
                if (layoutX == i)
                {
                    VectorXcs dfX1_0;
                    xSolver.dSolve(identity, 2*nx-1, 0, dfX1_0);
                    dfX1 += dfX1_0;
                }
                else if (layoutX == i+1)
                {
                    VectorXcs dfX1_1;
                    xSolver.dSolve(identity, 2*nx-1, 1, dfX1_1);
                    dfX1 += dfX1_1;
                }
            }

            for (int ii = 0; ii < nx; ++ii)
            {
                for (int jj = 0; jj < ny; ++jj)
                {
                    int m = ii - (nx - 1) / 2;
                    int n = jj - (ny - 1) / 2;

                    for (int kk = 0; kk < nx; ++kk)
                    {
                        for (int ll = 0; ll < ny; ++ll)
                        {
                            int j = kk - (nx - 1) / 2;
                            int l = ll - (ny - 1) / 2;
                            Ff(ii * ny + jj, kk * ny + ll) += 
                            fX1(m-j + nx - 1) * tFFi(n +(ny-1)/2, l+(ny-1)/2);

                            if (isValidX(layoutX, dFfx))
                            {
                                (*dFfx)(ii * ny + jj, kk * ny + ll) += 
                                dfX1(m-j + nx - 1) * tFFi(n +(ny-1)/2, l+(ny-1)/2);
                            }

                            if (isValidY(layoutY, dFfy))
                            {
                                (*dFfy)(ii * ny + jj, kk * ny + ll) += 
                                fX1(m-j + nx - 1) * dtFFi(n +(ny-1)/2, l+(ny-1)/2);
                            }

                        }
                    }
                }
            }           
        }       
    }
    else if (opt == ContinuousXY)
    {
        MatrixXcs fourier2D(2*nx-1, 2*ny-1);
        fourier2D.setZero();

        MatrixXcs dfourier2Dx, dfourier2Dy;
        if (isValidX(layoutX, dFfx))
        {
            dfourier2Dx.resize(2*nx-1, 2*ny-1);
            dfourier2Dx.setZero();
        }
        if (isValidY(layoutY, dFfy))
        {
            dfourier2Dy.resize(2*nx-1, 2*ny-1);
            dfourier2Dy.setZero();
        }


        FourierSolver1D ySolver(coordY_, dy_);
        for (int i = 0; i < coordX_.size()-1; ++i)
        {
            VectorXcs Fi = f.col(i);
            VectorXcs fFi;
            ySolver.solve(Fi, 2*ny-1, fFi);

            VectorXcs dfFi;
            if (isValidY(layoutY, dFfy))
            {
                ySolver.dSolve(Fi, 2*ny-1, layoutY, dfFi);
            }

            VectorXs tempCoordX(2);
            tempCoordX(0) = coordX_(i);
            tempCoordX(1) = coordX_(i+1);
            FourierSolver1D xSolver(tempCoordX, dx_);
            VectorXcs fX1;
            VectorXcs identity(1);
            identity(0) = scalex(1.0, 0.0);
            xSolver.solve(identity, 2*nx-1, fX1);

            VectorXcs dfX1(2*nx-1);
            dfX1.setZero();

            if (isValidX(layoutX, dFfx))
            {
                if (layoutX == i)
                {
                    VectorXcs dfX1_0;
                    xSolver.dSolve(identity, 2*nx-1, 0, dfX1_0);
                    dfX1 += dfX1_0;
                }
                else if (layoutX == i+1)
                {
                    VectorXcs dfX1_1;
                    xSolver.dSolve(identity, 2*nx-1, 1, dfX1_1);
                    dfX1 += dfX1_1;
                }               
            }

            fourier2D += fX1 * fFi.transpose();

            if (isValidX(layoutX, dFfx))
            {
                dfourier2Dx += dfX1 * fFi.transpose();
            }

            if (isValidY(layoutY, dFfy))
            {
                dfourier2Dy += fX1 * dfFi.transpose();
            }
        }

        Ff = BlockToeplitzMatrixXcs(fourier2D).toDense();


        if(isValidX(layoutX, dFfx))
        {
            (*dFfx) = BlockToeplitzMatrixXcs(dfourier2Dx).toDense();
        }

        if (isValidY(layoutY, dFfy))
        {
            (*dFfy) = BlockToeplitzMatrixXcs(dfourier2Dy).toDense();
        }
    }

}
