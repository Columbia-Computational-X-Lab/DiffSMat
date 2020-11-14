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

#pragma once

#include <Eigen/Core>
#include <math.h>
#include <complex>
#include <tuple>

typedef double scalar;
typedef std::complex<scalar> scalex;

namespace Eigen
{
    typedef Eigen::Matrix<scalar, 2, 1> Vector2s;
    typedef Eigen::Matrix<scalar, 3, 1> Vector3s;
    typedef Eigen::Matrix<scalar, Eigen::Dynamic, 1> VectorXs;
    typedef Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;

    typedef Eigen::Matrix<scalar, Eigen::Dynamic, 2> MatrixX2s;
    typedef Eigen::Matrix<scalar, Eigen::Dynamic, 3> MatrixX3s;
    typedef Eigen::Matrix<scalar, 2, Eigen::Dynamic> Matrix2Xs;
    typedef Eigen::Matrix<scalar, 3, Eigen::Dynamic> Matrix3Xs;
    
    typedef Eigen::DiagonalMatrix<scalex, Eigen::Dynamic> DiagonalMatrixXcs;
    typedef Eigen::Matrix<scalex, Eigen::Dynamic, 1> VectorXcs;
    typedef Eigen::Matrix<scalex, Eigen::Dynamic, Eigen::Dynamic> MatrixXcs;

    typedef Eigen::Matrix<scalex, Eigen::Dynamic, 2> MatrixX2cs;
    typedef Eigen::Matrix<scalex, Eigen::Dynamic, 3> MatrixX3cs;
    typedef Eigen::Matrix<scalex, 2, Eigen::Dynamic> Matrix2Xcs;
    typedef Eigen::Matrix<scalex, 3, Eigen::Dynamic> Matrix3Xcs;

    typedef std::tuple<Eigen::MatrixXcs, Eigen::MatrixXcs, Eigen::MatrixXcs> MatrixDerivativeXcs;

}

constexpr scalar Pi = M_PI;
constexpr scalex e_Vacuum = scalex(1., .0);
