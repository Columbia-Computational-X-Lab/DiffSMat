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

#include "defns.h"
#include <Eigen/Core>


class BlockToeplitzMatrixXcs
{
private:
    // size of matrix should be (nx_ * ny_) * (nx_ * ny_)
    int nx_;
    int ny_;

    // real storage only need (2 * nx_ - 1) * (2 * ny_ - 1)
    // ====================================================
    // Storage
    // ====================================================
    // * * * * * * * * * * * * * * * * 
    // *       *       *       *
    // *       *       *       *
    // *       *       *       *
    // * * * *
    // * 
    // *
    // *
    // * * * *
    // * 
    // *
    // *
    // * * * *
    // * 
    // *
    // *    
    // ====================================================

    Eigen::MatrixXcs data_;

public:
    BlockToeplitzMatrixXcs(int nx, int ny): nx_(nx), ny_(ny)
    {
        data_.resize(2 * nx_ - 1, 2 * ny_ - 1);
    }
    
    BlockToeplitzMatrixXcs(const Eigen::MatrixXcs& mat);
    scalex operator()(int row, int col) const;
    Eigen::MatrixXcs toDense() const;

    Eigen::MatrixXcs operator*(const Eigen::DiagonalMatrixXcs& rhs) const;
};
