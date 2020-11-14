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


class ToeplitzMatrixXcs
{
private:
    // size of matrix should be (n_, n_)
    int n_;

    // real storage only need (2 * n_ - 1)
    // ====================================================
    // Storage
    // ====================================================
    // * * * * 
    // *      
    // *       
    // *      
    // ====================================================
    Eigen::VectorXcs data_;

public:
    ToeplitzMatrixXcs(int n): n_(n)
    {
        data_.resize(2 * n_ - 1);
    }
    
    ToeplitzMatrixXcs(const Eigen::VectorXcs& vec);
    scalex operator()(int row, int col) const;
    Eigen::MatrixXcs toDense() const;
    Eigen::MatrixXcs inverse() const;
};