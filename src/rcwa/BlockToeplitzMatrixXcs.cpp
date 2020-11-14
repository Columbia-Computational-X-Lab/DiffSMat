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

#include "BlockToeplitzMatrixXcs.h"

using namespace Eigen;

BlockToeplitzMatrixXcs::BlockToeplitzMatrixXcs(const MatrixXcs& mat)
{
    nx_ = (mat.rows() + 1) / 2;
    ny_ = (mat.cols() + 1) / 2;

    data_.resize((2 * nx_ - 1), (2 * ny_ - 1));

    for (int i = 0; i < 2 * nx_ - 1; ++i)  // block iteration
    {
        for (int j = 0; j < 2 * ny_ - 1; ++j) // in-block iteration
        {
            int diffi = nx_ - i - 1; // m - j
            int diffj = ny_ - j - 1; // n - l
            data_(i, j) = mat(nx_ + diffi - 1, ny_ + diffj - 1);
        }
    }
}

scalex BlockToeplitzMatrixXcs::operator()(int row, int col) const
{
    int m = row / ny_;
    int n = row % ny_;
    int j = col / ny_;
    int l = col % ny_;

    int diffi = m - j;
    int diffj = n - l;

    return data_(nx_ - 1 - diffi, ny_ - 1 - diffj);
}

MatrixXcs BlockToeplitzMatrixXcs::toDense() const
{
    MatrixXcs result(nx_ * ny_, nx_ * ny_);

    for (int i = 0; i < nx_ * ny_; ++i)
    {
        for (int j = 0; j < nx_ * ny_; ++j)
        {
            result(i, j) = (*this)(i, j);
        }
    }

    return result;
}



Eigen::MatrixXcs BlockToeplitzMatrixXcs::operator*(const Eigen::DiagonalMatrixXcs& rhs) const
{
    MatrixXcs result(this->toDense());
    for (int i = 0; i < result.cols(); ++i)
    {
        result.col(i) *= rhs.diagonal()(i);
    }

    return result;
}