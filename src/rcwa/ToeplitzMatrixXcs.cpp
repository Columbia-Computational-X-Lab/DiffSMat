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

#include "ToeplitzMatrixXcs.h"

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

ToeplitzMatrixXcs::ToeplitzMatrixXcs(const Eigen::VectorXcs& vec)
{
    n_ = (vec.size() + 1) / 2;
    data_.resize(2 * n_ - 1);
    for (int i = 0; i < 2 * n_ - 1; ++i)
    {
        int diff = n_ - i - 1;
        data_(i) = vec(n_ + diff - 1);
    }
}

scalex ToeplitzMatrixXcs::operator()(int row, int col) const
{
    return data_(n_ - 1 - (row - col));
}


Eigen::MatrixXcs ToeplitzMatrixXcs::toDense() const
{
    MatrixXcs result(n_, n_);

    for (int i = 0; i < n_; ++i)
    {
        for (int j = 0; j < n_; ++j)
        {
            result(i, j) = (*this)(i, j);
        }
    }

    return result;
}

// need optimization
Eigen::MatrixXcs ToeplitzMatrixXcs::inverse() const
{
    return this->toDense().inverse();
}