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
#include <Eigen/Core>
#include <Eigen/Dense>
#include "core/MatrixDerivative.hpp"

using namespace std;
using namespace dtmm;


static void test_inverse_deriv()
{
    using namespace Eigen;

    MatrixXd m = MatrixXd::Random(11,11);
    MatrixXd mprime = MatrixXd::Random(11,11);

    const double h = 1E-8;
    MatrixXd inv1 = (m + mprime*h).inverse();
    MatrixXd inv2 = (m - mprime*h).inverse();
    MatrixXd dEst = (inv1 - inv2)/(h*2);

    MatrixXd ret = mat_inverse_deriv(m, mprime);
    cout << (ret - dEst).cwiseAbs().maxCoeff() << endl;
}

int main()
{
    test_inverse_deriv();
    return 0;
}
