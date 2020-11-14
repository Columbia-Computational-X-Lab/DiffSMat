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
#include "core/MatrixDerivative.hpp"
#include "dcomplex/MatrixExponential.hpp"
#include <unsupported/Eigen/MatrixFunctions>

using namespace std;
using namespace dtmm;

static void test1()
{
    using namespace Eigen;

    typedef complex<double> CompD;

    Matrix< CompD, 2, 2 > m1;
    m1(0, 0) = CompD(1,2)/10.;
    m1(0, 1) = CompD(4,1)/10.;
    m1(1, 0) = CompD(5,-3)/10.;
    m1(1, 1) = CompD(-3,7)/10.;
    
    Matrix< CompD, 2, 2 > delta;
    delta(0, 0) = CompD(1,0);
    delta(0, 1) = CompD(0,0);
    delta(1, 0) = CompD(0,1);
    delta(1, 1) = CompD(1,1);

    Matrix< CompD, 4, 4 > m2, U, V;
    m2.setZero();
    m2.topLeftCorner(2,2) = m1;
    m2.topRightCorner(2,2) = delta;
    m2.bottomRightCorner(2,2) = m1;

    Eigen::internal::matrix_exp_pade13(m2, U, V);

    Matrix< CompD, 2, 2 > Ua, Ub, Va, Vb;

    dtmm::internal::matrix_exp_pade13(m1, delta, Ua, Ub, Va, Vb);
    cout << (U.topLeftCorner(2,2) - Ua).cwiseAbs().maxCoeff() << endl;
    cout << (U.topRightCorner(2,2) - Ub).cwiseAbs().maxCoeff() << endl;
    cout << (V.topLeftCorner(2,2) - Va).cwiseAbs().maxCoeff() << endl;
    cout << (V.topRightCorner(2,2) - Vb).cwiseAbs().maxCoeff() << endl;

    /*
    cout << U.topLeftCorner(2,2) << endl;
    cout << " ---------" << endl;
    cout << U.topRightCorner(2,2) << endl;
    cout << " ------------------- " << endl;
    cout << V.topLeftCorner(2,2) << endl;
    cout << " ---------" << endl;
    cout << V.topRightCorner(2,2) << endl;

    cout << endl << " ~~~~ " << endl << endl;
    cout << Ua << endl;
    cout << " ---------" << endl;
    cout << Ub << endl;
    cout << " ------------------- " << endl;
    cout << Va << endl;
    cout << " ---------" << endl;
    cout << Vb << endl;
    */
}

static void test2()
{
    using namespace Eigen;

    typedef complex<double> CompD;

    Matrix< CompD, 2, 2 > m1;
    m1(0, 0) = CompD(1,2)/1.;
    m1(0, 1) = CompD(4,1)/1.;
    m1(1, 0) = CompD(5,-3)/1.;
    m1(1, 1) = CompD(-3,7)/1.;
    
    const double h0 = 1E-7;
    Matrix< CompD, 2, 2 > delta;
    delta(0, 0) = CompD(2,0);
    delta(0, 1) = CompD(0,0);
    delta(1, 0) = CompD(0,1);
    delta(1, 1) = CompD(1,4);
    //cout << delta.cwiseAbs().colwise().sum() << endl;
    
    Matrix< CompD, 4, 4 > m2, U, V;
    m2.setZero();
    m2.topLeftCorner(2,2) = m1;
    m2.topRightCorner(2,2) = delta;
    m2.bottomRightCorner(2,2) = m1;

    cout << ((m1+delta*h0).exp() - (m1-delta*h0).exp())/(h0*2) << endl;
    cout << " ------ " << endl;

    Matrix< CompD, 4, 4 > ret = m2.exp();
    cout << ret.topRightCorner(2,2) << endl;

    cout << " ------ " << endl;
    Matrix< CompD, 2, 2 > deriv = mat_exp_deriv(m1, delta);
    cout << deriv << endl;
}

int main()
{
    test1();
    test2();
}
