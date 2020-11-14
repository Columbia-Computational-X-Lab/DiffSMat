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
