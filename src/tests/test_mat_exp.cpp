#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include "core/MatrixDerivative.hpp"
#include "dcomplex/DComplexMath.hpp"
#include "dcomplex/EigenNumTraits.hpp"
#include "dcomplex/MatrixExponential.hpp"

using namespace std;
using namespace dtmm;

static void test1()
{
    Eigen::Matrix2d mat;
    mat(0, 0) = 1;
    mat(1, 0) = -2;
    mat(0, 1) = -0.3;
    mat(1, 1) = 2;
    //cout << mat << endl;
    //cout << "   exp(mat) = " << endl;

    //Eigen::Matrix2d s = mat.sqrt();
    cout << mat.exp() << endl;

    cout << "------" << endl;
    
    //cout << s*s << endl;
    
    // dtmm::ExpmWithSqrReturnValue<Eigen::Matrix2d> ret = dtmm::expm_with_sqr(mat, mat);
    Eigen::Matrix2d m2 = mat*mat;
    Eigen::MatrixXd ret = dtmm::expm_with_sqr(mat, m2);
    cout << ret << endl;
    // cout << ret.cols() << " , " << ret.rows() << endl;

    cout << endl << " --------------------------------------------- " << endl;
}

static void test1_1()
{
    Eigen::Matrix2d m;
    m(0, 0) = 1;
    m(1, 0) = -2;
    m(0, 1) = -0.3;
    m(1, 1) = 2;
    Eigen::Matrix2d m2 = m*m;
    
    Eigen::Matrix2d aa = m.unaryExpr(dtmm::internal::MatrixExponentialScalingOp<double>(3));
    Eigen::Matrix2d bb = m2.unaryExpr(dtmm::internal::MatrixExponentialScalingOp<double>(6));
    
    cout << aa*aa << endl;
    cout << "------" << endl;
    cout << bb << endl;
    // cout << ret.cols() << " , " << ret.rows() << endl;

    cout << endl << " --------------------------------------------- " << endl;
}

static void test1_2()
{
    Eigen::Matrix3d mat;
    mat(0, 0) = 1;
    mat(0, 1) = -0.3;
    mat(0, 2) = -0.1;
    mat(1, 0) = -2;
    mat(1, 1) = 2;
    mat(1, 2) = 1;
    mat(2, 0) = -0.5;
    mat(2, 1) = 1;
    mat(2, 2) = 2;

    //Eigen::Matrix2d s = mat.sqrt();
    cout << mat.exp() << endl;

    cout << "------" << endl;
    
    //cout << s*s << endl;
    
    // dtmm::ExpmWithSqrReturnValue<Eigen::Matrix2d> ret = dtmm::expm_with_sqr(mat, mat);
    Eigen::Matrix3d m2 = mat*mat;
    Eigen::MatrixXd ret = dtmm::expm_with_sqr(mat, m2);
    cout << ret << endl;
    // cout << ret.cols() << " , " << ret.rows() << endl;

    cout << endl << " --------------------------------------------- " << endl;
}

static void test1_3()
{
    srand((unsigned int) time(0));
    using namespace Eigen;
    MatrixXd m = MatrixXd::Random(21,21);

    //Eigen::Matrix2d s = mat.sqrt();
    MatrixXd r = m.exp();

    //cout << s*s << endl;
    
    // dtmm::ExpmWithSqrReturnValue<Eigen::Matrix2d> ret = dtmm::expm_with_sqr(mat, mat);
    MatrixXd m2 = m*m;
    MatrixXd ret = dtmm::expm_with_sqr(m, m2);

    cout << (ret-r).norm() << " / " << m.norm() << endl;

    cout << endl << " --------------------------------------------- " << endl;
}

static void test2()
{
    typedef Eigen::Matrix<DComplexd, 2, 2 > DCMat;
    DCMat mat;

    mat(0, 0) = DComplexd(1,2, 4, 1);
    mat(0, 1) = DComplexd(4,1, 2,-7);
    mat(1, 0) = DComplexd(5,-3,3,-2);
    mat(1, 1) = DComplexd(-3,7,4,-5);

    //Eigen::Matrix<DComplexd, Eigen::Dynamic, Eigen::Dynamic> ret = dtmm::expm_with_sqr(mat, mat);

    complex<double> D(1,2);
    cout << std::norm(D) << "   " << std::abs(D) << endl;
    //cout << sizeof(DCMat::Scalar) << endl;
    //cout << sizeof(Eigen::Matrix2d::Scalar) << endl;
    const double l1norm = mat.cwiseAbs().colwise().sum().maxCoeff();
    cout << l1norm << endl;
    static_assert( dtmm::is_dcomplex_type<DCMat::Scalar>::value, "Is Dcomplex Type" );

    cout << endl << " --------------------------------------------- " << endl;
}

static void test2_1()
{
    typedef complex<double> CompD;

    Eigen::Matrix< CompD, 2, 2 > m1;
    m1(0, 0) = CompD(1,2)/10.;
    m1(0, 1) = CompD(4,1)/10.;
    m1(1, 0) = CompD(5,-3)/10.;
    m1(1, 1) = CompD(-3,7)/10.;
    cout << m1.exp() << endl;

    cout << " ------ " << endl;

    typedef Eigen::Matrix<DComplexd, 2, 2 > DCMat;
    DCMat mat;
    mat(0, 0) = DComplexd(1,2,0,0)/10.;
    mat(0, 1) = DComplexd(4,1,0,0)/10.;
    mat(1, 0) = DComplexd(5,-3,0,0)/10.;
    mat(1, 1) = DComplexd(-3,7,0,0)/10.;
    DCMat mat2 = mat*mat;
    DCMat ret = dtmm::expm_with_sqr(mat, mat2);
    cout << ret << endl;
    
    cout << endl << " --------------------------------------------- " << endl;
}

static void test3()
{
    cout << " --- test D[exp(A)] ---" << endl;
    using namespace Eigen;

    typedef complex<double> CompD;

    Matrix< CompD, 2, 2 > m1;
    m1(0, 0) = CompD(1,2)/10.;
    m1(0, 1) = CompD(4,1)/10.;
    m1(1, 0) = CompD(5,-3)/10.;
    m1(1, 1) = CompD(-3,7)/10.;
    
    const double h0 = 1E-7;
    Matrix< CompD, 2, 2 > delta;
    delta(0, 0) = CompD(1,0);
    delta(0, 1) = CompD(0,0);
    delta(1, 0) = CompD(0,1);
    delta(1, 1) = CompD(1,1);

    cout << ((m1+delta*h0).exp() - (m1-delta*h0).exp())/(h0*2) << endl;
    cout << " ------ " << endl;

    typedef Matrix<DComplexd, 2, 2 > DCMat;

    const double h = 1E-10;
    DCMat mat;
    mat(0, 0) = DComplexd( 0.1, 0.2, h,0);
    mat(0, 1) = DComplexd( 0.4, 0.1, 0,0);
    mat(1, 0) = DComplexd( 0.5,-0.3, 0,h);
    mat(1, 1) = DComplexd(-0.3, 0.7, h,h);
    DCMat mat2 = mat*mat;
    DCMat ret = dtmm::expm_with_sqr(mat, mat2);

    Matrix< CompD, 2, 2 > rr;
    for(int i = 0;i < 3;++ i)
        for(int j = 0;j < 3; ++ j)
            rr(i,j) = ret(i,j).imag_ / h;
    cout << rr << endl;

    cout << endl << " --------------------------------------------- " << endl;

}

static void test4()
{
    cout << " --- test4 D[exp(A)] ---" << endl;
    using namespace Eigen;

    typedef complex<double> CompD;

    Matrix< CompD, 2, 2 > m1;
    m1(0, 0) = CompD(1,2)/10.;
    m1(0, 1) = CompD(4,1)/10.;
    m1(1, 0) = CompD(5,-3)/10.;
    m1(1, 1) = CompD(-3,7)/10.;
    
    const double h0 = 1E-7;
    Matrix< CompD, 2, 2 > delta;
    delta(0, 0) = CompD(1,0);
    delta(0, 1) = CompD(0,0);
    delta(1, 0) = CompD(0,1);
    delta(1, 1) = CompD(1,1);

    cout << ((m1+delta*h0).exp() - (m1-delta*h0).exp())/(h0*2) << endl;
    cout << " ------ " << endl;

    Matrix< CompD, 4, 4 > m2;
    m2.setZero();
    m2.topLeftCorner(2,2) = m1;
    m2.topRightCorner(2,2) = delta;
    m2.bottomRightCorner(2,2) = m1;

    Matrix< CompD, 4, 4 > ret = m2.exp();
    cout << ret.topRightCorner(2,2) << endl;
}

static void test5()
{
    cout << " --- test4 D[exp(A)] ---" << endl;
    using namespace Eigen;

    typedef complex<double> CompD;

    Matrix< CompD, 2, 2 > m1;
    m1(0, 0) = CompD(1,2)/10.;
    m1(0, 1) = CompD(4,1)/10.;
    m1(1, 0) = CompD(5,-3)/10.;
    m1(1, 1) = CompD(-3,7)/10.;
    
    const double h0 = 1E-7;
    Matrix< CompD, 2, 2 > delta;
    delta(0, 0) = CompD(1,0);
    delta(0, 1) = CompD(0,0);
    delta(1, 0) = CompD(0,1);
    delta(1, 1) = CompD(1,1);
    //cout << delta.cwiseAbs().colwise().sum() << endl;

    cout << ((m1+delta*h0).exp() - (m1-delta*h0).exp())/(h0*2) << endl;
    cout << " ------ " << endl;

    Matrix< CompD, 2, 2 > deriv = mat_exp_deriv(m1, delta);
    cout << deriv << endl;
}

int main()
{
    test1();
    test1_1();
    test1_2();
    test1_3();
    test2();
    test2_1();
    test3();
    test4();
    test5();
    return 0;
}
