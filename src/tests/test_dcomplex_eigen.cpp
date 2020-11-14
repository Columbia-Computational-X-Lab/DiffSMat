#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "dcomplex/DComplexMath.hpp"
#include "dcomplex/EigenNumTraits.hpp"
#include "dcomplex/DComplexMatrixFunc.hpp"

using namespace std;
using namespace dtmm;

static void test1()
{
    Eigen::Matrix< DComplexd, 3, 1 > cVec;
    for(int i = 0;i < 3;++ i)
    {
        cVec[i] = DComplexd(0, i+1, 0, (i+1)*2);
        cout << cVec[i] << endl;
    }

    cout << "--------" << endl;
    cout << cVec.sum() << endl;
    cout << cVec.mean() << endl;

    cout << endl << " --------------------------------------------- " << endl;
}

static void test2()
{
    Eigen::Matrix< DComplexd, 3, 1 > cVec;
    Eigen::Matrix< DComplexd, 3, 1 > bVec;
    for(int i = 0;i < 3;++ i)
    {
        cVec[i] = DComplexd(0, i+1, 0, (i+1)*2);
        bVec[i] = DComplexd(i, i+1, i*2, (i+1)*2);
        cout << cVec[i] << "  ,  " << bVec[i] << "  ,  " << cVec[i].conj() << " , " << cVec[i].conj()*bVec[i] << endl;
    }

    cout << cVec.dot(bVec) << endl;

    cout << endl << " --------------------------------------------- " << endl;
}

static void test3()
{
    typedef complex<double> CompD;

    cout << "--- matrix mul ---" << endl;
    Eigen::Matrix< DComplexd, 2, 2 > mat;
    mat(0, 0) = DComplexd(1,2, 4, 1);
    mat(0, 1) = DComplexd(4,1, 2,-7);
    mat(1, 0) = DComplexd(5,-3,3,-2);
    mat(1, 1) = DComplexd(-3,7,4,-5);
    cout << mat*mat << endl;
    cout << endl;

    // Test D[A*A] = A*A' + A'*A
    Eigen::Matrix< CompD, 2, 2 > m1, m2;
    m1(0, 0) = CompD(1,2);
    m1(0, 1) = CompD(4,1);
    m1(1, 0) = CompD(5,-3);
    m1(1, 1) = CompD(-3,7);

    m2(0, 0) = CompD(4,1);
    m2(0, 1) = CompD(2,-7);
    m2(1, 0) = CompD(3,-2);
    m2(1, 1) = CompD(4,-5);
    cout << m1*m2 + m2*m1 << endl;

    double h = 1E-10;
    mat(0, 0) = DComplexd(1,2, 4*h, 1*h);
    mat(0, 1) = DComplexd(4,1, 2*h,-7*h);
    mat(1, 0) = DComplexd(5,-3,3*h,-2*h);
    mat(1, 1) = DComplexd(-3,7,4*h,-5*h);
    Eigen::Matrix< DComplexd, 2, 2 > mret = mat*mat;
    cout << "----" << endl;
    cout << mret(0,0).imag_/h << " <---> " << mret(0,1).imag_/h << endl;
    cout << mret(1,0).imag_/h << " <---> " << mret(1,1).imag_/h << endl;
    
    cout << endl << " --------------------------------------------- " << endl;
}

static void test4()
{
    typedef complex<double> CompD;

    cout << " --- test LU factorization --- " << endl;
    Eigen::Matrix< CompD, 2, 2 > m1;
    m1(0, 0) = CompD(1,2);
    m1(0, 1) = CompD(4,1);
    m1(1, 0) = CompD(5,-3);
    m1(1, 1) = CompD(-3,7);
    //Eigen::ColPivHouseholderQR<Eigen::Matrix< CompD,2,2> > dec(m1);
    auto qr = m1.fullPivLu();
    Eigen::Matrix< CompD, 2, 2 > up1 = qr.matrixLU().triangularView<Eigen::Upper>();
    cout << up1 << endl;

    typedef Eigen::Matrix< DComplexd, 2, 2 > DCMat;
    DCMat mat;
    //double h = 1E-10;
    mat(0, 0) = DComplexd(1,2,0,0);
    mat(0, 1) = DComplexd(4,1,0,0);
    mat(1, 0) = DComplexd(5,-3,0,0);
    mat(1, 1) = DComplexd(-3,7,0,0);
    //auto qr2 = mat.colPivHouseholderQr();
    auto qr2 = mat.fullPivLu();
    DCMat up = qr2.matrixLU().triangularView<Eigen::Upper>();
    cout << " ---- UP TRI Mat ---- " << endl;
    cout << up << endl;
    cout << endl << " --------------------------------------------- " << endl;
}

// test 2x2 linear solve
static void test5()
{
    typedef complex<double> CompD;

    cout << " --- test 2x2 linear solve --- " << endl;
    typedef Eigen::Matrix< CompD, 2, 2 > ComplexMat;
    ComplexMat m1;
    m1(0, 0) = CompD(1,2);
    m1(0, 1) = CompD(4,1);
    m1(1, 0) = CompD(5,-3);
    m1(1, 1) = CompD(-3,7);
    ComplexMat m1Inv = m1.fullPivLu().solve(Eigen::Matrix< CompD, 2, 2 >::Identity());
    cout << m1Inv << endl;

    typedef Eigen::Matrix< DComplexd, 2, 2 > DCMat;
    DCMat mat;
    mat(0, 0) = DComplexd( 1, 2,0,0);
    mat(0, 1) = DComplexd( 4, 1,0,0);
    mat(1, 0) = DComplexd( 5,-3,0,0);
    mat(1, 1) = DComplexd(-3, 7,0,0);
    DCMat m2Inv = mat.fullPivLu().solve(DCMat::Identity());
    cout << " --- " << endl;
    cout << m2Inv << endl;

    double h0 = 1E-6;
    m1(0, 0) = CompD(1,2+h0);
    ComplexMat m1Plus  = m1.fullPivLu().solve(Eigen::Matrix< CompD, 2, 2 >::Identity());
    m1(0, 0) = CompD(1,2-h0);
    ComplexMat m1Minus = m1.fullPivLu().solve(Eigen::Matrix< CompD, 2, 2 >::Identity());
    ComplexMat derEst = (m1Plus - m1Minus)/(2.*h0);
    cout << "FD Est: " << endl << derEst << endl;

    const double h = 1E-10;
    mat(0, 0) = DComplexd( 1, 2,0,h);
    DCMat mmInv = mat.fullPivLu().solve(DCMat::Identity());
    cout << "Complex purturb Est: " << endl;
    //for(int i = 0;i < 2;++ i)
    //    for(int j = 0;j < 2; ++ j)
    //        m1(i,j) = mmInv(i,j).imag_ / h;
    //cout << m1 << endl;
    cout << mmInv(0,0).imag_/h << " , " << mmInv(0,1).imag_/h << endl;
    cout << mmInv(1,0).imag_/h << " , " << mmInv(1,1).imag_/h << endl;
    cout << endl << " --------------------------------------------- " << endl;
}

static void test6()
{
    typedef complex<double> CompD;

    cout << " --- test 3x3 linear solve --- " << endl;
    typedef Eigen::Matrix< CompD, 3, 3 > ComplexMat;
    ComplexMat m1;
    m1(0, 0) = CompD(1,2);
    m1(0, 1) = CompD(4,1);
    m1(0, 2) = CompD(2,-6);
    m1(1, 0) = CompD(5,-3);
    m1(1, 1) = CompD(-3,7);
    m1(1, 2) = CompD(4,-3);
    m1(2, 0) = CompD(4,2);
    m1(2, 1) = CompD(-6,7);
    m1(2, 2) = CompD(1, 3);

    double h0 = 1E-6;
    m1(1, 2) = CompD(4+h0,-3+h0);
    ComplexMat m1Plus  = m1.fullPivLu().solve(Eigen::Matrix< CompD, 3, 3 >::Identity());
    m1(1, 2) = CompD(4-h0,-3-h0);
    ComplexMat m1Minus = m1.fullPivLu().solve(Eigen::Matrix< CompD, 3, 3 >::Identity());
    ComplexMat derEst = (m1Plus - m1Minus)/(2.*h0);
    cout << "FD Est: " << endl << derEst << endl;

    // now estimate using complex perturbation
    typedef Eigen::Matrix< DComplexd, 3, 3 > DCMat;
    DCMat mat;
    mat(0, 0) = DComplexd( 1, 2, 0, 0);
    mat(0, 1) = DComplexd( 4, 1, 0, 0);
    mat(0, 2) = DComplexd( 2,-6, 0, 0);
    mat(1, 0) = DComplexd( 5,-3, 0, 0);
    mat(1, 1) = DComplexd(-3, 7, 0, 0);
    mat(1, 2) = DComplexd( 4,-3, 0, 0);
    mat(2, 0) = DComplexd( 4, 2, 0, 0);
    mat(2, 1) = DComplexd(-6, 7, 0, 0);
    mat(2, 2) = DComplexd( 1, 3, 0, 0);

    const double h = 1E-10;
    mat(1, 2) = DComplexd( 4,-3, h, h);
    DCMat mmInv = mat.fullPivLu().solve(DCMat::Identity());
    //DCMat mmInv = mat.partialPivLu().solve(DCMat::Identity());
    cout << "Complex purturb Est: " << endl;
    for(int i = 0;i < 3;++ i)
        for(int j = 0;j < 3; ++ j)
            m1(i,j) = mmInv(i,j).imag_ / h;
    cout << m1 << endl;
    cout << endl << " --------------------------------------------- " << endl;
}

static void test7()
{
    typedef complex<double> CompD;

    typedef Eigen::Matrix< CompD, 3, 3 > ComplexMat;
    ComplexMat m1;
    m1(0, 0) = CompD(1,2);
    m1(0, 1) = CompD(4,1);
    m1(0, 2) = CompD(2,-6);
    m1(1, 0) = CompD(5,-3);
    m1(1, 1) = CompD(-3,7);
    m1(1, 2) = CompD(4,-3);
    m1(2, 0) = CompD(4,2);
    m1(2, 1) = CompD(-6,7);
    m1(2, 2) = CompD(1, 3);

    typedef Eigen::Matrix< DComplexd, 3, 3 > DCMat;
    DCMat mat = dcomplex_mat_from_pair(m1, m1, std::complex<double>(1,0));

    cout << mat << endl;
    cout << endl << " --------------------------------------------- " << endl;
}

int main()
{
    test1();
    test2();
    test3();
    test4();
    test5();
    test6();
    test7();
    
    return 0;
}
