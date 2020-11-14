#include <iostream>
#include <Eigen/Dense>
#include "dcomplex/DComplexMatrixFunc.hpp"
#include "dcomplex/EigenNumTraits.hpp"
#include "dcomplex/DComplexMath.hpp"

using namespace std;

void test1()
{
    using namespace dtmm;

    //typedef complex<double> CompD;
    typedef Eigen::Matrix<DComplexd, Eigen::Dynamic, Eigen::Dynamic > DCMat;
    DCMat m1(13,8), m2(8,17), ret, ret2;

    for(int i = 0;i < 13;++ i)
    for(int j = 0;j < 8;++ j)
        m1(i,j).set( std::complex<double>(i,0), std::complex<double>(0,j) );

    for(int i = 0;i < 8;++ i)
    for(int j = 0;j < 17;++ j)
        m2(i,j).set( std::complex<double>(i,0), std::complex<double>(0,j) );

    ret = m1 * m2;
    dtmm::dcomplex_mat_product(m1, m2, ret2);
    for(int i = 0;i < ret.rows();++ i)
        for(int j = 0;j < ret.cols();++ j)
            if (dtmm::abs2(ret(i,j) - ret2(i,j))>1E-10)
            {
                cerr << i << ' ' << j << ":   " << ret(i,j) << " ~ " << ret2(i,j) << endl;
            }
}

void test2()
{
    using namespace dtmm;

    //typedef complex<double> CompD;
    typedef Eigen::Matrix<DComplexd, Eigen::Dynamic, Eigen::Dynamic > DCMat;
    DCMat m1 = DCMat::Random(110,110); 
    DCMat m2 = DCMat::Random(110,110); 
    DCMat ret, ret2;

    ret = m1 * m2;
    dtmm::dcomplex_mat_product(m1, m2, ret2);
    for(int i = 0;i < ret.rows();++ i)
        for(int j = 0;j < ret.cols();++ j)
            if (dtmm::abs2(ret(i,j) - ret2(i,j))>1E-10)
            {
                cerr << i << ' ' << j << ":   " << ret(i,j) << " ~ " << ret2(i,j) << endl;
            }
}

void test3()
{
    using namespace dtmm;

    //typedef complex<double> CompD;
    typedef Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> CMat;
    typedef Eigen::Matrix<DComplexd, Eigen::Dynamic, Eigen::Dynamic > DCMat;
    DCMat m1 = dtmm::dcomplex_mat_from_pair(CMat::Random(110,10), CMat::Random(110,10), complex<double>(1,0)); 
    DCMat m2 = dtmm::dcomplex_mat_from_pair(CMat::Random(10,100), CMat::Random(10,100), complex<double>(1,0)); 
    DCMat ret, ret2;
    /*
    DCMat m1(110,110); DCMat m2(110,110);
    DCMat ret, ret2;
    for(int i = 0;i < m1.rows();++ i)
    for(int j = 0;j < m1.cols();++ j)
    {
        m1(i,j) = DComplexd(rand(), rand(), rand(), rand());
        m1(i,j) /= RAND_MAX;
        m2(i,j) = DComplexd(rand(), rand(), rand(), rand());
        m2(i,j) /= RAND_MAX;
    }
    */

    ret = m1 * m2;
    dtmm::dcomplex_mat_product(m1, m2, ret2);
    for(int i = 0;i < ret.rows();++ i)
        for(int j = 0;j < ret.cols();++ j)
            if (dtmm::abs2(ret(i,j) - ret2(i,j))>1E-10)
            {
                cerr << i << ' ' << j << ":   " << ret(i,j) << " ~ " << ret2(i,j) << endl;
            }
}

int main()
{
    test3();
    return 0;
}

