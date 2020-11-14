#include <iostream>
#include <fstream>

#include "utils/timer.hpp"
#include "rcwa/WaveEquationCoeff.h"
#include "core/RCWAScatterMatrix.h"
#include "core/ScatterMatrixDerivative.hpp"

using namespace std;

constexpr scalex e_Si = scalex(3.48*3.48, 0.);
constexpr scalex e_SiO2 = scalex(1.445*1.445, 0.);

void test1()
{
    const int N = 7;
    const scalar lambda = 1.55;
    const scalar wx = 1.0;
    const scalar wy = 0.22;
    const scalar wz = 0.31;

    Eigen::MatrixXcs P, Q, dP, dQ;

    WaveEquationCoeff coeff(N, N, 5., 5.);
    dtmm::RCWAScatterMatrix rcwa(N, N, 5., 5.);

    coeff.solve(lambda, wx, wy, P, Q, &dP, nullptr, &dQ, nullptr);
    rcwa.compute(lambda, wz, P, Q);

    dtmm::ScatterMatrixDerivative<dtmm::RCWAScatterMatrix> deriv;
    deriv.compute(rcwa, P, Q, dP, dQ);

    // check Omega'
    Eigen::MatrixXcs dOmg = deriv.dOmega();
    cout << dOmg.norm() << endl;

    // FD estimation
    const double H = 1E-6;
    coeff.solve(lambda, wx+H, wy, P, Q, &dP, nullptr, &dQ, nullptr);
    rcwa.compute(lambda, wz, P, Q);
    deriv.compute(rcwa, P, Q, dP, dQ);
    Eigen::MatrixXcs omgP = deriv.Omega();

    coeff.solve(lambda, wx-H, wy, P, Q, &dP, nullptr, &dQ, nullptr);
    rcwa.compute(lambda, wz, P, Q);
    deriv.compute(rcwa, P, Q, dP, dQ);
    Eigen::MatrixXcs omgM = deriv.Omega();

    cout << "Omega Delta = " << (dOmg - (omgP-omgM)/(H*2)).cwiseAbs().maxCoeff() << endl;
    cout << endl << " --------------------------------------------- " << endl;
}

/*
void test2()
{
    const int N = 7;
    const scalar lambda = 1.55;
    const scalar wx = 1.0;
    const scalar wy = 0.22;
    const scalar wz = 0.31;

    Eigen::MatrixXcs P, Q, dP, dQ;

    WaveEquationCoeff coeff(N, N, 5., 5.);
    dtmm::RCWAScatterMatrix rcwa(N, N, 5., 5.);

    coeff.solve(lambda, wx, wy, P, Q, &dP, nullptr, &dQ, nullptr);
    rcwa.compute(lambda, wz, P, Q);

    dtmm::ScatterMatrixDerivative<dtmm::RCWAScatterMatrix> deriv;
    deriv.compute(rcwa, P, Q, dP, dQ);

    // FD estimation
    const double H = 1E-6;
    coeff.solve(lambda, wx+H, wy, P, Q, &dP, nullptr, &dQ, nullptr);
    rcwa.compute(lambda, wz, P, Q);
    Eigen::MatrixXcs abP = rcwa.invAB();

    coeff.solve(lambda, wx-H, wy, P, Q, &dP, nullptr, &dQ, nullptr);
    rcwa.compute(lambda, wz, P, Q);
    Eigen::MatrixXcs abM = rcwa.invAB();

    cout << "InvAB Delta = " << (deriv.dInvAB() - (abP-abM)/(H*2)).cwiseAbs().maxCoeff() << endl;
    cout << endl << " --------------------------------------------- " << endl;
}
*/

void test3()
{
    const int N = 7;
    const scalar lambda = 1.55;
    const scalar wx = 1.0;
    const scalar wy = 0.22;
    const scalar wz = 0.31;

    Eigen::MatrixXcs P, Q, dP, dQ;

    WaveEquationCoeff coeff(N, N, 5., 5.);
    dtmm::RCWAScatterMatrix rcwa(N, N, 5., 5.);
    coeff.solve(lambda, wx, wy, P, Q, &dP, nullptr, &dQ, nullptr);
    rcwa.compute(lambda, wz, P, Q);

    dtmm::ScatterMatrixDerivative<dtmm::RCWAScatterMatrix> deriv;
    deriv.compute(rcwa, P, Q, dP, dQ);

    // check Omega'
    Eigen::MatrixXcs dv = deriv.dWXinvW();
    cout << dv.norm() << endl;

    // FD estimation
    const double H = 1E-6;
    coeff.solve(lambda, wx+H, wy, P, Q, &dP, nullptr, &dQ, nullptr);
    rcwa.compute(lambda, wz, P, Q);
    deriv.compute(rcwa, P, Q, dP, dQ);
    Eigen::MatrixXcs vP = deriv.WXinvW();

    coeff.solve(lambda, wx-H, wy, P, Q, &dP, nullptr, &dQ, nullptr);
    rcwa.compute(lambda, wz, P, Q);
    deriv.compute(rcwa, P, Q, dP, dQ);
    Eigen::MatrixXcs vM = deriv.WXinvW();

    cout << "WXW^{-1} Delta = " << (dv - (vP-vM)/(H*2)).cwiseAbs().maxCoeff() << endl;
    cout << endl << " --------------------------------------------- " << endl;
}

#if 0
void test4()
{
    const int N = 9;
    const scalar lambda = 1.55;
    const scalar wx = 1.0;
    const scalar wy = 0.22;
    const scalar wz = 0.31;

    Eigen::MatrixXcs P, Q, dP, dQ;

    WaveEquationCoeff coeff(N, N, 5., 5.);
    dtmm::RCWAScatterMatrix rcwa(N, N, 5., 5.);

    coeff.solve(lambda, wx, wy, P, Q, &dP, nullptr, &dQ, nullptr);
    rcwa.compute(lambda, wz, P, Q);

    dtmm::ScatterMatrixDerivative<dtmm::RCWAScatterMatrix> deriv;
    deriv.compute(rcwa, P, Q, dP, dQ);

    // FD estimation
    const double H = 1E-6;
    coeff.solve(lambda, wx+H, wy, P, Q, &dP, nullptr, &dQ, nullptr);
    rcwa.compute(lambda, wz, P, Q);
    Eigen::MatrixXcs abP = rcwa.ret_;

    coeff.solve(lambda, wx-H, wy, P, Q, &dP, nullptr, &dQ, nullptr);
    rcwa.compute(lambda, wz, P, Q);
    Eigen::MatrixXcs abM = rcwa.ret_;

    cout << "InvAXB Delta = " << (deriv.ret_ - (abP-abM)/(H*2)).cwiseAbs().maxCoeff() << endl;
    cout << endl << " --------------------------------------------- " << endl;
}
#endif

void test5()
{
    const int N = 23;
    const scalar lambda = 1.55;
    const scalar wx = 1.0;
    const scalar wy = 0.22;
    const scalar wz = 0.31;

    Eigen::MatrixXcs P, Q, dP, dQ;

    auto t1 = sploosh::now();
    WaveEquationCoeff coeff(N, N, 5., 5.);
    dtmm::RCWAScatterMatrix rcwa(N, N, 5., 5.);

    coeff.solve(lambda, wx, wy, P, Q, &dP, nullptr, &dQ, nullptr);
    rcwa.compute(lambda, wz, P, Q);
    auto t2 = sploosh::now();
    cout << "T1: " << sploosh::duration_milli_d(t1, t2) << endl;

    dtmm::ScatterMatrixDerivative<dtmm::RCWAScatterMatrix> deriv;
    auto t3 = sploosh::now();
    deriv.compute(rcwa, P, Q, dP, dQ);
    auto t4 = sploosh::now();
    cout << "T2: " << sploosh::duration_milli_d(t3, t4) << endl;

    // FD estimation
    const double H = 1E-6;
    coeff.solve(lambda, wx+H, wy, P, Q, &dP, nullptr, &dQ, nullptr);
    rcwa.compute(lambda, wz, P, Q);
    Eigen::MatrixXcs abP = rcwa.Tuu();

    coeff.solve(lambda, wx-H, wy, P, Q, &dP, nullptr, &dQ, nullptr);
    rcwa.compute(lambda, wz, P, Q);
    Eigen::MatrixXcs abM = rcwa.Tuu();

    cout << "Scattering Matrix Delta = " << (deriv.dTuu() - (abP-abM)/(H*2)).cwiseAbs().maxCoeff() << endl;
    cout << endl << " --------------------------------------------- " << endl;
}

void test6()
{
    const int N = 23;
    const scalar lambda = 1.55;
    const scalar wx = 1.0;
    const scalar wy = 0.22;
    const scalar wz = 0.31;

    Eigen::MatrixXcs P, Q, dP, dQ;

    auto t1 = sploosh::now();
    WaveEquationCoeff coeff(N, N, 5., 5.);
    dtmm::RCWAScatterMatrix rcwa(N, N, 5., 5.);

    coeff.solve(lambda, wx, wy, P, Q, &dP, nullptr, &dQ, nullptr);
    rcwa.compute(lambda, wz, P, Q);
    auto t2 = sploosh::now();
    cout << "T1: " << sploosh::duration_milli_d(t1, t2) << endl;

    //dtmm::ScatterMatrixDerivative<dtmm::RCWAScatterMatrix> deriv;
    auto t3 = sploosh::now();
    rcwa.compute_dz(lambda, wz);
    //deriv.compute(rcwa, P, Q, dP, dQ);
    auto t4 = sploosh::now();
    cout << "T2: " << sploosh::duration_milli_d(t3, t4) << endl;

    // FD estimation
    const double H = 1E-6;
    coeff.solve(lambda, wx, wy, P, Q, &dP, nullptr, &dQ, nullptr);
    rcwa.compute(lambda, wz+H, P, Q);
    Eigen::MatrixXcs abP = rcwa.Tuu();

    coeff.solve(lambda, wx, wy, P, Q, &dP, nullptr, &dQ, nullptr);
    rcwa.compute(lambda, wz-H, P, Q);
    Eigen::MatrixXcs abM = rcwa.Tuu();

    cout << "Scattering Matrix Delta = " << (rcwa.dTuuz() - (abP-abM)/(H*2)).cwiseAbs().maxCoeff() << endl;
    cout << endl << " --------------------------------------------- " << endl;
}

int main()
{
    test1();
    //test2();
    test3();
    //test4();
    test5();
    test6();
    return 0;
}
