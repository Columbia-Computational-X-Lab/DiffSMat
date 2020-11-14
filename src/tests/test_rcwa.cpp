#include <iostream>
#include <fstream>

#include "utils/timer.hpp"
#include "rcwa/WaveEquationCoeff.h"
#include "rcwa/FourierSolver1D.h"
#include "rcwa/FourierSolver2D.h"
#include "rcwa/ScatterMatrixSolver.h"
#include "core/RCWAScatterMatrix.h"

using namespace std;
using namespace Eigen;
constexpr scalex e_Si = scalex(3.48*3.48, 0.);
constexpr scalex e_SiO2 = scalex(1.445*1.445, 0.);

void test_fourier_1d()
{
    VectorXs x(4);
    x << 0., 1.5, 2.5, 4.;
    FourierSolver1D solver(x, 4.);
    VectorXcs f(3);
    f << e_SiO2, e_Si, e_SiO2;

    int M = 21;
    VectorXcs Ff;

    solver.solve(f, M, Ff);
    cout << Ff << endl;

    VectorXcs dFf;
    solver.dSolve(f, M, 0, dFf);

    cout << "Ff(x) + 0.0001 * dFf(x)" << endl;
    cout << Ff + 0.0001 * dFf << endl;

    VectorXs x1(4);
    x1 << 0.0001, 1.5, 2.5, 4;

    FourierSolver1D solver1(x1, 4.);
    solver1.solve(f, M, Ff);

    cout << "Ff(x+0.0001)" << endl;
    cout << Ff << endl;
}

void test_fourier_2d()
{
    VectorXs x(4);
    x << 0., 1.5, 2.5, 4.;

    VectorXs y(4);
    y << 0., 1.89, 2.11, 4.;

    MatrixXcs eps;
    eps.resize(3, 3);
    eps(0, 0) = eps(0, 1) = eps(0, 2) = 
    eps(1, 0) =             eps(1, 2) =
    eps(2, 0) = eps(2, 1) = eps(2, 2) = e_SiO2;

                eps(1, 1) = e_Si;


    FourierSolver2D e11Solver(x, y);
    MatrixXcs fe11, dfe11x, dfe11y;
    e11Solver.solve(eps, 3, 3, ContinuousXY, fe11, 2, &dfe11x, 2, &dfe11y);

    ofstream fout0("fe11_0.txt");
    fout0 << fe11 << endl;

    ofstream fout("fe11.txt");
    fout << fe11 + 0.0001 * dfe11x + 0.0001 * dfe11y << endl;


    VectorXs x1(4);
    x1 << 0., 1.5, 2.5001, 4.;
    VectorXs y1(4);
    y1 << 0., 1.89, 2.1101, 4.;

    FourierSolver2D e11Solver_ref(x1, y1);
    MatrixXcs fe11ref;
    e11Solver_ref.solve(eps, 3, 3, ContinuousXY, fe11ref);
    ofstream foutref("fe11ref.txt");
    foutref << fe11ref << endl;    
}

void test_wave_equation_coeff()
{
    WaveEquationCoeff coeff(21, 21, 5.0, 5.0);

    scalar lambda = 1.55;
    scalar wx = 1.0;
    scalar wy = 0.22;

    MatrixXcs P, Q, dPx, dPy, dQx, dQy;
    coeff.solve(lambda, wx, wy, P, Q, &dPx, &dPy, &dQx, &dQy);

    ofstream fout("PQ0.txt");
    fout << P << endl;
    fout << endl;
    fout << Q << endl;

    ofstream fout1("PQ.txt");
    fout1 << P + 0.0001 * dPx << endl;
    fout1 << endl;
    fout1 << Q + 0.0001 * dQx << endl;

    ofstream foutref("PQref.txt");
    coeff.solve(lambda, wx+0.0001, wy, P, Q);
    foutref << P << endl;
    foutref << endl;
    foutref << Q << endl;
}

void test1()
{
    ScatterMatrixSolver solver(5, 5, 5.0, 5.0);

    scalar lambda = 1.55;
    scalar wx = 1.0;
    scalar wy = 0.22;
    scalar wz = 0.31;

    MatrixXcs Tuu, Rud, Rdu, Tdd;
        //auto t1 = sploosh::now();
    solver.solve(lambda, wx, wy, wz, Tuu, Rud, Rdu, Tdd);
        //auto t2 = sploosh::now();
        //cout << "TIME: " << sploosh::duration_milli_d(t1, t2) << endl;

    ofstream fout1("Tuu1.txt");
    fout1 << Tuu << endl;

    ofstream fout2("Rud1.txt");
    fout2 << Rud << endl;

    ofstream fout3("Rdu1.txt");
    fout3 << Rdu << endl;

    ofstream fout4("Tdd1.txt");
    fout4 << Tdd << endl;
}

void test2()
{
    WaveEquationCoeff coeff(9, 9, 5.0, 5.0);

    const scalar lambda = 1.55;
    const scalar wx = 1.0;
    const scalar wy = 0.42;

    MatrixXcs P, Q, dPx, dQx, dPy, dQy;
    coeff.solve(lambda, wx, wy, P, Q, &dPx, &dPy, &dQx, &dQy);

    MatrixXcs Pplus, Pminus, Qplus, Qminus;
    const scalar h = 1E-6;
    coeff.solve(lambda, wx+h, wy, Pplus, Qplus);
    coeff.solve(lambda, wx-h, wy, Pminus, Qminus);

    MatrixXcs dPEst = (Pplus - Pminus) / (2.*h);
    MatrixXcs dQEst = (Qplus - Qminus) / (2.*h);

    cout << "Delta = " << (dPEst - dPx).cwiseAbs().maxCoeff() << endl;
    cout << "Delta = " << (dQEst - dQx).cwiseAbs().maxCoeff() << endl;

    // -----------------
    
    coeff.solve(lambda, wx, wy+h, Pplus, Qplus);
    coeff.solve(lambda, wx, wy-h, Pminus, Qminus);

    dPEst = (Pplus - Pminus) / (2.*h);
    dQEst = (Qplus - Qminus) / (2.*h);

    cout << " -------------- " << endl;

    cout << "Delta = " << (dPEst - dPy).cwiseAbs().maxCoeff() << endl;
    cout << "Delta = " << (dQEst - dQy).cwiseAbs().maxCoeff() << endl;
}

void test3()
{
    cout << " -----------------------------" << endl;
    const int N = 21;
    ScatterMatrixSolver solver(N, N, 5.0, 5.0);

    scalar lambda = 1.55;
    scalar wx = 1.0;
    scalar wy = 0.22;
    scalar wz = 0.31;

    MatrixXcs Tuu, Rud, Rdu, Tdd;
    solver.solve(lambda, wx, wy, wz, Tuu, Rud, Rdu, Tdd);

    // ------------------------------------------

    WaveEquationCoeff coeff(N, N, 5., 5.);
    MatrixXcs P, Q;
    coeff.solve(lambda, wx, wy, P, Q);
    dtmm::RCWAScatterMatrix rcwa(N, N, 5., 5.);
    rcwa.compute(lambda, wz, P, Q);

    cout << "Delta: " << (rcwa.Tuu() - Tuu).cwiseAbs().maxCoeff() << endl;
    cout << "Delta: " << (rcwa.Rud() - Rud).cwiseAbs().maxCoeff() << endl;
}

int main()
{
    test1();
    test3();
    return 0;
}
