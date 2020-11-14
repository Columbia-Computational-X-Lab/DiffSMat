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
#include <iomanip>
#include "dcomplex/DComplexMath.hpp"

using namespace std;
using namespace dtmm;

static void test1()
{
    DComplexd dd(1,2,3,5);
    cout << dd << endl;

    DComplexd one(1,0,0,0);
    cout << (one/dd) << "    " << (one/dd)* dd << endl;

    DComplexd t2(0,1,0,0);
    cout << (t2/dd)* dd << endl;

    DComplexd t3(0,0,1,0);
    cout << (t3/dd)* dd << endl;

    DComplexd t4(0,0,0,1);
    cout << (t4/dd)* dd << endl;

    cout << endl << " --------------------------------- " << endl;
}

static void test2()
{
    complex<double> a(1, 2);
    cout << "--- D[z*z] --- " << endl;
    cout << a*2.0 << endl;
    double h = 1E-10;
    DComplexd b(a, h);
    DComplexd tt = b*b;
    cout << tt.imag()/h << endl;

    cout << "--- D[1/z] --- " << endl;
    cout << complex<double>(-1.,0)/(a*a) << endl;
    tt = DComplexd(1.)/b;
    cout << tt.imag()/h << endl;

    cout << endl << " --------------------------------- " << endl;
}

static void test3()
{
    cout << "--- D[exp(2*z)] --- " << endl;
    complex<double> a(1, 3);
    cout << exp(a*2.)*2. << endl;
    double h = 1E-10;
    DComplexd b(a, h);
    DComplexd tt = std::exp(b*2.);
    cout << tt.imag_/h << endl;

    cout << endl << " --------------------------------- " << endl;
}

static void test4()
{
    complex<double> a(0.1, 0.3);
    cout << "--- D[sin(z*z)] --- " << endl;
    cout << cos(a*a) * a * 2. << endl;

    double h = 1E-12;
    DComplexd b(a, h);
    DComplexd tt = sin(b*b);
    cout << tt.imag_/h << endl;

    cout << endl << "--- D[exp(cos(z*z)))] --- " << endl;
    cout << exp( cos(a*a) ) * sin(a*a) * a * (-2.) << endl;

    DComplexd t2 = std::exp( cos(b*b) );
    cout << t2.imag_/h << endl;
    cout << endl << " --------------------------------- " << endl;
}

int main()
{
    test1();
    cout << setprecision(16);
    test2();
    test3();
    test4();
    return 0;
}
