#include "MKLEigenSolver.h"
#include <mkl.h>
#include <stdio.h>
#include <stdlib.h>

void MKLEigenSolver::compute(const Eigen::MatrixXcs& A)
{
    int n = A.rows();
    lapack_complex_double * a = nullptr;
    lapack_complex_double * w = nullptr;
    lapack_complex_double * vr = nullptr;

    a = new lapack_complex_double[n*n];
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            a[i * n + j].real = A(i, j).real();
            a[i * n + j].imag = A(i, j).imag();
        }
    }
    w = new lapack_complex_double[n];
    vr = new lapack_complex_double[n*n];

    LAPACKE_zgeev(LAPACK_ROW_MAJOR,
        'N',
        'V',
        n,
        a,
        n,
        w,
        nullptr,
        n,
        vr,
        n);

    V_.resize(n , n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            V_(i, j) = scalex(vr[i * n + j].real, vr[i * n + j].imag);
        }
    }

    d_.resize(n);
    for (int i = 0; i < n; ++i)
    {
        d_(i) = scalex(w[i].real, w[i].imag);
    }

    delete[] a;
    delete[] w;
    delete[] vr;
}
