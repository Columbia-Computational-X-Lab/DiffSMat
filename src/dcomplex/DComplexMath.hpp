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

#ifndef DCOMPLEX_MATH_INC
#   define DCOMPLEX_MATH_INC

#include "DComplex.hpp"
#include <cmath>

namespace dtmm {

template <typename T>
inline DComplex<T> conj(const DComplex<T>& x)
{
    return x.conj();
}

template <typename T>
inline std::complex<T> real(const DComplex<T>& x)
{
    return x.real_;
}

template <typename T>
inline std::complex<T> imag(const DComplex<T>& x)
{
    return x.imag_;
}

template <typename T>
DComplex<T> exp(const DComplex<T>& x)
{
    return  DComplex<T>(std::cos(x.imag_), std::sin(x.imag_)) * std::exp(x.real_);
}

/*
 * sin(𝑎+𝑏𝑖)=sin𝑎cosh𝑏+𝑖cos𝑎sinh𝑏
 */
template <typename T>
DComplex<T> sin(const DComplex<T>& x)
{
    return DComplex<T>(
            std::sin(x.real_) * std::cosh(x.imag_),
            std::cos(x.real_) * std::sinh(x.imag_));
}

/*
 * cos(𝑎+𝑏𝑖)=cos𝑎cosh𝑏−𝑖sin𝑎sinh𝑏
 */
template <typename T>
DComplex<T> cos(const DComplex<T>& x)
{
    return DComplex<T>(
            std::cos(x.real_) * std::cosh(x.imag_),
           -std::sin(x.real_) * std::sinh(x.imag_));
}

template <typename T>
T abs(const DComplex<T>& x)
{
#ifdef FAST_DCOMPLEX_PRODUCT
    return std::sqrt(std::norm(x.real_));
#else
    return std::sqrt(std::norm(x.real_) + std::norm(x.imag_));
#endif   // --- end FAST_DCOMPLEX_PRODUCT ---
}

template <typename T>
T abs2(const DComplex<T>& x)
{
#ifdef FAST_DCOMPLEX_PRODUCT
    return std::norm(x.real_);
#else
    return std::norm(x.real_) + std::norm(x.imag_);
#endif   // --- end FAST_DCOMPLEX_PRODUCT ---
}

} // end dtmm

namespace std {

template <typename T>
inline dtmm::DComplex<T> exp(const dtmm::DComplex<T>& x)
{
    return dtmm::exp(x);
}

}

#endif
