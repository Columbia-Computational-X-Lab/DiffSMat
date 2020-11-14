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

#ifndef DCOMPLEX_INC
#   define DCOMPLEX_INC

#include <type_traits>
#include <complex>

namespace dtmm {

template <typename T>
class DComplex
{
    public:
        typedef T   Scalar;

        std::complex<T> real_;
        std::complex<T> imag_;

        // ========== Constructors ===========
        DComplex():real_((T)0), imag_((T)0) { }
        DComplex(T a):real_(a, (T)0), imag_((T)0, (T)0) { }
        DComplex(T a, T b, T c, T d):real_(a, b), imag_(c, d) { }
        DComplex(T a, T b):real_(a, (T)0), imag_(b, (T)0) { }
        DComplex(const std::complex<T>& r, const std::complex<T>& i):real_(r), imag_(i) { }
        DComplex(const std::complex<T>& r, T i):real_(r), imag_(i) { }

        void set(const std::complex<T>& r, const std::complex<T>& i)
        {
            real_ = r;
            imag_ = i;
        }

        /*! Copy operator */
        DComplex<T>& operator = (const DComplex<T>& rhs)
        {
            real_ = rhs.real_;
            imag_ = rhs.imag_;
            return *this;
        }
        DComplex<T>& operator = (T a)
        {
            real_ = std::complex<T>(a, (T)0);
            imag_ = std::complex<T>((T)0, (T)0);
            return *this;
        }

        /*! Addition operator */
        DComplex<T> operator + (const DComplex<T>& rhs) const 
        {
            return DComplex<T> (real_ + rhs.real_, imag_ + rhs.imag_);
        }

        friend DComplex<T> operator + (T lhs, const DComplex<T>& rhs) 
        { 
            return DComplex<T>(lhs + rhs.real_, rhs.imag_);
        }

        DComplex<T>& operator += (const DComplex<T>& rhs)
        {
            real_ += rhs.real_;
            imag_ += rhs.imag_;
            return *this;
        }

        DComplex<T>& operator += (const T rhs)
        {
            real_.real(real_.real() + rhs);
            return *this;
        }
       
        /*! Substraction operator */
        DComplex<T> operator - (const DComplex<T>& rhs) const 
        {
            return DComplex<T> (real_ - rhs.real_, imag_ - rhs.imag_);
        }

        friend DComplex<T> operator - (T lhs, const DComplex<T>& rhs) 
        { 
            return DComplex<T>(lhs - rhs.real_, rhs.imag_);
        }

        DComplex<T>& operator -= (const DComplex<T>& rhs)
        {
            real_ -= rhs.real_;
            imag_ -= rhs.imag_;
            return *this;
        }
       
        DComplex<T> operator - () const
        {
            return DComplex<T>(-real_, -imag_);
        }

        /*! Multiplication operator */
/*
 * When FAST_DCOMPLEX_PRODUCT is defined, we ignore all the terms that involve
 * high-order imaginary component to save a bit computation.
 */
#ifdef FAST_DCOMPLEX_PRODUCT
        DComplex<T> operator * (const DComplex<T>& rhs) const
        {
            return DComplex<T>(
                real_*rhs.real_, //// - imag_*rhs.imag_,
                real_*rhs.imag_ + imag_*rhs.real_);
        }

        DComplex<T>& operator *= (const DComplex<T>& rhs) 
        {
            const std::complex<T> rr = real_*rhs.real_; //// - imag_*rhs.imag_;
            const std::complex<T> ii = real_*rhs.imag_ + imag_*rhs.real_;

            real_ = rr;
            imag_ = ii;
            return *this;
        }
#else
        DComplex<T> operator * (const DComplex<T>& rhs) const
        {
            return DComplex<T>(
                real_*rhs.real_ - imag_*rhs.imag_,
                real_*rhs.imag_ + imag_*rhs.real_);
        }

        DComplex<T>& operator *= (const DComplex<T>& rhs) 
        {
            const std::complex<T> rr = real_*rhs.real_ - imag_*rhs.imag_;
            const std::complex<T> ii = real_*rhs.imag_ + imag_*rhs.real_;

            real_ = rr;
            imag_ = ii;
            return *this;
        }
#endif   // --- end FAST_DCOMPLEX_PRODUCT ---

        DComplex<T> operator * (T rhs) const 
        {
            return DComplex<T>(real_*rhs, imag_*rhs);
        }

        friend DComplex<T> operator * (T lhs, const DComplex<T>& rhs) 
        { 
            return DComplex<T>(rhs.real_ * lhs, rhs.imag_ * lhs);
        }

        DComplex<T>& operator *= (T rhs)
        {
            real_ *= rhs;
            imag_ *= rhs;
            return *this;
        }

        DComplex<T> operator * (const std::complex<T>& rhs) const 
        {
            return DComplex<T>(real_*rhs, imag_*rhs);
        }

        DComplex<T>& operator *= (const std::complex<T>& rhs)
        {
            real_ *= rhs;
            imag_ *= rhs;
            return *this;
        }

        bool operator != (const DComplex<T>& rhs) const
        {
            return (real_ != rhs.real_ || imag_ != rhs.imag_);
        }

        /*! Division operator */
#ifdef FAST_DCOMPLEX_PRODUCT
        DComplex<T> operator / (const DComplex<T>& rhs) const 
        {
            const std::complex<T> s = (T)1 / (rhs.real_*rhs.real_); //// + rhs.imag_*rhs.imag_);
            return DComplex<T>(
                real_ * rhs.real_ * s,
                (imag_*rhs.real_ - real_*rhs.imag_) * s);
        }

        DComplex<T>& operator /= (const DComplex<T>& rhs) 
        {
            const std::complex<T> s = (T)1 / (rhs.real_*rhs.real_); //// + rhs.imag_*rhs.imag_);
            const std::complex<T> rr = real_ * rhs.real_ * s;
            const std::complex<T> ii = (imag_*rhs.real_ - real_*rhs.imag_) * s;

            real_ = rr;
            imag_ = ii;
            return *this;
        }
#else
        DComplex<T> operator / (const DComplex<T>& rhs) const 
        {
            const std::complex<T> s = (T)1 / (rhs.real_*rhs.real_ + rhs.imag_*rhs.imag_);
            return DComplex<T>(
                (real_*rhs.real_ + imag_*rhs.imag_) * s,
                (imag_*rhs.real_ - real_*rhs.imag_) * s);
        }

        DComplex<T>& operator /= (const DComplex<T>& rhs) 
        {
            const std::complex<T> s = (T)1 / (rhs.real_*rhs.real_ + rhs.imag_*rhs.imag_);
            const std::complex<T> rr = (real_*rhs.real_ + imag_*rhs.imag_) * s;
            const std::complex<T> ii = (imag_*rhs.real_ - real_*rhs.imag_) * s;

            real_ = rr;
            imag_ = ii;
            return *this;
        }
#endif   // --- end FAST_DCOMPLEX_PRODUCT ---

        DComplex<T> operator / (T rhs) const 
        {
            return DComplex<T>(real_/rhs, imag_/rhs);
        }

        DComplex<T>& operator /= (T rhs)
        {
            real_ /= rhs;
            imag_ /= rhs;
            return *this;
        }

        DComplex<T> operator / (const std::complex<T>& rhs) const 
        {
            const std::complex<T> s = std::complex<T>((T)1,(T)0) / rhs;
            return DComplex<T>(real_ * s, imag_ * s);
        }


        DComplex<T>& operator /= (const std::complex<T>& rhs)
        {
            const std::complex<T> s = std::complex<T>((T)1,(T)0) / rhs;
            real_ *= s;
            imag_ *= s;
            return *this;
        }

        // ----------------------------------------------------------
        
        DComplex<T> conj() const
        {   
            return DComplex<T>(real_, -imag_);
        }

        const std::complex<T>& real() const
        {   return real_; }

        const std::complex<T>& imag() const
        {   return imag_; }

        //============== stream operator =================
        /*! Output to stream operator */
        friend std::ostream& operator<<(std::ostream& lhs, const DComplex<T>& rhs) 
        {
            lhs << '[' << rhs.real_ << ", " << rhs.imag_ << ']';
            return lhs;
        }
};

typedef class DComplex<double>      DComplexd;
typedef class DComplex<float>       DComplexf;

template<typename T> struct is_dcomplex_type : std::false_type {};
template<> struct is_dcomplex_type<DComplexd> : std::true_type {};
template<> struct is_dcomplex_type<DComplexf> : std::true_type {};

}

#endif
