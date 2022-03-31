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

#ifndef DCOMPLEX_MAT_FUNC_INC
#   define DCOMPLEX_MAT_FUNC_INC

#include <tbb/tbb.h>
#include <Eigen/Core>
#include "DComplex.hpp"
#include "utils/macros.h"

namespace dtmm {

namespace internal {

template<typename T> struct is_complex: std::false_type {};
template<> struct is_complex< std::complex<double> > : std::true_type {};
template<> struct is_complex< std::complex<float> > : std::true_type {};
}

template <typename Derived>
class DComplexMatFromPairReturnByValue : public Eigen::ReturnByValue< DComplexMatFromPairReturnByValue<Derived> >
{
    typedef typename Eigen::NumTraits<typename Derived::Scalar>::Real  Real;
    typedef typename Derived::Scalar    ComplexType;
    typedef typename Derived::Index     Index;

    public:
        DComplexMatFromPairReturnByValue(const Derived& m1, const Derived& m2, const ComplexType& s, Real h):
                m1_(m1), m2_(m2), s_(s), h_(h) {}

        Index rows() const { return m1_.rows(); }
        Index cols() const { return m1_.cols(); }

        template <typename ResultType>
        inline void evalTo(ResultType& result) const
        {
            result.resize(rows(), cols());
            tbb::parallel_for( tbb::blocked_range2d<Eigen::Index>(0, rows(), 0, cols()),
                    [&](const tbb::blocked_range2d<Eigen::Index>& r) {
                        for(Eigen::Index i = r.rows().begin();i != r.rows().end(); ++ i) {
                        for(Eigen::Index j = r.cols().begin();j != r.cols().end(); ++ j) {
                            result(i,j).set( m1_(i,j)*s_, m2_(i,j)*(s_*h_) );
                        } }
                    }
            );
        }

    private:
        const Derived& m1_;
        const Derived& m2_;
        const ComplexType&  s_;
        const Real          h_;
};

} // end of namespace dtmm

namespace Eigen {
namespace internal {

/** to support ReturnByValue used in Eigen */
template<typename Derived>
struct traits< dtmm::DComplexMatFromPairReturnByValue<Derived> >
{
    typedef typename NumTraits<typename Derived::Scalar>::Real  Real;
    typedef Matrix< dtmm::DComplex<Real>, Dynamic, Dynamic >    ReturnType;
};

}
} // end of namespace Eigen

// ---------------------------------------------------------------------------

namespace dtmm {

template <typename Derived> 
const DComplexMatFromPairReturnByValue<Derived>
dcomplex_mat_from_pair(const Eigen::MatrixBase<Derived>& m1,
                       const Eigen::MatrixBase<Derived>& m2,
                       const typename Derived::Scalar& s,
                       const typename Eigen::NumTraits<typename Derived::Scalar>::Real h=1)
{
    static_assert( internal::is_complex<typename Derived::Scalar>::value,
                   "The input must be complex matrices" );
    assert(m1.rows() == m2.rows() && m1.cols() == m2.cols());
    return DComplexMatFromPairReturnByValue<Derived>(m1.derived(), m2.derived(), s, h);
}

template <typename ArgType, typename RetType>
void dcomplex_mat_deriv(
        const ArgType& ret, 
        const typename Eigen::NumTraits<typename ArgType::Scalar>::Real h,
        RetType& out)
{
    static_assert( is_dcomplex_type<typename ArgType::Scalar>::value,
                   "The input must be dcomplex matrices" );
    out.resize(ret.rows(), ret.cols());
    
    tbb::parallel_for( tbb::blocked_range2d<Eigen::Index>(0, ret.rows(), 0, ret.cols()),
            [&](const tbb::blocked_range2d<Eigen::Index>& r) {
                for(Eigen::Index i = r.rows().begin();i != r.rows().end(); ++ i) {
                for(Eigen::Index j = r.cols().begin();j != r.cols().end(); ++ j) {
                    out(i,j) = ret(i,j).imag() / h;
                } }
            }
    );
}

/*
 * compute ret = A*B
 */
template <typename Derived, typename ResultType>
void dcomplex_mat_product(const Eigen::MatrixBase<Derived>& A,
                          const Eigen::MatrixBase<Derived>& B,
                          ResultType& ret)
{
    typedef typename Derived::Scalar    DComplexType;
    static_assert( is_dcomplex_type<typename Derived::Scalar>::value,
                   "The input must be dcomplex matrices" );
    assert(A.cols() == B.rows());

    const Eigen::Index M = A.rows();
    const Eigen::Index N = B.cols();
    const Eigen::Index L = A.cols();
    static tbb::affinity_partitioner ap;

    ret.resize(M, N);
    tbb::parallel_for( tbb::blocked_range2d<Eigen::Index>(0, M, 0, N),
            [&](const tbb::blocked_range2d<Eigen::Index>& r) {
                for(Eigen::Index i = r.rows().begin();i != r.rows().end(); ++ i) {
                for(Eigen::Index j = r.cols().begin();j != r.cols().end(); ++ j) {
                    DComplexType sum;
                    for(Eigen::Index k = 0;k < L;++ k)
                        sum += A(i,k) * B(k,j);
                    ret(i,j) = sum;
                } }
            }, ap
    );
}

/*
 * compute ret += a1*A1 + a2*A2 + a3*A3 + c*I
 */
template <typename T, typename Derived, typename ResultType>
void dcomplex_mat_axpy(T a1, const Eigen::MatrixBase<Derived>& A1,
                       T a2, const Eigen::MatrixBase<Derived>& A2,
                       T a3, const Eigen::MatrixBase<Derived>& A3,
                       T c, ResultType& ret)
{
    static_assert( is_dcomplex_type<typename Derived::Scalar>::value,
                   "The input must be dcomplex matrices" );
    assert(A1.rows() == A2.rows() && A1.cols() == A2.cols() &&
           A2.rows() == A3.rows() && A2.cols() == A3.cols() &&
           A3.rows() == ret.rows() && A3.cols() == ret.cols());

    tbb::parallel_for( tbb::blocked_range2d<Eigen::Index>(0, ret.rows(), 0, ret.cols()),
            [&](const tbb::blocked_range2d<Eigen::Index>& r) {
                for(Eigen::Index i = r.rows().begin();i != r.rows().end(); ++ i) {
                for(Eigen::Index j = r.cols().begin();j != r.cols().end(); ++ j) {
                    ret(i,j) += (A1(i,j)*a1 + A2(i,j)*a2 + A3(i,j)*a3);
                    if ( i == j ) ret(i,j) += c;
                } }
            }
    );
} // end

/*
 * compute ret = a1*A1 + a2*A2 + a3*A3
 */
template <typename T, typename Derived, typename ResultType>
void dcomplex_mat_sum(T a1, const Eigen::MatrixBase<Derived>& A1,
                      T a2, const Eigen::MatrixBase<Derived>& A2,
                      T a3, const Eigen::MatrixBase<Derived>& A3,
                      ResultType& ret)
{
    static_assert( is_dcomplex_type<typename Derived::Scalar>::value,
                   "The input must be dcomplex matrices" );
    assert(A1.rows() == A2.rows() && A1.cols() == A2.cols() &&
           A2.rows() == A3.rows() && A2.cols() == A3.cols()); 
    ret.resize(A1.rows(), A1.cols());

    tbb::parallel_for( tbb::blocked_range2d<Eigen::Index>(0, ret.rows(), 0, ret.cols()),
            [&](const tbb::blocked_range2d<Eigen::Index>& r) {
                for(Eigen::Index i = r.rows().begin();i != r.rows().end(); ++ i) {
                for(Eigen::Index j = r.cols().begin();j != r.cols().end(); ++ j) {
                    ret(i,j) = (A1(i,j)*a1 + A2(i,j)*a2 + A3(i,j)*a3);
                } }
            }
    );
} // end

} // end namespace dtmm

// ----------------------------------------------------------------------------------------

#if 1
namespace Eigen {
namespace internal {

/*
 * This is to allow overload Eigen's * operator for matrix<DComplexd> multiplication.
 *
 * We use TBB multi-threading to improve its performance.
 */
template<typename Index, int LhsStorageOrder, bool ConjugateLhs, int RhsStorageOrder, bool ConjugateRhs, int ResInnerStride>
struct general_matrix_matrix_product<Index, dtmm::DComplexd, LhsStorageOrder, ConjugateLhs, dtmm::DComplexd, RhsStorageOrder, ConjugateRhs, ColMajor,ResInnerStride>
{
    typedef gebp_traits<dtmm::DComplexd, dtmm::DComplexd> Traits;

    static void run(Index rows, Index cols, Index depth,
            const dtmm::DComplexd* _lhs, Index lhsStride,
            const dtmm::DComplexd* _rhs, Index rhsStride,
            dtmm::DComplexd* res, Index resIncr, Index resStride,
            dtmm::DComplexd alpha,
            level3_blocking<dtmm::DComplexd, dtmm::DComplexd>& /*blocking*/,
            GemmParallelInfo<Index>* /*info = 0*/)
    {
        UNSUPPORTED_FUNC_IMPL(-1);
    }
};

template<typename Index, int ResInnerStride>
struct general_matrix_matrix_product<Index, dtmm::DComplexd, ColMajor, false, dtmm::DComplexd, ColMajor, false, ColMajor, ResInnerStride>
{
    typedef gebp_traits<dtmm::DComplexd, dtmm::DComplexd> Traits;

    static void run(Index rows, Index cols, Index depth,
            const dtmm::DComplexd* _lhs, Index lhsStride,
            const dtmm::DComplexd* _rhs, Index rhsStride,
            dtmm::DComplexd* res, Index resIncr, Index resStride,
            dtmm::DComplexd alpha,
            level3_blocking<dtmm::DComplexd, dtmm::DComplexd>& /*blocking*/,
            GemmParallelInfo<Index>* /*info = 0*/)
    {
        static tbb::affinity_partitioner ap;

        tbb::parallel_for( tbb::blocked_range2d<Index>(0, rows, 0, cols), 
                [&](const tbb::blocked_range2d<Index>& r) {
                    for(Index i = r.rows().begin();i != r.rows().end(); ++ i) {
                    for(Index j = r.cols().begin();j != r.cols().end(); ++ j) {
                        dtmm::DComplexd sum;
                        for(Eigen::Index k = 0;k < depth;++ k)
                            sum += _lhs[k*lhsStride+i] * _rhs[j*rhsStride+k];  //A(i,k) * B(k,j);
                        res[j*resStride+i] += sum * alpha; 
                    } }
                }, ap
        );
    }
};

} } // end of namespace 
#endif

#endif
