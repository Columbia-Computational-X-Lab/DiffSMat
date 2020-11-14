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

/***********************************************************************
 * File              : EigenNumTraits.hpp
 * Author            : Changxi Zheng <cxz@cs.columbia.edu>
 * Date              : 01/07/2020
 * Last Modified Date: 01/07/2020
 * Last Modified By  : Changxi Zheng <cxz@cs.columbia.edu>
 ***********************************************************************/
#ifndef EIGEN_NUM_TRAITS_INC
#   define EIGEN_NUM_TRAITS_INC

#include "DComplexMath.hpp"

/*
 * Define NumTraits to allow DComplex to be used in Eigen
 */
namespace Eigen {

template <typename _Real>
struct NumTraits< dtmm::DComplex<_Real> > 
    : GenericNumTraits< dtmm::DComplex<_Real> >
{
    typedef _Real           Real;
    // type for DComplex operations resulting in non-integer values
    typedef dtmm::DComplex<_Real> NonInteger;
    typedef typename NumTraits<_Real>::Literal Literal;

    enum {
        IsComplex = 1,
        IsInteger = 0,
        IsSigned  = NumTraits<_Real>::IsSigned,
        RequireInitialization = NumTraits<_Real>::RequireInitialization,
        ReadCost = 2 * NumTraits< std::complex<_Real> >::ReadCost,
        AddCost = 2 * NumTraits< std::complex<_Real> >::AddCost,
        MulCost = 4 * NumTraits< std::complex<_Real> >::MulCost + 2 * NumTraits< std::complex<_Real> >::AddCost,
    };

    static inline Real epsilon() { return NumTraits<Real>::epsilon(); }

    static inline Real dummy_precision() { return NumTraits<Real>::dummy_precision(); }

    static inline int digits10() { return NumTraits<Real>::digits10(); }
};

} // end namespace Eigen

#endif

/*
namespace internal {

template<>
struct abs2_impl<dtmm::DComplexd>
{
  typedef typename NumTraits<dtmm::DComplexd>::Real RealScalar;
  EIGEN_DEVICE_FUNC
  static inline RealScalar run(const dtmm::DComplexd& x)
  {
    return dtmm::abs2(x);
  }
};

} // end namespace internal
*/

