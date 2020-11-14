#ifndef MATRIX_DERIVATIVE_INC
#   define MATRIX_DERIVATIVE_INC

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

namespace dtmm
{
namespace internal
{

/*
 * input is [ A dA ]
 *          [ 0  A ]
 */
template <typename ArgMat, typename RetMat>
void matrix_exp_pade3(const ArgMat& A, const ArgMat& dA, 
                      RetMat& Ua, RetMat& Ub, RetMat& Va, RetMat& Vb)
{
    typedef typename ArgMat::PlainObject MatrixType;
    typedef typename Eigen::NumTraits<typename MatrixType::Scalar>::Real RealScalar;
    const RealScalar b[] = {120.L, 60.L, 12.L, 1.L};

    const MatrixType A2_a = A * A;
    const MatrixType A2_b = A*dA + dA*A;

    const MatrixType tmp_a = b[3]*A2_a + b[1] * MatrixType::Identity(A.rows(), A.cols());
    const MatrixType tmp_b = b[3]*A2_b;

    Ua.noalias() = A * tmp_a;
    Ub.noalias() = A * tmp_b + dA * tmp_a;

    Va.noalias() = b[2]*A2_a + b[0]*MatrixType::Identity(A.rows(), A.cols());
    Vb.noalias() = b[2]*A2_b;
}

template <typename ArgMat, typename RetMat>
void matrix_exp_pade5(const ArgMat& A, const ArgMat& dA, 
                      RetMat& Ua, RetMat& Ub, RetMat& Va, RetMat& Vb)
{
    typedef typename ArgMat::PlainObject MatrixType;
    typedef typename Eigen::NumTraits<typename MatrixType::Scalar>::Real RealScalar;
    const RealScalar b[] = {30240.L, 15120.L, 3360.L, 420.L, 30.L, 1.L};

    const MatrixType A2_a = A * A;
    const MatrixType A2_b = A*dA + dA*A;

    const MatrixType A4_a = A2_a * A2_a;
    const MatrixType A4_b = A2_a * A2_b + A2_b * A2_a;

    const MatrixType tmp_a = b[5]*A4_a + b[3]*A2_a + b[1]*MatrixType::Identity(A.rows(), A.cols());
    const MatrixType tmp_b = b[5]*A4_b + b[3]*A2_b;

    Ua.noalias() = A * tmp_a;
    Ub.noalias() = A * tmp_b + dA * tmp_a;

    Va.noalias() = b[4]*A4_a + b[2]*A2_a + b[0]*MatrixType::Identity(A.rows(), A.cols());
    Vb.noalias() = b[4]*A4_b + b[2]*A2_b;
}

template <typename ArgMat, typename RetMat>
void matrix_exp_pade7(const ArgMat& A, const ArgMat& dA, 
                      RetMat& Ua, RetMat& Ub, RetMat& Va, RetMat& Vb)
{
    typedef typename ArgMat::PlainObject MatrixType;
    typedef typename Eigen::NumTraits<typename MatrixType::Scalar>::Real RealScalar;
    const RealScalar b[] = {17297280.L, 8648640.L, 1995840.L, 277200.L, 25200.L, 1512.L, 56.L, 1.L};

    const MatrixType A2_a = A * A;
    const MatrixType A2_b = A*dA + dA*A;

    const MatrixType A4_a = A2_a * A2_a;
    const MatrixType A4_b = A2_a * A2_b + A2_b * A2_a;

    const MatrixType A6_a = A4_a * A2_a;
    const MatrixType A6_b = A4_a * A2_b + A4_b * A2_a;

    const MatrixType tmp_a = b[7]*A6_a + b[5]*A4_a + b[3]*A2_a + b[1]*MatrixType::Identity(A.rows(), A.cols());
    const MatrixType tmp_b = b[7]*A6_b + b[5]*A4_b + b[3]*A2_b;

    Ua.noalias() = A * tmp_a;
    Ub.noalias() = A * tmp_b + dA * tmp_a;

    Va.noalias() = b[6]*A6_a + b[4]*A4_a + b[2]*A2_a + b[0]*MatrixType::Identity(A.rows(), A.cols());
    Vb.noalias() = b[6]*A6_b + b[4]*A4_b + b[2]*A2_b; 
}

template <typename ArgMat, typename RetMat>
void matrix_exp_pade9(const ArgMat& A, const ArgMat& dA, 
                      RetMat& Ua, RetMat& Ub, RetMat& Va, RetMat& Vb)
{
    typedef typename ArgMat::PlainObject MatrixType;
    typedef typename Eigen::NumTraits<typename MatrixType::Scalar>::Real RealScalar;

    const RealScalar b[] = {17643225600.L, 8821612800.L, 2075673600.L, 302702400.L, 30270240.L,
                            2162160.L, 110880.L, 3960.L, 90.L, 1.L};

    const MatrixType A2_a = A * A;
    const MatrixType A2_b = A*dA + dA*A;

    const MatrixType A4_a = A2_a * A2_a;
    const MatrixType A4_b = A2_a * A2_b + A2_b * A2_a;

    const MatrixType A6_a = A4_a * A2_a;
    const MatrixType A6_b = A4_a * A2_b + A4_b * A2_a;

    const MatrixType A8_a = A4_a * A4_a;
    const MatrixType A8_b = A4_a * A4_b + A4_b * A4_a;

    const MatrixType tmp_a = b[9]*A8_a + b[7]*A6_a + b[5]*A4_a + b[3]*A2_a + b[1]*MatrixType::Identity(A.rows(), A.cols());
    const MatrixType tmp_b = b[9]*A8_b + b[7]*A6_b + b[5]*A4_b + b[3]*A2_b;

    Ua.noalias() = A * tmp_a;
    Ub.noalias() = A * tmp_b + dA * tmp_a;
    Va.noalias() = b[8]*A8_a + b[6]*A6_a + b[4]*A4_a + b[2]*A2_a + b[0]*MatrixType::Identity(A.rows(), A.cols());
    Vb.noalias() = b[8]*A8_b + b[6]*A6_b + b[4]*A4_b + b[2]*A2_b;
}

template <typename ArgMat, typename RetMat>
void matrix_exp_pade13(const ArgMat& A, const ArgMat& dA, 
                      RetMat& Ua, RetMat& Ub, RetMat& Va, RetMat& Vb)
{
    typedef typename ArgMat::PlainObject MatrixType;
    typedef typename Eigen::NumTraits<typename MatrixType::Scalar>::Real RealScalar;
    const RealScalar b[] = {64764752532480000.L, 32382376266240000.L, 7771770303897600.L,
                            1187353796428800.L, 129060195264000.L, 10559470521600.L, 670442572800.L,
                            33522128640.L, 1323241920.L, 40840800.L, 960960.L, 16380.L, 182.L, 1.L};

    const MatrixType A2_a = A * A;
    const MatrixType A2_b = A*dA + dA*A;

    const MatrixType A4_a = A2_a * A2_a;
    const MatrixType A4_b = A2_a * A2_b + A2_b * A2_a;

    const MatrixType A6_a = A4_a * A2_a;
    const MatrixType A6_b = A4_a * A2_b + A4_b * A2_a;

    Va.noalias() = b[13]*A6_a + b[11]*A4_a + b[9]*A2_a; // used for temporary storage
    Vb.noalias() = b[13]*A6_b + b[11]*A4_b + b[9]*A2_b;

    MatrixType tmp_a = A6_a * Va;
    MatrixType tmp_b = A6_a * Vb + A6_b * Va;

    tmp_a += b[7]*A6_a + b[5]*A4_a + b[3]*A2_a + b[1]*MatrixType::Identity(A.rows(), A.cols());
    tmp_b += b[7]*A6_b + b[5]*A4_b + b[3]*A2_b;

    Ua.noalias() = A * tmp_a;
    Ub.noalias() = A * tmp_b + dA * tmp_a;

    tmp_a = b[12]*A6_a + b[10]*A4_a + b[8]*A2_a;
    tmp_b = b[12]*A6_b + b[10]*A4_b + b[8]*A2_b;

    Va.noalias() = A6_a * tmp_a;
    Vb.noalias() = A6_a * tmp_b + A6_b * tmp_a;

    Va += b[6]*A6_a + b[4]*A4_a + b[2]*A2_a + b[0]*MatrixType::Identity(A.rows(), A.cols());
    Vb += b[6]*A6_b + b[4]*A4_b + b[2]*A2_b;
}

template <typename MatrixType, 
          typename RealScalar = typename Eigen::NumTraits<typename MatrixType::Scalar>::Real>
struct matrix_exp_deriv_computeUV
{
    static int run(const MatrixType& A, const MatrixType& dA,
                   MatrixType& Ua, MatrixType& Ub, 
                   MatrixType& Va, MatrixType& Vb);
};

template <typename MatrixType>
struct matrix_exp_deriv_computeUV<MatrixType, double>
{
    typedef typename Eigen::NumTraits<typename MatrixType::Scalar>::Real RealScalar;

    template <typename ArgType>
    static int run(const ArgType& A, const ArgType& dA,
                   MatrixType& Ua, MatrixType& Ub, 
                   MatrixType& Va, MatrixType& Vb)
    {
        using std::frexp;

        const double l1norm = (A.cwiseAbs().colwise().sum() + 
                              dA.cwiseAbs().colwise().sum()).maxCoeff();
        int squarings = 0;

        if (l1norm < 1.495585217958292e-002) {
            matrix_exp_pade3(A, dA, Ua, Ub, Va, Vb);
        } else if (l1norm < 2.539398330063230e-001) {
            matrix_exp_pade5(A, dA, Ua, Ub, Va, Vb);
        } else if (l1norm < 9.504178996162932e-001) {
            matrix_exp_pade7(A, dA, Ua, Ub, Va, Vb);
        } else if (l1norm < 2.097847961257068e+000) {
            matrix_exp_pade9(A, dA, Ua, Ub, Va, Vb);
        } else {
            const RealScalar maxnorm = 5.371920351148152;
            frexp(l1norm / maxnorm, &squarings);
            if (squarings <= 0) 
            {
                squarings = 0;
                matrix_exp_pade13(A, dA, Ua, Ub, Va, Vb);
            }
            else
            {
                MatrixType As  = A.unaryExpr(Eigen::internal::MatrixExponentialScalingOp<double>(squarings));
                MatrixType dAs = dA.unaryExpr(Eigen::internal::MatrixExponentialScalingOp<double>(squarings));
                matrix_exp_pade13(As, dAs, Ua, Ub, Va, Vb);
            }
        } // end if
        return squarings;
    } // end run
};

template <typename ArgType, typename ResultType>
void mat_exp_deriv_compute(const ArgType& A, const ArgType& dA, ResultType& ret)
{
    typedef typename ArgType::PlainObject MatrixType;
    MatrixType Ua, Ub, Va, Vb;

    int squarings = matrix_exp_deriv_computeUV<MatrixType>::run(A, dA, Ua, Ub, Va, Vb);

    MatrixType numer_a = Ua + Va;
    MatrixType numer_b = Ub + Vb;
    MatrixType denom_a = Va - Ua;
    MatrixType denom_b = Vb - Ub;

    // compute the matrix division
    auto lu = denom_a.partialPivLu();
    MatrixType ret_a = lu.solve(numer_a);   // A_1^{-1}A_2
    MatrixType ret_b = lu.solve(denom_b);   // A_1^{-1}V_1
    ret = lu.solve(numer_b);
    ret -= ret_b * ret_a;

    for (int i = 0;i < squarings;++ i)
    {
        ret = ret_a*ret + ret*ret_a;
        ret_a *= ret_a;
    }
}

} // end of namespace internal

// ----------------------------------------------------------------------

template <typename Derived>
class MatInverseDerivReturnByValue : public Eigen::ReturnByValue< MatInverseDerivReturnByValue<Derived> >
{
    typedef typename Derived::Index Index;

    public:
        MatInverseDerivReturnByValue(const Derived& m, const Derived& mprime):m_(m), mprime_(mprime) { }

        Index rows() const { return m_.rows(); }
        Index cols() const { return m_.cols(); }

        template <typename ResultType>
        inline void evalTo(ResultType& result) const
        {
            Eigen::PartialPivLU<Derived> inv(m_);
            Derived mInv = inv.inverse();
            result = -mInv * mprime_ * mInv;
        }

    private:
        const Derived& m_;
        const Derived& mprime_;
};

template <typename Derived>
class MatExpDerivReturnByValue : public Eigen::ReturnByValue< MatExpDerivReturnByValue<Derived> >
{
    typedef typename Derived::Index Index;

    public:
        MatExpDerivReturnByValue(const Derived& m, const Derived& mprime):m_(m), mprime_(mprime) { }

        Index rows() const { return m_.rows(); }
        Index cols() const { return m_.cols(); }

        template <typename ResultType>
        inline void evalTo(ResultType& result) const
        {
            internal::mat_exp_deriv_compute(m_, mprime_, result);
        }

    private:
        const Derived& m_;
        const Derived& mprime_;
};

} // end of namespace dtmm

namespace Eigen {
namespace internal {

/** to support ReturnByValue used in Eigen */
template<typename Derived>
struct traits< dtmm::MatInverseDerivReturnByValue<Derived> >
{
    typedef typename Derived::PlainObject ReturnType;
};

/** to support ReturnByValue used in Eigen */
template<typename Derived>
struct traits< dtmm::MatExpDerivReturnByValue<Derived> >
{
    typedef typename Derived::PlainObject ReturnType;
};

}
}

// ------------------------------------------------------------------

namespace dtmm
{
/*
 * Given M and M', compute the derivative of M^{-1}
 */
template <typename Derived>
const MatInverseDerivReturnByValue<Derived>
mat_inverse_deriv(const Eigen::MatrixBase<Derived>& M,
                  const Eigen::MatrixBase<Derived>& Mprime)
{
    assert(M.rows() == Mprime.rows() && M.cols() == Mprime.cols() && M.rows() == M.cols());
    return MatInverseDerivReturnByValue<Derived>(M.derived(), Mprime.derived());
}

/*
 * Given M and M', compute the derivative of exp(M)
 */
template <typename Derived>
const MatExpDerivReturnByValue<Derived>
mat_exp_deriv(const Eigen::MatrixBase<Derived>& M,
              const Eigen::MatrixBase<Derived>& Mprime)
{
    assert(M.rows() == Mprime.rows() && M.cols() == Mprime.cols() && M.rows() == M.cols());
    return MatExpDerivReturnByValue<Derived>(M.derived(), Mprime.derived());
}

} // end of namespace dtmm

#endif
