#ifndef DCOMPLEX_MATRIX_EXPONENTIAL_INC
#   define DCOMPLEX_MATRIX_EXPONENTIAL_INC

#include <Eigen/Core>
#include "DComplex.hpp"
#include "DComplexMatrixFunc.hpp"
#include "utils/timer.hpp"

namespace dtmm
{

namespace internal
{

/** Scaling operator.
 *
 * This struct is used by CwiseUnaryOp to scale a matrix by 2^{-s}.
 */
template <typename RealScalar>
struct MatrixExponentialScalingOp
{
    /** \brief Constructor.
     *
     * \param[in] squarings  The integer \f$ s \f$ in this document.
     */
    MatrixExponentialScalingOp(int squarings) : squarings_(squarings) { }


    /** \brief Scale a matrix coefficient.
     *
     * \param[in,out] x  The scalar to be scaled, becoming \f$ 2^{-s} x \f$.
     */
    inline const RealScalar operator() (const RealScalar& x) const
    {
        using std::ldexp;
        return ldexp(x, -squarings_);
    }

    // -------------------------------------------------------------------

    typedef std::complex<RealScalar> ComplexScalar;

    /** \brief Scale a matrix coefficient.
     *
     * \param[in,out] x  The scalar to be scaled, becoming \f$ 2^{-s} x \f$.
     */
    inline const ComplexScalar operator() (const ComplexScalar& x) const
    {
        using std::ldexp;
        return ComplexScalar(ldexp(x.real(), -squarings_), ldexp(x.imag(), -squarings_));
    }

    // -------------------------------------------------------------------

    typedef DComplex<RealScalar> DComplexScalar;

    inline const DComplexScalar operator() (const DComplexScalar& x) const
    {
        using std::ldexp;
        return DComplexScalar(
            ldexp(x.real_.real(), -squarings_), 
            ldexp(x.real_.imag(), -squarings_),
            ldexp(x.imag_.real(), -squarings_),
            ldexp(x.imag_.imag(), -squarings_));
    }

    private:
        const int squarings_;
};


/** \brief Compute the (3,3)-Pad&eacute; approximant to the exponential.
 *
 *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
 *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
 */
template <typename ArgMat, typename RetMat>
void matrix_exp_pade3(const ArgMat& A, const ArgMat& A2, RetMat& U, RetMat& V)
{
    typedef typename ArgMat::PlainObject MatrixType;
    typedef typename Eigen::NumTraits<typename MatrixType::Scalar>::Real RealScalar;

    const RealScalar b[] = {120.L, 60.L, 12.L, 1.L};
    const MatrixType tmp = b[3] * A2 + b[1] * MatrixType::Identity(A.rows(), A.cols());
    U.noalias() = A * tmp;
    V = b[2] * A2 + b[0] * MatrixType::Identity(A.rows(), A.cols());
}

/** \brief Compute the (5,5)-Pad&eacute; approximant to the exponential.
 *
 *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
 *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
 */
template <typename ArgMat, typename RetMat>
void matrix_exp_pade5(const ArgMat& A, const ArgMat& A2, RetMat& U, RetMat& V)
{
    typedef typename ArgMat::PlainObject MatrixType;
    typedef typename Eigen::NumTraits<typename MatrixType::Scalar>::Real RealScalar;

    const RealScalar b[] = {30240.L, 15120.L, 3360.L, 420.L, 30.L, 1.L};
    const MatrixType A4 = A2 * A2;
    const MatrixType tmp = b[5] * A4 + b[3] * A2 + b[1] * MatrixType::Identity(A.rows(), A.cols());
    U.noalias() = A * tmp;
    V = b[4] * A4 + b[2] * A2 + b[0] * MatrixType::Identity(A.rows(), A.cols());
}

/** \brief Compute the (7,7)-Pad&eacute; approximant to the exponential.
 *
 *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
 *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
 */
template <typename ArgMat, typename RetMat>
void matrix_exp_pade7(const ArgMat& A, const ArgMat& A2, RetMat& U, RetMat& V)
{
    typedef typename ArgMat::PlainObject MatrixType;
    typedef typename Eigen::NumTraits<typename MatrixType::Scalar>::Real RealScalar;

    const RealScalar b[] = {17297280.L, 8648640.L, 1995840.L, 277200.L, 25200.L, 1512.L, 56.L, 1.L};
    const MatrixType A4 = A2 * A2;
    const MatrixType A6 = A4 * A2;
    const MatrixType tmp = b[7] * A6 + b[5] * A4 + b[3] * A2 + b[1] * MatrixType::Identity(A.rows(), A.cols());
    U.noalias() = A * tmp;
    V = b[6] * A6 + b[4] * A4 + b[2] * A2 + b[0] * MatrixType::Identity(A.rows(), A.cols());
}

/** \brief Compute the (9,9)-Pad&eacute; approximant to the exponential.
 *
 *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
 *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
 */
template <typename ArgMat, typename RetMat>
void matrix_exp_pade9(const ArgMat& A, const ArgMat& A2, RetMat& U, RetMat& V)
{
    typedef typename ArgMat::PlainObject MatrixType;
    typedef typename Eigen::NumTraits<typename MatrixType::Scalar>::Real RealScalar;

    const RealScalar b[] = {17643225600.L, 8821612800.L, 2075673600.L, 302702400.L, 30270240.L,
                            2162160.L, 110880.L, 3960.L, 90.L, 1.L};
    const MatrixType A4 = A2 * A2;
    const MatrixType A6 = A4 * A2;
    const MatrixType A8 = A4 * A4;
    const MatrixType tmp = b[9] * A8 + b[7] * A6 + b[5] * A4 + b[3] * A2 + b[1] * MatrixType::Identity(A.rows(), A.cols());
    U.noalias() = A * tmp;
    V = b[8] * A8 + b[6] * A6 + b[4] * A4 + b[2] * A2 + b[0] * MatrixType::Identity(A.rows(), A.cols());
}

/** \brief Compute the (13,13)-Pad&eacute; approximant to the exponential.
 *
 *  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pad&eacute;
 *  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
 */
template <typename ArgMat, typename RetMat>
void matrix_exp_pade13(const ArgMat& A, const ArgMat& A2, RetMat& U, RetMat& V, std::false_type)
{
    typedef typename ArgMat::PlainObject MatrixType;
    typedef typename Eigen::NumTraits<typename MatrixType::Scalar>::Real RealScalar;

    const RealScalar b[] = {64764752532480000.L, 32382376266240000.L, 7771770303897600.L,
                            1187353796428800.L, 129060195264000.L, 10559470521600.L, 670442572800.L,
                            33522128640.L, 1323241920.L, 40840800.L, 960960.L, 16380.L, 182.L, 1.L};
    const MatrixType A4 = A2 * A2;
    const MatrixType A6 = A4 * A2;
    V = b[13] * A6 + b[11] * A4 + b[9] * A2;    // used for temporary storage
    MatrixType tmp = A6 * V;                    // b13*A12 + b11*A10 + b9*A8
    tmp += b[7] * A6 + b[5] * A4 + b[3] * A2 + b[1] * MatrixType::Identity(A.rows(), A.cols());
    U.noalias() = A * tmp;
    tmp = b[12] * A6 + b[10] * A4 + b[8] * A2;
    V.noalias() = A6 * tmp;
    V += b[6] * A6 + b[4] * A4 + b[2] * A2 + b[0] * MatrixType::Identity(A.rows(), A.cols());
}

/**
 * Specific implementation for DComplex matrix.
 */
template <typename ArgMat, typename RetMat>
void matrix_exp_pade13(const ArgMat& A, const ArgMat& A2, RetMat& U, RetMat& V, std::true_type)
{
    typedef typename ArgMat::PlainObject MatrixType;
    typedef typename Eigen::NumTraits<typename MatrixType::Scalar>::Real RealScalar;

    const RealScalar b[] = {64764752532480000.L, 32382376266240000.L, 7771770303897600.L,
                            1187353796428800.L, 129060195264000.L, 10559470521600.L, 670442572800.L,
                            33522128640.L, 1323241920.L, 40840800.L, 960960.L, 16380.L, 182.L, 1.L};
    MatrixType A4, A6, tmp;
    dcomplex_mat_product(A2, A2, A4);
    dcomplex_mat_product(A2, A4, A6);

    //V = b[13] * A6 + b[11] * A4 + b[9] * A2;    // used for temporary storage
    dcomplex_mat_sum(b[13], A6, b[11], A4, b[9], A2, V);
    
    dcomplex_mat_product(A6, V, tmp);

    //tmp += b[7] * A6 + b[5] * A4 + b[3] * A2 + b[1] * MatrixType::Identity(A.rows(), A.cols());
    dcomplex_mat_axpy(b[7], A6, b[5], A4, b[3], A2, b[1], tmp);

    dcomplex_mat_product(A, tmp, U);

    //tmp = b[12] * A6 + b[10] * A4 + b[8] * A2;
    dcomplex_mat_sum(b[12], A6, b[10], A4, b[8], A2, tmp);
    
    dcomplex_mat_product(A6, tmp, V);

    //V += b[6] * A6 + b[4] * A4 + b[2] * A2 + b[0] * MatrixType::Identity(A.rows(), A.cols());
    dcomplex_mat_axpy(b[6], A6, b[4], A4, b[2], A2, b[0], V);
}

template <typename MatrixType, 
          typename ScalarType = typename Eigen::NumTraits<typename MatrixType::Scalar>::Real>
struct matrix_exp_sqr_computeUV
{
    /*
     * Implement the paper: "THE SCALING AND SQUARING METHOD FOR THE MATRIX EXPONENTIAL REVISITED" by Higham 05
     *
     * Compute U and V such that (V+U)(V-U)^{-1} is a Pade approximation of exp(M). 
     * Return an integer of squarings used for scale the matrix
     *
     * M:  the matrix for which we comptue exp(M)
     * M2: the matrix that holds the matrix square, M*M
     */
    static int run(const MatrixType& M, const MatrixType& M2, MatrixType& U, MatrixType& V);
};

/*
 * Implementation with double precision
 */
template <typename MatrixType>
struct matrix_exp_sqr_computeUV<MatrixType, double>
{
    template <typename ArgType>
    static int run(const ArgType& M, const ArgType& M2, MatrixType& U, MatrixType& V)
    {
        using std::frexp;

        const double l1norm = M.cwiseAbs().colwise().sum().maxCoeff();
        int squarings = 0;
        
        if (l1norm < 1.495585217958292e-002) {
            matrix_exp_pade3(M, M2, U, V);
        } else if (l1norm < 2.539398330063230e-001) {
            matrix_exp_pade5(M, M2, U, V);
        } else if (l1norm < 9.504178996162932e-001) {
            matrix_exp_pade7(M, M2, U, V);
        } else if (l1norm < 2.097847961257068e+000) {
            matrix_exp_pade9(M, M2, U, V);
        } else {
            const double maxnorm = 5.371920351148152;
            frexp(l1norm / maxnorm, &squarings);

            if (squarings <= 0) 
            {
                squarings = 0;
                matrix_exp_pade13(M, M2, U, V, is_dcomplex_type<typename ArgType::Scalar>());
                //matrix_exp_pade13(M, M2, U, V, std::false_type());
            } else {
                MatrixType A  = M.unaryExpr(MatrixExponentialScalingOp<double>(squarings));
                MatrixType A2 = M2.unaryExpr(MatrixExponentialScalingOp<double>(squarings*2));
                matrix_exp_pade13(A, A2, U, V, is_dcomplex_type<typename ArgType::Scalar>());
                //matrix_exp_pade13(A, A2, U, V, std::false_type());
            }
        } // end if

        return squarings;
    } // end run
};

#if 0
template <typename ArgType, typename ResultType>
void compute_expm_with_sqr(const ArgType& m, const ArgType& m2, ResultType& ret)
{
    typedef typename ArgType::PlainObject MatrixType;

    typename MatrixType::Scalar mu = m.trace() / m.cols();

    MatrixType newM  = m - MatrixType::Identity(m.rows(), m.cols())*mu;
    MatrixType newM2 = m2 - m*(mu*2.) +  MatrixType::Identity(m.rows(), m.cols())*(mu*mu);
    MatrixType U, V;
    int squarings = matrix_exp_sqr_computeUV<MatrixType>::run(newM, newM2, U, V);

    MatrixType numer = V + U;
    MatrixType denom = V - U;

    ret = denom.partialPivLu().solve(numer);
    for (int i = 0;i < squarings;++ i)
        ret *= ret;
    ret *= std::exp(mu);
}
#endif

template <typename ArgType, typename ResultType>
void compute_expm_with_sqr(const ArgType& m, const ArgType& m2, ResultType& ret, std::true_type)
{
    typedef typename ArgType::PlainObject MatrixType;

    MatrixType U, V;
    
    //auto t1 = sploosh::now();
    int squarings = matrix_exp_sqr_computeUV<MatrixType>::run(m, m2, U, V);
    //auto t2 = sploosh::now();

    MatrixType numer = V + U;
    MatrixType denom = V - U;

    ret = denom.partialPivLu().solve(numer);

    //auto t3 = sploosh::now();
    for (int i = 0;i < squarings;++ i)
    {
        if ( i % 2 )
            dcomplex_mat_product(denom, denom, ret);
        else
            dcomplex_mat_product(ret, ret, denom);
    }
    if ( squarings % 2 == 1 ) ret = denom;

    //auto t4 = sploosh::now();
    //std::cout << "expm Timings: (" << squarings << ")  " 
    //          << sploosh::duration_milli_d(t1, t2) << "ms  " 
    //          << sploosh::duration_milli_d(t2, t3) << "ms  "  
    //          << sploosh::duration_milli_d(t3, t4) << "ms  "  
    //          << std::endl;
}

template <typename ArgType, typename ResultType>
void compute_expm_with_sqr(const ArgType& m, const ArgType& m2, ResultType& ret, std::false_type)
{
    typedef typename ArgType::PlainObject MatrixType;

    MatrixType U, V;
    
    int squarings = matrix_exp_sqr_computeUV<MatrixType>::run(m, m2, U, V);

    MatrixType numer = V + U;
    MatrixType denom = V - U;

    ret = denom.partialPivLu().solve(numer);

    for (int i = 0;i < squarings;++ i)
        ret *= ret;
}

} // end of namespace dtmm::internal 

template <typename Derived>
class ExpmWithSqrReturnValue : public Eigen::ReturnByValue< ExpmWithSqrReturnValue<Derived> >
{
    typedef typename Derived::Index Index;

    public:
        ExpmWithSqrReturnValue(const Derived& m, const Derived& m2):m_(m), m2_(m2) { }

        Index rows() const { return m_.rows(); }
        Index cols() const { return m_.cols(); }

        template <typename ResultType>
        inline void evalTo(ResultType& result) const
        {
            internal::compute_expm_with_sqr(m_, m2_, result, 
                    is_dcomplex_type<typename Derived::Scalar>());
        }

    private:
        const Derived& m_;
        const Derived& m2_;
};

} // end of namespace dtmm

// ------------------------------------------------------------------

namespace Eigen {
namespace internal {

/** to support ReturnByValue used in Eigen */
template<typename Derived>
struct traits< dtmm::ExpmWithSqrReturnValue<Derived> >
{
    typedef typename Derived::PlainObject ReturnType;
};

}
}

// ------------------------------------------------------------------

namespace dtmm {

/*
 * Computes exp(M) with the M^2 explicitly provided.
 */
template <typename Derived>
const ExpmWithSqrReturnValue<Derived>
expm_with_sqr(const Eigen::MatrixBase<Derived>& m, const Eigen::MatrixBase<Derived>& m2)
{
    assert(m.rows() == m.cols() && m2.rows() == m2.cols() && m.rows() == m2.rows());
    return ExpmWithSqrReturnValue<Derived>(m.derived(), m2.derived());
}

template <typename Derived, typename ResultType>
void expm_with_sqr(const Eigen::MatrixBase<Derived>& m, 
        const Eigen::MatrixBase<Derived>& m2,
        ResultType& result)
{
    assert(m.rows() == m.cols() && m2.rows() == m2.cols() && m.rows() == m2.rows());

    result.resize(m.rows(), m.cols());
    internal::compute_expm_with_sqr(m, m2, result,
            is_dcomplex_type<typename Derived::Scalar>());
}

}

#endif
