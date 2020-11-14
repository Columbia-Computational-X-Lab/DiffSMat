#ifndef SCATTERING_MATRIX_DERIVATIVE
#   define SCATTERING_MATRIX_DERIVATIVE

#include <tbb/tbb.h>
#include <Eigen/Core>
#include "MatrixDerivative.hpp"
#include "dcomplex/DComplexMath.hpp"
#include "dcomplex/EigenNumTraits.hpp"
#include "dcomplex/DComplexMatrixFunc.hpp"
#include "dcomplex/MatrixExponential.hpp"
#include "utils/timer.hpp"

namespace dtmm
{

/*
 * ScatterMatrixSolver provide the access to the scattering matrix 
 * and its intermediate matrix results.
 */
template <class ScatterMatrixSolver>
class ScatterMatrixDerivative
{
    typedef typename ScatterMatrixSolver::MatrixType MatrixType;
    typedef typename ScatterMatrixSolver::VectorType VectorType;
    typedef typename ScatterMatrixSolver::SolverType SolverType;

    public:
        /*
         * Usage:
         * - use WaveEquationCoeff to compute P, Q, dP, and dQ
         * - compute Scattering Matrix to set sMatSolver
         * - call this method
         */
        void compute(const ScatterMatrixSolver& sMatSolver,
                     const MatrixType& P, const MatrixType& Q,
                     const MatrixType& dP, const MatrixType& dQ);

        /* retrieve the resulting derivative */
        const MatrixType& dTuu() const
        {   return dTuu_; }
        const MatrixType& dTdd() const
        {   return dTuu_; }
        const MatrixType& dRud() const
        {   return dRud_; }
        const MatrixType& dRdu() const
        {   return dRud_; }

        // --- for sanity check ---
        const MatrixType& Omega() const
        {   return Omega_; }
        const MatrixType& dOmega() const
        {   return dOmega_; }
        const MatrixType& WXinvW() const
        {   return WXinvW_; }
        const MatrixType& dWXinvW() const
        {   return dWXinvW_; }

    private:
        void deriv_Omega(const ScatterMatrixSolver& sMatSolver,
                         const MatrixType& P, const MatrixType& Q,
                         const MatrixType& dP, const MatrixType& dQ);

        void deriv_T(const ScatterMatrixSolver& sMatSolver,
                     const MatrixType& Q, const MatrixType& dQ,
                     MatrixType& TV0, MatrixType& dTV0);

    private:
        MatrixType  Omega_, dOmega_;
        MatrixType  WXinvW_, dWXinvW_;  // WXinvW = exp^{j*Omega*L/k_0}
        MatrixType  dTuu_, dRud_;

    //public:
    //    MatrixType  ret_;
};

// -----------------------------------------------------------------------

template <class ScatterMatrixSolver>
void ScatterMatrixDerivative<ScatterMatrixSolver>::compute(
        const ScatterMatrixSolver& sMatSolver,
        const ScatterMatrixDerivative<ScatterMatrixSolver>::MatrixType& P, 
        const ScatterMatrixDerivative<ScatterMatrixSolver>::MatrixType& Q,
        const ScatterMatrixDerivative<ScatterMatrixSolver>::MatrixType& dP, 
        const ScatterMatrixDerivative<ScatterMatrixSolver>::MatrixType& dQ)
{
    MatrixType TV0, dTV0;
    // compute Omega and Omega'
    //auto t1 = sploosh::now();
    deriv_Omega(sMatSolver, P, Q, dP, dQ);
    // compute TV0 and dTV0
    //auto t2 = sploosh::now();
    deriv_T(sMatSolver, Q, dQ, TV0, dTV0);

    //auto t3 = sploosh::now();
    // --- compute (A^{-1}B)' ---
    const MatrixType W0_plus_TV0 = sMatSolver.W0() + TV0;
    const MatrixType W0_minus_TV0 = sMatSolver.W0() - TV0;
    const MatrixType dInv_W0_plus_TV0 = mat_inverse_deriv(W0_plus_TV0, dTV0);

    SolverType inv_W0_plus_TV0( W0_plus_TV0 );

    const MatrixType dInvAB = dInv_W0_plus_TV0 * W0_minus_TV0 - inv_W0_plus_TV0.solve(dTV0);
    // --------------------------

    //auto t4 = sploosh::now();
    // --- compute (A^{-1}XA)' and (A^{-1}XB)' ---
    MatrixType tmpT1 = dInv_W0_plus_TV0 * WXinvW_;
    MatrixType tmpT2 = inv_W0_plus_TV0.solve(dWXinvW_);
    MatrixType tmpT3 = inv_W0_plus_TV0.solve(WXinvW_ * dTV0);
    const MatrixType dInvAXA = (tmpT1 + tmpT2) * W0_plus_TV0 + tmpT3;
    const MatrixType dInvAXB = (tmpT1 + tmpT2) * W0_minus_TV0 - tmpT3;

    //auto t5 = sploosh::now();
    // === Finally, compute the derivative of scattering sub-matrices
    tmpT1 = - sMatSolver.invAXB()*dInvAXB - dInvAXB*sMatSolver.invAXB();
    tmpT2 = mat_inverse_deriv(sMatSolver.I_minus_invAXB2(), tmpT1);
    dRud_ = tmpT2 * sMatSolver.invAXB_invAXA_invAB() +
            sMatSolver.inv_I_minus_invAXB2_sol().solve(
                    dInvAXB*sMatSolver.invAXA() + sMatSolver.invAXB()*dInvAXA - dInvAB);
    dTuu_ = tmpT2 * sMatSolver.invAXA_invAXB_invAB() + 
            sMatSolver.inv_I_minus_invAXB2_sol().solve(
                    dInvAXA - dInvAXB*sMatSolver.invAB() - sMatSolver.invAXB()*dInvAB);
    //auto t6 = sploosh::now();
    //std::cout << "Timings: " 
    //          << sploosh::duration_milli_d(t1, t2) << "ms  " 
    //          << sploosh::duration_milli_d(t2, t3) << "ms  "  
    //          << sploosh::duration_milli_d(t3, t4) << "ms  "  
    //          << sploosh::duration_milli_d(t4, t5) << "ms  "  
    //          << sploosh::duration_milli_d(t5, t6) << "ms  "  
    //          << std::endl;
}

template <class ScatterMatrixSolver>
void ScatterMatrixDerivative<ScatterMatrixSolver>::deriv_Omega(
        const ScatterMatrixSolver& sMatSolver,
        const ScatterMatrixDerivative<ScatterMatrixSolver>::MatrixType& P, 
        const ScatterMatrixDerivative<ScatterMatrixSolver>::MatrixType& Q,
        const ScatterMatrixDerivative<ScatterMatrixSolver>::MatrixType& dP, 
        const ScatterMatrixDerivative<ScatterMatrixSolver>::MatrixType& dQ)
{
    //auto t1 = sploosh::now();
    const VectorType& gamma = sMatSolver.Gamma();
    const MatrixType& W = sMatSolver.W();
    const MatrixType invW = sMatSolver.invW_sol().inverse();

    // --- compute Omega = W*\Gamma*W^{-1} ---
    Omega_ = W * gamma.asDiagonal() * invW;
    //auto t2 = sploosh::now();

    // --- compute Omega' ---
    const MatrixType dPQ = dP*Q + P*dQ;
    MatrixType C = invW * dPQ * W;
    tbb::parallel_for( tbb::blocked_range2d<Eigen::Index>(0, C.rows(), 0, C.cols()),
            [&](const tbb::blocked_range2d<Eigen::Index>& r) {
                for(Eigen::Index i = r.rows().begin();i != r.rows().end();++ i) {
                    for(Eigen::Index j = r.cols().begin();j != r.cols().end(); ++ j) {
                        C(i,j) /= (gamma(i) + gamma(j));
                    }
                } // end for
            } // end lambda func
    );

    dOmega_ = W * C * invW;     // compute Omega' = W * Y * W^{-1} (here C stores Y matrix)
    //auto t3 = sploosh::now();

    // -----------------------------------------------------------

    // --- compute exp^{j*Omega*L/k_0} ---
    WXinvW_ = W * sMatSolver.X().asDiagonal() * invW;

    //auto t4 = sploosh::now();

    // --- compute the derivative of exp^{j*Omega*L/k_0} ---
    typedef typename Eigen::NumTraits<typename MatrixType::Scalar>::Real Real;
    const Real s = sMatSolver.Lz_normalized();

    //static_assert(sizeof(Real) == 8, "Only double is considered for now");
    //auto t5 = sploosh::now();
#ifdef USE_DCOMPLEX_DERIVATIVE
    std::cout << "use dcomplex deriv" << std::endl;
    typedef Eigen::Matrix< DComplex<Real>, Eigen::Dynamic, Eigen::Dynamic > DCMat;

    const Real H = 1E-10;
    DCMat dcOmega = dcomplex_mat_from_pair(Omega_, dOmega_, std::complex<Real>(0, s), H);
    DCMat dcOmegaSqr = dcomplex_mat_from_pair(sMatSolver.PQ(), dPQ, std::complex<Real>(-s*s), H);
    DCMat ret = expm_with_sqr(dcOmega, dcOmegaSqr);
    dcomplex_mat_deriv(ret, H, dWXinvW_);
#else
    std::cout << "use mat exp deriv" << std::endl;
    MatrixType omega  = Omega_ * std::complex<Real>(0, s);
    MatrixType domega = dOmega_ * std::complex<Real>(0, s);

    dWXinvW_ = mat_exp_deriv(omega, domega);
#endif
    //auto t6 = sploosh::now();

    //std::cout << "=>> Timings: " 
    //          << sploosh::duration_milli_d(t1, t2) << "ms  " 
    //          << sploosh::duration_milli_d(t2, t3) << "ms  "  
    //          << sploosh::duration_milli_d(t3, t4) << "ms  "  
    //          << sploosh::duration_milli_d(t4, t5) << "ms  "  
    //          << sploosh::duration_milli_d(t5, t6) << "ms  "  
    //          << std::endl;
} // end

template <class ScatterMatrixSolver>
void ScatterMatrixDerivative<ScatterMatrixSolver>::deriv_T(
        const ScatterMatrixSolver& sMatSolver,
        const ScatterMatrixDerivative<ScatterMatrixSolver>::MatrixType& Q,
        const ScatterMatrixDerivative<ScatterMatrixSolver>::MatrixType& dQ,
        ScatterMatrixDerivative<ScatterMatrixSolver>::MatrixType& TV0,
        ScatterMatrixDerivative<ScatterMatrixSolver>::MatrixType& dTV0)
{
    typename ScatterMatrixSolver::SolverType invQ(Q);
    const MatrixType invQV0 = invQ.solve(sMatSolver.V0());

    // --- compute T ---
    TV0 = Omega_ * invQV0; // Omega * Q^{-1} * V0

    // --- compute dTV0 ---
    // Omega'*Q^{-1}*V0 - Omega * Q^{-1} * Q' * Q^{-1} * V0
    dTV0 = dOmega_ * invQV0 - Omega_ * invQ.solve(dQ * invQV0);
}

} // end of namespace dtmm

#endif
