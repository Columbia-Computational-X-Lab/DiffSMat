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

#ifndef RCWA_SCATTERING_MATRIX_INC
#   define RCWA_SCATTERING_MATRIX_INC

#include <Eigen/Dense>
#include "rcwa/WaveEquationSolver.h"

namespace dtmm 
{

class RCWAScatterMatrix
{
    public:
        typedef Eigen::MatrixXcs                MatrixType;
        typedef Eigen::VectorXcs                VectorType;
        typedef Eigen::PartialPivLU<MatrixType> SolverType;

        RCWAScatterMatrix(int nx, int ny, scalar Lx, scalar Ly):
                nx_(nx), ny_(ny), Lx_(Lx), Ly_(Ly), isComputed_(false) { }

        /*
         * Compute scattering matrix and store all the intermediate 
         * matrices needed for derivative computation
         */
        void compute(scalar lambda, scalar wz, 
                     const MatrixType& P,
                     const MatrixType& Q);

        /*
         * Compute scattering matrix derivative along z direction
         */
        void compute_dz(scalar lambda, scalar wz);

        inline const MatrixType& Tuu() const
        {   return Tuu_; }
        inline const MatrixType& Rud() const
        {   return Rud_; }
        inline const MatrixType& Rdu() const
        {   return Rud_; }
        inline const MatrixType& Tdd() const
        {   return Tuu_; }
        inline const MatrixType& dTuuz() const
        {   return dTuuz_; }
        inline const MatrixType& dRudz() const
        {   return dRudz_; }
        inline const MatrixType& dRduz() const
        {   return dRudz_; }
        inline const MatrixType& dTddz() const
        {   return dTuuz_; }
        inline const MatrixType& W() const
        {   return W_; }
        inline const MatrixType& PQ() const
        {   return PQ_; }
        inline const MatrixType& W0() const
        {   return W0_; }
        inline const MatrixType& V0() const
        {   return V0_; }
        inline scalar Lz_normalized() const
        {   return normalizedLz_; }

        inline const VectorType& Gamma() const
        {   return gamma_; }
        inline const VectorType& X() const
        {   return X_; }

        inline const MatrixType& invAB() const
        {   return invAB_; }
        inline const MatrixType& invAXA() const
        {   return invAXA_; }
        inline const MatrixType& invAXB() const
        {   return invAXB_; }
        inline const MatrixType& I_minus_invAXB2() const
        {   return I_minus_invAXB2_; }
        inline const MatrixType& invAXB_invAXA_invAB() const
        {   return invAXB_invAXA_invAB_; }
        inline const MatrixType& invAXA_invAXB_invAB() const
        {   return invAXA_invAXB_invAB_; }

        inline const SolverType& invW_sol() const  // W^{-1}
        {   return invWSol_; }
        inline const SolverType& inv_I_minus_invAXB2_sol() const
        {   return inv_I_minus_invAXB2_Sol_; }

    private:
    // assume no PML layer
    void evaluateHomogeneousLayer(
        scalar lambda, scalex eps,
        MatrixType& eigvecE,
        MatrixType& eigvecH,
        VectorType& keff);

    private:
        int     nx_, ny_;
        scalar  Lx_, Ly_;
        scalar  normalizedLz_;

        MatrixType Tuu_, Rud_;
        MatrixType W0_, V0_; 
        MatrixType PQ_, W_, V_;
        MatrixType A_, B_;
        MatrixType invAB_, invAXB_, invAXA_;
        MatrixType I_minus_invAXB2_;
        MatrixType invAXB_invAXA_invAB_;
        MatrixType invAXA_invAXB_invAB_;

        MatrixType dTuuz_, dRudz_;


        VectorType gamma_, X_;

        SolverType invWSol_;
        SolverType inv_I_minus_invAXB2_Sol_;   // I - (invA*X*B)^2
        SolverType invASol_;

        bool isComputed_;
    //public:
    //    MatrixType ret_;    // for sanity check
};

}

#endif
