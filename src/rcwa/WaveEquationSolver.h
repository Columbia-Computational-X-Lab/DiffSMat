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

#include "defns.h"

class WaveEquationSolver
{
    private:
        int nx_; // number of different Fx harmonics
        int ny_; // number of different Fy harmonics
    
        scalar Lx_; // the spatial period along x direction
        scalar Ly_; // the spatial period along y direction
    
    public:
        WaveEquationSolver(int nx, int ny, scalar Lx, scalar Ly):
            nx_(nx), ny_(ny), Lx_(Lx), Ly_(Ly) {}

        void solve(
            scalar lambda,
            scalar wx,
            scalar wy,
            Eigen::MatrixXcs& eigvecE,
            Eigen::MatrixXcs& eigvecH,
            Eigen::VectorXcs& keff) const;
    
        void solve(
            scalar lambda,
            const Eigen::MatrixXcs& P, 
            const Eigen::MatrixXcs& Q,
            Eigen::MatrixXcs& PQ,
            Eigen::MatrixXcs& eigvecE,
            Eigen::MatrixXcs& eigvecH,
            Eigen::VectorXcs& keff) const;
};
