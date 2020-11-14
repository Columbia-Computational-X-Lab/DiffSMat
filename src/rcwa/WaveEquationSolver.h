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
