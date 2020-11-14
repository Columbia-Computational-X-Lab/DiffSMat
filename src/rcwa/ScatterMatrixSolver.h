#include "defns.h"

class ScatterMatrixSolver
{
private:
	int nx_; // number of different Fx harmonics
	int ny_; // number of different Fy harmonics

	scalar Lx_; // the spatial period along x direction
	scalar Ly_; // the spatial period along y direction
public:
	ScatterMatrixSolver(int nx, int ny, scalar Lx, scalar Ly):
		nx_(nx), ny_(ny), Lx_(Lx), Ly_(Ly) {}

	void solve(
	scalar lambda,
	scalar wx,
	scalar wy,
	scalar wz,
	Eigen::MatrixXcs& Tuu,
	Eigen::MatrixXcs& Rud,
	Eigen::MatrixXcs& Rdu,
	Eigen::MatrixXcs& Tdd);


	// assume no PML layer
	void evaluateHomogeneousLayer(
		scalar lambda,
		scalex eps,
		Eigen::MatrixXcs& eigvecE,
		Eigen::MatrixXcs& eigvecH,
		Eigen::VectorXcs& keff);

};
