#include "defns.h"
#include "LayerStructure.h"

class WaveEquationCoeff
{
private:
	int nx_; // number of different Fx harmonics
	int ny_; // number of different Fy harmonics

	scalar Lx_; // the spatial period along x direction
	scalar Ly_; // the spatial period along y direction

public:
	WaveEquationCoeff(int nx, int ny, scalar Lx, scalar Ly):
		nx_(nx), ny_(ny), Lx_(Lx), Ly_(Ly) {}

	void solve(
		scalar lambda,
		scalar wx, // width of scatterer along x direction
		scalar wy, // width of scatterer along y direction
		Eigen::MatrixXcs& P,
		Eigen::MatrixXcs& Q,
		Eigen::MatrixXcs* dPx = nullptr, // partial P / partial wx
		Eigen::MatrixXcs* dPy = nullptr, // partial P / partial wy
		Eigen::MatrixXcs* dQx = nullptr, // partial Q / partial wx
		Eigen::MatrixXcs* dQy = nullptr);// partial Q / partial wy

	void solve(
		scalar lambda,
		const struct LayerStructure& layer,
		int layoutX,
		int layoutY,
		Eigen::MatrixXcs& P,
		Eigen::MatrixXcs& Q,
		Eigen::MatrixXcs* dPx = nullptr, // partial P / partial x_i
		Eigen::MatrixXcs* dPy = nullptr, // partial P / partial y_i
		Eigen::MatrixXcs* dQx = nullptr, // partial Q / partial x_i
		Eigen::MatrixXcs* dQy = nullptr);// partial Q / partial y_i)

	void solve(
		scalar lambda,
		const Eigen::VectorXs& alphas,
		Eigen::MatrixXcs& P,
		Eigen::MatrixXcs& Q,		
		const std::vector<int>* layouts = nullptr,
		std::vector<Eigen::MatrixXcs>* dPs = nullptr,
		std::vector<Eigen::MatrixXcs>* dQs = nullptr);
        
	void evaluateKMatrices(
		Eigen::MatrixXcs& Kx,
		Eigen::MatrixXcs& Ky);

private:
	void evaluateAlphaBeta_(
		Eigen::DiagonalMatrixXcs& alpha,
		Eigen::DiagonalMatrixXcs& beta);

	void constructLayer_(
        scalar lambda,
		Eigen::VectorXs& coordX,
		Eigen::VectorXs& coordY,
		Eigen::MatrixXcs& eps,
		scalar wx, scalar wy);
};