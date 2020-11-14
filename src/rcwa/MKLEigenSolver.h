#pragma once
#include "defns.h"
#include <Eigen/Core>


class MKLEigenSolver
{
private:
	Eigen::MatrixXcs V_;
	Eigen::VectorXcs d_;

public:
	MKLEigenSolver() {}
	void compute(const Eigen::MatrixXcs& A);
	const Eigen::MatrixXcs& eigenvectors() const { return V_; }
	const Eigen::VectorXcs& eigenvalues() const { return d_; }
};