#pragma once
#include "defns.h"
#include <Eigen/Core>
#include <vector>

class FourierSolverPolygon
{
private:
	Eigen::VectorXs alphas_; // the radius from center to each vertex

	scalar dx_;
	scalar dy_;

	scalar offset_; // used for debugging, set the offset of the first radius angle
public:
	FourierSolverPolygon(const Eigen::VectorXs& alphas,
		scalar dx, scalar dy);

	void setOffset(scalar offset) { offset_ = offset; }
	
	void solve(scalex f1, // filling
		scalex f0, // background
		int nx, int ny,
		Eigen::MatrixXcs& Ff,
		const std::vector<int>* layouts = nullptr,
		std::vector<Eigen::MatrixXcs>* dFfs = nullptr);
};