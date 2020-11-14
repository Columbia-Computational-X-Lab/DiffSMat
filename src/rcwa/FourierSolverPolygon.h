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