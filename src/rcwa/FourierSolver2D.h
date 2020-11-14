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

enum FourierSolverOption
{
    ContinuousXY = 0,
    DiscontinuousX = 1,
    DiscontinuousY = 2
};


class FourierSolver2D
{
private:
    Eigen::VectorXs coordX_;
    Eigen::VectorXs coordY_;

    scalar dx_;
    scalar dy_;

public:
    FourierSolver2D(const Eigen::VectorXs& coordX,
        const Eigen::VectorXs& coordY):
        coordX_(coordX), coordY_(coordY) {
        dx_ = coordX_(coordX_.size() - 1) - coordX_(0);
        dy_ = coordY_(coordY_.size() - 1) - coordY_(0);
    }

    void solve(const Eigen::MatrixXcs& f,
        int nx, int ny,
        FourierSolverOption opt,
        Eigen::MatrixXcs& Ff,
        int layoutX = -1,
        Eigen::MatrixXcs* dFfx = nullptr,
        int layoutY = -1,
        Eigen::MatrixXcs* dFfy = nullptr);

private:
    bool isValidX(int layoutX, Eigen::MatrixXcs* dFfx) {
        return layoutX > 0 && layoutX < coordX_.size() - 1 && dFfx != nullptr;
    }
    bool isValidY(int layoutY, Eigen::MatrixXcs* dFfy) {
        return layoutY > 0 && layoutY < coordY_.size() - 1 && dFfy != nullptr;
    }
};
