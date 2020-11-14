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
