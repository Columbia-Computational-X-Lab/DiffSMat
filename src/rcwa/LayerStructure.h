#pragma once

#include "defns.h"
#include <Eigen/Core>

struct LayerStructure
{
	Eigen::VectorXs coordX;
	Eigen::VectorXs coordY;
	Eigen::MatrixXcs eps;
};