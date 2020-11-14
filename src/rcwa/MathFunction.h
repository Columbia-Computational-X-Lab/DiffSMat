#pragma once
#include "defns.h"
#include <cmath>


inline scalex sinc(scalar x)
{
	return x == .0 ? scalex(1.) : scalex(sin(Pi * x) / (Pi * x));
}

inline scalar cross2(const Eigen::Vector2s& v1, const Eigen::Vector2s& v2)
{
	return -v1.y() * v2.x() + v1.x() * v2.y();
}