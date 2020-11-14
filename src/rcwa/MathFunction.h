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
#include <cmath>


inline scalex sinc(scalar x)
{
    return x == .0 ? scalex(1.) : scalex(sin(Pi * x) / (Pi * x));
}

inline scalar cross2(const Eigen::Vector2s& v1, const Eigen::Vector2s& v2)
{
    return -v1.y() * v2.x() + v1.x() * v2.y();
}