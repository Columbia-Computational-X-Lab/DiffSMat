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

#ifndef SPLOOSH_UTILS_TIMER
#   define SPLOOSH_UTILS_TIMER

#include <chrono>

namespace sploosh
{

inline auto now() 
{
    return std::chrono::high_resolution_clock::now();
}

template <class Clock>
inline double duration_milli_d(const std::chrono::time_point<Clock>& t1, 
                               const std::chrono::time_point<Clock>& t2)
{
    return std::chrono::duration<double, std::milli>(t2 - t1).count();
}

}

#endif
