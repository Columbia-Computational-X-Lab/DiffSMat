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

#include "dispersion.h"

scalex permSi(scalar lambda)
{
    return scalex(1+10.6684293*lambda*lambda/(lambda*lambda-pow(.301516485, 2)) + 0.003043478 * lambda * lambda/(lambda*lambda-pow(1.13475115, 2))+1.54133408*lambda * lambda / (lambda * lambda-1104*1104));
}