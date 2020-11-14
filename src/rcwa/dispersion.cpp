#include "dispersion.h"

scalex permSi(scalar lambda)
{
    return scalex(1+10.6684293*lambda*lambda/(lambda*lambda-pow(.301516485, 2)) + 0.003043478 * lambda * lambda/(lambda*lambda-pow(1.13475115, 2))+1.54133408*lambda * lambda / (lambda * lambda-1104*1104));
}