#ifndef MEDIUMSAMPLE_HPP_
#define MEDIUMSAMPLE_HPP_

#include "math/Vec.hpp"

namespace Tungsten {

class PhaseFunction;

struct MediumSample
{
    PhaseFunction *phase;
    Vec3f p;
    float continuedT;
    Vec3f continuedWeight;
    float t;
    Vec3f weight;	// weight is sampled to update the path throughput up to the surface or medium scattering event
    float pdf;	// chances of a photon of light being scattered in a particular direction
    bool exited;
};

}

#endif /* MEDIUMSAMPLE_HPP_ */
