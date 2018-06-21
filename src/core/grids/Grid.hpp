#ifndef GRID_HPP_
#define GRID_HPP_

#include "math/Mat4f.hpp"
#include "math/Box.hpp"

#include "io/JsonSerializable.hpp"

namespace Tungsten {

class PathSampleGenerator;

class Grid : public JsonSerializable
{
public:
    virtual ~Grid() {}

	/* naturalTransform() is worldToLocal transform combine **user-define transform**. */
	virtual Mat4f scaleTransform() const;
	virtual Mat4f invScaleTransform() const;
	virtual Mat4f configTransform() const;
	virtual Mat4f invConfigTransform() const;
    virtual Mat4f naturalTransform() const;
    virtual Mat4f invNaturalTransform() const;
	//virtual Mat4f worldToLocalNaturalTransform() const;
	//virtual Mat4f localToWorldNaturalTransform() const;
    virtual Box3f bounds() const;
	virtual float stepSize() const = 0;

    virtual float density(Vec3f p) const = 0;
	virtual Vec3f gradient(Vec3f p) const = 0;
	virtual float maxDensity() = 0;
	virtual Vec3f maxGradient() = 0;

	/* Compute transparency, i.e. e^-(integral of absorption between [t0, t1], along Ray(p, w))
	 * i.e. how much gets transmitted (pass through)
	 */
    virtual Vec3f transmittance(PathSampleGenerator &sampler, Vec3f p, Vec3f w, float t0, float t1, Vec3f sigmaT) const = 0;

    virtual Vec2f inverseOpticalDepth(PathSampleGenerator &sampler, Vec3f p, Vec3f w, float t0, float t1,
            float sigmaT, float xi) const = 0;
};

}

#endif /* GRID_HPP_ */
