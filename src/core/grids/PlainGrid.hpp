#ifndef PLAINGRID_HPP_
#define PLAINGRID_HPP_

#include "Grid.hpp"

#include "io/FileUtils.hpp"

namespace Tungsten {

/*
 * PlainGrid class do not support user-define transform now,
 * assuming read-in value is already transformed by user,
 * thus this class only do world-to-local transform.
 */
class PlainGrid : public Grid {
private:
	PathPtr _path;	// data file path
	std::string _gridName;

	int size_x, size_y, size_z;
	float *_data;

	/// scale factor
	Mat4f _transform;
	Mat4f _invTransform;

	float _stepSize;
	Box3f _bounds;	// bbox

public:
	/// Constructor
	PlainGrid();

	/// Override Function
	virtual void fromJson(JsonPtr value, const Scene &scene) override;
	virtual void loadResources() override;

	// virtual Mat4f naturalTransform() const;
	// virtual Mat4f invNaturalTransform() const;
	Mat4f worldToLocalTransform() const;
	Mat4f localToWorldTransform() const;
	int indexAt(Vec3f p) const;
	
	Box3f bounds() const override;

	virtual float density(Vec3f p) const override;
	virtual Vec3f transmittance(PathSampleGenerator &sampler, Vec3f p, Vec3f w/*ray dir*/, float t0, float t1, Vec3f sigmaT) const override;
	virtual Vec2f inverseOpticalDepth(PathSampleGenerator &sampler, Vec3f p, Vec3f w, float t0, float t1,
		float sigmaT, float xi) const override;
};

}


#endif