#ifndef PLAINGRID_HPP_
#define PLAINGRID_HPP_

#include "Grid.hpp"

#include "io/FileUtils.hpp"
#include "../common.hpp"

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
	Vec3f *_gradient;

	/// scale factor
	Mat4f _scaleTransform;
	Mat4f _invScaleTransform;

	Mat4f _configTransform;
	Mat4f _invConfigTransform;

	float _stepSize;
	Box3f _bounds;	// bbox

	float _maxDensity;
	Vec3f _maxGradient;
public:
	/// Constructor
	PlainGrid();

	/// Override Function
	virtual void fromJson(JsonPtr value, const Scene &scene) override;
	virtual void loadResources() override;

	Mat4f scaleTransform() const override;	// naturalTransform() is worldToLocalTransform()
	Mat4f invScaleTransform() const override;	// invNaturalTransform is localToWorldTransform(0
	Mat4f configTransform() const override;
	Mat4f invConfigTransform() const override;
	Mat4f naturalTransform() const override;
	Mat4f invNaturalTransform() const override;
	//virtual Mat4f worldToLocalNaturalTransform() const;
	//virtual Mat4f localToWorldNaturalTransform() const;

	int indexAt(Vec3f p) const;
	int indexAt(int x, int y, int z) const;
	
	Box3f bounds() const override;

	virtual float density(Vec3f p) const override;
	virtual Vec3f gradient(Vec3f p) const override;

	virtual Vec3f transmittance(PathSampleGenerator &sampler, Vec3f p, Vec3f w/*ray dir*/, float t0, float t1, Vec3f sigmaT) const override;
	virtual Vec2f inverseOpticalDepth(PathSampleGenerator &sampler, Vec3f p, Vec3f w, float t0, float t1,
		float sigmaT, float xi) const override;

	virtual float maxDensity() override;
	virtual Vec3f maxGradient() override;
	float precomputeMaxDensity();
	void precomputeGradient();

	virtual float stepSize() const;
};

}


#endif