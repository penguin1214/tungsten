#ifndef VOLUME_HPP_
#define VOLUME_HPP_

#include "Primitive.hpp"

class Medium;

namespace Tungsten {

/*
 * Volume is an abstract class,
 * Only used as bounding box, to illustrate volume scale.
 */
class Volume : public Primitive {

private:
	PathPtr _path;	// data file ptr
	Vec3f _pos;
	Vec3f _scale;
	// area?

	// _proxy?
	std::shared_ptr<Medium> _medium_ptr;	// unique_ptr?
	std::shared_ptr<Bsdf> _bsdf;

	Box3f _bounds;	// ??
public:
	Volume();
	
	/// Override Functions
	virtual void fromJson(JsonPtr value, const Scene &scene) override;
	// virtual rapidjson::Value toJson(Allocator &allocator) const override;

	virtual bool intersect(Ray &ray, IntersectionTemporary &data) const override;
	virtual bool occluded(const Ray &ray) const override;
	virtual bool hitBackside(const IntersectionTemporary &data) const override;
	virtual void intersectionInfo(const IntersectionTemporary &data, IntersectionInfo &info) const override;
	virtual bool tangentSpace(const IntersectionTemporary &data, const IntersectionInfo &info,
		Vec3f &T, Vec3f &B) const override;

	virtual bool isSamplable() const override;
	virtual void makeSamplable(const TraceableScene &scene, uint32 threadIndex);

	virtual bool invertParametrization(Vec2f uv, Vec3f &pos) const override;

	virtual bool isDirac() const override;
	virtual bool isInfinite() const override;

	virtual float approximateRadiance(uint32 threadIndex, const Vec3f &p) const override;

	virtual Box3f bounds() const override;

	virtual const TriangleMesh &asTriangleMesh();

	virtual int numBsdfs() const override;
	virtual std::shared_ptr<Bsdf> &bsdf(int index);
	virtual void setBsdf(int index, std::shared_ptr<Bsdf> &bsdf);

	virtual Primitive *clone();
};

}

#endif
