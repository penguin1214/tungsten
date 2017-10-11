#include "Volume.hpp"

#include "io/JsonObject.hpp"
#include "io/Scene.hpp"

namespace Tungsten {

struct VolumeIntersection {
	Vec3f p;
	float tmin, tmax;
};

/// Constructor
Volume::Volume()
:	_pos(0.0f),
	_scale(1.0f)	/// TODO
	{}

/// Functions
void Volume::fromJson(JsonPtr value, const Scene &scene) {
	Primitive::fromJson(value, scene);
	//TODO: attach data
	if (auto path = value["data_file"]) _path = scene.fetchResource(path);

	if (auto medium = value["medium"]) {
		_medium_ptr = scene.fetchMedium(medium);
	}

	// bsdf
	if (auto bsdf = value["bsdf"]) _bsdf = scene.fetchBsdf(bsdf);
}

// TODO: check
bool Volume::intersect(Ray &ray, IntersectionTemporary &data) const {
	Vec3f p = ray.pos() - _pos;
	Vec3f d = ray.dir();

	Vec3f invD = 1.0f / d;
	Vec3f relMin((-_scale - p));
	Vec3f relMax((_scale - p));

	float ttMin = ray.nearT(), ttMax = ray.farT();
	for (int i = 0; i < 3; ++i) {
		if (invD[i] >= 0.0f) {
			ttMin = max(ttMin, relMin[i] * invD[i]);
			ttMax = min(ttMax, relMax[i] * invD[i]);
		}
		else {
			ttMax = min(ttMax, relMin[i] * invD[i]);
			ttMin = max(ttMin, relMax[i] * invD[i]);
		}
	}

	/// TODO
	if (ttMin <= ttMax) {
		data.primitive = this;
		// no need?
		// VolumeIntersection *isect = data.as<VolumeIntersection>();
		/// ???
		if (ttMin > ray.nearT() && ttMin < ray.farT()) {
			ray.setFarT(ttMin);
			ray.setNearT(ttMax);
		}
		else if (ttMax > ray.nearT() && ttMax < ray.farT()) {
			ray.setFarT(ttMax);
			ray.setNearT(ttMin);
		}
		return true;
	}
	return false;
}

/// TODO
bool Volume::occluded(const Ray &ray) const {
	return false;
}

bool Volume::hitBackside(const IntersectionTemporary &data) const { return false; }
void Volume::intersectionInfo(const IntersectionTemporary &data, IntersectionInfo &info) const {}
/// TODO
bool Volume::tangentSpace(const IntersectionTemporary &data, const IntersectionInfo &info,
	Vec3f &T, Vec3f &B) const { return false; }

bool Volume::isSamplable() const { return true; }
void Volume::makeSamplable(const TraceableScene &scene, uint32 threadIndex) {}

bool Volume::invertParametrization(Vec2f uv, Vec3f &pos) const { return false; }

bool Volume::isDirac() const { return false; }
bool Volume::isInfinite() const { return false; }

float Volume::approximateRadiance(uint32 threadIndex, const Vec3f &p) const { return 0.0f; }

Box3f Volume::bounds() const {
	Box3f box;	// 3d box represented in float
	for (int i = 0; i < 8; i++) {
		box.grow(_pos + Vec3f(
			(i & 1 ? _scale.x() : -_scale.x()),
			(i & 2 ? _scale.y() : -_scale.y()),
			(i & 4 ? _scale.z() : -_scale.z())));
	}
	return box;
}

/// Volume is an abstract class!
int Volume::numBsdfs() const { return 0; }

std::shared_ptr<Bsdf> &Volume::bsdf(int index) { return _bsdf; }

void Volume::setBsdf(int index, std::shared_ptr<Bsdf> &bsdf) { _bsdf = bsdf; }

Primitive *Volume::clone() {
	return new Volume(*this);	// placement new
}


}