#include "VoxelMedium.hpp"

#include <fstream>

#include "sampling/PathSampleGenerator.hpp"

#include "math/TangentFrame.hpp"
#include "math/Ray.hpp"

#include "io/JsonObject.hpp"
#include "io/Scene.hpp"

#include "Common.h"

namespace Tungsten {

VoxelMedium::VoxelMedium()
: _sigmaA(0.0f),
  _sigmaS(0.0f)
{
}

void VoxelMedium::fromJson(JsonPtr value, const Scene &scene)
{
	Medium::fromJson(value, scene);
	value.getField("sigma_a", _sigmaA);
	value.getField("sigma_s", _sigmaS);
	_grid = scene.fetchGrid(value.getRequiredMember("grid"));
}

rapidjson::Value VoxelMedium::toJson(Allocator &allocator) const
{
	return JsonObject{Medium::toJson(allocator), allocator,
		"type", "voxel",
		"sigma_a", _sigmaA,
		"sigma_s", _sigmaS,
		"grid", *_grid
	};
}

void VoxelMedium::loadResources()
{
	_grid->loadResources();
}

///TODO
bool VoxelMedium::isHomogeneous() const
{
	return true;	// should be hetergeneous
}

void VoxelMedium::prepareForRender()
{
	_sigmaT = _sigmaA + _sigmaS;
	_absorptionOnly = _sigmaS == 0.0f;

	_worldToGrid = _grid->invNaturalTransform();
	_gridBounds = _grid->bounds();

	std::cout << "world to grid: " << _worldToGrid << std::endl;
	//std::cout << _gridBounds << std::endl;
	//std::cout << Box3f(_grid->naturalTransform()*_gridBounds.min(),
					   //_grid->naturalTransform()*_gridBounds.max()) << std::endl;
}

static inline bool bboxIntersection(const Box3f &box, const Vec3f &o, const Vec3f &d,
		float &tMin, float &tMax)
{
	Vec3f invD = 1.0f/d;
	Vec3f relMin((box.min() - o));
	Vec3f relMax((box.max() - o));

	float ttMin = tMin, ttMax = tMax;
	for (int i = 0; i < 3; ++i) {
		if (invD[i] >= 0.0f) {
			ttMin = max(ttMin, relMin[i]*invD[i]);
			ttMax = min(ttMax, relMax[i]*invD[i]);
		} else {
			ttMax = min(ttMax, relMin[i]*invD[i]);
			ttMin = max(ttMin, relMax[i]*invD[i]);
		}
	}

	if (ttMin <= ttMax) {
		tMin = ttMin;
		tMax = ttMax;
		return true;
	}
	return false;
}

Vec3f VoxelMedium::sigmaA(Vec3f p) const
{
#if USE_DENSITY
	return _grid->density(p)*_sigmaA;
#else
	return _grid->gradient(p)*_sigmaA;
#endif
}

Vec3f VoxelMedium::sigmaS(Vec3f p) const
{
#if USE_DENSITY
	return _grid->density(p)*_sigmaS;
#else
	return _grid->gradient(p)*_sigmaS;
#endif
}

Vec3f VoxelMedium::sigmaT(Vec3f p) const
{
#if USE_DENSITY
	return _grid->density(p)*_sigmaT;
#else
	return _grid->gradient(p)*_sigmaT;
#endif
}

/* samples a medium scattering interaction along the ray between [0, tmax] 
 * sample the integral form of the equation of transfer
 * pbrt-v3 section 15.2
 * 
 * ! use DELTA-TRACKING	for heterogeneous medium
 *
 * should set     t, weight, pdf
 */
bool VoxelMedium::sampleDistance(PathSampleGenerator &sampler, const Ray &ray,
		MediumState &state, MediumSample &sample) const
{
	if (state.bounce > _maxBounce)
		return false;

#if GRID_NORMAL_TRANSFORM
	// transform the ray into medium coordinate system
	// and normalizing the ray direction [pbrt-v3 p896]
	float maxT = ray.farT();	// as tmax for bboxIntersection
	Vec3f p = _worldToGrid*ray.pos();	// local space point coordinate
	Vec3f w = _worldToGrid.transformVector(ray.dir());	// transform ray direction to grid coordinate
	float wPrime = w.length();
	w /= wPrime;	// normalize ray dir
	float t0 = 0.0f, t1 = maxT*wPrime;
#else
	float maxT = ray.farT();
	Vec3f p = _worldToGrid*ray.pos();
	Vec3f w = _worldToGrid.transformVector(ray.dir());
	float wPrime = 1.0f;
	float t0 = 0.0f; float t1 = maxT;
#endif

	// compute parametric range
	if (!bboxIntersection(_gridBounds, p/*ray origin*/, w/*ray dir*/, t0, t1)) {	// t0, t1 get updated to box intersection point
		/// do not intersect volume, no absorption and scatter and emission
		sample.t = maxT;
		sample.weight = Vec3f(1.0f);
		sample.pdf = 1.0f;
		sample.exited = true;
		return true;
	}


	/*std::ofstream f("intersection", std::ios::out | std::ios::app);
	if (f.good()) {
		f << _grid->naturalTransform()*(p + t0*w);
		f << " ";
		f << _grid->naturalTransform()*(p + t1*w);
		f << " ";
	}
	f.close();*/

	if (_absorptionOnly) {
		/// if absorption only, no integral in RTE need to be solved!
		/// we only need to solve the transmittance integral
		sample.t = maxT;
		sample.weight = _grid->transmittance(sampler, p, w, t0, t1, _sigmaT);	// transparency
		sample.pdf = 1.0f;
		sample.exited = true;

		/*float t = t0 + _grid->stepSize();
		if (t > t1) {
			sample.t = maxT;
			sample.weight = Vec3f(1.0);
			sample.pdf = 1.0f;
			sample.exited = true;
		}
		else {
			sample.weight = _grid->transmittance(sampler, p, w, sample.t, t, _sigmaT);
			sample.t = t;
			sample.pdf = 1.0 / t1 - t0;
			sample.exited = false;
		}*/

	} else {
		int component = sampler.nextDiscrete(3);
		float sigmaTc = _sigmaT[component];

		float gridValue;
		float maxGridValue;
#if USE_DENSITY
		maxGridValue = _grid->maxDensity();
#else
		maxGridValue = _grid->maxGradient()[component];
#endif

		float maxSigmaT = sigmaTc * maxGridValue;
		float sigmaN = maxSigmaT - sigmaTc;

		float t = t0;
		sample.pdf = 1.0f;
		sample.weight = Vec3f(1.0f);

		while (true) {
			t -= std::log(1 - sampler.next1D()) / maxSigmaT;
#if USE_DENSITY
			gridValue = _grid->density(p + t*w);
#else
#if GRID_NORMAL_TRANSFORM
			gridValue = _grid->gradient(p + t*w)[component];
#else
			gridValue = _grid->gradient(_grid->invNaturalTransform() * (p + t*w))[component];
#endif
#endif
			// std::cout << _grid->density(p + t*w) << ", " << _grid->gradient(p + t*w) << std::endl;

			if (t > t1) {
				sample.t = t1;
				sample.exited = true;
				sample.pdf *= _grid->transmittance(sampler, p, w, t0, t1, _sigmaT).avg(); // transmittance
				sample.weight = Vec3f(1.0f);	// beta
				break;
			}

			float zeta = sampler.next1D();
			if (zeta < (_sigmaA*gridValue / maxSigmaT).avg()) {
				// absorb
				sample.exited = false;
				sample.weight = _sigmaA / sigmaTc;
				sample.t = t;
				break;
			}
			else if (zeta < (gridValue / maxGridValue) ) {
				// scatter 
				sample.t = t;
				sample.exited = false;
				sample.weight = _sigmaS / _sigmaT;
				break;
			}

			/*if ((_grid->density(p + t*w) / _grid->maxDensity()) > sampler.next1D()) {
				sample.t = t;
				sample.exited = false;
				sample.weight = _sigmaS / _sigmaT;
				break;
			}*/
		}
		sample.t /= wPrime;
		state.advance();
	}
	sample.p = ray.pos() + sample.t*ray.dir();
	sample.phase = _phaseFunction.get();

	return true;
}

Vec3f VoxelMedium::transmittance(PathSampleGenerator &sampler, const Ray &ray) const
{
#if GRID_NORMAL_TRANSFORM
	Vec3f p = _worldToGrid*ray.pos();
	Vec3f w = _worldToGrid.transformVector(ray.dir());
	float wPrime = w.length();
	w /= wPrime;
	float t0 = 0.0f, t1 = ray.farT()*wPrime;
#else
	float maxT = ray.farT();
	Vec3f p = _worldToGrid*ray.pos();
	Vec3f w = _worldToGrid.transformVector(ray.dir());
	float wPrime = 1.0f;
	float t0 = 0.0f; float t1 = maxT;
#endif

	if (!bboxIntersection(_gridBounds, p, w, t0, t1))
		return Vec3f(1.0f);

	return _grid->transmittance(sampler, p, w, t0, t1, _sigmaT);
}

float VoxelMedium::pdf(PathSampleGenerator &sampler, const Ray &ray, bool onSurface) const
{
	if (_absorptionOnly) {
		return 1.0f;
	} else {
		Vec3f p = _worldToGrid*ray.pos();
		Vec3f w = _worldToGrid.transformVector(ray.dir());
		float wPrime = w.length();
		w /= wPrime;
		float t0 = 0.0f, t1 = ray.farT()*wPrime;
		if (!bboxIntersection(_gridBounds, p, w, t0, t1))
			return 1.0f;

		Vec3f transmittance = _grid->transmittance(sampler, p, w, t0, t1, _sigmaT/wPrime);
		if (onSurface) {
			return transmittance.avg();
		} else {
			return (_grid->density(p)*_sigmaT*transmittance).avg();
		}
	}
}

/* transmittance() & pdf()
 * called from: TraceBase::generalizedShadowRayImpl
*/
Vec3f VoxelMedium::transmittanceAndPdfs(PathSampleGenerator &sampler, const Ray &ray, bool startOnSurface,
		bool endOnSurface, float &pdfForward, float &pdfBackward) const
{
#if	GRID_NORMAL_TRANSFORM
	Vec3f p = _worldToGrid*ray.pos();
	Vec3f w = _worldToGrid.transformVector(ray.dir());
	float wPrime = w.length();
	w /= wPrime;
	float t0 = 0.0f, t1 = ray.farT()*wPrime;
#else
	float maxT = ray.farT();
	Vec3f p = _worldToGrid*ray.pos();
	Vec3f w = _worldToGrid.transformVector(ray.dir());
	float wPrime = 1.0f;
	float t0 = 0.0f; float t1 = maxT;
#endif

	if (!bboxIntersection(_gridBounds, p, w, t0, t1)) {
		pdfForward = pdfBackward = 1.0f;
		return Vec3f(1.0f);
	}

	Vec3f transmittance = _grid->transmittance(sampler, p, w, t0, t1, _sigmaT/wPrime);

	if (_absorptionOnly) {
		pdfForward = pdfBackward = 1.0f;
	} else {
		pdfForward  =   endOnSurface ? transmittance.avg() : (_grid->density(p)*_sigmaT*transmittance).avg();
		pdfBackward = startOnSurface ? transmittance.avg() : (_grid->density(p)*_sigmaT*transmittance).avg();
	}

	return transmittance;
}

}
