#include "PathTracer.hpp"

#include "bsdfs/TransparencyBsdf.hpp"

namespace Tungsten {

PathTracer::PathTracer(TraceableScene *scene, const PathTracerSettings &settings, uint32 threadId)
: TraceBase(scene, settings, threadId),
  _settings(settings),
  _trackOutputValues(!scene->rendererSettings().renderOutputs().empty())
{
}

Vec3f PathTracer::traceSample(Vec2u pixel, PathSampleGenerator &sampler)
{
	// TODO: Put diagnostic colors in JSON? -dev
	const Vec3f nanDirColor = Vec3f(0.0f);
	const Vec3f nanEnvDirColor = Vec3f(0.0f);
	const Vec3f nanBsdfColor = Vec3f(0.0f);

	try {

	PositionSample point;
	if (!_scene->cam().samplePosition(sampler, point))
		return Vec3f(0.0f);
	DirectionSample direction;
	if (!_scene->cam().sampleDirection(sampler, point, pixel, direction))
		return Vec3f(0.0f);

	Vec3f throughput/*transparency*/ = point.weight*direction.weight;	/// ???
	Ray ray(point.p, direction.d);
	ray.setPrimaryRay(true);	// primary ray, no bounce yet

	MediumSample mediumSample;
	SurfaceScatterEvent surfaceEvent;
	IntersectionTemporary data;
	Medium::MediumState state;
	state.reset();
	IntersectionInfo info;
	Vec3f emission(0.0f);
	/// What if multiple medium?
	const Medium *medium = _scene->cam().medium().get();	// what for????

	bool recordedOutputValues = false;
	float hitDistance = 0.0f;

	// cast ray
	int bounce = 0;
	bool didHit = _scene->intersect(ray, data, info);
	bool wasSpecular = true;
	while ((didHit || medium) && bounce < _settings.maxBounces) {
		bool hitSurface = true;	// whether exited volume
		if (medium) {
			if (!medium->sampleDistance(sampler, ray, state, mediumSample))
				return emission;
			throughput *= mediumSample.weight;
			hitSurface = mediumSample.exited;	/// !
			if (hitSurface && !didHit)
				break;
		}

		/// do scattering
		/// if not exited volume yet, do volume scattering, else do surface scattering
		if (hitSurface) {
			/// why farT()
			hitDistance += ray.farT();	// hitDistance is 0 for primary ray

			surfaceEvent = makeLocalScatterEvent(data, info, ray, &sampler);
			Vec3f transmittance(-1.0f);
			bool terminate = !handleSurface(surfaceEvent, data, info, medium, bounce, false,
					_settings.enableLightSampling, ray, throughput, emission, wasSpecular, state, &transmittance);

			/// ???
			if (_trackOutputValues && !recordedOutputValues && (!wasSpecular || terminate)) {
				if (_scene->cam().depthBuffer())
					_scene->cam().depthBuffer()->addSample(pixel, hitDistance);
				if (_scene->cam().normalBuffer())
					_scene->cam().normalBuffer()->addSample(pixel, info.Ns);
				if (_scene->cam().albedoBuffer()) {
					Vec3f albedo;
					if (const TransparencyBsdf *bsdf = dynamic_cast<const TransparencyBsdf *>(info.bsdf))
						albedo = (*bsdf->base()->albedo())[info];
					else
						albedo = (*info.bsdf->albedo())[info];
					if (info.primitive->isEmissive())
						albedo += info.primitive->evalDirect(data, info);
					_scene->cam().albedoBuffer()->addSample(pixel, albedo);
				}
				if (_scene->cam().visibilityBuffer() && transmittance != -1.0f)
					_scene->cam().visibilityBuffer()->addSample(pixel, transmittance.avg());
				recordedOutputValues = true;
			}

			if (terminate)
				return emission;
		} else {
			if (!handleVolume(sampler, mediumSample, medium, bounce, false,
					_settings.enableVolumeLightSampling, ray, throughput, emission, wasSpecular))
				return emission;
		}

		if (throughput.max() == 0.0f)	// nothing left after absorption
			break;

		float roulettePdf = std::abs(throughput).max();	/// ???
		if (bounce > 2 && roulettePdf < 0.1f) {
			if (sampler.nextBoolean(roulettePdf))
				throughput /= roulettePdf;	// bounce!
			else
				return emission;
		}

		if (std::isnan(ray.dir().sum() + ray.pos().sum()))
			return nanDirColor;
		if (std::isnan(throughput.sum() + emission.sum()))
			return nanBsdfColor;

		bounce++;
		if (bounce < _settings.maxBounces)
			didHit = _scene->intersect(ray, data, info);
		// recursive bounce
	}

	if (bounce >= _settings.minBounces && bounce < _settings.maxBounces)
		/// ???
		handleInfiniteLights(data, info, _settings.enableLightSampling, ray, throughput, wasSpecular, emission);
	if (std::isnan(throughput.sum() + emission.sum()))
		return nanEnvDirColor;

	if (_trackOutputValues && !recordedOutputValues) {
		if (_scene->cam().depthBuffer() && bounce == 0)
			_scene->cam().depthBuffer()->addSample(pixel, 0.0f);
		if (_scene->cam().normalBuffer())
			_scene->cam().normalBuffer()->addSample(pixel, -ray.dir());
		if (_scene->cam().albedoBuffer() && info.primitive && info.primitive->isInfinite())
			_scene->cam().albedoBuffer()->addSample(pixel, info.primitive->evalDirect(data, info));
	}

	return emission;

	} catch (std::runtime_error &e) {
		std::cout << tfm::format("Caught an internal error at pixel %s: %s", pixel, e.what()) << std::endl;

		return Vec3f(0.0f);
	}
}

}
