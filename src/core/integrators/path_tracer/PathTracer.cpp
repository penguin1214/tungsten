#include "PathTracer.hpp"

#include <fstream>
#include <iostream>

#include "bsdfs/TransparencyBsdf.hpp"

#define ENABLE_INTERSECTION_TEST 0
#define DEBUG_PATH_TRACER 0
#define OUTPUT_THROUGHPUT 0

namespace Tungsten {

	PathTracer::PathTracer(TraceableScene *scene, const PathTracerSettings &settings, uint32 threadId)
		: TraceBase(scene, settings, threadId),
		_settings(settings),
		_trackOutputValues(!scene->rendererSettings().renderOutputs().empty()) {
	}

	Vec3f PathTracer::traceSample(Vec2u pixel, PathSampleGenerator &sampler) {
		// TODO: Put diagnostic colors in JSON? -dev
		const Vec3f nanDirColor = Vec3f(0.0f);
		const Vec3f nanEnvDirColor = Vec3f(0.0f);
		const Vec3f nanBsdfColor = Vec3f(0.0f);

		std::fstream f;

		try {
			PositionSample point;
			if (!_scene->cam().samplePosition(sampler, point))
				return Vec3f(0.0f);
			DirectionSample direction;
			if (!_scene->cam().sampleDirection(sampler, point, pixel, direction))
				return Vec3f(0.0f);

			Vec3f throughput = point.weight*direction.weight;
			Ray ray(point.p, direction.d);

#if ENABLE_INTERSECTION_TEST
			saveIntersect(true, point.p);
#endif

			ray.setPrimaryRay(true);	// primary ray, no bounce yet

			MediumSample mediumSample;
			SurfaceScatterEvent surfaceEvent;
			IntersectionTemporary data;
			Medium::MediumState state;
			state.reset();
			IntersectionInfo info;
			Vec3f emission(0.0f);
			/// What if multiple medium?
			const Medium *medium = _scene->cam().medium().get();

			bool recordedOutputValues = false;
			float hitDistance = 0.0f;

			// cast ray
			int bounce = 0;
			bool didHit = _scene->intersect(ray, data, info);
			// if (didHit) std::cout << "intersect, bounce: " << bounce << std::endl;
			bool wasSpecular = true;

			Vec3f left(0.0f, 1.0f, 0.0f);
			// normal shading
			/*if (didHit)
				emission += info.Ns.dot(left) * 0.1;*/
				//emission += info.Ns.dot(_scene->cam().up()) * 0.1;
			//emission = Vec3f(0.1);
			/*if (didHit)
				emission = Vec3f(info.Ng.avg());*/

#if DEBUG_PATH_TRACER
			if (didHit) {
				std::cout << "=============start tracing==============" << std::endl;
			}
#endif
#if 1
			while ((didHit || medium) && bounce < _settings.maxBounces) {
#if DEBUG_PATH_TRACER
				std::cout << "bounce" << std::endl;
#endif
				bool hitSurface = true;	// whether exited volume
				if (medium) {
					if (!medium->sampleDistance(sampler, ray, state, mediumSample))
						return emission;

					throughput *= mediumSample.weight;
					hitSurface = mediumSample.exited;
					
					if (hitSurface && !didHit)
						break;
				}

				if (hitSurface) {
					hitDistance += ray.farT();	// hitDistance is 0 for primary ray

					surfaceEvent = makeLocalScatterEvent(data, info, ray, &sampler);
					Vec3f transmittance(-1.0f);

#if 0
					std::ofstream f("intersection", std::ios::out | std::ios::app);
					if (f.good()) {
						f << ray.hitpoint() << std::endl;
					}
					f.close();
#endif
#if ENABLE_INTERSECTION_TEST
					saveIntersect(false, ray.hitpoint());
#endif
#if DEBUG_PATH_TRACER
					std::cout << "ray.dir: " << ray.dir() << std::endl;
					std::cout << "hit point geo normal: " << info.Ng << std::endl;
					std::cout << "hit point shading normal: " << info.Ns << std::endl;
					std::cout << "frame Normal: " << surfaceEvent.frame.normal << std::endl;
					std::cout << "frame Tangent: " << surfaceEvent.frame.tangent << std::endl;
					std::cout << "frame Bitangent: " << surfaceEvent.frame.bitangent << std::endl;
					std::cout << "frame local wi: " << surfaceEvent.wi;
					std::cout << std::endl;
#endif
					bool terminate = !handleSurface(surfaceEvent, data, info, medium, bounce, false,
						_settings.enableLightSampling, ray, throughput, emission, wasSpecular, state, &transmittance);
#if OUTPUT_THROUGHPUT
						std::cout << "throughput after surface: " << throughput << std::endl;
#endif

#if 1
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
#endif
					if (terminate) return emission;
				}
				else {
#if ENABLE_INTERSECTION_TEST
					saveIntersect(false, mediumSample.p);
#endif
					if (!handleVolume(sampler, mediumSample, medium, bounce, false,
						_settings.enableVolumeLightSampling, ray, throughput, emission, wasSpecular)) {
#if OUTPUT_THROUGHPUT
						std::cout << "throughput after volume: " << throughput << std::endl;
#endif
						return emission;
					}
				}

				if (throughput.max() < 0.01f)	// nothing left after absorption
					break;

				// russian roulette to decide whether terminate
				float roulettePdf = std::abs(throughput).max();
				if (bounce > 3 && roulettePdf < 0.1f) {
					if (sampler.nextBoolean(roulettePdf))
						throughput /= roulettePdf;
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
			}

			// evaluate direct lighting of infinite lights
			if (bounce >= _settings.minBounces && bounce < _settings.maxBounces)
				handleInfiniteLights(data, info, _settings.enableLightSampling, ray, throughput, wasSpecular, emission);
			if (std::isnan(throughput.sum() + emission.sum()))
				return nanEnvDirColor;
#endif

#if 1
			if (_trackOutputValues && !recordedOutputValues) {
				if (_scene->cam().depthBuffer() && bounce == 0)
					_scene->cam().depthBuffer()->addSample(pixel, 0.0f);
				if (_scene->cam().normalBuffer())
					_scene->cam().normalBuffer()->addSample(pixel, -ray.dir());
				if (_scene->cam().albedoBuffer() && info.primitive && info.primitive->isInfinite())
					_scene->cam().albedoBuffer()->addSample(pixel, info.primitive->evalDirect(data, info));
			}
#endif
			return emission;

		}
		catch (std::runtime_error &e) {
			std::cout << tfm::format("Caught an internal error at pixel %s: %s", pixel, e.what()) << std::endl;

			return Vec3f(0.0f);
		}
	}

	void PathTracer::saveIntersect(bool is_start, Vec3f p) {
		std::ofstream f("intersection", std::ios::out|std::ios::app);
		if (f.good()) {
			if (is_start) {
				// start new line
				f << std::endl << p;
			}
			else {
				f << p;
			}
			f << " ";
		}
		f.close();
	}

	void PathTracer::saveIntersect(bool is_start, Vec3f p, float tr) {
		std::ofstream f("intersection", std::ios::out | std::ios::app);
		if (f.good()) {
			if (is_start) {
				// start new line
				f << std::endl << p;
			}
			else {
				f << p << " " << tr;
			}
			f << " ";
		}
		f.close();
	}

}
