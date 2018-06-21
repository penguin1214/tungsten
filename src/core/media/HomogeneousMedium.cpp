#include "HomogeneousMedium.hpp"

#include "sampling/PathSampleGenerator.hpp"

#include "math/TangentFrame.hpp"
#include "math/FastMath.hpp"
#include "math/Ray.hpp"

#include "io/JsonObject.hpp"

namespace Tungsten {

HomogeneousMedium::HomogeneousMedium()
: _materialSigmaA(0.0f),
  _materialSigmaS(0.0f),
  _density(1.0f)
{
}

void HomogeneousMedium::fromJson(JsonPtr value, const Scene &scene)
{
    Medium::fromJson(value, scene);
    value.getField("sigma_a", _materialSigmaA);
    value.getField("sigma_s", _materialSigmaS);
    value.getField("density", _density);
}

rapidjson::Value HomogeneousMedium::toJson(Allocator &allocator) const
{
    return JsonObject{Medium::toJson(allocator), allocator,
        "type", "homogeneous",
        "sigma_a", _materialSigmaA,
        "sigma_s", _materialSigmaS,
        "density", _density
    };
}

bool HomogeneousMedium::isHomogeneous() const
{
    return true;
}

void HomogeneousMedium::prepareForRender()
{
    _sigmaA = _materialSigmaA*_density;
    _sigmaS = _materialSigmaS*_density;
    _sigmaT = _sigmaA + _sigmaS;
    _absorptionOnly = _sigmaS == 0.0f;
}

Vec3f HomogeneousMedium::sigmaA(Vec3f /*p*/) const
{
    return _sigmaA;
}

Vec3f HomogeneousMedium::sigmaS(Vec3f /*p*/) const
{
    return _sigmaS;
}

Vec3f HomogeneousMedium::sigmaT(Vec3f /*p*/) const
{
    return _sigmaT;
}

bool HomogeneousMedium::sampleDistance(PathSampleGenerator &sampler, const Ray &ray,
        MediumState &state, MediumSample &sample) const
{
    if (state.bounce > _maxBounce)
        return false;

    float maxT = ray.farT();
    if (_absorptionOnly) {
		/*if (maxT == Ray::infinity())
			return false;
		sample.t = maxT;
		sample.weight = FastMath::exp(-_sigmaT*maxT);
		sample.pdf = 1.0f;
		sample.exited = true;*/

		float t = 0;
		bool terminate = false;

		do {
			float zeta = sampler.next1D();
			t -= std::log(1 - zeta);
			sample.exited = (t > maxT);

			if (t > maxT) {
				terminate = true;
				sample.t = maxT;
				sample.weight = FastMath::exp(-_sigmaT*maxT);
				sample.pdf = 1.0f;
				break;
			}
			else {
				// absorb
				terminate = true;
				sample.t = t;
				sample.weight = _sigmaT * FastMath::exp(-_sigmaT*t);
				sample.pdf *= sample.weight.avg();
			}
			
		} while (!terminate);

		state.advance();

    } else {
		// sample a spectrum channel to pick a particular extinction coefficient
        int component = sampler.nextDiscrete(3);
        float sigmaTc = _sigmaT[component];

		// importance sampling the exponential distribution (of Tr),
		// as the interaction point
        float t = -std::log(1.0f - sampler.next1D())/sigmaTc;
        sample.t = min(t, maxT);	// test whether the sampled point is in the interval or on surface

        sample.continuedT = t;

        sample.weight = FastMath::exp(-sample.t*_sigmaT);	// Tr
        sample.continuedWeight = FastMath::exp(-sample.continuedT*_sigmaT);

        sample.exited = (t >= maxT);
        if (sample.exited) {
            sample.pdf = sample.weight.avg();
        } else {
            sample.pdf = (_sigmaT*sample.weight).avg();	// _sigmaT is due to MC integration
            sample.weight *= _sigmaS;
        }
        sample.weight /= sample.pdf;
        sample.continuedWeight = _sigmaS*sample.continuedWeight/(_sigmaT*sample.continuedWeight).avg();

        state.advance();
    }
    sample.p = ray.pos() + sample.t*ray.dir();
    sample.phase = _phaseFunction.get();

    return true;
}

bool HomogeneousMedium::invertDistance(WritablePathSampleGenerator &sampler, const Ray &ray, bool onSurface) const
{
    float maxT = ray.farT();
    if (_absorptionOnly) {
        if (maxT == Ray::infinity())
            return false;
        return true;
    } else {
        Vec3f Tr = std::exp(-_sigmaT*maxT);
        Vec3f pdfs = onSurface ? Tr : _sigmaT*Tr;
        Vec3f cdf(pdfs.x(), pdfs.x() + pdfs.y(), pdfs.sum());

        float target = sampler.untracked1D()*cdf.z();
        int component = target < cdf.x() ? 0 : (target < cdf.y() ? 1 : 2);
        float Trc = Tr[component];

        float xi = 1.0f - Trc;
        if (onSurface)
            xi -= (1.0f - Trc)*sampler.next1D();

        sampler.putDiscrete(3, component);
        sampler.put1D(xi);
    }
    return true;
}

Vec3f HomogeneousMedium::transmittance(PathSampleGenerator &/*sampler*/, const Ray &ray) const
{
    if (ray.farT() == Ray::infinity())
        return Vec3f(0.0f);
    else {
        return FastMath::exp(-_sigmaT*ray.farT());
    }
}

float HomogeneousMedium::pdf(PathSampleGenerator &/*sampler*/, const Ray &ray, bool onSurface) const
{
    if (_absorptionOnly) {
        return 1.0f;
    } else {
        if (onSurface)
            return FastMath::exp(-ray.farT()*_sigmaT).avg();
        else
            return (_sigmaT*FastMath::exp(-ray.farT()*_sigmaT)).avg();
    }
}

Vec3f HomogeneousMedium::transmittanceAndPdfs(PathSampleGenerator &/*sampler*/, const Ray &ray, bool startOnSurface,
        bool endOnSurface, float &pdfForward, float &pdfBackward) const
{
    if (ray.farT() == Ray::infinity()) {
        pdfForward = pdfBackward = 0.0f;
        return Vec3f(0.0f);
    } else if (_absorptionOnly) {
        pdfForward = pdfBackward = 1.0f;
        return FastMath::exp(-_sigmaT*ray.farT());
    } else {
        Vec3f weight = FastMath::exp(-_sigmaT*ray.farT());
        pdfForward  =   endOnSurface ? weight.avg() : (_sigmaT*weight).avg();
        pdfBackward = startOnSurface ? weight.avg() : (_sigmaT*weight).avg();
        return weight;
    }
}

bool HomogeneousMedium::invert(WritablePathSampleGenerator &sampler, const Ray &ray, bool onSurface) const
{
    if (_absorptionOnly)
        return true;

    Vec3f transmittance = std::exp(-ray.farT()*_sigmaT);
    Vec3f pdfs = _sigmaT*transmittance;
    float target = sampler.untracked1D()*pdfs.sum();
    int component = (target < pdfs.x() ? 0 : (target < pdfs.x() + pdfs.y() ? 1 : 2));

    float xi = 1.0f - transmittance[component];
    if (onSurface)
        xi += (1.0f - xi)*sampler.untracked1D();
    sampler.putDiscrete(3, component);
    sampler.put1D(xi);

    return true;
}

}
