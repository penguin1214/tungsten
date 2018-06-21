#include "DielectricBsdf.hpp"
#include "Fresnel.hpp"

#include "samplerecords/SurfaceScatterEvent.hpp"

#include "sampling/PathSampleGenerator.hpp"
#include "sampling/SampleWarp.hpp"

#include "math/MathUtil.hpp"
#include "math/Angle.hpp"
#include "math/Vec.hpp"

#include "io/JsonObject.hpp"

#include <rapidjson/document.h>

#define DEBUG_DIELECTRIC_BSDF 0

namespace Tungsten {

DielectricBsdf::DielectricBsdf()
: _ior(1.5f),
  _enableT(true)
{
    _lobes = BsdfLobes(BsdfLobes::SpecularReflectionLobe | BsdfLobes::SpecularTransmissionLobe);
}

DielectricBsdf::DielectricBsdf(float ior)
: _ior(ior),
  _enableT(true)
{
    _lobes = BsdfLobes(BsdfLobes::SpecularReflectionLobe | BsdfLobes::SpecularTransmissionLobe);
}

void DielectricBsdf::fromJson(JsonPtr value, const Scene &scene)
{
    Bsdf::fromJson(value, scene);
    value.getField("ior", _ior);
    value.getField("enable_refraction", _enableT);
}

rapidjson::Value DielectricBsdf::toJson(Allocator &allocator) const
{
    return JsonObject{Bsdf::toJson(allocator), allocator,
        "type", "dielectric",
        "ior", _ior,
        "enable_refraction", _enableT
    };
}

/* sample output direction. */
bool DielectricBsdf::sample(SurfaceScatterEvent &event) const
{
    bool sampleR = event.requestedLobe.test(BsdfLobes::SpecularReflectionLobe);
    bool sampleT = event.requestedLobe.test(BsdfLobes::SpecularTransmissionLobe) && _enableT;

    float eta = event.wi.z() < 0.0f ? _ior : _invIor;	// event.wi.z() < 0 => outgoing
#if DEBUG_DIELECTRIC_BSDF
	std::cout << "hitpoint frame normal: " << event.frame.normal << std::endl;
	std::cout << "ior: " << eta << std::endl;
#endif

    float cosThetaT = 0.0f;
#if DEBUG_DIELECTRIC_BSDF
	std::cout << "event.wi.z(): " << std::abs(event.wi.z()) << std::endl;
#endif
    float F = Fresnel::dielectricReflectance(eta, std::abs(event.wi.z())/*cosThetaI*/, cosThetaT);
#if DEBUG_DIELECTRIC_BSDF
	std::cout << "F: " << F << std::endl;
#endif
    float reflectionProbability;
    if (sampleR && sampleT)
        reflectionProbability = F;
    else if (sampleR)
        reflectionProbability = 1.0f;
    else if (sampleT)
        reflectionProbability = 0.0f;
    else
        return false;

    if (event.sampler->nextBoolean(reflectionProbability)) {	// reflection
#if DEBUG_DIELECTRIC_BSDF
		std::cout << "reflect" << std::endl;
#endif
        event.wo = Vec3f(-event.wi.x(), -event.wi.y(), event.wi.z());
        event.pdf = reflectionProbability;
        event.sampledLobe = BsdfLobes::SpecularReflectionLobe;
        event.weight = sampleT ? Vec3f(1.0f) : Vec3f(F);
    } else {	// refraction
#if DEBUG_DIELECTRIC_BSDF
		std::cout << "refract" << std::endl;
#endif
        if (F == 1.0f)
            return false;
		
        event.wo = Vec3f(-event.wi.x()*eta, -event.wi.y()*eta, -std::copysign(cosThetaT, event.wi.z()));	// refraction direction
#if DEBUG_DIELECTRIC_BSDF
		std::cout << "frame local wo: " << event.wo << std::endl;
#endif
		/*std::cout << "wi: " << event.wi << std::endl;
		std::cout << "wo: " << event.wo << std::endl;
		std::cout << std::abs(event.wo.x() / event.wo.z()) << std::endl;
		std::cout << cosThetaT << std::endl;*/
        event.pdf = 1.0f - reflectionProbability;
        event.sampledLobe = BsdfLobes::SpecularTransmissionLobe;
        event.weight = sampleR ? Vec3f(1.0f) : Vec3f(1.0f - F);
    }
    event.weight *= albedo(event.info);

    return true;
}

/* evaluate transparency (not transmittance!) */
Vec3f DielectricBsdf::eval(const SurfaceScatterEvent &event) const
{
    bool evalR = event.requestedLobe.test(BsdfLobes::SpecularReflectionLobe);	// test if reflective
    bool evalT = event.requestedLobe.test(BsdfLobes::SpecularTransmissionLobe) && _enableT;	// test if transmissive


    float eta = event.wi.z() < 0.0f ? _ior : _invIor;
    float cosThetaT = 0.0f;
    float F = Fresnel::dielectricReflectance(eta, std::abs(event.wi.z()), cosThetaT);

    if (event.wi.z()*event.wo.z() >= 0.0f) {	// event.wo = -wi
		if (evalR && checkReflectionConstraint(event.wi, event.wo)) {
			//std::cout << "checkreflectionconstraint" << std::endl;
            return F*albedo(event.info);
		}
		else {
			return Vec3f(0.0f);
		}
    } else {
		if (evalT && checkRefractionConstraint(event.wi, event.wo, eta, cosThetaT)) {
			//std::cout << "check refraction constraint" << std::endl;
			return (1.0f - F)*albedo(event.info);
		}
		else {
			//std::cout << "transmit no constraint" << std::endl;
			return Vec3f(0.0f);
		}
    }
}

bool DielectricBsdf::invert(WritablePathSampleGenerator &sampler, const SurfaceScatterEvent &event) const
{
    bool evalR = event.requestedLobe.test(BsdfLobes::SpecularReflectionLobe);
    bool evalT = event.requestedLobe.test(BsdfLobes::SpecularTransmissionLobe) && _enableT;

    float eta = event.wi.z() < 0.0f ? _ior : _invIor;
    float cosThetaT = 0.0f;
    float F = Fresnel::dielectricReflectance(eta, std::abs(event.wi.z()), cosThetaT);

    float reflectionProbability;
    if (evalR && evalT)
        reflectionProbability = F;
    else if (evalR)
        reflectionProbability = 1.0f;
    else if (evalT)
        reflectionProbability = 0.0f;
    else
        return false;

    if (event.wi.z()*event.wo.z() >= 0.0f) {
        if (evalR && checkReflectionConstraint(event.wi, event.wo)) {
            sampler.putBoolean(reflectionProbability, true);
            return true;
        } else {
            return false;
        }
    } else {
        if (evalT && checkRefractionConstraint(event.wi, event.wo, eta, cosThetaT)) {
            sampler.putBoolean(reflectionProbability, false);
            return true;
        } else {
            return false;
        }
    }
}

float DielectricBsdf::pdf(const SurfaceScatterEvent &event) const
{
    bool sampleR = event.requestedLobe.test(BsdfLobes::SpecularReflectionLobe);
    bool sampleT = event.requestedLobe.test(BsdfLobes::SpecularTransmissionLobe) && _enableT;

    float eta = event.wi.z() < 0.0f ? _ior : _invIor;
    float cosThetaT = 0.0f;
    float F = Fresnel::dielectricReflectance(eta, std::abs(event.wi.z()), cosThetaT);

    if (event.wi.z()*event.wo.z() >= 0.0f) {
        if (sampleR && checkReflectionConstraint(event.wi, event.wo))
            return sampleT ? F : 1.0f;
        else
            return 0.0f;
    } else {
        if (sampleT && checkRefractionConstraint(event.wi, event.wo, eta, cosThetaT))
            return sampleR ? 1.0f - F : 1.0f;
        else
            return 0.0f;
    }
}

float DielectricBsdf::eta(const SurfaceScatterEvent &event) const
{
    if (event.wi.z()*event.wo.z() >= 0.0f)
        return 1.0f;
    else
        return event.wi.z() < 0.0f ? _ior : _invIor;
}

void DielectricBsdf::prepareForRender()
{
    if (_enableT)
        _lobes = BsdfLobes(BsdfLobes::SpecularReflectionLobe | BsdfLobes::SpecularTransmissionLobe);
    else
        _lobes = BsdfLobes(BsdfLobes::SpecularReflectionLobe);

    _invIor = 1.0f/_ior;
}

}
