#include "SmokeMedium.hpp"

#include "io/JsonObject.hpp"

#include "math/TangentFrame.hpp"
#include "math/Ray.hpp"

#include "sampling/PathSampleGenerator.hpp"

namespace Tungsten {

/// Constructor
SmokeMedium::SmokeMedium() : _data(nullptr), _absorptionOnly(true) {}


/// Member Function
void SmokeMedium::bindData(float *d) { _data = d; }


/// Override Function
void SmokeMedium::fromJson(JsonPtr value, const Scene &scene) {
	/// Medium::fromJson
	/// fetch [name], [phase_function], [max_bounce]
	Medium::fromJson(value, scene);
	value.getField("sigma_a", _materialSigmaA);
	value.getField("sigma_s", _materialSigmaS);
}

// rapidjson::Value SmokeMedium::toJson(Allocator &allocator) const {}

bool SmokeMedium::isHomogeneous() const { return false; }

void SmokeMedium::prepareForRender() {
	/// called after data is binded!
}


Vec3f SmokeMedium::sigmaA(Vec3f p) const { return density(p) * _sigmaA; }
Vec3f SmokeMedium::sigmaS(Vec3f p) const { return density(p) * _sigmaS; }
Vec3f SmokeMedium::sigmaT(Vec3f p) const { return density(p) * _sigmaT; }

bool SmokeMedium::sampleDistance(PathSampleGenerator &sampler, const Ray &ray,
	MediumState &state, MediumSample &sample) const { return true; }
bool SmokeMedium::invertDistance(WritablePathSampleGenerator &sampler, const Ray &ray, bool onSurface) const { return true; }
Vec3f SmokeMedium::transmittance(PathSampleGenerator &sampler, const Ray &ray) const { return Vec3f(0.0); }
float SmokeMedium::pdf(PathSampleGenerator &sampler, const Ray &ray, bool onSurface) const { return 0.0; }
Vec3f SmokeMedium::transmittanceAndPdfs(PathSampleGenerator &sampler, const Ray &ray, bool startOnSurface,
	bool endOnSurface, float &pdfForward, float &pdfBackward) const { return Vec3f(0.0); }
// const PhaseFunction *SmokeMedium::phaseFunction(const Vec3f &p) const {}

/// Private Function
float SmokeMedium::density(Vec3f p) const { return 0.0; }

}