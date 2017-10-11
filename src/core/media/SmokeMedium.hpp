#ifndef SMOKEMEDIUM_HPP_
#define SMOKEMEDIUM_HPP_

#include "Medium.hpp"

namespace Tungsten {

/// inherent
/// @ std::shared_ptr<PhaseFunction> _phaseFunction;
/// @ int _maxBounce;
class SmokeMedium : public Medium {
private:
	float *_data;	// density value
	Vec3f _materialSigmaA, _materialSigmaS;
	Vec3f _sigmaA, _sigmaS;	// density weighted coefficient
	Vec3f _sigmaT;
	bool _absorptionOnly;
	/// TODO : scale factors?

	inline float density(Vec3f p) const;

public:

	/// Constructor
	SmokeMedium();

	/// Member Function
	void bindData(float *d);

	/// Override Function
	virtual void fromJson(JsonPtr value, const Scene &scene) override;
	// virtual rapidjson::Value toJson(Allocator &allocator) const override;

	virtual bool isHomogeneous() const override;

	virtual void prepareForRender() override;
	// virtual void teardownAfterRender() override;

	virtual Vec3f sigmaA(Vec3f p) const override;
	virtual Vec3f sigmaS(Vec3f p) const override;
	virtual Vec3f sigmaT(Vec3f p) const override;

	virtual bool sampleDistance(PathSampleGenerator &sampler, const Ray &ray,
		MediumState &state, MediumSample &sample) const override;
	virtual bool invertDistance(WritablePathSampleGenerator &sampler, const Ray &ray, bool onSurface) const override;
	virtual Vec3f transmittance(PathSampleGenerator &sampler, const Ray &ray) const override;
	virtual float pdf(PathSampleGenerator &sampler, const Ray &ray, bool onSurface) const override;
	virtual Vec3f transmittanceAndPdfs(PathSampleGenerator &sampler, const Ray &ray, bool startOnSurface,
		bool endOnSurface, float &pdfForward, float &pdfBackward) const override;
	// virtual const PhaseFunction *phaseFunction(const Vec3f &p) const override;
};

}

#endif
