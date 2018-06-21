#ifndef MEDIUM_HPP_
#define MEDIUM_HPP_

#include "phasefunctions/PhaseFunction.hpp"

#include "samplerecords/MediumSample.hpp"

#include "sampling/WritablePathSampleGenerator.hpp"

#include "math/Ray.hpp"

#include "io/JsonSerializable.hpp"

#include "grids/Grid.hpp"
#include "grids/PlainGrid.hpp"

#include <memory>

namespace Tungsten {

class Scene;

/*
 * Distance sampling according to a single color channel of the extinction coefficient
 * can be seen as just one sampling strategy of the integrand.
 * There are several sampling strategies available,
 * however - one for each color channel.
 * This means that we can randomly pick one of the color channels to sample a distance from,
 * and then weight the result with multiple importance sampling using the pdfs of all the other color channels.
 * The results achieved using this method are magnitudes better in noise for strongly spectrally varying media than just sampling with one strategy.
 */
class Medium : public JsonSerializable
{
protected:
    std::shared_ptr<PhaseFunction> _phaseFunction;
    int _maxBounce;

public:
    struct MediumState
    {
        bool firstScatter;
        int component;
        int bounce;

        void reset()
        {
            firstScatter = true;
            bounce = 0;
        }

        void advance()
        {
            firstScatter = false;
            bounce++;
        }
    };

    Medium();

    virtual void fromJson(JsonPtr value, const Scene &scene) override;
    virtual rapidjson::Value toJson(Allocator &allocator) const override;

    virtual bool isHomogeneous() const = 0;

    virtual void prepareForRender() {}
    virtual void teardownAfterRender() {}

    virtual Vec3f sigmaA(Vec3f p) const = 0;
    virtual Vec3f sigmaS(Vec3f p) const = 0;
    virtual Vec3f sigmaT(Vec3f p) const = 0;	// extinction
	virtual std::shared_ptr<PlainGrid> grid();

	/* pbrt-v3 section 15.2 */
    virtual bool sampleDistance(PathSampleGenerator &sampler, const Ray &ray,
            MediumState &state, MediumSample &sample) const = 0;
    virtual bool invertDistance(WritablePathSampleGenerator &sampler, const Ray &ray, bool onSurface) const;
    virtual Vec3f transmittance(PathSampleGenerator &sampler, const Ray &ray) const = 0;
    virtual float pdf(PathSampleGenerator &sampler, const Ray &ray, bool onSurface) const = 0;
    virtual Vec3f transmittanceAndPdfs(PathSampleGenerator &sampler, const Ray &ray, bool startOnSurface,
            bool endOnSurface, float &pdfForward, float &pdfBackward) const;
    virtual const PhaseFunction *phaseFunction(const Vec3f &p) const;
};

}



#endif /* MEDIUM_HPP_ */
