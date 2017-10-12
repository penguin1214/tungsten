#include "PlainGrid.hpp"

#include <iostream>

#include "io/JsonObject.hpp"
#include "io/Scene.hpp"

namespace Tungsten {

/// Constructor
PlainGrid::PlainGrid()
:_stepSize(1.0f){}

/// Function

void PlainGrid::fromJson(JsonPtr value, const Scene &scene) {
	if (auto path = value["data_file"]) _path = scene.fetchResource(path);	// get file name?
	value.getField("grid_name", _gridName);
	value.getField("step_size", _stepSize);
	// value.getField("transform"), _transform;
}

template<class T>
T* loadBinary(const char * fn, int &gridx, int &gridy, int &gridz) {
	if (sizeof(T) != sizeof(float))return nullptr;
	FILE *fp = fopen(fn, "rb");
	fread(&gridx, 1, sizeof(int), fp);
	fread(&gridy, 1, sizeof(int), fp);
	fread(&gridz, 1, sizeof(int), fp);

	int total = gridx * gridy * gridz;
	T *data = new T[total];

	fread(data, sizeof(T), total, fp);
	/*for (int i = 0; i < total; i++) {
	if (data[i] > 1) {
	std::cout << data[i] << ", ";
	}
	}*/
	fclose(fp);
	printf("loaded %s <%d,%d,%d>\n", fn, gridx, gridy, gridz);
	return data;
}

/* load binary file */
void PlainGrid::loadResources() {
	const std::string str_filename = _path->absolute().asString();
	const char* filename = str_filename.c_str();
	_data = loadBinary<float>(filename, size_x, size_y, size_z);

	/// set transform
	/// TODO
	// for test, the single scale grid starts from (0,0,0)
	Vec3i minP(0, 0, 0); Vec3i maxP(size_x-1, size_y-1, size_z-1);
	Vec3f diag = Vec3f(maxP - minP);	// diagonal
	float scale = 1.0f / diag.max();	// chose greatest one as reference
	diag *= scale;
	// center is center in bottom surface
	// Vec3f center = Vec3f(minP)*scale + Vec3f(diag.x(), 0.0f, diag.z())*0.5f;
	// use grid center
	Vec3f center = Vec3f(minP)*scale + Vec3f(diag.x(), diag.y(), diag.z())*0.5f;
	std::cout << minP << " -> " << maxP << std::endl;

	/// **ATTENTION** NOT CONSISTENT WITH PLAINGRID NOW!
	// transform to [-0.5, 0.5] now
	_transform = Mat4f::translate(-center)*Mat4f::scale(Vec3f(scale));
	_invTransform = Mat4f::scale(Vec3f(1.0f / scale))*Mat4f::translate(center);
	_bounds = Box3f(Vec3f(minP), Vec3f(maxP));
}

/// Getter
Mat4f PlainGrid::worldToLocalTransform() const {
	return _transform;
}

Mat4f PlainGrid::localToWorldTransform() const {
	return _invTransform;
}

Box3f PlainGrid::bounds() const {
	return _bounds;
}

/// TODO
/* return 1D index of position p */
int PlainGrid::indexAt(Vec3f p) const { return p.z()*size_x*size_y + p.y()*size_x + p.x(); }

float PlainGrid::density(Vec3f p) const {
	// p -> index
	// return _data[index]
	/// TODO: check bound efficiently
	if (p.x() < 0.0 || p.x() > size_x || p.y() < 0.0 || p.y() > size_y || p.z() < 0.0 || p.z() > size_z) return 0.0;
	// if (_data[indexAt(p)] > 0.0) std::cout << _data[indexAt(p)] << std::endl;
	return _data[indexAt(p)];
}

Vec3f PlainGrid::transmittance(PathSampleGenerator &sampler, Vec3f p, Vec3f w, float t0, float t1, Vec3f sigmaT) const {
	// ray marching
	// T(s_i, s_i+1) = exp(-(k(t)*dt), where k(t) is absorption coefficient, assuming absorption only model.
	// k(t)*dt is approximation of integral of k(t) from s_i to s_i+1, i.e. optical depth

	float ta = t0;
	float fa = density(p + w*t0);// interpolation
	float integral = 0.0f;
	float dT = sampler.next1D()*_stepSize;
	do {
		float tb = min(ta + dT, t1);
		float fb = density(p + w*tb);
		integral += (fa + fb)*0.5f*(tb - ta);	// ?
		ta = tb;
		fa = fb;
		dT = _stepSize;
	} while (ta < t1);
	return std::exp(-integral*sigmaT);
}

Vec2f PlainGrid::inverseOpticalDepth(PathSampleGenerator &sampler, Vec3f p, Vec3f w, float t0, float t1,
	float sigmaT, float xi) const {
	return Vec2f(0.0);
}




}