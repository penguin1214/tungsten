#include "PlainGrid.hpp"

#include <iostream>
#include <fstream>

#include "Common.h"

#include "io/JsonObject.hpp"
#include "io/Scene.hpp"

namespace Tungsten {

	/// Constructor
	PlainGrid::PlainGrid()
		:_stepSize(1.0f) {
	}

	/// Function

	void PlainGrid::fromJson(JsonPtr value, const Scene &scene) {
		if (auto path = value["data_file"]) _path = scene.fetchResource(path);
		value.getField("grid_name", _gridName);
		value.getField("step_size", _stepSize);
		value.getField("transform", _configTransform);	// parse natural transform
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

#if 0
		for (int z = 0; z < size_z; z++) {
			for (int y = 0; y < size_y; y++) {
				for (int x = 0; x < size_x; x++) {
					if (z == 20) {}
					else {
						int idx = size_x*size_y*z + size_x*y + x;
						_data[idx] = 0.0;
					}
				}
			}
		}
#endif

#if GRID_NORMAL_TRANSFORM
		/// Transform grid to local space for intersection
		// do natural transforms first
		// for test, the single scale grid starts from (0,0,0)
		Vec3i minP(0, 0, 0); Vec3i maxP(size_x - 1, size_y - 1, size_z - 1);
		Vec3f diag = Vec3f(maxP - minP);	// diagonal
		float scale = 1.0f / diag.max();	// chooe greatest one as reference
		diag *= scale;
		// center is center in bottom surface
		// Vec3f center = Vec3f(minP)*scale + Vec3f(diag.x(), 0.0f, diag.z())*0.5f;
		// use grid center
		Vec3f center = Vec3f(minP)*scale + Vec3f(diag.x(), diag.y(), diag.z())*0.5f;
		std::cout << minP << " -> " << maxP << std::endl;

		/// **ATTENTION** NOT CONSISTENT WITH PLAINGRID NOW!
		// transform to [-0.5, 0.5] now
		float o = 1.0f;
		float z = 0.0f;
		Mat4f m(o, z, z, z,
			z, o, z, z,
			z, z, -o, z,
			z, z, z, o);

		_scaleTransform = m * Mat4f::translate(-center)*Mat4f::scale(Vec3f(scale));
		_invScaleTransform = Mat4f::scale(Vec3f(1.0f / scale))*Mat4f::translate(center) * m;
		
		_invConfigTransform = _configTransform.invert();

		_bounds = Box3f(Vec3f(minP), Vec3f(maxP));
		std::cout << _bounds << std::endl;
#else
		// do not do normalization
		Vec3f minP(0.0f, 0.0f, 0.0f); Vec3f maxP(float(size_x - 1), float(size_y - 1), float(size_z - 1));
		_invConfigTransform = _configTransform.invert();
		_bounds = Box3f(_configTransform * Vec3f(minP), _configTransform * Vec3f(maxP));
		std::cout << "transformed bound: " << _bounds << std::endl;
#endif		

		/// precompute max density
		_maxDensity = precomputeMaxDensity();
		precomputeGradient();
	}

	/// Getter
	/* localScale to worldScale transform */
	Mat4f PlainGrid::scaleTransform() const {
		return _scaleTransform;
	}

	/* worldScale to localScale transform */
	Mat4f PlainGrid::invScaleTransform() const {
		return _invScaleTransform;
	}

	Mat4f PlainGrid::configTransform() const {
		return _configTransform;
	}

	Mat4f PlainGrid::invConfigTransform() const {
		return _invConfigTransform;
	}

	Mat4f PlainGrid::naturalTransform() const {
#if GRID_NORMAL_TRANSFORM
		return _configTransform*_scaleTransform;
#else
		return _configTransform;
#endif
	}

	Mat4f PlainGrid::invNaturalTransform() const {
#if GRID_NORMAL_TRANSFORM
		return _invScaleTransform*_invConfigTransform;
#else
		return _invConfigTransform;
#endif
	}

	Box3f PlainGrid::bounds() const {
		return _bounds;
	}

	/* return 1D index of position p */
	int PlainGrid::indexAt(Vec3f p) const {
		return int(p.z())*size_x*size_y + int(p.y())*size_x + int(p.x());
	}

	int PlainGrid::indexAt(int x, int y, int z) const {
		return z*size_x*size_y + y*size_x + x;
	}

	/*double interpolate(double t, double f1, double f2) {
		return (double)(f1 + (f2 - f1)*t);
	}

	double linearInterpolate(Vec3f p, Vec3f o, int i, int j, int k, double x, double y, double z) {
		int i = int(p.x()); int j = int(p.y()); int k = int(p.z());
		double t1 = interpolate(x, g[indexAt(i, j, k)], g[indexAt(i + 1, j, k)]);
		double t2 = interpolate(x, g[indexAt(i, j + 1, k)], g[indexAt(i + 1, j + 1, k)]);
		double t3 = interpolate(x, g[indexAt(i, j, k + 1)], g[indexAt(i + 1, j, k + 1)]);
		double t4 = interpolate(x, g[indexAt(i, j + 1, k + 1)], g[indexAt(i + 1, j + 1, k + 1)]);

		double t5 = interpolate(y, t1, t2);
		double t6 = interpolate(y, t3, t4);

		double t7 = interpolate(z, t5, t6);
		return t7;
	}*/

	float PlainGrid::density(Vec3f p) const {
		if (p.x() < 0.0 || p.x() > size_x - 1 || p.y() < 0.0 || p.y() > size_y - 1 || p.z() < 0.0 || p.z() > size_z - 1) return 0.0;
		return _data[indexAt(p)];
	}

	Vec3f PlainGrid::gradient(Vec3f p) const {
		if (p.x() < 0.0 || p.x() > size_x - 1 || p.y() < 0.0 || p.y() > size_y - 1 || p.z() < 0.0 || p.z() > size_z - 1) return Vec3f(0.0);
		return _gradient[indexAt(p)];
	}

	Vec3f PlainGrid::transmittance(PathSampleGenerator &sampler, Vec3f p, Vec3f w, float t0, float t1, Vec3f sigmaT) const {
#if 0
		float ta = t0;
		float fa = density(p + w*t0);// interpolation
		float integral = 0.0f;
		float dT = sampler.next1D()*_stepSize;
		do {
			float tb = min(ta + dT, t1);
			float fb = density(p + w*tb);
			integral += (fa + fb)*0.5f*(tb - ta);
			ta = tb;
			fa = fb;
			dT = _stepSize;
		} while (ta < t1);
		return std::exp(-integral*sigmaT);
#endif
		/// ray marching
		Vec3f T(0.0);
		float t;
		float d;

		for (t = t0; t < t1; t += _stepSize) {
#if USE_DENSITY
			d = density(p + t*w);
#else
			d = gradient(p + t*w).max();
#endif
			T += d *_stepSize;
		}

		T -= (t - t1) * d;

		return std::exp(-T * sigmaT);
	}

	/* delta tracking iteration */
	Vec2f PlainGrid::inverseOpticalDepth(PathSampleGenerator &sampler, Vec3f p, Vec3f w, float t0, float t1,
		float sigmaT, float xi) const {
		float ta = t0;
		float fa = density(p + w*t0);
		float integral = 0.0f;
		float dT = sampler.next1D()*_stepSize;
		do {
			float tb = min(ta + dT, t1);
			float fb = density(p + w*tb);
			float delta = (fa + fb)*sigmaT*0.5f*(tb - ta);	// k(x)*dx	optical depth

			if (integral + delta >= xi) {
				float a = (fb - fa)*sigmaT;
				float b = fa*sigmaT;
				float c = (integral - xi) / (tb - ta);
				float mantissa = max(b*b - 2.0f*a*c, 0.0f);
				float x1 = (-b + std::sqrt(mantissa)) / a;
				return Vec2f(ta + (tb - ta)*x1, fa + (fb - fa)*x1);
			}

			integral += delta;
			ta = tb;
			fa = fb;
			dT = _stepSize;
		} while (ta < t1);
		return Vec2f(t1, integral);
	}

	float PlainGrid::maxDensity() {
		return _maxDensity;
	}

	Vec3f PlainGrid::maxGradient() {
		return _maxGradient;
	}

	float PlainGrid::precomputeMaxDensity() {
		float m = 0.0f;
		for (int i = 0; i < size_x*size_y*size_z; i++) {
			if (_data[i] > m) m = _data[i];
		}
		return m;
	}

	float PlainGrid::stepSize() const {
		return _stepSize;
	}

	float curve(float x, float k) {
		// "half sigmoid"
		// f(x) = kx / k-x+1
		float ret;
		if (x < 0.5) ret = x;
		else ret = 1.0f / (1.0f + std::exp(-k*(x - 0.5)));
		/*std::ofstream f;
		f.open("curve", std::ios::out | std::ios::app);
		if (f.good()) {
			f << x << "," << ret << std::endl;
		}
		f.close();*/
		//return k*x / (k - x + 1.0);

		/// https://en.wikipedia.org/wiki/Logistic_function
		// 1 / 1+e^(-k*(x-0.5))
		return ret;
	}
	void PlainGrid::precomputeGradient() {
		_gradient = new Vec3f[size_x*size_y*size_z];

		float inv_d = 1.0 / _stepSize;
		float fx, fy, fz;
		float max_x = 0.0, max_y = 0.0, max_z = 0.0;
		float min_x = 0.0, min_y = 0.0, min_z = 0.0;
#if 0
		for (int x = 1; x < size_x - 1; x++) {
			for (int y = 1; y < size_y; y++) {
				for (int z = 1; z < size_z - 1; z++) {
					int index = size_x * size_y * z + size_x * y + x;
					if ((x - 1) < 0 || (y - 1) < 0 || (z - 1) < 0) {
						// at edge, use forward difference
						fx = curve(_data[indexAt(x + 1, y, z)] - _data[indexAt(x, y, z)] * inv_d, 6.5);
						fy = curve(_data[indexAt(x, y + 1, z)] - _data[indexAt(x, y, z)] * inv_d, 6.5);
						fz = curve(_data[indexAt(x, y, z + 1)] - _data[indexAt(x, y, z)] * inv_d, 6.5);
					}
					else if (x == size_x || y == size_y || z == size_z) {
						// at edge, use backward difference
						fx = curve(_data[indexAt(x, y, z)] - _data[indexAt(x - 1, y, z)] * inv_d, 6.5);
						fy = curve(_data[indexAt(x, y, z)] - _data[indexAt(x, y - 1, z)] * inv_d, 6.5);
						fz = curve(_data[indexAt(x, y, z)] - _data[indexAt(x, y, z - 1)] * inv_d, 6.5);
					}
					else {
						// inside, use central difference
						// 1 / 2h *[f(x + h) - f(x - h)]
						fx = curve(_data[indexAt(x + 1, y, z)] - _data[indexAt(x - 1, y, z)] * inv_d * 0.5, 6.5);
						fy = curve(_data[indexAt(x, y + 1, z)] - _data[indexAt(x, y - 1, z)] * inv_d * 0.5, 6.5);
						fz = curve(_data[indexAt(x, y, z + 1)] - _data[indexAt(x, y, z - 1)] * inv_d * 0.5, 6.5);
					}
					_gradient[index] = Vec3f(fx, fy, fz);
					if (fx > max_x) max_x = fx; if (fx < min_x) min_x = fx;
					if (fy > max_y) max_y = fy; if (fy < min_y) min_y = fy;
					if (fz > max_z) max_z = fz; if (fz < min_z) min_z = fz;
				}
			}
		}

#endif

		for (int x = 1; x < size_x - 1; x++) {
			for (int y = 1; y < size_y; y++) {
				for (int z = 1; z < size_z - 1; z++) {
					int index = size_x * size_y * z + size_x * y + x;
					if ((x - 1) < 0 || (y - 1) < 0 || (z - 1) < 0) {
						// at edge, use forward difference
						fx = _data[indexAt(x + 1, y, z)] - _data[indexAt(x, y, z)] * inv_d;
						fy = _data[indexAt(x, y + 1, z)] - _data[indexAt(x, y, z)] * inv_d;
						fz = _data[indexAt(x, y, z + 1)] - _data[indexAt(x, y, z)] * inv_d;
					}
					else if (x == size_x || y == size_y || z == size_z) {
						// at edge, use backward difference
						fx = _data[indexAt(x, y, z)] - _data[indexAt(x - 1, y, z)] * inv_d;
						fy = _data[indexAt(x, y, z)] - _data[indexAt(x, y - 1, z)] * inv_d;
						fz = _data[indexAt(x, y, z)] - _data[indexAt(x, y, z - 1)] * inv_d;
					}
					else {
						// inside, use central difference
						// 1 / 2h *[f(x + h) - f(x - h)]
						fx = _data[indexAt(x + 1, y, z)] - _data[indexAt(x - 1, y, z)] * inv_d * 0.5;
						fy = _data[indexAt(x, y + 1, z)] - _data[indexAt(x, y - 1, z)] * inv_d * 0.5;
						fz = _data[indexAt(x, y, z + 1)] - _data[indexAt(x, y, z - 1)] * inv_d * 0.5;
					}
					_gradient[index] = Vec3f(fx, fy, fz);
					if (fx > max_x) max_x = fx; if (fx < min_x) min_x = fx;
					if (fy > max_y) max_y = fy; if (fy < min_y) min_y = fy;
					if (fz > max_z) max_z = fz; if (fz < min_z) min_z = fz;
				}
			}
		}


		float denom_x = 1.0f / (max_x - min_x); float denom_y = 1.0f / (max_y - min_y); float denom_z = 1.0f / (max_z - min_z);
		float mx = 0.0; float my = 0.0; float mz = 0.0;
		for (int x = 1; x < size_x - 1; x++) {
			for (int y = 1; y < size_y; y++) {
				for (int z = 1; z < size_z - 1; z++) {
					int index = size_x * size_y * z + size_x * y + x;
					fx = _gradient[index].x(); fy = _gradient[index].y(); fz = _gradient[index].z();
					fx *= denom_x; fy *= denom_y; fz *= denom_z;
					if (fx > mx) mx = fx; if (fy > my) my = fy; if (fz > mz) mz = fz;
					_gradient[index] = Vec3f(fx, fy, fz);
				}
			}
		}


		// compute max gradient
		Vec3f m(0.0);
		_maxGradient = Vec3f(max_x, max_y, max_z);
		for (int i = 0; i < size_x*size_y*size_z; i++) {
			if (_gradient[i].max() > m.max()) m = _gradient[i];
		}
		//_maxGradient = Vec3f(max_x, max_y, max_z);
		_maxGradient = Vec3f(mx, my, mz);
	}
}
