#include "Grid.hpp"

namespace Tungsten {

	Mat4f Grid::scaleTransform() const {
		return Mat4f();
	}

	Mat4f Grid::invScaleTransform() const {
		return Mat4f();
	}

	Mat4f Grid::configTransform() const {
		return Mat4f();
	}

	Mat4f Grid::invConfigTransform() const {
		return Mat4f();
	}

	Mat4f Grid::naturalTransform() const {
		return Mat4f();
	}

	Mat4f Grid::invNaturalTransform() const {
		return Mat4f();
	}

	//Mat4f Grid::worldToLocalNaturalTransform() const {
	//	return Mat4f();
	//}
	//
	//Mat4f Grid::localToWorldNaturalTransform() const {
	//	return Mat4f();
	//}

	Box3f Grid::bounds() const {
		return Box3f();
	}

}
