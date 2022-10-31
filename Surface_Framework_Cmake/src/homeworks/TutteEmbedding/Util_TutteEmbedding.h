#include <vector>
#include <glm\glm.hpp>

#include "../PolyMesh/include/PolyMesh/PolyMesh.h"

namespace tutte
{
	enum class UVBoundaryType {
		POLYGON_CIRCLE,
		POLYGON_TRIANGLE,
		POLYGON_SQUARE,
		POLYGON_PENTAGON,
		// TODO : uncomment this to implement more complex boundaries, for later experiments
		// POLYGON_STAR,
		// POLYGON_CROSS
	};

	std::vector<glm::dvec2> GetBoundaryUVs(size_t numVerts, UVBoundaryType type);

	bool IsOinTriangle(glm::dvec2 P1, glm::dvec2 P2, glm::dvec2 P3);

	auto FloaterParam_I_a(acamcad::polymesh::MVert* v,
		const std::vector<acamcad::polymesh::MVert*>& adjVerts);
	auto FloaterParam_I_b(double thetaSum, const std::vector<glm::dvec2>& neighborAngleInfo);
	auto FloaterParam_II(const std::vector<glm::dvec2>& adjUVs);

	void AverageParam(acamcad::polymesh::MVert* v,
		const std::vector<acamcad::polymesh::MVert*>& adjVerts,
		std::vector<double>& weights);

	void FloaterParam(acamcad::polymesh::MVert* v,
		const std::vector<acamcad::polymesh::MVert*>& adjVerts,
		std::vector<double>& weights);
}