#include <vector>

#include <Eigen3\Dense>

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

	std::vector<Eigen::Vector2d> GetBoundaryUVs(size_t numVerts, UVBoundaryType type);

	bool IsOinTriangle(Eigen::Vector2d P1, Eigen::Vector2d P2, Eigen::Vector2d P3);

	auto FloaterParam_I_a(acamcad::polymesh::MVert* v,
		const std::vector<acamcad::polymesh::MVert*>& adjVerts);
	auto FloaterParam_I_b(double thetaSum, const std::vector<Eigen::Vector2d>& neighborAngleInfo);
	auto FloaterParam_II(const std::vector<Eigen::Vector2d>& adjUVs);

	void AverageParam(acamcad::polymesh::MVert* v,
		const std::vector<acamcad::polymesh::MVert*>& adjVerts,
		std::vector<double>& weights);

	void FloaterParam(acamcad::polymesh::MVert* v,
		const std::vector<acamcad::polymesh::MVert*>& adjVerts,
		std::vector<double>& weights);
}