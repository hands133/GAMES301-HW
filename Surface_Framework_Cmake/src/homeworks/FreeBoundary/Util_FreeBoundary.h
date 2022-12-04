#pragma once

#include <queue>
#include <vector>
#include <Eigen/Sparse>

#include "../PolyMesh/include/PolyMesh/PolyMesh.h"

namespace freeb
{
	class FreeBoundarySolver
	{
	public:
		void Solve(acamcad::polymesh::PolyMesh* mesh);

	private:
		void LSCMUV(acamcad::polymesh::PolyMesh* mesh);

		void TutteWithMeanValue(acamcad::polymesh::PolyMesh* mesh);
		void ConstrainUV(acamcad::polymesh::PolyMesh* mesh) const;

		std::vector<double> boundaryEdgesLength(const std::vector<acamcad::polymesh::MVert*>& boundaries) const;
		std::vector<double> MeanValueWeight(acamcad::polymesh::MVert* v,
			const std::vector<acamcad::polymesh::MVert*>& adjVerts) const;

	private:
		// data member
		Eigen::VectorXd m_UVList;

		// statistics
		size_t m_Iters{ 0 };
	};
}