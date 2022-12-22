#pragma once

#include <unordered_map>

#include <Eigen\Dense>
#include <Eigen\Sparse>

#include "../PolyMesh/include/PolyMesh/PolyMesh.h"

namespace bff
{
	class BFFSolver
	{
	public:
		void Solve(acamcad::polymesh::PolyMesh* mesh);

	private:
		// mainstream methods
		Eigen::VectorXd CalculateBoundaryLengths(const std::vector<acamcad::polymesh::MVert*>& bVertices) const;
		std::tuple<Eigen::VectorXd, Eigen::VectorXd> CalculateDiscreteCurvatures(acamcad::polymesh::PolyMesh* mesh) const;
		Eigen::SparseMatrix<double> BuildLaplacian(acamcad::polymesh::PolyMesh* mesh);

		
		
		
		// utility methods
		void Util_RearrangeVertIndex(acamcad::polymesh::PolyMesh* mesh);
		double Util_VertAngleDefect(acamcad::polymesh::MVert* vert, acamcad::polymesh::PolyMesh* mesh) const;
		double Util_VertExteriorAngle(acamcad::polymesh::MVert* vert, acamcad::polymesh::PolyMesh* mesh) const;
		double Util_HalfEdgeAngleCotan(acamcad::polymesh::MHalfedge* he, acamcad::polymesh::PolyMesh* mesh) const;

	private:
		// data member
		Eigen::VectorXd m_UVList;

		std::unordered_map<acamcad::polymesh::MVert*, size_t> m_VertIndexMap;

		std::vector<acamcad::polymesh::MVert*> m_BoundaryVertices;

		// statistics
		size_t m_Iters{ 0 };
	};
}