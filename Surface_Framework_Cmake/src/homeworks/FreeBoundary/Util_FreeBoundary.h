#pragma once

#include <queue>
#include <vector>
#include <Eigen/Sparse>

#include "../PolyMesh/include/PolyMesh/PolyMesh.h"

namespace freeb
{
	// return the edge is or not the cutting edge
	// only work for genus 1 and higher mesh, and genus 0 mesh with boundary
	std::vector<char> CutEdge(acamcad::polymesh::PolyMesh* mesh);

	class FreeBoundarySolver
	{
	public:
		void PresetMeshUV(acamcad::polymesh::PolyMesh* mesh);
		void UpdateMeshUV(acamcad::polymesh::PolyMesh* mesh);

	private:
		void ConstrainUV(acamcad::polymesh::PolyMesh* mesh) const;
		bool CheckPlanar(acamcad::polymesh::PolyMesh* mesh,
			const std::vector<acamcad::polymesh::MVert*>& boundaries) const;

		std::vector<double> boundaryEdgesLength(const std::vector<acamcad::polymesh::MVert*>& boundaries) const;

		std::vector<double> ConformWeight(acamcad::polymesh::MVert* v,
			const std::vector<acamcad::polymesh::MVert*>& adjVerts) const;
		std::vector<double> MeanValueWeight(acamcad::polymesh::MVert* v,
			const std::vector<acamcad::polymesh::MVert*>& adjVerts) const;

	private:
		// data member
		Eigen::SparseMatrix<double> m_Coeffs;
		Eigen::VectorXd m_UVList;

		// statistics
		size_t m_Iters{ 0 };
	};
}