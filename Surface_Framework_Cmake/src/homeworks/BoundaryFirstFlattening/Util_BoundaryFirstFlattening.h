#pragma once

#include <unordered_map>

#include <Eigen\Dense>
#include <Eigen\Sparse>

#include "../PolyMesh/include/PolyMesh/PolyMesh.h"

namespace tutte
{
	enum class UVBoundaryType;
}

namespace bff
{
	enum class FlattenType { FREE, DISK, FIXED };
	class BFFSolver
	{
	public:
		void Solve(acamcad::polymesh::PolyMesh* mesh, FlattenType type);
		void SetExteriorAngle(size_t numBV, tutte::UVBoundaryType shape);

	private:
		// procedural methods
		Eigen::VectorXd CalculateBoundaryLengths(const std::vector<acamcad::polymesh::MVert*>& bVertices) const;
		std::tuple<Eigen::VectorXd, Eigen::VectorXd> CalculateDiscreteCurvatures() const;
		Eigen::SparseMatrix<double> BuildLaplacian();

		void InitParam();
		void Flatten(const Eigen::VectorXd& givenTarget, bool givenScaleFactor = false);
		void Flatten(const Eigen::VectorXd& scaleFactor,
					 const Eigen::VectorXd& exteriorAngle,
					 bool conjugate = true);
		void FlattenToDisk();
		Eigen::MatrixXd ConstructBestFitCurve(const Eigen::VectorXd& l,
										      const Eigen::VectorXd& lstar,
										      const Eigen::VectorXd& ktilde) const;
		std::tuple<Eigen::VectorXd, Eigen::VectorXd> ExtendCurve(const Eigen::MatrixXd& gamma_tilde, bool freeBoundary);

		// utility methods
		void Util_RearrangeVertIndex();
		size_t Util_VertIndex(acamcad::polymesh::MVert* vert) const { return m_VertIndexMap.find(vert)->second; }
		double Util_VertAngleDefect(acamcad::polymesh::MVert* vert) const;
		double Util_VertExteriorAngle(acamcad::polymesh::MVert* vert) const;
		double Util_HalfEdgeAngleCotan(acamcad::polymesh::MHalfedge* he) const;
		Eigen::VectorXd Util_CvtDirichletToNeumann(const Eigen::VectorXd& phi,
												   const Eigen::VectorXd& g) const;
		Eigen::VectorXd Util_CvtNeumannToDirichlet(const Eigen::VectorXd& phi,
												   const Eigen::VectorXd& h);

		Eigen::VectorXd Util_VecVertConcat(const Eigen::VectorXd& A, const Eigen::VectorXd& B) const;
		double Util_TargetBoundaryLengths(Eigen::VectorXd& lstar) const;
		double Util_TargetDualBoundaryLengths(Eigen::VectorXd& lstar, Eigen::VectorXd& ldual) const;
		void Util_ConstrainUV();

	private:
		// data member
		acamcad::polymesh::PolyMesh* m_Mesh = nullptr;
		std::unordered_map<acamcad::polymesh::MVert*, size_t> m_VertIndexMap;
		bool m_FlattenToDisk = false;

		size_t m_NumV{ 0 };
		size_t m_NumBV{ 0 };
		size_t m_NumIV{ 0 };

		std::vector<acamcad::polymesh::MVert*> m_BoundaryVertices;
		Eigen::VectorXd m_BoundaryEdgeLengths;
		Eigen::VectorXd m_TargetEdgeLengths;
		Eigen::VectorXd m_BoundaryScaleFactor;
		Eigen::VectorXd m_TargetExteriorAngle;
		
		Eigen::VectorXd m_GaussianCurv, m_GeodesicCurv;
		Eigen::SparseMatrix<double> m_LapMat, m_LapII, m_LapIB, m_LapBB;
		Eigen::MatrixXd m_UVList;

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> m_LUSolver;
	};
}