#pragma once

#include <array>

#include <Eigen\Dense>

#include "../PolyMesh/include/PolyMesh/PolyMesh.h"

namespace eigensys
{
	inline Eigen::Vector4d opVEC(const Eigen::Matrix2d& M)
	{
		return Eigen::Vector4d{ M(0, 0), M(1, 0), M(0, 1), M(1, 1) };
	}

	// return vertex->index(), multiply it by 2 for UV
	// We use idxs as indexes, which refers to the index of vertice in mesh->vertices()
	//		and use ids as indices, which refers to the index of UV, and indices = indexes * 2
	inline Eigen::Vector3i PolyFaceVertIdxs(acamcad::polymesh::PolyMesh* mesh, size_t faceID)
	{
		auto verts = mesh->polygonVertices(mesh->polyface(faceID));
		return Eigen::Vector3i{ verts[0]->index(), verts[1]->index(), verts[2]->index() };
	}

	struct QPW_Invariables
	{
		void CalculateInvariants(const Eigen::Matrix2d& S);
		double I1{ 2.0 };	// I1 = sum(sigma_i)
		double I2{ 2.0 };	// I2 = sum(sigma_i^2)
		double I3{ 1.0 };	// I3 = prod(sigma_i)
	};

	struct QPW_Decomposition
	{
		void CalculateSVDnPolar(const Eigen::Matrix2d& F);
		Eigen::Vector2d sigma;	// { U, Sigma, V } = SVD(F)
		Eigen::Matrix2d U, V;
		Eigen::Matrix2d S, R;	// polar decomposition
	};

	struct QPW_DeformVectors
	{
		QPW_DeformVectors(const Eigen::Matrix2d& F, const QPW_Decomposition& decomp);

		Eigen::Vector4d f{ Eigen::Vector4d::Zero() };		// vec(F)
		Eigen::Vector4d g{ Eigen::Vector4d::Zero() };		// vec(R)
		Eigen::Vector4d r{ Eigen::Vector4d::Zero() };		// r / sqrt(2)
		Eigen::Vector4d t{ Eigen::Vector4d::Zero() };		// vec(U  twist VT) / sqrt(2)
		Eigen::Vector4d p{ Eigen::Vector4d::Zero() };		// vec(U [1 -1] VT) / sqrt(2)
		Eigen::Vector4d l{ Eigen::Vector4d::Zero() };		// vec(U   flip VT) / sqrt(2)
		Eigen::Vector4d d1{ Eigen::Vector4d::Zero() };		// vec(U {1  0} VT)
		Eigen::Vector4d d2{ Eigen::Vector4d::Zero() };		// vec(U {0  1} VT)
	};

	// cell, storing Dm and pFq_pfq
	class QPW_Cell
	{
		friend struct QPW_EigenSystem2D;

	public:
		using vec6d = Eigen::Vector<double, 6>;
		using mat6d = Eigen::Matrix<double, 6, 6>;

		QPW_Cell(const Eigen::Matrix2d& dm);

		double GetVolumeWeight() const { return Dm.determinant() / 2.0; }
		Eigen::Matrix<double, 4, 6> GetDeformGradDerivative() const { return pfq_pxq; }

		std::pair<vec6d, mat6d> CalculateGradNHess(const Eigen::Matrix2d& Ds);
		double CalculateEnergy(const Eigen::Matrix2d& Ds);
		double GetLastEnergy() const { return CalculateEnergy(m_Invariables); }

	private:
		double CalculateEnergy(const QPW_Invariables& invars) const;
		QPW_DeformVectors CalculateVars(const Eigen::Matrix2d& F);
		Eigen::Vector4d CalculateEnergyDeformGrad(const QPW_Invariables& invars, const QPW_DeformVectors& defvecs) const;
		Eigen::Matrix4d CalculateEnergyDeformHess(const QPW_EigenSystem2D& eigensys) const;

	private:
		Eigen::Matrix2d Dm;	// F = Ds / Dm
		Eigen::Matrix<double, 4, 6> pfq_pxq{ Eigen::Matrix<double, 4, 6>::Zero() };	// pf / px per quad-point

		QPW_Decomposition	m_Decomposition;
		QPW_Invariables		m_Invariables;
	};

	// eigen system
	struct QPW_EigenValVecPair
	{
		double l{ 0.0 };
		Eigen::Vector4d e{ Eigen::Vector4d::Zero() };
	};

	struct QPW_EigenSystem2D
	{
		QPW_EigenSystem2D(const QPW_Cell& cell, const QPW_DeformVectors& defvecs);
		std::array<QPW_EigenValVecPair, 4> pairs;
	};
}
