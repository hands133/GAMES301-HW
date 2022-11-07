#pragma once

#include <vector>

#include <Eigen\Dense>

#include "../PolyMesh/include/PolyMesh/PolyMesh.h"

namespace eigensys
{
	// axuliary functions
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
		// Initialized with unit matrix S
		double I1{ 2.0 };	// I1 = sum(sigma_i)
		double I2{ 2.0 };	// I2 = sum(sigma_i^2)
		double I3{ 1.0 };	// I3 = prod(sigma_i)
	};

	struct QPW_Decomposition
	{
		void CalculateSVDnPolar(const Eigen::Matrix2d& F);
		Eigen::Matrix2d U, Sigma, V;	// { U, Sigma, V } = SVD(F)
		Eigen::Matrix2d S, R;			// polar decomposition
	};

	struct QPW_DeformGradient
	{
		void CalculateDeformGradient(const Eigen::Matrix2d& Ds, const Eigen::Matrix2d& DmInv);
		Eigen::Matrix2d F{ Eigen::Matrix2d::Identity() };	// deformation gradient F
	};

	struct QPW_DeformVectors
	{
		void CalculateVectors(const Eigen::Matrix2d& F, const QPW_Decomposition& decomp);

		Eigen::Vector4d f{ Eigen::Vector4d::Zero() };		// vec(F)
		Eigen::Vector4d g{ Eigen::Vector4d::Zero() };		// vec(R)
		Eigen::Vector4d r{ Eigen::Vector4d::Zero() };		// r / sqrt(2)
		Eigen::Vector4d t{ Eigen::Vector4d::Zero() };		// vec(U  twist VT) / sqrt(2)
		Eigen::Vector4d p{ Eigen::Vector4d::Zero() };		// vec(U [1 -1] VT) / sqrt(2)
		Eigen::Vector4d l{ Eigen::Vector4d::Zero() };		// vec(U   flip VT) / sqrt(2)
		Eigen::Vector4d d1{ Eigen::Vector4d::Zero() };		// vec(U {1  0} VT)
		Eigen::Vector4d d2{ Eigen::Vector4d::Zero() };		// vec(U {0  1} VT)
	};

	struct QPW_EigenValVecPair
	{
		double l{ 0.0 };
		Eigen::Vector4d e{ Eigen::Vector4d::Zero() };
	};

	struct QPW_EigenSystem2D
	{
		std::array<QPW_EigenValVecPair, 4> valvecpairs;
	};

	// global Energy Gradient
	//struct EnergySD_2DGradient
	//{
	//	Eigen::MatrixXd operator()(
	//		acamcad::polymesh::PolyMesh* mesh,
	//		const Eigen::VectorXd& UVs,
	//		std::vector<QPW_DataPack>& packs);

	//	Eigen::MatrixXd operator()(
	//		acamcad::polymesh::PolyMesh* mesh,
	//		const Eigen::VectorXd& UVs,
	//		const std::vector<Eigen::Matrix2>& DmList);
	//};

	// F Dm = Ds
	//class QPW_DataPack
	//{
	//	friend struct QPW_EnergySD_2D;
	//	friend class ProjectNewtonSolver;
	//public:
	//	QPW_DataPack(const Eigen::Matrix2d& Dm);

	//	void CalculatepF_pDs(const Eigen::Matrix2d& Ds);

	//	// call CalculateGradient ahead!
	//	Eigen::Matrix<double, 4, 6> GetLocalpfq_pxq() const { return pfq_pxq; }
	//	Eigen::MatrixXd GetGlobalpf_px(acamcad::polymesh::PolyMesh* mesh, size_t faceID) const;

	//private:
	//	void Calculatepfq_pxq(const Eigen::Matrix2d& DmInv);

	//private:
	//	Eigen::Matrix2d DmINV;	// F = Ds / Dm
	//	Eigen::Matrix<double, 4, 6> pfq_pxq{ Eigen::Matrix<double, 4, 6>::Zero() };	// pf / px per quad-point

	//	QPW_DeformGradient	m_DeformGradient;
	//	QPW_Decomposition	m_Decomposition;
	//	QPW_Invariables		m_Invariables;
	//	QPW_DeformVectors	m_DeformVectors;
	//};
}
