#pragma once

#include <tuple>
#include <vector>
#include <functional>
#include <Eigen\Dense>
#include <Eigen\Sparse>

#include <TinyAD\Scalar.hh>

#include "../PolyMesh/include/PolyMesh/PolyMesh.h"
#include "Util_DataStructure.h"

// 2D, Symmetric Dirichlet Energy Only
namespace eigensys
{
	class ProjectNewtonSolver
	{
	public:
		void PresetMeshUV(acamcad::polymesh::PolyMesh* mesh, const Eigen::SparseMatrix<double>& UVmat);
		bool UpdateMeshUV(acamcad::polymesh::PolyMesh* mesh);

	private:
		Eigen::VectorXd	SetupUVs(acamcad::polymesh::PolyMesh* mesh);
		
		//QPW_EigenSystem2D Eval_Energy_EigenSystem(const QPW_DataPack& pfq_pxq) const;
		
		//Eigen::Matrix4d Calculatep2PSIq_pfq2(const QPW_EigenSystem2D& eigensys) const;
		//Eigen::Vector4d CalculatepPSIq_pfq(const QPW_DataPack& pack);

		std::pair<double, Eigen::VectorXd> Line_Search(
			acamcad::polymesh::PolyMesh* mesh,
			const Eigen::VectorXd& d,
			const Eigen::VectorXd& grad,
			const Eigen::VectorXd& UVs,
			double lastEnergy,
			double gamma, double c);

		double QPW_CalculateEnergySD_2D(const Eigen::Matrix2d& DmINV, const Eigen::Matrix2d& Ds) const;

		double CalculateEnergySD_2D(acamcad::polymesh::PolyMesh* polymesh, const Eigen::VectorXd& UVs) const;

		std::tuple<Eigen::VectorXd, Eigen::SparseMatrix<double>> CalculateGlobalEnergyDerivative(acamcad::polymesh::PolyMesh* mesh, const Eigen::VectorXd& UVs);
		std::tuple<Eigen::Vector<double, 6>, Eigen::Matrix<double, 6, 6>> CalculateLocalEnergyDerivative(acamcad::polymesh::PolyMesh* mesh, const Eigen::Matrix2d& Dm, const Eigen::Matrix2d& Ds);

	private:
		std::vector<Eigen::Matrix2d> m_DmList;
		Eigen::VectorXd m_UVList;

		size_t m_Iters = 0;
		double m_Energy = 0.0;

		bool m_FirstUpdate = true;
	};
}