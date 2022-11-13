#pragma once

#include <tuple>
#include <vector>
#include <filesystem>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <Eigen/SparseLU>

#include "../PolyMesh/include/PolyMesh/PolyMesh.h"

#include "Util_DataStructure.h"

// 2D, Symmetric Dirichlet Energy Only
namespace eigensys
{
	class ProjectNewtonSolver
	{
	public:
		ProjectNewtonSolver();
		void PresetMeshUV(acamcad::polymesh::PolyMesh* mesh, const Eigen::SparseMatrix<double>& UVmat);
		bool UpdateMeshUV(acamcad::polymesh::PolyMesh* mesh);
		void SaveEnergies(const std::filesystem::path& filePath) const;

	private:
		void ConstrainUV(acamcad::polymesh::PolyMesh* mesh, const Eigen::VectorXd& UVs);
		
		std::pair<double, Eigen::VectorXd> Line_Search(
			acamcad::polymesh::PolyMesh* mesh,
			const Eigen::VectorXd& d,
			const Eigen::VectorXd& grad,
			double lastEnergy,
			double gamma, double c, double a0);

		// procedure function
		double CalculateEnergySD_2D(acamcad::polymesh::PolyMesh* mesh, const Eigen::VectorXd& UVs) const;
		std::tuple<Eigen::VectorXd, Eigen::SparseMatrix<double>> CalculateEnergyDerivative(acamcad::polymesh::PolyMesh* mesh, const Eigen::VectorXd& UVs);

		double CalculateNoFlipoverStep(acamcad::polymesh::PolyMesh* mesh, const Eigen::VectorXd& grad) const;

	private:
		// data meber
		std::vector<QPW_Cell> m_CellList;
		Eigen::VectorXd m_UVList;

		// statistics
		size_t m_Iters{ 0 };
		double m_Energy{ 0.0 };

		std::vector<double> m_EnergyRecord;

		// utility
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> m_LDLTSolver;

		// control
		bool m_FirstUpdate{ true };
	};
}