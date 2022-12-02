#include "Util_EnergyEigensystems.h"

#include <numeric>
#include <algorithm>
#include <fstream>

namespace eigensys
{
	// constrain UV in range [0, 1]x[0, 1]
	void ProjectNewtonSolver::ConstrainUV(acamcad::polymesh::PolyMesh* mesh, const Eigen::VectorXd& UVs)
	{
		double xmin = std::numeric_limits<double>::max();
		double xmax = -std::numeric_limits<double>::max();
		double ymin = std::numeric_limits<double>::max();
		double ymax = -std::numeric_limits<double>::max();

		size_t numV = mesh->vertices().size();
		for (size_t vertID = 0; vertID < numV; ++vertID)
		{
			double x = UVs(vertID * 2 + 0);
			double y = UVs(vertID * 2 + 1);
			
			xmin = std::min(xmin, x);
			xmax = std::max(xmax, x);
			ymin = std::min(ymin, y);
			ymax = std::max(ymax, y);
		}

		double midx = (xmin + xmax) * 0.5;
		double midy = (ymin + ymax) * 0.5;

		for (size_t vertID = 0; vertID < numV; ++vertID)
		{
			auto* vert = mesh->vert(vertID);
			float UVx = static_cast<float>((UVs(2 * vertID + 0) - midx) / (xmax - xmin));
			float UVy = static_cast<float>((UVs(2 * vertID + 1) - midy) / (ymax - ymin));
			vert->setTexture(UVx * 2.0, UVy * 2.0, 0.0);
		}
	}

	ProjectNewtonSolver::ProjectNewtonSolver()
	{
		m_UVList.setZero();

		m_EnergyRecord.reserve(1000);

		m_LDLTSolver.setShift(std::numeric_limits<float>::epsilon());	// necessary
	}

	// x comes from tutte's embedding results
	void ProjectNewtonSolver::PresetMeshUV(acamcad::polymesh::PolyMesh* mesh, const Eigen::SparseMatrix<double>& UVmat)
	{
		assert(mesh != nullptr);
		m_Iters = 0;
		m_FirstUpdate = true;

		size_t numF = mesh->polyfaces().size();
		size_t numV = mesh->vertices().size();

		m_UVList = Eigen::VectorXd(2 * numV);
		m_UVList.setZero();
		for (size_t vertID = 0; vertID < numV; ++vertID)
		{
			m_UVList(2 * vertID + 0) = UVmat.coeff(vertID, 0);
			m_UVList(2 * vertID + 1) = UVmat.coeff(vertID, 1);
		}

		m_CellList.clear();
		m_CellList.reserve(numF);
		for (size_t faceID = 0; faceID < numF; ++faceID)
		{
			auto faceVerts = mesh->polygonVertices(mesh->polyface(faceID));
			auto* v0 = faceVerts[0];
			auto* v1 = faceVerts[1];
			auto* v2 = faceVerts[2];

			auto v10 = v1->position() - v0->position();
			auto v20 = v2->position() - v0->position();

			// y
			// ¦«
			// |    v2
			// .<---.
			// ¦«   /|\
			// |  / | \
			// | /  |  \
			// |/   V   \
			// .===>.--->.-> x
			// v0   vp   v1
			Eigen::Vector2d P0{ 0, 0 };
			Eigen::Vector2d P1{ v10.norm(), 0 };
			Eigen::Vector2d P2{ v10.dot(v20) / v10.norm(), v10.cross(v20).norm() / v10.norm() };

			Eigen::Matrix2d Dm;
			Dm << P1 - P0, P2 - P0;
			m_CellList.emplace_back(Dm);
		}

		m_EnergyRecord.clear();
	}

	// Please make sure that the UV point are stored in the texture coordinate
	bool ProjectNewtonSolver::UpdateMeshUV(acamcad::polymesh::PolyMesh* mesh)
	{
		if (mesh == nullptr)
		{
			std::cout << "the mesh pointer is NULL!\n";
			return false;
		}

		size_t numV = mesh->vertices().size();
		size_t numF = mesh->polyfaces().size();

		std::cout << "[ProjectNewton] ";

		auto [G, H] = CalculateEnergyDerivative(mesh, m_UVList);

		std::cout << "max |Gi| = " << G.cwiseAbs().maxCoeff() << "\t";

		if (G.cwiseAbs().maxCoeff() < std::numeric_limits<float>::epsilon())
		{
			std::cout << "================CONVERGE!!!================\n";
			m_Iters = 0;
			return false;
		}

		if (m_FirstUpdate)
		{
			m_Energy = CalculateEnergySD_2D(mesh, m_UVList);
			m_EnergyRecord.emplace_back(m_Energy);
			m_FirstUpdate = false;
		}

		m_LDLTSolver.compute(H);
		Eigen::VectorXd d = m_LDLTSolver.solve(-G);

		double step = CalculateNoFlipoverStep(mesh, d);
		std::cout << "A = " << step << "\t";

		auto [energy, updatedUV] = Line_Search(mesh, d, G, m_Energy, 0.8, 1.0e-4, std::min(step * 0.8, 1.0));
		m_EnergyRecord.emplace_back(energy);

		m_Energy = energy;
		m_UVList = updatedUV;

		ConstrainUV(mesh, updatedUV);

		return true;
	}

	void ProjectNewtonSolver::SaveEnergies(const std::filesystem::path& filePath) const
	{
		std::ofstream file(filePath);
		if (!file.is_open()) {
			file.close();
			assert(false, "file open error!");
		}

		for (double energy : m_EnergyRecord)
			file << energy << ", ";

		file.close();
	}

	std::pair<double, Eigen::VectorXd> ProjectNewtonSolver::Line_Search(
		acamcad::polymesh::PolyMesh* mesh,
		const Eigen::VectorXd& d,
		const Eigen::VectorXd& grad,
		double lastEnergy,
		double gamma, double c, double a0)
	{
		double alpha = a0;
		double currentEnergy = lastEnergy;

		const int MAX_LOOP_SIZE = 64;

		auto updatedUVs = m_UVList;
		size_t iter = 0;
		for (iter = 0; iter < MAX_LOOP_SIZE; ++iter)
		{
			updatedUVs = m_UVList + alpha * d;
			currentEnergy = CalculateEnergySD_2D(mesh, updatedUVs);

			if (currentEnergy <= lastEnergy + c * alpha * d.dot(grad))	break;
			alpha *= gamma;
		}
		std::cout << "Iter = " << ++m_Iters << "\tEo = " << lastEnergy
			<< "\tEn = " << currentEnergy << "\talpha = " << alpha << "\n";

		return std::make_pair(currentEnergy, updatedUVs);
	}

	double ProjectNewtonSolver::CalculateEnergySD_2D(
		acamcad::polymesh::PolyMesh* mesh,
		const Eigen::VectorXd& UVs) const
	{
		double energy = 0.0;
		size_t numF = mesh->polyfaces().size();
		
		for (size_t faceID = 0; faceID < numF; ++faceID)
		{
			auto ids = PolyFaceVertIdxs(mesh, faceID) * 2;

			Eigen::Matrix2d Ds = Eigen::Matrix2d::Zero();
			Ds << UVs.segment(ids(1), 2) - UVs.segment(ids(0), 2),
				  UVs.segment(ids(2), 2) - UVs.segment(ids(0), 2);

			auto& cell = m_CellList[faceID];
			energy += cell.GetVolumeWeight() * cell.GetLastEnergy();
		}
		
		return energy;
	}

	std::tuple<Eigen::VectorXd, Eigen::SparseMatrix<double>>
		ProjectNewtonSolver::CalculateEnergyDerivative(
		acamcad::polymesh::PolyMesh* mesh, const Eigen::VectorXd& UVs)
	{
		size_t numV = mesh->vertices().size();
		size_t numF = mesh->polyfaces().size();

		Eigen::VectorXd GRAD(2 * numV);
		GRAD.setZero();

		Eigen::SparseMatrix<double> HESS(2 * numV, 2 * numV);
		HESS.setZero();
		std::vector<Eigen::Triplet<double>> HESStrips;
		HESStrips.reserve(6 * numF);

		for (size_t faceID = 0; faceID < numF; ++faceID)
		{
			auto ids = PolyFaceVertIdxs(mesh, faceID) * 2;
			int GlobalCoord[6] = { ids(0), ids(0) + 1, ids(1), ids(1) + 1, ids(2), ids(2) + 1 };

			Eigen::Matrix2d Ds;
			Ds << UVs.segment(ids(1), 2) - UVs.segment(ids(0), 2),
				  UVs.segment(ids(2), 2) - UVs.segment(ids(0), 2);

			auto [gradq, hessq] = m_CellList[faceID].CalculateGradNHess(Ds);
			
			for (int n = 0; n < 6; ++n)
			{
				GRAD(GlobalCoord[n]) += gradq(n);
				for (int m = 0; m < 6; ++m)
					HESStrips.emplace_back(GlobalCoord[n], GlobalCoord[m], hessq(n, m));
			}
			assert(!std::isnan(GRAD.norm()));
		}
		
		HESS.setFromTriplets(HESStrips.begin(), HESStrips.end());
		HESS.makeCompressed();
		return std::make_tuple(GRAD, HESS);
	}

	double ProjectNewtonSolver::CalculateNoFlipoverStep(acamcad::polymesh::PolyMesh* mesh, const Eigen::VectorXd& grad) const
	{
		size_t numF = mesh->polyfaces().size();

		const auto& UVs = m_UVList;
		double minStep = 1.0;

		for (size_t faceID = 0; faceID < numF; ++faceID)
		{
			auto ids = PolyFaceVertIdxs(mesh, faceID) * 2;
			Eigen::Matrix2d Ds = Eigen::Matrix2d::Zero();
			Ds << UVs.segment(ids(1), 2) - UVs.segment(ids(0), 2),
				  UVs.segment(ids(2), 2) - UVs.segment(ids(0), 2);

			Eigen::Matrix2d D = Eigen::Matrix2d::Zero();
			D << grad.segment(ids(1), 2) - grad.segment(ids(0), 2),
				 grad.segment(ids(2), 2) - grad.segment(ids(0), 2);

			double A = D.determinant();
			double B = D.determinant() + Ds.determinant() - (Ds - D).determinant();
			double C = Ds.determinant();

			double t1 = 0.0;
			double t2 = 0.0;
			if (std::abs(A) > 1.0e-10)
			{
				double Delta = B * B - 4 * A * C;
				if (Delta <= 0)	continue;

				double delta = std::sqrt(Delta); // delta >= 0
				if (B >= 0)
				{
					double bd = -B - delta;
					t1 = 2 * C / bd;
					t2 = bd / (2 * A);
				}
				else
				{
					double bd = -B + delta;
					t1 = bd / (2 * A);
					t2 = (2 * C) / bd;
				}

				assert(std::isfinite(t1));
				assert(std::isfinite(t2));

				if (A < 0) std::swap(t1, t2); // make t1 > t2
				minStep = (t1 > 0) ? (std::min(minStep, t2 > 0 ? t2 : t1)) : minStep;
			}
			else
			{
				t1 = -C / B;
				//    avoid divide-by-zero
				minStep = (B == 0) ? minStep : ((t1 > 0) ? std::min(t1, minStep) : minStep);
			}
		}

		return minStep;
	}
}