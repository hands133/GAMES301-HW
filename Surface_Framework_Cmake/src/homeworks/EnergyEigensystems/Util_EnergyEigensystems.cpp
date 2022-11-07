#include "Util_EnergyEigensystems.h"

#include <numeric>

#include <TinyAD/Scalar.hh>
#include <TinyAD/Operations/SVD.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/HessianProjection.hh>

#include <Eigen\SparseLU>
#include <Eigen\IterativeLinearSolvers>

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

		for (size_t vertID = 0; vertID < numV; ++vertID)
		{
			auto* vert = mesh->vert(vertID);
			float UVx = static_cast<float>(UVs(2 * vertID + 0));
			float UVy = static_cast<float>(UVs(2 * vertID + 1));
			UVx = (UVx - xmin) / (xmax - xmin);
			UVy = (UVy - ymin) / (ymax - ymin);
			vert->setTexture(UVx, UVy, 0.0);
			vert->setPosition(UVx, UVy, 0.0);
		}
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

		m_DmList.clear();
		m_DmList.resize(numF);
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

			m_DmList[faceID] << P1 - P0, P2 - P0;
		}
	}

	// Please guarantee that the UV point are stored in the textore coordinate
	bool ProjectNewtonSolver::UpdateMeshUV(acamcad::polymesh::PolyMesh* mesh)
	{
		size_t numV = mesh->vertices().size();
		size_t numF = mesh->polyfaces().size();

		std::cout << "[ProjectNewton] ";

		auto [G, H] = CalculateEnergyDerivative(mesh, m_UVList);

		std::cout << "max |Gi| = " << G.cwiseAbs().maxCoeff() << "\t";

		if (G.cwiseAbs().maxCoeff() < 1.0e-4)
		{
			std::cout << "================CONVERGE!!!================\n";
			m_Iters = 0;
			return false;
		}

		if (m_FirstUpdate)
		{
			m_Energy = CalculateEnergySD_2D(mesh, m_UVList);
			m_FirstUpdate = false;
		}
		
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		solver.setShift(std::numeric_limits<float>::epsilon());		// necessary
		solver.analyzePattern(H);
		solver.factorize(H);
		Eigen::VectorXd d = solver.solve(-G);

		auto [energy, updatedUV] = Line_Search(mesh, d, G, m_UVList, m_Energy, 0.8, 1.0e-4);

		m_Energy = energy;
		m_UVList = updatedUV;

		ConstrainUV(mesh, updatedUV);

		return true;
	}

	//Eigen::Matrix4d ProjectNewtonSolver::Calculatep2PSIq_pfq2(const QPW_EigenSystem2D& eigensys) const
	//{
	//	return std::accumulate(eigensys.valvecpairs.begin(), eigensys.valvecpairs.end(), Eigen::Matrix4d{ Eigen::Matrix4d::Zero() },
	//		[](Eigen::Matrix4d sum, const QPW_EigenValVecPair& pair)
	//		{	return sum + std::max(pair.l, 0.0) * (pair.e * pair.e.transpose());	});
	//}

	//// for debug, equal to sum(max(lambda, 0) ei eiT)
	//Eigen::Matrix4d ProjectNewtonSolver::CalculateHqAnother(const QPW_pfq_pxq& pfq_pxq) const
	//{
	//	const auto& defvecs = pfq_pxq.m_DeformVectors;
	//	const auto& invar = pfq_pxq.m_Invariables;

	//	auto rhat = defvecs.r / std::sqrt(2.0);

	//	Eigen::Matrix4d P1 = (1.0 + 1.0 / (invar.I3 * invar.I3)) * Eigen::Matrix4d::Identity();
	//	Eigen::Matrix4d P2 = -invar.I2 / std::pow(invar.I3, 3.0) * (
	//		rhat * rhat.transpose() + defvecs.t * defvecs.t.transpose() -
	//		defvecs.p * defvecs.p.transpose() - defvecs.l * defvecs.l.transpose());
	//	Eigen::Matrix4d P3 = -2.0 / std::pow(invar.I3, 3.0) * (
	//		defvecs.g * defvecs.f.transpose() + defvecs.f * defvecs.g.transpose());

	//	return P1 + P2 + P3;
	//}

	//Eigen::Vector4d ProjectNewtonSolver::CalculatepPSIq_pfq(const QPW_DataPack& pack)
	//{
	//	const auto& defvecs = pack.m_DeformVectors;
	//	const auto& invar = pack.m_Invariables;

	//	double pPSI_pI[3] = { 0.0,												// pPSI / pI1 = 0
	//						(1.0 + 1.0 / (invar.I3 * invar.I3)) / 2.0,			// pPSI / pI2 = (1 + 1 / I3^2) / 2
	//						- invar.I2 / (invar.I3 * invar.I3 * invar.I3) };	// pPSI / pI3 = -I2 / I3^3
	//		
	//	Eigen::Vector4d pI_pFq[3] = { defvecs.r,		// pI1 / pfq = r
	//								  defvecs.f * 2.0,	// pI2 / pfq = 2f
	//								  defvecs.g };		// pI3 / pfq = g
	//		
	//	return std::inner_product(pPSI_pI, pPSI_pI + 3, pI_pFq, Eigen::Vector4d{ Eigen::Vector4d::Zero() });
	//}

	std::pair<double, Eigen::VectorXd> ProjectNewtonSolver::Line_Search(
		acamcad::polymesh::PolyMesh* mesh,
		const Eigen::VectorXd& d,
		const Eigen::VectorXd& grad,
		const Eigen::VectorXd& UVs,
		double lastEnergy,
		double gamma, double c)
	{
		double alpha = 1.0;
		double currentEnergy = lastEnergy;

		const int MAX_LOOP_SIZE = 64;

		auto updatedUVs = UVs;
		size_t iter = 0;
		for (iter = 0; iter < MAX_LOOP_SIZE; ++iter)
		{
			updatedUVs = UVs + alpha * d;
			currentEnergy = CalculateEnergySD_2D(mesh, updatedUVs);

			if (currentEnergy <= lastEnergy + c * alpha * d.dot(grad))
			{
				std::cout << "Iter = " << ++m_Iters << "\tlast E = " << lastEnergy
					<< "\tcurrent E = " << currentEnergy << "\talpha = " << alpha << "\n";
				break;
			}
			alpha *= gamma;
		}
		return std::make_pair(currentEnergy, updatedUVs);
	}

	// for 2D symmetric energy Psi = (I2 + I2 / I3^2) / 2
	double ProjectNewtonSolver::QPW_CalculateEnergySD_2D(const Eigen::Matrix2d& DmINV, const Eigen::Matrix2d& Ds) const
	{
		Eigen::Matrix2d F = Ds * DmINV;
		return 0.5 * (F.squaredNorm() + F.inverse().squaredNorm());
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

			const auto& Dm = m_DmList[faceID];
			double volWeight = Dm.determinant() / 2.0;
			energy += volWeight * QPW_CalculateEnergySD_2D(Dm.inverse(), Ds);
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

			Eigen::Matrix2d Ds;
			Ds << UVs.segment(ids(1), 2) - UVs.segment(ids(0), 2),
				  UVs.segment(ids(2), 2) - UVs.segment(ids(0), 2);

			auto [gradq, hessq] = QPW_CalculateEnergyDerivative(mesh, m_DmList[faceID], Ds);

			int GlobalCoord[6] = { ids(0), ids(0) + 1, ids(1), ids(1) + 1, ids(2), ids(2) + 1 };
			
			for (int n = 0; n < 6; ++n)
			{
				GRAD(GlobalCoord[n]) += gradq(n);
				for (int m = 0; m < 6; ++m)
					HESStrips.emplace_back(GlobalCoord[n], GlobalCoord[m], hessq(n, m));
			}
			assert(!std::isnan(GRAD.norm()));
		}

		HESS.setFromTriplets(HESStrips.begin(), HESStrips.end());
		return std::make_tuple(GRAD, HESS);
	}

	std::tuple<Eigen::Vector<double, 6>, Eigen::Matrix<double, 6, 6>> 
		ProjectNewtonSolver::QPW_CalculateEnergyDerivative(
		acamcad::polymesh::PolyMesh* mesh, const Eigen::Matrix2d& Dm, const Eigen::Matrix2d& Ds)
	{
		Eigen::Matrix2d DmInv = Dm.inverse();
		Eigen::Matrix2d F = Ds * DmInv;

		Eigen::Matrix<double, 4, 6> pfq_pxq;
		pfq_pxq.setZero();

		{	// derivative
			double u = DmInv(0, 0);
			double v = DmInv(0, 1);
			double w = DmInv(1, 0);
			double t = DmInv(1, 1);
			double s1 = u + w;
			double s2 = v + t;

			//			pfq  pfq  pfq  pfq  pfq  pfq
			//			---  ---  ---  ---  ---  ---
			//		   px1x px1y px2x px2y px3x px3y
			pfq_pxq << -s1, 0.0,   u, 0.0,   w, 0.0,
					   0.0, -s1, 0.0,   u, 0.0,   w,
					   -s2, 0.0,   v, 0.0,   t, 0.0,
					   0.0, -s2, 0.0,   v, 0.0,   t;
		}

		Eigen::Vector4d pPSIq_pfq;
		Eigen::Matrix4d p2PSIq_pfq2;

		Eigen::Matrix2d U, V;
		Eigen::Vector2d sigma;

		{	// PSI(x) derivatives
			Eigen::JacobiSVD<Eigen::Matrix2d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
			svd.computeU();
			svd.computeV();
			U = svd.matrixU();
			V = svd.matrixV();
			sigma = svd.singularValues();

			{	// flip correction
				Eigen::Matrix2d L = Eigen::Matrix2d::Identity();
				if ((U * V.transpose()).determinant() < 0)	L(1, 1) = -1;

				double detU = U.determinant();
				double detV = V.determinant();

				if (detU < 0 && detV > 0)	U = U * L;
				if (detU > 0 && detV < 0)	V = V * L;
				sigma = (sigma.asDiagonal() * L).diagonal();
			}

			Eigen::Matrix2d R = U * V.transpose();
			Eigen::Matrix2d S = V * sigma.asDiagonal() * V.transpose();

			double sigma1 = sigma(0);
			double sigma2 = sigma(1);

			Eigen::Matrix2d twist;
			twist << 0, -1, 1, 0;
			Eigen::Matrix2d flip;
			flip << 0, 1, 1, 0;

			double I2 = S.squaredNorm();
			double I3 = S.determinant();

			double lam1 = 1 + 3 / std::pow(sigma1, 4.0);
			double lam2 = 1 + 3 / std::pow(sigma2, 4.0);
			double lam3 = 1 + 1 / std::pow(I3, 2.0) + I2 / std::pow(I3, 3.0);
			double lam4 = 1 + 1 / std::pow(I3, 2.0) - I2 / std::pow(I3, 3.0);

			lam3 = std::max(lam3, 0.0);
			lam4 = std::max(lam4, 0.0);

			Eigen::Matrix2d D1 = U * Eigen::Vector2d(1, 0).asDiagonal() * V.transpose();
			Eigen::Matrix2d D2 = U * Eigen::Vector2d(0, 1).asDiagonal() * V.transpose();
			Eigen::Matrix2d L = 1 / sqrt(2) * U * flip * V.transpose();
			Eigen::Matrix2d T = 1 / sqrt(2) * U * twist * V.transpose();

			//Eigen::Matrix2d G = twist * F * twist.transpose();
			Eigen::Matrix2d G = U * sigma.reverse().asDiagonal() * V.transpose();
			pPSIq_pfq = (1 + 1 / (I3 * I3)) * opVEC(F) - I2 / (I3 * I3 * I3) * opVEC(G);
			p2PSIq_pfq2 = lam1 * opVEC(D1) * opVEC(D1).transpose() +
						  lam2 * opVEC(D2) * opVEC(D2).transpose() +
						  lam3 * opVEC(L) * opVEC(L).transpose() +
						  lam4 * opVEC(T) * opVEC(T).transpose();
		}

		double volumeWeight = Dm.determinant();
		
		Eigen::Vector<double, 6> GRADq = volumeWeight * pfq_pxq.transpose() * pPSIq_pfq;
		assert(!std::isnan(GRADq.norm()));

		Eigen::Matrix<double, 6, 6> HESSq = volumeWeight * pfq_pxq.transpose() * p2PSIq_pfq2 * pfq_pxq;

		return std::make_pair(GRADq, HESSq);
	}
}