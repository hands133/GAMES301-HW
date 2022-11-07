#include "Util_EnergyEigensystems.h"

#include <numeric>

#include <TinyAD/Scalar.hh>
#include <TinyAD/Operations/SVD.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/HessianProjection.hh>

#include <Eigen\SparseLU>
#include <Eigen\IterativeLinearSolvers>

//#include <fstream>
//#include <string>

namespace eigensys
{
	// only used for UpdateMeshUV
	Eigen::VectorXd ProjectNewtonSolver::SetupUVs(acamcad::polymesh::PolyMesh* mesh)
	{
		size_t numV = mesh->vertices().size();

		Eigen::VectorXd UVs(2 * numV);
		UVs.setZero();
		for (size_t vertID = 0; vertID < numV; ++vertID)
		{
			acamcad::Texcoord uvw = mesh->vert(vertID)->getTexture();
			UVs[2 * vertID + 0] = uvw[0];
			UVs[2 * vertID + 1] = uvw[1];
		}
		return UVs;
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

			//m_DmList[faceID] << Eigen::Vector2d(v10.norm(), 0),
			//	Eigen::Vector2d(v10.dot(v20), v10.cross(v20).norm()) / v10.norm();
		}
	}

	// Please guarantee that the UV point are stored in the textore coordinate
	bool ProjectNewtonSolver::UpdateMeshUV(acamcad::polymesh::PolyMesh* mesh)
	{
		size_t numV = mesh->vertices().size();
		size_t numF = mesh->polyfaces().size();

		std::cout << "[ProjectNewton] ";

		//m_UVList = SetupUVs(mesh);

		auto [G, H] = CalculateGlobalEnergyDerivative(mesh, m_UVList);

		//std::ofstream Hfile("H.txt");
		//std::ofstream bfile("b.txt");

		//for (int k = 0; k < H.outerSize(); ++k)
		//	for (Eigen::SparseMatrix<double>::InnerIterator it(H, k); it; ++it)
		//	{
		//		it.value();
		//		it.row();   // row index
		//		it.col();   // col index (here it is equal to k)
		//		it.index(); // inner index, here it is equal to it.row()
		//		Hfile << it.row() << "\t" << it.col() << "\t" << it.value() << "\n";
		//	}

		//for (size_t vertID = 0; vertID < numV; ++vertID)
		//	bfile << G(vertID * 2) << "\n" << G(vertID * 2 + 1) << "\n";

		//Hfile.close();
		//bfile.close();

		//std::cout << "max |Gi| = " << G.cwiseAbs().maxCoeff() << "\t";

		//if (G.cwiseAbs().maxCoeff() < 1.0e-4)
		//{
		//	std::cout << "================CONVERGE!!!================\n";
		//	Iters = 0;
		//	return false;
		//}

		if (m_FirstUpdate)
		{
			m_Energy = CalculateEnergySD_2D(mesh, m_UVList);
			m_FirstUpdate = false;
		}
		
		//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
		//Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
		//Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteCholesky<double>> solver;
		//solver.setTolerance(std::numeric_limits<float>::epsilon());
		//solver.setMaxIterations(2000);
		//Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		solver.setShift(std::numeric_limits<float>::epsilon());
		solver.analyzePattern(H);
		solver.factorize(H);
		Eigen::VectorXd d = solver.solve(-G);


		std::cout << "|G| = " << G.norm() << "\t|d| = " << d.norm() << "\t";

		double dTG = -0.5 * d.dot(G);
		std::cout << "-0.5dot(d,b) = " << dTG << "\t";

		if (dTG < 1.0e-4)
		{
			std::cout << "================CONVERGE!!!================\n";
			m_Iters = 0;
			return false;
		}

		auto [energy, updatedUV] = Line_Search(mesh, d, G, m_UVList, m_Energy, 0.8, 1.0e-4);
		assert(!std::isnan(energy) && !std::isinf(energy));

		m_Energy = energy;
		m_UVList = updatedUV;

		for (size_t vertID = 0; vertID < numV; ++vertID)
		{
			auto* vert = mesh->vert(vertID);
			float UVx = static_cast<float>(updatedUV(2 * vertID + 0));
			float UVy = static_cast<float>(updatedUV(2 * vertID + 1));
			vert->setTexture(UVx, UVy, 0.0);
			vert->setPosition(UVx, UVy, 0.0);
		}

		return true;
	}

	//QPW_EigenSystem2D ProjectNewtonSolver::Eval_Energy_EigenSystem(const QPW_DataPack& pack) const
	//{
	//	const auto& decomp = pack.m_Decomposition;
	//	const auto& invars = pack.m_Invariables;

	//	const auto& U = decomp.U;
	//	const auto& Sigma = decomp.Sigma;
	//	const auto& V = decomp.V;

	//	double I1 = invars.I1;
	//	double I2 = invars.I2;
	//	double I3 = invars.I3;

	//	Eigen::Matrix2d twist;
	//	twist << 0, -1, 1, 0;

	//	Eigen::Matrix2d flip;
	//	flip << 0, 1, 1, 0;

	//	QPW_EigenSystem2D eigensys;

	//	Eigen::Matrix2d D1 = U * Eigen::Vector2d(1, 0).asDiagonal() * V.transpose();
	//	eigensys.valvecpairs[0] = { 1.0 + 3.0 / std::pow(Sigma(0, 0), 4.0), opVEC(D1) };

	//	Eigen::Matrix2d D2 = U * Eigen::Vector2d(0, 1).asDiagonal() * V.transpose();
	//	eigensys.valvecpairs[1] = { 1.0 + 3.0 / std::pow(Sigma(1, 1), 4.0), opVEC(D2) };

	//	Eigen::Matrix2d L = 1.0 / std::sqrt(2.0) * U * flip * V.transpose();
	//	eigensys.valvecpairs[2] = { 1.0 + 1.0 / std::pow(I3, 2.0) + I2 / std::pow(I3, 3.0), opVEC(L) };

	//	Eigen::Matrix2d T = 1.0 / std::sqrt(2.0) * U * twist * V.transpose();
	//	eigensys.valvecpairs[3] = { 1.0 + 1.0 / std::pow(I3, 2.0) - I2 / std::pow(I3, 3.0), opVEC(T) };

	//	return eigensys;
	//}

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
		ProjectNewtonSolver::CalculateGlobalEnergyDerivative(
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

			Eigen::Vector2d P0 = UVs.segment(ids(0), 2);
			Eigen::Vector2d P1 = UVs.segment(ids(1), 2);
			Eigen::Vector2d P2 = UVs.segment(ids(2), 2);

			acamcad::Texcoord C0 = mesh->vert(ids(0) / 2)->getTexture();
			acamcad::Texcoord C1 = mesh->vert(ids(1) / 2)->getTexture();
			acamcad::Texcoord C2 = mesh->vert(ids(2) / 2)->getTexture();

			Eigen::Matrix2d Ds;
			//Ds << UVs.segment(ids(1), 2) - UVs.segment(ids(0), 2),
			//	  UVs.segment(ids(2), 2) - UVs.segment(ids(0), 2);
			Ds << P1 - P0, P2 - P0;

			auto [gradq, hessq] = CalculateLocalEnergyDerivative(mesh, m_DmList[faceID], Ds);

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

	// return gradient grad(PSIq) and Hq
	//std::tuple<Eigen::Vector<double, 6>, Eigen::Matrix<double, 6, 6>> 
	//	ProjectNewtonSolver::CalculateLocalEnergyDerivativeNoArea(
	//		acamcad::polymesh::PolyMesh* mesh,
	//		const Eigen::Matrix2d& Dm,
	//		const Eigen::Matrix2d& Ds,
	//		const Eigen::VectorXd& UVs,
	//		size_t faceID)
	//{
	//	QPW_DataPack pack(Dm);
	//	pack.CalculatepF_pDs(Ds);

	//	auto pfq_pxq = pack.GetLocalpfq_pxq();
	//	auto pPSIq_pfq = CalculatepPSIq_pfq(pack);

	//	// gradient
	//	auto GRADq = pfq_pxq.transpose() * pPSIq_pfq;

	//	// hessian
	//	auto eigensys = Eval_Energy_EigenSystem(pack);
	//	auto p2PSIq_pfq2 = Calculatep2PSIq_pfq2(eigensys);

	//	auto HESSq = pfq_pxq.transpose() * p2PSIq_pfq2 * pfq_pxq;

	//	// NOTE p2PSIq_pfq2 positive-semi definite, 
	//	//{
	//	//	using LocalDataType = TinyAD::Double<6>;
	//	//	// use tinyad to make sure the result are correct
	//	//	auto ids = PolyFaceVertIdxs(mesh, faceID) * 2;
	//	//	Eigen::Vector<LocalDataType, 6> vvec = LocalDataType::make_active(
	//	//		{
	//	//			UVs(ids(0)), UVs(ids(0) + 1),
	//	//			UVs(ids(1)), UVs(ids(1) + 1),
	//	//			UVs(ids(2)), UVs(ids(2) + 1)
	//	//		});
	//	//	Eigen::Vector2<LocalDataType> P0;
	//	//	Eigen::Vector2<LocalDataType> P1;
	//	//	Eigen::Vector2<LocalDataType> P2;

	//	//	P0 << vvec[0], vvec[1];
	//	//	P1 << vvec[2], vvec[3];
	//	//	P2 << vvec[4], vvec[5];

	//	//	Eigen::Matrix2<LocalDataType> vDs;
	//	//	vDs << P1 - P0, P2 - P0;
	//	//	Eigen::Matrix2<LocalDataType> vDm = Dm;

	//	//	Eigen::Matrix2<LocalDataType> FF = vDs * vDm.inverse();
	//	//	auto cost = 0.5 * (FF.squaredNorm() + FF.inverse().squaredNorm());
	//	//	std::cout << "|cost.grad - GRADq| = " << (cost.grad - GRADq).transpose().norm() << std::endl;
	//	//	return std::make_tuple(cost.grad, cost.Hess);
	//	//}

	//	return std::make_tuple(GRADq, HESSq);
	//}

	std::tuple<Eigen::Vector<double, 6>, Eigen::Matrix<double, 6, 6>> 
		ProjectNewtonSolver::CalculateLocalEnergyDerivative(
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
				Eigen::Matrix2d transsigma = sigma.asDiagonal() * L;
				sigma = transsigma.diagonal();
			}

			Eigen::Matrix2d twist;
			twist << 0, -1, 1, 0;
			Eigen::Matrix2d flip;
			flip << 0, 1, 1, 0;

			double I2 = F.squaredNorm();
			double I3 = F.determinant();
			double lam1 = 1 + 3 / (sigma[0] * sigma[0] * sigma[0] * sigma[0]);
			double lam2 = 1 + 3 / (sigma[1] * sigma[1] * sigma[1] * sigma[1]);
			double lam3 = 1 + 1 / (I3 * I3) + I2 / (I3 * I3 * I3);
			double lam4 = 1 + 1 / (I3 * I3) - I2 / (I3 * I3 * I3);

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