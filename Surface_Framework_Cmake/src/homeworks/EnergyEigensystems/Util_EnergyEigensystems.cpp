#include "Util_EnergyEigensystems.h"

#include <numeric>

#include <TinyAD/Scalar.hh>
#include <TinyAD/Operations/SVD.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/HessianProjection.hh>

#include <Eigen\SparseLU>

namespace eigensys
{
	// axuliary functions
	Eigen::Vector4d opVEC(const Eigen::Matrix2d& M)
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

	inline double UTriangleArea(const Eigen::Vector3i& vertIds, const Eigen::VectorXd& UVs)
	{
		Eigen::Matrix2d Ds;
		Ds << UVs.segment(vertIds(1), 2) - UVs.segment(vertIds(0), 2),
			UVs.segment(vertIds(2), 2) - UVs.segment(vertIds(0), 2);
		return std::abs(Ds.determinant()) / 2.0;
	}

	// member functions
	void QPW_Invariables::CalculateInvariants(const Eigen::Matrix2d& S)
	{
		I1 = S.trace();
		I2 = S.squaredNorm();
		I3 = S.determinant();
	}

	void QPW_Decomposition::CalculateSVDnPolar(const Eigen::Matrix2d& F)
	{
		Eigen::JacobiSVD<Eigen::Matrix2d> SVD(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
		SVD.computeU();
		SVD.computeV();
		U = SVD.matrixU();
		V = SVD.matrixV();
		Sigma = SVD.singularValues().asDiagonal();

		Eigen::Matrix2d L = Eigen::Matrix2d::Identity();
		if ((U * V.transpose()).determinant() < 0)	L(1, 1) = -1;

		double detU = U.determinant();
		double detV = V.determinant();

		if (detU < 0 && detV > 0)	U = U * L;
		if (detU > 0 && detV < 0)	V = V * L;
		Sigma = Sigma * L;

		R = U * V.transpose();
		S = V * Sigma * V.transpose();
	}

	void QPW_DeformGradient::CalculateDeformGradient(const Eigen::Matrix2d& Ds, const Eigen::Matrix2d& DmInv)
	{
		F = Ds * DmInv;
	}

	// DmINV = [ [ u, w ]T [ v, t ]T ], pfq_pxq has nothing to do with Ds
	void QPW_DataPack::Calculatepfq_pxq(const Eigen::Matrix2d& DmInv)
	{
		double u = DmInv(0, 0);
		double v = DmInv(0, 1);
		double w = DmInv(1, 0);
		double t = DmInv(1, 1);
		double s1 = u + w;
		double s2 = v + t;

		//			pfq  pfq  pfq  pfq  pfq  pfq
		//			---  ---  ---  ---  ---  ---
		//		   px1x px1y px2x px2y px3x px3y
		pfq_pxq << -s1, 0.0, u, 0.0, w, 0.0,
				   0.0, -s1, 0.0, u, 0.0, w,
				   -s2, 0.0, v, 0.0, t, 0.0,
				   0.0, -s2, 0.0, v, 0.0, t;
	}

	void QPW_DeformVectors::CalculateVectors(const Eigen::Matrix2d& F, const QPW_Decomposition& decomp)
	{
		Eigen::Matrix2d twist;
		twist << 0, -1, 1, 0;

		Eigen::Matrix2d flip;
		flip << 0, 1, 1, 0;

		auto& U = decomp.U;
		auto& SI = decomp.Sigma;
		auto& VT = decomp.V.transpose();

		double sqrt2 = std::sqrt(2.0);

		//Eigen::Matrix2d G = twist * F * twist.transpose();
		Eigen::Matrix2d G = U * SI.diagonal().reverse().asDiagonal() * VT;
		Eigen::Matrix2d T = U * twist * VT / sqrt2;
		Eigen::Matrix2d P = U * Eigen::Vector2d(1, -1).asDiagonal() * VT / sqrt2;
		Eigen::Matrix2d L = U * flip * VT / sqrt2;
		Eigen::Matrix2d D1 = U * Eigen::Vector2d(1, 0).asDiagonal() * VT;
		Eigen::Matrix2d D2 = U * Eigen::Vector2d(0, 1).asDiagonal() * VT;

		f = opVEC(F);
		g = opVEC(G);
		r = opVEC(decomp.R);
		t = opVEC(T);
		p = opVEC(P);
		l = opVEC(L);
		d1 = opVEC(D1);
		d2 = opVEC(D2);
	}
}

namespace eigensys
{
	QPW_DataPack::QPW_DataPack(const Eigen::Matrix2d& Dm)
	{
		DmINV = Dm.inverse();
		Calculatepfq_pxq(DmINV);
	}

	void QPW_DataPack::CalculatepF_pDs(const Eigen::Matrix2d& Ds)
	{
		m_DeformGradient.CalculateDeformGradient(Ds, DmINV);
		m_Decomposition.CalculateSVDnPolar(m_DeformGradient.F);
		m_Invariables.CalculateInvariants(m_Decomposition.S);
		m_DeformVectors.CalculateVectors(m_DeformGradient.F, m_Decomposition);
	}

	Eigen::MatrixXd QPW_DataPack::GetGlobalpf_px(acamcad::polymesh::PolyMesh* mesh, size_t faceID) const
	{
		size_t numV = mesh->vertices().size();

		Eigen::MatrixXd pf_px(2 * 2, 2 * numV);
		pf_px.setZero();
		auto ids = PolyFaceVertIdxs(mesh, faceID) * 2;

		pf_px.col(ids(0) + 0) = pfq_pxq.col(0);
		pf_px.col(ids(0) + 1) = pfq_pxq.col(1);
		pf_px.col(ids(1) + 0) = pfq_pxq.col(2);
		pf_px.col(ids(1) + 1) = pfq_pxq.col(3);
		pf_px.col(ids(2) + 0) = pfq_pxq.col(4);
		pf_px.col(ids(2) + 1) = pfq_pxq.col(5);

		return pf_px;
	}
}

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
	void ProjectNewtonSolver::PresetMeshUV(acamcad::polymesh::PolyMesh* mesh)
	{
		assert(mesh != nullptr);
		Iters = 0;
		m_FirstUpdate = true;

		size_t numF = mesh->polyfaces().size();
		size_t numV = mesh->vertices().size();

		// parameterization UV
		m_UVList = Eigen::VectorXd(2 * numV);
		m_UVList.setZero();
		for (size_t vertID = 0; vertID < numV; ++vertID)
		{
			acamcad::Texcoord uvw = mesh->vert(vertID)->getTexture();
			m_UVList[2 * vertID + 0] = uvw[0];
			m_UVList[2 * vertID + 1] = uvw[1];
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
			if (v10.cross(v20).dot({ 0, 0, 1 }) < 0)	P2.y() *= -1;

			m_DmList[faceID] << P1 - P0, P2 - P0;
		}
	}

	// Please guarantee that the UV point are stored in the textore coordinate
	bool ProjectNewtonSolver::UpdateMeshUV(acamcad::polymesh::PolyMesh* mesh)
	{
		size_t numV = mesh->vertices().size();
		size_t numF = mesh->polyfaces().size();

		std::cout << "[ProjectNewton] ";

		m_UVList = SetupUVs(mesh);

		auto [G, H] = CalculateGlobalEnergyDerivative(mesh, m_UVList);
		//double GNorm = G.norm();

		//std::cout << "|G| = " << G.norm() << "\t";

		//if (G.norm() < 1.0e-4)
		//{
		//	std::cout << "GMax = " << G.maxCoeff() << " < 1e-4\n";
		//	std::cout << "================CONVERGE!!!================\n";
		//	Iters = 0;
		//	return false;
		//}

		if (m_FirstUpdate)
		{
			m_Energy = CalculateEnergySD_2D(mesh, m_UVList, m_DmList);
			m_FirstUpdate = false;
		}
		
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
		solver.analyzePattern(H);
		solver.factorize(H);
		Eigen::VectorXd d = solver.solve(-G);

		std::cout << "|G| = " << G.norm() << "|d| = " << d.norm() << "\t";

		double dTG = -0.5 * d.dot(G);
		std::cout << "-0.5dot(d,b) = " << dTG << "\t";

		if (dTG < 1.0e-4)
		{
			std::cout << "-0.5dot(d,b) = " << dTG << " < 1.0e-2";
			std::cout << "================CONVERGE!!!================\n";
			Iters = 0;
			return false;
		}

		auto [alpha, energy] = Line_Search(mesh, d, G, m_UVList, m_Energy, 0.9, 1.0e-4);
		assert(!std::isnan(energy) && !std::isinf(energy));

		//double alpha = 1.0;
		//double energy = CalculateEnergySD_2D(mesh, m_UVList, m_DmList);

		std::cout << "Iter = " << ++Iters << " last energy = " << m_Energy 
				  << "\t current energy = " << energy << "\t alpha = " << alpha << "\n";
		m_Energy = energy;

		Eigen::VectorXd uvList = m_UVList + alpha * d;

		//double UVx_min = DBL_MAX;
		//double UVx_max = -DBL_MAX;
		//double UVy_min = DBL_MAX;
		//double UVy_max = -DBL_MAX;

		//for (size_t vertID = 0; vertID < numV; ++vertID)
		//{ 
		//	double x = uvList(2 * vertID);
		//	double y = uvList(2 * vertID + 1);

		//	UVx_min = std::min(UVx_min, x);
		//	UVx_max = std::max(UVx_max, x);

		//	UVy_min = std::min(UVy_min, y);
		//	UVy_max = std::max(UVy_max, y);
		//}

		//double x_scale = UVx_max - UVx_min;
		//double y_scale = UVy_max - UVy_min;

		//double x_scale = 1.0;
		//double y_scale = 1.0;

		for (size_t vertID = 0; vertID < numV; ++vertID)
		{
			auto* vert = mesh->vert(vertID);
			//double UVx = uvList(2 * vertID + 0) / x_scale;
			//double UVy = uvList(2 * vertID + 1) / y_scale;
			double UVx = uvList(2 * vertID + 0);
			double UVy = uvList(2 * vertID + 1);
			vert->setTexture(UVx, UVy, 0.0);
			vert->setPosition(UVx, UVy, 0.0);
		}

		m_UVList = uvList;

		return true;
	}

	QPW_EigenSystem2D ProjectNewtonSolver::Eval_Energy_EigenSystem(const QPW_DataPack& pack) const
	{
		const auto& decomp = pack.m_Decomposition;
		const auto& invars = pack.m_Invariables;

		const auto& U = decomp.U;
		const auto& Sigma = decomp.Sigma;
		const auto& V = decomp.V;

		double I1 = invars.I1;
		double I2 = invars.I2;
		double I3 = invars.I3;

		Eigen::Matrix2d twist;
		twist << 0, -1, 1, 0;

		Eigen::Matrix2d flip;
		flip << 0, 1, 1, 0;

		QPW_EigenSystem2D eigensys;

		Eigen::Matrix2d D1 = U * Eigen::Vector2d(1, 0).asDiagonal() * V.transpose();
		eigensys.valvecpairs[0] = { 1.0 + 3.0 / std::pow(Sigma(0, 0), 4.0), opVEC(D1) };

		Eigen::Matrix2d D2 = U * Eigen::Vector2d(0, 1).asDiagonal() * V.transpose();
		eigensys.valvecpairs[1] = { 1.0 + 3.0 / std::pow(Sigma(1, 1), 4.0), opVEC(D2) };

		Eigen::Matrix2d L = 1.0 / std::sqrt(2.0) * U * flip * V.transpose();
		eigensys.valvecpairs[2] = { 1.0 + 1.0 / std::pow(I3, 2.0) + I2 / std::pow(I3, 3.0), opVEC(L) };

		Eigen::Matrix2d T = 1.0 / std::sqrt(2.0) * U * twist * V.transpose();
		eigensys.valvecpairs[3] = { 1.0 + 1.0 / std::pow(I3, 2.0) - I2 / std::pow(I3, 3.0), opVEC(T) };

		return eigensys;
	}

	Eigen::Matrix4d ProjectNewtonSolver::Calculatep2PSIq_pfq2(const QPW_EigenSystem2D& eigensys) const
	{
		return std::accumulate(eigensys.valvecpairs.begin(), eigensys.valvecpairs.end(), Eigen::Matrix4d{ Eigen::Matrix4d::Zero() },
			[](Eigen::Matrix4d sum, const QPW_EigenValVecPair& pair)
			{	return sum + std::max(pair.l, 0.0) * (pair.e * pair.e.transpose());	});
	}

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

	Eigen::Vector4d ProjectNewtonSolver::CalculatepPSIq_pfq(const QPW_DataPack& pack)
	{
		const auto& defvecs = pack.m_DeformVectors;
		const auto& invar = pack.m_Invariables;

		double pPSI_pI[3] = { 0.0,												// pPSI / pI1 = 0
							(1.0 + 1.0 / (invar.I3 * invar.I3)) / 2.0,			// pPSI / pI2 = (1 + 1 / I3^2) / 2
							- invar.I2 / (invar.I3 * invar.I3 * invar.I3) };	// pPSI / pI3 = -I2 / I3^3
			
		Eigen::Vector4d pI_pFq[3] = { defvecs.r,		// pI1 / pfq = r
									  defvecs.f * 2.0,	// pI2 / pfq = 2f
									  defvecs.g };		// pI3 / pfq = g
			
		return std::inner_product(pPSI_pI, pPSI_pI + 3, pI_pFq, Eigen::Vector4d{ Eigen::Vector4d::Zero() });
	}

	std::pair<double, double> ProjectNewtonSolver::Line_Search(
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

		for (size_t i = 0; i < MAX_LOOP_SIZE; ++i)
		{
			auto updatedUVs = UVs + alpha * d;
			currentEnergy = CalculateEnergySD_2D(mesh, updatedUVs, m_DmList);

			if (currentEnergy <= lastEnergy + c * alpha * d.squaredNorm())	break;
			alpha *= gamma;
		}

		return std::make_pair(alpha, currentEnergy);
	}

	// for 2D symmetric energy Psi = (I2 + I2 / I3^2) / 2
	double ProjectNewtonSolver::CalculateQPWEnergySD_2DNoArea(const Eigen::Matrix2d& DmINV, const Eigen::Matrix2d& Ds) const
	{
		Eigen::Matrix2d F = Ds * DmINV;
		return 0.5 * (F.squaredNorm() + F.inverse().squaredNorm());
	}

	double ProjectNewtonSolver::CalculateQPWEnergySD_2DNoArea(const QPW_Invariables& invars) const
	{
		return (invars.I2 + invars.I2 / (invars.I3 * invars.I3)) / 2.0;
	}

	double ProjectNewtonSolver::CalculateEnergySD_2D(
		acamcad::polymesh::PolyMesh* mesh,
		const Eigen::VectorXd& UVs,
		const std::vector<Eigen::Matrix2d>& DmList) const
	{
		double energy = 0.0;
		size_t numF = mesh->polyfaces().size();
		
		for (size_t faceID = 0; faceID < numF; ++faceID)
		{
			auto ids = PolyFaceVertIdxs(mesh, faceID) * 2;

			Eigen::Matrix2d Ds = Eigen::Matrix2d::Zero();
			Ds << UVs.segment(ids(1), 2) - UVs.segment(ids(0), 2),
				UVs.segment(ids(2), 2) - UVs.segment(ids(0), 2);

			energy += UTriangleArea(ids, UVs) * CalculateQPWEnergySD_2DNoArea(DmList[faceID].inverse(), Ds);
			
			//QPW_DataPack pack(DmList[faceID]);
			//pack.CalculatepF_pDs(Ds);
			//energy += UTriangleArea(ids, UVs) * CalculateQPWEnergySD_2DNoArea(pack.m_Invariables);
		}
		
		return energy;
	}

	std::tuple<Eigen::VectorXd, Eigen::SparseMatrix<double>>
		ProjectNewtonSolver::CalculateGlobalEnergyDerivative(
		acamcad::polymesh::PolyMesh* mesh, const Eigen::VectorXd& UVs)
	{
		int numV = mesh->vertices().size();
		int numF = mesh->polyfaces().size();

		Eigen::VectorXd GRAD(2 * numV);
		GRAD.setZero();

		Eigen::SparseMatrix<double> HESS(2 * numV, 2 * numV);
		HESS.setZero();
		std::vector<Eigen::Triplet<double>> HESStrips;
		HESStrips.reserve(6 * numF);

		for (size_t faceID = 0; faceID < numF; ++faceID)
		{
			auto ids = PolyFaceVertIdxs(mesh, faceID) * 2;
			double volWeight = UTriangleArea(ids, UVs);

			Eigen::Matrix2d Ds;
			Ds << UVs.segment(ids(1), 2) - UVs.segment(ids(0), 2),
				  UVs.segment(ids(2), 2) - UVs.segment(ids(0), 2);

			auto [gradq, hessq] = CalculateLocalEnergyDerivativeNoArea(mesh, m_DmList[faceID], Ds, UVs, faceID);
			int GlobalCoord[6] = { ids(0), ids(0) + 1, ids(1), ids(1) + 1, ids(2), ids(2) + 1 };
			
			for (int m = 0; m < 6; ++m)
			{
				GRAD(GlobalCoord[m]) += volWeight * gradq(m);
				for (int n = 0; n < 6; ++n)
					HESStrips.emplace_back(GlobalCoord[m], GlobalCoord[n], hessq(m, n));
			}
		}

		HESS.setFromTriplets(HESStrips.begin(), HESStrips.end());
		return std::make_tuple(GRAD, HESS);
	}

	// return gradient grad(PSIq) and Hq
	std::tuple<Eigen::Vector<double, 6>, Eigen::Matrix<double, 6, 6>> 
		ProjectNewtonSolver::CalculateLocalEnergyDerivativeNoArea(
			acamcad::polymesh::PolyMesh* mesh,
			const Eigen::Matrix2d& Dm,
			const Eigen::Matrix2d& Ds,
			const Eigen::VectorXd& UVs,
			size_t faceID)
	{
		QPW_DataPack pack(Dm);
		pack.CalculatepF_pDs(Ds);

		auto pfq_pxq = pack.GetLocalpfq_pxq();
		auto pPSIq_pfq = CalculatepPSIq_pfq(pack);

		// gradient
		auto GRADq = pfq_pxq.transpose() * pPSIq_pfq;

		// hessian
		auto eigensys = Eval_Energy_EigenSystem(pack);
		auto p2PSIq_pfq2 = Calculatep2PSIq_pfq2(eigensys);

		auto HESSq = pfq_pxq.transpose() * p2PSIq_pfq2 * pfq_pxq;

		// NOTE p2PSIq_pfq2 positive-semi definite, 
		//{
		//	using LocalDataType = TinyAD::Double<6>;
		//	// use tinyad to make sure the result are correct
		//	auto ids = PolyFaceVertIdxs(mesh, faceID) * 2;
		//	Eigen::Vector<LocalDataType, 6> vvec = LocalDataType::make_active(
		//		{
		//			UVs(ids(0)), UVs(ids(0) + 1),
		//			UVs(ids(1)), UVs(ids(1) + 1),
		//			UVs(ids(2)), UVs(ids(2) + 1)
		//		});
		//	Eigen::Vector2<LocalDataType> P0;
		//	Eigen::Vector2<LocalDataType> P1;
		//	Eigen::Vector2<LocalDataType> P2;

		//	P0 << vvec[0], vvec[1];
		//	P1 << vvec[2], vvec[3];
		//	P2 << vvec[4], vvec[5];

		//	Eigen::Matrix2<LocalDataType> vDs;
		//	vDs << P1 - P0, P2 - P0;
		//	Eigen::Matrix2<LocalDataType> vDm = Dm;

		//	Eigen::Matrix2<LocalDataType> FF = vDs * vDm.inverse();
		//	auto cost = 0.5 * (FF.squaredNorm() + FF.inverse().squaredNorm());
		//	std::cout << "|cost.grad - GRADq| = " << (cost.grad - GRADq).transpose().norm() << std::endl;
		//	return std::make_tuple(cost.grad, cost.Hess);
		//}

		return std::make_tuple(GRADq, HESSq);
	}
}