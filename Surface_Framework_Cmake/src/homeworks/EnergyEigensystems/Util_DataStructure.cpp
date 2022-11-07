#include "Util_DataStructure.h"

namespace eigensys
{
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
	//void QPW_DataPack::Calculatepfq_pxq(const Eigen::Matrix2d& DmInv)
	//{
	//	double u = DmInv(0, 0);
	//	double v = DmInv(0, 1);
	//	double w = DmInv(1, 0);
	//	double t = DmInv(1, 1);
	//	double s1 = u + w;
	//	double s2 = v + t;

	//	//			pfq  pfq  pfq  pfq  pfq  pfq
	//	//			---  ---  ---  ---  ---  ---
	//	//		   px1x px1y px2x px2y px3x px3y
	//	pfq_pxq << -s1, 0.0, u, 0.0, w, 0.0,
	//		0.0, -s1, 0.0, u, 0.0, w,
	//		-s2, 0.0, v, 0.0, t, 0.0,
	//		0.0, -s2, 0.0, v, 0.0, t;
	//}

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
	/*
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
	}*/
}
