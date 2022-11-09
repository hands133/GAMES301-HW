#include "Util_DataStructure.h"

#include <numeric>
#include <algorithm>

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
		sigma = SVD.singularValues();

		Eigen::Matrix2d L = Eigen::Matrix2d::Identity();
		if ((U * V.transpose()).determinant() < 0)	L(1, 1) = -1;

		double detU = U.determinant();
		double detV = V.determinant();

		if (detU < 0 && detV > 0)	U = U * L;
		if (detU > 0 && detV < 0)	V = V * L;
		sigma = (sigma.asDiagonal() * L).diagonal();

		R = U * V.transpose();
		S = V * sigma.asDiagonal() * V.transpose();
	}

	QPW_DeformVectors::QPW_DeformVectors(const Eigen::Matrix2d& F, const QPW_Decomposition& decomp)
	{
		Eigen::Matrix2d twist;
		twist << 0, -1, 1, 0;

		Eigen::Matrix2d flip;
		flip << 0, 1, 1, 0;

		auto& U = decomp.U;
		auto& sigma = decomp.sigma;
		auto& VT = decomp.V.transpose();

		double sqrt2 = std::sqrt(2.0);

		//Eigen::Matrix2d G = twist * F * twist.transpose();
		Eigen::Matrix2d G = U * sigma.reverse().asDiagonal() * VT;
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

	QPW_EigenSystem2D::QPW_EigenSystem2D(const QPW_Cell& cell, const QPW_DeformVectors& defvecs)
	{
		const auto& decomp = cell.m_Decomposition;
		const auto& invar = cell.m_Invariables;

		pairs[0].l = 1 + 3 / std::pow(decomp.sigma(0), 4.0);
		pairs[0].e = defvecs.d1;

		pairs[1].l = 1 + 3 / std::pow(decomp.sigma(1), 4.0);
		pairs[1].e = defvecs.d2;
		
		pairs[2].l = 1 + 1 / std::pow(invar.I3, 2.0) + invar.I2 / std::pow(invar.I3, 3.0);
		pairs[2].e = defvecs.l;

		pairs[3].l = 1 + 1 / std::pow(invar.I3, 2.0) - invar.I2 / std::pow(invar.I3, 3.0);
		pairs[3].e = defvecs.t;
	}
}

namespace eigensys
{
	QPW_Cell::QPW_Cell(const Eigen::Matrix2d& dm) : Dm(dm)
	{
		// DmINV = [ [ u, w ]T [ v, t ]T ], pfq_pxq has nothing to do with Ds
		Eigen::Matrix2d DmINV = Dm.inverse();
		double u = DmINV(0, 0);
		double v = DmINV(0, 1);
		double w = DmINV(1, 0);
		double t = DmINV(1, 1);
		double s1 = u + w;
		double s2 = v + t;

		//			pfq  pfq  pfq  pfq  pfq  pfq
		//			---  ---  ---  ---  ---  ---
		//		   px1x px1y px2x px2y px3x px3y
		pfq_pxq << -s1, 0.0,   u, 0.0,    w, 0.0,
				   0.0, -s1, 0.0,    u, 0.0,   w,
				   -s2, 0.0,   v, 0.0,    t, 0.0,
				   0.0, -s2, 0.0,    v, 0.0,   t;
	}

	std::pair<QPW_Cell::vec6d, QPW_Cell::mat6d> QPW_Cell::CalculateGradNHess(const Eigen::Matrix2d& Ds)
	{
		// volume weight
		double volumeWeight = GetVolumeWeight();

		// 1. Prepare F
		Eigen::Matrix2d F = Ds * Dm.inverse();
		
		// 2. Variables (decomposition, invariables, vectors)
		auto defvecs = CalculateVars(F);

		// 3. energy-deform gradient
		auto pPSIq_pfq = CalculateEnergyDeformGrad(m_Invariables, defvecs);
		vec6d GRAD = volumeWeight * pfq_pxq.transpose() * pPSIq_pfq;

		// 4. energy-deform hessian
		QPW_EigenSystem2D eigensys(*this, defvecs);
		auto p2PSIq_pfq2 = CalculateEnergyDeformHess(eigensys);
		mat6d HESS = volumeWeight * pfq_pxq.transpose() * p2PSIq_pfq2 * pfq_pxq;

		return std::make_pair(GRAD, HESS);
	}

	double QPW_Cell::CalculateEnergy(const Eigen::Matrix2d& Ds)
	{
		// 1. Prepare F
		Eigen::Matrix2d F = Ds * Dm.inverse();

		// 2. Variables (decomposition, invariables)
		CalculateVars(F);

		// 3. energy
		return GetVolumeWeight() * CalculateEnergy(m_Invariables);
	}

	double QPW_Cell::CalculateEnergy(const QPW_Invariables& invars) const
	{
		return (invars.I2 + invars.I2 / std::pow(invars.I3, 2.0)) / 2.0;
	}

	QPW_DeformVectors QPW_Cell::CalculateVars(const Eigen::Matrix2d& F)
	{
		m_Decomposition.CalculateSVDnPolar(F);
		m_Invariables.CalculateInvariants(m_Decomposition.S);
		return QPW_DeformVectors(F, m_Decomposition);
	}

	Eigen::Vector4d QPW_Cell::CalculateEnergyDeformGrad(const QPW_Invariables& invars, const QPW_DeformVectors& defvecs) const
	{
		Eigen::Vector4d pPSIq_pfq = Eigen::Vector4d::Zero();

		const auto& invar = m_Invariables;

		double pPSI_pIs[3] = { 0.0,												// pPSI / pI1 = 0
							   (1.0 + 1.0 / std::pow(invar.I3, 2.0)) / 2.0,		// pPSI / pI2 = (1 + 1 / I3^2) / 2
							   - invar.I2 / std::pow(invar.I3, 3.0) };			// pPSI / pI3 = -I2 / I3^3
				
		Eigen::Vector4d pI_pFqs[3] = { defvecs.r,								// pI1 / pfq = r
									   defvecs.f * 2.0,							// pI2 / pfq = 2f
									   defvecs.g };								// pI3 / pfq = g
				
		return std::inner_product(pPSI_pIs, pPSI_pIs + 3, pI_pFqs, Eigen::Vector4d{ Eigen::Vector4d::Zero() });
	}

	Eigen::Matrix4d QPW_Cell::CalculateEnergyDeformHess(const QPW_EigenSystem2D& eigensys) const
	{
		return std::accumulate(eigensys.pairs.begin(), eigensys.pairs.end(), Eigen::Matrix4d{ Eigen::Matrix4d::Zero() },
			[](Eigen::Matrix4d sum, const QPW_EigenValVecPair& pair)
			{	return sum + std::max(pair.l, 0.0) * (pair.e * pair.e.transpose());	});
	}
}