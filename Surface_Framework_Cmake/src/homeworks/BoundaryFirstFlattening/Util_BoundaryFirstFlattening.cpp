#include "Util_BoundaryFirstFlattening.h"

#include "TutteEmbedding/Util_TutteEmbedding.h"

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

namespace bff
{
	void BFFSolver::Solve(acamcad::polymesh::PolyMesh* mesh, FlattenType type)
	{
		assert(!mesh->isEmpty(), "mesh is empty!");
		m_Mesh = mesh;

		InitParam();
		
		switch (type)
		{
		case bff::FlattenType::FREE:
		{
			m_BoundaryScaleFactor = Eigen::VectorXd::Zero(m_NumBV);
			Flatten(m_BoundaryScaleFactor, true);
		}
			break;
		case bff::FlattenType::DISK:
		{
			m_FlattenToDisk = true;
			FlattenToDisk();
		}
			break;
		case bff::FlattenType::FIXED:
		{
			assert(m_TargetExteriorAngle.rows() == m_NumBV, "Require Input Exterior Angle!");
			// rearrange vert index
			Eigen::VectorXd tmpBoundaryExteriorAngle = m_TargetExteriorAngle;
			m_TargetExteriorAngle = Eigen::VectorXd::Zero(m_NumBV);
			for (size_t bVertID = 0; bVertID < m_NumBV; ++bVertID)
			{
				size_t bID = Util_VertIndex(m_BoundaryVertices[bVertID]) - m_NumIV;
				m_TargetExteriorAngle(bID) = tmpBoundaryExteriorAngle(bVertID);
			}
			Flatten(m_TargetExteriorAngle, false);
		}
			break;
		}
	}

	void BFFSolver::SetExteriorAngle(size_t numBV, tutte::UVBoundaryType shape)
	{
		m_TargetExteriorAngle = Eigen::VectorXd::Zero(numBV);
		auto targetExteriorBoundaryUVs = tutte::GetBoundaryUVs(numBV, shape);
		for (size_t j = 0; j < numBV; ++j)
		{
			size_t i = (j + numBV - 1) % numBV;
			size_t k = (j + 1) % numBV;

			auto vecij = targetExteriorBoundaryUVs[j] - targetExteriorBoundaryUVs[i];
			auto vecjk = targetExteriorBoundaryUVs[k] - targetExteriorBoundaryUVs[j];

			m_TargetExteriorAngle(j) = acamcad::vectorAngle(acamcad::MVector3(vecij.x(), vecij.y(), 0.0),
															acamcad::MVector3(vecjk.x(), vecjk.y(), 0.0));
		}
	}

	Eigen::VectorXd BFFSolver::CalculateBoundaryLengths(const std::vector<acamcad::polymesh::MVert*>& bVertices) const
	{
		size_t numBV = bVertices.size();
		Eigen::VectorXd edges = Eigen::VectorXd(numBV);
		edges.setZero();

		for (size_t bVertID = 0; bVertID < numBV; ++bVertID)
		{
			auto* pBVert0 = bVertices[bVertID];
			auto* pBVert1 = bVertices[(bVertID + 1) % m_NumBV];
			size_t bID = Util_VertIndex(pBVert0) - m_NumIV;
			edges(bID) = acamcad::distance(pBVert0->position(), pBVert1->position());
		}
		return edges;
	}

	std::tuple<Eigen::VectorXd, Eigen::VectorXd> BFFSolver::CalculateDiscreteCurvatures() const
	{
		Eigen::VectorXd gaussianCurv = Eigen::VectorXd::Zero(m_NumIV);
		Eigen::VectorXd geodesicCurv = Eigen::VectorXd::Zero(m_NumBV);

		for (auto* pVert : m_Mesh->vertices())
		{
			size_t iVertID = Util_VertIndex(pVert);
			if (!m_Mesh->isBoundary(pVert))
				gaussianCurv(iVertID) = Util_VertAngleDefect(pVert);
		}

		for (auto* pBVert : m_BoundaryVertices)
		{
			size_t bID = Util_VertIndex(pBVert) - m_NumIV;
			geodesicCurv(bID) = Util_VertExteriorAngle(pBVert);
		}

		return std::make_tuple(gaussianCurv, geodesicCurv);
	}

	Eigen::SparseMatrix<double> BFFSolver::BuildLaplacian()
	{
		std::vector<Eigen::Triplet<double>> ATriplets;
		ATriplets.reserve(m_Mesh->halfEdges().size());

		for (auto* pFace : m_Mesh->polyfaces())
		{
			acamcad::polymesh::MHalfedge* he = pFace->halfEdge();
			do
			{
				size_t i = Util_VertIndex(he->fromVertex());
				size_t j = Util_VertIndex(he->toVertex());
				double w = 0.5 * Util_HalfEdgeAngleCotan(he);

				ATriplets.emplace_back(i, i, w);
				ATriplets.emplace_back(j, j, w);
				ATriplets.emplace_back(i, j, -w);
				ATriplets.emplace_back(j, i, -w);

				he = he->next();
			} while (he != pFace->halfEdge());
		}

		Eigen::SparseMatrix<double> A = Eigen::SparseMatrix<double>(m_NumV, m_NumV);
		A.setFromTriplets(ATriplets.begin(), ATriplets.end());
		A.makeCompressed();

		auto AShift = Eigen::SparseMatrix<double>(m_NumV, m_NumV);
		AShift.setIdentity();
		A += AShift * std::numeric_limits<float>::epsilon();

		return A;
	}

	Eigen::MatrixXd BFFSolver::ConstructBestFitCurve(const Eigen::VectorXd& l,
												     const Eigen::VectorXd& lstar,
												     const Eigen::VectorXd& ktilde) const
	{
		Eigen::MatrixXd T_tilde = Eigen::MatrixXd::Zero(2, m_NumBV);
		Eigen::VectorXd Ndiag = Eigen::VectorXd::Zero(m_NumBV);
		double angle = 0.0;

		for (size_t bVertID = 0; bVertID < m_NumBV; ++bVertID)
		{
			size_t bIDI = Util_VertIndex(m_BoundaryVertices[bVertID]) - m_NumIV;

			angle += ktilde(bIDI);
			T_tilde.col(bIDI) = Eigen::Vector2d(std::cos(angle), std::sin(angle));
			Ndiag(bIDI) = 1.0 / l(bIDI);
		}
		
		Eigen::MatrixXd N = Ndiag.asDiagonal();
		Eigen::VectorXd ltilde = lstar - N.inverse() * T_tilde.transpose() * (T_tilde * N.inverse() * T_tilde.transpose()).inverse() * T_tilde * lstar;

		Eigen::MatrixXd gamma_tilde = Eigen::MatrixXd::Zero(2, m_NumBV);

		for (size_t bVertID = 1; bVertID < m_NumBV; ++bVertID)
		{
			size_t ID1 = Util_VertIndex(m_BoundaryVertices[bVertID]) - m_NumIV;
			size_t ID0 = Util_VertIndex(m_BoundaryVertices[bVertID - 1]) - m_NumIV;

			gamma_tilde.col(ID1) = gamma_tilde.col(ID0) + ltilde(ID0) * T_tilde.col(ID0);
		}

		return gamma_tilde;
	}

	void BFFSolver::InitParam()
	{
		m_BoundaryVertices = m_Mesh->boundaryVertices();

		m_NumV = m_Mesh->vertices().size();
		m_NumBV = m_BoundaryVertices.size();
		m_NumIV = m_NumV - m_NumBV;

		Util_RearrangeVertIndex();

		m_BoundaryEdgeLengths = CalculateBoundaryLengths(m_BoundaryVertices);

		{
			auto& [Gk, Ck] = CalculateDiscreteCurvatures();
			m_GaussianCurv = Gk;
			m_GeodesicCurv = Ck;
		}

		m_LapMat = BuildLaplacian();
		m_LapII = m_LapMat.block(0, 0, m_NumIV, m_NumIV);
		m_LapIB = m_LapMat.block(0, m_NumIV, m_NumIV, m_NumBV);
		m_LapBB = m_LapMat.block(m_NumIV, m_NumIV, m_NumBV, m_NumBV);

		m_LUSolver.compute(m_LapMat);
	}

	void BFFSolver::Flatten(const Eigen::VectorXd& givenTarget,
							bool givenScaleFactor)
	{
		if (givenScaleFactor)
		{
			Eigen::VectorXd dudn = Util_CvtDirichletToNeumann(-m_GaussianCurv, givenTarget);
			Eigen::VectorXd tExteriorAngle = m_GeodesicCurv - dudn;
			double EPS = std::abs(tExteriorAngle.sum() - M_PI * 2);
			std::cout << "Exterior Angle EPS = " << EPS << "\n";

			Flatten(givenTarget, tExteriorAngle, true);
		}
		else
		{
			m_BoundaryScaleFactor = Util_CvtNeumannToDirichlet(-m_GaussianCurv, m_GeodesicCurv - givenTarget);
			Flatten(m_BoundaryScaleFactor, givenTarget, false);
		}

		Util_ConstrainUV();
	}

	void BFFSolver::Flatten(const Eigen::VectorXd& scaleFactor,
							const Eigen::VectorXd& exteriorAngle,
							bool freeBoundary)
	{
		Eigen::VectorXd tmpBoundaryLengths = Eigen::VectorXd::Zero(m_NumBV);
		Util_TargetBoundaryLengths(tmpBoundaryLengths);

		Eigen::MatrixXd gamma_tilde = ConstructBestFitCurve(m_BoundaryEdgeLengths, tmpBoundaryLengths, exteriorAngle);

		auto& [a, b] = ExtendCurve(gamma_tilde, freeBoundary);

		m_UVList = Eigen::MatrixXd(m_NumV, 2);
		const auto& vertices = m_Mesh->vertices();
		for (size_t vertID = 0; vertID < m_NumV; ++vertID)
		{
			size_t orderedVertID = Util_VertIndex(m_Mesh->vert(vertID));
			m_UVList(vertID, 0) = a(orderedVertID);
			m_UVList(vertID, 1) = -b(orderedVertID);
		}
	}

	void BFFSolver::FlattenToDisk()
	{
		m_BoundaryScaleFactor = Eigen::VectorXd::Zero(m_NumBV);
		m_TargetExteriorAngle = Eigen::VectorXd::Zero(m_NumBV);
		
		m_TargetEdgeLengths = Eigen::VectorXd::Zero(m_NumBV);
		
		Eigen::VectorXd tDualLength = Eigen::VectorXd::Zero(m_NumBV);
		for (size_t repeat = 0; repeat < 20; ++repeat)
		{
			Util_TargetBoundaryLengths(m_TargetEdgeLengths);
			double L = Util_TargetDualBoundaryLengths(m_TargetEdgeLengths, tDualLength);

			for (auto* pBVert : m_BoundaryVertices)
			{
				size_t bID = Util_VertIndex(pBVert) - m_NumIV;
				m_TargetExteriorAngle(bID) = 2 * M_PI * tDualLength(bID) / L;
			}
			m_BoundaryScaleFactor = Util_CvtNeumannToDirichlet(-m_GaussianCurv, m_GeodesicCurv - m_TargetExteriorAngle);
		}

		Flatten(m_TargetExteriorAngle, false);
	}

	void BFFSolver::Util_RearrangeVertIndex()
	{
		m_VertIndexMap.clear();
		m_VertIndexMap.reserve(m_NumV);

		size_t ivertCount = 0;

		for (auto* pVert : m_Mesh->vertices())
		{
			if (m_Mesh->isBoundary(pVert))	continue;
			m_VertIndexMap.insert({ pVert, ivertCount++ });
		}

		std::for_each(m_BoundaryVertices.begin(), m_BoundaryVertices.end(), 
			[&](acamcad::polymesh::MVert* v) { m_VertIndexMap.insert({ v, ivertCount++ }); });
	}

	// Integrated gaussian curvature for inner vertices
	double BFFSolver::Util_VertAngleDefect(acamcad::polymesh::MVert* vert) const
	{
		const auto& neighbors = m_Mesh->vertAdjacentVertices(vert);
		size_t numNV = neighbors.size();

		double angleSum{ 0.0 };
		for (size_t vertID = 0; vertID < numNV; ++vertID)
		{
			size_t i = vertID;
			size_t j = (i + 1) % numNV;
			angleSum += std::abs(acamcad::vectorAngle(
				neighbors[i]->position() - vert->position(),
				neighbors[j]->position() - vert->position()));
		}
		return 2 * M_PI - angleSum;
	}

	double BFFSolver::Util_VertExteriorAngle(acamcad::polymesh::MVert* vert) const
	{
		const auto& neighbors = m_Mesh->vertAdjacentVertices(vert);
		size_t numNV = neighbors.size();

		double angleSum{ 0.0 };
		for (size_t vertID = 0; vertID < numNV - 1; ++vertID)
		{
			size_t i = vertID;
			size_t j = i + 1;
			angleSum += std::abs(acamcad::vectorAngle(
				neighbors[i]->position() - vert->position(),
				neighbors[j]->position() - vert->position()));
		}
		return M_PI - angleSum;
	}

	double BFFSolver::Util_HalfEdgeAngleCotan(acamcad::polymesh::MHalfedge* he) const
	{
		auto* pVertA = he->fromVertex();
		auto* pVertB = he->next()->fromVertex();
		auto* pVertC = he->prev()->fromVertex();

		auto vecU = pVertA->position() - pVertC->position();
		auto vecV = pVertB->position() - pVertC->position();

		return acamcad::dot(vecU, vecV) / acamcad::cross(vecU, vecV).norm();
	}

	Eigen::VectorXd BFFSolver::Util_CvtDirichletToNeumann(const Eigen::VectorXd& phi,
														  const Eigen::VectorXd& g) const
	{
		Eigen::VectorXd b = phi - m_LapIB * g;
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldltSolver;
		ldltSolver.compute(m_LapII);
		Eigen::VectorXd a = ldltSolver.solve(b);

		return - (m_LapIB.transpose() * a - m_LapBB * g);
	}
	
	Eigen::VectorXd BFFSolver::Util_CvtNeumannToDirichlet(const Eigen::VectorXd& phi,
														  const Eigen::VectorXd& h)
	{
		Eigen::VectorXd b = Util_VecVertConcat(phi, -h);
		Eigen::VectorXd a = m_LUSolver.solve(b);

		return a.segment(m_NumIV, m_NumBV);
	}

	Eigen::VectorXd BFFSolver::Util_VecVertConcat(const Eigen::VectorXd& A, const Eigen::VectorXd& B) const
	{
		size_t mA = A.size();
		size_t mB = B.size();

		Eigen::VectorXd concatVec = Eigen::VectorXd(mA + mB);
		concatVec << A, B;

		return concatVec;
	}

	double BFFSolver::Util_TargetBoundaryLengths(Eigen::VectorXd& tBoundaryLength) const
	{
		tBoundaryLength = Eigen::VectorXd::Zero(m_NumBV);
		for (size_t bVertID = 0; bVertID < m_NumBV; ++bVertID)
		{
			size_t bID0 = Util_VertIndex(m_BoundaryVertices[bVertID]) - m_NumIV;
			size_t bID1 = Util_VertIndex(m_BoundaryVertices[(bVertID + 1) % m_NumBV]) - m_NumIV;

			double oBoundaryLength = m_BoundaryEdgeLengths[bID0];

			double u0 = m_BoundaryScaleFactor(bID0);
			double u1 = m_BoundaryScaleFactor(bID1);

			tBoundaryLength(bID0) = std::exp(0.5 * (u0 + u1)) * oBoundaryLength;
		}
		return tBoundaryLength.sum();
	}

	double BFFSolver::Util_TargetDualBoundaryLengths(Eigen::VectorXd& tBoundaryLength, Eigen::VectorXd& tDualLength) const
	{
		tDualLength = Eigen::VectorXd::Zero(m_NumBV);
		for (size_t bVertID = 0; bVertID < m_NumBV; ++bVertID)
		{
			size_t bID0 = Util_VertIndex(m_BoundaryVertices[bVertID]) - m_NumIV;
			size_t bID1 = Util_VertIndex(m_BoundaryVertices[(bVertID + 1) % m_NumBV]) - m_NumIV;

			tDualLength(bID1) = 0.5 * (tBoundaryLength(bID0) + tBoundaryLength(bID1));
		}
		return tDualLength.sum();
	}

	std::tuple<Eigen::VectorXd, Eigen::VectorXd> BFFSolver::ExtendCurve(const Eigen::MatrixXd& gamma_tilde, bool freeBoundary)
	{
		Eigen::VectorXd Re_gamma = gamma_tilde.row(0).transpose();
		Eigen::VectorXd rightB = -m_LapIB * Re_gamma;

		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldltSolver;
		ldltSolver.compute(m_LapII);
		Eigen::VectorXd aI = ldltSolver.solve(rightB);
		Eigen::VectorXd a = Util_VecVertConcat(aI, Re_gamma);
		Eigen::VectorXd b = Eigen::VectorXd::Zero(m_NumIV);

		if (freeBoundary)
		{
			Eigen::VectorXd h = Eigen::VectorXd::Zero(m_NumV);

			for (size_t bVertID = 0; bVertID < m_NumBV; ++bVertID)
			{
				size_t beforeID = (bVertID + m_NumBV - 1) % m_NumBV;
				size_t afterID = (bVertID + 1) % m_NumBV;

				auto* pBVert = m_BoundaryVertices[bVertID];
				auto* pBeforeBVert = m_BoundaryVertices[beforeID];
				auto* pAfterBVert = m_BoundaryVertices[afterID];

				h(Util_VertIndex(pBVert)) = 0.5 * (a(Util_VertIndex(pBeforeBVert)) - a(Util_VertIndex(pAfterBVert)));
			}

			ldltSolver.compute(m_LapMat);
			b = ldltSolver.solve(h);
		}
		else
		{
			Eigen::VectorXd Im_gamma = gamma_tilde.row(1).transpose();
			rightB = -m_LapIB * Im_gamma;
			Eigen::VectorXd bI = ldltSolver.solve(rightB);
			b = Util_VecVertConcat(bI, Im_gamma);
		}

		return std::make_tuple(a, b);
	}

	// constrain UV in range [0, 1]x[0, 1]
	void BFFSolver::Util_ConstrainUV()
	{
		double xmin = std::numeric_limits<double>::max();
		double xmax = -std::numeric_limits<double>::max();
		double ymin = std::numeric_limits<double>::max();
		double ymax = -std::numeric_limits<double>::max();

		for (size_t vertID = 0; vertID < m_NumV; ++vertID)
		{
			if (!m_Mesh->isBoundary(m_Mesh->vert(vertID)))	continue;
			double x = m_UVList(vertID, 0);
			double y = m_UVList(vertID, 1);

			xmin = std::min(xmin, x);
			xmax = std::max(xmax, x);
			ymin = std::min(ymin, y);
			ymax = std::max(ymax, y);
		}
		
		double midx = (xmin + xmax) * 0.5;
		double midy = (ymin + ymax) * 0.5;

		double maxScale = std::max(xmax - xmin, ymax - ymin);

		for (size_t vertID = 0; vertID < m_NumV; ++vertID)
		{
			auto* vert = m_Mesh->vert(vertID);
			double UVx = (m_UVList(vertID, 0) - midx) / maxScale;
			double UVy = (m_UVList(vertID, 1) - midy) / maxScale;

			vert->setTexture(UVy * 2.0f, -UVx * 2.0f, 0.0);
		}
	}
}
