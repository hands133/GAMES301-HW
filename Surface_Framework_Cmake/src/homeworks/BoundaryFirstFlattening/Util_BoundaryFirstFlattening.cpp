#include "Util_BoundaryFirstFlattening.h"


namespace bff
{
	void BFFSolver::Solve(acamcad::polymesh::PolyMesh* mesh)
	{
		// step 0. prepare : u, boundary lengths l
		//					 scale factor u or exterior angles k_tilde
		m_BoundaryVertices = mesh->boundaryVertices();
		auto bEdgeLengths = CalculateBoundaryLengths(m_BoundaryVertices);

		Util_RearrangeVertIndex(mesh);

		// step 1. discrete curvatures
		auto& [curvK, geodK] = CalculateDiscreteCurvatures(mesh);

		// step 2. build Laplacian matrix
		Eigen::SparseMatrix<double> A = BuildLaplacian(mesh);

		size_t numV = mesh->vertices().size();
		size_t numBV = mesh->boundaryVertices().size();
		size_t numIV = numV - numBV;

		Eigen::SparseMatrix<double> AII = A.block(0, 0, numIV, numIV);
		Eigen::SparseMatrix<double> AIB = A.block(0, numIV, numIV, numBV);
		Eigen::SparseMatrix<double> ABB = A.block(numIV, numIV, numBV, numBV);

	}

	Eigen::VectorXd BFFSolver::CalculateBoundaryLengths(const std::vector<acamcad::polymesh::MVert*>& bVertices) const
	{
		size_t numBV = bVertices.size();
		Eigen::VectorXd edges = Eigen::VectorXd(numBV);
		edges.setZero();

		for (size_t bVertID = 0; bVertID < numBV; ++bVertID)
			edges[bVertID] = acamcad::distance(bVertices[(bVertID + 1) % numBV]->position(),
											   bVertices[bVertID]->position());
		return edges;
	}

	std::tuple<Eigen::VectorXd, Eigen::VectorXd> BFFSolver::CalculateDiscreteCurvatures(acamcad::polymesh::PolyMesh* mesh) const
	{
		// Integrated gaussian curvature vs integrated geodesic curvature
		size_t numV = mesh->vertices().size();
		size_t numBV = m_BoundaryVertices.size();
		size_t numIV = numV - numBV;
		Eigen::VectorXd curvK = Eigen::VectorXd(numIV);
		Eigen::VectorXd geodk = Eigen::VectorXd(numBV);
		curvK.setZero();
		geodk.setZero();

		// Integrated gaussian curvature for interior vertices
		size_t innerVertCount = 0;
		for (auto* vert : mesh->vertices())
			if (!mesh->isBoundary(vert))
				curvK[innerVertCount++] = Util_VertAngleDefect(vert, mesh);

		// Integrated geodesic curvature for boundary vertices
		size_t boundVertCount = 0;
		for (auto* vert : m_BoundaryVertices)
			geodk[boundVertCount++] = Util_VertExteriorAngle(vert, mesh);

		return std::make_tuple(curvK, geodk);
	}

	Eigen::SparseMatrix<double> BFFSolver::BuildLaplacian(acamcad::polymesh::PolyMesh* mesh)
	{
		size_t numV = mesh->vertices().size();
		size_t numF = mesh->polyfaces().size();

		Eigen::SparseMatrix<double> A;
		std::vector<Eigen::Triplet<double>> ATriplets;

		for (auto* pFace : mesh->polyfaces())
		{
			acamcad::polymesh::MHalfedge* he = pFace->halfEdge();
			do
			{
				size_t i = m_VertIndexMap[he->fromVertex()];
				size_t j = m_VertIndexMap[he->toVertex()];
				double w = 0.5 * Util_HalfEdgeAngleCotan(he, mesh);

				ATriplets.emplace_back(i, i, w);
				ATriplets.emplace_back(j, j, w);
				ATriplets.emplace_back(i, j, -w);
				ATriplets.emplace_back(j, i, -w);

				he = he->next();
			} while (he != pFace->halfEdge());
		}

		A = Eigen::SparseMatrix<double>(numV, numV);
		A.setFromTriplets(ATriplets.begin(), ATriplets.end());

		auto AShift = Eigen::SparseMatrix<double>(numV, numV);
		AShift.setIdentity();
		A += AShift * std::numeric_limits<float>::epsilon();

		return A;
	}

	void BFFSolver::Util_RearrangeVertIndex(acamcad::polymesh::PolyMesh* mesh)
	{
		m_VertIndexMap.clear();
		m_VertIndexMap.reserve(mesh->vertices().size());

		size_t vertCount = 0;

		for (auto* pVert : mesh->vertices())
			if (!mesh->isBoundary(pVert))
				m_VertIndexMap.insert({ pVert, vertCount++ });

		for (auto* pBVert : mesh->boundaryVertices())
			m_VertIndexMap.insert({ pBVert, vertCount++ });
	}

	// Integrated gaussian curvature for inner vertices
	double BFFSolver::Util_VertAngleDefect(acamcad::polymesh::MVert* vert, acamcad::polymesh::PolyMesh* mesh) const
	{
		auto& neighbors = mesh->vertAdjacentVertices(vert);
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

	double BFFSolver::Util_VertExteriorAngle(acamcad::polymesh::MVert* vert, acamcad::polymesh::PolyMesh* mesh) const
	{
		assert(mesh->isBoundary(vert));

		auto& neighbors = mesh->vertAdjacentVertices(vert);
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
		return M_PI - angleSum;
	}

	double BFFSolver::Util_HalfEdgeAngleCotan(acamcad::polymesh::MHalfedge* he, acamcad::polymesh::PolyMesh* mesh) const
	{
		if (mesh->isBoundary(he->edge()))	return 0.0;

		auto* pVertA = he->fromVertex();
		auto* pVertB = he->next()->fromVertex();
		auto* pVertC = he->prev()->toVertex();

		auto vecU = pVertA->position() - pVertC->position();
		auto vecV = pVertB->position() - pVertC->position();

		double w = acamcad::dot(vecU, vecV) / acamcad::cross(vecU, vecV).norm();
		if (std::isinf(w) || std::isnan(w))	w = 0.0;

		return w;
	}
}
