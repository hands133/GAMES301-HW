#include "Util_FreeBoundary.h"

#include <algorithm>
#include <numeric>

#include <Eigen\Sparse>
//#include <Eigen\SparseQR>
#include <Eigen\src\QR\FullPivHouseholderQR.h>

#include <Eigen\OrderingMethods>


namespace freeb
{
	// CutGraph Stage 1. calculate the dual spanning tree of the original graph
	void CalculateDualSpanningTree(acamcad::polymesh::PolyMesh* mesh, size_t chosenFaceIdx,
		std::vector<char>& edgeSharped, std::vector<char>& faceTouched)
	{
		size_t numE = mesh->edges().size();
		size_t numF = mesh->polyfaces().size();

		std::queue<size_t> faceIdxQueue;
		faceIdxQueue.push(chosenFaceIdx);
		faceTouched[chosenFaceIdx] = 0x01;

		while (!faceIdxQueue.empty())
		{
			size_t faceIdx = faceIdxQueue.front();
			auto* pFace = mesh->polyface(faceIdx);
			faceIdxQueue.pop();

			for (auto* pHalfEdge : mesh->polygonHalfedges(pFace))
			{
				auto* pEdge = pHalfEdge->edge();
				size_t edgeIdx = pEdge->index();

				if (!mesh->isBoundary(pEdge))
				{
					auto* pTwinFace = pHalfEdge->pair()->polygon();
					size_t twinFaceIdx = pTwinFace->index();
					if (faceTouched[twinFaceIdx] == 0x00)
					{
						faceTouched[twinFaceIdx] = 0x01;
						edgeSharped[edgeIdx] = 0x01;

						faceIdxQueue.push(twinFaceIdx);
					}
				}
			}
		}
	}

	void PruneSpanningTreeBranches(acamcad::polymesh::PolyMesh* mesh, std::vector<char>& edgeSharped)
	{
		size_t numV = mesh->vertices().size();
		size_t numE = mesh->edges().size();

		std::vector<uint16_t> vertValences(numV, 0);

		std::vector<size_t> edgesTobeRemoved;
		edgesTobeRemoved.reserve(std::accumulate(edgeSharped.begin(), edgeSharped.end(), size_t{ 0 },
			[](size_t sum, char sharped) { return sum + static_cast<size_t>(sharped); }));

		while (true)
		{
			bool prune = false;

			for (size_t vertID = 0; vertID < numV; ++vertID)
			{
				acamcad::polymesh::MEdge* removeEdge = nullptr;
				vertValences[vertID] = 0;

				for (auto* pvAdjEdge : mesh->vertAdjacentEdge(mesh->vert(vertID)))
					if (edgeSharped[pvAdjEdge->index()] == 0x01)
					{
						vertValences[vertID]++;
						removeEdge = pvAdjEdge;
					}

				if (vertValences[vertID] == 1)
				{
					edgesTobeRemoved.emplace_back(removeEdge->index());
					prune = true;
				}
			}

			std::for_each(edgesTobeRemoved.begin(), edgesTobeRemoved.end(),
				[&](size_t edgeID)
				{
					edgeSharped[edgeID] = 0x00;
			vertValences[mesh->edge(edgeID)->getVert(0)->index()]--;
			vertValences[mesh->edge(edgeID)->getVert(1)->index()]--;
				});

			edgesTobeRemoved.clear();
			if (!prune)	return;
		}
	}

	std::vector<char> CutEdge(acamcad::polymesh::PolyMesh* mesh)
	{
		size_t numE = mesh->edges().size();
		size_t numF = mesh->polyfaces().size();

		std::vector<char> edgeSharped(numE, 0x00);
		std::vector<char> faceTouched(numF, 0x00);

		CalculateDualSpanningTree(mesh, 0, edgeSharped, faceTouched);

		for (size_t edgeID = 0; edgeID < numE; ++edgeID)
			edgeSharped[edgeID] = 1 - edgeSharped[edgeID];

		// Stage 2. prune vertices which valences == 1 till terminated
		PruneSpanningTreeBranches(mesh, edgeSharped);

		return edgeSharped;
	}
}

namespace freeb
{
	// constrain UV in range [0, 1]x[0, 1]
	void FreeBoundarySolver::ConstrainUV(acamcad::polymesh::PolyMesh* mesh) const
	{
		double xmin = std::numeric_limits<double>::max();
		double xmax = -std::numeric_limits<double>::max();
		double ymin = std::numeric_limits<double>::max();
		double ymax = -std::numeric_limits<double>::max();

		size_t numV = mesh->vertices().size();
		for (size_t vertID = 0; vertID < numV; ++vertID)
		{
			double x = m_UVList(vertID * 2 + 0);
			double y = m_UVList(vertID * 2 + 1);

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
			float UVx = static_cast<float>((m_UVList(2 * vertID + 0) - midx) / (xmax - xmin));
			float UVy = static_cast<float>((m_UVList(2 * vertID + 1) - midy) / (ymax - ymin));
			vert->setTexture(UVx * 2.0, UVy * 2.0, 0.0);
		}
	}

	// Check if the planar m_UVList has no flip-over triangles
	bool FreeBoundarySolver::CheckPlanar(acamcad::polymesh::PolyMesh* mesh,
		const std::vector<acamcad::polymesh::MVert*>& boundaries) const
	{
		double angle = 0.0;

		size_t nB = boundaries.size();
		for (size_t b = 0; b < nB; ++b)
		{
			size_t n = (b + 1) % nB;
			size_t p = (b + nB - 1) % nB;

			size_t vpID = boundaries[p]->index();
			size_t v0ID = boundaries[b]->index();
			size_t vnID = boundaries[n]->index();

			acamcad::MPoint3 vp(m_UVList[vpID * 2], m_UVList[vpID * 2 + 1], 0.0);
			acamcad::MPoint3 v0(m_UVList[v0ID * 2], m_UVList[v0ID * 2 + 1], 0.0);
			acamcad::MPoint3 vn(m_UVList[vnID * 2], m_UVList[vnID * 2 + 1], 0.0);

			angle += acamcad::signedAngle(v0 - vp, vn - v0, { 0, 0, 1 });
		}

		return std::abs(std::abs(angle) - 2.0 * M_PI) < std::numeric_limits<float>::epsilon();
	}

	std::vector<double> FreeBoundarySolver::boundaryEdgesLength(const std::vector<acamcad::polymesh::MVert*>& boundaries) const
	{
		size_t nB = boundaries.size();
		std::vector<double> lengths(nB, 0.0);

		for (size_t b = 0; b < nB; ++b)
		{
			size_t n = (b + 1) % nB;
			size_t viID = boundaries[b]->index();
			size_t vjID = boundaries[n]->index();

			Eigen::Vector2d vi;
			vi << m_UVList(viID * 2), m_UVList(viID * 2 + 1);
			Eigen::Vector2d vj;
			vj << m_UVList(vjID * 2), m_UVList(vjID * 2 + 1);

			lengths[b] = (vi - vj).norm();
		}

		return lengths;
	}

	void FreeBoundarySolver::PresetMeshUV(acamcad::polymesh::PolyMesh* mesh)
	{
		size_t numV = mesh->vertices().size();
		size_t numE = mesh->polyfaces().size();

		std::vector<Eigen::Triplet<double>> coeffTriplets;
		coeffTriplets.reserve(6 * numV);

		m_Coeffs = Eigen::SparseMatrix<double>(2 * numV + 4, 2 * numV);

		// step 1. inner nodes
		for (auto* pVert : mesh->vertices())
		{
			auto vNeighbors = mesh->vertAdjacentVertices(pVert);
			auto nWeights = ConformWeight(pVert, vNeighbors);

			for (size_t n = 0; n < vNeighbors.size(); ++n)
			{
				auto* nVert = vNeighbors[n];
				coeffTriplets.emplace_back(pVert->index() * 2, nVert->index() * 2, nWeights[n]);
				coeffTriplets.emplace_back(pVert->index() * 2 + 1, nVert->index() * 2 + 1, nWeights[n]);
			}

			double sumWeights = std::accumulate(nWeights.begin(), nWeights.end(), 0.0);
			coeffTriplets.emplace_back(pVert->index() * 2, pVert->index() * 2, -sumWeights);
			coeffTriplets.emplace_back(pVert->index() * 2 + 1, pVert->index() * 2 + 1, -sumWeights);
		}

		// step 2. boundary nodes
		auto boundaries = mesh->boundaryVertices();
		size_t nB = boundaries.size();
		for (size_t b = 0; b < nB; ++b)
		{
			size_t b1 = (b + 1) % nB;
			size_t b2 = (b + nB - 1) % nB;

			auto* pv1 = boundaries[b1];
			auto* pv2 = boundaries[b2];

			size_t vID = boundaries[b]->index();
			// equation for x
			coeffTriplets.emplace_back(vID * 2, pv1->index() * 2 + 1, -1.0);
			coeffTriplets.emplace_back(vID * 2, pv2->index() * 2 + 1,  1.0);
			// equation for y
			coeffTriplets.emplace_back(vID * 2 + 1, pv1->index() * 2,  1.0);
			coeffTriplets.emplace_back(vID * 2 + 1, pv2->index() * 2, -1.0);
		}

		// step 3. fix two points
		Eigen::SparseVector<double> b = Eigen::SparseVector<double>(2 * numV + 4);
		b.setZero();

		size_t tv0ID = 0;			// fixed at (0,  0.5)
		size_t tv1ID = numV - 1;	// fixed at (0, -0.5)
		//coeffTriplets.emplace_back(tv0ID * 2 + 0, tv0ID * 2 + 0, 1.0);
		//coeffTriplets.emplace_back(tv0ID * 2 + 1, tv0ID * 2 + 1, 1.0);
		//coeffTriplets.emplace_back(tv1ID * 2 + 0, tv1ID * 2 + 0, 1.0);
		//coeffTriplets.emplace_back(tv1ID * 2 + 1, tv1ID * 2 + 1, 1.0);

		//b.insert(tv0ID * 2 + 0) =  0.0;
		//b.insert(tv0ID * 2 + 1) =  0.5;
		//b.insert(tv1ID * 2 + 0) =  0.0;
		//b.insert(tv1ID * 2 + 1) = -0.5;

		coeffTriplets.emplace_back(numV * 2 + 0, tv0ID * 2 + 0, 1.0);
		coeffTriplets.emplace_back(numV * 2 + 1, tv0ID * 2 + 1, 1.0);
		coeffTriplets.emplace_back(numV * 2 + 2, tv1ID * 2 + 0, 1.0);
		coeffTriplets.emplace_back(numV * 2 + 3, tv1ID * 2 + 1, 1.0);

		b.insert(numV * 2 + 0) =  0.0;
		b.insert(numV * 2 + 1) =  0.5;
		b.insert(numV * 2 + 2) =  0.0;
		b.insert(numV * 2 + 3) = -0.5;


		m_Coeffs.setFromTriplets(coeffTriplets.begin(), coeffTriplets.end());
		m_Coeffs.makeCompressed();

		Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::NaturalOrdering<int>> solver;
		solver.compute(m_Coeffs);
		m_UVList = solver.solve(b);
		
		std::cout << "Iter = " << ++m_Iters << "\t||Ax-b||=" << (m_Coeffs * m_UVList - b).norm() << "\n";
		//std::cout << m_Coeffs << "\n\n" << b.transpose() << "\n" << m_UVList.transpose() << "\n\n";

		ConstrainUV(mesh);
	}

	void FreeBoundarySolver::UpdateMeshUV(acamcad::polymesh::PolyMesh* mesh)
	{
		size_t numV = mesh->vertices().size();
		size_t numE = mesh->polyfaces().size();

		auto boundaries = mesh->boundaryVertices();

		Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;

		std::vector<Eigen::Triplet<double>> coeffTriplets;
		coeffTriplets.reserve(6 * numV);

		while (!CheckPlanar(mesh, boundaries))
		{
			coeffTriplets.clear();

			m_Coeffs.setZero();
			m_Coeffs = Eigen::SparseMatrix<double>(2 * numV + 4, 2 * numV);

			// step 1. inner nodes
			for (auto* pVert : mesh->vertices())
			{
				auto vNeighbors = mesh->vertAdjacentVertices(pVert);
				auto nWeights = MeanValueWeight(pVert, vNeighbors);

				for (size_t n = 0; n < vNeighbors.size(); ++n)
				{
					auto* nVert = vNeighbors[n];
					coeffTriplets.emplace_back(pVert->index() * 2, nVert->index() * 2, nWeights[n]);
					coeffTriplets.emplace_back(pVert->index() * 2 + 1, nVert->index() * 2 + 1, nWeights[n]);
				}

				double sumWeights = std::accumulate(nWeights.begin(), nWeights.end(), 0.0);
				coeffTriplets.emplace_back(pVert->index() * 2, pVert->index() * 2, -sumWeights);
				coeffTriplets.emplace_back(pVert->index() * 2 + 1, pVert->index() * 2 + 1, -sumWeights);
			}

			// step 2. boundary nodes
			auto boundaries = mesh->boundaryVertices();
			auto boundaryLengths = boundaryEdgesLength(boundaries);

			size_t nB = boundaries.size();
			for (size_t b = 0; b < nB; ++b)
			{
				size_t b1 = (b + 1) % nB;
				size_t b2 = (b + nB - 1) % nB;

				auto* pv1 = boundaries[b1];
				auto* pv2 = boundaries[b2];

				double inv_r1 = 1.0 / boundaryLengths[b1];
				double inv_r2 = 1.0 / boundaryLengths[b2];

				size_t vID = boundaries[b]->index();
				// equation for x
				coeffTriplets.emplace_back(vID * 2, vID * 2 + 1, inv_r1 - inv_r2);
				coeffTriplets.emplace_back(vID * 2, pv1->index() * 2 + 1, -inv_r1);
				coeffTriplets.emplace_back(vID * 2, pv2->index() * 2 + 1, inv_r2);
				// equation for y
				coeffTriplets.emplace_back(vID * 2 + 1, vID * 2, inv_r2 - inv_r1);
				coeffTriplets.emplace_back(vID * 2 + 1, pv1->index() * 2, inv_r1);
				coeffTriplets.emplace_back(vID * 2 + 1, pv2->index() * 2, -inv_r2);
			}

			// step 3. fix two points
			Eigen::SparseVector<double> b = Eigen::SparseVector<double>(2 * numV + 4, 1);

			size_t tv0ID = 0;			// fixed at (0,  0.5)
			size_t tv1ID = numV - 1;	// fixed at (0, -0.5)
			coeffTriplets.emplace_back(numV * 2 + 0, tv0ID * 2 + 0, 1.0);
			coeffTriplets.emplace_back(numV * 2 + 1, tv0ID * 2 + 1, 1.0);
			coeffTriplets.emplace_back(numV * 2 + 2, tv1ID * 2 + 0, 1.0);
			coeffTriplets.emplace_back(numV * 2 + 3, tv1ID * 2 + 1, 1.0);

			//b.insert(tv0ID * 2 + 0) =  0.0;
			b.insert(tv0ID * 2 + 1) = 0.5;
			//b.insert(tv1ID * 2 + 0) =  0.0;
			b.insert(tv1ID * 2 + 1) = -0.5;

			m_Coeffs.setFromTriplets(coeffTriplets.begin(), coeffTriplets.end());
			m_Coeffs.makeCompressed();

			solver.compute(m_Coeffs);
			m_UVList = solver.solve(b);

			//std::cout << "Iter = " << ++m_Iters << "\t" << m_UVList.transpose() << std::endl;
			std::cout << "Iter = " << ++m_Iters << "\t||Ax-b||=" << (m_Coeffs * m_UVList - b).norm() << "\n";
		}

		ConstrainUV(mesh);
	}

	std::vector<double> FreeBoundarySolver::ConformWeight(acamcad::polymesh::MVert* v,
		const std::vector<acamcad::polymesh::MVert*>& adjVerts) const
	{
		size_t nN = adjVerts.size();	// number of neighbors

		std::vector<double> nWeights(nN, 0.0);
		std::vector<Eigen::Vector2d> cosAngles(nN, Eigen::Vector2d::Zero());

		for (size_t i = 0; i < nN; ++i)
		{
			size_t j = (i + 1) % nN;

			auto* pVi = adjVerts[i];
			auto* pVj = adjVerts[j];

			double a2 = (pVi->position() - pVj->position()).normSq();
			double b2 = (pVi->position() - v->position()).normSq();
			double c2 = (pVj->position() - v->position()).normSq();

			double cosB = (a2 + c2 - b2) / (2.0 * std::sqrt(a2) * std::sqrt(c2));
			double cosC = (a2 + b2 - c2) / (2.0 * std::sqrt(a2) * std::sqrt(b2));

			cosAngles[i] << cosB, cosC;
		}

		for (size_t i = 0; i < nN; ++i)
		{
			size_t p = (i + nN - 1) % nN;

			double cosBi = cosAngles[i].x();
			double cosCp = cosAngles[p].y();

			double sinBi = std::sqrt(1.0 - cosBi * cosBi);
			double sinCp = std::sqrt(1.0 - cosCp * cosCp);

			double cotBi = cosBi / sinBi;
			double cotCp = cosCp / sinCp;

			nWeights[i] = cotBi + cotCp;
		}

		return nWeights;
	}

	std::vector<double> FreeBoundarySolver::MeanValueWeight(acamcad::polymesh::MVert* v, const std::vector<acamcad::polymesh::MVert*>& adjVerts) const
	{
		size_t nN = adjVerts.size();	// number of neighbors

		std::vector<double> nWeights(nN, 0.0);
		std::vector<Eigen::Vector2d> cosAnglesNDistanceSquare(nN, Eigen::Vector2d::Zero());

		size_t v0ID = v->index();

		for (size_t i = 0; i < nN; ++i)
		{
			size_t j = (i + 1) % nN;

			size_t viID = adjVerts[i]->index();
			size_t vjID = adjVerts[j]->index();

			Eigen::Vector2d vecA = m_UVList.segment(2 * viID, 2) - m_UVList.segment(2 * vjID, 2);
			Eigen::Vector2d vecB = m_UVList.segment(2 * viID, 2) - m_UVList.segment(2 * v0ID, 2);
			Eigen::Vector2d vecC = m_UVList.segment(2 * vjID, 2) - m_UVList.segment(2 * v0ID, 2);

			double a2 = vecA.dot(vecA);
			double b2 = vecA.dot(vecA);
			double c2 = vecB.dot(vecB);

			double cosA = (b2 + c2 - a2) / (2.0 * std::sqrt(b2) * std::sqrt(c2));
			cosAnglesNDistanceSquare[i] << cosA, b2;
		}

		for (size_t i = 0; i < nN; ++i)
		{
			size_t p = (i + nN - 1) % nN;

			double ri = cosAnglesNDistanceSquare[i].y();

			double cosAi = cosAnglesNDistanceSquare[i].x();
			double cosAp = cosAnglesNDistanceSquare[p].x();

			double sinAi = std::sqrt(1.0 - cosAi * cosAi);
			double sinAp = std::sqrt(1.0 - cosAp * cosAp);

			double tanAi_2 = sinAi / (1.0 + cosAi);
			double tanAp_2 = sinAp / (1.0 + cosAp);

			nWeights[i] = (tanAi_2 + tanAp_2) / ri;
		}
		
		return nWeights;
	}

}