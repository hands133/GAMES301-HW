#include "Util_FreeBoundary.h"

#include <algorithm>
#include <iostream>
#include <numeric>
#include <chrono>

#include <Eigen\Dense>
#include <Eigen\Sparse>
#include <Eigen\SparseQR>
#include <Eigen\IterativeLinearSolvers>

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
	void FreeBoundarySolver::Solve(acamcad::polymesh::PolyMesh* mesh)
	{
		m_UVList = Eigen::VectorXd(mesh->vertices().size() * 2);
		m_UVList.setZero();

		// step 1. get coarse boundary with LSCM
		LSCMUV(mesh);

		// step 2. calculate tutte's parameterization with mean-value coordinate
		TutteWithMeanValue(mesh);

		// step 3. write back
		ConstrainUV(mesh);
	}

	void FreeBoundarySolver::LSCMUV(acamcad::polymesh::PolyMesh* mesh)
	{
		size_t numV = mesh->vertices().size();
		size_t numF = mesh->polyfaces().size();

		// boundaries
		auto meshBoundaries = mesh->boundaryVertices();
		size_t pID0 = meshBoundaries[0]->index();
		size_t pID1 = meshBoundaries[meshBoundaries.size() / 2]->index();
		size_t numP = 2;
		std::cout << "pinned index = [" << pID0 << ", " << pID1 << "]\n";

		std::vector<acamcad::polymesh::MVert*> pinnedVert(numP, nullptr);

		// A = [ Re(Mf) -Im(Mf) ]
		//     [ Im(Mf)  Re(Mf) ](2F * 2(V-P))
		//		Mf(F * (V-P))
		Eigen::SparseMatrix<double> matA = Eigen::SparseMatrix<double>(2 * numF, 2 * (numV - numP));
		// B = [ Re(Mp) -Im(Mp) ]
		//     [ Im(Mp)  Re(Mp) ](2F * 2P)
		//		Mp(F * P)
		Eigen::SparseMatrix<double> matB = Eigen::SparseMatrix<double>(2 * numF, 2 * numP);

		std::vector<Eigen::Triplet<double>> matAtriples;
		std::vector<Eigen::Triplet<double>> matBtriples;

		std::cout << "Free Boundary with LSCM\n";
		auto timeStart = std::chrono::steady_clock::now();

		std::unordered_map<acamcad::polymesh::MVert*, size_t> freeVertIndexMap;
		{
			std::vector<acamcad::polymesh::MVert*> freeVert;
			freeVert.reserve(numV - numP);
			for (auto* pV : mesh->vertices())
			{
				if (pV->index() == pID0 || pV->index() == pID1)	 continue;
				freeVert.emplace_back(pV);
			}

			for (size_t i = 0; i < numV - numP; ++i)	freeVertIndexMap.insert({ freeVert[i], i });
		}

		// step 1. position on local orthonormal bases
		for (size_t faceID = 0; faceID < numF; ++faceID)
		{
			auto faceVerts = mesh->polygonVertices(mesh->polyface(faceID));
			auto v21 = faceVerts[1]->position() - faceVerts[0]->position();
			auto v31 = faceVerts[2]->position() - faceVerts[0]->position();
			// y
			// ¦«
			// |    v3
			// .<---.
			// ¦«   /|\
			// |  / | \
			// | /  |  \
			// |/   V   \
			// .===>.--->.-> x
			// v1   vp   v2
			Eigen::Vector2d P1{ 0, 0 };
			Eigen::Vector2d P2{ v21.norm(), 0 };
			Eigen::Vector2d P3{ v21.dot(v31) / v21.norm(), v21.cross(v31).norm() / v21.norm() };

			double faceArea2 = (P1.x() * P2.y() - P1.y() * P2.x()) +
							   (P2.x() * P3.y() - P2.y() * P3.x()) +
							   (P3.x() * P1.y() - P3.y() * P1.x());

			std::array<Eigen::Vector2d, 3> Ws = {	Eigen::Vector2d{ P3.x() - P2.x(), P3.y() - P2.y() } / std::sqrt(faceArea2),
													Eigen::Vector2d{ P1.x() - P3.x(), P1.y() - P3.y() } / std::sqrt(faceArea2),
													Eigen::Vector2d{ P2.x() - P1.x(), P2.y() - P1.y() } / std::sqrt(faceArea2) };

			// step 2. calculate matrix A and matrix B
			for (size_t vID = 0; vID < 3; ++vID)
			{
				size_t vertID = faceVerts[vID]->index();
				auto W = Ws[vID];

				if (vertID == pID0 || vertID == pID1)
				{
					size_t targetVertID = (vertID == pID1);
					matBtriples.emplace_back(faceID, targetVertID, W.x());
					matBtriples.emplace_back(faceID + numF, targetVertID + numP, W.x());
					matBtriples.emplace_back(faceID, targetVertID + numP, -W.y());
					matBtriples.emplace_back(faceID + numF, targetVertID, W.y());
				} else {
					size_t targetVertID = freeVertIndexMap[faceVerts[vID]];
					matAtriples.emplace_back(faceID, targetVertID, W.x());
					matAtriples.emplace_back(faceID + numF, targetVertID + numV - numP, W.x());
					matAtriples.emplace_back(faceID, targetVertID + numV - numP, -W.y());
					matAtriples.emplace_back(faceID + numF, targetVertID, W.y());
				}
			}
		}

		matA.setFromTriplets(matAtriples.begin(), matAtriples.end());
		matA.makeCompressed();
		matB.setFromTriplets(matBtriples.begin(), matBtriples.end());

		// step 3. calculate matrix Mf and Mp
		{
			Eigen::Vector4d vecU;
			vecU << 0.0, 0.0, 0.0, 1.0;

			// b = -B [ Re(Up) Im(Up) ]T
			Eigen::VectorXd vecb = -matB * vecU;
			
			Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver;
			solver.compute(matA);
			Eigen::VectorXd X = solver.solve(vecb);

			size_t bIndex = 0;
			for (size_t vertID = 0; vertID < numV; ++vertID)
			{
				if (vertID == pID0 || vertID == pID1)	continue;
				m_UVList(vertID * 2 + 0) = X(bIndex);
				m_UVList(vertID * 2 + 1) = X(bIndex + numV - numP);
				bIndex++;
			}
			m_UVList.segment(pID0 * 2, 2) = vecU.segment(0, 2);
			m_UVList.segment(pID1 * 2, 2) = vecU.segment(2, 2);
		}

		auto timeEnd = std::chrono::steady_clock::now();
		auto ms = std::chrono::duration_cast<std::chrono::microseconds>(timeEnd - timeStart).count() / 1000.0;

		std::cout << "Free Boundary with LSCM in " << ms << " ms\n";
	}

	void FreeBoundarySolver::TutteWithMeanValue(acamcad::polymesh::PolyMesh* mesh)
	{

	}

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

		double maxScale = std::max(xmax - xmin, ymax - ymin);

		for (size_t vertID = 0; vertID < numV; ++vertID)
		{
			auto* vert = mesh->vert(vertID);
			float UVx = static_cast<float>((m_UVList(2 * vertID + 0) - midx) / maxScale);
			float UVy = static_cast<float>((m_UVList(2 * vertID + 1) - midy) / maxScale);
			vert->setTexture(UVx * 2.0f, UVy * 2.0f, 0.0);
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