#include "Util_FreeBoundary.h"

#include <algorithm>
#include <numeric>

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