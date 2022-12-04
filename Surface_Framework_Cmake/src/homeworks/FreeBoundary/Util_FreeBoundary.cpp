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
		
		std::cout << "[Free Boundary] pinned index = [" << pID0 << ", " << pID1 << "]\n";

		Eigen::Vector2d fixP0;
		fixP0 << -0.5, 0.0;
		Eigen::Vector2d fixP1;
		fixP1 << 0.5, 0.0;
		m_UVList.segment(pID0 * 2, 2) = fixP0;
		m_UVList.segment(pID1 * 2, 2) = fixP1;

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

		std::cout << "[Free Boundary] LSCM Start\n";
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

		{	// step 3. calculate matrix Mf and Mp
			Eigen::Vector4d vecU;
			vecU << fixP0.x(), fixP1.x(), fixP0.y(), fixP1.y();

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
		}

		auto timeEnd = std::chrono::steady_clock::now();
		auto ms = std::chrono::duration_cast<std::chrono::microseconds>(timeEnd - timeStart).count() / 1000.0;

		std::cout << "[Free Boundary] LSCM End in " << ms << " ms\n";
	}

	void FreeBoundarySolver::TutteWithMeanValue(acamcad::polymesh::PolyMesh* mesh)
	{
		using acamcad::polymesh::MVert;

		auto timeStart = std::chrono::steady_clock::now();

		auto boundaryVertLists = mesh->boundaryVertices();
		std::cout << "[Free Boundary] Boundary Points [" << boundaryVertLists.size() << "]\n";

		// 1. prepare convex polygon
		size_t numB = boundaryVertLists.size();
		size_t numV = mesh->vertices().size();

		// 2. prepare matrix (Eigen::Sparse)
		Eigen::SparseMatrix<double> A(numV, numV);	A.setZero();
		Eigen::SparseMatrix<double> x(numV, 2);		x.setZero();
		Eigen::SparseMatrix<double> b(numV, 2);		b.setZero();

		std::vector<Eigen::Triplet<double>> ATriplets;
		std::vector<Eigen::Triplet<double>> bTriplets;

		// a) boundary vertices
		for (auto* pBVert : boundaryVertLists)
		{
			size_t vertID = pBVert->index();
			ATriplets.emplace_back(vertID, vertID, 1.0);
			bTriplets.emplace_back(vertID, 0, m_UVList(vertID * 2 + 0));
			bTriplets.emplace_back(vertID, 1, m_UVList(vertID * 2 + 1));
		}

		// b) inner vertices
		for (auto* pVert : mesh->vertices())
		{
			if (mesh->isBoundary(pVert))	continue;

			size_t vertID = pVert->index();

			auto adjVerts = mesh->vertAdjacentVertices(pVert);
			auto adjWeights = MeanValueWeight(pVert, adjVerts);
			double weightSUM = std::accumulate(adjWeights.begin(), adjWeights.end(), 0.0);

			for (int j = 0; j < adjVerts.size(); ++j)
				ATriplets.emplace_back(vertID, adjVerts[j]->index(), adjWeights[j]);
			ATriplets.emplace_back(vertID, vertID, -weightSUM);
		}

		A.setFromTriplets(ATriplets.begin(), ATriplets.end());
		b.setFromTriplets(bTriplets.begin(), bTriplets.end());

		// 3. solve the equation Ax = b
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
		solver.compute(A);
		x = solver.solve(b);

		auto timeEnd = std::chrono::steady_clock::now();
		auto ms = std::chrono::duration_cast<std::chrono::microseconds>(timeEnd - timeStart).count() / 1000.0;

		std::cout << "[Free Boundary] Tutte's parameterization end in " << ms << "ms\n";
		
		for (auto* pVert : mesh->vertices())
		{
			if (mesh->isBoundary(pVert)) continue;
			size_t vertID = pVert->index();
			m_UVList.segment(vertID * 2, 2) = x.row(vertID);
		}
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
			double b2 = vecB.dot(vecB);
			double c2 = vecC.dot(vecC);

			double cosA = (b2 + c2 - a2) / (2.0 * std::sqrt(b2) * std::sqrt(c2));
			cosAnglesNDistanceSquare[i] << cosA, b2;
		}

		for (size_t i = 0; i < nN; ++i)
		{
			size_t p = (i + nN - 1) % nN;

			double ri = std::sqrt(cosAnglesNDistanceSquare[i].y());

			double cosAi = cosAnglesNDistanceSquare[i].x();
			double cosAp = cosAnglesNDistanceSquare[p].x();
			double sinAi = std::sqrt(1.0 - cosAi * cosAi);
			double sinAp = std::sqrt(1.0 - cosAp * cosAp);
			double tanAi_2 = sinAi / (1.0 + cosAi);
			double tanAp_2 = sinAp / (1.0 + cosAp);

			nWeights[i] = (tanAi_2 + tanAp_2) / ri;

			assert(nWeights[i] > 0.0, "negative!");
		}
		
		return nWeights;
	}

}