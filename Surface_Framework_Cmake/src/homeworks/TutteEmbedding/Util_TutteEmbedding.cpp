#include "Util_TutteEmbedding.h"

#include <algorithm>
#include <numeric>

#include <Eigen\Eigen>

namespace tutte
{
	std::vector<glm::dvec2> GetBoundaryUVs(size_t numVerts, UVBoundaryType type)
	{
		const double CIRCLE_RADIUS = 0.5;
		std::vector<glm::dvec2> boundaryUVs(numVerts, { CIRCLE_RADIUS, CIRCLE_RADIUS });

		// 1. polygon edges
		size_t polygonEdges = 3;
		switch (type)
		{
		case UVBoundaryType::POLYGON_CIRCLE:	polygonEdges = numVerts;	break;
		case UVBoundaryType::POLYGON_TRIANGLE:	polygonEdges = 3;			break;
		case UVBoundaryType::POLYGON_SQUARE:	polygonEdges = 4;			break;
		case UVBoundaryType::POLYGON_PENTAGON:	polygonEdges = 5;			break;
		}

		polygonEdges = std::min(polygonEdges, numVerts);

		// 2. initialization
		std::vector<int> polygonEdgePoints(polygonEdges, 0);
		{
			int V = numVerts;
			int i = 0;
			while (V-- > 0)
			{
				polygonEdgePoints[i]++;
				i = (i + 1) % polygonEdges;
			}
		}
		assert(std::accumulate(polygonEdgePoints.begin(), polygonEdgePoints.end(), 0) == numVerts);

		double angleGap = 2.0 * M_PI / static_cast<double>(polygonEdges);
		double edgeLength = 2.0 * CIRCLE_RADIUS * glm::sin(angleGap / 2.0);

		// rotate matrix
		glm::dmat2x2 polygonRotateMat = glm::dmat2({ 0, 1 }, { -1, 0 });

		// 3. calculation
		auto boundaryUVIter = boundaryUVs.begin();
		for (size_t edgeIdx = 0; edgeIdx < polygonEdges; ++edgeIdx)
		{
			glm::dvec2 basePoint = CIRCLE_RADIUS * glm::dvec2{ glm::cos(edgeIdx * angleGap), glm::sin(edgeIdx * angleGap) };
			glm::dvec2 edgeDirection = CIRCLE_RADIUS * glm::dvec2(glm::cos(((edgeIdx + 1) % polygonEdges) * angleGap),
				glm::sin(((edgeIdx + 1) % polygonEdges) * angleGap)) - basePoint;

			int edgePoints = polygonEdgePoints[edgeIdx];
			double edgeDisDelta = edgeLength / static_cast<double>(edgePoints);

			int i = 0;
			while (i++ < edgePoints) {
				glm::dvec2 pointUV = basePoint + (i / static_cast<double>(edgePoints)) * edgeDirection;
				*boundaryUVIter += (polygonRotateMat * pointUV);
				boundaryUVIter++;
			}
		}

		return boundaryUVs;
	}

	auto FloaterParam_I_a(acamcad::polymesh::MVert* v,
		const std::vector<acamcad::polymesh::MVert*>& adjVerts)
	{
		double thetaI{ 0.0 };
		int N = adjVerts.size();

		std::vector<glm::dvec2> disANDangles(N, { 0.0, 0.0 });

		for (int i = 0; i < N; ++i)
		{
			int j = (i + 1) % N;
			auto vecI = adjVerts[i]->position() - v->position();
			auto vecJ = adjVerts[j]->position() - v->position();

			glm::dvec3 vi(vecI.x(), vecI.y(), vecI.z());
			glm::dvec3 vj(vecJ.x(), vecJ.y(), vecJ.z());

			double di = glm::sqrt(glm::dot(vi, vi));
			double dj = glm::sqrt(glm::dot(vj, vj));

			disANDangles[i].x = di;
			disANDangles[i].y = glm::acos(glm::dot(vi, vj) / (di * dj));
			thetaI += disANDangles[i].y;
		}
		assert(thetaI > 0);
		return std::make_tuple(thetaI, disANDangles);
	}

	auto FloaterParam_I_b(double thetaSum, const std::vector<glm::dvec2>& neighborAngleInfo)
	{
		int N = neighborAngleInfo.size();

		std::vector<glm::dvec2> adjUVs(N, { 0.0, 0.0 });
		adjUVs[0] = { neighborAngleInfo[0].x, 0.0 };

		for (int j = 1; j < N; ++j)
		{
			// 1. rotate factor
			double thetaI = neighborAngleInfo[j - 1].y * (2.0 * M_PI / thetaSum);
			double cosThetaI = glm::cos(thetaI);
			double sinThetaI = glm::sin(thetaI);
			glm::dmat2x2 rotateMat2({ cosThetaI, sinThetaI }, { -sinThetaI, cosThetaI });

			// 2. scale factor
			double scale = neighborAngleInfo[j].x / neighborAngleInfo[j - 1].x;
			glm::dmat2x2 scaleMat2({ scale, 0.0 }, { 0.0, scale });

			adjUVs[j] = scaleMat2 * rotateMat2 * adjUVs[j - 1];
		}

		return adjUVs;
	}

	// If triangle P1P2P3 centering the original point (0, 0)
	bool IsOinTriangle(glm::dvec2 P1, glm::dvec2 P2, glm::dvec2 P3)
	{
		double t1 = P1.x * P2.y - P1.y * P2.x;
		double t2 = P2.x * P3.y - P2.y * P3.x;
		double t3 = P3.x * P1.y - P3.y * P1.x;

		return t1 * t2 >= 0 && t1 * t3 >= 0 && t2 * t3 >= 0;
	}

	auto FloaterParam_II(const std::vector<glm::dvec2>& adjUVs)
	{
		int N = adjUVs.size();
		std::vector<double> mu_ls(N, 0.0);

		for (int i = 0; i < N; ++i) {

			auto Pl = adjUVs.begin() + i;

			auto Pr = adjUVs.begin() + (i + 1) % N;
			auto Pn = adjUVs.begin() + (i + 2) % N;

			auto IterAscend = [&]() {
				Pr++;
				Pn++;
				if (Pr == adjUVs.end())	Pr = adjUVs.begin();
				if (Pn == adjUVs.end()) Pn = adjUVs.begin();
			};

			while (Pn != Pl)
			{
				if (!IsOinTriangle(*Pl, *Pr, *Pn))
				{
					IterAscend();
					continue;
				};

				Eigen::Matrix3d A;
				A(0, 0) = Pl->x;
				A(1, 0) = Pl->y;
				A(2, 0) = 1.0;
				A(0, 1) = Pr->x;
				A(1, 1) = Pr->y;
				A(2, 1) = 1.0;
				A(0, 2) = Pn->x;
				A(1, 2) = Pn->y;
				A(2, 2) = 1.0;
				Eigen::Vector3d b(0.0, 0.0, 1.0);

				auto x = A.lu().solve(b);

				mu_ls[Pl - adjUVs.begin()] += x(0);
				mu_ls[Pr - adjUVs.begin()] += x(1);
				mu_ls[Pn - adjUVs.begin()] += x(2);
				mu_ls[i] += x(0, 0);

				IterAscend();
			}
		}

		std::for_each(mu_ls.begin(), mu_ls.end(), [N](double& v) { return v / static_cast<double>(N); });
		return mu_ls;
	}

	void AverageParam(acamcad::polymesh::MVert* v,
		const std::vector<acamcad::polymesh::MVert*>& adjVerts,
		std::vector<double>& weights)
	{
		weights.clear();
		weights.resize(adjVerts.size(), 1.0);
	}

	void FloaterParam(acamcad::polymesh::MVert* v,
		const std::vector<acamcad::polymesh::MVert*>& adjVerts,
		std::vector<double>& weights)
	{
		int N = adjVerts.size();
		weights.clear();

		/// ====== Stage 1 : Initialize (u, v) list ======
		// a) thetaI
		auto [thetaI, disANDangles] = FloaterParam_I_a(v, adjVerts);

		// b) UVs
		auto UVs = FloaterParam_I_b(thetaI, disANDangles);

		/// ====== Stage 2 : Calculate lambda_[i, jk] ======
		weights = FloaterParam_II(UVs);
	}
}