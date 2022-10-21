#include <QtCore>
#include "MeshViewerWidget.h"

#include <tuple>
#include <Eigen3\Dense>
#include <glm\glm.hpp>

MeshViewerWidget::MeshViewerWidget(QWidget* parent)
	: QGLViewerWidget(parent),
	ptMin(0,0,0),
	ptMax(0,0,0),
	isEnableLighting(true),
	isTwoSideLighting(false),
	isDrawBoundingBox(false),
	isDrawBoundary(false)
{
}

MeshViewerWidget::~MeshViewerWidget(void)
{
}

bool MeshViewerWidget::LoadMesh(const std::string & filename)
{
	Clear();

	bool read_OK = acamcad::polymesh::loadMesh(filename, polyMesh);
	std::cout << "Load mesh from file " << filename << std::endl;
	if (read_OK)
	{
		strMeshFileName = QString::fromStdString(filename);
		QFileInfo fi(strMeshFileName);
		strMeshPath = fi.path();
		strMeshBaseName = fi.baseName();
		UpdateMesh();
		update();
		return true;
	}
	return false;
}

void MeshViewerWidget::Clear(void)
{
	polyMesh->clear();
}

void MeshViewerWidget::UpdateMesh(void)
{
	polyMesh->updateFacesNormal();
	polyMesh->updateMeshNormal();
	polyMesh->updateVerticesNormal();
	if (polyMesh->numVertices() == 0)
	{
		std::cerr << "ERROR: UpdateMesh() No vertices!" << std::endl;
		return;
	}
	ptMin[0] = ptMin[1] = ptMin[2] = DBL_MAX;
	ptMax[0] = ptMax[1] = ptMax[2] = -DBL_MAX;

	for (const auto& vh : polyMesh->vertices())
	{
		auto p = vh->position();
		for (size_t i = 0; i < 3; i++)
		{
			ptMin[i] = ptMin[i] < p[i] ? ptMin[i] : p[i];
			ptMax[i] = ptMax[i] > p[i] ? ptMax[i] : p[i];
		}
	}

	double avelen = 0.0;
	double maxlen = 0.0;
	double minlen = DBL_MAX;
	for (const auto& eh : polyMesh->edges()) {
		double len = eh->length();
		maxlen = len > maxlen ? len : maxlen;
		minlen = len < minlen ? len : minlen;
		avelen += len;
	}

	SetScenePosition((ptMin + ptMax)*0.5, (ptMin - ptMax).norm()*0.5);
	std::cout << "Information of the input mesh:" << std::endl;
	std::cout << "  [V, E, F] = [" << polyMesh->numVertices()<< ", " << polyMesh->numEdges() << ", " << polyMesh->numPolygons() << "]\n";
	std::cout << "  BoundingBox:\n";
	std::cout << "  X: [" << ptMin[0] << ", " << ptMax[0] << "]\n";
	std::cout << "  Y: [" << ptMin[1] << ", " << ptMax[1] << "]\n";
	std::cout << "  Z: [" << ptMin[2] << ", " << ptMax[2] << "]\n";
	std::cout << "  Diag length of BBox: " << (ptMax - ptMin).norm() << std::endl;
	std::cout << "  Edge Length: [" << minlen << ", " << maxlen << "]; AVG: " << avelen / polyMesh->numEdges() << std::endl;
}

bool MeshViewerWidget::SaveMesh(const std::string & filename)
{
	return acamcad::polymesh::writeMesh(filename, polyMesh);
}

bool MeshViewerWidget::ScreenShot()
{
	update();
	QString filename = strMeshPath + "/" + QDateTime::currentDateTime().toString("yyyyMMddHHmmsszzz") + QString(".png");
	QImage image = grabFramebuffer();
	image.save(filename);
	std::cout << "Save screen shot to " << filename.toStdString() << std::endl;
	return true;
}

void MeshViewerWidget::SetDrawBoundingBox(bool b)
{
	isDrawBoundingBox = b;
	update();
}
void MeshViewerWidget::SetDrawBoundary(bool b)
{
	isDrawBoundary = b;
	update();
}
void MeshViewerWidget::EnableLighting(bool b)
{
	isEnableLighting = b;
	update();
}
void MeshViewerWidget::EnableDoubleSide(bool b)
{
	isTwoSideLighting = b;
	update();
}

void MeshViewerWidget::ResetView(void)
{
	ResetModelviewMatrix();
	ViewCenter();
	update();
}

void MeshViewerWidget::ViewCenter(void)
{
	if (polyMesh->numVertices()!=0)
	{
		UpdateMesh();
	}
	update();
}

void MeshViewerWidget::CopyRotation(void)
{
	CopyModelViewMatrix();
}

void MeshViewerWidget::LoadRotation(void)
{
	LoadCopyModelViewMatrix();
	update();
}

void MeshViewerWidget::PrintMeshInfo(void)
{
	std::cout << "Mesh Info:\n";
	std::cout << "  [V, E, F] = [" << polyMesh->numVertices() << ", " << polyMesh->numEdges() << ", " << polyMesh->numPolygons() << "]\n";
	std::cout << "  BoundingBox:\n";
	std::cout << "  X: [" << ptMin[0] << ", " << ptMax[0] << "]\n";
	std::cout << "  Y: [" << ptMin[1] << ", " << ptMax[1] << "]\n";
	std::cout << "  Z: [" << ptMin[2] << ", " << ptMax[2] << "]\n";
	std::cout << "  Diag length of BBox: " << (ptMax - ptMin).norm() << std::endl;
	
}

void MeshViewerWidget::DrawScene(void)
{
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixd(&projectionmatrix[0]);
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixd(&modelviewmatrix[0]);
	//DrawAxis();
	if (isDrawBoundingBox) DrawBoundingBox();
	if (isDrawBoundary) DrawBoundary();
	if (isEnableLighting) glEnable(GL_LIGHTING);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, isTwoSideLighting);
	DrawSceneMesh();
	if (isEnableLighting) glDisable(GL_LIGHTING);
}

void MeshViewerWidget::DrawSceneMesh(void)
{
	if (polyMesh->numVertices() == 0) { return; }
	SetMaterial();
	switch (drawmode)
	{
	case POINTS:
		DrawPoints();
		break;
	case WIREFRAME:
		DrawWireframe();
		break;
	case HIDDENLINES:
		DrawHiddenLines();
		break;
	case FLATLINES:
		DrawFlatLines();
		break;
	case FLAT:
		glColor3d(0.8, 0.8, 0.8);
		DrawFlat();
		break;
	case SMOOTH:
		DrawSmooth();
		break;
	default:
		break;
	}
}

void MeshViewerWidget::DrawPoints(void) const
{
	glColor3d(1.0, 0.5, 0.5);
	glPointSize(5);
	glBegin(GL_POINTS);
	for (const auto& vh : polyMesh->vertices()) {
		glNormal3dv(vh->normal().data());
		glVertex3dv(vh->position().data());
	}
	glEnd();
}

void MeshViewerWidget::DrawWireframe(void) const
{
	glColor3d(0.2, 0.2, 0.2);
	glBegin(GL_LINES);
	for (const auto& eh : polyMesh->edges()) {
		auto heh = eh->halfEdge();
		auto v0 = heh->fromVertex();
		auto v1 = heh->toVertex();
		glNormal3dv(v0->normal().data());
		glVertex3dv(v0->position().data());
		glNormal3dv(v1->normal().data());
		glVertex3dv(v1->position().data());
	}
	glEnd();
}

void MeshViewerWidget::DrawHiddenLines() const
{
	glLineWidth(1.0);
	float backcolor[4];
	glGetFloatv(GL_COLOR_CLEAR_VALUE, backcolor);
	glColor4fv(backcolor);
	glDepthRange(0.01, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	if (glIsEnabled(GL_LIGHTING))
	{
		glDisable(GL_LIGHTING);
		DrawFlat();
		glEnable(GL_LIGHTING);
	}
	else
	{
		DrawFlat();
	}
	glDepthRange(0.0, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glColor3d(.3, .3, .3);
	DrawFlat();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void MeshViewerWidget::DrawFlatLines(void) const
{
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.5f, 2.0f);
	glShadeModel(GL_FLAT);
	//glColor3d(0.8, 0.8, 0.8);
	glColor3d(1.0, 1.0, 1.0);
	DrawFlat();
	glDisable(GL_POLYGON_OFFSET_FILL);
	if (glIsEnabled(GL_LIGHTING))
	{
		glDisable(GL_LIGHTING);
		DrawWireframe();
		glEnable(GL_LIGHTING);
	}
	else
	{
		DrawWireframe();
	}
}

void MeshViewerWidget::DrawFlat(void) const
{
	glBegin(GL_TRIANGLES);
	for (const auto& fh : polyMesh->polyfaces())
	{
		glNormal3dv(fh->normal().data());
		for (const auto& fvh :polyMesh->polygonVertices(fh))
		{
			glVertex3dv(fvh->position().data());
		}
	}
	glEnd();
}

// my own-shading code
void MeshViewerWidget::DrawSmooth() const
{
	glShadeModel(GL_SMOOTH);

	glBindTexture(GL_TEXTURE_2D, glTextureID);
	glEnable(GL_TEXTURE_2D);

	glBegin(GL_TRIANGLES);

	for (const auto& fh : polyMesh->polyfaces())
	{
		for (const auto& fvh : polyMesh->polygonVertices(fh))
		{
			glNormal3dv(fvh->normal().data());
			auto uv = fvh->getTextureUVW().uv;
			float uvScale = 10.0f;
			glTexCoord2f(uv[0] * uvScale, uv[1] * uvScale);
			glVertex3dv(fvh->position().data());
		}
	}

	glEnd();
}

void MeshViewerWidget::DrawBoundingBox(void) const
{
	float linewidth;
	glGetFloatv(GL_LINE_WIDTH, &linewidth);
	glLineWidth(2.0f);
	glColor3d(.3, .7, .3);
	glBegin(GL_LINES);
	for (const auto& i : { 0, 1 })
	{
		for (const auto& j : { 0, 1 })
		{
			for (const auto& k : { 0, 1 })
			{
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(~i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], ~j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], k ? ptMin[2] : ptMax[2]);
				glVertex3d(i ? ptMin[0] : ptMax[0], j ? ptMin[1] : ptMax[1], ~k ? ptMin[2] : ptMax[2]);
			}
		}
	}
	glEnd();
	glLineWidth(linewidth);
}

void MeshViewerWidget::DrawBoundary(void) const
{
	float linewidth;
	glGetFloatv(GL_LINE_WIDTH, &linewidth);
	glLineWidth(2.0f);
	glColor3d(0.1, 0.1, 0.1);
	glBegin(GL_LINES);

	for (const auto& eh : polyMesh->edges()) {
		if (polyMesh->isBoundary(eh)) {
			auto heh = eh->halfEdge();
			auto v0 = heh->fromVertex();
			auto v1 = heh->toVertex();
			glNormal3dv(v0->normal().data());
			glVertex3dv(v0->position().data());
			glNormal3dv(v1->normal().data());
			glVertex3dv(v1->position().data());
		}
	}
	
	glEnd();
	glLineWidth(linewidth);
}

/// ==================== my parameterization functions ====================

enum class UVBoundaryType { 
	POLYGON_CIRCLE, 
	POLYGON_TRIANGLE,
	POLYGON_SQUARE,
	POLYGON_PENTAGON,
	// TODO : uncomment this to implement more complex boundaries, for later experiments
	// POLYGON_STAR,
	// POLYGON_CROSS
} ;

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
	glm::dmat2x2 polygonRotateMat = glm::dmat2({ 0, 1 }, { -1, 0});

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

// return sum theta and { ||xj - xi||, theta_j } of each xj respectively
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

// return UVs of neighbor vertices of v
auto FloaterParam_I_b(double thetaSum, const std::vector<glm::dvec2>& neighborAngleInfo)
{
	int N = neighborAngleInfo.size();

	glm::dvec2 P(0.0, 0.0);
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
		double dP = scale * glm::sqrt(glm::dot(adjUVs[j - 1] - P, adjUVs[j - 1] - P));
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

			assert(!std::isnan(x(0) * x(1) * x(2)));

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

std::vector<double> CalAdjectWeight(acamcad::polymesh::MVert* v, 
	const std::vector<acamcad::polymesh::MVert*>& adjVerts, 
	MeshViewerWidget::TutteParamType type)
{
	std::vector<double> weights;

	switch (type)
	{
	case MeshViewerWidget::TutteParamType::AVERAGE_WEIGHTED:
		AverageParam(v, adjVerts, weights);
		break;
	case MeshViewerWidget::TutteParamType::FLOATER_WEIGHTED:
		FloaterParam(v, adjVerts, weights);
		break;
	}
	return weights;
}

void MeshViewerWidget::TutteParam(TutteParamType type)
{
	using mat = Eigen::MatrixXd;
	using acamcad::polymesh::MVert;
	std::cout << "Tutte's Parameterization\n";

	if (polyMesh->numVertices() == 0)
	{
		std::cerr << "ERROR: TutteParam() No vertices!" << std::endl;
		return;
	}

	/// ====== calculate boundary vertices ======
	auto boundaryVertLists = polyMesh->boundaryVertices();
	std::cout << "Boundary Points [" << boundaryVertLists.size() << "]\n";

	// 1. prepare convex polygon
	int M = boundaryVertLists.size();
	int N = polyMesh->numVertices();

	// TODO : Better to add a function to control the boundary polygon, or a slider?
	auto boundaryUVs = GetBoundaryUVs(M, UVBoundaryType::POLYGON_CIRCLE);
	//auto boundaryUVs = GetBoundaryUVs(M, UVBoundaryType::POLYGON_TRIANGLE);
	//auto boundaryUVs = GetBoundaryUVs(M, UVBoundaryType::POLYGON_SQUARE);
	//auto boundaryUVs = GetBoundaryUVs(M, UVBoundaryType::POLYGON_PENTAGON);

	// 2. prepare matrix (Eigen::Dense)
	mat A = mat::Zero(N, N);
	mat x = mat::Zero(N, 2);
	mat b = mat::Zero(N, 2);

	// a) boundary vertices
	for (int i = 0; i < M; ++i)
	{
		int vertID = boundaryVertLists[i]->index();

		A(vertID, vertID) = 1.0;
		b(vertID, 0) = boundaryUVs[i].x;
		b(vertID, 1) = boundaryUVs[i].y;
	}

	// b) inner vertices
	std::unordered_set<MVert*> boundaryVertsDict(boundaryVertLists.begin(), boundaryVertLists.end());

	for (int i = 0; i < N; ++i)
	{
		auto pVert = polyMesh->vert(i);
		if (boundaryVertsDict.count(pVert))	continue;

		int vID = pVert->index();

		auto adjVerts = polyMesh->vertAdjacentVertices(pVert);
		auto weights = CalAdjectWeight(pVert, adjVerts, type);
		double weightSUM = std::accumulate(weights.begin(), weights.end(), 0.0);
		for (int j = 0; j < adjVerts.size(); ++j) {
			auto adjID = adjVerts[j]->index();
			A(vID, adjID) = weights[j];
			A(vID, vID) = -weightSUM;
		}
	}

	/// 3. solve the equation Ax = b
	x = A.lu().solve(b);

	std::cout << "Tutte's parameterization finished, updating mesh...\n";

	for (int i = 0; i < N; ++i)
	{
		auto pVert = polyMesh->vert(i);
		auto vID = pVert->index();
	
		// TODO : Better to add a function to set vertex position
		//pVert->setPosition(x(vID, 0), x(vID, 1), 0.0);
		pVert->setTexture(x(vID, 0), x(vID, 1));
	}

	std::cout << "Done\n";

	UpdateMesh();
	update();
}
