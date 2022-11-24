#include <QtCore>
#include "MeshViewerWidget.h"

#include <Eigen\Dense>
#include <Eigen\Sparse>
#include <Eigen\SparseLU>

#include "TutteEmbedding\Util_TutteEmbedding.h"
#include "FreeBoundary\Util_FreeBoundary.h"

MeshViewerWidget::MeshViewerWidget(QWidget* parent)
	: QGLViewerWidget(parent),
	ptMin(0,0,0),
	ptMax(0,0,0),
	isEnableLighting(true),
	isTwoSideLighting(false),
	isDrawBoundingBox(false),
	isDrawBoundary(false),
	isDrawMeshUV(false)
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

		isParameterized = false;
		isProjNewtonSolver = false; 

		paramUVs.setZero();

		edgeCut = freeb::CutEdge(polyMesh);

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

void MeshViewerWidget::ChangeUVState(void)
{
	isDrawMeshUV = !isDrawMeshUV;
	update();
}

void MeshViewerWidget::SetUVScaleIndex(int v)
{
	double scaleUVMin = 1.0;
	double scaleUVMax = 10.0;
	m_UVscale = scaleUVMin + static_cast<double>(v) / 100.0 * (scaleUVMax - scaleUVMin);
	update();
}

void MeshViewerWidget::DrawScene(void)
{
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixd(&projectionmatrix[0]);
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixd(&modelviewmatrix[0]);
	//DrawAxis();
	if (isDrawBoundingBox)	DrawBoundingBox();
	if (isDrawBoundary)		DrawBoundary();
	if (isEnableLighting)	glEnable(GL_LIGHTING);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, isTwoSideLighting);

	glBindTexture(GL_TEXTURE_2D, 0);
	glDisable(GL_TEXTURE_2D);
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
	case UV_EMBEDDING:
		DrawUVEmbedding();
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

// ============ hw1: Smooth shade model  ============
void MeshViewerWidget::DrawSmooth() const
{
	glShadeModel(GL_SMOOTH);

	if (isDrawMeshUV)
	{
		glBindTexture(GL_TEXTURE_2D, glTextureID);
		glEnable(GL_TEXTURE_2D);
	}

	glBegin(GL_TRIANGLES);

	for (const auto& fh : polyMesh->polyfaces())
	{
		for (const auto& fvh : polyMesh->polygonVertices(fh))
		{
			glNormal3dv(fvh->normal().data());
			auto uv = fvh->getTextureUVW().uv;
			glTexCoord2f(uv[0] * m_UVscale, uv[1] * m_UVscale);
			glVertex3dv(fvh->position().data());
		}
	}

	glEnd();
}

void MeshViewerWidget::DrawUVEmbedding() const
{
	glShadeModel(GL_SMOOTH);

	if (isDrawMeshUV)
	{
		glBindTexture(GL_TEXTURE_2D, glTextureID);
		glEnable(GL_TEXTURE_2D);
	}

	auto transfer = (ptMin + ptMax) * 0.5;
	double radius = (ptMin - ptMax).norm() * 0.5 / std::sqrt(3.0);

	{
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0f, 1.0f);

		glBegin(GL_TRIANGLES);
		for (const auto& fh : polyMesh->polyfaces())
		{
			for (const auto& fvh : polyMesh->polygonVertices(fh))
			{
				glNormal3dv(fvh->normal().data());
				auto uv = fvh->getTextureUVW().uv;
				glTexCoord2f(uv[0] * m_UVscale, uv[1] * m_UVscale);
				glVertex3f(uv[0] * radius + transfer.x(), uv[1] * radius + transfer.y(), transfer.z());
			}
		}
		glEnd();

		glDisable(GL_POLYGON_OFFSET_FILL);
	}
	{
		glDisable(GL_LIGHTING);
		glColor3d(0.2, 0.2, 0.2);

		glBegin(GL_LINES);
		for (const auto& eh : polyMesh->edges()) {
			auto heh = eh->halfEdge();
			auto v0 = heh->fromVertex();
			auto v1 = heh->toVertex();
			auto uv0 = v0->getTextureUVW().uv;
			auto uv1 = v1->getTextureUVW().uv;

			glNormal3dv(v0->normal().data());
			glVertex3d(uv0[0] * radius + transfer.x(), uv0[1] * radius + transfer.y(), transfer.z());
			glNormal3dv(v1->normal().data());
			glVertex3d(uv1[0] * radius + transfer.x(), uv1[1] * radius + transfer.y(), transfer.z());
		}
		glEnd();

		glEnable(GL_LIGHTING);
	}
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
	{
		//glLineWidth(2.0f);
		//glColor3d(0.0, 0.0, 1.0);
		//glBegin(GL_LINES);

		//for (const auto& eh : polyMesh->edges()) {
		//	if (polyMesh->isBoundary(eh)) {
		//		auto heh = eh->halfEdge();
		//		auto v0 = heh->fromVertex();
		//		auto v1 = heh->toVertex();
		//		glNormal3dv(v0->normal().data());
		//		glVertex3dv(v0->position().data());
		//		glNormal3dv(v1->normal().data());
		//		glVertex3dv(v1->position().data());
		//	}
		//}
		//glEnd();
	}

	{
		auto transfer = (ptMin + ptMax) * 0.5;
		double radius = (ptMin - ptMax).norm() * 0.5 / std::sqrt(3.0);

		glLineWidth(2.0f);
		glColor3d(1.0, 0.2, 0.3);
		if (drawmode == DrawMode::UV_EMBEDDING)
		{
			glBegin(GL_LINES);
			for (const auto& eh : polyMesh->edges()) {
				size_t edgeIdx = eh->index();
				if (edgeCut[edgeIdx] == 0x01)
				{
					auto heh = eh->halfEdge();
					auto v0 = heh->fromVertex();
					auto v1 = heh->toVertex();

					auto uv0 = v0->getTextureUVW();
					auto uv1 = v1->getTextureUVW();

					glNormal3dv(v0->normal().data());
					glVertex3d(uv0[0] * radius + transfer.x(), uv0[1] * radius + transfer.y(), transfer.z());
					glNormal3dv(v1->normal().data());
					glVertex3d(uv1[0] * radius + transfer.x(), uv1[1] * radius + transfer.y(), transfer.z());
				}
			}
			glEnd();
		} else {
			glBegin(GL_LINES);
			for (const auto& eh : polyMesh->edges()) {
				size_t edgeIdx = eh->index();
				if (edgeCut[edgeIdx] == 0x01)
				{
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
		}
	}

	glLineWidth(linewidth);
}

/// ==================== my parameterization functions ====================
std::vector<double> MeshViewerWidget::CalAdjectWeight(acamcad::polymesh::MVert* v, 
	const std::vector<acamcad::polymesh::MVert*>& adjVerts, 
	MeshViewerWidget::TutteParamType type)
{
	std::vector<double> weights;

	switch (type)
	{
	case TutteParamType::AVERAGE_WEIGHTED:
		tutte::AverageParam(v, adjVerts, weights);
		break;
	case TutteParamType::FLOATER_WEIGHTED:
		tutte::FloaterParam(v, adjVerts, weights);
		break;
	}
	return weights;
}

Eigen::SparseMatrix<double> MeshViewerWidget::TutteParam(TutteParamType type)
{
	// time
	auto timeStart = std::chrono::steady_clock::now();

	using acamcad::polymesh::MVert;
	std::cout << "Tutte's Parameterization\n";

	if (polyMesh->numVertices() == 0)
	{
		std::cerr << "ERROR: TutteParam() No vertices!" << std::endl;
		return Eigen::SparseMatrix<double>();
	}

	/// ====== calculate boundary vertices ======
	auto boundaryVertLists = polyMesh->boundaryVertices();
	std::cout << "Boundary Points [" << boundaryVertLists.size() << "]\n";

	// 1. prepare convex polygon
	int M = boundaryVertLists.size();
	int N = polyMesh->numVertices();


	auto boundaryUVs = tutte::GetBoundaryUVs(M, tutte::UVBoundaryType::POLYGON_CIRCLE);
	//auto boundaryUVs = tutte::GetBoundaryUVs(M, tutte::UVBoundaryType::POLYGON_TRIANGLE);
	//auto boundaryUVs = tutte::GetBoundaryUVs(M, tutte::UVBoundaryType::POLYGON_SQUARE);
	//auto boundaryUVs = tutte::GetBoundaryUVs(M, tutte::UVBoundaryType::POLYGON_PENTAGON);

	// 2. prepare matrix (Eigen::Sparse)
	Eigen::SparseMatrix<double> A(N, N);	A.setZero();
	Eigen::SparseMatrix<double> x(N, 2);	x.setZero();
	Eigen::SparseMatrix<double> b(N, 2);	b.setZero();

	std::vector<Eigen::Triplet<double>> ATriplets;
	ATriplets.reserve(N * 6);

	std::vector<Eigen::Triplet<double>> bTriplets;
	bTriplets.reserve(N);

	// a) boundary vertices
	for (int i = 0; i < M; ++i)
	{
		int vertID = boundaryVertLists[i]->index();

		ATriplets.emplace_back(vertID, vertID, 1.0);
		bTriplets.emplace_back(vertID, 0, boundaryUVs[i].x());
		bTriplets.emplace_back(vertID, 1, boundaryUVs[i].y());
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
		
		for (int j = 0; j < adjVerts.size(); ++j)
			ATriplets.emplace_back(vID, adjVerts[j]->index(), weights[j]);
		ATriplets.emplace_back(vID, vID, -weightSUM);
	}

	A.setFromTriplets(ATriplets.begin(), ATriplets.end());
	b.setFromTriplets(bTriplets.begin(), bTriplets.end());
	
	/// 3. solve the equation Ax = b
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);
	x = solver.solve(b);

	auto timeEnd = std::chrono::steady_clock::now();
	auto ms = std::chrono::duration_cast<std::chrono::microseconds>(timeEnd - timeStart).count() / 1000.0;

	std::cout << "Tutte's parameterization finished with " << ms << "ms, updating mesh...\n";

	for (int i = 0; i < N; ++i)
	{
		auto pVert = polyMesh->vert(i);
		auto vID = pVert->index();
		pVert->setTexture(x.coeff(vID, 0), x.coeff(vID, 1));
	}

	isParameterized = true;

	std::cout << "Done\n";

	UpdateMesh();
	update();

	paramUVs = x;

	return paramUVs;
}

void MeshViewerWidget::ProjNewtonSolver()
{
	if (polyMesh->numVertices() == 0)
	{
		std::cerr << "ERROR: ProjNewtonSolver() No vertices!" << std::endl;
		return;
	}

	if (!isParameterized)
	{
		std::cout << "The mesh hasn't been parameterized yet, run with Tutte's Embedding...\n";
		paramUVs = this->TutteParam(TutteParamType::FLOATER_WEIGHTED);
		isParameterized = true;
	}

	if (!isProjNewtonSolver)
	{
		m_ProjNewtonSolver.PresetMeshUV(polyMesh, paramUVs);
		isProjNewtonSolver = true;
	}

	auto timeStart = std::chrono::steady_clock::now();

	while (m_ProjNewtonSolver.UpdateMeshUV(polyMesh));

	auto timeEnd = std::chrono::steady_clock::now();
	auto ms = std::chrono::duration_cast<std::chrono::microseconds>(timeEnd - timeStart).count() / 1000.0;

	std::cout << "Project Newton Solver finished with " << ms << " ms\n";

	//m_ProjNewtonSolver.SaveEnergies("energies.txt");

	isProjNewtonSolver = false;

	update();
}
