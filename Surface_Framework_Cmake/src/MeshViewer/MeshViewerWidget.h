#pragma once
#include <QString>
#include "QGLViewerWidget.h"
#include "../PolyMesh/include/PolyMesh/IOManager.h"

#include "EnergyEigensystems\Util_EnergyEigensystems.h"
#include "FreeBoundary\Util_FreeBoundary.h"
#include "BoundaryFirstFlattening\Util_BoundaryFirstFlattening.h"

#include <Eigen\Sparse>

class MeshViewerWidget : public QGLViewerWidget
{
	Q_OBJECT
public:
	MeshViewerWidget(QWidget* parent = 0);
	virtual ~MeshViewerWidget(void);
	bool LoadMesh(const std::string & filename);
	void Clear(void);
	void UpdateMesh(void);
	bool SaveMesh(const std::string & filename);
	bool ScreenShot(void);
	void SetDrawBoundingBox(bool b);
	void SetDrawBoundary(bool b);
	void EnableLighting(bool b);
	void EnableDoubleSide(bool b);
	void ResetView(void);
	void ViewCenter(void);
	void CopyRotation(void);
	void LoadRotation(void);
	// ============ hw1: parameterization enums ============
	enum class TutteParamType { AVERAGE_WEIGHTED, FLOATER_WEIGHTED };
	// ============ hw1: parameterization function entrance ============
	Eigen::SparseMatrix<double> TutteParam(TutteParamType type);
	// ============ hw2: Project Newton Method ============
	void ProjNewtonSolver();
	// ============ hw3: Boundary Free Method ============
	void FreeBoundarySolver();
	// ============ hw4: Boundary First Flattening Method ============
	void BFFSolver();

signals:
	void LoadMeshOKSignal(bool, QString);
public slots:
	void PrintMeshInfo(void);
	void ChangeUVState(void);
	void SetUVScaleIndex(int);

protected:
	virtual void DrawScene(void) override;
	void DrawSceneMesh(void);

private:
	void DrawPoints(void) const;
	void DrawWireframe(void) const;
	void DrawHiddenLines(void) const;
	void DrawFlatLines(void) const;
	void DrawFlat(void) const;
	void DrawSmooth() const;
	void DrawUVEmbedding() const;
	void DrawBoundingBox(void) const;
	void DrawBoundary(void) const;

	// Tutte embedding auxiliary function
	std::vector<double> CalAdjectWeight(acamcad::polymesh::MVert* v,
		const std::vector<acamcad::polymesh::MVert*>& adjVerts,
		MeshViewerWidget::TutteParamType type);

protected:
	acamcad::polymesh::PolyMesh* polyMesh = new acamcad::polymesh::PolyMesh();
	QString strMeshFileName;
	QString strMeshBaseName;
	QString strMeshPath;
	acamcad::MPoint3 ptMin;
	acamcad::MPoint3 ptMax;
	
	bool isEnableLighting;
	bool isTwoSideLighting;
	bool isDrawBoundingBox;
	bool isDrawBoundary;

	bool isDrawMeshUV;
	double m_UVscale = 2.0f;

	// ============ hw2: Project Newton Method member ============
	bool isParameterized = false;
	bool isProjNewtonSolver = false;

	Eigen::SparseMatrix<double> paramUVs;

	eigensys::ProjectNewtonSolver m_ProjNewtonSolver;

	// ============ hw3: Free Boundary member ============
	bool isFreeBoundarySolver = false;
	freeb::FreeBoundarySolver m_FreeBoundarySolver;

	// ============ hw4: Boundary First Flattening member ============
	bff::BFFSolver m_BFFSolver;
};
