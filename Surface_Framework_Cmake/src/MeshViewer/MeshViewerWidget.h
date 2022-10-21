#pragma once
#include <QString>
#include "QGLViewerWidget.h"
#include "../PolyMesh/include/PolyMesh/IOManager.h"

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
	// my parameterization enums
	enum class TutteParamType { AVERAGE_WEIGHTED, FLOATER_WEIGHTED };
	// my parameterization functions
	void TutteParam(TutteParamType type);

signals:
	void LoadMeshOKSignal(bool, QString);
public slots:
	void PrintMeshInfo(void);
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
	void DrawBoundingBox(void) const;
	void DrawBoundary(void) const;

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
};
