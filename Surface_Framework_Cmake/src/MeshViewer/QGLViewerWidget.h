#pragma once
#include <QOpenGLWidget>
#include "../PolyMesh/include/PolyMesh/IOManager.h"

class QOpenGLTexture;

class QGLViewerWidget : public QOpenGLWidget
{
	Q_OBJECT
public:
	QGLViewerWidget(QWidget* _parent = 0);
	virtual ~QGLViewerWidget(void);
private:
	void Init(void);
public:
	QSize minimumSizeHint(void) const override;
	QSize sizeHint(void) const override;
	const double & Radius(void) const;
	const acamcad::MVector3& Center(void) const;
	const double* GetModelviewMatrix(void) const;
	void ResetModelviewMatrix(void);
	void CopyModelViewMatrix(void);
	void LoadCopyModelViewMatrix(void);
	const double* GetProjectionMatrix(void) const;

	enum ProjectionMode{ PERSPECTIVE, ORTHOGRAPHIC };
	void SetProjectionMode(const ProjectionMode &pm);
	const ProjectionMode & GetProjectionMode(void) const;

	enum DrawMode{ POINTS, WIREFRAME, HIDDENLINES, FLATLINES, FLAT, SMOOTH, UV_EMBEDDING };
	void SetDrawMode(const DrawMode &dm);
	const DrawMode& GetDrawMode(void) const;

protected:
	enum MaterialType { MaterialDefault, MaterialGold, MaterialSilver, MaterialEmerald, MaterialTin };
	void SetMaterial(const MaterialType & mattype = MaterialGold) const;
	void SetDefaultLight(void) const;
	void initializeGL(void) override;
	void resizeGL(int w, int h) override;
	void paintGL(void) override;
	virtual void DrawScene(void);
	virtual void mousePressEvent(QMouseEvent*) override;
	virtual void mouseMoveEvent(QMouseEvent*) override;
	virtual void mouseReleaseEvent(QMouseEvent*) override;
	virtual void wheelEvent(QWheelEvent*) override;
	virtual void keyPressEvent(QKeyEvent*) override;
	virtual void keyReleaseEvent(QKeyEvent*) override;
private:
	void Translation(const QPoint & p);
	void Translate(const acamcad::MVector3& trans);

	void Rotation(const QPoint & p);
	void Rotate(const acamcad::MVector3& axis, const double & angle);
	bool MapToSphere(const QPoint & point, acamcad::MVector3& result);
	void UpdateProjectionMatrix(void);
	// ============ hw1: load uv mapping texture ============
	void LoadTexture();
public:
	void SetScenePosition(const acamcad::MVector3& c, const double & r);
	void ViewAll(void);
protected:
	DrawMode drawmode;
	ProjectionMode projectionmode;
	double windowleft;
	double windowright;
	double windowtop;
	double windowbottom;
	Qt::MouseButton mousemode;
	acamcad::MVector3 center;
	double radius;
	std::vector<double> projectionmatrix;
	std::vector<double> modelviewmatrix;
	std::vector<double> copymodelviewmatrix;
	QPoint lastpoint2;
	acamcad::MVector3 lastpoint3;
	bool lastpointok;
	// ============ hw1: store texture ID ============
	unsigned int glTextureID;
private:
	static const double trackballradius;
};
