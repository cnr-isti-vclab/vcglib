#include <GL\glew.h>
#include <qgl.h>
#include <wrap\gl\trimesh.h>
#include <wrap\gui\trackball.h>
#include <segmentator.h>
#include <sim/tri_pde_integrator.h>
#include <vcg/complex/trimesh/update/bounding.h>
#include <vcg/complex/trimesh/update/bounding.h>
#include <wrap/io_trimesh/export_ply.h>
#include <segmentform.h>

class SimpleGLWidget: public QGLWidget{

Q_OBJECT

private :
	int _H;
	int _W;
	vcg::Trackball Track;
	double zoom;
	GLdouble projection[16];
	GLdouble modelMatrix[16];
	GLint viewport[4];
	bool _showslides;
	int _numslide;
	bool wire;
	bool blocked;
	bool extForces;
	bool intForces;
	bool resultForces;
	bool continue_int;
	
	//vcg::GlTrimesh<Segmentator::MyTriMesh> *Wrap;
	
public:
	SegmentForm *w;
	SimpleGLWidget( QWidget * parent = 0, const char * name = 0, const QGLWidget * shareWidget = 0, WFlags f = 0 );

	virtual void glDraw();
	//virtual void paintEvent ( QPaintEvent * ) ;
    void resizeGL( int w, int h );
    virtual void mousePressEvent ( QMouseEvent * e );
	virtual void mouseReleaseEvent(QMouseEvent * e );
	virtual void mouseMoveEvent ( QMouseEvent * e );
	virtual void wheelEvent ( QWheelEvent * e );
	//virtual void keyPressEvent(QKeyEvent *k);
	virtual void initializeGL();
	virtual void SaveMatrix();
	virtual void Save();
	void LoadMatrix();
	void drawSlide();
	void SmoothMesh();
	void Step();
	void SetExtractionParameters();
	void WriteInfo();
	void ClearMesh();
	void OpenDirectory();

	
	//virtual void keyPressEvent(QKeyEvent *qk);

	public slots:

	void Open()
	{
		OpenDirectory();
	}

	void ShowSlides()
	{
		_showslides=!_showslides;
		repaint();
	}

	void SetWire()
	{
		wire=!wire;
		repaint();
	}

	void SetShowBlocked()
	{
		blocked=!blocked;
		repaint();
	}

	void ShowExternalForces()
	{
		extForces=!extForces;
		repaint();
	}

	void ShowInternalForces()
	{
		intForces=!intForces;
		repaint();
	}

	void ShowResultForces()
	{
		resultForces=!resultForces;
		repaint();
	}
	
	void Smooth()
	{
		SmoothMesh();
		repaint();
	}

	void SavePly()
	{
		Save();
	}

	void Apply()
	{
		SetExtractionParameters();
	}

	void Extract()
	{
		continue_int=!continue_int;
	}

	void Update()
	{
		Step();
	}
	
	void Clear()
	{
		ClearMesh();
	}
};