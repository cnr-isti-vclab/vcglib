#include <GL/glew.h>
#include <qgl.h>
#include <wrap/gl/glwraptetra.h>
#include <wrap/gl/pick.h>
#include <wrap/gui/trackball.h>
#include "myclasses.h"
#include "tetrastats.h"



class MyGLWidget: public QGLWidget{

Q_OBJECT

private :
	int _H;
	int _W;
	vcg::Trackball Track;
	vcg::GLWrapTetra<MyTetraMesh::TetraContainer > *WT;

	GLdouble projection[16];
	GLdouble modelMatrix[16];
	GLint viewport[4];

	
	int modality;//rendering modality
	enum mousemod {MMTrackball, MMSection,MMNavigateSection};//modality of using mouse
	mousemod mouse_modality;
	

	vcg::Trackball TrackClip;
	
/// This are the flags pf info of the mesh that we want to show
	int  _ShowBar;

  enum {
	SIMPLEX        = 0x00000001,	 // show vertex number and tetrahedrons number
	PHYSICS		   = 0x00000002,	 // show also physical information about the mesh
	QUALITY        = 0x00000004,	 // show informations about aspect ratio
	};

	
public:
	
	MyGLWidget( QWidget * parent = 0, const char * name = 0, const QGLWidget * shareWidget = 0, WFlags f = 0 );

	virtual void glDraw();
    void resizeGL( int w, int h );
    virtual void mousePressEvent ( QMouseEvent * e );
	virtual void mouseReleaseEvent(QMouseEvent * e );
	virtual void mouseMoveEvent ( QMouseEvent * e );
	virtual void wheelEvent ( QWheelEvent * e );
	virtual void keyPressEvent(QKeyEvent *k);
	virtual void initializeGL();
	virtual void SaveMatrix();
	void DrawTetraMesh();
	void DrawBox();
	bool ShowTextSimplex();
	bool ShowTextPhysics();
	bool ShowTextQuality();
	void DrawTextInfo();
	void LoadMatrix();

public slots:

	///bounding box visualization modality
	void setBox(){
		modality=0;
		repaint();
	};

	///wireframe modality
	void setWire(){
		modality=1;
		repaint();
	};

	///hiddenlines modality
	void setHidden(){
		modality=2;
		repaint();
	};

	///alternate wire visualization
	void setFlat(){
		modality=3;
		repaint();
	};

	///alternate wire visualization
	void setFlatWire(){
		modality=4;
		repaint();
	};

	///alternate wire visualization
	void setSmooth(){
		modality=5;
		repaint();
	};

	///alternate wire visualization
	void setSmallTetra()
	{
		modality=6;
		repaint();
	};

	//set trackball modality
	void TrackMouseModality()
	{
		mouse_modality=MMTrackball;
	};

	//set trackball modality
	void SectionMouseModality()
	{
		mouse_modality=MMSection;
	};

	///switching to modality of viewing txt info on simplex
	void SwitchTextSimplex()
	{
		if (ShowTextSimplex())
			_ShowBar&=~SIMPLEX;
		else
			_ShowBar|=SIMPLEX;
		repaint();
	};

	///switching to modality of viewing txt info on physics
	void SwitchTextPhysics()
	{
		if (ShowTextPhysics())
			_ShowBar&=~PHYSICS;
		else
		_ShowBar|=PHYSICS;
		repaint();
	};

	///switching to modality of viewing txt info on quality
	void SwitchTextQuality()
	{
		if (ShowTextQuality())
			_ShowBar&=~QUALITY;
		else
		_ShowBar|=QUALITY;
		repaint();
	};

};