//#include <wrap\gl\trimesh.h>
#include <vcg\physics\methods\lem\interface_lem_vertex.h>
//#include <vcg\physics\methods\lem\interface_lem_face.h>
#include <vcg\physics\methods\lem\interface_lem_remesher.h>
//#include <qapplication.h>
//#include <qgl.h>

//#include <bardrawer.h>

#include <simplex\vertex\with\afvn.h>

#include <simplex\face\with\afav.h>
#include <complex\trimesh\base.h>
#include <wrap\io_trimesh\import_ply.h>
#include <complex\trimesh\update\topology.h>
#include <complex\trimesh\update\bounding.h>



//#include "form1.h"


class MyFace;

class MyVertex: public vcg::VertexAFVNd<vcg::DUMMYEDGETYPE,MyFace,vcg::DUMMYTETRATYPE>
{
  public:
	  ///we suppose at maximum 4 bars for each direction
	  vcg::Bar<MyVertex>* B[6];
	  double Delta[6];
	  double Num[6];
	  //int Num;
};

class MyFace:public vcg::FaceAFAV<MyVertex,vcg::DUMMYEDGETYPE,MyFace>
{
	 public:
	  ///we suppose at maximum 4 bars for each direction
	  vcg::Bar<MyFace>* B[6];
};

typedef  vcg::tri::TriMesh< std::vector<MyVertex> ,std::vector<MyFace> > MyTriMesh;

//typedef  vcg::Lem_Face<MyTriMesh> LemType;
typedef  vcg::Lem_Vertex<MyTriMesh> LemType;
typedef vcg::Lem_Remesher<MyTriMesh> LemRemesherType;
LemType LS;
LemRemesherType LR;
vcg::tri::UpdateTopology<MyTriMesh> UT;
MyTriMesh *tm;
vcg::tri::UpdateBounding<MyTriMesh> UB;

//vcg::GLWrapBar<LemType::LemModel::vectBar> *WB;
//vcg::GlTrimesh<MyTriMesh> *glT;

//struct MyGl: public QGLWidget{
//			MyGl( QWidget * parent = 0, const char * name = 0, const QGLWidget * shareWidget = 0, WFlags f = 0 )
//				:QGLWidget(parent,name){}
//			//void QGLWidget::paintEvent ( QPaintEvent * ) [virtual protected]
//				double lr,ud,tz;
//				int cx,cy,z;
//
//			virtual void glDraw(){
//
//        glClearColor(0.2,0.2,0.2,1);
//				glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
//				glMatrixMode(GL_PROJECTION);
//				glLoadIdentity();
//				gluPerspective(45,1,0.01,20);
//				glMatrixMode(GL_MODELVIEW);
//				glLoadIdentity();
//				gluLookAt(0,0,1,0,0,0,0,10,0);
//        
//
//				glTranslatef(0,0,tz);
//				glRotatef(lr,0,1,0);
//				glRotatef(ud,1,0,0);
//
//        glScalef(1/tm->bbox.Diag(),1/tm->bbox.Diag(),1/tm->bbox.Diag());
//        vcg::Point3d p=tm->bbox.Center();
//				glTranslate(-p);
//
//        WB->Draw();
//        WB->DrawMesh<MyTriMesh>(tm);
//
//        /*glT->Draw<vcg::GLW:: DMFlatWire,vcg::GLW:: CMNone,vcg::GLW:: TMNone> ();*/
//        QGLWidget::glDraw();
//				}
//
//      void resizeGL( int w, int h )
//        {
//            //// setup viewport, projection etc.:
//            glViewport( 0, 0, (GLint)w, (GLint)h );
//        }
//		
//	  virtual void mousePressEvent ( QMouseEvent * e ){
//			cx = e->x();
//			cy = e->y();
//			
//			//tr.MouseDown(e->x(),e->y(),(Trackball::Button)(int)(e->button()));
//			//QWidget::mousePressEvent(e);
//			}
//
//		virtual void mouseMoveEvent ( QMouseEvent * e ){
//			//tr.MouseMove(e->x(),e->y());
//			//QWidget::mouseMoveEvent(e);
//			lr+=e->x()-cx;
//			ud-=e->y()-cy;
//			cx = e->x();
//			cy = e->y();
//			repaint();
//			}
//
//	virtual void wheelEvent ( QWheelEvent * e ){
//			tz +=e->delta()/360.f;
//					repaint();
//			QWidget::wheelEvent(e);
//		}
//
//			virtual void initializeGL(){
//
//					
//				GLfloat f[4]={0.2,0.2,0.2,1.f};
//				GLfloat p[4]={3,3,5,0};
//				glLightfv(GL_LIGHT0, GL_AMBIENT,f);
//				glLightfv(GL_LIGHT1, GL_POSITION,p);
//				glLightfv(GL_LIGHT1, GL_DIFFUSE,f);
//				glLightfv(GL_LIGHT1, GL_SPECULAR,f);
//
//				glEnable(GL_LIGHT0);
//				glEnable(GL_LIGHT1);
//				glEnable(GL_LIGHTING);	
//				glEnable(GL_DEPTH_TEST);
//				glDepthFunc(GL_LESS);
//				glPolygonMode(GL_FRONT,GL_FILL);
//				glEnable(GL_BACK);
//				glCullFace(GL_BACK);
//				}
//
//};


int main( int argc, char ** argv )
{    LS=LemType(0.01,0.1,0.01);

	//LS=LemType(0.01,200,0.01);
	LS._SetDir(0);
	LS._SetDir(1);
	LS._SetDir(2);
	LS._SetDir(3);
	LS._SetDir(4);
	LS._SetDir(5);

    tm=new MyTriMesh();
    vcg::tri::io::ImporterPLY<MyTriMesh> Imp=vcg::tri::io::ImporterPLY<MyTriMesh>();
    char *name="cube.ply";
    Imp.Open((*tm),name);

    UT.FaceFace(*tm);
    UT.VertexFace(*tm);
    UB.Box(*tm);
	
	LR._SetDir(0);
	LR._SetDir(1);
	LR._SetDir(2);
	LR._SetDir(3);
	LR._SetDir(4);
	LR._SetDir(5);
	LR.Remesh((*tm),0.5,0.2);
	LS.Init(tm,0.5);
    
	/*WB=new vcg::GLWrapBar<std::vector<LemType::BarType> >(LS.LEM.Bars);*/

	LS.SetTouchedBar(&LS.LEM.Bars[9],1);
    LS.ComputeStep(tm);

   /* glT=new vcg::GlTrimesh<MyTriMesh>();
    glT->m=tm;*/
    /*QApplication a( argc, argv );*/
   /* Form1 w;*/
    /*MyGl *gl = new MyGl(&w);
	  gl->setMinimumSize(800,800);
    w.show();
    a.connect( &a, SIGNAL( lastWindowClosed() ), &a, SLOT( quit() ) );
    return a.exec();*/
}
