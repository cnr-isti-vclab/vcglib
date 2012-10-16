#include "myglwidget.h"
#include <vcg\space\tetra3.h>
#include <vcg\space\point3.h>
#include "mainframe.h"

extern MyTetraMesh *tm;
extern TetraStats<MyTetraMesh> Stats;
//extern MainFrame *wp;

bool MyGLWidget::ShowTextSimplex()
{
	return (_ShowBar & SIMPLEX);
}
	
bool MyGLWidget::ShowTextPhysics()
{
	return (_ShowBar & PHYSICS);
}

bool MyGLWidget::ShowTextQuality()
{
	return (_ShowBar & QUALITY);
}


MyGLWidget::MyGLWidget( QWidget * parent, const char * name, const QGLWidget * shareWidget, WFlags f ):
QGLWidget(parent, name)
{
	Track.Reset();
	Track.radius= 1;
	WT=0;
	modality=3;
	mouse_modality=MMTrackball;
	_ShowBar=SIMPLEX;
	grabKeyboard();
}


void MyGLWidget::DrawTextInfo()
	{
		glPushAttrib(0xffffffff);
		
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_SRC_ALPHA);
		glEnable(GL_LIGHTING);
		glEnable(GL_NORMALIZE);
		glEnable(GL_COLOR_MATERIAL);
		glDisable(GL_CLIP_PLANE0);
		glColor4d(0.7,0,0.7,0.6);

		glDepthRange(0.0,0.1);

		glBegin(GL_QUADS);
			glVertex3d(-0.5,-0.5,0);
			glVertex3d(-0.5,-0.3,0);
			glVertex3d(0.5,-0.3,0);
			glVertex3d(0.5,-0.5,0);
		glEnd();
		
		if (Stats.TCurrent()!=0)
		{
			glBegin(GL_QUADS);
				glVertex3d(0.25,0.5,0);
				glVertex3d(0.5,0.5,0);
				glVertex3d(0.5,0.2,0);
				glVertex3d(0.25,0.2,0);
			glEnd();
		}
		

		renderText( (width() - 10) / 2, 15, "a" );

		QFont f( "arial", 12 );
		QFontMetrics fmc( f );
		glColor3d(1,1,1);
		
		QString str="";
		int level=0;
		
		glDisable( GL_LIGHTING );
		glDisable( GL_TEXTURE_2D );

		if (ShowTextSimplex())
		{
			level++;
			str.sprintf( "Tetrahedrons :  %i  Vertex: %i ",tm->tn,tm->vn);
			renderText( 20, height() - level*20, str, f );
		}
		if (ShowTextPhysics())
		{
			level++;
			str.sprintf( "Volume :  %03f   ",Stats.volume);
			renderText( 20, height() - level*20, str, f );
		}
		if (ShowTextQuality())
		{
			level++;
			str.sprintf( "Aspect Ratio :  %03f  ",Stats.ratio);
			renderText( 20, height() - level*20, str, f );
		}

		//at the end i draw the window for informations about a tretrahedron
		if (Stats.TCurrent()!=0)
		{
			str="";
			str.sprintf( "Volume :  %03f  ",Stats.TCurrent()->ComputeVolume());
			renderText( width()-150, 30, str, f );
			str.sprintf( "Aspect Ratio :  %03f  ",Stats.TCurrent()->AspectRatio());
			renderText( width()-150, 50, str, f );
		
		LoadMatrix();
		glColor3d(1,0,0);
		
		glDisable(GL_BLEND);
		//write values of the tetrahedron
		for (int i=0;i<4;i++)
		{
			double x=Stats.TCurrent()->V(i)->P().V(0);//x value of vertex
			double y=Stats.TCurrent()->V(i)->P().V(1);//y value of vertex
			double z=Stats.TCurrent()->V(i)->P().V(2);//z value of vertex
			str.sprintf("%i",i);
			renderText(x,y,z,str,f);
		}
		Stats.TCurrent()->SetS();
		}
		glPopAttrib();
	}


void MyGLWidget::DrawBox()
{
	glPushAttrib(0xffffffff);
	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHTING);
	glDisable(GL_NORMALIZE);
	glColor3d(1,1,1);

	glBegin(GL_LINE_LOOP);
		glVertex(tm->bbox.P(0));
		glVertex(tm->bbox.P(1));
		glVertex(tm->bbox.P(3));
		glVertex(tm->bbox.P(2));
	glEnd();

	glBegin(GL_LINE_LOOP);
		glVertex(tm->bbox.P(4));
		glVertex(tm->bbox.P(5));
		glVertex(tm->bbox.P(7));
		glVertex(tm->bbox.P(6));
	glEnd();

	glBegin(GL_LINE_LOOP);
		glVertex(tm->bbox.P(0));
		glVertex(tm->bbox.P(1));
		glVertex(tm->bbox.P(5));
		glVertex(tm->bbox.P(4));
	glEnd();
	
	glBegin(GL_LINE_LOOP);
		glVertex(tm->bbox.P(3));
		glVertex(tm->bbox.P(2));
		glVertex(tm->bbox.P(6));
		glVertex(tm->bbox.P(7));
	glEnd();

	glPopAttrib();
}

void MyGLWidget::DrawTetraMesh()
{

switch (modality)
{
	case 0:DrawBox();break;
	case 1:WT->Draw<vcg::GLW::DMWire,vcg::GLW::NMFlat,vcg::GLW::CMNone>();break;
	case 2:WT->Draw<vcg::GLW::DMHidden,vcg::GLW::NMFlat,vcg::GLW::CMNone>();break;
	case 3:WT->Draw<vcg::GLW::DMFlat,vcg::GLW::NMFlat,vcg::GLW::CMNone>();break;
	case 4:WT->Draw<vcg::GLW::DMFlatWire,vcg::GLW::NMFlat,vcg::GLW::CMNone>();break;
	case 5:WT->Draw<vcg::GLW::DMFlat,vcg::GLW::NMSmooth,vcg::GLW::CMNone>();break;
	case 6:WT->Draw<vcg::GLW::DMSmallTetra,vcg::GLW::NMFlat,vcg::GLW::CMNone>();break;
}
}

void MyGLWidget::SaveMatrix()
{
	glGetDoublev(GL_PROJECTION_MATRIX ,projection);
	glGetDoublev(GL_MODELVIEW_MATRIX ,modelMatrix);
}

void MyGLWidget::LoadMatrix()
{
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixd(projection);
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixd(modelMatrix);
}

void MyGLWidget::glDraw(){

	glClearColor(0.2,0.2,0.2,1);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45,1,0.01,20);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0,0,1,0,0,0,0,10,0);

	if (tm!=0){
		glPushMatrix();
		
		glScalef(1/tm->bbox.Diag(),1/tm->bbox.Diag(),1/tm->bbox.Diag());

		Track.GetView();
		Track.Apply();
		Track.Draw();

		vcg::Point3d p=tm->bbox.Center();
		glTranslate(-p);
		

		//if not exist crete an instance of wrapper
		if (WT==0)
		{
			WT= new vcg::GLWrapTetra<std::vector<MyTetrahedron> >(tm->tetra);
			WT->SetHint(vcg::GLW::HShrinkFactor, 0.8);
		}
		
		/*glGetDoublev(GL_PROJECTION_MATRIX,proj);
		glGetDoublev(GL_mode_MATRIX,mod);
		glGetDoublev(GL_PROJECTION_MATRIX,);*/
		SaveMatrix();
		DrawTetraMesh();
		glPopMatrix();
		DrawTextInfo();
	}
	QGLWidget::glDraw();
}

void MyGLWidget::resizeGL( int w, int h )
  {
	  // setup viewport, projection etc.:
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	float ratio=(float)w/(float)h;
	gluPerspective(45,ratio,1,20);
	_W=w;
	_H=h;
	glViewport (0, 0, (GLsizei) w, (GLsizei) h);
	glMatrixMode (GL_MODELVIEW);
	repaint();

  }

void MyGLWidget::mousePressEvent ( QMouseEvent * e )
{
if (e->button()==Qt::LeftButton)
{
	MyTetraMesh::TetraIterator ti;
	int face;
	switch(mouse_modality)
	{				
	case MMTrackball:
		Track.MouseDown(e->x(),_H-e->y(),vcg::Trackball::BUTTON_LEFT);
		break;

	case MMSection:
		LoadMatrix();
		vcg::GLPickTetra<MyTetraMesh>::PickNearestTetraFace(e->x(),_H-e->y(),*tm,ti,face);
		if (ti!=0)
		{
			///find external face
			
			/*while (!ti->IsBorderF(face))
			face++;*/
			/*ti->SetS();*/
			vcg::Point3d p0=ti->V(vcg::Tetra::VofF(face,0))->P();
			vcg::Point3d p1=ti->V(vcg::Tetra::VofF(face,1))->P();
			vcg::Point3d p2=ti->V(vcg::Tetra::VofF(face,2))->P();

			//put the trackball on the barycenter of the face
			MyTetraMesh::VertexType::CoordType b=(p0+p1+p2)/3.f;

			WT->AddClipSection(p0,p1,p2);
			TrackClip.Reset();
			TrackClip.radius=1;
			TrackClip.center.V(0)=(float)b.V(0);
			TrackClip.center.V(1)=(float)b.V(1);
			TrackClip.center.V(2)=(float)b.V(2);
			mouse_modality=MMNavigateSection;
			TrackClip.MouseDown(e->x(),_H-e->y(),vcg::Trackball::BUTTON_LEFT);
		}
	break;

	case MMNavigateSection:
		TrackClip.MouseDown(e->x(),_H-e->y(),vcg::Trackball::BUTTON_LEFT);
	break;
	
	}	
}

else if (e->button()==Qt::RightButton)
{
	MyTetraMesh::TetraIterator ti;
	LoadMatrix();
	WT->section.GlClip();
	vcg::GLPickTetra<MyTetraMesh>::PickNearestTetra(e->x(),_H-e->y(),*tm,ti);
	if (ti!=0)
	{
		Stats.SetTetraInfo(&*ti);
	}
}

repaint();
}

void MyGLWidget::mouseReleaseEvent(QMouseEvent * e )
  {
	  Track.MouseUp(e->x(),_H-e->y(),vcg::Trackball::BUTTON_LEFT);
	  repaint();
  }

void MyGLWidget::mouseMoveEvent ( QMouseEvent * e )
	{	
		if (mouse_modality==MMTrackball)
		{ 
			Track.MouseMove(e->x(),_H-e->y());
			repaint();
		}
		else
		if ((mouse_modality==MMNavigateSection)&&(e->state() & Qt::LeftButton))
		{
			LoadMatrix();
			TrackClip.MouseMove(e->x(),_H-e->y());
			TrackClip.GetView();
			TrackClip.Apply();
			WT->section.Transform(TrackClip.track.Matrix());
			repaint();
		}
	}

void MyGLWidget::wheelEvent ( QWheelEvent * e ){
		/*
		if (mouse_modality==MMTrackball)
		{ 
		const int WHEEL_DELTA =120;
		Track.MouseWheel( e->delta()/ float(WHEEL_DELTA) );
		repaint();	
		}else*/
	if (mouse_modality==MMNavigateSection)
		{
			const int WHEEL_DELTA =120;
			float delta= e->delta()/ float(WHEEL_DELTA);
			WT->section.Translate(delta/10);

			///for casting from double to float
			TrackClip.center.V(0)=(float)WT->section.P.V(0);
			TrackClip.center.V(1)=(float)WT->section.P.V(1);
			TrackClip.center.V(2)=(float)WT->section.P.V(2);
			
			repaint();
		}
	}


void MyGLWidget::keyPressEvent(QKeyEvent *k)
{
	if (k->key()==Qt::Key_Escape)//&&((mouse_modality==MMNavigateSection)||(mouse_modality==MMSection)))
	{
		mouse_modality=MMTrackball;
		WT->ClearClipSection();
	}
}

void MyGLWidget::initializeGL(){

	GLfloat f[4]={0.2,0.2,0.2,1.f};
	GLfloat p[4]={3,3,5,0};
	glLightfv(GL_LIGHT0, GL_AMBIENT,f);
	glLightfv(GL_LIGHT1, GL_POSITION,p);
	glLightfv(GL_LIGHT1, GL_DIFFUSE,f);
	glLightfv(GL_LIGHT1, GL_SPECULAR,f);

	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHTING);	
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glPolygonMode(GL_FRONT,GL_FILL);
	glEnable(GL_BACK);
	glCullFace(GL_BACK);

	}
