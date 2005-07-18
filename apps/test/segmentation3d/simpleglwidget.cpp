#include <SimpleGLWidget.h>
#include <qlineedit.h> 
#include <qfiledialog.h> 
#include <qimage.h>
#include <qdir.h> 
#include <qmessagebox.h> 
#include <qslider.h>

extern Segmentator *s;

SimpleGLWidget::SimpleGLWidget( QWidget * parent, const char * name, const QGLWidget * shareWidget, WFlags f ):
QGLWidget(parent, name)
{

	//grabKeyboard();
	_showslides=false;
	blocked=false;
	wire=true;
	extForces=false;
	intForces=false;
	resultForces=false;
	continue_int=false;
	_numslide=0;

	TrackM.center=Point3f(0,0,0);
	TrackM.Reset();
	TrackM.radius= 100.f;

	TrackS.center=Point3f(0,0,0);
	TrackS.Reset();
	TrackS.radius= 100.f;
	
	zoom=1;
	path="";
	/*s=new Segmentator();*/
	//timer = new QTimer(this );
	//QTimer::connect( timer, SIGNAL(timeout()), this, SLOT(Update()) );
 //   timer->start(0); // 2 seconds single-shot timer
	
}

void SimpleGLWidget::SaveMatrix()
{
	glGetDoublev(GL_PROJECTION_MATRIX ,projection);
	glGetDoublev(GL_MODELVIEW_MATRIX ,modelMatrix);
	glGetIntegerv(GL_VIEWPORT,viewport);
}

void SimpleGLWidget::LoadMatrix()
{
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixd(projection);
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixd(modelMatrix);
}

//load as texture the i jpg of the directory
void SimpleGLWidget::LoadTextureJpg(QString p,int level)
{
	int e=glGetError();

	QImage qI=QImage();
	QDir Qd=QDir(p);
	QString qformat;
	QString Path=QString(p);
	Qd.setNameFilter("*.jpg");
	//Qd.setNameFilter("*.bmp");

	Qd.setSorting(QDir::Name);
	QString PathFile=Path;
	PathFile.append(Qd[level]);

	qI.load(PathFile,qformat);
	QImage tx = QGLWidget::convertToGLFormat (qI);

	
	glGenTextures(1, &texName);
	glBindTexture(GL_TEXTURE_2D, texName);
	glTexImage2D( GL_TEXTURE_2D, 0, 3, tx.width(), tx.height(), 0,GL_RGBA, GL_UNSIGNED_BYTE, tx.bits() );
	glEnable(GL_TEXTURE_2D);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);

	e=glGetError();
}


//load the texture corrisponding to a level
void SimpleGLWidget::LoadTextureRaw(int level)
{	
	//glGenTextures(1, &texName);
	//glBindTexture(GL_TEXTURE_2D, texName);
	//glTexImage2D( GL_TEXTURE_2D, 0,GL_LUMINANCE, 512, 512, 0,GL_LUMINANCE, GL_SHORT, &s->V.Data[level*512*512] );
	//glEnable(GL_TEXTURE_2D);
	//glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL );
	//glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	//glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
}

void SimpleGLWidget::drawSlide()
{	
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glPolygonMode(GL_FRONT,GL_FILL);
	glMatrixMode (GL_TEXTURE);
	glLoadIdentity();
	glEnable(GL_TEXTURE_2D);
	/*glEnable(GL_NORMALIZE);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	*/
	glDisable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glEnable(GL_COLOR_MATERIAL);

	float dx=s->BBox().DimX();
	float dy=s->BBox().DimY();

	//to change take scale from segmentator.. better!
	float n=atof(w->S_dist->text());
	float level=((float)_numslide)*n;

	Point3f p0=s->BBox().min;
	p0+=Point3f(0.f,0.f,level);

	Point3f p1=p0+Point3f(dx,0.f,0.f);
	Point3f p2=p0+Point3f(dx,dy,0.f);
	Point3f p3=p0+Point3f(0.f,dy,0.f);

	glColor3d(0.8,0.8,0.8);

	///texture
	//_numslide
	glBegin(GL_QUADS);
		glNormal(Point3d(0,0,1));
		glTexCoord(Point3f(0,1,0));
		glVertex(p0);
		glTexCoord(Point3f(1,1,0));
		glVertex(p1);
		glTexCoord(Point3f(1,0,0));
		glVertex(p2);
		glTexCoord(Point3f(0,0,0));
		glVertex(p3);
	glEnd();
	glPopAttrib();

//glBegin(GL_QUADS);
//for (int x=0;x<((s->V.dimX())-1);x++)
//	for (int y=0;y<((s->V.dimY())-1);y++)
//	{
//		glNormal(Point3d(0,0,1));
//		Point3<int> p0=Point3<int>(x,y,_numslide);
//		double color=((double)s->V.getAt(p0))/256.f;
//		glColor3d(color,color,color);
//		glVertex(p0);
//
//		Point3<int> p1=Point3<int>(x+1,y,_numslide);
//		color=((double)s->V.getAt(p1))/256.f;
//		glColor3d(color,color,color);
//		glVertex(p1);
//
//		Point3<int> p2=Point3<int>(x+1,y+1,_numslide);
//		color=((double)s->V.getAt(p2))/256.f;
//		glColor3d(color,color,color);
//		glVertex(p2);
//
//		Point3<int> p3=Point3<int>(x,y+1,_numslide);
//		color=((double)s->V.getAt(p3))/256.f;
//		glColor3d(color,color,color);
//		glVertex(p3);
//	}
//glEnd();		
}

void drawForces(Segmentator::Part_VertexContainer *pv,int typeForce)
{
Segmentator::Part_VertexContainer::iterator vi;
glPushAttrib(GL_ALL_ATTRIB_BITS);
glLineWidth(0.3f);
glDisable(GL_NORMALIZE);
glDisable(GL_LIGHTING);
glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);

if (typeForce==0)
	glColor3d(1,0,0);
else
if (typeForce==1)
	glColor3d(0,0,1);
else
if (typeForce==2)
	glColor3d(0,1,0);

for (vi=pv->begin();vi<pv->end();++vi)
{
	glBegin(GL_LINE_LOOP);
	vcg::glVertex((*vi).P());
	if (typeForce==0)
		vcg::glVertex((*vi).P()+((*vi).IntForce()*4.f));
	else
	if (typeForce==1)
		vcg::glVertex((*vi).P()+((*vi).ExtForce()*4.f));
	else
	if (typeForce==2)
		vcg::glVertex((*vi).P()+(((*vi).ExtForce()+(*vi).IntForce())*4.f));
	glEnd();
}
glPopAttrib();
}

void SimpleGLWidget::Save()
{
	QString filename = QFileDialog::getSaveFileName("prova.ply",
													"Ply files (*.ply)",
													this,
													"save file dialog",
													"Choose a filename to save under" );
	if (filename!=NULL)
	{
		const char *path_save=filename.ascii();
		vcg::tri::io::ExporterPLY<Segmentator::MyTriMesh>::Save((*s->m),path_save);	
	}
	
}

bool TimeSelfIntersection()
{
	static clock_t time=0;
	clock_t elapsedsecs=abs(time-clock());
	if (elapsedsecs>500)
	{
		time=clock();
		return true;
	}
	return false;
}

//void SimpleGLWidget::WriteInfo()
//{
//	
//	if (s!=0)
//	{
//		
//		glPushAttrib(0xffffffff);
//		
//		glEnable(GL_BLEND);
//		glBlendFunc(GL_SRC_ALPHA, GL_SRC_ALPHA);
//		glEnable(GL_LIGHTING);
//		glEnable(GL_NORMALIZE);
//		glEnable(GL_COLOR_MATERIAL);
//		glDisable(GL_CLIP_PLANE0);
//		glColor4d(0.7,0,0.7,0.6);
//
//		glDepthRange(0.0,0.1);
//
//		glBegin(GL_QUADS);
//			glVertex3d(-0.5,-0.5,0);
//			glVertex3d(-0.5,-0.3,0);
//			glVertex3d(0.5,-0.3,0);
//			glVertex3d(0.5,-0.5,0);
//		glEnd();
//		
//		renderText( (width() - 10) / 2, 15, "a" );
//
//		QFont f( "arial", 12 );
//		QFontMetrics fmc( f );
//		glColor3d(1,1,1);
//		
//		QString str="";
//		int level=0;
//		
//		glDisable( GL_LIGHTING );
//		glDisable( GL_TEXTURE_2D );
//
//		level++;
//		str.sprintf( "Triangles :  %i  Vertex: %i ",s->m.fn,s->m.vn);
//		renderText( 20, height() - level*20, str, f );
//	}
//}

void SimpleGLWidget::SmoothMesh()
{
	s->Smooth();
}

void SimpleGLWidget::glDraw(){
	//glClearColor(0.2f,0.2f,0.2f,1.f);
	glClearColor(1.f,1.f,1.f,1.f);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45,1,0.01,20);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0,0,1,0,0,0,0,10,0);

	if (s!=0){
		
		glPushMatrix();
		if (_showslides)
		{
			vcg::Point3f p=Point3f((float)s->BBox().Center().V(0),(float)s->BBox().Center().V(1),(float)s->BBox().Center().V(2));
			TrackS.radius=s->BBox().Diag();
			TrackS.GetView();
			TrackS.Apply();
			TrackS.Draw();
			glScalef(1.f/s->BBox().Diag(),1.f/s->BBox().Diag(),1.f/s->BBox().Diag());
			//glScalef(GLfloat(zoom),GLfloat(zoom),GLfloat(zoom));
			glTranslate(-p);
			//save transformation matrixes
			SaveMatrix();
		}
		else
		{
			if (s->m!=NULL)
			{
				//vcg::tri::UpdateBounding<Segmentator::MyTriMesh>::Box(*(s->m));
#ifndef _TORUS
				vcg::Point3f p=s->m->bbox.Center();
#endif
				TrackM.GetView();
				TrackM.Apply();
				TrackM.Draw();
				//glScalef(1.f/s->m->bbox.Diag(),1.f/s->m->bbox.Diag(),1.f/s->m->bbox.Diag());
#ifndef _TORUS
				glTranslate(-p);
#endif
				//glTranslate(-CenterExtraction);
				//glScalef(1.f/s->m->bbox.Diag(),1.f/s->m->bbox.Diag(),1.f/s->m->bbox.Diag());
			}
		}
		
		glEnable(GL_NORMALIZE);
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_COLOR_MATERIAL);
		if (_showslides)
			drawSlide();

		if (intForces)
			drawForces(&s->P_Vertex,0);
		if (extForces)
			drawForces(&s->P_Vertex,1);
		if (resultForces)
			drawForces(&s->P_Vertex,2);

		//draw faces
		glEnable(GL_NORMALIZE);
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_COLOR_MATERIAL);
		glDisable(GL_TEXTURE_2D);

		Segmentator::MyTriMesh::FaceIterator fi;
		

		if (wire)
		{
			glDisable(GL_NORMALIZE);
			glDisable(GL_LIGHTING);
			glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
		}
		else
			glPolygonMode(GL_FRONT,GL_FILL);

		if (s->m!=NULL)
		{
			glBegin(GL_TRIANGLES);
			for (fi=s->m->face.begin();fi<s->m->face.end();fi++)
			{
				if(!(*fi).IsD())
				{
					/*if ((!fi->V(0)->IsR())||(!fi->V(1)->IsR())||(!fi->V(2)->IsR()))
					assert(0);
					*/

					//glColor3d(0.4,0.8,0.8);
					if (wire)
					{
						glColor3d(0,0,0);
						glLineWidth(2);
					}
					else 
						glColor3d(0.4,0.8,0.8);
					if ((fi->intersected)&&(wire))
						glColor3d(1.f,0,0);

					if (((blocked)&&(!fi->IsBlocked()))||(!blocked))
					{
						if (!wire)
						{
							vcg::glNormal(fi->V(0)->N());
							vcg::glVertex(fi->V(0)->P());
							vcg::glNormal(fi->V(1)->N());
							vcg::glVertex(fi->V(1)->P());
							vcg::glNormal(fi->V(2)->N());
							vcg::glVertex(fi->V(2)->P());
						}
						else
						{
							vcg::glVertex(fi->V(0)->P());
							vcg::glVertex(fi->V(1)->P());
							vcg::glVertex(fi->V(2)->P());
						}

					}
				}
			}
		}
		glEnd();
	glPopMatrix();
	//WriteInfo();
		
	}
	QGLWidget::glDraw();
}

///open the directiry and initialize dataset
void SimpleGLWidget::OpenDirectory()
{
	QString filename = QFileDialog::getExistingDirectory(
              ".",
              this,
              "open file dialog"
              "Choose a Directory" );
	if (filename!=NULL)
	{
		rawImage=false;
		filename+="/";
		path=filename;
		const char *pa=filename.ascii();
		char *p=(char*)pa;
		s->LoadFromDir(p,"prova.txt");
		LoadTextureJpg(p,0);
		_showslides=true;
		w->SlidesButton->setOn(true);
		repaint();
	}
}

///open the directiry and initialize dataset
void SimpleGLWidget::OpenRaw()
{
	QString filename = QFileDialog::getOpenFileName(
			  "D:/Endocas/",
              "Raw Files(*.raw)",
              this,
              "open file dialog"
              "Choose a Raw file" );

	if (filename!=NULL)
	{
		rawImage=true;
		const char *pa=filename.ascii();
		char *p=(char*)pa;
		s->LoadFromRaw(p);
		LoadTextureRaw(0);
		_showslides=true;
		w->SlidesButton->setOn(true);
		repaint();
	}
}

void SimpleGLWidget::resizeGL( int w, int h )
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

void SimpleGLWidget::ClearMesh()
{
	s->m->Clear();	
	repaint();
}

void SimpleGLWidget::UpdateBBMesh()
{
	vcg::tri::UpdateBounding<Segmentator::MyTriMesh>::Box(*(s->m));
}

void SimpleGLWidget::mousePressEvent ( QMouseEvent * e )
{
	if (e->button()==Qt::LeftButton )
	{

		if ((!_showslides)&&(s->m!=NULL))
		{
			/*vcg::tri::UpdateBounding<Segmentator::MyTriMesh>::Box(*(s->m));
			TrackM.radius=s->m->bbox.Diag();
			TrackM.center=s->m->bbox.Center();
			*/
			UpdateBBMesh();
			TrackM.MouseDown(e->x(),_H-e->y(),vcg::Trackball::BUTTON_LEFT);
#ifndef _TORUS
			CenterExtraction=s->m->bbox.Center();
#endif
		}
		else if (_showslides)
			TrackS.MouseDown(e->x(),_H-e->y(),vcg::Trackball::BUTTON_LEFT);
	}
	else
	//test mass spring model
	if ((e->button()==Qt::RightButton)&&(_showslides))
	{
		float winz;
		double x;
		double y;
		double z;
		//LoadMatrix();
		glReadPixels(e->x(),_H-e->y(),1,1,GL_DEPTH_COMPONENT,GL_FLOAT,&winz);
		gluUnProject(e->x(),_H-e->y(),winz,modelMatrix,projection,viewport,&x,&y,&z);

		SetExtractionParameters();
		s->InitSegmentation(Point3f(x,y,z));
		w->slider1->setValue(s->gray_init);

		UpdateBBMesh();
		//vcg::tri::UpdateBounding<Segmentator::MyTriMesh>::Box(*(s->m));
		repaint();
	}
	//vcg::tri::UpdateBounding<Segmentator::MyTriMesh>::Box(s->m);
}

void SimpleGLWidget::setColor()
{
	s->gray_init=w->slider1->value();
}

///only for debugghing
void SimpleGLWidget::keyPressEvent(QKeyEvent *k)
{
	s->AutoStep();
	repaint();
}

void SimpleGLWidget::wheelEvent(QWheelEvent *e)
{
	if (!_showslides)
	{
	/*	zoom+=e->delta()/120.f;
		repaint();*/
		TrackM.MouseWheel(e->delta()/120.f);
		repaint();
	}
	else
	{
		int oldnum=_numslide;
		_numslide+=e->delta()/120.f;
		if ((_numslide<0)||(_numslide>=s->V.dimZ()))
				_numslide=oldnum;
		if (s!=0)
		{
			if (!rawImage)
				LoadTextureJpg(path,_numslide);
			else
				LoadTextureRaw(_numslide);
		}
	}

	repaint();
}

void SimpleGLWidget::mouseReleaseEvent(QMouseEvent * e )
  {
	  if (!_showslides)
		TrackM.MouseUp(e->x(),_H-e->y(),vcg::Trackball::BUTTON_LEFT);
	  else
		TrackS.MouseUp(e->x(),_H-e->y(),vcg::Trackball::BUTTON_LEFT);

	  repaint();
  }

void SimpleGLWidget::Step()
	{
		if ((s!=0)&&(continue_int))
			{
				s->AutoStep();
				repaint();
			}
	}

//void SimpleGLWidget::MarchingCubesExtraction()
//	{
//		if (s!=0)
//			{
//				s->MarchingCubeExtraction();
//				repaint();
//			}
//	}

void SimpleGLWidget::MarchingCube()
{
	if (s->m->fn>0)
		s->Resample();

	//vcg::tri::io::ExporterPLY<Segmentator::MyTriMesh>::Save((*s->new_m),"d:/march.ply");	
	//vcg::tri::io::ExporterPLY<Segmentator::MyTriMesh>::Save((*s->m),"d:/march.ply");	
}

void SimpleGLWidget::SetExtractionParameters()
	{
		float mass=atof(w->M_particles->text());
		float k_elanst=atof(w->K_elanst->text());
		float dihedral=atof(w->D_angle->text());
		float timestep=atof(w->T_step->text());
		float edge=atof(w->E_size->text());
		float tolerance=atof(w->Tolerance->text());
		float dinstance=atof(w->S_dist->text());
		//int color =atoi(w->Color->text());
		///to modify

		s->SetSegmentParameters(tolerance,mass,k_elanst,dihedral,timestep,edge,Point3f(1.f,1.f,dinstance));
	}

void SimpleGLWidget::mouseMoveEvent ( QMouseEvent * e )
	{
		if (_showslides)
			TrackS.MouseMove(e->x(),_H-e->y());
		else
			TrackM.MouseMove(e->x(),_H-e->y());

		repaint();
	}


void SimpleGLWidget::initializeGL(){

	GLfloat f[4]={0.2f,0.2f,0.2f,1.f};
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
