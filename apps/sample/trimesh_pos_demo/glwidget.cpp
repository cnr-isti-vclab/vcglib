/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *   
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
/****************************************************************************
  History

$Log: not supported by cvs2svn $
*/

#include <QtGui>
#include <GL/glew.h>
#include <QtOpenGL>

#include <math.h>
#include <wrap/io_trimesh/import_PLY.h>
#include <wrap/gl/picking.h>
#include <wrap/gl/space.h>
#include <wrap/gl/pos.h>
#include <vcg/complex/trimesh/update/bounding.h>
#include <vcg/complex/trimesh/update/normal.h>
#include <vcg/complex/trimesh/update/topology.h>

#include "glwidget.h"

GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(parent)
{
    object = 0;

    trolltechGreen = QColor::fromCmykF(0.40, 0.0, 1.0, 0.0);
    trolltechPurple = QColor::fromCmykF(0.39, 0.39, 0.0, 0.0);
		track.SetPosition(vcg::Point3f(0.0,0.0,0.0));
		track.SetIdentity();
		track.radius = 0.4;
		pos.f=NULL;
}

GLWidget::~GLWidget()
{
    makeCurrent();
    glDeleteLists(object, 1);
}

QSize GLWidget::minimumSizeHint() const
{
    return QSize(50, 50);
}

QSize GLWidget::sizeHint() const
{
    return QSize(400, 400);
}

void GLWidget::LoadTriMesh(QString &namefile)
{
	vcg::tri::io::ImporterPLY<MyStraightMesh>::Open(mesh,namefile.toAscii());
	vcg::tri::UpdateBounding<MyStraightMesh>::Box(mesh);
	vcg::tri::UpdateNormals<MyStraightMesh>::PerFace(mesh);
	vcg::tri::UpdateNormals<MyStraightMesh>::PerVertex(mesh);
	vcg::tri::UpdateTopology<MyStraightMesh>::FaceFace(mesh);
	pos.f=0;
}

void GLWidget::OpenFile(){
	QStringList filters;
	

	QString	fileName = QFileDialog::getOpenFileName(this,tr("Open File"),".", filters.join("\n"));
	
	if (fileName.isEmpty())	return;
	else
	 LoadTriMesh( fileName );

}

void GLWidget::flipV( ){
	if(pos.f) pos.FlipV();
	repaint();
}
void GLWidget::flipE( ){
	if(pos.f) pos.FlipE();
	repaint();
}
void GLWidget::flipF( ){
	if(pos.f) pos.FlipF();
	repaint();
}
void GLWidget::nextE( ){
	if(pos.f) pos.NextE();
	repaint();
}

void GLWidget::initializeGL()
{
    qglClearColor(trolltechPurple.dark());
    object = makeObject();
    glShadeModel(GL_FLAT);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
}

void GLWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,2,   0,0,0,   0,1,0);    
    
 		track.GetView();
    track.Apply();

		
		glScalef(1/glWrap.m->bbox.Diag(),1/glWrap.m->bbox.Diag(),1/glWrap.m->bbox.Diag());
		glTranslate(-glWrap.m->bbox.Center());	

		// to do some picking
		 MyStraightMesh::FaceType* fp=NULL;
		 	if(doPick)
		 	{
				std::vector<MyStraightMesh::FaceType*> res;
				int yes = vcg::Pick<MyStraightMesh::FaceContainer>(pic_x,ScreenH-pic_y+1,mesh.face,res,vcg::glTriangle3<MyStraightMesh::FaceType>,1,1);
				if(yes) 
					{fp = res[0];
						pos.Set(fp,0,fp->V(0));
					}
				doPick=false;
		  	}
		
		glWrap.Draw<vcg::GLW::DMFlatWire,vcg::GLW::CMNone,vcg::GLW::TMNone> ();
 
		if(pos.f!=NULL) {
			glPushAttrib(0xffffffff);
			glDisable(GL_LIGHTING);
			glColor3f(0.0,1.0,0.0);
			glDepthRange(0.0,0.999);
			vcg::GlPos<vcg::face::Pos<MyStraightMesh::FaceType> >::Draw(pos);
			glPopAttrib();
		}
}





void GLWidget::mouseMoveEvent(QMouseEvent *e)
{
			track.MouseMove(e->x(),ScreenH-e->y()+1);
			repaint();

    //if (event->buttons() & Qt::LeftButton) {
    //    setXRotation(xRot + 8 * dy);
    //} else if (event->buttons() & Qt::RightButton) {
    //    setXRotation(xRot + 8 * dy);
    //}
 //   lastPos = event->pos();
}
 void GLWidget::keyPressEvent ( QKeyEvent * e ) {
		if((keypress == Qt::Key_Control)&&(e->key()==Qt::Key_Control))
			keypress = -1;
		else
			keypress = e->key();
 }

GLuint GLWidget::makeObject()
{
    GLuint list = glGenLists(1);
    glNewList(list, GL_COMPILE);

    glBegin(GL_QUADS);

    GLdouble x1 = +0.06;
    GLdouble y1 = -0.14;
    GLdouble x2 = +0.14;
    GLdouble y2 = -0.06;
    GLdouble x3 = +0.08;
    GLdouble y3 = +0.00;
    GLdouble x4 = +0.30;
    GLdouble y4 = +0.22;

    quad(x1, y1, x2, y2, y2, x2, y1, x1);
    quad(x3, y3, x4, y4, y4, x4, y3, x3);

    extrude(x1, y1, x2, y2);
    extrude(x2, y2, y2, x2);
    extrude(y2, x2, y1, x1);
    extrude(y1, x1, x1, y1);
    extrude(x3, y3, x4, y4);
    extrude(x4, y4, y4, x4);
    extrude(y4, x4, y3, x3);

    const double Pi = 3.14159265358979323846;
    const int NumSectors = 200;

    for (int i = 0; i < NumSectors; ++i) {
        double angle1 = (i * 2 * Pi) / NumSectors;
        GLdouble x5 = 0.30 * sin(angle1);
        GLdouble y5 = 0.30 * cos(angle1);
        GLdouble x6 = 0.20 * sin(angle1);
        GLdouble y6 = 0.20 * cos(angle1);

        double angle2 = ((i + 1) * 2 * Pi) / NumSectors;
        GLdouble x7 = 0.20 * sin(angle2);
        GLdouble y7 = 0.20 * cos(angle2);
        GLdouble x8 = 0.30 * sin(angle2);
        GLdouble y8 = 0.30 * cos(angle2);

        quad(x5, y5, x6, y6, x7, y7, x8, y8);

        extrude(x6, y6, x7, y7);
        extrude(x8, y8, x5, y5);
    }

    glEnd();

    glEndList();
    return list;
}

void GLWidget::quad(GLdouble x1, GLdouble y1, GLdouble x2, GLdouble y2,
                    GLdouble x3, GLdouble y3, GLdouble x4, GLdouble y4)
{
    qglColor(trolltechGreen);

    glVertex3d(x1, y1, -0.05);
    glVertex3d(x2, y2, -0.05);
    glVertex3d(x3, y3, -0.05);
    glVertex3d(x4, y4, -0.05);

    glVertex3d(x4, y4, +0.05);
    glVertex3d(x3, y3, +0.05);
    glVertex3d(x2, y2, +0.05);
    glVertex3d(x1, y1, +0.05);
}

void GLWidget::extrude(GLdouble x1, GLdouble y1, GLdouble x2, GLdouble y2)
{
    qglColor(trolltechGreen.dark(250 + int(100 * x1)));

    glVertex3d(x1, y1, +0.05);
    glVertex3d(x2, y2, +0.05);
    glVertex3d(x2, y2, -0.05);
    glVertex3d(x1, y1, -0.05);
}

void GLWidget::normalizeAngle(int *angle)
{
    while (*angle < 0)
        *angle += 360 * 16;
    while (*angle > 360 * 16)
        *angle -= 360 * 16;
}

 void GLWidget:: mousePressEvent(QMouseEvent *e)
{
		if( (keypress==Qt::Key_Control) && (e->button() == Qt::LeftButton) )
					track.MouseDown(e->x(),ScreenH-e->y()+1,vcg::Trackball::KEY_CTRL|vcg::Trackball::BUTTON_LEFT );
			else
				if(e->button() == Qt::LeftButton )
					track.MouseDown(e->x(),ScreenH-e->y()+1,vcg::Trackball::BUTTON_LEFT);
				else
					if(e->button() == Qt::RightButton)
					{
						doPick=true;
						pic_x = e->x();
						pic_y = e->y();
					}
			repaint();
		}

 void GLWidget::mouseReleaseEvent ( QMouseEvent * e ){
			
			if( (keypress==Qt::Key_Control) && (e->button() == Qt::LeftButton) )
					track.MouseUp(e->x(),ScreenH-e->y()+1,vcg::Trackball::KEY_CTRL );
			else
				if(e->button() == Qt::LeftButton )
					track.MouseUp(e->x(),ScreenH-e->y()+1,vcg::Trackball::BUTTON_LEFT);
			repaint();
		}


 void GLWidget::resizeGL(int w,int h){
			ScreenW=w; ScreenH=h;
			glViewport(0,0,w,h);
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			gluPerspective(45,ScreenW/(float)ScreenH,0.01,5);
		}
 void GLWidget::wheelEvent ( QWheelEvent * e ){
			int v =  e->delta()/(float) 120;
			track.MouseWheel(v);
			repaint();
		}
