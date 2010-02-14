/****************************************************************************
 * VCGLib                                                            o o     *
 * Visual and Computer Graphics Library                            o     o   *
 *                                                                _   O  _   *
 * Copyright(C) 2007                                                \/)\/    *
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
Revision 1.1  2007/10/18 08:52:06  benedetti
Initial release.


****************************************************************************/

#include <QtGui>
#include "glarea.h"
#include <wrap/qt/trackball.h>

GLArea::GLArea (QWidget * parent)
          :QGLWidget (parent)
{
	drawmode= SMOOTH;
	GLArea::loadTetrahedron();
}

void GLArea::selectDrawMode(int mode){
	drawmode=DrawMode(mode);
	updateGL();
}

void GLArea::loadMesh(QString fileName)
{	
   int err=vcg::tri::io::ImporterPLY<CMesh>::Open(mesh,(fileName.toStdString()).c_str());
	if(err!=0){
	  const char* errmsg=vcg::tri::io::ImporterPLY<CMesh>::ErrorMsg(err);
          QMessageBox::warning(this,tr("Error Loading Mesh"),QString(errmsg));
	}
	initMesh("Loaded \""+fileName+"\".");
}

void GLArea::loadTetrahedron(){
	vcg::tri::Tetrahedron(mesh);
	initMesh(tr("Tethraedron [builtin]"));
}

void GLArea::loadDodecahedron(){
	vcg::tri::Dodecahedron(mesh);
	initMesh(tr("Dodecahedron [builtin]"));
}

void GLArea::initMesh(QString message)
{
	// update bounding box
	vcg::tri::UpdateBounding<CMesh>::Box(mesh);
	// update Normals
        vcg::tri::UpdateNormals<CMesh>::PerVertexNormalizedPerFace(mesh);
        vcg::tri::UpdateNormals<CMesh>::PerFaceNormalized(mesh);
	// Initialize the opengl wrapper
 	glWrap.m = &mesh;
  	glWrap.Update();
  	updateGL();
  	emit setStatusBar(message);	
}

void GLArea::initializeGL ()
{
  glClearColor(0, 0, 0, 0); 
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_NORMALIZE);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
}

void GLArea::resizeGL (int w, int h)
{
  glViewport (0, 0, (GLsizei) w, (GLsizei) h); 
  initializeGL();
 }

vcg::Point3f EP_right_shoulder(float s_dx, float r, float r1, unsigned int N ,unsigned int i){
    float delta = (r-r1)*0.5/N;
    if(i<N/2){
        float x = s_dx-delta*i;
        return vcg::Point3f(x,0,-sqrt(1.f-x*x/(r*r*0.25))*r1*0.5);}
    else{
        float x = -r/2+delta*i;
        return vcg::Point3f(x,0,sqrt(1.f-x*x/(r*r*0.25))*r1*0.5);}
    }
}

vcg::Point3f CP_right_shoulder(float r, unsigned int N ,unsigned int i){
    float ang = M_PI/N*i;
    return vcg::Point3f(r*vcg::math::Cos(ang),0,-r*vcg::math::Sin(ang))
}

void DrawHE(float r, float r1, float sg=1.0){
   const float delta = 100.f;
    glBegin(GL_LINE_STRIP);
    for(float x = -r/2.f; x < r/2.f+r/(delta-1.f); x+=r/delta)
        glVertex3f(x,0,sg*sqrt(1.f-x*x/(r*r*0.25))*r1*0.5);
    glEnd();
}
void DrawShoulder_d(float dx,float dy,float rad){
    glPushMatrix();
    glTranslatef(dx,dy,0.f);
    gluSphere(gluNewQuadric(),rad,100,100);
    glPopMatrix();
}

void DrawE(float n_r1, float n_r2){
    DrawHE(n_r1,n_r2,1.0);
    DrawHE(n_r1,n_r2,-1.0);

}
void  DrawNeck(float ne_x, float ne_y, float angle, float n_r1, float n_r2){

    glPushMatrix();
    glTranslatef(ne_x,ne_y,0.0);
    glRotatef(angle,1,0,0);

    DrawE(n_r1,n_r2);
    glTranslatef(0,3,0);
    DrawE(n_r1,n_r2);
    glPopMatrix();
}

void DrawBody(){
    DrawHE(50,40,1.0);      // petto basso
    DrawHE(50,20,-1.0);     // schiena basso

    float rad = 4;
    float     s_dx = -15.f,s_dy = 10.f;
    float     s_sx =  15.f,s_sy = 10.f;
    DrawShoulder_d(s_dx,s_dy,rad);
    DrawShoulder_d(s_sx,s_sy,rad);

    float delta_neck = 3.f;
    float ne_x = (s_dx+s_sx)*0.5;
    float ne_y = (s_dy+s_sy+2*rad)*0.5+delta_neck;
    float angle = 10.f;
    float n_r1 = 10.f;
    float n_r2 = 12.f;

    DrawNeck(ne_x,ne_y, angle,n_r1,n_r2);

    // FRR
    float delta_low = (s_dx-r/2)/20;
    for(float x = -r/2; r <= s_dx;r=r+delta_low){
        glVertex(EP(r,r1,1.0,x));
        glVertex(EP(r,r1,1.0,x+delta_low));
        glVertex(CP(r,1.0,x+delta_high));

    }







}
void GLArea::paintGL ()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40, GLArea::width()/(float)GLArea::height(), 0.1, 200);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,100,   0,0,0,   0,1,0);
    track.center=vcg::Point3f(0, 0, 0);
    track.radius= 30;
	track.GetView();
    track.Apply(false);
    glPushMatrix();
    glPolygonMode(GL_FRONT,GL_LINE);


    // disegna il modello
    glDisable(GL_LIGHTING);
    DrawBody();

    glPopMatrix();
   // track.DrawPostApply();
} 

void GLArea::keyReleaseEvent (QKeyEvent * e)
{
  e->ignore ();
  if (e->key () == Qt::Key_Control)
    track.ButtonUp (QT2VCG (Qt::NoButton, Qt::ControlModifier));
  if (e->key () == Qt::Key_Shift)
    track.ButtonUp (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
  if (e->key () == Qt::Key_Alt)
    track.ButtonUp (QT2VCG (Qt::NoButton, Qt::AltModifier));
  updateGL ();
}

void GLArea::keyPressEvent (QKeyEvent * e)
{
  e->ignore ();
  if (e->key () == Qt::Key_Control)
    track.ButtonDown (QT2VCG (Qt::NoButton, Qt::ControlModifier));
  if (e->key () == Qt::Key_Shift)
    track.ButtonDown (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
  if (e->key () == Qt::Key_Alt)
    track.ButtonDown (QT2VCG (Qt::NoButton, Qt::AltModifier));
  updateGL ();
}

void GLArea::mousePressEvent (QMouseEvent * e)
{
  e->accept ();
  setFocus ();
  track.MouseDown (e->x (), height () - e->y (), QT2VCG (e->button (), e->modifiers ()));
  updateGL ();
}

void GLArea::mouseMoveEvent (QMouseEvent * e)
{
  if (e->buttons ()) {
    track.MouseMove (e->x (), height () - e->y ());
    updateGL ();
  }
}

void GLArea::mouseReleaseEvent (QMouseEvent * e)
{
  track.MouseUp (e->x (), height () - e->y (), QT2VCG (e->button (), e->modifiers ()));
  updateGL ();
}

void GLArea::wheelEvent (QWheelEvent * e)
{
  const int WHEEL_STEP = 120;
  track.MouseWheel (e->delta () / float (WHEEL_STEP), QTWheel2VCG (e->modifiers ()));
  updateGL ();
}
