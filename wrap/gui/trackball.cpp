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
Revision 1.12  2004/12/17 10:28:10  ricciodimare
*** empty log message ***

Revision 1.11  2004/09/28 15:30:12  ponchio
Added a 'else'.

Revision 1.10  2004/09/09 14:38:52  ponchio
#include <gl... -> #include <GL...

Revision 1.9  2004/07/11 22:06:55  cignoni
Added scaling by wheel

Revision 1.8  2004/06/09 14:01:13  cignoni
Heavily restructured. To be completed only rotation works...

Revision 1.7  2004/05/14 03:15:09  ponchio
Redesigned partial version.

Revision 1.6  2004/05/12 20:55:18  ponchio
*** empty log message ***

Revision 1.5  2004/05/07 12:46:08  cignoni
Restructured and adapted in a better way to opengl

Revision 1.4  2004/04/07 10:54:10  cignoni
Commented out unused parameter names and other minor warning related issues

Revision 1.3  2004/03/31 15:08:03  ponchio
Fixed current_action initialization.

Revision 1.2  2004/03/25 14:55:25  ponchio
Adding copyright.


****************************************************************************/
#include <GL/glew.h>
#include "trackball.h"

#include <wrap/gl/math.h>
#include <wrap/gl/space.h>
using namespace vcg;

#include <iostream>  //debug!
using namespace std;

Transform::Transform() {
  track.SetIdentity();
  radius=1.0f;
  center=Point3f(0,0,0);
}

Trackball::Trackball(): current_button(0), current_mode(NULL),
			dragging(false), spinnable(true), spinning(false), 
			history_size(10) {
  //here we add mode
  modes[0]                       = NULL;
  modes[BUTTON_LEFT]             = new SphereMode();
  modes[BUTTON_LEFT | KEY_CTRL]  = new PlaneMode(Plane3f(0, Point3f(1, 0, 0)));
  modes[BUTTON_LEFT | KEY_SHIFT] = new ScaleMode();
  modes[WHEEL]                   = new ScaleMode();
  SetCurrentAction();
}

Trackball::~Trackball() {
  map<int, TrackMode *>::iterator i;
  //for(i = modes.begin(); i != modes.end(); i++)
  //  delete (*i).second;
}

void Trackball::SetIdentity() {
  track.SetIdentity();
  Reset();
}
void Trackball::SetPosition(const Point3f &c, int /* millisec */) {
  center = c;
}

void Trackball::GetView() {
  camera.GetView();
    
  /*  //lets get view matrix  
  Similarityf m = last_track;
  Point3f c_obj = m*center;   //coordinate of the center of the trackball in obj coords
  Point3f ScreenCenter = camera.Project(c_obj); //center of the trackball in screen coords.	
  Point3f ScreenRadius = 10.0f/Distance(center, camera.UnProject(Point3f(ScreenCenter[0] + 10, ScreenCenter[1],     ScreenCenter[2])));
  
  Point3f X, Y, Z, C;
  X = camera.UnProject(Point3f(ScreenCenter[0] + 100, ScreenCenter[1],     ScreenCenter[2]));
  Y = camera.UnProject(Point3f(ScreenCenter[0],     ScreenCenter[1] - 100, ScreenCenter[2]));
  Z = camera.UnProject(Point3f(ScreenCenter[0],     ScreenCenter[1],     ScreenCenter[2] + 0.1f));	
  C = c_obj;
  X = X - C; X.Normalize();
  Y = Y - C; Y.Normalize();
  Z = X ^ Y;	
  
  Matrix44f view_axis; //this is before applying last (or track...)
  view_axis.SetIdentity();
  view_axis.element(0, 0) = X[0]; view_axis.element(0, 1) = X[1]; view_axis.element(0, 2) = X[2];
  view_axis.element(1, 0) = Y[0]; view_axis.element(1, 1) = Y[1]; view_axis.element(1, 2) = Y[2];
  view_axis.element(2, 0) = Z[0]; view_axis.element(2, 1) = Z[1]; view_axis.element(2, 2) = Z[2];
  view.SetIdentity();
  view.FromMatrix(view_axis);	*/
}

void Trackball::Apply() { 
  glTranslate(center);
  glMultMatrix(track.Matrix());
  glTranslate(-center);
}

void Trackball::ApplyInverse() { 
  glTranslate(center);
  glMultMatrix(track.InverseMatrix());
  glTranslate(-center);
}

/***************************************************************/

void Trackball::DrawCircle() {
  int nside=DH.CircleStep;
  const double pi2=3.14159265*2.0;
  glBegin(GL_LINE_STRIP);
  for(double i=0;i<=nside;i++){
    glNormal3d(cos(i*pi2/nside), sin(i*pi2/nside),  0.0);
    glVertex3d(cos(i*pi2/nside), sin(i*pi2/nside),  0.0);
  }
  glEnd();
  DrawPlaneHandle();
}

void Trackball::DrawPlane() {
  const int nl=10;
  float w=5.0f/3.0f;
  float u;
  glBegin(GL_LINES);
  glNormal3f(0.0,0.0,1.0);
  for( u=-w; u<=w+0.01f; u+=2*w/nl){
    glVertex3f(-w,	+u,	0);
    glVertex3f(+w,	+u,	0);
    glVertex3f(+u,	-w, 0);
    glVertex3f(+u,	+w, 0);
  }
  glEnd();
}

void Trackball::DrawPlaneHandle() {
  float r=1.0;
  float dr=r/10.0f;
  glBegin(GL_LINE_STRIP);
  glVertex3f(+r+dr,   +r,   0.0);
  glVertex3f(+r   ,   +r+dr,0.0);
  glVertex3f(+r-dr,   +r,   0.0);
  glVertex3f(+r   ,   +r-dr,0.0);
  glVertex3f(+r+dr,   +r,   0.0);
  glEnd();
  glBegin(GL_LINE_STRIP);
  glVertex3f(-r+dr,   -r,   0.0);
  glVertex3f(-r   ,   -r+dr,0.0);
  glVertex3f(-r-dr,   -r,   0.0);
  glVertex3f(-r   ,   -r-dr,0.0);
  glVertex3f(-r+dr,   -r,   0.0);
  glEnd();
}

void Trackball::Draw() {
  
  glPushMatrix();
  ApplyInverse();
/*  glBegin(GL_POINTS);
  for(vector<Point3f>::iterator vi=Hits.begin();vi!=Hits.end();++vi)
    glVertex(*vi);
  glEnd()*/;
  glPopMatrix();

  glPushMatrix();
  
  glTranslate(center);
  glScalef(radius,radius,radius);
  glScalef(1.0f/track.sca,1.0f/track.sca,1.0f/track.sca);


  /// Here start the real drawing stuff
  
  float amb[4] ={.3f,.3f,.3f,1.0f};
  float col[4] ={.5f,.5f,.8f,1.0f};
  //float col2[4]={.9f,.9f,1.0f,1.0f};
  glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT | GL_CURRENT_BIT | GL_LIGHTING_BIT);
  glLineWidth(2.0);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_BLEND);

  glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,amb);
  glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,col);
  glPushMatrix();
    DrawCircle();
    glPushMatrix();
  
      glRotatef(90,1,0,0);
      DrawCircle();
      glRotatef(90,0,1,0);
      DrawCircle();
      
    glPopMatrix();
  glPopMatrix();
		
  glColor4f(1.0,.8f,.8f,1.0f);

  /*  switch(current_action) {
  case TRACK_ROTATE_X:
  case TRACK_ROTATE_Y:
  case TRACK_ROTATE_Z:
    //Point3d raxt(0,0,0),rax;// compute the rotation axis
    //raxt[current_action.motion-ROTATE_X]=1;
    //RotM.Apply(rax,raxt);
    //
    //glDisable(GL_LIGHTING);
    //glBegin(GL_LINE_STRIP);
    //		glVertex(Manip.c-raxt*TrackballRadius*1.3);
    //		glVertex(Manip.c+raxt*TrackballRadius*1.3);
    //glEnd();
    break;

  case DRAG_XY:
    glPushMatrix();
    //glTranslate(Manip.c);
    DrawPlane();
    glPopMatrix();
    break;

  case DRAG_XY:
    glPushMatrix();
    //glTranslate(Manip.c);
    glRotatef(90,1,0,0);
    DrawPlane();
    glPopMatrix();
    break;

  case DRAG_XY:
    glPushMatrix();
    //glTranslate(Manip.c);
    glRotatef(90,0,1,0);
    DrawPlane();
    glPopMatrix();
    break;
  default:
    break;
    }*/

  glPopAttrib();

  glPopMatrix();
}

void Trackball::Reset() {
  track.SetIdentity();
}

//interface
void Trackball::MouseDown(int x, int y, int button) {
  current_button |= button;  
  SetCurrentAction();
  last_point = Point3f((float)x, (float)y, 0);
  Hits.clear();
}

void Trackball::MouseMove(int x, int y) {  
  if(current_mode == NULL) return;  
  if(last_point[2] == -1) { //changed mode in the middle of moving
    last_point = Point3f((float)x, (float)y, 0);
    return;
  }
  current_mode->Apply(this, Point3f(float(x), float(y), 0));
} 

void Trackball::MouseUp(int /* x */, int /* y */, int button) { 
  current_button &= (~button);
  SetCurrentAction();
} 
 // it assumes that a notch of 1.0 is a single step of the wheel
void Trackball::MouseWheel(float notch  ) {
  if(current_mode == NULL)
  {
    SphereMode tm;
    tm.TrackMode::Apply(this, notch);
  } else
  current_mode->Apply(this, notch);
}

void Trackball::ButtonDown(Trackball::Button button) {
  current_button |= button;  
  SetCurrentAction();
}

void Trackball::ButtonUp(Trackball::Button button) { 
  current_button &= (~button);  
  SetCurrentAction();
}



//spinning interface
void Trackball::SetSpinnable(bool /* on*/ ){}
bool Trackball::IsSpinnable() {
  return spinnable;
}  
void Trackball::SetSpinning(Quaternionf &/* spin*/){}
void Trackball::StopSpinning(){}
bool Trackball::IsSpinning() {
  return spinning;
}  

//interfaccia navigation:
void Trackball::Back(){}
void Trackball::Forward(){}
void Trackball::Home(){}
void Trackball::HistorySize(int /* lenght */){}

void Trackball::SetCurrentAction() {  
  //I use strict matching.
  assert(modes.count(0));
  if(!modes.count(current_button))
    current_mode = NULL;
  else
    current_mode = modes[current_button];

  last_point = Point3f(0, 0, -1);
  last_track = track;
  //  last_view = view;
}

////return center of trackball in Window coordinates.
//Point3f Trackball::ScreenOrigin() { 
//  return camera.Project(ModelOrigin());	
//}


//return center of trackball in Model coordinates
//Point3f Trackball::ModelOrigin() {
//  return center;     
//}

//Matrix44f Trackball::ScreenToModel() {
//  return camera.inverse;
//}
//
//Similarityf Trackball::ModelToLocal() {
//  Similarityf m = local * last_track;
//  return m;
//}


