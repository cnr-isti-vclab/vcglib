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
Revision 1.4  2004/04/07 10:54:10  cignoni
Commented out unused parameter names and other minor warning related issues

Revision 1.3  2004/03/31 15:08:03  ponchio
Fixed current_action initialization.

Revision 1.2  2004/03/25 14:55:25  ponchio
Adding copyright.


****************************************************************************/

#include <gl/glew.h>

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

Trackball::Trackball(): current_button(0), last_x(-1), last_y(-1), dragging(false), 
                        spinnable(true), spinning(false), history_size(10) {
  //qui si aggiungono tutte le actions
  actions[0] = Action(LOCAL, NONE);       
  actions[BUTTON_LEFT] = Action(VIEW, ROTATE);
  actions[BUTTON_LEFT | KEY_CTRL] = Action(VIEW, DRAG_XY);
  actions[BUTTON_LEFT | KEY_SHIFT] = Action(VIEW, SCALE);
  actions[WHEEL] = Action(SCREEN, SCALE);  
  SetCurrentAction();
}

void Trackball::SetIdentity() {
  track.SetIdentity();
  Reset();
}
void Trackball::SetPosition(const Point3f &c, int /* millisec */) {
  center = c;
}

//operating
void Trackball::GetView() {
  camera.GetView();
    
  //lets get view matrix  
  Similarityf m = last_track;
	Point3f c_obj = m*center;   //coordinate of the center of the trackball in obj coords
	ScreenCenter = camera.Project(c_obj); //center of the trackball in screen coords.	
	ScreenRadius = 10.0f/Distance(center, camera.UnProject(Point3f(ScreenCenter[0] + 10, ScreenCenter[1],     ScreenCenter[2])));

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
  view.FromMatrix(view_axis);	
  //view = view * Inverse(Similarityf(track.rot));   	  
  //spinning ignored
}
void Trackball::Apply() { 
  glMultMatrix(track.Matrix());
}

/***************************************************************/

void Trackball::DrawCircle()
{
	const int nside=18;
	const double pi2=3.14159265*2.0;
	glBegin(GL_LINE_STRIP);
		for(double i=0;i<=nside;i++){
			glNormal3d(cos(i*pi2/nside), sin(i*pi2/nside),  0.0);
			glVertex3d(cos(i*pi2/nside), sin(i*pi2/nside),  0.0);
		}
	glEnd();
	DrawPlaneHandle();
}

void Trackball::DrawPlane()
{
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

void Trackball::DrawPlaneHandle()
{
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
  glTranslate(center);
  glScalef(radius,radius,radius);

  /// Here start the real drawing stuff
  
    float amb[4] ={.3f,.3f,.3f,1.0f};
		float col[4] ={.5f,.5f,.8f,1.0f};
		float col2[4]={.9f,.9f,1.0f,1.0f};
		glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT | GL_CURRENT_BIT | GL_LIGHTING_BIT);
		glLineWidth(2.0);
		glEnable(GL_LIGHTING);
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
		if(current_action.motion==ROTATE_X || current_action.motion==ROTATE_Y || current_action.motion==ROTATE_Z ){
			//Point3d raxt(0,0,0),rax;// compute the rotation axis
			//raxt[current_action.motion-ROTATE_X]=1;
			//RotM.Apply(rax,raxt);
			//
			//glDisable(GL_LIGHTING);
			//glBegin(GL_LINE_STRIP);
			//		glVertex(Manip.c-raxt*TrackballRadius*1.3);
			//		glVertex(Manip.c+raxt*TrackballRadius*1.3);
			//glEnd();
		}
		glPopAttrib();
	
	if(current_action.motion==DRAG_XY)	{
		glPushMatrix();
		//glTranslate(Manip.c);
		DrawPlane();
		glPopMatrix();
		}
	if(current_action.motion==DRAG_XZ)	{	
		glPushMatrix();
		//glTranslate(Manip.c);
		glRotatef(90,1,0,0);
		DrawPlane();
		glPopMatrix();
	}
	if(current_action.motion==DRAG_YZ)	{	
		glPushMatrix();
		//glTranslate(Manip.c);
		glRotatef(90,0,1,0);
		DrawPlane();
		glPopMatrix();
	}

  glPopMatrix();
}

void Trackball::Reset() {
  track.SetIdentity();
}

//interface
void Trackball::MouseDown(int x, int y, Trackball::Button button) {
  current_button |= button;  
  SetCurrentAction();
  last_x = x;
	last_y = y;	     
}

void Trackball::MouseMove(int x, int y) {  
  if(current_action.motion == NONE) return;  
  if(last_x == -1 && last_y == -1) { //changed mode in the middle of moving
    last_x = x;
    last_y = y;
    return;
  }

  TrackMode *mode = CurrentMode();
  mode->tb=this;

  Point3f new_point = Point3f(float(x), float(y), 0);
  Point3f old_point = Point3f(float(last_x), float(last_y), 0);
  Similarityf diff=mode->ComputeFromWindow(old_point,new_point);

  track = last_track*diff;  
} 

void Trackball::MouseUp(int /* x */, int /* y */, Trackball::Button button) { 
  current_button &= (~button);
  SetCurrentAction();
} 

void Trackball::MouseWheel(Trackball::Button /* notch */ ) {
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
  assert(actions.count(0));
  if(!actions.count(current_button))
    current_action = actions[0];
  else
    current_action = actions[current_button];

  last_x = -1;
  last_y = -1;
  last_track = track;
  last_view = view;
}

TrackMode *Trackball::CurrentMode() {
  Point3f x(1, 0, 0), y(0, 1, 0), z(0, 0, 1);
  TrackMode *mode = NULL;
  switch(current_action.motion) {
    case NONE: mode = new TrackMode(); break;
    case ROTATE: mode = new SphereMode(); break;
    case ROTATE_DUMMY: mode = new GravityMode(); break;
    case ROTATE_X: mode = new CylinderMode(x); break;
    case ROTATE_Y: mode = new CylinderMode(y); break;
    case ROTATE_Z: mode = new CylinderMode(z); break;
    case DRAG_X: mode = new LineMode(x); break;
    case DRAG_Y: mode = new LineMode(y); break;
    case DRAG_Z: mode = new LineMode(z); break;
    case DRAG_XY: mode = new PlaneMode(x, y); break; 
    case DRAG_YZ: mode = new PlaneMode(y, z); break;
    case DRAG_XZ: mode = new PlaneMode(z, x); break;
    case SCALE: mode = new ScaleMode(); break;
    default: break;
  }
  return mode;
}
////return center of trackball in Window coordinates.
//Point3f Trackball::ScreenOrigin() { 
//  return camera.Project(ModelOrigin());	
//}


//return center of trackball in Model coordinates
Point3f Trackball::ModelOrigin() {
  return center;     
}

//Matrix44f Trackball::ScreenToModel() {
//  return camera.inverse;
//}
//
//Similarityf Trackball::ModelToLocal() {
//  Similarityf m = local * last_track;
//  return m;
//}


