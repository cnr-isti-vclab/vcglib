#include "trackball.h"
#include <wrap/gl/wrap.h>
using namespace vcg;

#include <iostream>  //debug!
using namespace std;

Transform::Transform() {
  local.SetIdentity();
  track.SetIdentity();
}

Trackball::Trackball(): current_button(0), last_x(-1), last_y(-1), dragging(false), 
                        spinnable(true), spinning(false), history_size(10) {
  //qui si aggiungono tutte le actions
  actions[0] = Action(LOCAL, NONE);       
  actions[BUTTON_LEFT] = Action(VIEW, ROTATE);
  actions[BUTTON_LEFT | KEY_CTRL] = Action(VIEW, DRAG_XY);
  actions[BUTTON_LEFT | KEY_SHIFT] = Action(VIEW, SCALE);
  actions[WHEEL] = Action(SCREEN, SCALE);
  current_action = Action(LOCAL, NONE);
}

void Trackball::SetIdentity() {
  local.SetIdentity();
  Reset();
}
void Trackball::SetPosition(const Similarityf &m, int millisec) {
  local = m;
  //millisec ignored at the moment. 
}

//operating
void Trackball::GetView() {
  camera.GetView();
    
  //lets get view matrix  
  Similarityf m = local * last_track;
	Point3f c_obj = Point3f(0, 0, 0) * m;   //coordinate of the center of the trackball in obj coords
	Point3f c_view = camera.Project(c_obj); //center of the trackball in screen coords.	
	Point3f X, Y, Z, C;
	X = camera.UnProject(Point3f(c_view[0] + 100, c_view[1],     c_view[2]));
	Y = camera.UnProject(Point3f(c_view[0],     c_view[1] - 100, c_view[2]));
	Z = camera.UnProject(Point3f(c_view[0],     c_view[1],     c_view[2] + 0.1));	
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

void Trackball::Draw() {
  glPushMatrix();
  glMultMatrix(local.Matrix());
  //Disegnamo un cubo:
  glColor3f(1, 0, 0);
  glScalef(0.5, 0.5, 0.5);
  glBegin(GL_LINE_STRIP);
  glVertex3f(-1, -1, -1);
  glVertex3f( 1, -1, -1);
  glVertex3f( 1,  1, -1);
  glVertex3f(-1,  1, -1);
  glVertex3f(-1,  1,  1);
  glVertex3f( 1,  1,  1);
  glVertex3f( 1, -1,  1);
  glVertex3f(-1, -1,  1);
  glVertex3f(-1, -1,  1);
  glEnd();
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

  Point3f origin = camera.ViewportToScreen(ScreenOrigin());
  Point3f new_point = camera.ViewportToScreen(Point3f(x, y, 0)) - origin;
  Point3f old_point = camera.ViewportToScreen(Point3f(last_x, last_y, 0)) - origin;
  new_point *= 2;
  old_point *= 2;    

  Similarityf u;
  TrackMode *mode = CurrentMode();
  Similarityf new_track = mode->Apply(new_point, u);
  Similarityf old_track = mode->Apply(old_point, u);
  delete mode;

  Invert(old_track);
  new_track = old_track * new_track;

  Similarityf diff;
  switch(current_action.system) {
    case VIEW:            
      u = last_view * Similarityf((-local.tra) * Similarityf(last_track.rot)) * Similarityf(-last_track.tra);
      diff = Inverse(u) * new_track * u;
      break;  
    default: break;
  }
  
  track = diff * last_track;  

} 

void Trackball::MouseUp(int x, int y, Trackball::Button button) { 
  current_button &= (~button);
  SetCurrentAction();
} 

void Trackball::MouseWheel(Trackball::Button notch) {
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
void Trackball::SetSpinnable(bool on){}
bool Trackball::IsSpinnable() {
  return spinnable;
}  
void Trackball::SetSpinning(Quaternionf &spin){}
void Trackball::StopSpinning(){}
bool Trackball::IsSpinning() {
  return spinning;
}  

//interfaccia navigation:
void Trackball::Back(){}
void Trackball::Forward(){}
void Trackball::Home(){}
void Trackball::HistorySize(int lenght){}

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
//return center of trackball in Screen coordinates.
Point3f Trackball::ScreenOrigin() { 
  return camera.Project(ModelOrigin());	
}

//return center of trackball in Model coordinates
Point3f Trackball::ModelOrigin() {
  Similarityf m = local * last_track;
  return Point3f(0, 0, 0) * m;     
}

Matrix44f Trackball::ScreenToModel() {
  return camera.inverse;
}

Similarityf Trackball::ModelToLocal() {
  Similarityf m = local * last_track;
  return m;
}


