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
Revision 1.26  2008/02/24 14:37:00  ponchio
Restored trackball functionality. Not very much tested, and code will need some
cleanup.

Revision 1.25  2008/02/22 18:57:46  benedetti
first attempt to correct after quaternion ToMatrix() inversion (does not work yet)

Revision 1.24  2007/12/21 12:29:34  tarini
*** empty log message ***

Revision 1.23  2007/10/12 14:02:39  corsini
solve memory leak in dtor

Revision 1.22  2007/07/09 22:47:18  benedetti
Removed using namespace std and modified accordingly.

Revision 1.21  2007/06/20 12:59:43  corsini
adjust wheel back-compatibility

Revision 1.20  2007/06/13 17:15:08  benedetti
Added one-level undo system and sticky trackmodes.

Revision 1.19  2007/05/15 15:00:27  benedetti
Moved the drawing code to trackmodes, some other minor changes

Revision 1.18  2007/02/26 01:30:02  cignoni
Added reflection Name

Revision 1.17  2007/01/15 15:04:15  tarini
added "ToAscii" and "SetFromAscii" methods to load/store current trackball status from/to ascii strings
(intended uses: clipboard operations and comments inside png snapshots!)

Revision 1.16  2006/07/26 13:54:45  cignoni
Reversed the direction of wheel scaling and added middle mouse panning

Revision 1.15  2006/02/13 13:15:52  cignoni
Added Scale and Translate methods.
Added many drawing hints and raised the default num. of steps when drawing circles.
Added MouseDown without coords (for remembering changes of keys modifiers)
Added ZMode to the default modes under Alt+left
Added DrawPostApply (to be completed)

Revision 1.14  2005/10/17 01:29:46  cignoni
Main restructuring. Removed the Draw function and slightly changed the meaning of the trackball itself.
See the notes at the beginning of trackball.h

Revision 1.13  2005/04/17 17:48:24  ganovelli
modes deallocation commented (quick and dirty solution..to debug)

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
#include<set>

#include <wrap/gl/math.h>
#include <wrap/gl/space.h>

using namespace vcg;

Transform::Transform() {
  track.SetIdentity();
  radius=1.0f;
  center=Point3f(0,0,0);
}

Trackball::Trackball(): current_button(0), current_mode(NULL), inactive_mode(NULL),
			dragging(false), last_time(0), spinnable(true), spinning(false),
			history_size(10), fixedTimestepMode(false) {
  setDefaultMapping ();
}

Trackball::~Trackball()
{
	ClearModes();
	delete  inactive_mode;
}

void Trackball::ClearModes()
{
	// Note: people ofter maps different keys to the same modes.
	// so we should avoid double deletion of these double referenced modes.
	std::set<TrackMode *> goodModes;
	std::map<int, TrackMode *>::iterator it;
  for(it = modes.begin(); it != modes.end(); it++)
		if ((*it).second) goodModes.insert( (*it).second);

	std::set<TrackMode *>::iterator its;
	for(its = goodModes.begin(); its != goodModes.end(); its++)
			delete *its;

	modes.clear();
}

void Trackball::setDefaultMapping () {
  idle_and_keys_mode = NULL;

  inactive_mode = new InactiveMode ();
	ClearModes();
  modes[0] = NULL;

  modes[BUTTON_MIDDLE | KEY_ALT] =
  modes[BUTTON_LEFT] = new SphereMode ();

  modes[BUTTON_LEFT | KEY_CTRL] = new PanMode ();

  modes[BUTTON_MIDDLE] = new PanMode ();

  modes[WHEEL] =
  modes[BUTTON_LEFT | KEY_SHIFT] = new ScaleMode ();

  modes[BUTTON_LEFT | KEY_ALT] = new ZMode ();

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
}

// the drawing code has been moved to the trackmodes
void Trackball::DrawPostApply() {
	if(current_mode !=NULL){
		current_mode->Draw(this);
	}else{
		if (inactive_mode != NULL) inactive_mode->Draw(this);
  }
}

void Trackball::Apply () {
  glTranslate (center);
  glMultMatrix (track.Matrix());
  glTranslate (-center);
}

void Trackball::Apply(bool ToDraw) {
  Apply();
  if(ToDraw){
    DrawPostApply();
  }
}


void Trackball::ApplyInverse() {
  glTranslate(center);
  glMultMatrix(track.InverseMatrix());
  glTranslate(-center);
}

// T(c) S R T(t) T(-c) => S R T(S^(-1) R^(-1)(c) + t - c)
Matrix44f Trackball::Matrix() const{
  #ifndef VCG_USE_EIGEN
  Matrix44f r; track.rot.ToMatrix(r);
  Matrix44f sr    = Matrix44f().SetScale(track.sca, track.sca, track.sca) * r;
  Matrix44f s_inv = Matrix44f().SetScale(1/track.sca, 1/track.sca, 1/track.sca);
  Matrix44f t     = Matrix44f().SetTranslate(s_inv*r.transpose()*center + track.tra - center);

  return Matrix44f(sr*t);
  #else
  Eigen::Quaternionf rot(track.rot);
  Eigen::Translation3f tr( (1/track.sca) * (rot.inverse() * center) + track.tra - center );
  return ( Eigen::Scaling3f(track.sca) * (rot * tr) ).matrix();
  #endif
}

Matrix44f Trackball::InverseMatrix() const{
  return Inverse(Matrix());
}

void Trackball::Scale(const float s)
{
  track.sca*=s;
}

void Trackball::Translate(Point3f tr)
{
  Quaternionf irot = track.rot;
  irot.Invert();
  track.tra = last_track.tra + irot.Rotate(tr)/track.sca;
}

/***************************************************************/
// DrawCircle () e DrawPlane() have been moved to trackutils.h
// the drawing code has been moved to the trackmodes
/*
void Trackball::DrawCircle() {
  int nside=DH.CircleStep;
  const double pi2=3.14159265*2.0;
  glBegin(GL_LINE_LOOP);
  for(double i=0;i<nside;i++){
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
*/

void Trackball::ToAscii(char* result){
  float * f = (float*) &track;
  sprintf(result, "trackball(%f,%f,%f,%f,%f,%f,%f,%f)",
                  f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7] );
}

bool Trackball::SetFromAscii(const char * st){
  float * f = (float*) &track;
  int res=  sscanf(st, "trackball(%f,%f,%f,%f,%f,%f,%f,%f)",
                  f+0,f+1,f+2,f+3,f+4,f+5,f+6,f+7 );

  return (res==8);
}

// DrawPlaneHandle() e DrawIcon() have been moved to trackutils.h
// the drawing code has been moved to the trackmodes
/*
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

void Trackball::DrawIcon() {
  glPushMatrix();

  glScale(radius);
  /// Here start the real drawing stuff
  float amb[4] ={.3f,.3f,.3f,1.0f};
  float col[4] ={.5f,.5f,.8f,1.0f};
  //float col2[4]={.9f,.9f,1.0f,1.0f};
  glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT | GL_CURRENT_BIT | GL_LIGHTING_BIT);


  if(current_mode == NULL ) glLineWidth(DH.LineWidthStill);
                       else glLineWidth(DH.LineWidthMoving);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glColor(DH.color);

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

  //glColor4f(1.0,.8f,.8f,1.0f);

  glPopAttrib();

  glPopMatrix();
}
*/

void Trackball::Reset() {
  track.SetIdentity();
  undo_track = track;
  std::map<int, TrackMode *>::iterator i;
  for(i = modes.begin(); i != modes.end(); i++){
   TrackMode * mode=(*i).second;
   if(mode!=NULL)
     mode->Reset();
  }
  if (inactive_mode != NULL) inactive_mode->Reset();
 }

//interface
void Trackball::MouseDown(int button) {
  undo_track = track;
  current_button |= button;
  SetCurrentAction();
  Hits.clear();
}
void Trackball::MouseDown(int x, int y, int button) {
  undo_track = track;
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
  undo_track = track;
  current_mode->Apply(this, Point3f(float(x), float(y), 0));
}

bool Trackball::IsAnimating(unsigned int msec){
	bool res;
	if(idle_and_keys_mode == NULL) res=false; else res=idle_and_keys_mode->IsAnimating(this);

	if (!fixedTimestepMode)	{
  	if (msec==0) msec = clock()*1000/CLOCKS_PER_SEC;
	  if (!res) {
  	  last_time = msec;
	  }
	}
	return res;
}

void Trackball::Sync(unsigned int msec) {
	if (!fixedTimestepMode)	Animate(msec);
}

void Trackball::Animate(unsigned int msec){
	unsigned int delta;
	if (fixedTimestepMode) delta=msec;
	else {
		if (msec==0) msec = clock()*1000/CLOCKS_PER_SEC;
    delta = msec -last_time;
	  last_time = msec;
	}
	if(idle_and_keys_mode == NULL) return;
	idle_and_keys_mode->Animate(delta,this);
}

void Trackball::MouseUp(int /* x */, int /* y */, int button) {
  undo_track = track;
	ButtonUp(vcg::Trackball::Button(button));
  //current_button &= (~button);
  //SetCurrentAction();
}

// it assumes that a notch of 1.0 is a single step of the wheel
void Trackball::MouseWheel(float notch)
{
  undo_track = track;
	int buttons = current_button;
	current_button = WHEEL | (buttons&(KEY_SHIFT|KEY_CTRL|KEY_ALT));
	SetCurrentAction();
  if (current_mode == NULL)
  {
    //ScaleMode scalemode;
    //scalemode.Apply (this, notch);
  }
	else
	{
    current_mode->Apply(this, notch);
  }
	current_button = buttons;
	SetCurrentAction();
}

void Trackball::MouseWheel(float notch, int button)
{
  undo_track = track;
  current_button |= button;
  SetCurrentAction();
  if (current_mode == NULL) {
    ScaleMode scalemode;
    scalemode.Apply (this, notch);
  } else {
    current_mode->Apply (this, notch);
  }
  current_button &= (~button);
  SetCurrentAction ();
}

void Trackball::ButtonDown(Trackball::Button button, unsigned int msec) {
	Sync(msec);
  bool old_sticky=false, new_sticky=false;
  assert (modes.count (0));

	Button b=Button(current_button & MODIFIER_MASK);
  if ( ( modes.count (b) ) && ( modes[b] != NULL ) ) old_sticky = modes[b]->isSticky();

  current_button |= button;
	b=Button(current_button & MODIFIER_MASK);
	if ( ( modes.count (b) ) && ( modes[b] != NULL ) ) new_sticky = modes[b]->isSticky();

  if ( !old_sticky && !new_sticky) SetCurrentAction();

}

void Trackball::ButtonUp(Trackball::Button button) {
  bool old_sticky=false, new_sticky=false;
  assert (modes.count (0));

	Button b=Button(current_button & MODIFIER_MASK);
  if ( ( modes.count (b) ) && ( modes[b] != NULL ) ) old_sticky = modes[b]->isSticky();

  current_button &= (~button);
	b=Button(current_button & MODIFIER_MASK);
	if ( ( modes.count (b) ) && ( modes[b] != NULL ) ) new_sticky = modes[b]->isSticky();

  if ( !old_sticky && !new_sticky) SetCurrentAction();
}

void Trackball::Undo(){
  track = undo_track;
  if(current_mode != NULL)
    current_mode->Undo();
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

//navigation interface:
void Trackball::Back(){}
void Trackball::Forward(){}
void Trackball::Home(){}
void Trackball::HistorySize(int /* length */){}

void Trackball::SetCurrentAction ()
{
  //I use strict matching.
  assert (modes.count (0));
  if (!modes.count (current_button & MODIFIER_MASK)) {
    current_mode = NULL;
  } else {
    current_mode = modes[current_button & MODIFIER_MASK];
    if(current_mode != NULL)
      current_mode->SetAction();
  }
  last_point = Point3f (0, 0, -1);
  last_track = track;
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

