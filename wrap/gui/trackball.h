/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                      
\/)\/    *
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
Revision 1.4  2004/05/07 12:46:08  cignoni
Restructured and adapted in a better way to opengl

Revision 1.3  2004/04/07 10:54:10  cignoni
Commented out unused parameter names and other minor warning related issues

Revision 1.2  2004/03/25 14:55:25  ponchio
Adding copyright.


****************************************************************************/

#ifndef TRACKBALL_H
#define TRACKBALL_H

#include <vcg/math/similarity.h>
#include <wrap/gui/view.h>
#include <wrap/gui/trackmode.h>
#include <list>
#include <map>

namespace vcg {
/* A trackball stores two transformations
the first one, local, is the one 'placing' the trackball somewhere, 
the second one, track, is the one that effectively rotate the object.

The 'local' transformation is the one that contains information about
the rotation center, the size and the orientation of the trackball.
the 'track' is the one implementing the effective transformation.

*/
class Transform {
public: 
  Transform();
  Similarityf track;

  Point3f center;  // la posizione della trackball nello spazio di modello. il defgault e' 000
  float   radius; /// size of the widget in spazio di modello.
  //Quaternion orientation;
};

Transform interpolate(const Transform &a, const Transform &b, float t);
class TrackMode;
class Trackball: public Transform {
public:
  enum Button { BUTTON_NONE   = 0x0000, 
                BUTTON_LEFT   = 0x0001, 
                BUTTON_MIDDLE = 0x0002, 
                BUTTON_RIGHT  = 0x0004, 
                WHEEL         = 0x0008,
		            KEY_SHIFT     = 0x0010, 
                KEY_CTRL      = 0x0020, 
                KEY_ALT       = 0x0040, 
                HANDLE        = 0x0080 };

  Trackball();
  void SetIdentity();
  void SetPosition(const Point3f &c, int millisec = 0);
  void SetScale(const float s) {radius=s;};
  void SetTransform(const Transform &transform, int miilisec = 0);
//  float ScreenScale();

  //operating
  void GetView();
  void Apply();
  void Draw();
  void ApplynDraw() { Apply(); Draw(); }
  void Reset();

  // Internal Drawing stuff
  void DrawCircle ();
  void DrawPlane();
  void DrawPlaneHandle();


  //interface
  void MouseDown(int x, int y, Button button);
	void MouseMove(int x, int y); 
	void MouseUp(int x, int y, Button button); 
	void MouseWheel(Button notch);
  void ButtonUp(Button button);
  void ButtonDown(Button button);

  //default sensitivity 1
  void SetSensitivity(float s);

  //spinning interface
  void SetSpinnable(bool on);
  bool IsSpinnable();  
  void SetSpinning(Quaternionf &spin);
  void StopSpinning();
  bool IsSpinning();  

  //interfaccia navigation:
  void Back();
  void Forward();
  void Home();
  void Store();
  void HistorySize(int lenght);

  //internals

  enum System { LOCAL, VIEW, SCREEN };
  
  enum Motion { NONE = 0, ROTATE = 1, ROTATE_DUMMY = 2,         //really makes sense only in VIEW system
		            ROTATE_X = 3, ROTATE_Y = 4, ROTATE_Z = 5,		    // Axis Constrained Rotation 
		            DRAG_X = 6,   DRAG_Y = 7,   DRAG_Z = 8,					// Drag constrained to an axis (trackball axis)
		            DRAG_XY = 9,  DRAG_YZ = 10,  DRAG_XZ = 11,  		// Drag constrained to a plane
		            SCALE = 12                                      //scale respect to center of trackball		            
  };  

	

  struct Action {
    System system;
    Motion motion;
    Action() {}
    Action(System s, Motion m): system(s), motion(m) {}
  };
  ///Find the current action ussing the current button
  void SetCurrentAction();

//protected:
  View<float> camera;
  Similarityf view;            //Rotate LOCAL coordinate into VIEW coordinates

  float ScreenRadius;
  Point3f ScreenCenter;
  
  int current_button;
  Action current_action;

  TrackMode *CurrentMode();
  std::map<int, Action> actions;

  Similarityf last_track;
  Similarityf last_view;
  int last_x, last_y;
  bool dragging;
  int button_mask;

  Quaternionf spin;
  bool spinnable;
  bool spinning;
  
  std::list<Transform> history;
  int history_size;

   
  //Point3f ScreenOrigin();      //center of trackball in Screen coord   
  Point3f ModelOrigin();       //center of trackball in Model coord
  
  Matrix44f ScreenToModel();  //forse non serve.....
  Similarityf ModelToLocal();
  
  //Point3f ScreenToLocal(const Point3f &p);
  //Point3f LocalToScreen(const Point3f &p);  
};


}//namespace

#endif