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
Revision 1.6  2004/05/14 03:15:09  ponchio
Redesigned partial version.

Revision 1.5  2004/05/12 20:55:18  ponchio
*** empty log message ***

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
#include <vector>
#include <map>

namespace vcg {
  /* A trackball stores a transformation called 'track' that effectively rotate the object.
     the rotation center, and size are kept in center and radius.
   
  */

  class Transform {
  public: 
    Transform();
    Similarityf track;
    
    /// la posizione della track nello spazio di modello. il defgault e' 000
    Point3f center; 
    /// size of the widget in spazio di modello. 
    float   radius; 
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
    ~Trackball();
    void SetIdentity();
    void SetPosition(const Point3f &c, int millisec = 0);
    void SetScale(const float s) {radius=s;};
    void SetTransform(const Transform &transform, int miilisec = 0);

    //operating
    void GetView();
    void Apply();
    void ApplyInverse();
    void Draw();
    void ApplynDraw() { Apply(); Draw(); }
    void Reset();

    // Internal Drawing stuff
    static void DrawCircle ();
    static void DrawPlane();
    static void DrawPlaneHandle();

    //interface
    void MouseDown(int x, int y, /*Button*/ int button);
    void MouseMove(int x, int y); 
    void MouseUp(int x, int y, /*Button */ int button); 
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
    enum Action { NONE = 0,
		  VIEW_ROTATE = 1,
		  // Axis Constrained Rotation 
		  TRACK_ROTATE_X = 3, TRACK_ROTATE_Y = 4, TRACK_ROTATE_Z = 5,
		  // Drag constrained to an axis (trackball axis)
		  DRAG_X = 6,   DRAG_Y = 7,   DRAG_Z = 8,
		  // Drag constrained to a plane
		  DRAG_XY = 9,  DRAG_YZ = 10,  DRAG_XZ = 11,
		  //scale model respect to center of trackball
		  VIEW_SCALE = 12,
		  //scale trackball and model
		  TRACK_SCALE = 13
    };

	
    //protected:
    View<float> camera;

    void SetCurrentAction();
  
    int current_button;
    TrackMode *current_mode;

    std::map<int, TrackMode *> modes;

    Similarityf last_track;
    Similarityf last_view;
    Point3f last_point;
    std::vector<Point3f> Hits;
    bool dragging;
    int button_mask;

    Quaternionf spin;
    bool spinnable;
    bool spinning;
  
    std::list<Transform> history;
    int history_size;

    friend class TrackMode;
  };


}//namespace

#endif
