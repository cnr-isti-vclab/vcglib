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
Revision 1.14  2007/05/15 15:00:47  benedetti
Moved the drawing code to trackmodes, some other minor changes

Revision 1.13  2007/02/26 01:30:02  cignoni
Added reflection Name

Revision 1.12  2007/01/15 15:04:15  tarini
added "ToAscii" and "SetFromAscii" methods to load/store current trackball status from/to ascii strings
(intended uses: clipboard operations and comments inside png snapshots!)

Revision 1.11  2006/08/23 15:40:57  marfr960
*** empty log message ***

Revision 1.10  2006/02/13 13:15:52  cignoni
Added Scale and Translate methods.
Added many drawing hints and raised the default num. of steps when drawing circles.
Added MouseDown without coords (for remembering changes of keys modifiers)
Added ZMode to the default modes under Alt+left
Added DrawPostApply (to be completed)

Revision 1.9  2005/10/17 01:29:46  cignoni
Main restructuring. Removed the Draw function and slightly changed the meaning of the trackball itself.
See the notes at the beginning of trackball.h

Revision 1.8  2004/07/11 22:06:56  cignoni
Added scaling by wheel

Revision 1.7  2004/06/09 14:01:13  cignoni
Heavily restructured. To be completed only rotation works...

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
/****************************************************************************
Short usage note:

The trackball is a manipulator of an object

Center specify the center of rotation and scaling of the trackball and usually 
is set by the program and do not interactively change
Radius specify the radius of the interactive ball shaped icon to specify rotation.
It is in absolute unit but it should be in screen related units like the previoous 
one it is not changed during interaction.

When you specify a traslation with the trackball the trackball center remain UNCHANGED, 
in other words it means that the object move out of the trackball icon. 
Similarly when you apply a scaling the size of the iconshaped ball do not change.


Typical use:

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60, float(width())/float(height()), 1, 100);
	  glMatrixMode(GL_MODELVIEW);
	  glLoadIdentity();
	  gluLookAt(0,0,3,   0,0,0,   0,1,0);        
    
    trackball.center=Point3f(0, 0, 0);
    trackball.radius= 1;
    
    trackball.GetView();
    trackball.Apply();
        
    float d=1.0f/mesh.bbox.Diag();
    glScale(d);
    glTranslate(-mesh.bbox.Center());
    mesh->Render();

Note on the typical use:
Perspective and gllookat are choosed to frame the origin centered 1-radius 
trackball.
The final scale and translate are just to fit a generic mesh to the 1sized 
origin centered where the trackball stays box.
The trackball works also on Orthographic projections 
BUT that are not centered around origin (just move it back along the Z)

****************************************************************************/

#ifndef TRACKBALL_H
#define TRACKBALL_H

#include <vcg/math/similarity.h>
#include <vcg/space/color4.h>
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
    
    /// track position in model space. default is 0,0,0
    Point3f center; 
    /// size of the widget in model space.
    float   radius; 
  };

  Transform interpolate(const Transform &a, const Transform &b, float t);

  class TrackMode;
   class Trackball: public Transform {
  public:
// the drawing code has been moved to the trackmodes
//  class DrawingHint {

// DrawingHint DH;

  
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
    void SetTransform(const Transform &transform, int millisec = 0);
    void Translate(Point3f tr);
    void Scale(const float f);


    //operating
    void GetView();
    void Apply(bool Draw); 
    void Apply (); 
    void DrawPostApply();
    void ApplyInverse();
    // DrawIcon() has been moved to trackutils.h
    //void DrawIcon();
    void Reset();

    // DrawCircle (), DrawPlane(), DrawPlaneHandle() has been moved to trackutils.h
    // the drawing code has been moved to the trackmodes
    // void DrawCircle ();
    // void DrawPlane();
    // void DrawPlaneHandle();

    //interface
    void MouseDown(/*Button*/ int button);
    void MouseDown(int x, int y, /*Button*/ int button);
    void MouseMove(int x, int y); 
    void MouseUp(int x, int y, /*Button */ int button); 
    void MouseWheel(float notch);  // it assumes that a notch of 1.0 is a single step of the wheel
    void MouseWheel (float notch, /*Button */ int button);
    void ButtonUp(Button button);
    void ButtonDown(Button button);
    void Undo();

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

/*    //internals  // commented out no more used this stuff!
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
*/
    // loads/stores current status from/to ascii stings
    void ToAscii(char * st);
    bool SetFromAscii(char * st);
	
    //protected:
    View<float> camera;

    void SetCurrentAction();
  
    int current_button;
    TrackMode *current_mode;

	// inactive_mode is used to draw the inactive trackball
    // can be assigned, for example, to draw an area or a path
    // even when the user is not interacting with it
    TrackMode *inactive_mode;
   
    // reset modes to default mapping.
    void setDefaultMapping ();

    std::map<int, TrackMode *> modes;

    Similarityf last_track;
    
    // undo_track and last_track have different meanings..
    Similarityf undo_track; 
	
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
