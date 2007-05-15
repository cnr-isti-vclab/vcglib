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
Revision 1.10  2007/02/26 01:30:02  cignoni
Added reflection Name

Revision 1.9  2006/02/13 13:10:27  cignoni
Added Zmode for moving objects along the perpendicular to the viewplane

Revision 1.8  2004/07/18 06:54:08  cignoni
Added Scaling

Revision 1.7  2004/07/11 22:06:56  cignoni
Added scaling by wheel

Revision 1.6  2004/06/09 14:01:13  cignoni
Heavily restructured. To be completed only rotation works...

Revision 1.5  2004/05/14 03:15:09  ponchio
Redesigned partial version.

Revision 1.4  2004/05/07 12:46:08  cignoni
Restructured and adapted in a better way to opengl

Revision 1.3  2004/04/07 10:54:11  cignoni
Commented out unused parameter names and other minor warning related issues

Revision 1.2  2004/03/25 14:55:25  ponchio
Adding copyright.


****************************************************************************/

#ifndef TRACKMODE_H
#define TRACKMODE_H

#include <vcg/space/line3.h>
#include <vcg/space/plane3.h>
#include <vcg/space/segment3.h>
#include <vcg/space/ray3.h>
#include <wrap/gui/view.h>

using namespace std;

namespace vcg {
  
class Trackball;

// Base class for all the track modes.
// This class' functions does nothing.
class TrackMode {
public:
  virtual ~TrackMode () {
  } 
  virtual void Apply (Trackball * trackball, Point3f new_point);
  virtual void Apply (Trackball * trackball, float WheelNotch);
  virtual void SetAction ();
  virtual void Reset ();  
  virtual const char *Name (){
    return "TrackMode";
  };
  virtual void Draw (Trackball * trackball);
}; 

// Inactive mode.
// useful only for drawing the inactive trackball
class InactiveMode:public TrackMode {
public:
  const char *Name () {
    return "InactiveMode";
  };
  void Draw (Trackball * trackball);
};

/* View space modes */
 
// old interfaces
/*
class SphereMode: public TrackMode {
}
class CylinderMode: public TrackMode {
}
class PlaneMode: public TrackMode {
}
class ZMode: public TrackMode {
}
class LineMode: public TrackMode {
}
class LineMode: public TrackMode {
}
class ScaleMode: public TrackMode {

*/

// Sphere mode.
// The classic trackball.
class SphereMode:public TrackMode {
public:
  void Apply (Trackball * trackball, Point3f new_point);
  const char *Name () {
    return "SphereMode";
  };
  void Draw (Trackball * trackball);
};

// Panning mode.
// The user can drag the model on the view plane.
class PanMode:public TrackMode {
public:
  void Apply (Trackball * trackball, Point3f new_point);
  const char *Name () {
    return "PanMode";
  };
  void Draw (Trackball * trackball);
};

// Z mode.
// Dragging the mouse up and down or scrolling the 
// mouse wheel will move the object along the Z of the camera.
class ZMode:public TrackMode {
public:
  const char *Name () {
    return "ZMode";
  };
  void Apply (Trackball * trackball, Point3f new_point);
  void Apply (Trackball * trackball, float WheelNotch);
  void Draw (Trackball * trackball);
};

// Scale Mode.
// Dragging the mouse up and down or scrolling the 
// mouse wheel will scale the object.
class ScaleMode:public TrackMode {
public:
  const char *Name () {
    return "ScaleMode";
  };
  void Apply (Trackball * trackball, Point3f new_point);
  void Apply (Trackball * trackball, float WheelNotch);
  void Draw (Trackball * trackball);
};

// Axis mode.
// Moves the object in a costrained direction.
// The user can either drag the mouse or scroll the wheel.
// The direction can be specified either with a line
// or a origin and a direction.
// the object posistion is not needed to be on the line.
class AxisMode:public TrackMode {
public:
  AxisMode (const Line3f & ln)
    : axis (ln) {
  } 
  AxisMode (const Point3f & origin, const Point3f & direction) {
    axis = Line3fN (origin, direction);
  }
  const char *Name () {
    return "AxisMode";
  };
  void Apply (Trackball * trackball, Point3f new_point);
  void Apply (Trackball * trackball, float WheelNotch);
  void Draw (Trackball * trackball);
private:
    Line3fN axis;
};

// Plane mode.
// The user can drag the object in a costrained plane.
// The plane can be specified either with a plane
// or the plane's equation parameters
// the object posistion is not needed to be on the plane.
class PlaneMode:public TrackMode {
public:
  PlaneMode (float a, float b, float c, float d)
    : plane(Plane3f(d,Point3f(a,b,c))){
  }  
  PlaneMode (Plane3f & pl)
    : plane(pl) {
  }
  const char *Name () {
    return "PlaneMode";
  };
  void Apply (Trackball * trackball, Point3f new_point);
  void Draw (Trackball * trackball);
private:
  Plane3f plane;
};

// Cylinder mode.
// Rotates the object along a fixed axis
// The user can either drag the mouse or scroll the wheel,
// in either cases the rotation's angle is influenced by 
// the radius of the trackball.
// The axis can be specified either with a line
// or a origin and a direction 
// when the user drags the mouse, if the axis is too
// perpendicular to view plane, the angle is specified
// only by the vertical component of the mouse drag and the radius.
class CylinderMode:public TrackMode {
public:
  CylinderMode (Line3fN & ln)
    : axis (ln){
  }
  CylinderMode (const Point3f & origin, const Point3f & direction)
    : axis (Line3fN(origin,direction)){  
  }
  const char *Name () {
    return "CylinderMode";
  };
  void Apply (Trackball * trackball, Point3f new_point);
  void Apply (Trackball * trackball, float WheelNotch);
  void Draw (Trackball * trackball);
private:
  Line3fN axis;
};

// Path mode.
// move the object along an eventually closed path.
// The user can either drag the mouse or scroll the wheel,
// when the user drags the mouse, the object tries to slide toward it.
// if the path is a simple segment, it can be specified just with the endpoints,
// otherwise it's specified with a point vector and, eventually, a boolean value used for closing the path.
// the object is assumed to initially be on the same position of the first point on the path.
// you can try to set the starting point calling SetStartNear(Point3f)
// the path is NOT assumed to have 0-length segments, so, if you want to close the path, please DO NOT add
// a copy of the first point on the end of the vector...
// the vector passed to build the path is copied locally.
class PathMode:public TrackMode {
public:
  PathMode ( const vector < Point3f > &pts, bool w = false)
    : points(), wrap(w), current_state(0), initial_state(0), old_hitpoint()
  {
    Init(pts);
	assert(min_seg_length > 0.0f);
  }
  PathMode ( const Point3f &start, const Point3f &end )
    : points(), wrap(false), current_state(0), initial_state(0), old_hitpoint()
  {
	points.push_back(start); 	
	points.push_back(end); 
	path_length=Distance(start,end);
	min_seg_length=path_length;
	assert(min_seg_length > 0.0f);
  }  
  const char *Name () {
    return "PathMode";
  };
  void Apply (Trackball * trackball, Point3f new_point);
  void Apply (Trackball * trackball, float WheelNotch);
  void Draw (Trackball * trackball);
  void SetAction ();
  void Reset (); 
  Point3f SetStartNear(Point3f p);
private:
  void Init(const vector < Point3f > &points);
  void GetPoints(float state, Point3f & point, Point3f & prev_point, Point3f & next_point);
  float Normalize(float state);
  float HitPoint(float state, Ray3fN ray, Point3f &hit_point);
  int Verse(Point3f reference_point,Point3f current_point,Point3f prev_point,Point3f next_point);

  vector < Point3f > points;
  bool wrap;
  float current_state;
  float initial_state;
  float path_length;
  float min_seg_length;
  Point3f old_hitpoint;

};

// Area mode.
// The user can drag the object inside a planar area, defined by a polygon.
// the polygon can be non convex, and is specified with a vector of vertexes
// if the object's trajectory intersects some poligon side, it tries to slide 
// around it, in a "rubber band flavoured" way.
// for the vertexes vector its calculated the plane of the polygon, and then
// every point in the vector is projected on this plane.
// the object is assumed to initially be on the same position of the first vertex.
// you can try to set the starting point calling SetStartNear(Point3f)
// the vector is assumed to be formed of NON collinear points
// the polygon is NOT assumed to have 0-length sides, so please DO NOT add
// a copy of the first point on the end of the vector...
// the vector passed to build the polygon is copied locally.

class AreaMode:public TrackMode {
public:
  AreaMode (const vector < Point3f > &pts)
  {
    Init(pts);
    assert(min_side_length > 0.0f);
  }
  const char *Name () {
    return "AreaMode";
  };
  void Apply (Trackball * trackball, Point3f new_point);
  void Draw (Trackball * trackball);
  void SetAction ();
  void Reset (); 
  Point3f SetStartNear(Point3f p);
private:
  void Init(const vector < Point3f > &pts);
  bool Inside(Point3f point);
  Point3f Move(Point3f start,Point3f end);

  vector < Point3f > points;
  bool begin_action;
  int first_coord_kept;
  int second_coord_kept;
  float min_side_length;
  Point3f status,delta_mouse,old_status,initial_status;
  
  Plane3f plane;   
  Point3f rubberband_handle ;
  vector < Point3f > path;
  
};

}//namespace 

#endif
