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
#include <wrap/gui/view.h>

namespace vcg {
  
class Trackball;

class TrackMode {
public:
  virtual ~TrackMode() {}
  virtual void Apply(Trackball *trackball, Point3f new_point) = 0;
  virtual void Apply(Trackball *trackball, float WheelNotch);

  virtual void Draw() {}
 protected:
  Plane3f GetViewPlane(const View<float> &view, const Point3f &center);
  Point3f HitViewPlane(Trackball *trackball, const Point3f &p);
};
 
/* View space modes */
 
class SphereMode: public TrackMode {
 public:  
  void Apply(Trackball *trackball, Point3f new_point);
 protected:
  Point3f Hit(Trackball *trackball, const Point3f &p);
  bool HitHyper(Point3f center,  float radius, Point3f viewpoint, Plane3f vp, Point3f hitplane, Point3f &hit) ;

};

class CylinderMode: public TrackMode {
public:
  CylinderMode(const Line3f &/*line*/, float /*radius = 1*/) {}
  void Apply(Trackball * /*trackball*/, Point3f /*new_point*/) {}
protected:
  Line3f line;
  float radius;
};

class PlaneMode: public TrackMode {
public:
  PlaneMode(const Plane3f &pl): plane(pl) {}
  void Apply(Trackball *trackball, Point3f new_point);
protected:
  Plane3f plane;
};

class LineMode: public TrackMode {
public:
  LineMode(const Line3f &/*line*/) {}
  void Apply(Trackball * /*trackball*/, Point3f /*new_point*/) {}
protected:
  Line3f line;
};

class ScaleMode: public TrackMode {
public:
  void Apply(Trackball * /*trackball*/, Point3f /*new_point*/) {}
};

}//namespace 

#endif
