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
Revision 1.3  2004/04/07 10:54:11  cignoni
Commented out unused parameter names and other minor warning related issues

Revision 1.2  2004/03/25 14:55:25  ponchio
Adding copyright.


****************************************************************************/

#ifndef TRACKMODE_H
#define TRACKMODE_H

#include <vcg/space/point3.h>
#include <vcg/math/similarity.h>
#include <wrap/gui/trackball.h>

namespace vcg {
class Trackball;
class TrackMode {
public:
  virtual ~TrackMode() {}
  //virtual void Draw() {}
  virtual Similarityf ComputeFromWindow(const Point3f &/* oldp */, const Point3f &/* newp */) { return Similarityf().SetIdentity(); }
  Point3f Hit(const Point3f &p);
  Trackball *tb;
  
};

class SphereMode: public TrackMode {
public:  
  Similarityf ComputeFromWindow(const Point3f &oldP, const Point3f &newP);
  //Plane3f SetViewPlane();
  Point3f Hit(const Point3f &p);
  Plane3f GetViewPlane();
//  Line3f GetViewLine(const Point3f &p);

};

class GravityMode: public TrackMode {
public:
};

class CylinderMode: public TrackMode {
public:
  CylinderMode(const Point3f _axis): axis(_axis) { axis.Normalize(); }
protected:
  Point3f axis;
};

class PlaneMode: public TrackMode {
public:
  PlaneMode(const Point3f _x, const Point3f _y): x(_x), y(_y) { x.Normalize(); y.Normalize(); }
  Similarityf ComputeFromWindow(const Point3f &oldP, const Point3f &newP);
protected:
  Point3f x;
  Point3f y;
};

class LineMode: public TrackMode {
public:
  LineMode(const Point3f _axis): axis(_axis) { axis.Normalize();}
protected:
  Point3f axis;
};

class ScaleMode: public TrackMode {
public:
};

}//namespace 

#endif