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
Revision 1.4  2004/05/07 12:46:08  cignoni
Restructured and adapted in a better way to opengl

Revision 1.3  2004/04/07 10:54:11  cignoni
Commented out unused parameter names and other minor warning related issues

Revision 1.2  2004/03/25 14:55:25  ponchio
Adding copyright.


****************************************************************************/


#include <wrap/gui/trackmode.h>
#include <wrap/gui/trackball.h>
#include <vcg/space/intersection3.h>
#include <vcg/math/similarity.h>
#include <iostream>
using namespace std;

using namespace vcg;

Plane3f TrackMode::GetViewPlane(const View<float> &camera, Point3f center) {
  Point3f vp = camera.ViewPoint();
  
  Plane3f pl; //plane perpedicular to view dir and passing through manip center
  pl.Set(vp - center, (vp - center)*center);
  return pl;
}

Point3f TrackMode::HitViewPlane(Trackball *tb, const Point3f &p) {
  // plane perpedicular to view direction and passing through manip center
  Plane3f vp = GetViewPlane(tb->camera, tb->center); 

  Line3fN ln= tb->camera.ViewLineFromWindow(Point3f(p[0],p[1],0));

  Point3f PonVP;
  bool res = Intersection<float>(vp,ln,PonVP);
  return PonVP;
}


void SphereMode::Apply(Trackball *tb, Point3f new_point) {
  Point3f hitOld=Hit(tb, tb->last_point);
  Point3f hitNew=Hit(tb, new_point);

  Point3f ref = (tb->camera.ViewPoint() - tb->center).Normalize();


  Point3f axis = hitNew^ref; 
  axis.Normalize();
  float dist = (hitNew - ref).Norm()/2;
  float phi = 2 * math::Asin(dist);

  Point3f oaxis = hitOld^ref; 
  oaxis.Normalize();
  float odist = (hitOld - ref).Norm()/2;
  float ophi = 2 * math::Asin(odist);

  
  Quaternionf r = tb->last_track.rot;
  Quaternionf diff = r * Quaternionf(phi, axis) * 
                     Quaternionf(-ophi, oaxis) * Inverse(r);

  tb->track = Similarityf().SetRotate(diff) * tb->last_track;
}

/* dato un punto in coordinate di schermo e.g. in pixel stile opengl 
   restituisce un punto in coordinate di mondo sulla superficie 
   della trackball.
   La superficie della trackball e' data da una sfera + una porzione 
   di iperboloide di rotazione.
   Assumiamo la sfera di raggio unitario e centrata sull'origine e 
   di guardare lungo la y negativa.

                                       X   0   sqrt(1/2)  1  
   eq sfera:              y=sqrt(1-x*x);   1   sqrt(1/2)  0   
   eq iperboloide :       y=1/2x;         inf  sqrt(1/2)  1/2
   eq cono                y=x+sqrt(2);

   */

Point3f SphereMode::Hit(Trackball *tb, const Point3f &p) {
  const float Thr = tb->radius/math::Sqrt(2.0f);
  Line3fN vn = tb->camera.ViewLineFromModel(tb->center);
  Line3fN ln = tb->camera.ViewLineFromWindow(Point3f(p[0],p[1],0));
  Point3f viewpoint = tb->camera.ViewPoint();

  Plane3f vp = GetViewPlane(tb->camera, tb->center); 
  vp.SetOffset(vp.Offset() + Thr);

  Point3f hit;
  bool res = Intersection<float>(vp, ln, hit);
  float d = Distance(tb->center - vn.Direction()*Thr, hit);
  if(d < Thr) {
    Point3f hit2;
    Sphere3f sphere(tb->center, tb->radius);
    bool res = Intersection<float>(sphere, ln, hit, hit2);

    //find closest intersection to sphere
    float d = (hit - viewpoint).Norm();
    float d2 = (hit2 - viewpoint).Norm();
    if(d > d2) hit = hit2;
    hit -= tb->center;
  } else {
    if(d > 2.99 * Thr) 
      d = 2.99 * Thr;
    Point3f norm = (hit - tb->center)^(viewpoint - tb->center);
    norm.Normalize();
    float phi = -M_PI/4 - 3*M_PI/8 *(d - Thr)/Thr;

    Quaternionf q(phi, norm);
    hit = q.Rotate((viewpoint - tb->center).Normalize() * tb->radius);
  }
  hit.Normalize();
  return hit;
}

void PlaneMode::Apply(Trackball *tb, Point3f new_point) {
  Point3f hitOld=HitViewPlane(tb, tb->last_point);
  Point3f hitNew=HitViewPlane(tb, new_point);

}
