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


#include <wrap/gui/trackmode.h>
#include <wrap/gui/trackball.h>
#include <vcg/space/intersection3.h>
#include <vcg/math/similarity.h>
#include <iostream>

using namespace std;
using namespace vcg;


void TrackMode::Apply(Trackball *trackball, float WheelNotch) {
    trackball->track.sca*=pow(1.2f,WheelNotch);
 }

void ScaleMode::Apply(Trackball *tb, Point3f new_point) {
  float ScreenHeight= tb->camera.viewport[3]-tb->camera.viewport[1];
  float dist=(new_point[1]-tb->last_point[1])/ScreenHeight;
  tb->track.sca= tb->last_track.sca*pow(3.0f,-dist);
}

/// Compute the plane plane perpedicular to view dir and passing through manip center
Plane3f TrackMode::GetViewPlane(const View<float> &camera, const Point3f &center) {
  Point3f vp = camera.ViewPoint();
  
  Plane3f pl; 
  Point3f plnorm= vp - center;
  plnorm.Normalize();
  pl.Set(plnorm, plnorm*center);
  return pl;
}

/// Given a point p in window coordinate it compute the point where the lie p 
/// over the plane paralell the viewplane and passing through the center of the trackball

Point3f TrackMode::HitViewPlane(Trackball *tb, const Point3f &p) {
  // plane perpedicular to view direction and passing through manip center
  Plane3f vp = GetViewPlane(tb->camera, tb->center); 

  Line3fN ln= tb->camera.ViewLineFromWindow(Point3f(p[0],p[1],0));

  Point3f PonVP;
  bool res = Intersection<float>(vp,ln,PonVP);
  return PonVP;
}

// the most important function; given a new point in window coord, it update the transformation computed by the trackball.
// General scheme : the transformation is a function of just the begin and current mouse positions, with greater precision is function of just two 3d points over the manipulator.
void SphereMode::Apply(Trackball *tb, Point3f new_point) {
  Point3f hitOld=Hit(tb, tb->last_point);
  Point3f hitNew=Hit(tb, new_point);
  tb->Hits.push_back(hitNew);

  Point3f axis = (hitNew- tb->center)^(hitOld- tb->center); 

  //  Figure out how much to rotate around that axis.
  //float phi=Angle((hitNew- tb->center),(hitOld- tb->center));
  float phi = Distance(hitNew,hitOld) / tb->radius;
    
  tb->track.rot = tb->last_track.rot * Quaternionf(phi,axis);
 /* Codice Originale Ponchio 

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

  tb->track = Similarityf().SetRotate(diff) * tb->last_track;*/

}

/*
dato un punto in coordinate di schermo e.g. in pixel stile opengl 
calcola il punto di intersezione tra la viewline  che passa per viewpoint e per hitplane e l'iperboloide.
l'iperboloide si assume essere quello di rotazione attorno alla retta viewpoint-center e di raggio rad
si assume come sistema di riferimento quello con l'origine su center ecome x la retta center-viewpoint

eq linea
       hitplane.y
y = - ----------- * x + hitplane.y 
      viewpoint.x

eq hiperboloide di raggio r (e.g. che passa per (r/sqrt2,r/sqrt2) 

     1
y = --- * (r^2 /2.0)
     x 

 hitplane.y
 ----------- * x^2 - hitplane.y *x + (r^2/2.0) == 0
 viewpoint.x

*/
bool SphereMode::HitHyper(Point3f center,  float radius, Point3f viewpoint, Plane3f vp, Point3f hitplane, Point3f &hit) 
{
  float hitplaney = Distance(center,hitplane);
  float viewpointx= Distance(center,viewpoint);

  float a =  hitplaney/viewpointx;
  float b = -hitplaney;
  float c = radius*radius/2.0f;
  float delta = b*b - 4*a*c;
  float x1,x2,xval,yval;
  if(delta>0)
  {
    x1= (- b - sqrt(delta))/(2.0f*a);
    x2= (- b + sqrt(delta))/(2.0f*a);
    
    xval=x1; // always take the minimum value solution
    yval=c/xval;     //  alternatively it also oould be the other part of the equation yval=-(hitplaney/viewpointx)*xval+hitplaney;
  }
  else 
  {
    return false;
  }
  // Computing the result in 3d space;
  Point3f dirRadial=hitplane-center;
  dirRadial.Normalize();
  Point3f dirView=vp.Direction();
  dirView.Normalize();
  hit= center +dirRadial*yval+dirView*xval;
  return true;
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
   eq iperboloide :       y=1/2*x;         inf  sqrt(1/2)  1/2
   eq cono                y=x+sqrt(2);

   */

Point3f SphereMode::Hit(Trackball *tb, const Point3f &p) {
  const float Thr = tb->radius/math::Sqrt(2.0f);

  Line3fN vn = tb->camera.ViewLineFromModel(tb->center);
  Line3fN ln = tb->camera.ViewLineFromWindow(Point3f(p[0],p[1],0));
  Point3f viewpoint = tb->camera.ViewPoint();
  Plane3f vp = GetViewPlane(tb->camera, tb->center); 
  Point3f hit,hitPlane,hitSphere,hitSphere1,hitSphere2,hitHyper;
  Intersection<float>(vp, ln, hitPlane);
  Sphere3f sphere(tb->center,tb->radius);
  bool resSp = Intersection<float>(sphere, ln, hitSphere1, hitSphere2);
  if(Distance(viewpoint,hitSphere1)<Distance(viewpoint,hitSphere2)) 
        hitSphere=hitSphere1;
  else  hitSphere=hitSphere2;

  float dl=Distance(ln,tb->center);
  bool resHp = HitHyper(tb->center, tb->radius, viewpoint, vp, hitPlane, hitHyper) ;

  // four cases 
  
  // 1) Degenerate line tangent to both sphere and hyperboloid! 
  if((!resSp && !resHp) )
  {
    hit=ClosestPoint(ln,tb->center);
    //printf("closest point to line %f\n",Distance(hit,tb->center));
    return hit;
  }
  if((resSp && !resHp) ) return hitSphere; // 2) line cross only the sphere
  if((!resSp && resHp) ) return hitHyper;  // 3) line cross only the hyperboloid

  // 4) line cross both sphere and hyperboloid: choose according angle.
  float angleDeg=math::ToDeg(Angle((viewpoint-tb->center),(hitSphere-tb->center)));
  //printf("Angle %f (%5.2f %5.2f %5.2f) (%5.2f %5.2f %5.2f)\n",angleDeg,hitSphere[0],hitSphere[1],hitSphere[2],hitHyper[0],hitHyper[1],hitHyper[2]);
  if(angleDeg<45) return hitSphere;
          else    return hitHyper;

  // 
  // Codice ORIGINALE PONCHIO
  //vp.SetOffset(vp.Offset() + Thr);

  //Point3f hit;
  //bool res = Intersection<float>(vp, ln, hit);
  //float d = Distance(tb->center - vn.Direction()*Thr, hit);
  //if(d < Thr) {
  //  Point3f hit2;
  //  Sphere3f sphere(tb->center, tb->radius);
  //  bool res = Intersection<float>(sphere, ln, hit, hit2);

  //  //find closest intersection to sphere
  //  float d = (hit - viewpoint).Norm();
  //  float d2 = (hit2 - viewpoint).Norm();
  //  if(d > d2) hit = hit2;
  //  hit -= tb->center;
  //} else {
  //  if(d > 2.99 * Thr) 
  //    d = 2.99 * Thr;
  //  Point3f norm = (hit - tb->center)^(viewpoint - tb->center);
  //  norm.Normalize();
  //  float phi = -M_PI/4 - 3*M_PI/8 *(d - Thr)/Thr;

  //  Quaternionf q(phi, norm);
  //  hit = q.Rotate((viewpoint - tb->center).Normalize() * tb->radius);
  //}
  // hit.Normalize();
  // return hit;
}

void PlaneMode::Apply(Trackball *tb, Point3f new_point) {
  Point3f hitOld=HitViewPlane(tb, tb->last_point);
  Point3f hitNew=HitViewPlane(tb, new_point);

}
