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


****************************************************************************/

#ifndef TRACKUTILS_H
#define TRACKUTILS_H

#include <assert.h>
#include <vcg/math/base.h>
#include <vcg/math/similarity.h>
#include <vcg/space/intersection3.h>
#include <vcg/space/line3.h>
#include <vcg/space/plane3.h>
#include <wrap/gl/math.h>
#include <wrap/gl/space.h>
#include <vector>

using namespace std;

namespace vcg {
	
namespace trackutils {

/// Compute the plane perpedicular to view dir and passing through manip center
Plane3f GetViewPlane (const View < float >&camera, const Point3f & center)
{
  Point3f vp = camera.ViewPoint ();
  Plane3f pl;
  Point3f plnorm = vp - center;
  plnorm.Normalize ();
  pl.Set (plnorm, plnorm * center);
  return pl;
}

Ray3f line2ray(const Line3f &l){
  Ray3f r(l.Origin(),l.Direction());
  r.Normalize();
  return r;
}

// Given a point p in window coordinate it compute the point where the lie p 
// over the plane paralell the viewplane and passing through the center of the trackball
Point3f HitViewPlane (Trackball * tb, const Point3f & p)
{
  // plane perpedicular to view direction and passing through manip center
  Plane3f vp = GetViewPlane (tb->camera, tb->center);
  Line3fN ln = tb->camera.ViewLineFromWindow (Point3f (p[0], p[1], 0));
  Point3f PonVP;
  /*bool res = */ IntersectionLinePlane < float >(vp, ln, PonVP);
  return PonVP;
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
bool HitHyper (Point3f center, float radius, Point3f viewpoint, Plane3f vp,
                           Point3f hitplane, Point3f & hit)
{
  float hitplaney = Distance (center, hitplane);
  float viewpointx = Distance (center, viewpoint);

  float a = hitplaney / viewpointx;
  float b = -hitplaney;
  float c = radius * radius / 2.0f;
  float delta = b * b - 4 * a * c;
  float x1, x2, xval, yval;

  if (delta > 0) {
    x1 = (-b - sqrt (delta)) / (2.0f * a);
    x2 = (-b + sqrt (delta)) / (2.0f * a);

    xval = x1;                  // always take the minimum value solution
    yval = c / xval;            //  alternatively it also could be the other part of the equation yval=-(hitplaney/viewpointx)*xval+hitplaney;
  }
  else {
    return false;
  }
  // Computing the result in 3d space;
  Point3f dirRadial = hitplane - center;
  dirRadial.Normalize ();
  Point3f dirView = vp.Direction ();
  dirView.Normalize ();
  hit = center + dirRadial * yval + dirView * xval;
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

Point3f HitSphere (Trackball * tb, const Point3f & p)
{
  // const float Thr = tb->radius/math::Sqrt(2.0f);
  Line3fN vn = tb->camera.ViewLineFromModel (tb->center);
  Line3fN ln = tb->camera.ViewLineFromWindow (Point3f (p[0], p[1], 0));
  Point3f viewpoint = tb->camera.ViewPoint ();
  Plane3f vp = GetViewPlane (tb->camera, tb->center);
  Point3f hit, hitPlane, hitSphere, hitSphere1, hitSphere2, hitHyper;
  IntersectionLinePlane < float >(vp, ln, hitPlane);

  Sphere3f sphere (tb->center, tb->radius);
  bool resSp = IntersectionLineSphere < float >(sphere, ln, hitSphere1, hitSphere2);

  if (resSp == true) {
    if (Distance (viewpoint, hitSphere1) < Distance (viewpoint, hitSphere2))
      hitSphere = hitSphere1;
    else
      hitSphere = hitSphere2;
  }

  /*float dl= */ Distance (ln, tb->center);
  bool resHp = HitHyper (tb->center, tb->radius, viewpoint, vp, hitPlane, hitHyper);

  // four cases

  // 1) Degenerate line tangent to both sphere and hyperboloid!
  if ((!resSp && !resHp)) {
    hit = ClosestPoint (ln, tb->center);
    //printf("closest point to line %f\n",Distance(hit,tb->center));
    return hit;
  }
  if ((resSp && !resHp))
    return hitSphere;           // 2) line cross only the sphere
  if ((!resSp && resHp))
    return hitHyper;            // 3) line cross only the hyperboloid

  // 4) line cross both sphere and hyperboloid: choose according angle.
  float angleDeg = math::ToDeg (Angle ((viewpoint - tb->center), (hitSphere - tb->center)));

  //printf("Angle %f (%5.2f %5.2f %5.2f) (%5.2f %5.2f %5.2f)\n",angleDeg,hitSphere[0],hitSphere[1],hitSphere[2],hitHyper[0],hitHyper[1],hitHyper[2]);
  if (angleDeg < 45)
    return hitSphere;
  else
    return hitHyper;

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

/*
 calculates the minimal distance between 2 lines.
 P and Q are the lines, P_s and Q_t are set to be the closest points on these lines.
 it's returned the distance from P_s and Q_t, and a boolean value which is true
 if the lines are parallel enough.
 if P and Q are parallel P_s and Q_t aren't set.
 the formula is taken from pages 81-83 of
 "Eric Lengyel - Mathematics for 3D Game Programming & Computer Graphics"
*/
pair< float, bool > LineLineDistance(const Line3f & P,const Line3f & Q,Point3f & P_s, Point3f & Q_t){
  Point3f p0 = P.Origin (), Vp = P.Direction ();
  Point3f q0 = Q.Origin (), Vq = Q.Direction ();
  float VPVP = Vp * Vp;
  float VQVQ = Vq * Vq;
  float VPVQ = Vp * Vq;
  const float det = ( VPVP * VQVQ ) - ( VPVQ * VPVQ );
  const float EPSILON = 0.00001f;
  if ( fabs(det) < EPSILON ) {
  	return make_pair(Distance(P,q0), true);
  }
  float b1= (q0 - p0) * Vp;
  float b2= (p0 - q0) * Vq;
  float s = ( (VQVQ * b1) + (VPVQ * b2) ) / det;
  float t = ( (VPVQ * b1) + (VPVP * b2) ) / det;
  P_s = p0 + (Vp * s);
  Q_t = q0 + (Vq * t);
  return make_pair(Distance(P_s,Q_t),false);
}

/*
 calculates the minimal distance between a ray and a line
 R is the ray and Q is the line, R_s and Q_t are set to be the closest points on 
 the ray and the line.
 it's returned the distance from R_s and Q_t, and a boolean value which is true
 if the ray and the line are parallel enough.
 if R and Q are parallel R_s and Q_t aren't set.
*/
pair< float, bool > RayLineDistance(const Ray3f & R,const Line3f & Q,Point3f & R_s, Point3f & Q_t){
  Point3f r0 = R.Origin (), Vr = R.Direction ();
  Point3f q0 = Q.Origin (), Vq = Q.Direction ();
  float VRVR = Vr * Vr;
  float VQVQ = Vq * Vq;
  float VRVQ = Vr * Vq;
  const float det = ( VRVR * VQVQ ) - ( VRVQ * VRVQ );
  const float EPSILON = 0.00001f;
  if ( ( det >= 0.0f ? det : -det) < EPSILON ) {
  	return make_pair(Distance(Q,r0), true);
  }
  float b1= (q0 - r0) * Vr;
  float b2= (r0 - q0) * Vq;
  float s = ( (VQVQ * b1) + (VRVQ * b2) ) / det;
  float t = ( (VRVQ * b1) + (VRVR * b2) ) / det;
  if(s<0){
    R_s = r0;
    Q_t = ClosestPoint(Q,R_s);    
  }else {
    R_s = r0 + (Vr * s);
    Q_t = q0 + (Vq * t);
  }
  return make_pair(Distance(R_s,Q_t),false);
}

/*
 calculates the minimal distance between 2 segments
 R e Q are the segments, R_s and Q_t are set to be the closest points on 
 the segments
 it's returned the distance from R_s and Q_t, and a boolean value which is true
 if the segments are parallel enough.
*/
pair< float, bool > SegmentSegmentDistance(const Segment3f & R, const Segment3f & Q, Point3f & R_s, Point3f & Q_t)
{
  float R_len=Distance(R.P0(),R.P1());
  float Q_len=Distance(Q.P0(),Q.P1());
  const float EPSILON_LENGTH = max(R_len,Q_len)*0.0001f; 
  if(R_len < EPSILON_LENGTH){
  	R_s=R.P0();
  	Q_t=ClosestPoint(Q,R_s);
  	return make_pair(Distance(R_s,Q_t),true);
  }
  if( Q_len < EPSILON_LENGTH){
  	Q_t=Q.P0();
  	R_s=ClosestPoint(R,Q_t);
  	return make_pair(Distance(R_s,Q_t),true);
  }  
  Point3f r0 = R.P0(), Vr = (R.P1()-R.P0()).Normalize();
  Point3f q0 = Q.P0(), Vq = (Q.P1()-Q.P0()).Normalize();
  float VRVR = Vr * Vr;
  float VQVQ = Vq * Vq;
  float VRVQ = Vr * Vq;
  const float det = ( VRVR * VQVQ ) - ( VRVQ * VRVQ );
  const float EPSILON = 0.00001f;
  if ( ( det >= 0.0f ? det : -det) < EPSILON ) {
  	Line3f lR(R.P0(),R.P1());
  	float qa=lR.Projection(Q.P0());
  	float qb=lR.Projection(Q.P1());  	
  	if( (qa<=0.0f) && qb<=(0.0f)){
      R_s=R.P0();
      Q_t=ClosestPoint(Q,R_s);
  	} else if ( (qa >= 1.0f) && (qb >= 1.0f) ){
      R_s=R.P1();
      Q_t=ClosestPoint(Q,R_s);
  	} else {
      if( (qa >= 0.0f) && (qa <= 1.0f) ){
        Q_t=Q.P0();
		R_s=ClosestPoint(R,Q_t);
	  } else if((qb >= 0.0f) && (qb <= 1.0f) ){
        Q_t=Q.P1();
        R_s=ClosestPoint(R,Q_t);
      } else {
        if( ((qa<=0.0f)&&(qb>=1.0f) ||((qb<=0.0f)&&(qa>=1.0f)))){
           R_s=R.P0();
           Q_t=ClosestPoint(Q,R_s);
        }else{
           assert(0);
        }
      }
  	}  	
	return make_pair(Distance(R_s,Q_t),true);
  }
  float b1= (q0 - r0) * Vr;
  float b2= (r0 - q0) * Vq;
  float s = ( (VQVQ * b1) + (VRVQ * b2) ) / det;
  float t = ( (VRVQ * b1) + (VRVR * b2) ) / det;
  if( s < 0 ){
    R_s = R.P0();
  }else if ( s > R_len ){
    R_s = R.P1();
  } else {
    R_s = r0 + (Vr * s);
  }
  if( t < 0){
  	Q_t = Q.P0(); 
  }else if ( t > Q_len ){
    Q_t = Q.P1();
  }else{
    Q_t = q0 + (Vq * t);
  }
  return make_pair(Distance(R_s,Q_t),false);
}


pair< Point3f,bool > HitNearestPointOnAxis (Trackball * tb,Line3f axis, Point3f point)
{
  Ray3fN ray = line2ray(tb->camera.ViewLineFromWindow (point));
  Point3f axis_p(0,0,0), ray_p(0,0,0);  
  pair< float, bool > resp=RayLineDistance(ray,axis,ray_p,axis_p);
  if(resp.second || (ray_p == ray.Origin())){
  	return make_pair(Point3f(0,0,0),false);
  }
  return make_pair(axis_p,true);
}


Line3f ProjectLineOnPlane(const Line3f & ln, const Plane3f & pl)
{
  Point3f l0=ln.Origin();
  Point3f l1=l0+ln.Direction();
  Point3f p1,p2;
  p1=pl.Projection(l0);
  p2=pl.Projection(l1);
  Line3f res(p1,p2-p1);
  return res;
}

float signedDistance(Line3f line,Point3f pt,Point3f positive_dir)
{
  return Distance(line,pt) * ((((pt-ClosestPoint(line,pt)) * positive_dir) >= 0.0 )? 1.0: -1.0);
}

float getDeltaY(Trackball * tb, Point3f new_point)
{
  float ScreenHeight = float (tb->camera.viewport[3] - tb->camera.viewport[1]);
  return (new_point[1] - tb->last_point[1]) / ScreenHeight;
}

/// intersection between RAY and plane
template<class T>
  inline bool IntersectionRayPlane( const Plane3<T> & pl, const Ray3<T> & ray, Point3<T> &po){
  const T epsilon = T(1e-8);

  T k = pl.Direction() * ray.Direction(); // Compute 'k' factor
  if( (k > -epsilon) && (k < epsilon))
    return false;
  T r = (pl.Offset() - pl.Direction()*ray.Origin())/k;  // Compute ray distance
  if (r < 0)
    return false;
  po = ray.Origin() + ray.Direction()*r;
  return true;
}

pair< Point3f, bool > HitPlane (Trackball * tb, Point3f point, Plane3f plane)
{
  Ray3fN ray = line2ray(tb->camera.ViewLineFromWindow (point));
  Point3f p(0,0,0);
  bool res = IntersectionRayPlane < float >(plane, ray, p);
  return make_pair(p,res);
}


// drawing section

class DrawingHint {
public:
  DrawingHint () {
    CircleStep = 64;
    HideStill = false;
    DrawTrack = false;
    LineWidthStill = 0.5f;
    LineWidthMoving = 1.5f;
    color = Color4b::LightBlue;
  }
  int CircleStep;
  bool HideStill, DrawTrack;
  Color4b color;
  float LineWidthStill;
  float LineWidthMoving;
};

DrawingHint DH;

void DrawPlaneHandle ()
{
  float r = 1.0;
  float dr = r / 10.0f;

  glBegin (GL_LINE_STRIP);
  glVertex3f (+r + dr, +r, 0.0);
  glVertex3f (+r, +r + dr, 0.0);
  glVertex3f (+r - dr, +r, 0.0);
  glVertex3f (+r, +r - dr, 0.0);
  glVertex3f (+r + dr, +r, 0.0);
  glEnd ();
  glBegin (GL_LINE_STRIP);
  glVertex3f (-r + dr, -r, 0.0);
  glVertex3f (-r, -r + dr, 0.0);
  glVertex3f (-r - dr, -r, 0.0);
  glVertex3f (-r, -r - dr, 0.0);
  glVertex3f (-r + dr, -r, 0.0);
  glEnd ();
}

void DrawCircle ()
{
  int nside = DH.CircleStep;
  const double pi2 = 3.14159265 * 2.0;
  glBegin (GL_LINE_LOOP);
  for (double i = 0; i < nside; i++) {
    glNormal3d (cos (i * pi2 / nside), sin (i * pi2 / nside), 0.0);
    glVertex3d (cos (i * pi2 / nside), sin (i * pi2 / nside), 0.0);
  }
  glEnd ();
  DrawPlaneHandle ();
}

void DrawSphereIcon (Trackball * tb,bool active)
{  
  glPushMatrix ();
  glTranslate (tb->center);
  glMultMatrix (tb->track.InverseMatrix ());
  glScale (tb->radius);
  Matrix44f r;
  tb->track.rot.ToMatrix (r);
  glMultMatrix (r);
  
  glPushMatrix ();
  float amb[4] = { .3f, .3f, .3f, 1.0f };
  float col[4] = { .5f, .5f, .8f, 1.0f };
  glPushAttrib (GL_ENABLE_BIT | GL_LINE_BIT | GL_CURRENT_BIT | GL_LIGHTING_BIT);
  if (active)
    glLineWidth (DH.LineWidthMoving);
  else
    glLineWidth (DH.LineWidthStill);
  glEnable (GL_LIGHTING);
  glEnable (GL_LIGHT0);
  glEnable (GL_LINE_SMOOTH);
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glColor (DH.color);
  glMaterialfv (GL_FRONT_AND_BACK, GL_EMISSION, amb);
  glMaterialfv (GL_FRONT_AND_BACK, GL_DIFFUSE, col);
  glPushMatrix ();
    DrawCircle ();
    glPushMatrix ();
      glRotatef (90, 1, 0, 0);
      DrawCircle ();
      glRotatef (90, 0, 1, 0);
      DrawCircle ();
    glPopMatrix ();
  glPopMatrix ();
  glPopAttrib ();
  glPopMatrix ();
  glPopMatrix ();
  
}

// TEMPORARY drawing section
// Disclaimer: the following code is of VERY POOR quality
// feel free to delete and rewrite everything

void prepara_attrib()
{
  float amb[4] = { .3f, .3f, .3f, 1.0f };
  float col[4] = { .5f, .5f, .8f, 1.0f };
  glEnable (GL_LIGHTING);
  glEnable (GL_LIGHT0);
  glEnable (GL_LINE_SMOOTH);
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glMaterialfv (GL_FRONT_AND_BACK, GL_EMISSION, amb);
  glMaterialfv (GL_FRONT_AND_BACK, GL_DIFFUSE, col);
}

void DrawUglyLetter(Trackball * tb,vector<Point3f> ugly_letter)
{
  Point3f center=tb->camera.Project(tb->center);
  float offset=0;  
  offset=max(offset,Distance(center,tb->camera.Project(tb->center+(Point3f(1,0,0) * tb->radius))));
  offset=max(offset,Distance(center,tb->camera.Project(tb->center+(Point3f(0,1,0) * tb->radius))));
  offset=max(offset,Distance(center,tb->camera.Project(tb->center+(Point3f(0,0,1) * tb->radius))));
  glPushMatrix();
  glPushAttrib (GL_ALL_ATTRIB_BITS);
   // go to world coords
  glTranslate (tb->center);
  glMultMatrix (tb->track.InverseMatrix ());
  glTranslate (-tb->center);
  prepara_attrib();
  glColor3f(1,1,1);
  glLineWidth(4.0);
  
  glBegin(GL_LINE_STRIP);
    for(unsigned int i=0;i<ugly_letter.size();i++){
  	  glVertex(tb->camera.UnProject(center+(ugly_letter[i] * offset * 0.25)
  	           +Point3f(-offset,-offset,0)));
    }  
  glEnd();
  glPopAttrib ();
  glPopMatrix();

}

void DrawUglyPanMode(Trackball * tb)
{
  vector<Point3f> ugly_p;
  ugly_p.push_back(Point3f(-1,-1,0));
  ugly_p.push_back(Point3f(-1,1,0));
  ugly_p.push_back(Point3f(1,1,0));
  ugly_p.push_back(Point3f(1,0,0));  
  ugly_p.push_back(Point3f(-1,0,0)); 
  
  DrawUglyLetter(tb,ugly_p);
}

void DrawUglyZMode(Trackball * tb)
{
  vector<Point3f> ugly_z;
  ugly_z.push_back(Point3f(-1,1,0));
  ugly_z.push_back(Point3f(1,1,0));
  ugly_z.push_back(Point3f(-1,-1,0));
  ugly_z.push_back(Point3f(1,-1,0));
  DrawUglyLetter(tb,ugly_z);
}

void DrawUglyScaleMode(Trackball * tb)
{
  vector<Point3f> ugly_s;
  ugly_s.push_back(Point3f(1,1,0));
  ugly_s.push_back(Point3f(-1,1,0));
  ugly_s.push_back(Point3f(-1,0,0));
  ugly_s.push_back(Point3f(1,0,0));
  ugly_s.push_back(Point3f(1,-1,0));
  ugly_s.push_back(Point3f(-1,-1,0));
  DrawUglyLetter(tb,ugly_s);
}

void DrawUglyAxisMode(Trackball * tb,Line3f axis)
{
  glPushMatrix();
  glPushAttrib (GL_ALL_ATTRIB_BITS);
  // go to world coords
  glTranslate (tb->center);
  glMultMatrix (tb->track.InverseMatrix ());
  glTranslate (-tb->center);
  prepara_attrib();
  glColor3f(0.9,0.9,0.2);
  glLineWidth(2.0);
   glBegin(GL_LINES);
    glVertex(axis.Origin()+(axis.Direction()*100));
    glVertex(axis.Origin()-(axis.Direction()*100));
  glEnd();
  glPointSize(8.0);
  glColor3f(0.2,0.2,0.9);
  glBegin(GL_POINTS);
    glVertex(axis.Origin());
  glEnd();
  glPopAttrib ();
  glPopMatrix();
}

void DrawUglyPlaneMode(Trackball * tb,Plane3f plane)
{
  glPushMatrix();
  glPushAttrib (GL_ALL_ATTRIB_BITS);
  // go to world coords
  glTranslate (tb->center);
  glMultMatrix (tb->track.InverseMatrix ());
  glTranslate (-tb->center);
  prepara_attrib();
  Point3f p0,d1,d2,norm;
  norm=plane.Direction();
  p0=plane.Projection(Point3f(0,0,0));
  d1=Point3f(0,1,0);
  if(norm == d1 || norm == -d1)
    d1 = Point3f(1,0,0);
  d2=plane.Projection(d1);
  d1=(d2 - p0).Normalize();  
  d2=(d1 ^ norm).Normalize();  
  glLineWidth(3.0);
  glColor3f(0.2,0.2,0.9);
  glBegin(GL_LINES);
    glVertex(p0);
    glVertex(p0+norm);
  glEnd();
  glLineWidth(1.0);
  for(float i=0.5;i<100.0; i+=0.7){
    glBegin(GL_LINE_LOOP);
    for(int a=0;a<360;a+=10){
      float f0=i*cos((M_PI*a)/180);
      float f1=i*sin((M_PI*a)/180);
      glVertex(p0+(d1*f0)+(d2*f1));
    }
    glEnd();
  }  
  glColor3f(0.9,0.9,0.2);
  glPointSize(8.0);
  glBegin(GL_POINTS);
    glVertex(p0);
  glEnd();
  glColor3f(0.7,0.7,0);
  glPointSize(6.0);
  glBegin(GL_POINTS);
    glVertex(p0+norm);
  glEnd();
  glPopAttrib ();
  glPopMatrix();
}

void DrawUglyCylinderMode(Trackball * tb,Line3f axis)
{
  glPushMatrix();
  glPushAttrib (GL_ALL_ATTRIB_BITS);
  // go to world coords
  glTranslate (tb->center);
  glMultMatrix (tb->track.InverseMatrix ());
  glTranslate (-tb->center);
  prepara_attrib();
  Plane3f plane;
  plane.Init(axis.Origin(),axis.Direction());
  Point3f p0,d1,d2,norm;
  norm=plane.Direction();
  p0=plane.Projection(Point3f(0,0,0));
  d1=Point3f(0,1,0);
  if(norm == d1 || norm == -d1)
    d1 = Point3f(1,0,0);
  d2=plane.Projection(d1);
  d1=(d2 - p0).Normalize();  
  d2=(d1 ^ norm).Normalize();
  glLineWidth(1.0);
  glColor3f(0.2,0.2,0.9);
  for(int i=-100;i<100;i++){
    glBegin(GL_LINE_LOOP);
    for(int a=0;a<360;a+=10){
      float f0=(tb->radius)*cos((M_PI*a)/180);
      float f1=(tb->radius)*sin((M_PI*a)/180);
      glVertex(p0+(norm*i)+(d1*f0)+(d2*f1));
    }
    glEnd();
  }  
  glLineWidth(3.0);
  glColor3f(0.2,0.2,0.9);
  glBegin(GL_LINES);
     glVertex(axis.Origin());
     glVertex(axis.Origin()+(axis.Direction()*100));
  glEnd();
  glLineWidth(1.5);
  glColor3f(0.9,0.2,0.9);
  glBegin(GL_LINES);
    glVertex(axis.Origin());
    glVertex(axis.Origin()-(axis.Direction()*100));
  glEnd();
  glColor3f(0.9,0.9,0.2);
  glPointSize(8.0);
  glBegin(GL_POINTS);
    glVertex(axis.Origin());
  glEnd();
  glPopAttrib ();
  glPopMatrix();
}

void DrawUglyPathMode(Trackball * tb,const vector < Point3f > &points,
                      Point3f current_point,Point3f prev_point,
                      Point3f next_point,Point3f old_hitpoint,bool wrap)
{
  glPushMatrix();
  glPushAttrib (GL_ALL_ATTRIB_BITS);
  // go to world coords
  glTranslate (tb->center);
  glMultMatrix (tb->track.InverseMatrix ());
  glTranslate (-tb->center);
  prepara_attrib();
  glColor3f(0.9,0.9,0.2);
  glLineWidth(2.0);
  if(wrap)
    glBegin(GL_LINE_LOOP);
  else
    glBegin(GL_LINE_STRIP);
  for (vector < Point3f >::const_iterator i = points.begin (); i != points.end (); ++i){
    glVertex(*i);
  }
  glEnd();
  glColor3f(1,0,1);
  glPointSize(8.0);
  glBegin(GL_POINTS);
    glVertex(current_point);
  glEnd();
  glColor3f(0.6,0,0.6);
  glPointSize(7.0);
  glBegin(GL_POINTS);
    glVertex(old_hitpoint);
  glEnd();
  glColor3f(0.7,0.7,0.7);
  glPointSize(6.5);
  glBegin(GL_POINTS);
    glVertex(prev_point);
    glVertex(next_point);
  glEnd();
  glPopAttrib ();
  glPopMatrix();
}

void DrawUglyAreaMode(Trackball * tb,const vector < Point3f > &points,
                      Point3f status,Point3f old_status,Plane3f plane,
                      const vector < Point3f > &path,Point3f rubberband_handle)
{
  glPushMatrix();
  glPushAttrib (GL_ALL_ATTRIB_BITS);
  // go to world coords
  glTranslate (tb->center);
  glMultMatrix (tb->track.InverseMatrix ());
  glTranslate (-tb->center);
  prepara_attrib();
  glColor3f(0.9,0.9,0.2);
  glLineWidth(2.0);
  glBegin(GL_LINE_LOOP);
  for (vector < Point3f >::const_iterator i = points.begin (); i != points.end (); ++i){
    glVertex(*i);
  }
  glEnd();
  glColor3f(0.0,0.9,0.2);
  glLineWidth(1.2);
  glBegin(GL_LINE_STRIP);
  for (vector < Point3f >::const_iterator i = path.begin (); i != path.end (); ++i){
    glVertex(*i);
  }
  glEnd();
   glColor3f(1,0,1);
  glPointSize(8.0);
  glBegin(GL_POINTS);
    glVertex(status);
  glEnd();
  glColor3f(0.6,0,0.6);
  glPointSize(7.0);
  glBegin(GL_POINTS);
    glVertex(old_status);
  glEnd();
  glColor3f(0.6,0,0.0);
  glPointSize(6.0);
  glBegin(GL_POINTS);
    glVertex(rubberband_handle);
  glEnd();
  glLineWidth(1.0);
  glBegin(GL_LINES);
    glVertex(rubberband_handle);
    glVertex(status);
  glEnd();
  Point3f p0,d1,d2,norm;
  norm=plane.Direction();
  p0=plane.Projection(Point3f(0,0,0));
  d1=Point3f(0,1,0);
  if(norm == d1 || norm == -d1)
    d1 = Point3f(1,0,0);
  d2=plane.Projection(d1);
  d1=(d2 - p0).Normalize();  
  d2=(d1 ^ norm).Normalize();  
  glLineWidth(3.0);
  glColor3f(0.2,0.2,0.9);
  glBegin(GL_LINES);
    glVertex(p0);
    glVertex(p0+norm);
  glEnd();
  glLineWidth(0.1);
  for(float i=0.5;i<100.0; i+=0.7){
    glBegin(GL_LINE_LOOP);
    for(int a=0;a<360;a+=10){
      float f0=i*cos((M_PI*a)/180);
      float f1=i*sin((M_PI*a)/180);
      glVertex(p0+(d1*f0)+(d2*f1));
    }
    glEnd();
  }  
  
  glPopAttrib ();
  glPopMatrix();
}


} //end namespace trackutils
	
} //end namespace vcg

#endif //TRACKUTILS_H
