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

#include <vcg/space/point2.h>
#include <vcg/space/point3.h>
#include <vcg/space/plane3.h>
#include <vcg/space/intersection3.h>
#include "trackmode.h"
#include <vcg/math/similarity.h>

using namespace vcg;




/* 
Le varie Apply prendono una coppia di punti in screen space e restituiscono una trasformazione
*/

Plane3f SphereMode::GetViewPlane()
{
  Point3f mo=tb->ModelOrigin();
  Point3f vp; vp.Import(tb->camera.ViewPoint());
	Plane3f pl;  // plane perpedicular to view direction and passing through manip center
	pl.Set(vp-tb->center, (vp-tb->center)*tb->center);
  return pl;
}



/* dato un punto in coordinate di schermo e.g. in pixel stile opengl 
   restituisce un punto in coordinate di mondo sulla superficie della trackball 
   La superficie della trackball e' data da una sfera + una porzione di iperboloide di rotazione
   assumiamo la sfera di raggio unitario e centrata sull'origine e di guardare lungo la y negativa.
                                       X   0   sqrt(1/2)  1  
   eq sfera:              y=sqrt(1-x*x);   1   sqrt(1/2)  0   
   eq iperboloide :       y=1/2x;         inf  sqrt(1/2)  1/2

   */

Point3f SphereMode::Hit(const Point3f &p)
{
  printf("Hit in screen space at %5.3f %5.3f %5.3f\n",p[0],p[1],p[2]); 

  Plane3f vp=GetViewPlane(); // plane perpedicular to view direction and passing through manip center
  Line3fN ln= tb->camera.ViewLineFromWindow(Point3f(p[0],p[1],0));
  //Point3f P0,P1;
  Point3f PonVP;
  bool res=Intersection<float>(vp,ln,PonVP);
  const float Thr=tb->radius/math::Sqrt(2.0f);

  Point3f HitPoint;
  float dd=Distance(tb->center,PonVP);
  if(dd<Thr)
  { // First case: We hit the sphere so, we must set the z accordingly 
    float hh=math::Sqrt(tb->radius*tb->radius - dd*dd);
    HitPoint=PonVP+ln.Direction()*hh;
    printf("Hit the sphere point on plane is %5.3f %5.3f %5.3f\n",PonVP[0],PonVP[1],PonVP[2]);
    printf(" Distance from center is %5.3f \n",dd);
    printf(" Heigth for view plane should be %5.3f \n",hh );
  }
  else
  { // Second Case we hit the hyperboloid
    float hh=tb->radius/(2.0f*(dd/tb->radius));
    HitPoint=PonVP+ln.Direction()*hh;

    printf("Hit the hiperboloid at %5.3f %5.3f %5.3f\n",PonVP[0],PonVP[1],PonVP[2]); 
    printf(" Distance from center is %5.3f \n",dd);
    printf(" Heigth for view plane should be %5.3f \n",hh );
    
  }
 return HitPoint;
}

/* 
  Nella trackball classica si considera 

*/
Similarityf PlaneMode::ComputeFromWindow(const Point3f &oldP, const Point3f &newP) 
{ 
  return Similarityf().SetIdentity();
}
/*
Restituisce la trasformazione originata dal drag in window coord da oldp a newp.
*/
Similarityf SphereMode::ComputeFromWindow(const Point3f &oldP, const Point3f &newP)
{
 Point3f hitOld=Hit(oldP);
 Point3f hitNew=Hit(newP);
 // Now compute the rotation defined on the sphere...

 Point3f norm;

 return Similarityf().SetIdentity();

	//Point3f axis = Point3f(0, 0, 1)^result; /* Axis of rotation */
	//axis.Normalize();
	//Point3f d = result - Point3f(0, 0, 1);
	//float t = d.Norm() / 2.0f;
	//if(t > thr)
	//	t += (t - thr) * 0.7f;
	//if (t > 1.0f) t = 1.0f;
	//if (t < -1.0f) t = -1.0f;
 // float phi = 2 * math::Asin(t);
	////return Similarityf().SetRotate(phi * 180/(float)M_PI, axis);
 // return Similarityf().SetRotate(phi, axis);
}
