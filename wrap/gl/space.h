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
Revision 1.2  2004/05/13 23:44:47  ponchio
<GL/GL.h>  -->   <GL/gl.h>

Revision 1.1  2004/04/05 11:56:14  cignoni
First working version!


****************************************************************************/

#ifndef VCG_GL_SPACE_H
#define VCG_GL_SPACE_H

// Please note that this file assume that you have already included your 
// gl-extension wrapping utility, and that therefore all the extension symbol are already defined.

#include <vcg/space/point2.h>
#include <vcg/space/point3.h>
#include <vcg/space/color4.h>

namespace vcg {

	inline void glVertex(Point2<int> const & p)   { glVertex2iv(p.V());}
	inline void glVertex(Point2<short> const & p) { glVertex2sv(p.V());}
	inline void glVertex(Point2<float> const & p) { glVertex2fv(p.V());}
	inline void glVertex(Point2<double> const & p){ glVertex2dv(p.V());}
	inline void glTexCoord(Point2<int> const & p)   { glTexCoord2iv(p.V());}
	inline void glTexCoord(Point2<short> const & p) { glTexCoord2sv(p.V());}
	inline void glTexCoord(Point2<float> const & p) { glTexCoord2fv(p.V());}
	inline void glTexCoord(Point2<double> const & p){ glTexCoord2dv(p.V());}
	inline void glTranslate(Point2<float> const & p) { glTranslatef(p.X(),p.Y(),0);}
	inline void glTranslate(Point2<double> const & p){ glTranslated(p.X(),p.Y(),0);}
	inline void glScale(Point2<float> const & p) { glScalef(p.X(),p.Y(),0);}
	inline void glScale(Point2<double> const & p){ glScaled(p.X(),p.Y(),0);}

  inline void glVertex(Point3<int> const & p)   { glVertex3iv(p.V());}
	inline void glVertex(Point3<short> const & p) { glVertex3sv(p.V());}
	inline void glVertex(Point3<float> const & p) { glVertex3fv(p.V());}
	inline void glVertex(Point3<double> const & p){ glVertex3dv(p.V());}
	inline void glNormal(Point3<int> const & p)   { glNormal3iv(p.V());}
	inline void glNormal(Point3<short> const & p) { glNormal3sv(p.V());}
	inline void glNormal(Point3<float> const & p) { glNormal3fv(p.V());}
	inline void glNormal(Point3<double> const & p){ glNormal3dv(p.V());}
	inline void glTexCoord(Point3<int> const & p)   { glTexCoord3iv(p.V());}
	inline void glTexCoord(Point3<short> const & p) { glTexCoord3sv(p.V());}
	inline void glTexCoord(Point3<float> const & p) { glTexCoord3fv(p.V());}
	inline void glTexCoord(Point3<double> const & p){ glTexCoord3dv(p.V());}
	inline void glTranslate(Point3<float> const & p) { glTranslatef(p.X(),p.Y(),p.Z());}
	inline void glTranslate(Point3<double> const & p){ glTranslated(p.X(),p.Y(),p.Z());}
	inline void glScale(Point3<float> const & p) { glScalef(p.X(),p.Y(),p.Z());}
	inline void glScale(Point3<double> const & p){ glScaled(p.X(),p.Y(),p.Z());}

  inline void glColor(Color4b const & c)   { glColor4ubv(c.V());}
  inline void glClearColor(Color4b const &c) {	::glClearColor(float(c[0])/255.0f,float(c[1])/255.0f,float(c[2])/255.0f,1.0f);}
  inline void glLight(GLenum light, GLenum pname,  Color4b const & c)   { 
    static float cf[4]; 
    cf[0]=float(cf[0]/255.0); cf[1]=float(c[1]/255.0); cf[2]=float(c[2]/255.0); cf[3]=float(c[3]/255.0); 
    glLightfv(light,pname,cf);
  }

}//namespace
#endif
