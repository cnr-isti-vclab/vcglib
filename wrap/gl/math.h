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
Revision 1.5  2004/05/12 20:54:55  ponchio
*** empty log message ***

Revision 1.4  2004/05/12 13:07:47  ponchio
Added #include <glew.h>

Revision 1.3  2004/05/04 23:36:23  cignoni
remove include of gl and added glextgension exploiting,

Revision 1.2  2004/04/07 10:47:03  cignoni
inlined functions for avoid multiple linking errors

Revision 1.1  2004/03/31 15:27:17  ponchio
*** empty log message ***


****************************************************************************/

#ifndef VCG_GL_MATH_H
#define VCG_GL_MATH_H

// Please note that this file assume that you have already included your 
// gl-extension wrapping utility, and that therefore all the extension symbol are already defined.

#include <vcg/math/matrix44.h>
#include <vcg/math/similarity.h>

namespace vcg {

inline void glMultMatrix(const Matrix44f &matrix) {
  //glMultMatrixf((const GLfloat *)(matrix[0]));  
  glMultTransposeMatrixf((const GLfloat *)(matrix[0])); 
}

inline void glMultMatrix(const Matrix44d &matrix) {
//  glMultMatrixd((const GLdouble *)(matrix[0]));
  glMultTransposeMatrixd((const GLdouble *)(matrix[0])); 
}

inline void glMultMatrix(const Similarityf &s) {
  glTranslatef(s.tra[0], s.tra[1], s.tra[2]);
  glScalef(s.sca, s.sca, s.sca);
  float alpha;
  Point3f axis;
  s.rot.ToAxis(alpha, axis);    
  glRotatef(math::ToDeg(alpha), axis[0], axis[1], axis[2]);    
  
}

inline void glMultMatrix(const Similarityd &s) {
  glTranslated(s.tra[0], s.tra[1], s.tra[2]);
  double alpha;
  Point3d axis;
  s.rot.ToAxis(alpha, axis);
  glRotated(math::ToDeg(alpha), axis[0], axis[1], axis[2]);
  glScaled(s.sca, s.sca, s.sca);
}

}//namespace
#endif
