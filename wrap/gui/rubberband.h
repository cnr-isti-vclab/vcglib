/****************************************************************************
 * MeshLab                                                           o o     *
 * A versatile mesh processing toolbox                             o     o   *
 *                                                                _   O  _   *
 * Copyright(C) 2008                                                \/)\/    *
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
#ifndef RUBBERBAND_H
#define RUBBERBAND_H

#include <vcg/space/color4.h>
#include <QGLWidget>

namespace vcg {

class Rubberband
{
public:
  //data:
  Color4b color;  
  // functions:  
  Rubberband(Color4b);
  virtual ~Rubberband() {}
  void Render(QGLWidget*);
  void Drag(QPoint);
  void Pin(QPoint);
  void Reset();
  bool IsReady();  
  void GetPoints(Point3f &,Point3f &);
  void RenderLabel(QString text,QGLWidget* gla);
private:
  // types:
  typedef enum { RUBBER_BEGIN = 0,
	             RUBBER_DRAGGING = 1,
	             RUBBER_DRAGGED = 2,
	           } RubberPhase;
  // data:  
  RubberPhase currentphase;
  QPoint qt_cursor;
  Point3f start, end;
  bool have_to_pick;
  QFont font;
  // functions:
  Point3f PixelConvert(const Point3f);
  
};

}//namespace

#endif /*RUBBERBAND_H*/
