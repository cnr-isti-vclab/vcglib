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
#ifndef ACTIVECOORDINATEFRAME_H
#define ACTIVECOORDINATEFRAME_H

#include "coordinateframe.h"
#include <wrap/gui/trackball.h>
#include <QGLWidget>

namespace vcg {

class ActiveCoordinateFrame: public MovableCoordinateFrame
{
public:
  ActiveCoordinateFrame(float);
  virtual ~ActiveCoordinateFrame();
  virtual void Render(QGLWidget*);
  virtual void Reset(bool, bool);
  virtual void SetPosition(const Point3f);
  virtual void SetRotation(const Quaternionf);
  virtual void AlignWith(const Point3f, const Point3f,const char,const char);
  void MouseDown(int, int, int);
  void MouseMove(int, int); 
  void MouseUp(int, int, int); 
  void ButtonUp(int);
  void ButtonDown(int);
  void SetSnap(float);

  Trackball *manipulator;
  bool drawmoves;
  bool drawrotations;
protected:
  // data:
  const int move_button,rotate_button;
  const int x_modifier,y_modifier,z_modifier;
  Point3f x_axis,y_axis,z_axis;
  float rot_snap_rad,mov_snap;
  // functions: 
  virtual void Move(const Similarityf);
  void Update();
private:
  int movx,movy,movz,rotx,roty,rotz;
};

}//namespace
#endif /*ACTIVECOORDINATEFRAME_H*/
