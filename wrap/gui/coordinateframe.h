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
#ifndef COORDINATEFRAME_H
#define COORDINATEFRAME_H


#include <vcg/math/similarity.h>
#include <vcg/space/color4.h>
#include <wrap/gui/trackball.h>

#include <QGLWidget>

namespace vcg {

class CoordinateFrame
{
public:
  // functions:  
  CoordinateFrame(float);
  virtual ~CoordinateFrame() {}
  virtual void Render(QGLWidget*);
  // data
  Color4b basecolor;
  Color4b xcolor;
  Color4b ycolor;
  Color4b zcolor;
  float size;
  float linewidth;
  QFont font;
  bool drawaxis;
  bool drawlabels;
  bool drawvalues;
  
protected:
  // functions:  
  void drawTickedLine(const Point3d &, const Point3d &, float, float,float);
  float calcSlope(const Point3d &, const Point3d &, float, int , double *, double *, GLint *);
  float niceRound(float);
};

class MovableCoordinateFrame: public CoordinateFrame
{

public:
  MovableCoordinateFrame(float);
  virtual ~MovableCoordinateFrame(){}
  virtual void Render(QGLWidget*);
  virtual void Reset(bool ,bool);
  virtual void SetPosition(const Point3f);
  virtual void SetRotation(const Quaternionf);
  virtual Point3f GetPosition();
  virtual Quaternionf GetRotation();
  virtual void GetTransform(Matrix44f &);
  virtual void Flip(const Point3f);
  virtual void AlignWith(const Point3f, const Point3f);

protected:
  // data:
  Point3f position;
  Quaternionf rotation;
  
  // functions: 
  virtual void Move(const Similarityf);
  void RotateToAlign(const Point3f, const Point3f);

};

class ActiveCoordinateFrame: public MovableCoordinateFrame
{
public:
  ActiveCoordinateFrame(float);
  virtual ~ActiveCoordinateFrame();
  virtual void Render(QGLWidget*);
  virtual void Reset(bool, bool);
  virtual void SetPosition(const Point3f);
  virtual void SetRotation(const Quaternionf);
  virtual void AlignWith(const Point3f, const Point3f);
  void MouseDown(QPoint,int, int, int);
  void MouseMove(QPoint,int, int); 
  void MouseUp(int, int, int); 
  void ButtonUp(int);
  void ButtonDown(int);
  void SetSnap(float);

  Trackball *manipulator;
  bool drawmoves;
  bool drawrotations;
protected:
  // data:
  const int movx,movy,movz,rotx,roty,rotz;
  Point3f x_axis,y_axis,z_axis;
  float rot_snap_rad,mov_snap;
  QPoint cursor;
  
  // functions: 
  virtual void Move(const Similarityf);
  void Update();

};

}//namespace
#endif /*COORDINATEFRAME_H*/
