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

#include <GL/glew.h>
#include <wrap/gl/math.h>
#include <wrap/gl/space.h>
#include <wrap/gl/addons.h>

#include "coordinateframe.h"

using namespace vcg;



CoordinateFrame::CoordinateFrame(float s)
://QObject(),
basecolor(Color4b::White),xcolor(Color4b::Red)
,ycolor(Color4b::Green),zcolor(Color4b::Blue),size(s),linewidth(2.0)
,font(),drawaxis(true),drawlabels(true),drawvalues(false)
{
  font.setFamily("Helvetica");
}

void CoordinateFrame::Render(QGLWidget* glw)
{
  assert( glw!= NULL);
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_2D);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_POINT_SMOOTH);
  glLineWidth(linewidth);
  glPointSize(linewidth*1.5);

  Point3d o(0,0,0);
  Point3d a(size,0,0);
  Point3d b(0,size,0);
  Point3d c(0,0,size);
  // Get gl state values
  double mm[16],mp[16];
  GLint vp[4];
  glGetDoublev(GL_MODELVIEW_MATRIX,mm);
  glGetDoublev(GL_PROJECTION_MATRIX,mp);
  glGetIntegerv(GL_VIEWPORT,vp);
  float slope_a=calcSlope(-a,a,2*size,10,mm,mp,vp);
  float slope_b=calcSlope(-b,b,2*size,10,mm,mp,vp);
  float slope_c=calcSlope(-c,c,2*size,10,mm,mp,vp);
  float scalefactor = size*0.02f;
  if(drawaxis){
    glBegin(GL_LINES);
      glColor(xcolor);
      glVertex(-a); glVertex(a);
      glColor(ycolor);
      glVertex(-b); glVertex(b);
      glColor(zcolor);
      glVertex(-c); glVertex(c);
    glEnd();
    glColor(basecolor);
    // assi positivi
    drawTickedLine(o,a,size,slope_a,linewidth);  // Draws x axis
    drawTickedLine(o,b,size,slope_b,linewidth);  // Draws y axis
    drawTickedLine(o,c,size,slope_c,linewidth);  // Draws z axis
    //assi negativi
    drawTickedLine(o,-a,size,slope_a,linewidth);  // Draws x axis
    drawTickedLine(o,-b,size,slope_b,linewidth);  // Draws y axis
    drawTickedLine(o,-c,size,slope_c,linewidth);  // Draws z axis
    glPushMatrix();
      glTranslate(a);
      glScalef(scalefactor,scalefactor,scalefactor);
      Add_Ons::Cone(10,linewidth*1.5,linewidth*0.5,true);
    glPopMatrix();
    glPushMatrix();
      glTranslate(b);
      glRotatef(90,0,0,1);
      glScalef(scalefactor,scalefactor,scalefactor);
      Add_Ons::Cone(10,linewidth*1.5,linewidth*0.5,true);
    glPopMatrix();
    glPushMatrix();
      glTranslate(c);
      glRotatef(-90,0,1,0);
      glScalef(scalefactor,scalefactor,scalefactor);
      Add_Ons::Cone(10,linewidth*1.5,linewidth*0.5,true);
    glPopMatrix();
  }
  if(drawlabels){
  	font.setBold(true);
    font.setPixelSize(10);
    glColor(xcolor);
    glw->renderText(size+(scalefactor*3),0,0,QString("X"),font);
    glColor(ycolor);
    glw->renderText(0,size+(scalefactor*3),0,QString("Y"),font);
    glColor(zcolor);
    glw->renderText(0,0,size+(scalefactor*3),QString("Z"),font);
  }  
  if(drawvalues){
  	font.setBold(false);  	
    font.setPixelSize(8);
    float i;
    glColor(Color4b::LightGray);
    for(i=slope_a;i<size;i+=slope_a){
      glw->renderText( i,0,0,QString(" %1").arg(i,3,'f',1),font);
      glw->renderText(-i,0,0,QString("-%1").arg(i,3,'f',1),font);
    }
    for(i=slope_b;i<size;i+=slope_b){
      glw->renderText(0, i,0,QString(" %1").arg(i,3,'f',1),font);
      glw->renderText(0,-i,0,QString("-%1").arg(i,3,'f',1),font);
    }
    for(i=slope_c;i<size;i+=slope_c){
      glw->renderText(0,0, i,QString(" %1").arg(i,3,'f',1),font);
      glw->renderText(0,0,-i,QString("-%1").arg(i,3,'f',1),font);
    }
  }
  
  glPopAttrib();
  assert(!glGetError());  
}

void CoordinateFrame::drawTickedLine(const Point3d &a,const Point3d &b, float dim,float tickDist,float linewidth)
{
  Point3d v(b-a);
  v = v /dim; // normalize without computing square roots and powers

  glBegin(GL_POINTS);
  float i;
  for(i=tickDist;i<dim;i+=tickDist)
    glVertex3f(a[0] + i*v[0],a[1] + i*v[1],a[2] + i*v[2]);
  glEnd();

  glPushAttrib(GL_POINT_BIT);
  glPointSize(linewidth*3);  
  glBegin(GL_POINTS);
       glVertex3f(a[0] + dim*v[0],a[1] + dim*v[1],a[2] + dim*v[2]);
  glEnd();

  glPopAttrib();
}

float CoordinateFrame::calcSlope(const Point3d &a,const Point3d &b,float dim,int spacing,double *mm,double *mp,GLint *vp)
{
   Point3d p1,p2;

  gluProject(a[0],a[1],a[2],mm,mp,vp,&p1[0],&p1[1],&p1[2]);
  gluProject(b[0],b[1],b[2],mm,mp,vp,&p2[0],&p2[1],&p2[2]);
  p1[2]=p2[2]=0;

  float tickNum = spacing/Distance(p2,p1);// pxl spacing
  float slope = dim*tickNum;
  float nslope = math::Min(
          math::Min(niceRound(slope), 0.5f*niceRound(2.0f*slope)), 
                                      0.2f*niceRound(5.0f*slope));
  nslope = math::Max<float>(niceRound(dim*.001f),nslope); // prevent too small slope
  return nslope;
}

float CoordinateFrame::niceRound(float val)
{
  return powf(10.f,ceil(log10(val)));
}

MovableCoordinateFrame::MovableCoordinateFrame(float size)
:CoordinateFrame(size),position(0,0,0),rotation(0,Point3f(1,0,0))
{
  // nothing here
}

void MovableCoordinateFrame::Render(QGLWidget* gla)
{
  glPushMatrix();
  
  glTranslate(position);  
  Matrix44f mrot; 
  rotation.ToMatrix(mrot);
  glMultMatrix(mrot);
  
  CoordinateFrame::Render(gla);
  
  glPopMatrix();
}

void MovableCoordinateFrame::GetTransform(Matrix44f & transform)
{
  // costruisco la matrice che porta le coordinate in spazio di mondo

  // resetto la trasf
  transform.SetIdentity();

  // ruoto
  Matrix44f rot;
  rotation.ToMatrix(rot);
  
  transform = rot * transform ;
  
  // sposto in posizione
  Matrix44f pos;
  pos.SetTranslate(position);
  
  transform = pos * transform;
  
}

void MovableCoordinateFrame::Reset(bool reset_position,bool reset_alignment)
{
  if(reset_position)
    position = Point3f(0,0,0);
  if(reset_alignment)
    rotation = Quaternionf(0,Point3f(1,0,0));
}

void MovableCoordinateFrame::SetPosition(const Point3f newpos)
{
  position = newpos;
}

void MovableCoordinateFrame::SetRotation(const Quaternionf newrot)
{
  rotation = newrot;
}

Point3f MovableCoordinateFrame::GetPosition()
{
  return position;
}

Quaternionf MovableCoordinateFrame::GetRotation()
{
  return rotation;
}

void MovableCoordinateFrame::Flip(const Point3f axis)
{
  Similarityf s;
  s.SetRotate(M_PI,Inverse(rotation).Rotate(axis));
  Move(s);
}

void MovableCoordinateFrame::AlignWith(const Point3f pri,const Point3f secondary)
{
  const float EPSILON=1e-6;
  Point3f primary=pri;

  if( primary.Norm() < EPSILON*size )
    return;

  primary.Normalize(); // ho l'asse primario, lo normalizzo
  Plane3f plane(0,primary); // piano di proiezione per la seconda rotazione 
  Point3f old_z = Inverse(rotation).Rotate(Point3f(0,0,1)); // l'asse z
  Point3f old_y_pro = plane.Projection(Inverse(rotation).Rotate(Point3f(0,1,0))); //la proiezione dell'asse y
  Point3f old_x_pro = plane.Projection(Inverse(rotation).Rotate(Point3f(1,0,0))); //la proiezione dell'asse x

  // allinea l'asse z corrente all'asse primary
  RotateToAlign(old_z,primary); // prima rotazione

  Point3f secondary_pro = plane.Projection(secondary); // la proiezione di secondary
  Point3f new_y_pro = plane.Projection(Inverse(rotation).Rotate(Point3f(0,1,0))); // la proiezione dell'asse y dopo la prima rotazione

  // se c'e` un asse secondary e la sua proiezione non e` 0
  if( secondary.Norm() > EPSILON*size && secondary_pro.Norm() > EPSILON ){
    // allinea la proiezione dell'asse y dopo la prima rotazione alla proiezione dell'asse secondary
    secondary_pro.Normalize();
    RotateToAlign(new_y_pro,secondary_pro);
    return;
  }
  // creco di riallineare la y
  if ( old_y_pro.Norm() > EPSILON ) {
    // allinea la proiezione dell'asse y dopo la prima rotazione alla proiezione dell'asse y 
    old_y_pro.Normalize();
    RotateToAlign(new_y_pro,old_y_pro);
    return;
  }
  // cerco di riallineare la x
  Point3f new_x_pro = plane.Projection(Inverse(rotation).Rotate(Point3f(1,0,0))); //la proiezione dell'asse x dopo la prima rotazione
  assert(old_x_pro.Norm() > EPSILON ); // la proiezione dell'asse x non dovrebbe essere 0
  // allinea la proiezione dell'asse x dopo la prima rotazione alla proiezione dell'asse x 
  old_x_pro.Normalize();
  RotateToAlign(new_x_pro,old_x_pro);
}

void MovableCoordinateFrame::Move(const Similarityf track)
{
  position = position + track.tra;
  rotation = rotation * track.rot;
}

void MovableCoordinateFrame::RotateToAlign(const Point3f source, const Point3f dest)
{
  const float EPSILON=1e-6;
  // source e dest devono essere versori
  assert( math::Abs(source.Norm() - 1) < EPSILON);
  assert( math::Abs(dest.Norm() - 1) < EPSILON);
  
  Point3f axis = dest ^ source;
  float sinangle = axis.Norm();
  float cosangle = dest * source;
  float angle = math::Atan2(sinangle,cosangle);  
  
  if( math::Abs(angle) < EPSILON )    
    return; // angolo ~ 0, annullo
  
  if( math::Abs(math::Abs(angle)-M_PI) < EPSILON){
    // devo trovare un asse su cui flippare
    Plane3f plane(0,source);
    axis=plane.Projection(Point3f(1,0,0)); // proietto un punto a caso sul piano normale a source
  	if(axis.Norm() < EPSILON){ // source era ~ [1,0,0]...
  	  axis=plane.Projection(Point3f(0,1,0)); 
      assert(axis.Norm() > EPSILON); // quest'altro punto deve andare bene
  	}
  }  
  rotation = rotation * Quaternionf(angle,axis);   
}

ActiveCoordinateFrame::ActiveCoordinateFrame(float size)
:MovableCoordinateFrame(size),manipulator(NULL),drawmoves(true),
drawrotations(true),
movx(Trackball::BUTTON_RIGHT ),
movy(Trackball::BUTTON_RIGHT | Trackball::KEY_CTRL),
movz(Trackball::BUTTON_RIGHT | Trackball::KEY_SHIFT),
rotx(Trackball::BUTTON_LEFT),
roty(Trackball::BUTTON_LEFT | Trackball::KEY_CTRL),
rotz(Trackball::BUTTON_LEFT | Trackball::KEY_SHIFT),
x_axis(1,0,0),y_axis(0,1,0),z_axis(0,0,1),rot_snap_rad(0.0f),mov_snap(0.0f)
{
  manipulator=new Trackball();
  std::map<int, TrackMode *>::iterator it;
  for(it = manipulator->modes.begin(); it != manipulator->modes.end(); it++)
  {
    if ((*it).second)
      delete (*it).second;
  }
  manipulator->modes.clear();
  manipulator->modes[0] = NULL;    
  Update();
 }

ActiveCoordinateFrame::~ActiveCoordinateFrame()
{
   if(manipulator!=NULL) {
     delete manipulator;
     manipulator=NULL;
  }
}

void ActiveCoordinateFrame::Render(QGLWidget* glw)
{
  glPushMatrix(); //occhio

  manipulator->radius=size;
  manipulator->center=position;
  manipulator->GetView();
  manipulator->Apply(false);  
   
  MovableCoordinateFrame::Render(glw);
  
  // non devo disegnare
  if(!drawmoves && !drawrotations){
    glPopMatrix(); //occhio
    return;  
  }

  int current_mode=manipulator->current_button;  
  bool rotating=(current_mode==rotx)||(current_mode==roty)||(current_mode==rotz);
  bool moving=(current_mode==movx)||(current_mode==movy)||(current_mode==movz);

  // non sto draggando o sto draggando quello che non devo disengnare
  if( (!rotating && !moving)||
      (!((rotating && drawrotations)||(moving && drawmoves)))
    ){
    glPopMatrix(); //occhio
    return;  
  }
  
  // sto draggando e devo disegnare qualcosa
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_2D);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_POINT_SMOOTH);
  
  Color4b color;  
  QString message;
  char axis_name;
  float verse;
  if(rotating && drawrotations){ // devo disegnare una rotazione
    Point3f axis, arc_point;
    float angle;
    manipulator->track.rot.ToAxis(angle,axis);
    if(current_mode==rotx){
      verse=((axis+x_axis).Norm()<1?-1:1);
      color=xcolor; axis_name='x'; arc_point=y_axis*(size*0.8);
    } else if(current_mode==roty) {
      verse=((axis+y_axis).Norm()<1?-1:1);
      color=ycolor; axis_name='y'; arc_point=z_axis*(size*0.8);
    } else if(current_mode==rotz) {
      verse=((axis+z_axis).Norm()<1?-1:1);
      color=zcolor; axis_name='z'; arc_point=x_axis*(size*0.8);
    } else assert(0); // doveva essere una rotazione
    // normalizzo la rotazione a [-180,180]
    float sign = ((angle*verse)<0) ? -1 : 1;
    float abs_angle = (angle<0) ? -angle : angle;
    angle = sign * ( (abs_angle>M_PI) ? 2*M_PI-abs_angle : abs_angle );
    axis = axis * verse;    
    message = QString("rotated %1 deg around %2")
                      .arg(((angle*180.0)/M_PI),5,'f',3)
                      .arg(axis_name);
    Quaternionf arc_rot;
    arc_rot.FromAxis(angle/18.0,axis);
    glColor(color);
    glBegin(GL_POLYGON);   
      glVertex(position);
      glVertex(position+arc_point);
      for(int i=0;i<18;i++){
      	 arc_point = arc_rot.Rotate(arc_point);
         glVertex(position+arc_point);
      }
    glEnd(); 
  } else if(moving && drawmoves){// devo disegnare una traslazione
    Point3f ntra=manipulator->track.tra;
    ntra.Normalize();
    if(current_mode==movx){
      verse=((ntra+x_axis).Norm()<1?-1:1);      
      color=xcolor; axis_name='x';
    }else if(current_mode==movy){
      verse=((ntra+y_axis).Norm()<1?-1:1);      
      color=ycolor; axis_name='y';
    }else if(current_mode==movz){
      verse=((ntra+z_axis).Norm()<1?-1:1);      
      color=zcolor; axis_name='z';
    }else assert(0); // doveva essere una traslazione
    message = QString("moved %1 units along %2")
                      .arg(verse*manipulator->track.tra.Norm(),5,'f',3)
                      .arg(axis_name);
    Point3f old_pos = position-manipulator->track.tra;
    glLineWidth(2*linewidth);
    glPointSize(4*linewidth);
    glColor(color);
    glBegin(GL_LINES);
      glVertex(position);
      glVertex(old_pos);
    glEnd();
    glBegin(GL_POINTS);
      glVertex(old_pos);
    glEnd();    
  } else assert(0); // qualcosa lo dovevo disegnare
  // disegno la stringa  
  font.setBold(true);
  font.setPixelSize(12);
  glw->renderText(cursor.x()+16,cursor.y()+16,message,font);

  glPopAttrib();
  glPopMatrix(); //occhio
}


void ActiveCoordinateFrame::Reset(bool reset_position,bool reset_alignment)
{
  MovableCoordinateFrame::Reset(reset_position, reset_alignment);
  Update();
  manipulator->Reset();  
}

void ActiveCoordinateFrame::SetPosition(const Point3f newpos)
{
  MovableCoordinateFrame::SetPosition(newpos);
  Update();
  manipulator->Reset();
}

void ActiveCoordinateFrame::SetRotation(const Quaternionf newrot)
{
  MovableCoordinateFrame::SetRotation(newrot);
  Update();
  manipulator->Reset();
}

void ActiveCoordinateFrame::AlignWith(const Point3f primary,const Point3f secondary=Point3f(0,0,0))
{
  MovableCoordinateFrame::AlignWith(primary,secondary);
  Update();
  manipulator->Reset();
}

void ActiveCoordinateFrame::MouseDown(QPoint c,int x, int y, /*Button*/ int button)
{
  cursor=c;
  Move(manipulator->track);
  manipulator->Reset();
  manipulator->MouseDown(x,y,button);
}

void ActiveCoordinateFrame::MouseMove(QPoint c,int x, int y)
{
  cursor=c;
  manipulator->MouseMove(x,y);
}

void ActiveCoordinateFrame::MouseUp(int x, int y, /*Button */ int button) 
{
  Move(manipulator->track);
  manipulator->Reset();
  manipulator->MouseUp(x, y, button);
}

void ActiveCoordinateFrame::ButtonUp(int button)
{
  Move(manipulator->track);
  manipulator->Reset();
  manipulator->ButtonUp((Trackball::Button) button);
}

void ActiveCoordinateFrame::ButtonDown(int button)
{
  Move(manipulator->track);
  manipulator->Reset();
  manipulator->ButtonDown((Trackball::Button) button);
}

void ActiveCoordinateFrame::SetSnap(float rot_deg)
{
  assert((rot_deg>=0.0)&&(rot_deg<=180));
  rot_snap_rad=rot_deg*M_PI/180.0;
  Update();
}

void ActiveCoordinateFrame::Move(const Similarityf track)
{
  MovableCoordinateFrame::Move(track);
  Update();
}

void ActiveCoordinateFrame::Update()
{
  Point3f p=position;
  Quaternionf r=Inverse(rotation);
  x_axis=r.Rotate(Point3f(1,0,0));
  y_axis=r.Rotate(Point3f(0,1,0));
  z_axis=r.Rotate(Point3f(0,0,1));
  
  if(manipulator->modes[movx]!=NULL) 
    delete manipulator->modes[movx];
  manipulator->modes[movx] = new AxisMode(p,x_axis);
  
  if(manipulator->modes[movy]!=NULL) 
    delete manipulator->modes[movy];
  manipulator->modes[movy] = new AxisMode(p,y_axis);
  
  if(manipulator->modes[movz]!=NULL)
    delete manipulator->modes[movz];
  manipulator->modes[movz] = new AxisMode(p,z_axis);
  
  if(manipulator->modes[rotx]!=NULL)
    delete manipulator->modes[rotx];
  manipulator->modes[rotx] = new CylinderMode(p,x_axis,rot_snap_rad);
  
  if(manipulator->modes[roty]!=NULL)
    delete manipulator->modes[roty];
  manipulator->modes[roty] = new CylinderMode(p,y_axis,rot_snap_rad);
  
  if(manipulator->modes[rotz]!=NULL)
    delete manipulator->modes[rotz];
  manipulator->modes[rotz] = new CylinderMode(p,z_axis,rot_snap_rad);
  
  manipulator->SetCurrentAction();
}
