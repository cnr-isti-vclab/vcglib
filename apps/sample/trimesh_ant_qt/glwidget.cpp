/****************************************************************************
**
** Copyright (C) 2011 Nokia Corporation and/or its subsidiary(-ies).
** All rights reserved.
** Contact: Nokia Corporation (qt-info@nokia.com)
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** You may use this file under the terms of the BSD license as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of Nokia Corporation and its Subsidiary(-ies) nor
**     the names of its contributors may be used to endorse or promote
**     products derived from this software without specific prior written
**     permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
** $QT_END_LICENSE$
**
****************************************************************************/

#include <QtGui>
#include <QtOpenGL>
#include <math.h>


#include "glwidget.h"
#include <wrap/qt/trackball.h>

enum DrawMode{SMOOTH=0,POINTS,WIRE,FLATWIRE,HIDDEN,FLAT};

TwBar *bar;
char * filename;/// filename of the mesh to load
CMesh mesh;     /// the active mesh instance
vcg::GlTrimesh<CMesh> glWrap;    /// the active mesh opengl wrapper
vcg::Trackball track;     /// the active manipulator
DrawMode drawmode;     /// the current drawmode
void initMesh(QString message);     /// mesh data structure initializer

void TW_CALL CopyCDStringToClient(char **destPtr, const char *src)
{
    size_t srcLen = (src!=NULL) ? strlen(src) : 0;
    size_t destLen = (*destPtr!=NULL) ? strlen(*destPtr) : 0;

    // Alloc or realloc dest memory block if needed
    if( *destPtr==NULL )
        *destPtr = (char *)malloc(srcLen+1);
    else if( srcLen>destLen )
        *destPtr = (char *)realloc(*destPtr, srcLen+1);

    // Copy src
    if( srcLen>0 )
        strncpy(*destPtr, src, srcLen);
    (*destPtr)[srcLen] = '\0'; // null-terminated string
}

void  TW_CALL loadTetrahedron(void *){
	vcg::tri::Tetrahedron(mesh);
	glWrap.m = &mesh;
  	glWrap.Update();

}

void TW_CALL loadMesh(void *)
{
  if(filename==0) return;
  int err=vcg::tri::io::ImporterPLY<CMesh>::Open(mesh,(char*)filename);
  if(err==ply::E_NOERROR)
  {
    vcg::tri::UpdateBounding<CMesh>::Box(mesh);
    vcg::tri::UpdateNormals<CMesh>::PerVertexNormalizedPerFace(mesh);
    vcg::tri::UpdateNormals<CMesh>::PerFaceNormalized(mesh);
    // Initialize the opengl wrapper
    glWrap.m = &mesh;
    glWrap.Update();
  }
}

GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
{
  filename=0;
     setWindowTitle(tr("Hello GL"));
     bar = TwNewBar("TweakBar");
     TwCopyCDStringToClientFunc (CopyCDStringToClient);

     TwAddVarRW(bar,"Input",TW_TYPE_CDSTRING, &filename," label='Filepath' group=SetMesh help=` Name of the file to load` ");
     TwAddButton(bar,"Load from file",loadMesh,0,	" label='Load Mesh' group=SetMesh help=`load the mesh` ");
     TwAddButton(bar,"Use tetrahedron",loadTetrahedron,0,	" label='Make Tetrahedron' group=SetMesh help=`use tetrahedron.` ");

     // ShapeEV associates Shape enum values with labels that will be displayed instead of enum values
     TwEnumVal drawmodes[6] = { {SMOOTH, "Smooth"}, {POINTS, "Per Points"}, {WIRE, "Wire"}, {FLATWIRE, "FlatWire"},{HIDDEN, "Hidden"},{FLAT, "Flat"}};
     // Create a type for the enum shapeEV
     TwType drawMode = TwDefineEnum("DrawMode", drawmodes, 6);
     // add 'g_CurrentShape' to 'bar': this is a variable of type ShapeType. Its key shortcuts are [<] and [>].
     TwAddVarRW(bar, "Draw Mode", drawMode, &drawmode, " keyIncr='<' keyDecr='>' help='Change draw mode.' ");
}
GLWidget::~GLWidget() { }
QSize GLWidget::sizeHint() const
{
    return QSize(800, 600);
}

void GLWidget::initializeGL ()
{
  glClearColor(0, 0, 0, 0);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_NORMALIZE);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
}

void GLWidget::resizeGL (int w, int h)
{
  glViewport (0, 0, (GLsizei) w, (GLsizei) h);
  TwWindowSize(w, h);
  initializeGL();
 }

void GLWidget::paintGL ()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40, GLWidget::width()/(float)GLWidget::height(), 0.1, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,5,   0,0,0,   0,1,0);
    track.center=vcg::Point3f(0, 0, 0);
    track.radius= 1;
    track.GetView();
    track.Apply(false);
    glPushMatrix();
    if(mesh.vert.size()>0)
    {
      float d=1.0f/mesh.bbox.Diag();
      vcg::glScale(d);
      glTranslate(-mesh.bbox.Center());
      // the trimesh drawing calls
      switch(drawmode)
      {
      case SMOOTH:    glWrap.Draw<vcg::GLW::DMSmooth,   vcg::GLW::CMNone,vcg::GLW::TMNone> ();  break;
      case POINTS:    glWrap.Draw<vcg::GLW::DMPoints,   vcg::GLW::CMNone,vcg::GLW::TMNone> (); break;
      case WIRE:      glWrap.Draw<vcg::GLW::DMWire,     vcg::GLW::CMNone,vcg::GLW::TMNone> ();  break;
      case FLATWIRE:  glWrap.Draw<vcg::GLW::DMFlatWire, vcg::GLW::CMNone,vcg::GLW::TMNone> (); break;
      case HIDDEN:    glWrap.Draw<vcg::GLW::DMHidden,   vcg::GLW::CMNone,vcg::GLW::TMNone> ();  break;
      case FLAT:      glWrap.Draw<vcg::GLW::DMFlat,     vcg::GLW::CMNone,vcg::GLW::TMNone> (); break;
      default:        break;
      }
    }
    glPopMatrix();
    track.DrawPostApply();
    TwDraw();
}

void GLWidget::keyReleaseEvent (QKeyEvent * e)
{
  e->ignore ();
  if (e->key () == Qt::Key_Control)  track.ButtonUp (QT2VCG (Qt::NoButton, Qt::ControlModifier));
  if (e->key () == Qt::Key_Shift)  track.ButtonUp (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
  if (e->key () == Qt::Key_Alt) track.ButtonUp (QT2VCG (Qt::NoButton, Qt::AltModifier));
  updateGL ();
}

int TwKeyPressedQt(QKeyEvent *e)
{
  int kmod = 0;
  if(e->modifiers() & Qt::ShiftModifier )  kmod |= TW_KMOD_SHIFT;
  if(e->modifiers() & Qt::ControlModifier )  kmod |= TW_KMOD_CTRL;
  if(e->modifiers() & Qt::AltModifier )  kmod |= TW_KMOD_ALT;
  int key = e->key();
  int k = 0;

  if( key>0 && key<0x7e ) k=key; // plain ascii codes

  if( key>=Qt::Key_F1 && key<=Qt::Key_F12  )
      k = TW_KEY_F1 + (key-Qt::Key_F1 );
  else
  {
      switch( key )
      {
      case Qt::Key_Left:      k = TW_KEY_LEFT;  break;
      case Qt::Key_Up:        k = TW_KEY_UP; break;
      case Qt::Key_Right:     k = TW_KEY_RIGHT;  break;
      case Qt::Key_Down:      k = TW_KEY_DOWN;   break;
      case Qt::Key_PageUp:    k = TW_KEY_PAGE_UP;  break;
      case Qt::Key_PageDown:  k = TW_KEY_PAGE_DOWN; break;
      case Qt::Key_Home:      k = TW_KEY_HOME; break;
      case Qt::Key_End:       k = TW_KEY_END; break;
      case Qt::Key_Insert:    k = TW_KEY_INSERT; break;
      case Qt::Key_Backspace:    k = TW_KEY_BACKSPACE; break;
      case Qt::Key_Delete:    k = TW_KEY_DELETE; break;
      case Qt::Key_Return:    k = TW_KEY_RETURN; break;
      case Qt::Key_Enter:    k = TW_KEY_RETURN; break;
      case Qt::Key_Escape:    k = TW_KEY_ESCAPE; break;
      case Qt::Key_Tab:    k = TW_KEY_TAB; break;
      }
  }

  return TwKeyPressed(k, kmod);
}


void GLWidget::keyPressEvent (QKeyEvent * e)
{
  e->ignore ();
  if (e->key () == Qt::Key_Control) track.ButtonDown (QT2VCG (Qt::NoButton, Qt::ControlModifier));
  if (e->key () == Qt::Key_Shift)  track.ButtonDown (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
  if (e->key () == Qt::Key_Alt)  track.ButtonDown (QT2VCG (Qt::NoButton, Qt::AltModifier));

  TwKeyPressedQt(e);
  updateGL ();
}

TwMouseButtonID Qt2TwMouseButtonId(QMouseEvent *e)
{
  if(e->button() && Qt::LeftButton) return TW_MOUSE_LEFT;
  if(e->button() && Qt::MidButton) return TW_MOUSE_MIDDLE;
  if(e->button() && Qt::RightButton) return TW_MOUSE_RIGHT;
  assert(0);
}

void GLWidget::mousePressEvent (QMouseEvent * e)
{
  e->accept ();
  setFocus ();
  track.MouseDown (e->x (), height () - e->y (), QT2VCG (e->button (), e->modifiers ()));
  TwMouseMotion(e->x (), e->y ());
  TwMouseButton(TW_MOUSE_PRESSED, Qt2TwMouseButtonId(e));
  updateGL ();
}

void GLWidget::mouseMoveEvent (QMouseEvent * e)
{
  if (e->buttons ()) {
    track.MouseMove (e->x (), height () - e->y ());
    updateGL ();
  }
  TwMouseMotion(e->x (), e->y ());
}

void GLWidget::mouseReleaseEvent (QMouseEvent * e)
{
  track.MouseUp (e->x (), height () - e->y (), QT2VCG (e->button (), e->modifiers ()));
  TwMouseButton(TW_MOUSE_RELEASED, Qt2TwMouseButtonId(e));
  updateGL ();
}

void GLWidget::wheelEvent (QWheelEvent * e)
{
  const int WHEEL_STEP = 120;
  track.MouseWheel (e->delta () / float (WHEEL_STEP), QTWheel2VCG (e->modifiers ()));
  updateGL ();
}
