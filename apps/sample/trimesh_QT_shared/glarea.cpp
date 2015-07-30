/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2007                                                \/)\/    *
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
Revision 1.1  2007/10/18 08:52:06  benedetti
Initial release.


****************************************************************************/


#include "glarea.h"
#include <QKeyEvent>
#include <QKeyEvent>
#include <QWheelEvent>
#include <wrap/qt/trackball.h>
#include <cassert>

#ifdef Q_OS_MAC
#define glGenVertexArrays glGenVertexArraysAPPLE
#define glBindVertexArray glBindVertexArrayAPPLE
#define glDeleteVertexArrays glDeleteVertexArraysAPPLE
#endif

GLArea::GLArea (CMeshO& m, MLThreadSafeGLMeshAttributesFeeder& feed,QWidget* parent,QGLWidget* sharedcont)
    :QGLWidget (parent,sharedcont),mesh(m),feeder(feed),rq(),drawmode(MDM_SMOOTH)
{
}

GLArea::~GLArea()
{
}

void GLArea::initializeGL()
{
	makeCurrent();
	glewExperimental=GL_TRUE;
	glewInit();
	glClearColor(0, 0, 0, 0);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_NORMALIZE);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);
}

void GLArea::resizeGL (int w, int h)
{
	makeCurrent();
	glViewport (0, 0, (GLsizei) w, (GLsizei) h);
    //initializeGL();
}

void GLArea::paintGL ()
{
	makeCurrent();
	glPushAttrib(GL_ALL_ATTRIB_BITS);
    //GLenum err = glGetError();
    //assert(err == GL_NO_ERROR);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluPerspective(25, GLArea::width()/(float)GLArea::height(), 0.1, 100);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	gluLookAt(0,0,5,   0,0,0,   0,1,0);
	track.center=vcg::Point3f(0, 0, 0);
	track.radius= 1;
	track.GetView();
	track.Apply();
	glPushMatrix();
	float d=2.0f/mesh.bbox.Diag();
	vcg::glScale(d);
	glTranslate(-mesh.bbox.Center());

	if (mesh.VN() > 0)
	{
		switch(drawmode)
		{
		case MDM_SMOOTH:
		case MDM_FLAT:
			feeder.drawTriangles(rq);
		case MDM_WIRE:
			feeder.drawWire(rq);
			break;
		case MDM_POINTS:
			feeder.drawPoints(rq);
			break;
		case MDM_FLATWIRE:
			feeder.drawFlatWire(rq);
			break;
		default:
			break;
		}
	}
	glPopMatrix();
	
	track.DrawPostApply();
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glPopAttrib();
	//doneCurrent();
	GLenum err = glGetError();
	assert(err == GL_NO_ERROR);
}

void GLArea::keyReleaseEvent (QKeyEvent * e)
{
	e->ignore ();
	if (e->key () == Qt::Key_Control)
		track.ButtonUp (QT2VCG (Qt::NoButton, Qt::ControlModifier));
	if (e->key () == Qt::Key_Shift)
		track.ButtonUp (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
	if (e->key () == Qt::Key_Alt)
		track.ButtonUp (QT2VCG (Qt::NoButton, Qt::AltModifier));
	updateGL ();
}

void GLArea::keyPressEvent (QKeyEvent * e)
{
	e->ignore ();
	if (e->key () == Qt::Key_Control)
		track.ButtonDown (QT2VCG (Qt::NoButton, Qt::ControlModifier));
	if (e->key () == Qt::Key_Shift)
		track.ButtonDown (QT2VCG (Qt::NoButton, Qt::ShiftModifier));
	if (e->key () == Qt::Key_Alt)
		track.ButtonDown (QT2VCG (Qt::NoButton, Qt::AltModifier));
	updateGL ();
}

void GLArea::mousePressEvent (QMouseEvent * e)
{
	e->accept ();
	setFocus ();
	track.MouseDown (QT2VCG_X(this,e), QT2VCG_Y(this,e), QT2VCG (e->button (), e->modifiers ()));
	updateGL ();
}

void GLArea::mouseMoveEvent (QMouseEvent * e)
{
	if (e->buttons ()) {
		track.MouseMove (QT2VCG_X(this,e), QT2VCG_Y(this,e));
		updateGL ();
	}
}

void GLArea::mouseReleaseEvent (QMouseEvent * e)
{
	track.MouseUp (QT2VCG_X(this,e), QT2VCG_Y(this,e), QT2VCG (e->button (), e->modifiers ()));
	updateGL ();
}

void GLArea::wheelEvent (QWheelEvent * e)
{
	const int WHEEL_STEP = 120;
	track.MouseWheel (e->delta () / float (WHEEL_STEP), QTWheel2VCG (e->modifiers ()));
	updateGL ();
}

void GLArea::updateRequested(MyDrawMode md,vcg::GLFeederInfo::ReqAtts& reqatts)
{
	drawmode = md;
	rq = reqatts;
	updateGL();
}

void GLArea::resetTrackBall()
{
    makeCurrent();
	track.Reset();
	updateGL();
    doneCurrent();
}



SharedDataOpenGLContext::SharedDataOpenGLContext( CMeshO& mesh,MLThreadSafeMemoryInfo& mi,QWidget* parent /*= 0*/ )
	:QGLWidget(parent),feeder(mesh,mi,100000)
{
}

SharedDataOpenGLContext::~SharedDataOpenGLContext()
{
	deAllocateBO();
}

void SharedDataOpenGLContext::myInitGL()
{
    makeCurrent();
    glewInit();
	doneCurrent();
}

void SharedDataOpenGLContext::passInfoToOpenGL(int mode)
{
	makeCurrent();
	MyDrawMode drawmode = static_cast<MyDrawMode>(mode);
    //_tsbm.setUpBuffers();
	vcg::GLFeederInfo::ReqAtts req;
	bool allocated = false;
	switch(drawmode)
	{
	case MDM_SMOOTH:
	case MDM_WIRE:
		req[vcg::GLFeederInfo::ATT_VERTPOSITION] = true;
		req[vcg::GLFeederInfo::ATT_VERTNORMAL] = true;
		req[vcg::GLFeederInfo::ATT_VERTINDEX] = true;
		req.primitiveModality() = vcg::GLFeederInfo::PR_TRIANGLES;
		break;
	case MDM_POINTS:
		req[vcg::GLFeederInfo::ATT_VERTPOSITION] = true;
		req[vcg::GLFeederInfo::ATT_VERTNORMAL] = true;
		req.primitiveModality() = vcg::GLFeederInfo::PR_POINTS;
		break;
	case MDM_FLAT:
	case MDM_FLATWIRE:
		req[vcg::GLFeederInfo::ATT_VERTPOSITION] = true;
		req[vcg::GLFeederInfo::ATT_FACENORMAL] = true;
		req.primitiveModality() = vcg::GLFeederInfo::PR_TRIANGLES;
		break;
	default:
		break;
	}
	vcg::GLFeederInfo::ReqAtts rq = feeder.setupRequestedAttributes(req,allocated);
	doneCurrent();
	emit dataReadyToBeRead(drawmode,rq);

}

void SharedDataOpenGLContext::deAllocateBO()
{
	makeCurrent();
	feeder.deAllocateBO();
	doneCurrent();
}
