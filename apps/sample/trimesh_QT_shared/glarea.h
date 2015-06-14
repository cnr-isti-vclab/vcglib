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

#ifndef GLAREA_H_
#define GLAREA_H_

/// Opengl related imports
#include <GL/glew.h>
#include <QGLWidget>

/// wrapper imports
#include <wrap/gui/trackball.h>


#include "mesh.h"
#include "ml_scene_renderer.h"
#include "ml_atomic_guard.h"
/// declaring edge and face type


enum MyDrawMode{MDM_SMOOTH=0,MDM_POINTS,MDM_WIRE,MDM_FLATWIRE,MDM_FLAT};

class SharedDataOpenGLContext : public QGLWidget
{
	Q_OBJECT
public:
	SharedDataOpenGLContext(CMesh& mesh,MLThreadSafeMemoryInfo& mi,QWidget* parent = 0);
	~SharedDataOpenGLContext();

    void myInitGL();

	MLThreadSafeGLMeshAttributesFeeder feeder;

public slots:
	/// widget-based user interaction slots
	void passInfoToOpenGL(int mode);

signals:
	void dataReadyToBeRead(MyDrawMode mode);
private:
	GLuint vaohandlespecificicforglcontext;
	MyDrawMode drawmode;
};

class GLArea:public QGLWidget
{
	Q_OBJECT 
public:
    GLArea (CMesh& m,MLThreadSafeGLMeshAttributesFeeder& feed,QWidget* parent = NULL,QGLWidget* sharedcont = NULL);
	~GLArea();
	/// we choosed a subset of the avaible drawing modes
public slots:
	void setupEnvironment(MyDrawMode mode);

signals:
		/// signal for setting the statusbar message
		void setStatusBar(QString message);
protected:
	/// opengl initialization and drawing calls
	void initializeGL ();
	void resizeGL (int w, int h);
	void paintGL ();
	/// keyboard and mouse event callbacks
	void keyReleaseEvent(QKeyEvent * e);
	void keyPressEvent(QKeyEvent * e);
	void mousePressEvent(QMouseEvent*e);
	void mouseMoveEvent(QMouseEvent*e);
	void mouseReleaseEvent(QMouseEvent*e);
	void wheelEvent(QWheelEvent*e); 
private:
	MLAtomicGuard sem;

	GLuint vaohandlespecificicforglcontext;
	/// the active mesh instance
	CMesh& mesh;
	/// the active manipulator
	vcg::Trackball track;
	/// mesh data structure initializer
	void initMesh(QString message);
	//MLThreadSafeMemoryInfo& mi;
	MLThreadSafeGLMeshAttributesFeeder& feeder;

	MyDrawMode drawmode;
};

//class GLAreaEXT:public GLArea
//{
//	Q_OBJECT 
//public:
//	GLAreaEXT (CMesh& m,MLThreadSafeGLMeshAttributesFeeder& feed,QWidget * parent = NULL,QGLWidget* shared = NULL);
//	~GLAreaEXT();
//
//	void paintGL ();
//};

#endif /*GLAREA_H_ */
