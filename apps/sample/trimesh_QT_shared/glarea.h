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
/// declaring edge and face type

class GLArea:public QGLWidget
{
	Q_OBJECT 
public:
	GLArea (CMesh& m,MLThreadSafeGLMeshAttributesFeeder& feed,QWidget * parent = 0);
	~GLArea();
	/// we choosed a subset of the avaible drawing modes
	enum DrawMode{SMOOTH=0,POINTS,WIRE,FLATWIRE,FLAT};
	public slots:
		/// widget-based user interaction slots
		void selectDrawMode(int mode);
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
	GLuint vaohandlespecificicforglcontext;
	/// the active mesh instance
	CMesh& mesh;
	/// the active manipulator
	vcg::Trackball track;
	/// the current drawmode
	DrawMode drawmode;
	/// mesh data structure initializer
	void initMesh(QString message);
	//MLThreadSafeMemoryInfo& mi;
	MLThreadSafeGLMeshAttributesFeeder& feeder;
};

#endif /*GLAREA_H_ */
