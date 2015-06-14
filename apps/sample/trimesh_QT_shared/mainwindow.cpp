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

****************************************************************************/


#include "mainwindow.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QDebug>

MainWindow::MainWindow (QWidget * parent)
	:QMainWindow (parent),mi(1000000000),mesh()
{
	ui.setupUi (this);
	QLayout* tmp = ui.glFrame->layout();

    //parent is set to NULL in order to avoid QT bug on MAC (business as usual...).
    //The QGLWidget are destroyed by hand in the MainWindow destructor...
    shared = new SharedDataOpenGLContext(mesh,mi,NULL);
	shared->setHidden(true);
	shared->myInitGL();
	connect (ui.drawModeComboBox, SIGNAL (currentIndexChanged(int)),shared, SLOT (passInfoToOpenGL(int)));

	for(int ii = 0;ii < 2;++ii)
	{
        glar[ii] = new GLArea(mesh,shared->feeder,NULL,shared);
		connect (shared,SIGNAL(dataReadyToBeRead(MyDrawMode)),glar[ii], SLOT (setupEnvironment(MyDrawMode)));
		tmp->addWidget(glar[ii]);
	}

	

	connect (ui.loadMeshPushButton, SIGNAL (clicked()),this, SLOT (chooseMesh()));
	connect (ui.loadTetrahedronPushButton, SIGNAL (clicked()),this, SLOT (loadTetrahedron()));
	connect (ui.loadDodecahedronPushButton, SIGNAL (clicked()),this, SLOT (loadDodecahedron()));
	//from toolFrame to glArea through mainwindow

}

// mesh chooser file dialog
void MainWindow::chooseMesh()
{
	mesh.Clear();
	QString fileName = QFileDialog::getOpenFileName(this,
		tr("Open Mesh"), QDir::currentPath(),
		tr("Poly Model (*.ply)"));
	int err=vcg::tri::io::ImporterPLY<CMesh>::Open(mesh,(fileName.toStdString()).c_str());
	if(err!=0)
	{
		const char* errmsg=vcg::tri::io::ImporterPLY<CMesh>::ErrorMsg(err);
		QMessageBox::warning(this,tr("Error Loading Mesh"),QString(errmsg));
	}
	initMesh(fileName);
}

void MainWindow::loadTetrahedron()
{
	mesh.Clear();
	vcg::tri::Tetrahedron(mesh);
	initMesh(tr("Tethraedron [builtin]"));
}

void MainWindow::loadDodecahedron()
{
	mesh.Clear();
	vcg::tri::Dodecahedron(mesh);
	initMesh(tr("Dodecahedron [builtin]"));
}

void MainWindow::initMesh(QString message)
{
	ui.statusbar->showMessage(message);
	// update bounding box
	vcg::tri::UpdateBounding<CMesh>::Box(mesh);
	// update Normals
	vcg::tri::UpdateNormal<CMesh>::PerVertexNormalizedPerFaceNormalized(mesh);
	shared->feeder.update(vcg::GLMeshAttributesFeeder<CMesh>::ATT_ALL);
	shared->passInfoToOpenGL(ui.drawModeComboBox->currentIndex());
}

MainWindow::~MainWindow()
{
	for(int ii = 0;ii < 2;++ii)
		delete glar[ii];
	delete shared;
}
