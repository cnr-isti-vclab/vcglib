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
		connect (shared,SIGNAL(dataReadyToBeRead(MyDrawMode,vcg::GLFeederInfo::ReqAtts&)),glar[ii], SLOT (updateRequested(MyDrawMode,vcg::GLFeederInfo::ReqAtts&)));
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
	int err=vcg::tri::io::ImporterPLY<CMeshO>::Open(mesh,(fileName.toStdString()).c_str());
	if(err!=0)
	{
		const char* errmsg=vcg::tri::io::ImporterPLY<CMeshO>::ErrorMsg(err);
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
	if (shared != NULL)
		shared->deAllocateBO();
	// update bounding box
	vcg::tri::UpdateBounding<CMeshO>::Box(mesh);
	// update Normals
	vcg::tri::UpdateNormal<CMeshO>::PerVertexNormalizedPerFaceNormalized(mesh);
	shared->passInfoToOpenGL(ui.drawModeComboBox->currentIndex());
	for(size_t ii = 0;ii < 2;++ii)
		if (glar[ii] != NULL)
			glar[ii]->resetTrackBall();
    ui.statusbar->showMessage(message);
}

MainWindow::~MainWindow()
{
	for(int ii = 0;ii < 2;++ii)
		delete glar[ii];
	delete shared;
}
