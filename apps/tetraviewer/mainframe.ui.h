/****************************************************************************
** ui.h extension file, included from the uic-generated form implementation.
**
** If you want to add, delete, or rename functions or slots, use
** Qt Designer to update this file, preserving your code.
**
** You should not define a constructor or destructor in this file.
** Instead, write your code in functions called init() and destroy().
** These will automatically be called by the form's constructor and
** destructor.
*****************************************************************************/


#include <qfiledialog.h> 
#include <qdir.h> 
#include <qmessagebox.h> 

extern void openTetraMesh(const char *);

void MainFrame::fileNew()
{

}


void MainFrame::fileOpen()
{
	 QString filename = QFileDialog::getOpenFileName(
              "",
              "Tetrahedral Meshes File (*.ts *.ply)",
              this,
              "open file dialog"
              "Choose a TS Tetrahedral mesh file" );
	if (filename!=NULL)
	{
		const char *path=filename.ascii();
		openTetraMesh(path);
	}
}


void MainFrame::fileSave()
{

}


void MainFrame::fileSaveAs()
{

}



void MainFrame::fileExit()
{
	
}


void MainFrame::helpIndex()
{

}


void MainFrame::helpContents()
{

}


void MainFrame::helpAbout()
{

}


void MainFrame::setWire()
{

}
