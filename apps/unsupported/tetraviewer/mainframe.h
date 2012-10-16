/****************************************************************************
** Form interface generated from reading ui file 'MainFrame.ui'
**
** Created: Mon Oct 4 19:00:57 2004
**      by: The User Interface Compiler ($Id: mainframe.h,v 1.1 2004-10-04 18:01:36 ganovelli Exp $)
**
** WARNING! All changes made in this file will be lost!
****************************************************************************/

#ifndef MAINFRAME_H
#define MAINFRAME_H

#include <qvariant.h>
#include <qpixmap.h>
#include <qmainwindow.h>

class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QSpacerItem;
class QAction;
class QActionGroup;
class QToolBar;
class QPopupMenu;
class MyGLWidget;
class QGroupBox;
class QPushButton;
class QButtonGroup;

class MainFrame : public QMainWindow
{
    Q_OBJECT

public:
    MainFrame( QWidget* parent = 0, const char* name = 0, WFlags fl = WType_TopLevel );
    ~MainFrame();

    QGroupBox* file;
    QPushButton* OpenButton;
    QButtonGroup* buttonGroup1;
    QPushButton* BoxButton;
    QPushButton* WireButton;
    QPushButton* HiddenButton;
    QPushButton* FlatWireButton;
    QPushButton* SmoothButton;
    QPushButton* SmallTetraButton;
    QPushButton* FlatButton;
    MyGLWidget* myGLWidget;
    QButtonGroup* buttonGroup2;
    QPushButton* TrackButton;
    QPushButton* SectionButton;
    QMenuBar *MenuBar;
    QPopupMenu *File;
    QPopupMenu *Help;
    QPopupMenu *Info_2;
    QAction* fileNewAction;
    QAction* fileOpenAction;
    QAction* fileSaveAction;
    QAction* fileSaveAsAction;
    QAction* filePrintAction;
    QAction* fileExitAction;
    QAction* helpContentsAction;
    QAction* helpIndexAction;
    QAction* helpAboutAction;
    QAction* new_menunew_itemAction;
    QAction* infonew_itemAction;
    QAction* infoSimplexAction;
    QAction* infoQualityAction;
    QAction* infoPhysicsAction;

public slots:
    virtual void fileNew();
    virtual void fileOpen();
    virtual void fileSave();
    virtual void fileSaveAs();
    virtual void fileExit();
    virtual void helpIndex();
    virtual void helpContents();
    virtual void helpAbout();
    void setWire();

protected:

protected slots:
    virtual void languageChange();

private:
    QPixmap image0;

};

#endif // MAINFRAME_H
