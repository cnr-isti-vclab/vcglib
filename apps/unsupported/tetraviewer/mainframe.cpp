/****************************************************************************
** Form implementation generated from reading ui file 'MainFrame.ui'
**
** Created: Mon Oct 4 19:00:57 2004
**      by: The User Interface Compiler ($Id: mainframe.cpp,v 1.2 2004-10-04 18:45:48 ganovelli Exp $)
**
** WARNING! All changes made in this file will be lost!
****************************************************************************/

#include "mainframe.h"

#include <qvariant.h>
#include "myglwidget.h"
#include <qgroupbox.h>
#include <qpushbutton.h>
#include <qbuttongroup.h>
#include <qlayout.h>
#include <qtooltip.h>
#include <qwhatsthis.h>
#include <qaction.h>
#include <qmenubar.h>
#include <qpopupmenu.h>
#include <qtoolbar.h>
#include <qimage.h>
#include <qpixmap.h>

#include "MainFrame.ui.h"
/*
 *  Constructs a MainFrame as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 *
 */
MainFrame::MainFrame( QWidget* parent, const char* name, WFlags fl )
    : QMainWindow( parent, name, fl )
{
    (void)statusBar();
    if ( !name )
	setName( "MainFrame" );
    setCentralWidget( new QWidget( this, "qt_central_widget" ) );

    file = new QGroupBox( centralWidget(), "file" );
    file->setGeometry( QRect( 70, 0, 90, 80 ) );

    OpenButton = new QPushButton( file, "OpenButton" );
    OpenButton->setGeometry( QRect( 10, 10, 70, 60 ) );
    OpenButton->setPixmap( QPixmap::fromMimeSource( "Open64.png" ) );

    buttonGroup1 = new QButtonGroup( centralWidget(), "buttonGroup1" );
    buttonGroup1->setGeometry( QRect( 160, 0, 470, 80 ) );
    buttonGroup1->setExclusive( TRUE );

    BoxButton = new QPushButton( buttonGroup1, "BoxButton" );
    BoxButton->setGeometry( QRect( 10, 30, 51, 31 ) );
    BoxButton->setToggleButton( TRUE );

    WireButton = new QPushButton( buttonGroup1, "WireButton" );
    WireButton->setGeometry( QRect( 80, 30, 40, 30 ) );
    WireButton->setToggleButton( TRUE );

    HiddenButton = new QPushButton( buttonGroup1, "HiddenButton" );
    HiddenButton->setGeometry( QRect( 130, 30, 60, 31 ) );
    HiddenButton->setToggleButton( TRUE );

    FlatWireButton = new QPushButton( buttonGroup1, "FlatWireButton" );
    FlatWireButton->setGeometry( QRect( 260, 30, 61, 31 ) );
    FlatWireButton->setToggleButton( TRUE );

    SmoothButton = new QPushButton( buttonGroup1, "SmoothButton" );
    SmoothButton->setGeometry( QRect( 330, 30, 50, 30 ) );
    SmoothButton->setToggleButton( TRUE );

    SmallTetraButton = new QPushButton( buttonGroup1, "SmallTetraButton" );
    SmallTetraButton->setGeometry( QRect( 390, 30, 60, 30 ) );
    SmallTetraButton->setToggleButton( TRUE );

    FlatButton = new QPushButton( buttonGroup1, "FlatButton" );
    FlatButton->setGeometry( QRect( 200, 30, 50, 30 ) );
    FlatButton->setToggleButton( TRUE );
    FlatButton->setOn( TRUE );
    FlatButton->setAutoDefault( FALSE );
    FlatButton->setDefault( FALSE );

    myGLWidget = new MyGLWidget( centralWidget(), "myGLWidget" );
    myGLWidget->setGeometry( QRect( 70, 80, 790, 720 ) );

    buttonGroup2 = new QButtonGroup( centralWidget(), "buttonGroup2" );
    buttonGroup2->setGeometry( QRect( 630, 0, 230, 80 ) );
    buttonGroup2->setExclusive( TRUE );

    TrackButton = new QPushButton( buttonGroup2, "TrackButton" );
    TrackButton->setGeometry( QRect( 140, 20, 61, 31 ) );
    TrackButton->setToggleButton( TRUE );
    TrackButton->setOn( TRUE );
    TrackButton->setDefault( FALSE );

    SectionButton = new QPushButton( buttonGroup2, "SectionButton" );
    SectionButton->setGeometry( QRect( 40, 20, 71, 31 ) );
    SectionButton->setToggleButton( TRUE );

    // actions
    fileNewAction = new QAction( this, "fileNewAction" );
    fileNewAction->setIconSet( QIconSet( QPixmap::fromMimeSource( "" ) ) );
    fileOpenAction = new QAction( this, "fileOpenAction" );
    fileOpenAction->setToggleAction( FALSE );
    fileOpenAction->setOn( FALSE );
    fileOpenAction->setIconSet( QIconSet( QPixmap::fromMimeSource( "" ) ) );
    fileSaveAction = new QAction( this, "fileSaveAction" );
    fileSaveAction->setIconSet( QIconSet( QPixmap::fromMimeSource( "" ) ) );
    fileSaveAsAction = new QAction( this, "fileSaveAsAction" );
    filePrintAction = new QAction( this, "filePrintAction" );
    filePrintAction->setIconSet( QIconSet( QPixmap::fromMimeSource( "" ) ) );
    fileExitAction = new QAction( this, "fileExitAction" );
    helpContentsAction = new QAction( this, "helpContentsAction" );
    helpIndexAction = new QAction( this, "helpIndexAction" );
    helpAboutAction = new QAction( this, "helpAboutAction" );
    new_menunew_itemAction = new QAction( this, "new_menunew_itemAction" );
    infonew_itemAction = new QAction( this, "infonew_itemAction" );
    infoSimplexAction = new QAction( this, "infoSimplexAction" );
    infoSimplexAction->setToggleAction( TRUE );
    infoSimplexAction->setOn( TRUE );
    infoQualityAction = new QAction( this, "infoQualityAction" );
    infoQualityAction->setToggleAction( TRUE );
    infoPhysicsAction = new QAction( this, "infoPhysicsAction" );
    infoPhysicsAction->setToggleAction( TRUE );


    // toolbars


    // menubar
    MenuBar = new QMenuBar( this, "MenuBar" );


    File = new QPopupMenu( this );
    fileNewAction->addTo( File );
    fileOpenAction->addTo( File );
    fileSaveAction->addTo( File );
    fileSaveAsAction->addTo( File );
    File->insertSeparator();
    filePrintAction->addTo( File );
    File->insertSeparator();
    fileExitAction->addTo( File );
    MenuBar->insertItem( QString(""), File, 1 );

    Help = new QPopupMenu( this );
    helpContentsAction->addTo( Help );
    helpIndexAction->addTo( Help );
    Help->insertSeparator();
    helpAboutAction->addTo( Help );
    MenuBar->insertItem( QString(""), Help, 2 );

    Info_2 = new QPopupMenu( this );
    new_menunew_itemAction->addTo( Info_2 );
    infoSimplexAction->addTo( Info_2 );
    infoQualityAction->addTo( Info_2 );
    infoPhysicsAction->addTo( Info_2 );
    MenuBar->insertItem( QString(""), Info_2, 3 );

    languageChange();
    resize( QSize(908, 846).expandedTo(minimumSizeHint()) );
    clearWState( WState_Polished );

    // signals and slots connections
    connect( fileNewAction, SIGNAL( activated() ), this, SLOT( fileNew() ) );
    connect( fileOpenAction, SIGNAL( activated() ), this, SLOT( fileOpen() ) );
    connect( fileSaveAction, SIGNAL( activated() ), this, SLOT( fileSave() ) );
    connect( fileSaveAsAction, SIGNAL( activated() ), this, SLOT( fileSaveAs() ) );
    connect( helpAboutAction, SIGNAL( activated() ), this, SLOT( helpAbout() ) );
    connect( helpContentsAction, SIGNAL( activated() ), this, SLOT( helpContents() ) );
    connect( helpIndexAction, SIGNAL( activated() ), this, SLOT( helpIndex() ) );
    connect( BoxButton, SIGNAL( pressed() ), myGLWidget, SLOT( setBox() ) );
    connect( WireButton, SIGNAL( pressed() ), myGLWidget, SLOT( setWire() ) );
    connect( HiddenButton, SIGNAL( pressed() ), myGLWidget, SLOT( setHidden() ) );
    connect( FlatButton, SIGNAL( pressed() ), myGLWidget, SLOT( setFlat() ) );
    connect( FlatWireButton, SIGNAL( pressed() ), myGLWidget, SLOT( setFlatWire() ) );
    connect( SmoothButton, SIGNAL( pressed() ), myGLWidget, SLOT( setSmooth() ) );
    connect( SmallTetraButton, SIGNAL( pressed() ), myGLWidget, SLOT( setSmallTetra() ) );
    connect( OpenButton, SIGNAL( clicked() ), this, SLOT( fileOpen() ) );
    connect( SectionButton, SIGNAL( pressed() ), myGLWidget, SLOT( SectionMouseModality() ) );
    connect( TrackButton, SIGNAL( pressed() ), myGLWidget, SLOT( TrackMouseModality() ) );
    connect( infoPhysicsAction, SIGNAL( activated() ), myGLWidget, SLOT( SwitchTextPhysics() ) );
    connect( infoQualityAction, SIGNAL( activated() ), myGLWidget, SLOT( SwitchTextQuality() ) );
    connect( infoSimplexAction, SIGNAL( activated() ), myGLWidget, SLOT( SwitchTextSimplex() ) );
}

/*
 *  Destroys the object and frees any allocated resources
 */
MainFrame::~MainFrame()
{
    // no need to delete child widgets, Qt does it all for us
}

/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void MainFrame::languageChange()
{
    setCaption( tr( "TetraView" ) );
    file->setTitle( QString::null );
    OpenButton->setText( QString::null );
    buttonGroup1->setTitle( QString::null );
    BoxButton->setText( tr( "box" ) );
    WireButton->setText( tr( "Wire" ) );
    HiddenButton->setText( tr( "Hidden" ) );
    FlatWireButton->setText( tr( "FlatWire" ) );
    SmoothButton->setText( tr( "Smooth" ) );
    SmallTetraButton->setText( tr( "SmallTetra" ) );
    FlatButton->setText( tr( "Flat" ) );
    buttonGroup2->setTitle( QString::null );
    TrackButton->setText( tr( "Trackball" ) );
    SectionButton->setText( tr( "Section" ) );
    fileNewAction->setText( tr( "New" ) );
    fileNewAction->setMenuText( tr( "&New" ) );
    fileNewAction->setAccel( tr( "Ctrl+N" ) );
    fileOpenAction->setText( tr( "Open" ) );
    fileOpenAction->setMenuText( tr( "&Open..." ) );
    fileOpenAction->setAccel( tr( "Ctrl+O" ) );
    fileSaveAction->setText( tr( "Save" ) );
    fileSaveAction->setMenuText( tr( "&Save" ) );
    fileSaveAction->setAccel( tr( "Ctrl+S" ) );
    fileSaveAsAction->setText( tr( "Save As" ) );
    fileSaveAsAction->setMenuText( tr( "Save &As..." ) );
    fileSaveAsAction->setAccel( QString::null );
    filePrintAction->setText( tr( "Print" ) );
    filePrintAction->setMenuText( tr( "&Print..." ) );
    filePrintAction->setAccel( tr( "Ctrl+P" ) );
    fileExitAction->setText( tr( "Exit" ) );
    fileExitAction->setMenuText( tr( "E&xit" ) );
    fileExitAction->setAccel( QString::null );
    helpContentsAction->setText( tr( "Contents" ) );
    helpContentsAction->setMenuText( tr( "&Contents..." ) );
    helpContentsAction->setAccel( QString::null );
    helpIndexAction->setText( tr( "Index" ) );
    helpIndexAction->setMenuText( tr( "&Index..." ) );
    helpIndexAction->setAccel( QString::null );
    helpAboutAction->setText( tr( "About" ) );
    helpAboutAction->setMenuText( tr( "&About" ) );
    helpAboutAction->setAccel( QString::null );
    new_menunew_itemAction->setText( QString::null );
    new_menunew_itemAction->setMenuText( QString::null );
    infonew_itemAction->setText( tr( "new item" ) );
    infonew_itemAction->setMenuText( tr( "new item" ) );
    infoSimplexAction->setText( tr( "Simplex" ) );
    infoSimplexAction->setMenuText( tr( "Simplex" ) );
    infoQualityAction->setText( tr( "Quality" ) );
    infoQualityAction->setMenuText( tr( "Quality" ) );
    infoPhysicsAction->setText( tr( "Physics" ) );
    infoPhysicsAction->setMenuText( tr( "Physics" ) );
    if (MenuBar->findItem(1))
        MenuBar->findItem(1)->setText( tr( "&File" ) );
    if (MenuBar->findItem(2))
        MenuBar->findItem(2)->setText( tr( "&Help" ) );
    if (MenuBar->findItem(3))
        MenuBar->findItem(3)->setText( tr( "Info" ) );
}

