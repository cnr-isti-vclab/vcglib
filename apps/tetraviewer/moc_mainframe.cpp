/****************************************************************************\
** MainFrame meta object code from reading C++ file 'mainframe.h'
**
** Created: Mon Oct 4 19:00:57 2004
**      by: The Qt MOC ($Id: moc_mainframe.cpp,v 1.1 2004-10-04 18:01:36 ganovelli Exp $)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#undef QT_NO_COMPAT
#include "mainframe.h"
#include <qmetaobject.h>
#include <qapplication.h>

#include <private/qucomextra_p.h>
#if !defined(Q_MOC_OUTPUT_REVISION) || (Q_MOC_OUTPUT_REVISION != 26)
#error "This file was generated using the moc from 3.3.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

const char *MainFrame::className() const
{
    return "MainFrame";
}

QMetaObject *MainFrame::metaObj = 0;
static QMetaObjectCleanUp cleanUp_MainFrame( "MainFrame", &MainFrame::staticMetaObject );

#ifndef QT_NO_TRANSLATION
QString MainFrame::tr( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "MainFrame", s, c, QApplication::DefaultCodec );
    else
	return QString::fromLatin1( s );
}
#ifndef QT_NO_TRANSLATION_UTF8
QString MainFrame::trUtf8( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "MainFrame", s, c, QApplication::UnicodeUTF8 );
    else
	return QString::fromUtf8( s );
}
#endif // QT_NO_TRANSLATION_UTF8

#endif // QT_NO_TRANSLATION

QMetaObject* MainFrame::staticMetaObject()
{
    if ( metaObj )
	return metaObj;
    QMetaObject* parentObject = QMainWindow::staticMetaObject();
    static const QUMethod slot_0 = {"fileNew", 0, 0 };
    static const QUMethod slot_1 = {"fileOpen", 0, 0 };
    static const QUMethod slot_2 = {"fileSave", 0, 0 };
    static const QUMethod slot_3 = {"fileSaveAs", 0, 0 };
    static const QUMethod slot_4 = {"fileExit", 0, 0 };
    static const QUMethod slot_5 = {"helpIndex", 0, 0 };
    static const QUMethod slot_6 = {"helpContents", 0, 0 };
    static const QUMethod slot_7 = {"helpAbout", 0, 0 };
    static const QUMethod slot_8 = {"setWire", 0, 0 };
    static const QUMethod slot_9 = {"languageChange", 0, 0 };
    static const QMetaData slot_tbl[] = {
	{ "fileNew()", &slot_0, QMetaData::Public },
	{ "fileOpen()", &slot_1, QMetaData::Public },
	{ "fileSave()", &slot_2, QMetaData::Public },
	{ "fileSaveAs()", &slot_3, QMetaData::Public },
	{ "fileExit()", &slot_4, QMetaData::Public },
	{ "helpIndex()", &slot_5, QMetaData::Public },
	{ "helpContents()", &slot_6, QMetaData::Public },
	{ "helpAbout()", &slot_7, QMetaData::Public },
	{ "setWire()", &slot_8, QMetaData::Public },
	{ "languageChange()", &slot_9, QMetaData::Protected }
    };
    metaObj = QMetaObject::new_metaobject(
	"MainFrame", parentObject,
	slot_tbl, 10,
	0, 0,
#ifndef QT_NO_PROPERTIES
	0, 0,
	0, 0,
#endif // QT_NO_PROPERTIES
	0, 0 );
    cleanUp_MainFrame.setMetaObject( metaObj );
    return metaObj;
}

void* MainFrame::qt_cast( const char* clname )
{
    if ( !qstrcmp( clname, "MainFrame" ) )
	return this;
    return QMainWindow::qt_cast( clname );
}

bool MainFrame::qt_invoke( int _id, QUObject* _o )
{
    switch ( _id - staticMetaObject()->slotOffset() ) {
    case 0: fileNew(); break;
    case 1: fileOpen(); break;
    case 2: fileSave(); break;
    case 3: fileSaveAs(); break;
    case 4: fileExit(); break;
    case 5: helpIndex(); break;
    case 6: helpContents(); break;
    case 7: helpAbout(); break;
    case 8: setWire(); break;
    case 9: languageChange(); break;
    default:
	return QMainWindow::qt_invoke( _id, _o );
    }
    return TRUE;
}

bool MainFrame::qt_emit( int _id, QUObject* _o )
{
    return QMainWindow::qt_emit(_id,_o);
}
#ifndef QT_NO_PROPERTIES

bool MainFrame::qt_property( int id, int f, QVariant* v)
{
    return QMainWindow::qt_property( id, f, v);
}

bool MainFrame::qt_static_property( QObject* , int , int , QVariant* ){ return FALSE; }
#endif // QT_NO_PROPERTIES
