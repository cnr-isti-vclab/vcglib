/****************************************************************************
** SimpleGLWidget meta object code from reading C++ file 'simpleglwidget.h'
**
** Created: Sat Dec 18 11:09:46 2004
**      by: The Qt MOC ($Id: moc_simpleglwidget.cpp,v 1.2 2004-12-20 17:56:01 pietroni Exp $)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#undef QT_NO_COMPAT
#include "simpleglwidget.h"
#include <qmetaobject.h>
#include <qapplication.h>

#include <private/qucomextra_p.h>
#if !defined(Q_MOC_OUTPUT_REVISION) || (Q_MOC_OUTPUT_REVISION != 26)
#error "This file was generated using the moc from 3.3.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

const char *SimpleGLWidget::className() const
{
    return "SimpleGLWidget";
}

QMetaObject *SimpleGLWidget::metaObj = 0;
static QMetaObjectCleanUp cleanUp_SimpleGLWidget( "SimpleGLWidget", &SimpleGLWidget::staticMetaObject );

#ifndef QT_NO_TRANSLATION
QString SimpleGLWidget::tr( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "SimpleGLWidget", s, c, QApplication::DefaultCodec );
    else
	return QString::fromLatin1( s );
}
#ifndef QT_NO_TRANSLATION_UTF8
QString SimpleGLWidget::trUtf8( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "SimpleGLWidget", s, c, QApplication::UnicodeUTF8 );
    else
	return QString::fromUtf8( s );
}
#endif // QT_NO_TRANSLATION_UTF8

#endif // QT_NO_TRANSLATION

QMetaObject* SimpleGLWidget::staticMetaObject()
{
    if ( metaObj )
	return metaObj;
    QMetaObject* parentObject = QGLWidget::staticMetaObject();
    static const QUMethod slot_0 = {"Open", 0, 0 };
    static const QUMethod slot_1 = {"ShowSlides", 0, 0 };
    static const QUMethod slot_2 = {"SetWire", 0, 0 };
    static const QUMethod slot_3 = {"SetShowBlocked", 0, 0 };
    static const QUMethod slot_4 = {"ShowExternalForces", 0, 0 };
    static const QUMethod slot_5 = {"ShowInternalForces", 0, 0 };
    static const QUMethod slot_6 = {"ShowResultForces", 0, 0 };
    static const QUMethod slot_7 = {"Smooth", 0, 0 };
    static const QUMethod slot_8 = {"SavePly", 0, 0 };
    static const QUMethod slot_9 = {"Apply", 0, 0 };
    static const QUMethod slot_10 = {"Extract", 0, 0 };
    static const QUMethod slot_11 = {"Update", 0, 0 };
    static const QUMethod slot_12 = {"Clear", 0, 0 };
    static const QMetaData slot_tbl[] = {
	{ "Open()", &slot_0, QMetaData::Public },
	{ "ShowSlides()", &slot_1, QMetaData::Public },
	{ "SetWire()", &slot_2, QMetaData::Public },
	{ "SetShowBlocked()", &slot_3, QMetaData::Public },
	{ "ShowExternalForces()", &slot_4, QMetaData::Public },
	{ "ShowInternalForces()", &slot_5, QMetaData::Public },
	{ "ShowResultForces()", &slot_6, QMetaData::Public },
	{ "Smooth()", &slot_7, QMetaData::Public },
	{ "SavePly()", &slot_8, QMetaData::Public },
	{ "Apply()", &slot_9, QMetaData::Public },
	{ "Extract()", &slot_10, QMetaData::Public },
	{ "Update()", &slot_11, QMetaData::Public },
	{ "Clear()", &slot_12, QMetaData::Public }
    };
    metaObj = QMetaObject::new_metaobject(
	"SimpleGLWidget", parentObject,
	slot_tbl, 13,
	0, 0,
#ifndef QT_NO_PROPERTIES
	0, 0,
	0, 0,
#endif // QT_NO_PROPERTIES
	0, 0 );
    cleanUp_SimpleGLWidget.setMetaObject( metaObj );
    return metaObj;
}

void* SimpleGLWidget::qt_cast( const char* clname )
{
    if ( !qstrcmp( clname, "SimpleGLWidget" ) )
	return this;
    return QGLWidget::qt_cast( clname );
}

bool SimpleGLWidget::qt_invoke( int _id, QUObject* _o )
{
    switch ( _id - staticMetaObject()->slotOffset() ) {
    case 0: Open(); break;
    case 1: ShowSlides(); break;
    case 2: SetWire(); break;
    case 3: SetShowBlocked(); break;
    case 4: ShowExternalForces(); break;
    case 5: ShowInternalForces(); break;
    case 6: ShowResultForces(); break;
    case 7: Smooth(); break;
    case 8: SavePly(); break;
    case 9: Apply(); break;
    case 10: Extract(); break;
    case 11: Update(); break;
    case 12: Clear(); break;
    default:
	return QGLWidget::qt_invoke( _id, _o );
    }
    return TRUE;
}

bool SimpleGLWidget::qt_emit( int _id, QUObject* _o )
{
    return QGLWidget::qt_emit(_id,_o);
}
#ifndef QT_NO_PROPERTIES

bool SimpleGLWidget::qt_property( int id, int f, QVariant* v)
{
    return QGLWidget::qt_property( id, f, v);
}

bool SimpleGLWidget::qt_static_property( QObject* , int , int , QVariant* ){ return FALSE; }
#endif // QT_NO_PROPERTIES
