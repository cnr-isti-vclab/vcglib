/****************************************************************************
** MyGLWidget meta object code from reading C++ file 'myglwidget.h'
**
** Created: Mon Oct 4 19:00:57 2004
**      by: The Qt MOC ($Id: moc_myglwidget.cpp,v 1.1 2004-10-04 18:01:36 ganovelli Exp $)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#undef QT_NO_COMPAT
#include "myglwidget.h"
#include <qmetaobject.h>
#include <qapplication.h>

#include <private/qucomextra_p.h>
#if !defined(Q_MOC_OUTPUT_REVISION) || (Q_MOC_OUTPUT_REVISION != 26)
#error "This file was generated using the moc from 3.3.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

const char *MyGLWidget::className() const
{
    return "MyGLWidget";
}

QMetaObject *MyGLWidget::metaObj = 0;
static QMetaObjectCleanUp cleanUp_MyGLWidget( "MyGLWidget", &MyGLWidget::staticMetaObject );

#ifndef QT_NO_TRANSLATION
QString MyGLWidget::tr( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "MyGLWidget", s, c, QApplication::DefaultCodec );
    else
	return QString::fromLatin1( s );
}
#ifndef QT_NO_TRANSLATION_UTF8
QString MyGLWidget::trUtf8( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "MyGLWidget", s, c, QApplication::UnicodeUTF8 );
    else
	return QString::fromUtf8( s );
}
#endif // QT_NO_TRANSLATION_UTF8

#endif // QT_NO_TRANSLATION

QMetaObject* MyGLWidget::staticMetaObject()
{
    if ( metaObj )
	return metaObj;
    QMetaObject* parentObject = QGLWidget::staticMetaObject();
    static const QUMethod slot_0 = {"setBox", 0, 0 };
    static const QUMethod slot_1 = {"setWire", 0, 0 };
    static const QUMethod slot_2 = {"setHidden", 0, 0 };
    static const QUMethod slot_3 = {"setFlat", 0, 0 };
    static const QUMethod slot_4 = {"setFlatWire", 0, 0 };
    static const QUMethod slot_5 = {"setSmooth", 0, 0 };
    static const QUMethod slot_6 = {"setSmallTetra", 0, 0 };
    static const QUMethod slot_7 = {"TrackMouseModality", 0, 0 };
    static const QUMethod slot_8 = {"SectionMouseModality", 0, 0 };
    static const QUMethod slot_9 = {"SwitchTextSimplex", 0, 0 };
    static const QUMethod slot_10 = {"SwitchTextPhysics", 0, 0 };
    static const QUMethod slot_11 = {"SwitchTextQuality", 0, 0 };
    static const QMetaData slot_tbl[] = {
	{ "setBox()", &slot_0, QMetaData::Public },
	{ "setWire()", &slot_1, QMetaData::Public },
	{ "setHidden()", &slot_2, QMetaData::Public },
	{ "setFlat()", &slot_3, QMetaData::Public },
	{ "setFlatWire()", &slot_4, QMetaData::Public },
	{ "setSmooth()", &slot_5, QMetaData::Public },
	{ "setSmallTetra()", &slot_6, QMetaData::Public },
	{ "TrackMouseModality()", &slot_7, QMetaData::Public },
	{ "SectionMouseModality()", &slot_8, QMetaData::Public },
	{ "SwitchTextSimplex()", &slot_9, QMetaData::Public },
	{ "SwitchTextPhysics()", &slot_10, QMetaData::Public },
	{ "SwitchTextQuality()", &slot_11, QMetaData::Public }
    };
    metaObj = QMetaObject::new_metaobject(
	"MyGLWidget", parentObject,
	slot_tbl, 12,
	0, 0,
#ifndef QT_NO_PROPERTIES
	0, 0,
	0, 0,
#endif // QT_NO_PROPERTIES
	0, 0 );
    cleanUp_MyGLWidget.setMetaObject( metaObj );
    return metaObj;
}

void* MyGLWidget::qt_cast( const char* clname )
{
    if ( !qstrcmp( clname, "MyGLWidget" ) )
	return this;
    return QGLWidget::qt_cast( clname );
}

bool MyGLWidget::qt_invoke( int _id, QUObject* _o )
{
    switch ( _id - staticMetaObject()->slotOffset() ) {
    case 0: setBox(); break;
    case 1: setWire(); break;
    case 2: setHidden(); break;
    case 3: setFlat(); break;
    case 4: setFlatWire(); break;
    case 5: setSmooth(); break;
    case 6: setSmallTetra(); break;
    case 7: TrackMouseModality(); break;
    case 8: SectionMouseModality(); break;
    case 9: SwitchTextSimplex(); break;
    case 10: SwitchTextPhysics(); break;
    case 11: SwitchTextQuality(); break;
    default:
	return QGLWidget::qt_invoke( _id, _o );
    }
    return TRUE;
}

bool MyGLWidget::qt_emit( int _id, QUObject* _o )
{
    return QGLWidget::qt_emit(_id,_o);
}
#ifndef QT_NO_PROPERTIES

bool MyGLWidget::qt_property( int id, int f, QVariant* v)
{
    return QGLWidget::qt_property( id, f, v);
}

bool MyGLWidget::qt_static_property( QObject* , int , int , QVariant* ){ return FALSE; }
#endif // QT_NO_PROPERTIES
