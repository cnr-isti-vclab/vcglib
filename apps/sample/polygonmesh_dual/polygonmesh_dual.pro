\macx: QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.9
QMAKE_MAC_SDK = macosx10.9
CONFIG += c++11

include(../common.pri)
TARGET = polygonmesh_base
SOURCES += polygonmesh_dual.cpp ../../../wrap/ply/plylib.cpp
