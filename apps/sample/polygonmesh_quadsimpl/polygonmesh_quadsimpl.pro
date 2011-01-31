TARGET = polygonmesh_quadsimpl
LIBPATH += 
DEPENDPATH += . 
INCLUDEPATH += . ../../..
CONFIG += console stl 
TEMPLATE = app
SOURCES += ../../../wrap/ply/plylib.cpp \
    quadsimpl.cpp
