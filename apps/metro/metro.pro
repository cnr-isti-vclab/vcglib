
TARGET = metro
LIBPATH += 
DEPENDPATH += . 
INCLUDEPATH += . ../..
CONFIG += console stl
TEMPLATE = app
HEADERS += defs.h sampling.h mesh_type.h
SOURCES += metro.cpp ../../wrap/ply/plylib.cpp