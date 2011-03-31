TARGET = trimesh_intersection
DEPENDPATH += . 
INCLUDEPATH += . ../../..
CONFIG += console stl 
TEMPLATE = app
SOURCES += trimesh_intersection.cpp ../../../wrap/ply/plylib.cpp
