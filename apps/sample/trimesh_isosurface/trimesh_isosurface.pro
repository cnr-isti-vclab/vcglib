INCLUDEPATH += . ../../..
CONFIG += console stl
TEMPLATE = app
# Input
HEADERS += trivial_walker.h simple_volume.h
SOURCES += trimesh_isosurface.cpp ../../../wrap/ply/plylib.cpp
