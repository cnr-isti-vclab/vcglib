TARGET = trimesh_edge
DEPENDPATH += .  ../../..
INCLUDEPATH += . ../../..
CONFIG += console stl opengl
TEMPLATE = app
SOURCES += trimesh_edge.cpp ../../../wrap/ply/plylib.cpp
HEADERS +=  ../../../vcg/complex/algorithms/update/topology.h
HEADERS +=  ../../../wrap/gl/glu_tessellator_cap.h
# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle
