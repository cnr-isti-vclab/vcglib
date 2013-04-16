#DEFINES += VCG_USE_EIGEN

TARGET = trimesh_voroatlas
DEPENDPATH += . ../../..

INCLUDEPATH += . ../../..

CONFIG += console stl
TEMPLATE = app
SOURCES += trimesh_voroatlas.cpp ../../../wrap/ply/plylib.cpp
HEADERS += ../../../vcg/complex/algorithms/parametrization/voronoi_atlas.h

# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle


