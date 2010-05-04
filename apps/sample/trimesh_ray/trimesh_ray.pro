TARGET = trimesh_ray
LIBPATH += 
DEPENDPATH += . 
INCLUDEPATH += . ../../..
CONFIG += console stl 
TEMPLATE = app
SOURCES += trimesh_ray.cpp ../../../wrap/ply/plylib.cpp

# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle
