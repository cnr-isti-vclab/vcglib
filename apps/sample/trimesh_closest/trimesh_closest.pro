
TARGET = trimesh_closest
DEPENDPATH += ../../..
INCLUDEPATH += . ../../..
CONFIG += console stl
TEMPLATE = app
HEADERS += 
SOURCES += trimesh_closest.cpp ../../../wrap/ply/plylib.cpp

#DEFINES += N_DEBUG
# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle
