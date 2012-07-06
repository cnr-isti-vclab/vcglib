
TARGET = trimesh_closest
DEPENDPATH += ../../..
INCLUDEPATH += . ../../..
CONFIG += console stl
TEMPLATE = app
HEADERS += 
SOURCES += trimesh_closest.cpp ../../../wrap/ply/plylib.cpp

release {DEFINES += NDEBUG}
# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle
