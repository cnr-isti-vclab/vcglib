DEPENDPATH += . ../../../../..
INCLUDEPATH += . ../../../.. ../../../../eigenlib

CONFIG += console c++11
TEMPLATE = app
# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle

QMAKE_CXXFLAGS += -std=c++11
TARGET = trimesh_field_smoothing
SOURCES += trimesh_field_smoothing.cpp ../../../ply/plylib.cpp
