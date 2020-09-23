
TARGET = metro
DEPENDPATH += ../..
INCLUDEPATH += . ../.. ../../eigenlib
CONFIG += console stl  c++11 debug_and_release
TEMPLATE = app
SOURCES += metro.cpp ../../wrap/ply/plylib.cpp

# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle
