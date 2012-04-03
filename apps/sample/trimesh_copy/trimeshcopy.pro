
TARGET = trimeshcopy
DEPENDPATH += ../..
INCLUDEPATH += . ../..
CONFIG += console stl debug_and_release
TEMPLATE = app
HEADERS += 
SOURCES += trimeshcopy.cpp ../../wrap/ply/plylib.cpp


# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle
