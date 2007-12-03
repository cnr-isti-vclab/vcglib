
TARGET = tridecimator
LIBPATH += 
DEPENDPATH += . 
INCLUDEPATH += . ../..
CONFIG += console stl debug_and_release
TEMPLATE = app
HEADERS += 
SOURCES += tridecimator.cpp ../../wrap/ply/plylib.cpp


# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle