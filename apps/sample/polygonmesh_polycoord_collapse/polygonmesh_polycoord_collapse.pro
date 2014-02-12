TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp


VCGLIBDIR = ../../../
DEPENDPATH += $$VCGLIBDIR
INCLUDEPATH += $$VCGLIBDIR
SOURCES += $$VCGLIBDIR/wrap/ply/plylib.cpp
