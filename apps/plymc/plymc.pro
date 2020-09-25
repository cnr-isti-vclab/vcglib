TARGET = plymc
DEPENDPATH += .
INCLUDEPATH += ../.. 
CONFIG += console c++17
CONFIG -= app_bundle
TEMPLATE = app

SOURCES += ../../wrap/ply/plylib.cpp \
    plymc_main.cpp

