
TARGET = trimeshcopy
DEPENDPATH += ../../..
INCLUDEPATH += . ../../..
CONFIG += console stl
TEMPLATE = app
HEADERS += 
SOURCES += trimeshcopy.cpp ../../../wrap/ply/plylib.cpp

#DEFINES += N_DEBUG
# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle
