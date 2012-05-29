QT += opengl svg
TARGET = space_packer
DEPENDPATH += . ../../..
INCLUDEPATH += . ../../..
CONFIG += console stl
TEMPLATE = app
SOURCES += space_packer.cpp     ../../../../vcglib/wrap/qt/PolyToQImage.cpp
HEADERS += ../../../vcg/space/rect_packer.h \
    ../../../vcg/math/similarity2.h \
    ../../../vcg/space/poly_packer.h


# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle

# Awful problem with windows..
win32{
  DEFINES += NOMINMAX
}

QT           += opengl svg