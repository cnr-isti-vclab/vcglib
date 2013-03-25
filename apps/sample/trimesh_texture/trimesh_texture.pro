include(../common.pri)
QT += opengl svg
TARGET = trimesh_texture
SOURCES += trimesh_texture.cpp  ../../../../vcglib/wrap/qt/PolyToQImage.cpp ../../../wrap/ply/plylib.cpp
