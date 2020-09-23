include(../common.pri)
TARGET = trimesh_import_fbx
SOURCES += trimesh_import_fbx.cpp  \
  ../../../wrap/openfbx/src/ofbx.cpp \
  ../../../wrap/openfbx/src/miniz.c

CONFIG += c++11

