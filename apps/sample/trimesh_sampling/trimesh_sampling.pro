include(../common.pri)
TARGET = trimesh_sampling
SOURCES += trimesh_sampling.cpp

# Awful..
win32{
  DEFINES += NOMINMAX
}
