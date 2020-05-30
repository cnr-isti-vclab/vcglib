include(../common.pri)

TARGET = trimesh_align_pair

HEADERS += \
	my_mesh.h

SOURCES += \
	trimesh_align_pair.cpp \
	../../../wrap/ply/plylib.cpp
