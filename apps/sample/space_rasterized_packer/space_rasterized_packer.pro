include(../common.pri)
QT += opengl svg
TARGET = space_rasterized_packer
SOURCES += space_rasterized_packer.cpp  \
../../../../vcglib/wrap/qt/Outline2ToQImage.cpp \
../../../../vcglib/wrap/qt/outline2_rasterizer.cpp
