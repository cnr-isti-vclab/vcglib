INCLUDEPATH += . ../../.. ../../../../code/lib
HEADERS       = glwidget.h \
                window.h \
		mesh_type	
SOURCES       = glwidget.cpp \
                main.cpp \
                window.cpp\
		../../../wrap/ply/plylib.cpp\
		../../../wrap/gui/trackball.cpp\
		../../../wrap/gui/trackmode.cpp\
		../../../wrap/gui/trackball.cpp
QT           += opengl

# install
target.path = $$./debug
sources.files = $$SOURCES $$HEADERS $$RESOURCES $$FORMS trimesh_pos_demo.pro
sources.path =  ./
INSTALLS += target sources

