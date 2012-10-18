VCGLIBDIR = ../../../vcg/vcglib
GLEWDIR   = ../../../code/lib/glew
ANTDIR    = ../../../code/lib/AntTweakBar1.14
COMISODIR = ./CoMISo

# Glew settings
DEFINES += GLEW_STATIC
INCLUDEPATH += $$GLEWDIR/include
SOURCES += $$GLEWDIR/src/glew.c

#Anttweakbar stuff
mac{
  LIBS +=$$ANTDIR/lib/libAntTweakBar.dylib
}

#Comiso
mac{
LIBS += -L $$COMISODIR/build/Build/lib/CoMISo/ -lCoMISo
}

mac{
  QMAKE_POST_LINK +="cp -P ../../../code/lib/AntTweakBar1.14/lib/libAntTweakBar.dylib . ; "
  QMAKE_POST_LINK +="install_name_tool -change ../lib/libAntTweakBar.dylib ./libAntTweakBar.dylib $$TARGET ; "
  QMAKE_POST_LINK +="cp -P $$COMISODIR/build/Build/lib/CoMISo/libCoMISo.dylib . ; "
}
