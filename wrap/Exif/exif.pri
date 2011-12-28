#Please define the EXIF_DIR variable which locates the glUtils directory in your system.
#eg. EXIF_DIR = ../sandbox/brivio/Exif

!contains(DEFINES, EXIF_DIR){
  DEFINES += EXIF_DIR

  INCLUDEPATH += $$EXIF_DIR/include
  SOURCES += $$EXIF_DIR/src/exif.cpp
  SOURCES += $$EXIF_DIR/src/gpsinfo.cpp
  SOURCES += $$EXIF_DIR/src/iptc.cpp
  SOURCES += $$EXIF_DIR/src/jhead.cpp
  SOURCES += $$EXIF_DIR/src/jpgfile.cpp
  SOURCES += $$EXIF_DIR/src/makernote.cpp
  win32:SOURCES += $$EXIF_DIR/src/myglob.cpp
  SOURCES += $$EXIF_DIR/src/paths.cpp
}
