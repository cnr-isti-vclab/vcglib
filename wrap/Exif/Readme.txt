How to read exif metadata from a jpg file:
- add "echo INCLUDEPATH += %rootrel%/code/lib/exif/include>> %filename%"
  and "echo SOURCES += %rootrel%\code\lib\exif\*.cpp>> %filename%" to pmake.bat
- include "jhead.h";
- invoke ::ProcessFile(filename);
- get the available information from the "ImageInfo" global variable of type ImageInfo_t, which is defined in jhead.h.
