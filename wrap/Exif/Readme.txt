How to read exif metadata from a jpg file.
 = = = = = = = = = = = = = = = = = = = = =

1- In your .pro file, define the EXIF_DIR variable which locates the Exif directory in your system, and include the exif.pri file
e.g.
  EXIF_DIR = ./yourSandboxPath/Exif
  include($$EXIF_DIR/exif.pri)

2- In your source file (where you need the exif reader), add:
  #include "jhead.h"
Then, invoke:
  ::ProcessFile(filename);
and get the available exif information from the "ImageInfo" global variable of type ImageInfo_t, defined in jhead.h.
