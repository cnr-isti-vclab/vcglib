/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *   
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

/****************************************************************************
  History

$Log: not supported by cvs2svn $
Revision 1.2  2005/01/26 22:43:19  cignoni
Add std:: to stl containers

Revision 1.1  2004/11/29 08:12:10  cignoni
Initial Update


****************************************************************************/

#ifndef __VCGLIB_IMPORT
#define __VCGLIB_IMPORT

#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/import_stl.h>
#include <wrap/io_trimesh/import_off.h>

#include <locale>

namespace vcg {
namespace tri {
namespace io {

/** 
This class encapsulate a filter for automatically importing meshes by guessing 
the right filter according to the extension
*/

template <class OpenMeshType>
class Importer
{
private:
  enum KnownTypes { KT_UNKNOWN, KT_PLY, KT_STL, KT_OFF };
static int &LastType()
{
  static int lastType= KT_UNKNOWN;
return lastType;
}

public:
// simple aux function that returns true if a given file has a given extesnion
static bool FileExtension(std::string filename,  std::string extension)
{
  std::locale loc1 ;
  std::use_facet<std::ctype<char> > ( loc1 ).tolower(&*filename.begin(),&*filename.end());
  std::use_facet<std::ctype<char> > ( loc1 ).tolower(&*extension.begin(),&*extension.end());
  std::string end=filename.substr(filename.length()-extension.length(),extension.length());
  return end==extension;
}

// Open Mesh
static int Open(OpenMeshType &m, const char *filename, CallBackPos *cb=0)
{
  int err;
  if(FileExtension(filename,"ply"))
  {
    err = ImporterPLY<OpenMeshType>::Open(m,filename,cb);
    LastType()=KT_PLY;
  }
  else if(FileExtension(filename,"stl"))
  {
    err = ImporterSTL<OpenMeshType>::Open(m,filename,cb);
    LastType()=KT_STL;
  }
   else if(FileExtension(filename,"off"))
  {
    err = ImporterOFF<OpenMeshType>::Open(m,filename,cb);
    LastType()=KT_OFF;
  }
 else {
    err=1;
    LastType()=KT_UNKNOWN;
  }

  return err;
}

static const char *ErrorMsg(int error)
{
  switch(LastType())
  {
    case KT_PLY : return ImporterPLY<OpenMeshType>::ErrorMsg(error); break;
    case KT_STL : return ImporterSTL<OpenMeshType>::ErrorMsg(error); break;
    case KT_OFF : return ImporterOFF<OpenMeshType>::ErrorMsg(error); break;
  }
  return "Unknown type";  
}

}; // end class
} // end Namespace tri
} // end Namespace io
} // end Namespace vcg

#endif
