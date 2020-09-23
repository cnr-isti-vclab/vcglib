/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
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
****************************************************************************/

#ifndef __VCGLIB_TETRA_IMPORT
#define __VCGLIB_TETRA_IMPORT

#include <wrap/io_tetramesh/import_ply.h>
#include <wrap/io_tetramesh/import_msh.h>

#include <locale>

namespace vcg {
namespace tetra {
namespace io {

/**
This class encapsulate a filter for automatically importing meshes by guessing
the right filter according to the extension
*/

template <class OpenMeshType>
class Importer
{
private:
  enum KnownTypes { KT_UNKNOWN, KT_PLY, KT_MSH };

static int &LastType()
{
  static int lastType= KT_UNKNOWN;
return lastType;
}

public:
enum ImporterError {
  E_NOERROR =0 // No error =0 is the standard for ALL the imported files.
};
// simple aux function that returns true if a given file has a given extesnion
static bool FileExtension(std::string filename, std::string extension)
{
  std::transform(filename.begin(),   filename.end(),  filename.begin(), ::tolower);
  std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);

  std::string end=filename.substr(filename.length() - extension.length(), extension.length());
  
  return end == extension;
}
// Open Mesh, returns 0 on success.
static int Open(OpenMeshType &m, const char *filename, CallBackPos *cb=0)
{
  int dummymask = 0;
  return Open(m, filename, dummymask, cb);
}

/// Open Mesh and fills the load mask (the load mask must be initialized first); returns 0 on success.
static int Open(OpenMeshType &m, const char *filename, int &loadmask, CallBackPos *cb=0)
{
	int err;
	if(FileExtension(filename,"ply"))
	{
                err = tetra::io::ImporterPLY<OpenMeshType>::Open(m, filename, loadmask, cb);
		LastType()=KT_PLY;
	}
	else if(FileExtension(filename,"msh"))
	{
                err = tetra::io::ImporterMSH<OpenMeshType>::Open(m, filename, cb);
		LastType()=KT_MSH;
	}
    else
    {
		err=1;
		LastType()=KT_UNKNOWN;
	}

	return err;
}

static bool ErrorCritical(int error)
{
  switch(LastType())
  {
    case KT_PLY : return (error > 0); break;
    case KT_MSH : return (error > 0); break;
  }

  return true;
}

static const char *ErrorMsg(int error)
{
  switch(LastType())
  {
    // case KT_PLY : return ImporterPLY<OpenMeshType>::ErrorMsg(error); break;
    // case KT_STL : return ImporterSTL<OpenMeshType>::ErrorMsg(error); break;
    // case KT_OFF : return ImporterOFF<OpenMeshType>::ErrorMsg(error); break;
    // case KT_OBJ : return ImporterOBJ<OpenMeshType>::ErrorMsg(error); break;
    // case KT_VMI : return ImporterVMI<OpenMeshType>::ErrorMsg(error); break;
  }
  return "Unknown type";
}

static bool LoadMask(const char * filename, int &mask)
{
	bool err;

	if(FileExtension(filename, "ply"))
	{
                err = tetra::io::ImporterPLY<OpenMeshType>::LoadMask(filename, mask);
		LastType() = KT_PLY;
	}
	else if(FileExtension(filename, "msh"))
	{
		mask = Mask::IOM_VERTCOORD | Mask::IOM_TETRAINDEX;
		err = true;
		LastType() = KT_MSH;
	}
	else
	{
		err = false;
		LastType()=KT_UNKNOWN;
	}

	return err;
}
}; // end class
} // end Namespace tetra
} // end Namespace io
} // end Namespace vcg

#endif
