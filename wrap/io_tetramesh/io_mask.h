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

#ifndef __VCGLIB_IOTETRAMESH_IO_MASK
#define __VCGLIB_IOTETRAMESH_IO_MASK

namespace vcg {
namespace tetra {
namespace io {

/**
@name Input/output data mask
*/
//@{
  
class Mask
{
public:

  /*
  Bitmask for specifying what data has to be loaded or saved or it is present in a given plyfile;
*/

enum {
	IOM_NONE         = 0x00000,

    IOM_VERTCOORD    = 0x00001,
	IOM_VERTFLAGS    = 0x00002,
	IOM_VERTCOLOR    = 0x00004,
	IOM_VERTQUALITY  = 0x00008,
	IOM_VERTRADIUS   = 0x10000,

	IOM_EDGEINDEX    = 0x80000,

	IOM_TETRAINDEX   = 0x00010,
	IOM_TETRAFLAGS   = 0x00020,
	IOM_TETRACOLOR   = 0x00040,
	IOM_TETRAQUALITY = 0x00080,

	IOM_CAMERA       = 0x08000,

	IOM_FLAGS        = IOM_VERTFLAGS + IOM_TETRAFLAGS,

	IOM_ALL          = 0xFFFFF
};

template <class MeshType>
static void ClampMask(MeshType &m, int &mask)
{
  if( (mask & IOM_TETRACOLOR)    && !HasPerTetraColor(m) )   mask = mask & (~IOM_TETRACOLOR);
  if( (mask & IOM_VERTCOLOR)    && !HasPerVertexColor(m) )   mask = mask & (~IOM_VERTCOLOR);
}

}; // end class
//@}

} // end namespace tetra
} // end namespace io
} // end namespace vcg
#endif
