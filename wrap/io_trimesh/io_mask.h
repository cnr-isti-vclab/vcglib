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
Revision 1.3  2004/05/12 10:19:30  ganovelli
new line added at the end of file

Revision 1.2  2004/03/09 21:26:47  cignoni
cr lf mismatch

Revision 1.1  2004/03/03 15:00:51  cignoni
Initial commit

****************************************************************************/
#ifndef __VCGLIB_IOTRIMESH_IO_MASK
#define __VCGLIB_IOTRIMESH_IO_MASK

#include<wrap/callback.h>
#include<wrap/ply/plylib.h>

namespace vcg {
namespace tri {
namespace io {

/**
@name Load and Save in Ply format
*/
//@{
  
class PLYMask
{
public:

  /*
  Bitmask for specifying what data has to be loaded or saved or it is present in a given plyfile;
*/

enum {
	PM_NONE         = 0x0000,

  PM_VERTCOORD    = 0x0001,
	PM_VERTFLAGS    = 0x0002, 
	PM_VERTCOLOR    = 0x0004,
	PM_VERTQUALITY  = 0x0008,
	PM_VERTNORMAL   = 0x0010,
	PM_VERTTEXCOORD = 0x0020,

	PM_FACEINDEX    = 0x0040,
	PM_FACEFLAGS    = 0x0080,
	PM_FACECOLOR    = 0x0100,
	PM_FACEQUALITY  = 0x0200,
	PM_FACENORMAL   = 0x0400,
	PM_WEDGCOLOR    = 0x0800,
	PM_WEDGTEXCOORD = 0x1000,
	PM_WEDGTEXMULTI = 0x2000, // Se ha anche l'indice di texture esplicito
	PM_WEDGNORMAL   = 0x4000,

	PM_CAMERA       = 0x8000,

	PM_FLAGS        = PM_VERTFLAGS + PM_FACEFLAGS,

	PM_ALL          = 0xFFFF
};


static void SMFlags2String( int mask, char str[] )
{
	str[0] = 0;

	strcat(str,"V:");
	if( mask & PM_VERTFLAGS    ) strcat(str,"flag,");
	if( mask & PM_VERTCOLOR    ) strcat(str,"color,");
	if( mask & PM_VERTQUALITY  ) strcat(str,"quality,");
	if( mask & PM_VERTTEXCOORD ) strcat(str,"tcoord,");
	if( mask & PM_VERTNORMAL ) strcat(str,"normal,");

	strcat(str," F:");
	if( mask & PM_FACEFLAGS    ) strcat(str,"mask,");
	if( mask & PM_FACECOLOR    ) strcat(str,"color,");
	if( mask & PM_FACEQUALITY  ) strcat(str,"quality,");
	if( mask & PM_FACENORMAL   ) strcat(str,"normal,");

	strcat(str," W:");
	if( mask & PM_WEDGCOLOR    ) strcat(str,"color,");
	if( mask & PM_WEDGTEXCOORD ) strcat(str,"tcoord,");
	if( mask & PM_WEDGNORMAL  ) strcat(str,"normal,");

	if( mask & PM_CAMERA ) strcat(str," camera");
}

}; // end class
//@}
} // end namespace tri
} // end namespace io
} // end namespace vcg
#endif
