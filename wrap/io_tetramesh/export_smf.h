/***************************************************************************** VCGLib                                                            o o     ** Visual and Computer Graphics Library                            o     o   **                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    ** Visual Computing Lab                                            /\/|      ** ISTI - Italian National Research Council                           |      **                                                                    \      ** All rights reserved.                                                      **                                                                           ** This program is free software; you can redistribute it and/or modify      *   * it under the terms of the GNU General Public License as published by      ** the Free Software Foundation; either version 2 of the License, or         ** (at your option) any later version.                                       **                                                                           ** This program is distributed in the hope that it will be useful,           ** but WITHOUT ANY WARRANTY; without even the implied warranty of            ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          ** for more details.                                                         **                                                                           *****************************************************************************/
/****************************************************************************
  History

$Log: not supported by cvs2svn $
Revision 1.3  2004/05/10 13:14:28  ganovelli
converted to library style (namespaces etc..)


****************************************************************************/#ifndef __VCGLIB_IMPORTERSMF
#define __VCGLIB_IMPORTERSMF#include <vcg/space/point3.h>

namespace vcg {
namespace tetra {
namespace io {

template <typename  MESHTYPE>
class ExporterTS{

	typedef typename MESHTYPE::VertexPointer VertexPointer;
	typedef typename MESHTYPE::VertexType VertexType;
	typedef typename MESHTYPE::TetraType FaceType;
	typedef typename MESHTYPE::VertexIterator VertexIterator;
	typedef typename MESHTYPE::FaceIterator FaceIterator;
	typedef MESHTYPE::ScalarType ScalarType;
	typedef Point3<ScalarType> Point3x;


int Save( MESHTYPE & m, const char * filename ){	
	FILE *f;
	f = fopen(filename,"w");
	if(f == NULL ) 
		{
			printf( "The file could not be opened\n" );
			return -1;
		}
   else
   {
		fprintf(f, "%i", m.vn );
		fprintf(f, "%i", m.tn );
		VertexIterator vi;
		for (vi = m.vert.begin(); vi != m.vert.end();++vi)
			fprintf(f, "%f %f %f \n", (*vi).P()[0],(*vi).P()[1],(*vi).P()[2] );

		TetraIterator ti;
		fr( ti = m.tetra.begin(); ti != m.tetra.end(); ++ti)
			fprintf(f, "%d %d %d %d \n",
							(*ti)->V(0)-&*m.tetra.begin(),
							(*ti)->V(1)-&*m.tetra.begin(),
							(*ti)->V(2)-&*m.tetra.begin(),
							(*ti)->V(3)-&*m.tetra.begin());
	 } return 0;}	};// end class		};// end of io	};// end of tri};// end of vcg#endif