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
class ImporterTS{

	typedef typename MESHTYPE::VertexPointer VertexPointer;
	typedef typename MESHTYPE::VertexType VertexType;
	typedef typename MESHTYPE::TetraType FaceType;
	typedef typename MESHTYPE::VertexIterator VertexIterator;
	typedef typename MESHTYPE::FaceIterator FaceIterator;
	typedef MESHTYPE::ScalarType ScalarType;
	typedef Point3<ScalarType> Point3x;


int Load( MESHTYPE & m, const char * filename ){	
	int nvertex;
	int ntetra;
	float x;
	float y;
	float z;
	int tp0;
	int tp1;
	int tp2;
	int tp3;
	float mass;
	FILE *f;
	Tetramesh::VertexType p1;
	f = fopen(filename,"r");
	if(f == NULL ) 
		{
			printf( "The file was not opened\n" );
			return -1;
		}
   else
   {
		fscanf(f, "%i", &nvertex );
		fscanf(f, "%i", &ntetra );
		int j;
		for (j=0;j<nvertex;j++)
		{
			fscanf(f, "%f", &x );
			fscanf(f, "%f", &y );
			fscanf(f, "%f", &z );
		  //fscanf(f, "%f", &mass );
      p1.ClearFlags();
      p1.P()=Point3x(x, y,z );
			vert.push_back(p1);
		}
		tetra.reserve(ntetra*10);
    vert.reserve(nvertex*10);
		for (j=0;j<ntetra;j++)
		{
			fscanf(f, "%i", &tp0 );
			fscanf(f, "%i", &tp1 );
			fscanf(f, "%i", &tp2 );
			fscanf(f, "%i", &tp3 );
			
			Tetramesh::TetraType  newTetra;
			tetra.push_back(newTetra);
			tetra.back().Init(&vert[tp0],&vert[tp1],&vert[tp2],&vert[tp3]); 
		}
	 }
		return 0;	 }	};// end class		};// end of io	};// end of tri};// end of vcg#endif