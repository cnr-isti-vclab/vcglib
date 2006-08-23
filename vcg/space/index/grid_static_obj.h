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
Revision 1.2  2005/12/02 00:29:00  cignoni
updated the templates of BasicGrid

Revision 1.1  2005/07/28 08:41:00  cignoni
First working version

Revision 1.13  2005/04/14 17:23:08  ponchio
*** empty log message ***

****************************************************************************/

#ifndef __VCGLIB_UGRID_OBJ
#define __VCGLIB_UGRID_OBJ

#include <vector>
#include <algorithm>
#include <stdio.h>

#include <vcg/space/box3.h>
#include <vcg/space/line3.h>
#include <vcg/space/index/grid_util.h>
namespace vcg {
/** Static Uniform Grid
A simple Spatial grid of object. 
Kept in the most trivial way. Every cell is allocated 
and contains one istance of the template class.
*/

template < typename class ObjType, class FLT=float  >
class GridStaticObj : public BasicGrid<ObjType, FLT>
{
 public:

	 /// La matriciona della griglia
	 ObjType *grid;

	 int size() const { return siz[0]*siz[1]*siz[2];}

	 inline  GridStaticObj() { grid = 0; }
	 inline ~GridStaticObj() { if(grid) delete[] grid; }
	 inline Init(const ObjType &val)
	 {
		 fill(grid,grid+size(),val);
	 }


	 /// Date le coordinate ritorna la cella
	 inline ObjType & Grid( const int x, const int y, const int z ) {return grid[GridInd(Point3i(x,y,z))]; }

	 // Dato un punto ritorna la cella  
	 inline ObjType & Grid( const Point3<FLT> & p )                 {				return grid[GridInd(p)];		}

	 inline int GridInd( const Point3i & pi ) const
	 {
#ifndef NDEBUG
		 if ( pi[0]<0 || pi[0]>=siz[0] || pi[1]<0 || pi[1]>=siz[1] || pi[2]<0 || pi[2]>=siz[2] )
		 {	assert(0);
		 return 0;
		 } 
#endif
		 return pi[0]+siz[0]*(pi[1]+siz[1]*pi[2]);
	 }

	 // Dato un punto ritorna l'indice della cella
	 inline int GridInd( const Point3<FLT> & p ) const { return GridInd(GridP(p)); 	}
  
	void Create( Point3i &_siz, const ObjType & init )
	{
		siz=_siz;
	 	voxel[0] = dim[0]/siz[0];
		voxel[1] = dim[1]/siz[1];
		voxel[2] = dim[2]/siz[2];

		if(grid) delete[] grid;
		int n = siz[0]*siz[1]*siz[2];
		grid = new ObjType[n];
		fill(grid,grid+n,init);
	}

	/// Crea una griglia di un dato bbox e con un certo numero di elem.
	/// il bbox viene gonfiato appositamente.

	template<class FLT2>
	void Create(const Box3<FLT2> & b, int ncell, const ObjType & init, bool Inflate = true )
	{
		bbox.Import(b);
		if(Inflate) bbox.Offset(0.01*bbox.Diag());
		dim  = bbox.max - bbox.min;

		// Calcola la dimensione della griglia
		Point3i _siz;
		BestDim( ncell, dim, _siz );
		Create(_siz, init );
	}
};
//end class SGrid


}
#endif
