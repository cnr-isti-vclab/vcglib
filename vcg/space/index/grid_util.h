/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005                                                \/)\/    *
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

****************************************************************************/
#ifndef __VCGLIB_GRID_UTIL
#define __VCGLIB_GRID_UTIL


 /** Calcolo dimensioni griglia.
      Calcola la dimensione della griglia in funzione
      della ratio del bounding box e del numero di elementi
  */
 namespace vcg {
   
   template<class scalar_type>
   void BestDim( const int elems, const Point3<scalar_type> & size, Point3i & dim )
    {
      const int mincells   = 1;		// Numero minimo di celle
      const double GFactor = 1.0;	// GridEntry = NumElem*GFactor
      double diag = size.Norm();	// Diagonale del box
      double eps  = diag*1e-4;		// Fattore di tolleranza
      
      assert(elems>0);
      assert(size[0]>=0.0);
      assert(size[1]>=0.0);
      assert(size[2]>=0.0);
      

      int ncell = int(elems*GFactor);	// Calcolo numero di voxel
      if(ncell<mincells)
	ncell = mincells;
      
      dim[0] = 1;
      dim[1] = 1;
      dim[2] = 1;
      
      if(size[0]>eps)
	{
	  if(size[1]>eps)
	    {
	      if(size[2]>eps)
		{
		  double k = pow((double)(ncell/(size[0]*size[1]*size[2])),double(1.0/3.f));
		  dim[0] = int(size[0] * k);
		  dim[1] = int(size[1] * k);
		  dim[2] = int(size[2] * k);
		} 
	      else 
		{
		  dim[0] = int(::sqrt(ncell*size[0]/size[1]));
		  dim[1] = int(::sqrt(ncell*size[1]/size[0]));
		}
	    }
	  else
	    {
	      if(size[2]>eps)
		{
		  dim[0] = int(::sqrt(ncell*size[0]/size[2]));
		  dim[2] = int(::sqrt(ncell*size[2]/size[0]));
		}
	      else
		dim[0] = int(ncell);
	    }
	}
      else
	{
	  if(size[1]>eps)
	    {
	      if(size[2]>eps)
		{
		  dim[1] = int(::sqrt(ncell*size[1]/size[2]));
		  dim[2] = int(::sqrt(ncell*size[2]/size[1]));
		}
	      else
		dim[1] = int(ncell);
	    }
	  else if(size[2]>eps)
	    dim[2] = int(ncell);
	}
      dim[0] = math::Max(dim[0],1);
      dim[1] = math::Max(dim[1],1);
      dim[2] = math::Max(dim[2],1);
    }
 }
#endif