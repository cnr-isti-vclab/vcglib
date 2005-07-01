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
Revision 1.1  2005/03/15 11:43:18  cignoni
Removed BestDim function from the grid_static_ptr class and moved to a indipendent file (grid_util.h) for sake of generality.


****************************************************************************/
#ifndef __VCGLIB_GRID_UTIL
#define __VCGLIB_GRID_UTIL


 namespace vcg {
   
// classe di base per tutte le strutture dati di indexing spaziale basate su griglia.
// contiene tutte le funzioni solite per le conversioni tra point3f e point3i

   
   template <class ScalarType>
   class BasicGrid {
   public:
     Box3<ScalarType> bbox;
     /// Dimensione spaziale (lunghezza lati) del bbox
    Point3<ScalarType> dim;
		/// Dimensioni griglia in celle
    Point3i siz;
		/// Dimensioni di una cella
    Point3<ScalarType> voxel;

    

   // Dato un punto ritorna le coordinate della cella
    inline Point3i GridP( const Point3<ScalarType> & p ) const 
		{
      Point3i pi; return PToIP(p,pi);
		}

		/// Dato un punto 3d ritorna l'indice del box corrispondente
    inline void PToIP(const Point3<ScalarType> & p, Point3i &pi ) const
		{
			Point3<ScalarType> t = p - bbox.min;
			pi[0] = int( t[0]/voxel[0] );
			pi[1] = int( t[1]/voxel[1] );
			pi[2] = int( t[2]/voxel[2] );
		}

			/// Dato un box reale ritorna gli indici dei voxel compresi dentro un ibox
    void BoxToIBox( const Box3d & b, Box3i & ib ) const
		{
			Point3d t;

			t = (b.min - bbox.min);	// Traslo il box b;
			t[0] /= voxel[0];
			t[1] /= voxel[1];
			t[2] /= voxel[2];
			ib.min[0] = int( t[0] );	// Trasformazione in intero
			ib.min[1] = int( t[1] );
			ib.min[2] = int( t[2] );
			assert(ib.min[0]>=0 && ib.min[1]>=0 && ib.min[2]>=0);
			t = (b.max - bbox.min);	// Taslo il box b;
			t[0] /= voxel[0];
			t[1] /= voxel[1];
			t[2] /= voxel[2];
			ib.max[0] = int( ceil(t[0]) );	// Trasformazione in intero
			ib.max[1] = int( ceil(t[1]) );
			ib.max[2] = int( ceil(t[2]) );
			assert(ib.max[0]>=0 && ib.max[1]>=0 && ib.max[2]>=0);			
		}



   };

   /** Calcolo dimensioni griglia.
      Calcola la dimensione della griglia in funzione
      della ratio del bounding box e del numero di elementi
  */
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