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
Revision 1.1  2004/03/08 09:21:31  cignoni
Initial commit

****************************************************************************/

#ifndef __VCGLIB_UGRID
#define __VCGLIB_UGRID

#include <vector>
#include <algorithm>
#include <vcg/space/box3.h>
#include <vcg/space/line3.h>

namespace vcg {

	/** Static Uniform Grid
	    A spatial search structure for a accessing a container of objects. It is based on a uniform grid overlayed over a protion of space. 
      The grid partion the space into cells. Cells contains just pointers to the object that are stored elsewhere. 
      The set of object is meant to be static and pointer stable. 
      Useful for situation were many space related query are issued over the same dataset (ray tracing, measuring distances between meshes, re-detailing ecc.). Works well for distribution that ar reasonably uniform.
      How to use it:

     */
template < typename ContainerType >
class GridStaticPtr
{
public:

	/** Internal class for keeping the first pointer of object.
	    Definizione Link dentro la griglia. Classe di supporto per GridStaticObj.
	 */
		class Link
		{
		public:
				/// Costruttore di default
				inline Link(){};
				/// Costruttore con inizializzatori
				inline Link(  typename ContainerType::iterator const nt, const int ni ){
							assert(ni>=0);
							t = nt;
							i = ni;
				};

				
				inline bool operator <  ( const Link & l ) const{ 	return i <   l.i; } 
				inline bool operator <= ( const Link & l ) const{ 	return i <=  l.i; }
				inline bool operator >  ( const Link & l ) const{ 	return i >   l.i; }
				inline bool operator >= ( const Link & l ) const{ 	return i >=  l.i; }
				inline bool operator == ( const Link & l ) const{ 	return i ==  l.i; }
				inline bool operator != ( const Link & l ) const{ 	return i !=  l.i; }
			
				inline typename ContainerType::iterator & Elem() {
					return t;
				}
				inline int & Index() {
					return i;
				}

		private:
				/// Puntatore all'elemento T
				typename ContainerType::iterator t;
				/// Indirizzo del voxel dentro la griglia
				int i;
				

		};//end class Link

		typedef typename ContainerType::value_type ObjType;
		typedef ObjType* ObjPtr;
		typedef typename ObjType::ScalarType ScalarType;
		typedef Point3<ScalarType> Point3x;
		typedef Box3<ScalarType> Box3x;
		typedef Line3<ScalarType> Line3x;
		typedef Link* Cell;
		
    Box3x   bbox;
		Point3x dim;   /// Dimensione spaziale (lunghezza lati) del bbox
    Point3i siz;   /// Dimensioni griglia in celle
		Point3x voxel; /// Dimensioni di una cella
    
		/// Insieme di tutti i links
    std::vector<Link>   links;
		/// Griglia vera e propria
    std::vector<Cell> grid;

		/// Dato un punto, ritorna la cella che lo contiene
    inline Cell* Grid( const Point3d & p )
		{
			int x = int( (p[0]-bbox.min[0])/voxel[0] );
			int y = int( (p[1]-bbox.min[1])/voxel[1] );
			int z = int( (p[2]-bbox.min[2])/voxel[2] );

#ifndef NDEBUG
			if ( x<0 || x>=siz[0] || y<0 || y>=siz[1] || z<0 || z>=siz[2] )
				return NULL;
			else
#endif

			return grid.begin() + ( x+siz[0]*(y+siz[1]*z) );
		}
		/// Date le coordinate ritorna la cella
    inline Cell* Grid( const int x, const int y, const int z )
		{
#ifndef NDEBUG
			if ( x<0 || x>=siz[0] || y<0 || y>=siz[1] || z<0 || z>=siz[2] )
				assert(0);
				//return NULL;
			else
#endif
			assert(((unsigned int)x+siz[0]*y+siz[1]*z)<grid.size());
			return &*grid.begin() + ( x+siz[0]*(y+siz[1]*z) );
		}


		/// Date le coordinate di un grid point ritorna le celle che condividono
		/// l'edge cell che parte dal grid point in direzione axis
    inline void Grid( Point3i p, const int axis,
											std::vector<Cell*> & cl,
											std::vector<Point3i> &o) 
		{
#ifndef NDEBUG
			if ( p[0]<0 || p[0]>siz[0] || p[1]<0 || p[1]>siz[1] || p[2]<0 || p[2]>siz[2] )
				assert(0);
				//return NULL;
			else
#endif
			assert(((unsigned int) p[0]+siz[0]*p[1]+siz[1]*p[2])<grid.size());
			
			int axis0 = (axis+1)%3;
			int axis1 = (axis+2)%3;
			int i,j,x,y;
			x = p[axis0];
			y = p[axis1];
			for(i = max(x-1,0); i <= min( x,siz[axis0]-1);++i)	
				for(j = max(y-1,0); j <= min( y,siz[axis1]-1);++j){
					p[axis0]=i;
					p[axis1]=j;
					cl.push_back(Grid(p[0]+siz[0]*(p[1]+siz[1]*p[2]))); ;
					o.push_back(p);
					}
		}

	 Cell* Grid(const  int i) {
		return &grid[i];
		}
	void Grid( const Point3d & p, Cell & first, Cell & last )
	{
		Cell* g = Grid(s);
		
		first = *g;
		last  = *(g+1);
	}
	void Grid( const Cell* g, Cell & first, Cell & last )
	{
		first = *g;
		last  = *(g+1);
	}
	void Grid( const int x, const int y, const int z, Cell & first, Cell & last )
	{
		Cell* g = Grid(x,y,z);
		
		first = *g;
		last  = *(g+1);
	}

		/// Setta il bounding box della griglia
    void SetBBox( const Box3x & b )
		{
			bbox = b;
			dim  = b.max - b.min;
		}

    void SetSafeBBox( const Box3x & b )
		{
			Box3x btmp=b;
			btmp.InflateFix(0.01);
			bbox = btmp;
			dim  = bbox.max - bbox.min;
		}

		/// Dato un punto 3d ritorna l'indice del box corrispondente
    inline void PToIP(const Point3x & p, Point3i &pi ) const
		{
			Point3x t = p - bbox.min;
			pi[0] = int( t[0]/voxel[0] );
			pi[1] = int( t[1]/voxel[1] );
			pi[2] = int( t[2]/voxel[2] );
		}
		/// Dato un box reale ritorna gli indici dei voxel compresi dentro un ibox
    void BoxToIBox( const Box3x & b, Box3i & ib ) const
		{
			PToIP(b.min,ib.min);
			PToIP(b.max,ib.max);
		}

	
	void ShowStats(FILE *fp)
	{
				// Conto le entry
		//int nentry = 0;
		//Hist H;
		//H.SetRange(0,1000,1000);
		//int pg;
		//for(pg=0;pg<grid.size()-1;++pg)
		//	if( grid[pg]!=grid[pg+1] )
		//	{
		//		++nentry;
		//		H.Add(grid[pg+1]-grid[pg]);
		//	}

		//	fprintf(fp,"Uniform Grid: %d x %d x %d (%d voxels), %.1f%% full, %d links \nNon empty Cell Occupancy Distribution Avg: %f (%4.0f %4.0f %4.0f) \n",
		//	siz[0],siz[1],siz[2],grid.size()-1,
		//	double(nentry)*100.0/(grid.size()-1),links.size(),H.Avg(),H.Percentile(.25),H.Percentile(.5),H.Percentile(.75)
  //    
		//);
	}


	/** Returns the closest posistion of a point p and its distance
		@param p a 3d point
		@return The closest element
	*/
	ObjPtr  GetClosest( const Point3x & p, ScalarType & min_dist, Point3x  & res)
	{
		ScalarType dx = ( (p[0]-bbox.min[0])/voxel[0] );
		ScalarType dy = ( (p[1]-bbox.min[1])/voxel[1] );
		ScalarType dz = ( (p[2]-bbox.min[2])/voxel[2] );
 
		int ix = int( dx );
		int iy = int( dy );
		int iz = int( dz );

		double voxel_min=voxel[0];
		if (voxel_min<voxel[1]) voxel_min=voxel[1];
		if (voxel_min<voxel[2]) voxel_min=voxel[2];

		ScalarType radius=(dx-ScalarType(ix));
		if (radius>0.5)  radius=(1.0-radius);	radius*=voxel[0];

		ScalarType 
		tmp=dy-ScalarType(iy); if (tmp>0.5) tmp=1.0-tmp; tmp*=voxel[1]; if (radius>tmp) radius=tmp;
		tmp=dz-ScalarType(iz); if (tmp>0.5) tmp=1.0-tmp; tmp*=voxel[2]; if (radius>tmp) radius=tmp;

		Point3x t_res;
		//ScalarType min_dist=1e10;
		ObjPtr winner=NULL;

		Link  *first, *last;
		Link *l;
		if ((ix>=0) && (iy>=0) && (iz>=0) && (ix<siz[0]) && (iy<siz[1]) && (iz<siz[2])) {

			Grid( ix, iy, iz, first, last );
			for(l=first;l!=last;++l)
			{
				if (!l->Elem()->IsD() && l->Elem()->Dist(p,min_dist,t_res)) {
					winner=&*(l->Elem());
					res=t_res;

				}
			};
		};

		//return winner;

		Point3i done_min=Point3i(ix,iy,iz), done_max=Point3i(ix,iy,iz);

		//printf(".");

		while (min_dist>radius) {
			//if (dy-ScalarType(iy))
			done_min[0]--; if (done_min[0]<0) done_min[0]=0;
			done_min[1]--; if (done_min[1]<0) done_min[1]=0;
			done_min[2]--; if (done_min[2]<0) done_min[2]=0;
			done_max[0]++; if (done_max[0]>=siz[0]-1) done_max[0]=siz[0]-1;
			done_max[1]++; if (done_max[1]>=siz[1]-1) done_max[1]=siz[1]-1;
			done_max[2]++; if (done_max[2]>=siz[2]-1) done_max[2]=siz[2]-1;
			radius+=voxel_min;
			//printf("+");
			for (ix=done_min[0]; ix<=done_max[0]; ix++) 
			for (iy=done_min[1]; iy<=done_max[1]; iy++) 
			for (iz=done_min[2]; iz<=done_max[2]; iz++) 
			{
				Grid( ix, iy, iz, first, last );
				for(l=first;l!=last;++l)
				{
					if (!l->Elem()->IsD() && l->Elem()->Dist(p,min_dist,t_res)) {
						winner=&*(l->Elem());
						res=t_res;
					};
				};
			}
		};
		return winner;
	};

	/// Inserisce una mesh nella griglia. Nota: prima bisogna 
	/// chiamare SetBBox che setta dim in maniera corretta
	void Set( ContainerType & s )
	{
		Set(s,s.size());
	}


	/// Inserisce una mesh nella griglia. Nota: prima bisogna 
	/// chiamare SetBBox che setta dim in maniera corretta
	void Set( ContainerType & s,int _size )
	{
		Point3i _siz;

		BestDim( _size, dim, _siz );
		Set(s,_siz);
	}
	void Set(ContainerType & s, Point3i _siz)
	{	
		siz=_siz;
			// Calcola la dimensione della griglia
		voxel[0] = dim[0]/siz[0];
		voxel[1] = dim[1]/siz[1];
		voxel[2] = dim[2]/siz[2];

		// "Alloca" la griglia: +1 per la sentinella
		grid.resize( siz[0]*siz[1]*siz[2]+1 );

    		// Ciclo inserimento dei tetraedri: creazione link
		links.clear();
		typename ContainerType::iterator pt;
		for(pt=s.begin(); pt!=s.end(); ++pt)
		{
			Box3x bb;			// Boundig box del tetraedro corrente
				(*pt).GetBBox(bb);
				bb.Intersect(bbox);
				if(! bb.IsNull() )
				{
					Box3i ib;		// Boundig box in voxels
					BoxToIBox( bb,ib );

					int x,y,z;
					for(z=ib.min[2];z<=ib.max[2];++z)
					{
						 int bz = z*siz[1];
						 for(y=ib.min[1];y<=ib.max[1];++y)
						 {
							 int by = (y+bz)*siz[0];
							 for(x=ib.min[0];x<=ib.max[0];++x)
								// Inserire calcolo cella corrente
								// if( pt->Intersect( ... )
								links.push_back( Link(pt,by+x) );
						 }
					}
				}
		}
		// Push della sentinella
		links.push_back( Link(NULL,grid.size()-1) );

		// Ordinamento dei links
		sort( links.begin(), links.end() );

		// Creazione puntatori ai links
		typename std::vector<Link>::iterator pl;
		unsigned int pg;
		pl = links.begin();
		for(pg=0;pg<grid.size();++pg)
		{
			assert(pl!=links.end());

			grid[pg] = &*pl;
			while( pg == pl->Index() )	// Trovato inizio
			{
				++pl;		// Ricerca prossimo blocco
				if(pl==links.end())
				break;
			}
		}

	}

		/** Calcolo dimensioni griglia.
		 Calcola la dimensione della griglia in funzione
			 della ratio del bounding box e del numero di elementi
	 */
		static void BestDim( const int elems, const Point3x & size, Point3i & dim )
		{
				const int mincells   = 1;		// Numero minimo di celle
				const double GFactor = 1.0;		// GridEntry = NumElem*GFactor
				double diag = size.Norm();		// Diagonale del box
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
			dim[0] = math::Max(dim[0],0);
			dim[1] = math::Max(dim[1],1);
			dim[2] = math::Max(dim[2],2);
		}


	int MemUsed()
	{
		return sizeof(GridStaticObj)+ sizeof(Link)*links.size() + sizeof(Cell) * grid.size();
	}
}; //end class GridStaticObj

}; // end namespace

#endif
