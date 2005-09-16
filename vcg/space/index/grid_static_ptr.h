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
Revision 1.23  2005/09/15 13:16:42  spinelli
fixed bugs

Revision 1.22  2005/09/15 11:14:39  pietroni
minor changes

Revision 1.21  2005/09/14 13:27:38  spinelli
minor changes

Revision 1.20  2005/09/14 12:57:52  pietroni
canged template parameters for Closest Function (use of TempMark class)

Revision 1.19  2005/09/14 09:05:32  pietroni
added * operator to Link
modified getClosest in order to use Temporary mark
corrected bug on functor calling compilation

Revision 1.18  2005/09/09 11:29:21  m_di_benedetto
Modified old GetClosest() to respect old min_dist semantic (in/out) and removed #included <limits>

Revision 1.17  2005/09/09 11:11:15  m_di_benedetto
#included <limits> for std::numeric_limits<ScalarType>::max() and corrected parameters bug in old GetClosest();

Revision 1.16  2005/09/09 11:01:02  m_di_benedetto
Modified GetClosest(): now it uses a functor for distance calculation.
Added comments and a GetClosest() method with backward compatibility.

Revision 1.15  2005/08/26 09:27:58  cignoni
Added a templated version of SetBBox

Revision 1.14  2005/08/02 11:18:36  pietroni
exetended form BasicGrid, changed type of t in class Link (from Iterator to Pointer to the object)

Revision 1.13  2005/04/14 17:23:08  ponchio
*** empty log message ***

Revision 1.12  2005/03/15 11:43:18  cignoni
Removed BestDim function from the grid_static_ptr class and moved to a indipendent file (grid_util.h) for sake of generality.

Revision 1.11  2005/01/03 11:21:26  cignoni
Added some casts

Revision 1.10  2004/09/28 10:25:05  ponchio
SetBox minimal change.

Revision 1.9  2004/09/23 14:29:42  ponchio
Small bugs fixed.

Revision 1.8  2004/09/23 13:44:25  ponchio
Removed SetSafeBBox. SetBBox is now safe enough.

Revision 1.7  2004/09/09 12:44:39  fasano
included stdio.h

Revision 1.6  2004/09/09 08:39:29  ganovelli
minor changes for gcc

Revision 1.5  2004/06/25 18:34:23  ganovelli
added Grid to return all the cells sharing a specified edge

Revision 1.4  2004/06/23 15:49:03  ponchio
Added some help and inndentation

Revision 1.3  2004/05/12 18:50:58  ganovelli
changed calls to Dist

Revision 1.2  2004/05/11 14:33:46  ganovelli
changed to grid_static_obj to grid_static_ptr

Revision 1.1  2004/05/10 14:44:13  ganovelli
created

Revision 1.1  2004/03/08 09:21:31  cignoni
Initial commit

****************************************************************************/

#ifndef __VCGLIB_UGRID
#define __VCGLIB_UGRID

#include <vector>
#include <algorithm>
#include <stdio.h>

#include <vcg/space/box3.h>
#include <vcg/space/line3.h>
#include <vcg/space/index/grid_util.h>
#include <vcg/simplex/face/distance.h>
namespace vcg {

	/** Static Uniform Grid
	A spatial search structure for a accessing a container of objects. 
	It is based on a uniform grid overlayed over a protion of space. 
	The grid partion the space into cells. Cells contains just pointers 
	to the object that are stored elsewhere.
	The set of objects is meant to be static and pointer stable. 

	Useful for situation were many space related query are issued over 
	the same dataset (ray tracing, measuring distances between meshes, 
	re-detailing ecc.). 
	Works well for distribution that ar reasonably uniform.
	How to use it:
	ContainerType must have a 'value_type' typedef inside.
	(stl containers already have it)

	Objects pointed by cells (of kind 'value_type') must have
	a 'ScalarType' typedef (float or double usually)
	and a member function:

	void GetBBox(Box3<ScalarType> &b)
	which return the bounding box of the object

	When using the GetClosest() method, the user must supply a functor object
	(whose type is a method template argument) which expose the following
	operator ():

	bool operator () (const ObjType & obj, const Point3f & point, ScalarType & mindist, Point3f & result);
	which return true if the distance from point to the object 'obj' is < mindist
	and set mindist to said distance, and result must be set as the closest 
	point of the object to point)
	*/

	template < typename ContainerType,class FLT=float >
	class GridStaticPtr:public BasicGrid<FLT>
	{
	public:

		typedef typename ContainerType::value_type ObjType;
		typedef ObjType* ObjPtr;
		typedef typename ObjType::ScalarType ScalarType;
		typedef Point3<ScalarType> CoordType;
		typedef Box3<ScalarType> Box3x;
		typedef Line3<ScalarType> Line3x;


		/** Internal class for keeping the first pointer of object.
		Definizione Link dentro la griglia. Classe di supporto per GridStaticObj.
		*/
		class Link
		{
		public:
			/// Costruttore di default
			inline Link(){};
			/// Costruttore con inizializzatori
			inline Link(ObjPtr nt, const int ni ){
				assert(ni>=0);
				t = nt;
				i = ni;
			};


			inline bool operator <  ( const Link & l ) const{ return i <   l.i; } 
			inline bool operator <= ( const Link & l ) const{ return i <=  l.i; }
			inline bool operator >  ( const Link & l ) const{ return i >   l.i; }
			inline bool operator >= ( const Link & l ) const{ return i >=  l.i; }
			inline bool operator == ( const Link & l ) const{ return i ==  l.i; }
			inline bool operator != ( const Link & l ) const{ return i !=  l.i; }

			inline ObjPtr & Elem() {
				return t;
			}

			ObjType &operator *(){return *(t);}

			inline int & Index() {
				return i;
			}

		private:
			/// Puntatore all'elemento T
			ObjPtr t;
			/// Indirizzo del voxel dentro la griglia
			int i;


		};//end class Link

		typedef Link* Cell;
		typedef Cell CellIterator;

		std::vector<Link>   links;   /// Insieme di tutti i links

		std::vector<Cell> grid;   /// Griglia vera e propria


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


		/// Date le coordinate di un grid point (corner minx,miy,minz) ritorna le celle che condividono
		/// l'edge cell che parte dal grid point in direzione axis
		inline void Grid( Point3i p, const int axis,
			std::vector<Cell*> & cl) 
		{
#ifndef NDEBUG
			if ( p[0]<0 || p[0]>siz[0] || 
				p[1]<0 || p[1]>siz[1] || 
				p[2]<0 || p[2]>siz[2] )
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
					cl.push_back(Grid(p[0]+siz[0]*(p[1]+siz[1]*p[2])));
				}
		}

		Cell* Grid(const  int i) {
			return &grid[i];
		}
		void Grid( const Point3d & p, Cell & first, Cell & last )
		{
			Cell* g = Grid(p);

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

		/// Set the bounding box of the grid
		///We need some extra space for numerical precision.
		template <class Box3Type>
    void SetBBox( const Box3Type & b )
		{
			bbox.Import( b );
			ScalarType t = bbox.Diag()/100.0;
			if(t == 0) t = ScalarType(1e20);  // <--- Some doubts on this (Cigno 5/1/04)
			bbox.Offset(t);
			dim  = bbox.max - bbox.min;
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
		@param max_dist maximum distance not to search beyond.
		@param dist_funct (templated type) a functor object used to calculate distance from a grid object to the point p.
		@param min_dist the returned closest distance
		@param res the returned closest point
		@return The closest element
		*/
		/*
			A DISTFUNCT object must implement an operator () with signature:
				bool operator () (const ObjType& obj, const CoordType & p, ScalarType & min_dist, CoordType & res);
		*/
		template <class DISTFUNCTOR,class TMARKER>
			ObjPtr  GetClosest( const CoordType & p, const ScalarType & max_dist, DISTFUNCTOR & dist_funct, ScalarType & min_dist, CoordType & res,TMARKER tm)
		{
			// Initialize min_dist with max_dist to exploit early rejection test.
			min_dist = max_dist;

			ScalarType dx = ( (p[0]-bbox.min[0])/voxel[0] );
			ScalarType dy = ( (p[1]-bbox.min[1])/voxel[1] );
			ScalarType dz = ( (p[2]-bbox.min[2])/voxel[2] );

			int ix = int( dx );
			int iy = int( dy );
			int iz = int( dz );

			if (!bbox.IsIn(p))
				assert (0);///the grid has to be extended until the point

			double voxel_min=voxel[0];
			if (voxel_min<voxel[1]) voxel_min=voxel[1];
			if (voxel_min<voxel[2]) voxel_min=voxel[2];

			ScalarType radius=(dx-ScalarType(ix));
			if (radius>0.5)  radius=(1.0-radius);	radius*=voxel[0];

			ScalarType tmp=dy-ScalarType(iy); 
			if (tmp>0.5) tmp=1.0-tmp; 
			tmp*=voxel[1]; 
			if (radius>tmp) radius=tmp;
			tmp=dz-ScalarType(iz); 
			if (tmp>0.5) tmp=1.0-tmp; 
			tmp*=voxel[2]; 
			if (radius>tmp) radius=tmp;

			CoordType t_res;
			//ScalarType min_dist=1e10;
			ObjPtr winner=NULL;

			tm.UnMarkAll();

			Link  *first, *last;
			Link *l;
			if ((ix>=0) && (iy>=0) && (iz>=0) && 
				(ix<siz[0]) && (iy<siz[1]) && (iz<siz[2])) {

					Grid( ix, iy, iz, first, last );
					for(l=first;l!=last;++l)
						if (!(**l).IsD())
						{
							if( ! tm.IsMarked(l->Elem()))
							{
								//if (!l->Elem()->IsD() && l->Elem()->Dist(p,min_dist,t_res)) {
								//if (!l->Elem()->IsD() && dist_funct(*(l->Elem()), p, min_dist, t_res)) { // <-- NEW: use of distance functor
								if (dist_funct((**l), p, min_dist, t_res))  // <-- NEW: use of distance functor
								{
									winner=l->Elem();
									res=t_res;
								}
								tm.Mark(l->Elem());
							}
						}
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
								if (!(**l).IsD())
								{
									if( ! tm.IsMarked(l->Elem()))
									{
										//if (!l->Elem()->IsD() && l->Elem()->Dist(p,min_dist,t_res)) {
										if (dist_funct((**l), p, min_dist, t_res)) // <-- NEW: use of distance functor
										{
											winner=l->Elem();
											res=t_res;
										};
										tm.Mark(l->Elem());
									}
								}
							};
						}
			};
			return winner;
		};


		/** Returns the closest posistion of a point p and its distance (OLD VERSION)
		@param p a 3d point
		@param min_dist the returned closest distance
		@param res the returned closest point
		@return The closest element
		*/
		/*
			NOTE: kept for backward compatibility.
			Same as template <class DISTFUNCT> GetClosest() but without maximum distance rejection
			and distance functor and with the assumption that ObjType expose a Dist() method
			acting like a DISTFUNCT funcor;
		*/

		class BackCompDist {
				public:
					inline bool operator () (const ObjType & obj, const CoordType & pt, ScalarType & mindist, CoordType & result) {
						return (vcg::face::PointDistance<ObjType>(obj,pt, mindist, result));
					}
			};

		template <class TMARKER>
		ObjPtr  GetClosest( const CoordType & p, ScalarType & min_dist, CoordType & res,TMARKER &tm) {
			
			const ScalarType max_dist = min_dist;
			return (this->GetClosest<BackCompDist,TMARKER>(p, max_dist, BackCompDist(), min_dist, res,tm));
		}

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
								links.push_back( Link(&(*pt),by+x) );
						}
					}
				}
			}
			// Push della sentinella
			/*links.push_back( Link((typename ContainerType::iterator)NULL,
			(grid.size()-1)));*/

			links.push_back( Link(NULL,
				(grid.size()-1)));

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
				while( (int)pg == pl->Index() )	// Trovato inizio
				{
					++pl;		// Ricerca prossimo blocco
					if(pl==links.end())
						break;
				}
			}

		}


		int MemUsed()
		{
			return sizeof(GridStaticPtr)+ sizeof(Link)*links.size() + 
				sizeof(Cell) * grid.size();
		}
	}; //end class GridStaticObj

}; // end namespace

#endif
