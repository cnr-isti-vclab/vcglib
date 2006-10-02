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
Revision 1.20  2006/07/26 08:12:56  pietroni
added InitEmpty Function

Revision 1.19  2006/07/10 12:43:13  turini
explicit cast in _IsInHtable() to resolve a warning

Revision 1.18  2006/04/20 08:30:27  cignoni
small GCC compiling issues

Revision 1.17  2006/01/23 21:26:57  ponchio
gcc compatibility (templates mostly)
bbox -> this->bbox
More consistent use of Box3x and such.

Revision 1.16  2006/01/23 15:26:31  ponchio
P1 --> HASH_P1
Old definition was conflicting with functions in segment.h

Revision 1.15  2005/11/23 15:55:35  pietroni
1 warning corrected

Revision 1.14  2005/11/07 14:15:36  pietroni
added dynamic spatial hashing class for dynamic updating of entries (and relative functions)

Revision 1.13  2005/10/05 17:04:18  pietroni
corrected bug on Set Function .... bbox must be exetended in order to have'nt any object on his borde

Revision 1.12  2005/10/03 13:58:21  pietroni
added GetInSphere and GetInBox functions

Revision 1.11  2005/10/03 10:05:26  pietroni
changed Set functions, added possibility to pass the bbox as parameter

Revision 1.10  2005/09/30 13:14:59  pietroni
added wrapping to functions defined in GridClosest:
- GetClosest
- GetKClosest
- DoRay

Revision 1.9  2005/09/21 14:22:49  pietroni
Added DynamicSpatialHAshTable class

Revision 1.8  2005/09/19 13:35:45  pietroni
use of standard grid interface
use of vector instead of map inside the cell
removed closest iterator

Revision 1.7  2005/06/15 11:44:47  pietroni
minor changes

Revision 1.6  2005/06/01 13:47:59  pietroni
resolved hash code conflicts

Revision 1.5  2005/03/15 09:50:44  ganovelli
there was a debug line, now removed

Revision 1.4  2005/03/14 15:11:18  ganovelli
ClosestK added and other minor changes

Revision 1.3  2005/03/11 15:25:29  ganovelli
added ClosersIterator and other minor changes. Not compatible with the previous version.
Still other modifications to do (temporary commit)

Revision 1.2  2005/02/21 12:13:25  ganovelli
added vcg header



****************************************************************************/

#ifndef VCGLIB_SPATIAL_HASHING
#define VCGLIB_SPATIAL_HASHING


#define HASH_P0 73856093
#define HASH_P1 19349663
#define HASH_P2 83492791

#include <vcg/space/index/grid_util.h>
#include <vcg/space/index/grid_closest.h>
//#include <map>
#include <vector>
#include <algorithm>
#ifdef WIN32
 #ifndef __MINGW32__
  #include <hash_map>
  #define STDEXT stdext
 #else
  #include <ext/hash_map>
  #define STDEXT __gnu_cxx
 #endif
#else
 #include <ext/hash_map>
 #define STDEXT __gnu_cxx
#endif


namespace vcg{


	/** Spatial Hash Table
	Spatial Hashing as described in
	"Optimized Spatial Hashing for Coll	ision Detection of Deformable Objects",
	Matthias Teschner and Bruno Heidelberger and Matthias Muller and Danat Pomeranets and Markus Gross
	*/
	template < typename OBJTYPE,class FLT=double>
	class SpatialHashTable:public BasicGrid<OBJTYPE,FLT>
	{

	public:

		typedef OBJTYPE ObjType;
		typedef ObjType* ObjPtr;
		typedef typename ObjType::ScalarType ScalarType;
		typedef Point3<ScalarType> CoordType;
        typedef typename BasicGrid<OBJTYPE, FLT>::Box3x Box3x;

//		typedef typename SpatialHashTable<ObjType,FLT> SpatialHashType;
        typedef SpatialHashTable SpatialHashType;
		//typedef typename SpatialHashTable<ObjType,FLT> GridType;

		//type of container of pointer to object in a Cell
		//typedef typename std::pair<ObjType*,int> EntryType ;
		class  EntryType : public std::pair<ObjType*,int>
		{
		public:
			EntryType(ObjType* sim,const int  &_tempMark)
			{
				this->first=sim;
				this->second=_tempMark;
			}

			ObjType& operator *(){return (*this->first);}
		};

		typedef typename std::vector<EntryType> CellContainerType;
		typedef typename CellContainerType::iterator IteMap;
		typedef EntryType* CellIterator;

		//This Class Identify the cell
		struct Cell
		{
		protected:

			Point3i cell_n;//cell number

		public:

			//elements
			CellContainerType _entries;

			Cell()
			{}

			Cell(ObjType* sim,Point3i _cell,const int  &_tempMark)
			{
				_entries.push_back(EntryType(sim,_tempMark));
				cell_n=_cell;
			}


			///return the number of elements stored in the cell
			int Size()
			{return (int)(_entries.size());}

			///find the simplex into the cell
			bool Find(ObjType* sim,IteMap &I)
			{
				for (I=_entries.begin();I<_entries.end();I++)
					if ((*I).first==sim)
						return true;
				return false;
			}

			///update or insert an element into a cell
			void Update(ObjType* sim, const int & _tempMark)
			{
				IteMap I;
				if (Find(sim,I))
					(*I).second=_tempMark;
				else
					_entries.push_back(EntryType(sim,_tempMark));
			}

			Point3i CellN()
			{return cell_n;}

			bool operator ==(const Cell &h)
			{return (cell_n==h.CellN());}

			bool operator !=(const Cell &h)
			{return ((cell_n!=h.CellN()));}

		}; // end struct Cell

		//hash table definition
		typedef typename STDEXT::hash_multimap<int,Cell> Htable;
		//record of the hash table
		typedef typename std::pair<int,Cell> HRecord;
		//iterator to the hash table
		typedef typename Htable::iterator IteHtable;

		SpatialHashTable(){HashSpace=1000;};//default value for hash_space
		~SpatialHashTable(){};

		int tempMark;

	protected:

		Htable hash_table;

		///number of possible hash code [0...HashSpace]
		int HashSpace;

		///number of conflicts created
		int conflicts;

		///insert a new cell
		void _InsertNewHentry(ObjType* s,Point3i cell)
		{
			int h=Hash(cell);
			hash_table.insert(HRecord(h,Cell(s,cell,tempMark)));
		}

		///return true and return the iterator to the cell if exist
		bool _IsInHtable(Point3i cell,IteHtable &result)
		{
			int h=Hash(cell);
			int count = (int) hash_table.count(h);
			if (count==0)///in this case there is no entry for that key
				return false;
			else
			{
				std::pair<IteHtable, IteHtable> p =hash_table.equal_range(h);
				IteHtable i = p.first;

				while((i != p.second)&&((*i).second.CellN()!=cell))++i;

				if (i==p.second)///the scan is terminated and we have not found the right cell
				{
					conflicts++;
					return false;
				}
				else	///we have found the right cell
				{
					result=i;
					return true;
				}
			}
		}

		virtual void _UpdateHMark(ObjType* s){(void)s;}

		///insert an element in a specified cell if the cell doesn't exist than
		///create it.
		void _InsertInCell(ObjType* s,Point3i cell)
		{
			IteHtable I;
			if (!_IsInHtable(cell,I))
				_InsertNewHentry(s,cell);
			else///there is the entry specified by the iterator I so update only the temporary mark
				(*I).second.Update(s,tempMark);
		}

		// hashing
		const int Hash(Point3i p) const
		{
			return ((p.V(0)*HASH_P0 ^ p.V(1)*HASH_P1 ^ p.V(2)*HASH_P2)%HashSpace);
		}


	public:

		/////We need some extra space for numerical precision.
		//template <class Box3Type>
		// void SetBBox( const Box3Type & b )
		//{
		//	bbox.Import( b );
		//	ScalarType t = bbox.Diag()/100.0;
		//	if(t == 0) t = ScalarType(1e20);  // <--- Some doubts on this (Cigno 5/1/04)
		//	bbox.Offset(t);
		//	dim  = bbox.max - bbox.min;
		//}

		vcg::Box3i Add( ObjType* s)
		{
			Box3<ScalarType> b;
			s->GetBBox(b);
			vcg::Box3i bb;
			BoxToIBox(b,bb);
			//then insert all the cell of bb
			for (int i=bb.min.X();i<=bb.max.X();i++)
				for (int j=bb.min.Y();j<=bb.max.Y();j++)
					for (int k=bb.min.Z();k<=bb.max.Z();k++)
						_InsertInCell(s,vcg::Point3i(i,j,k));

			_UpdateHMark(s);
			return bb;
		}

		
		/// set an empty spatial hash table
		void InitEmpty(const Box3x &_bbox, vcg::Point3i grid_size)
		{
			Box3x b;
			Box3x &bbox = this->bbox;
			CoordType &dim = this->dim;
			Point3i &siz = this->siz;
			CoordType &voxel = this->voxel;

			assert(!_bbox.IsNull());
			bbox=_bbox;
			dim  = bbox.max - bbox.min;
			assert((grid_size.V(0)>0)&&(grid_size.V(1)>0)&&(grid_size.V(2)>0));
			siz=grid_size;

			voxel[0] = dim[0]/siz[0];
			voxel[1] = dim[1]/siz[1];
			voxel[2] = dim[2]/siz[2];
		}

		/// Insert a mesh in the grid.
		template <class OBJITER>
			void Set(const OBJITER & _oBegin, const OBJITER & _oEnd,const Box3x &_bbox=Box3x() )
		{

			OBJITER i;
			Box3x b;
            Box3x &bbox = this->bbox;
            CoordType &dim = this->dim;
            Point3i &siz = this->siz;
            CoordType &voxel = this->voxel;

			int _size=std::distance<OBJITER>(_oBegin,_oEnd);
			if(!_bbox.IsNull()) this->bbox=_bbox;
			else
			{
				for(i = _oBegin; i!= _oEnd; ++i)
				{
					(*i).GetBBox(b);
					this->bbox.Add(b);
				}
				///inflate the bb calculated
				ScalarType infl=bbox.Diag()/_size;
				bbox.min -= CoordType(infl,infl,infl);
				bbox.max += CoordType(infl,infl,infl);
			}

				dim  = bbox.max - bbox.min;
				BestDim( _size, dim, siz );
				// find voxel size
				voxel[0] = dim[0]/siz[0];
				voxel[1] = dim[1]/siz[1];
				voxel[2] = dim[2]/siz[2];

				for(i = _oBegin; i!= _oEnd; ++i)
					Add(&(*i));
		}


		///return the simplexes of the cell that contain p
		void Grid( const Point3d & p, CellIterator & first, CellIterator & last )
		{
			IteHtable I;
			vcg::Point3i _c;
			this->PToIP(p,_c);
			Grid(_c,first,last);
		}

		///return the simplexes on a specified cell
		void Grid( int x,int y,int z, CellIterator & first, CellIterator & last )
		{
			Grid(vcg::Point3i(x,y,z),first,last);
		}

		///return the simplexes on a specified cell
		void Grid( const Point3i & _c, CellIterator & first, CellIterator & last )
		{
			IteHtable I;
			if (_IsInHtable(_c,I))//if there is the cell then
			{	///return pointers to first and last element cell elems
				first= &*(*I).second._entries.begin();
				last=  &*(*I).second._entries.end();
			}
			else
			{	///return 2 equals pointers
				first=&*(*hash_table.begin()).second._entries.begin();
				last= &*(*hash_table.begin()).second._entries.begin();
			}
		}

		///return the number of elemnts in the cell and the iterator to the cell
		///if the cell exist
		int numElemCell(Point3i _c,IteHtable &I)
		{
			if (_IsInHtable(_c,I))
				return ((*I).second.Size());
			else
				return 0;
		}

		///return the number of cell created
		int CellNumber()
		{return (hash_table.size());}

		int Conflicts()
		{return conflicts;}

		void Clear()
		{hash_table.clear();}

		void SetHashKeySpace(int n)
		{HashSpace=n;}

		void UpdateTmark()
		{tempMark++;}


		template <class OBJPOINTDISTFUNCTOR, class OBJMARKER>
			ObjPtr  GetClosest(OBJPOINTDISTFUNCTOR & _getPointDistance, OBJMARKER & _marker,
			const CoordType & _p, const ScalarType & _maxDist,ScalarType & _minDist, CoordType & _closestPt)
		{
			return (vcg::GridClosest<SpatialHashType,OBJPOINTDISTFUNCTOR,OBJMARKER>(*this,_getPointDistance,_marker, _p,_maxDist,_minDist,_closestPt));
		}


		template <class OBJPOINTDISTFUNCTOR, class OBJMARKER, class OBJPTRCONTAINER,class DISTCONTAINER, class POINTCONTAINER>
			unsigned int GetKClosest(OBJPOINTDISTFUNCTOR & _getPointDistance,OBJMARKER & _marker,
			const unsigned int _k, const CoordType & _p, const ScalarType & _maxDist,OBJPTRCONTAINER & _objectPtrs,
			DISTCONTAINER & _distances, POINTCONTAINER & _points)
		{
			return (vcg::GridGetKClosest<SpatialHashType,
				OBJPOINTDISTFUNCTOR,OBJMARKER,OBJPTRCONTAINER,DISTCONTAINER,POINTCONTAINER>
				(*this,_getPointDistance,_marker,_k,_p,_maxDist,_objectPtrs,_distances,_points));
		}

		template <class OBJPOINTDISTFUNCTOR, class OBJMARKER, class OBJPTRCONTAINER, class DISTCONTAINER, class POINTCONTAINER>
		unsigned int GetInSphere(OBJPOINTDISTFUNCTOR & _getPointDistance,
		OBJMARKER & _marker,
		const CoordType & _p,
		const ScalarType & _r,
		OBJPTRCONTAINER & _objectPtrs,
		DISTCONTAINER & _distances,
		POINTCONTAINER & _points)
		{
			return(vcg::GridGetInSphere<SpatialHashType,
				OBJPOINTDISTFUNCTOR,OBJMARKER,OBJPTRCONTAINER,DISTCONTAINER,POINTCONTAINER>
				(*this,_getPointDistance,_marker,_p,_r,_objectPtrs,_distances,_points));
		}

		template <class OBJMARKER, class OBJPTRCONTAINER>
			unsigned int GetInBox(OBJMARKER & _marker,
            const Box3x _bbox,
			OBJPTRCONTAINER & _objectPtrs)
		{
			return(vcg::GridGetInBox<SpatialHashType,OBJMARKER,OBJPTRCONTAINER>
				  (*this,_marker,_bbox,_objectPtrs));
		}

		template <class OBJRAYISECTFUNCTOR, class OBJMARKER>
			ObjPtr DoRay(OBJRAYISECTFUNCTOR & _rayIntersector, OBJMARKER & _marker, const Ray3<ScalarType> & _ray, const ScalarType & _maxDist, ScalarType & _t)
		{
			return(vcg::GridDoRay<SpatialHashType,OBJRAYISECTFUNCTOR,OBJMARKER>
				  (*this,_rayIntersector,_marker,_ray,_maxDist,_t));
		}


	}; // end class

	/** Spatial Hash Table Dynamic
	Update the Hmark value on the simplex for dynamic updating of contents of the cell.
	The simplex must have the HMark() function.
	*/
	template < typename ContainerType,class FLT=double>
	class DynamicSpatialHashTable: public SpatialHashTable<ContainerType,FLT>
	{
	public:
        typedef typename SpatialHashTable<ContainerType,FLT>::CoordType CoordType;
        typedef typename SpatialHashTable<ContainerType,FLT>::ObjType ObjType;
        typedef typename SpatialHashTable<ContainerType,FLT>::ObjPtr ObjPtr;
        typedef typename SpatialHashTable<ContainerType,FLT>::Box3x Box3x;
        typedef typename SpatialHashTable<ContainerType,FLT>::CellIterator CellIterator;

		void _UpdateHMark(ObjType* s){ s->HMark() = this->tempMark;}

		void getInCellUpdated(vcg::Point3i cell,std::vector<ObjPtr> &elems)
		{
			CellIterator first,last,l;
			Grid(cell,first,last);
			for (l=first;l!=last;l++)
			{
				if ((l->second)>=(**l).HMark())
					elems.push_back(&(**l));
			}
		}

	};



}// end namespace

#undef HASH_P0
#undef HASH_P1
#undef HASH_P2

#endif
