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
Revision 1.2  2005/09/28 08:27:11  cignoni
Added a control to avoid multiple check of the same cells during radial expansion
Still miss some code to properly initialize when point is out of the BBox of the grid.

Revision 1.1  2005/09/27 15:09:38  cignoni
First Version


****************************************************************************/

///** Returns the closest posistion of a point _p and its distance
//@param _p a 3d point
//@param _maxDist maximum distance not to search beyond.
//@param _getPointDistance (templated type) a functor object used to calculate distance from a grid object to the point _p.
//@param _minDist the returned closest distance
//@param _closestPt the returned closest point
//@return The closest element
//*/
///*
//	A DISTFUNCT object must implement an operator () with signature:
//		bool operator () (const ObjType& obj, const CoordType & _p, ScalarType & min_dist, CoordType & _closestPt);
//*/
#ifndef __VCGLIB_GRID_CLOSEST
#define __VCGLIB_GRID_CLOSEST

#include <vcg/space/index/space_iterators.h>

namespace vcg{
	
	template <class SPATIAL_INDEX,class OBJPOINTDISTFUNCTOR, class OBJMARKER>
		typename SPATIAL_INDEX::ObjPtr  GridClosest(SPATIAL_INDEX &Si,OBJPOINTDISTFUNCTOR _getPointDistance,
		OBJMARKER & _marker,const typename SPATIAL_INDEX::CoordType & _p,
		const typename SPATIAL_INDEX::ScalarType & _maxDist,typename SPATIAL_INDEX::ScalarType & _minDist,
		typename SPATIAL_INDEX:: CoordType &_closestPt)
	{
		typedef SPATIAL_INDEX::ObjPtr ObjPtr;
		typedef SPATIAL_INDEX SpatialIndex;
		typedef SPATIAL_INDEX::CoordType CoordType;
		typedef SPATIAL_INDEX::ScalarType ScalarType;

		// Initialize min_dist with _maxDist to exploit early rejection test.
		_minDist = _maxDist;

		ScalarType dx = ( (_p[0]-Si.bbox.min[0])/Si.voxel[0] );
		ScalarType dy = ( (_p[1]-Si.bbox.min[1])/Si.voxel[1] );
		ScalarType dz = ( (_p[2]-Si.bbox.min[2])/Si.voxel[2] );

		int ix = int( dx );
		int iy = int( dy );
		int iz = int( dz );

		if (ix<0) ix=0;
		if (iy<0) iy=0;
		if (iz<0) iz=0;
		if (ix>=Si.siz[0]-1) ix=Si.siz[0]-1;
		if (iy>=Si.siz[1]-1) iy=Si.siz[1]-1;
		if (iz>=Si.siz[2]-1) iz=Si.siz[2]-1;

		if (!Si.bbox.IsIn(_p)){
			assert (0);///the grid has to be extended until the point

		}
		double voxel_min=Si.voxel[0];
		if (voxel_min<Si.voxel[1]) voxel_min=Si.voxel[1];
		if (voxel_min<Si.voxel[2]) voxel_min=Si.voxel[2];

		ScalarType radius=(dx-ScalarType(ix));
		if (radius>0.5)  radius=(1.0-radius);	radius*=Si.voxel[0];

		ScalarType tmp=dy-ScalarType(iy); 
		if (tmp>0.5) tmp=1.0-tmp; 
		tmp*=Si.voxel[1]; 
		if (radius>tmp) radius=tmp;
		tmp=dz-ScalarType(iz); 
		if (tmp>0.5) tmp=1.0-tmp; 
		tmp*=Si.voxel[2]; 
		if (radius>tmp) radius=tmp;

		CoordType t_res;
		//ScalarType min_dist=1e10;
		ObjPtr winner=NULL;

		_marker.UnMarkAll();

		SpatialIndex::CellIterator first,last;
		SpatialIndex::CellIterator l;

		if ((ix>=0) && (iy>=0) && (iz>=0) && 
			(ix<Si.siz[0]) && (iy<Si.siz[1]) && (iz<Si.siz[2])) {

				Si.Grid( ix, iy, iz, first, last );
				for(l=first;l!=last;++l)
					if (!(**l).IsD())
					{
						ObjPtr elem=&(**l);
						if(!_marker.IsMarked(elem))
						{
							//if (!l->Elem()->IsD() && l->Elem()->Dist(_p,min_dist,t_res)) {
							//if (!l->Elem()->IsD() && _getPointDistance(*(l->Elem()), _p, min_dist, t_res)) { // <-- NEW: use of distance functor
							if (_getPointDistance((**l), _p,_minDist, t_res))  // <-- NEW: use of distance functor
							{
								winner=elem;
								_closestPt=t_res;
							}
							_marker.Mark(elem);
						}
					}
			};

		//return winner;

		// the portion of the grid that have already been checked.
		Point3i done_min, done_max; 
		// the new box that we want to traverse: todo is a superset of done.
		Point3i todo_min=Point3i(ix,iy,iz), todo_max=Point3i(ix,iy,iz);


		// we should traverse only (todo - done).
		while (_minDist>radius) {
			done_min=todo_min; done_max=todo_max;
			todo_min[0]--; if (todo_min[0]<0) todo_min[0]=0;
			todo_min[1]--; if (todo_min[1]<0) todo_min[1]=0;
			todo_min[2]--; if (todo_min[2]<0) todo_min[2]=0;
			todo_max[0]++; if (todo_max[0]>=Si.siz[0]-1) todo_max[0]=Si.siz[0]-1;
			todo_max[1]++; if (todo_max[1]>=Si.siz[1]-1) todo_max[1]=Si.siz[1]-1;
			todo_max[2]++; if (todo_max[2]>=Si.siz[2]-1) todo_max[2]=Si.siz[2]-1;
			radius+=voxel_min;
			for (ix=todo_min[0]; ix<=todo_max[0]; ix++) 
				for (iy=todo_min[1]; iy<=todo_max[1]; iy++) 
					for (iz=todo_min[2]; iz<=todo_max[2]; iz++) 
						if(ix<done_min[0] || ix>done_max[0] ||  // this test is to avoid to re-process already analyzed cells.
							iy<done_min[1] || iy>done_max[1] ||
							iz<done_min[2] || iz>done_max[2] )
						{
							Si.Grid( ix, iy, iz, first, last );
							for(l=first;l!=last;++l)
							{
								if (!(**l).IsD())
								{
									ObjPtr elem=&(**l);
									if( ! _marker.IsMarked(elem))
									{
										//if (!l->Elem()->IsD() && l->Elem()->Dist(_p,min_dist,t_res)) {
										if (_getPointDistance((**l), _p, _minDist, t_res)) 
										{
											winner=elem;
											_closestPt=t_res;
										};
										_marker.Mark(elem);
									}
								}
							};
						}
		};
		return winner;
	};

	template <class SPATIALINDEXING,class OBJPOINTDISTFUNCTOR, class OBJMARKER, 
	class OBJPTRCONTAINER, class DISTCONTAINER, class POINTCONTAINER>
	unsigned int GridGetKClosest(SPATIALINDEXING &_Si,OBJPOINTDISTFUNCTOR & _getPointDistance,
	OBJMARKER & _marker, const unsigned int _k, const typename SPATIALINDEXING::CoordType & _p, 
	const typename SPATIALINDEXING::ScalarType & _maxDist,OBJPTRCONTAINER & _objectPtrs,
	DISTCONTAINER & _distances, POINTCONTAINER & _points)
	{
		typedef vcg::ClosestIterator<SPATIALINDEXING,OBJPOINTDISTFUNCTOR,OBJMARKER> ClosestIteratorType;
		ClosestIteratorType	Cli=ClosestIteratorType(_Si,_getPointDistance);
		Cli.SetMarker(_marker);
		Cli.Init(_p,_maxDist);
		unsigned int i=0;
		_objectPtrs.clear();
		_distances.clear();
		_points.clear();
		while ((!Cli.End())&&(i<_k))
		{
			_objectPtrs.push_back(&(*Cli));
			_distances.push_back(Cli.Dist());
			_points.push_back(Cli.NearestPoint());
			++Cli;
			i++;
		}
		return (i);
	};
		
	template <class SPATIALINDEXING,class OBJRAYISECTFUNCTOR, class OBJMARKER>
	typename SPATIALINDEXING::ObjPtr GridDoRay(SPATIALINDEXING &_Si,OBJRAYISECTFUNCTOR &_rayIntersector, 
	OBJMARKER &_marker, const Ray3<typename SPATIALINDEXING::ScalarType> & _ray, 
	const typename SPATIALINDEXING::ScalarType & _maxDist,typename SPATIALINDEXING::ScalarType & _t) 
	{
		typedef vcg::RayIterator<SPATIALINDEXING,OBJRAYISECTFUNCTOR,OBJMARKER> RayIteratorType;
		RayIteratorType RayIte=RayIteratorType(_Si,_rayIntersector);
		RayIte.SetMarker(_marker);
		RayIte.Init(_ray);

		if (!RayIte.End())
		{
			_t=RayIte.Dist();
			return(&(*RayIte));
		}
		return 0;
	}
}//end namespace vcg
#endif

