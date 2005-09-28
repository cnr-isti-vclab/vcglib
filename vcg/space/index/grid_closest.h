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
Revision 1.1  2005/09/27 15:09:38  cignoni
First Version


****************************************************************************/

///** Returns the closest posistion of a point p and its distance
//@param p a 3d point
//@param max_dist maximum distance not to search beyond.
//@param dist_funct (templated type) a functor object used to calculate distance from a grid object to the point p.
//@param dist the returned closest distance
//@param res the returned closest point
//@return The closest element
//*/
///*
//	A DISTFUNCT object must implement an operator () with signature:
//		bool operator () (const ObjType& obj, const CoordType & p, ScalarType & min_dist, CoordType & res);
//*/
#ifndef __VCGLIB_GRID_CLOSEST
#define __VCGLIB_GRID_CLOSEST


namespace vcg{

	template <class SPATIAL_INDEX,class DISTFUNCTOR, class TMARKER>
		typename SPATIAL_INDEX::ObjPtr  GetClosest( const typename SPATIAL_INDEX::CoordType & p, 
		                                            const typename SPATIAL_INDEX::ScalarType & max_dist, 
		                                            DISTFUNCTOR & dist_funct, 
                                                typename SPATIAL_INDEX::ScalarType & dist,
		                                            typename SPATIAL_INDEX:: CoordType & res,
                                                TMARKER tm,
                                                SPATIAL_INDEX &Si)
	{
		typedef SPATIAL_INDEX::ObjPtr ObjPtr;
		typedef SPATIAL_INDEX SpatialIndex;
		typedef SPATIAL_INDEX::CoordType CoordType;
		typedef SPATIAL_INDEX::ScalarType ScalarType;

		// Initialize min_dist with max_dist to exploit early rejection test.
		dist = max_dist;

		ScalarType dx = ( (p[0]-Si.bbox.min[0])/Si.voxel[0] );
		ScalarType dy = ( (p[1]-Si.bbox.min[1])/Si.voxel[1] );
		ScalarType dz = ( (p[2]-Si.bbox.min[2])/Si.voxel[2] );

		int ix = int( dx );
		int iy = int( dy );
		int iz = int( dz );

    if (ix<0) ix=0;
    if (iy<0) iy=0;
    if (iz<0) iz=0;
    if (ix>=Si.siz[0]-1) ix=Si.siz[0]-1;
    if (iy>=Si.siz[1]-1) iy=Si.siz[1]-1;
    if (iz>=Si.siz[2]-1) iz=Si.siz[2]-1;

    if (!Si.bbox.IsIn(p)){
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

		tm.UnMarkAll();

		SpatialIndex::CellIterator first,last;
		SpatialIndex::CellIterator l;

		if ((ix>=0) && (iy>=0) && (iz>=0) && 
			(ix<Si.siz[0]) && (iy<Si.siz[1]) && (iz<Si.siz[2])) {

				Si.Grid( ix, iy, iz, first, last );
				for(l=first;l!=last;++l)
					if (!(**l).IsD())
					{
						ObjPtr elem=&(**l);
						if(!tm.IsMarked(elem))
						{
							//if (!l->Elem()->IsD() && l->Elem()->Dist(p,min_dist,t_res)) {
							//if (!l->Elem()->IsD() && dist_funct(*(l->Elem()), p, min_dist, t_res)) { // <-- NEW: use of distance functor
							if (dist_funct((**l), p,dist, t_res))  // <-- NEW: use of distance functor
							{
								winner=elem;
								res=t_res;
							}
							tm.Mark(elem);
						}
					}
			};

		//return winner;

		// the portion of the grid that have already been checked.
    Point3i done_min, done_max; 
		// the new box that we want to traverse: todo is a superset of done.
    Point3i todo_min=Point3i(ix,iy,iz), todo_max=Point3i(ix,iy,iz);


		// we should traverse only (todo - done).
		while (dist>radius) {
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
								if( ! tm.IsMarked(elem))
								{
									//if (!l->Elem()->IsD() && l->Elem()->Dist(p,min_dist,t_res)) {
									if (dist_funct((**l), p, dist, t_res)) 
									{
										winner=elem;
										res=t_res;
									};
									tm.Mark(elem);
								}
							}
						};
					}
		};
		return winner;
	};

}//end namespace vcg
#endif

