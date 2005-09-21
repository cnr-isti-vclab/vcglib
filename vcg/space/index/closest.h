
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

namespace vcg{

	template <class SPATIAL_INDEX,class DISTFUNCTOR,class TMARKER>
		typename SPATIAL_INDEX::ObjPtr  GetClosest( const typename SPATIAL_INDEX::CoordType & p, 
		const typename SPATIAL_INDEX::ScalarType & max_dist, 
		DISTFUNCTOR & dist_funct, typename SPATIAL_INDEX::ScalarType & dist,
		typename SPATIAL_INDEX:: CoordType & res,TMARKER tm,SPATIAL_INDEX &Si)
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

		if (!Si.bbox.IsIn(p))
			assert (0);///the grid has to be extended until the point

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

		Point3i done_min=Point3i(ix,iy,iz), done_max=Point3i(ix,iy,iz);

		//printf(".");

		while (dist>radius) {
			//if (dy-ScalarType(iy))
			done_min[0]--; if (done_min[0]<0) done_min[0]=0;
			done_min[1]--; if (done_min[1]<0) done_min[1]=0;
			done_min[2]--; if (done_min[2]<0) done_min[2]=0;
			done_max[0]++; if (done_max[0]>=Si.siz[0]-1) done_max[0]=Si.siz[0]-1;
			done_max[1]++; if (done_max[1]>=Si.siz[1]-1) done_max[1]=Si.siz[1]-1;
			done_max[2]++; if (done_max[2]>=Si.siz[2]-1) done_max[2]=Si.siz[2]-1;
			radius+=voxel_min;
			//printf("+");
			for (ix=done_min[0]; ix<=done_max[0]; ix++) 
				for (iy=done_min[1]; iy<=done_max[1]; iy++) 
					for (iz=done_min[2]; iz<=done_max[2]; iz++) 
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
