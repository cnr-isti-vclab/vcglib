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
Revision 1.11  2007/05/02 13:25:45  zifnab1974
only use typename when necessary

Revision 1.10  2007/04/10 22:46:57  pietroni
- line 152 changed call intersection to IntersectionPlaneTriangle because changing in function's name

Revision 1.9  2007/01/03 15:51:28  pietroni
added initial define and included missing files

Revision 1.8  2006/01/19 14:06:37  spinelli
add std:: namespace...

Revision 1.7  2005/10/03 16:18:15  spinelli
add template parameter for spatialindexing struction

Revision 1.6  2005/05/30 09:11:20  ganovelli
header added, error in include

Revision 1.3  2005/05/17 21:19:37  ganovelli
some std::and typename  missing  (CRS4)

Revision 1.2  2005/03/08 14:42:22  ganovelli
added vcg header


****************************************************************************/

#include<vector>
#include <algorithm>
#include<vcg/space/point3.h>
#include<vcg/space/plane3.h>
#include<vcg/space/segment3.h>
#include<vcg/space/intersection3.h>
#include<vcg/complex/edgemesh/allocate.h>
#include<vcg/complex/trimesh/allocate.h>
#include<vcg/complex/trimesh/subset.h>
#include<vcg/complex/trimesh/closest.h>
#include<vcg/complex/trimesh/base.h>

#ifndef __VCGLIB_INTERSECTION_TRI_MESH
#define __VCGLIB_INTERSECTION_TRI_MESH

namespace vcg{

/** \addtogroup complex */
/*@{*/
/** 
    Function computing the intersection between  a grid and a plane. It returns all the cells intersected
*/
template < typename  GridType,typename ScalarType>
bool Intersect(   GridType & grid,Plane3<ScalarType> plane, std::vector<typename GridType::Cell *> &cells){					
	Point3d p,_d;
	Plane3d pl;
	_d.Import(plane.Direction());
	pl.SetDirection(_d);
	pl.SetOffset(plane.Offset());
	for( int ax = 0; ax <3; ++ax)
			{ int axis = ax;
				int axis0 = (axis+1)%3;
				int axis1 = (axis+2)%3;
				int i,j;
				Point3i pi;

				Segment3<double> seg;
				seg.P0().Import(grid.bbox.min);
				seg.P1().Import(grid.bbox.min);
				seg.P1()[axis] = grid.bbox.max[axis];

				for(i = 0 ; i <= grid.siz[axis0]; ++i){
					for(j = 0 ; j <= grid.siz[axis1]; ++j)
						{
							seg.P0()[axis0] = grid.bbox.min[axis0]+ (i+0.01) * grid.voxel[axis0] ;
							seg.P1()[axis0] = grid.bbox.min[axis0]+ (i+0.01) * grid.voxel[axis0];
							seg.P0()[axis1] = grid.bbox.min[axis1]+ (j+0.01) * grid.voxel[axis1];
							seg.P1()[axis1] = grid.bbox.min[axis1]+ (j+0.01) * grid.voxel[axis1];
              if ( IntersectionPlaneSegmentEpsilon(pl,seg,p))
								{
									pi[axis] =	std::min(std::max(0,(int)floor((p[axis ]-grid.bbox.min[axis])/grid.voxel[axis])),grid.siz[axis]);
									pi[axis0] = i;
									pi[axis1] = j;
									grid.Grid(pi,axis,cells);
								}
						}
					}
			}
		sort(cells.begin(),cells.end());
		cells.erase(unique(cells.begin(),cells.end()),cells.end());
		
		return false;
	}

/*@}*/
/** 
	Basic Function computing the intersection between  a trimesh and a plane, provided a pointer 
	to an space indexing data structure (e.g. a grid, an oct-tree..)
*/
	template < typename  TriMeshType, typename EdgeMeshType, class ScalarType, class IndexingType >
	bool Intersection(	/*TriMeshType & m, */
	Plane3<ScalarType>  pl,
	EdgeMeshType & em,
	double& ave_length,
	IndexingType *grid,
	typename std::vector< typename IndexingType::Cell* >& cells)
{
		typedef typename TriMeshType::FaceContainer FaceContainer;
		typedef IndexingType GridType;
		typename EdgeMeshType::VertexIterator vi;
		typename TriMeshType::FaceIterator fi;
		std::vector<typename TriMeshType::FaceType*> v;
		v.clear();
		Intersect(*grid,pl,cells);
		Segment3<ScalarType> seg;
		ave_length = 0.0;
		typename std::vector<typename GridType::Cell*>::iterator ic;
		typename GridType::Cell fs,ls;
		for(ic = cells.begin(); ic != cells.end();++ic)
		{
			grid->Grid(*ic,fs,ls);
			typename GridType::Link * lk = fs;
			while(lk != ls){
				typename TriMeshType::FaceType & face = *(lk->Elem());
				if(!face.IsS())
				{
					face.SetS();
					v.push_back(&face);
					if(vcg::IntersectionPlaneTriangle(pl,face,seg))// intersezione piano triangolo
					{
						face.SetS();
						// add to em
						ave_length+=seg.Length();
						vcg::edg::Allocator<EdgeMeshType>::AddEdges(em,1);
						vi = vcg::edg::Allocator<EdgeMeshType>::AddVertices(em,2);
						(*vi).P() = seg.P0();
						em.edges.back().V(0) = &(*vi); 
						vi++;
						(*vi).P() = seg.P1();
						em.edges.back().V(1) = &(*vi); 
					}
				}//endif 
				lk++;
			}//end while
		}
		ave_length/=em.en;
		typename std::vector<typename TriMeshType::FaceType*>::iterator v_i;
    for(v_i=v.begin(); v_i!=v.end(); ++v_i) (*v_i)->ClearS();
		
		return true;
}

/** \addtogroup complex */
/*@{*/
/** 
    Basic Function computing the intersection between  a trimesh and a plane. It returns an EdgeMesh without needing anything else.
		Note: This version always returns a segment for each triangle of the mesh which intersects with the plane. In other
		words there are 2*n vertices where n is the number of segments fo the mesh. You can run vcg::edge:Unify to unify
		the vertices closer that a given value epsilon. Note that, due to subtraction error during triangle plane intersection,
		it is not safe to put epsilon to 0. 
// TODO si dovrebbe considerare la topologia face-face della trimesh per derivare quella della edge mesh..
*/
template < typename  TriMeshType, typename EdgeMeshType, class ScalarType >
bool Intersection(TriMeshType & m,
									Plane3<ScalarType>  pl,
									EdgeMeshType & em)
{
	typename EdgeMeshType::VertexIterator vi;
	typename TriMeshType::FaceIterator fi;
  em.Clear();
	Segment3<ScalarType> seg;
	for(fi=m.face.begin();fi!=m.face.end();++fi) 
		if(!(*fi).IsD())
		{
			if(vcg::IntersectionPlaneTriangle(pl,*fi,seg))// intersezione piano triangolo
							{
								vcg::edg::Allocator<EdgeMeshType>::AddEdges(em,1);
								vi = vcg::edg::Allocator<EdgeMeshType>::AddVertices(em,2);
								(*vi).P() = seg.P0();
								em.edges.back().V(0) = &(*vi); 
								vi++;
								(*vi).P() = seg.P1();
								em.edges.back().V(1) = &(*vi); 
							}
		 }//end for
	
	return true;
}

/** \addtogroup complex */
/*@{*/
/** 
    Compute the intersection between a trimesh and a plane. 
		given a plane return the set of faces that are contained 
		into intersected cells.
*/
template < typename  TriMeshType, class ScalarType, class IndexingType >
bool Intersection(Plane3<ScalarType>  pl,
									IndexingType *grid,
									typename std::vector<typename TriMeshType::FaceType*> &v)
{
	typedef typename TriMeshType::FaceContainer FaceContainer;
	typedef IndexingType GridType;
	typename TriMeshType::FaceIterator fi;
	v.clear();
	typename std::vector< typename GridType::Cell* > cells;
	Intersect(*grid,pl,cells);
	typename std::vector<typename GridType::Cell*>::iterator ic;
	typename GridType::Cell fs,ls;

	for(ic = cells.begin(); ic != cells.end();++ic)
	{
		grid->Grid(*ic,fs,ls);
		typename GridType::Link * lk = fs;
		while(lk != ls){
			typename TriMeshType::FaceType & face = *(lk->Elem());
			v.push_back(&face);
			lk++;
		}//end while
	}//end for
	return true;
}

/** 
	 Computes the intersection between a Ray and a Mesh. Returns a 3D Pointset.  
*/
template < typename  TriMeshType, class ScalarType>
bool IntersectionRayMesh(	
	/* Input Mesh */		TriMeshType * m, 
	/* Ray */				const Line3<ScalarType> & ray,
	/* Intersect Point */	Point3<ScalarType> & hitPoint)
{
	//typedef typename TriMeshType::FaceContainer FaceContainer;
	typename TriMeshType::FaceIterator fi;
	bool hit=false;

	if(m==0) return false;

	//TriMeshType::FaceIterator fi;
	//std::vector<TriMeshType::FaceType*>::iterator fi;

	ScalarType bar1,bar2,dist;
	Point3<ScalarType> p1;
	Point3<ScalarType> p2;
	Point3<ScalarType> p3;
	for(fi = m->face.begin(); fi != m->face.end(); ++fi)
	{
		p1=vcg::Point3<ScalarType>( (*fi).P(0).X() ,(*fi).P(0).Y(),(*fi).P(0).Z() );
		p2=vcg::Point3<ScalarType>( (*fi).P(1).X() ,(*fi).P(1).Y(),(*fi).P(1).Z() );
		p3=vcg::Point3<ScalarType>( (*fi).P(2).X() ,(*fi).P(2).Y(),(*fi).P(2).Z() );
		if(IntersectionLineTriangle<ScalarType>(ray,p1,p2,p3,dist,bar1,bar2))
		{
			hitPoint= p1*(1-bar1-bar2) + p2*bar1 + p3*bar2;
			hit=true;
		}
	}

	return hit;
}

/** 
	 Computes the intersection between a Ray and a Mesh. Returns a 3D Pointset, baricentric's coordinates 
	 and a pointer of intersected face.
*/
template < typename  TriMeshType, class ScalarType>
bool IntersectionRayMesh(	
	/* Input Mesh */		TriMeshType * m, 
	/* Ray */				const Line3<ScalarType> & ray,
	/* Intersect Point */	Point3<ScalarType> & hitPoint,
	/* Baricentric coord 1*/ ScalarType &bar1,
	/* Baricentric coord 2*/ ScalarType &bar2,
	/* Baricentric coord 3*/ ScalarType &bar3,
	/* FacePointer */ typename TriMeshType::FacePointer fp
	)
{
	//typedef typename TriMeshType::FaceContainer FaceContainer;
	typename TriMeshType::FaceIterator fi;
	bool hit=false;

	if(m==0) return false;

	//TriMeshType::FaceIterator fi;
	//std::vector<TriMeshType::FaceType*>::iterator fi;

	ScalarType dist;
	Point3<ScalarType> p1;
	Point3<ScalarType> p2;
	Point3<ScalarType> p3;
	for(fi = m->face.begin(); fi != m->face.end(); ++fi)
	{
		p1=vcg::Point3<ScalarType>( (*fi).P(0).X() ,(*fi).P(0).Y(),(*fi).P(0).Z() );
		p2=vcg::Point3<ScalarType>( (*fi).P(1).X() ,(*fi).P(1).Y(),(*fi).P(1).Z() );
		p3=vcg::Point3<ScalarType>( (*fi).P(2).X() ,(*fi).P(2).Y(),(*fi).P(2).Z() );
		if(IntersectionLineTriangle<ScalarType>(ray,p1,p2,p3,dist,bar1,bar2))
		{
			bar3 = (1-bar1-bar2);
			hitPoint= p1*bar3 + p2*bar1 + p3*bar2;
			fp = &(*fi);
			hit=true;
		}
	}

	return hit;
}

/** 
    Compute the intersection between a mesh and a ball. 
		given a mesh return a new mesh made by a copy of all the faces entirely includeded in the ball plus
		new faces created by refining the ones intersected by the ball border.
		It works by recursively splitting the triangles that cross the border, as long as their area is greater than
		a given value tol. If no value is provided, 1/10^5*2*pi*radius is used 
		NOTE: the returned mesh is a triangle soup 
*/
template < typename  TriMeshType, class ScalarType>
void IntersectionBallMesh(	 TriMeshType & m, const vcg::Sphere3<ScalarType> &ball, TriMeshType & res,
													float tol = 0){

	typename TriMeshType::VertexIterator v0,v1,v2;
	typename TriMeshType::FaceIterator fi;
	std::vector<typename TriMeshType:: FaceType*> closests;
	vcg::Point3<ScalarType>	witness;
	std::pair<ScalarType, ScalarType> info;

	if(tol == 0) tol = M_PI * ball.Radius() * ball.Radius() / 100000;

	for(fi = m.face.begin(); fi != m.face.end(); ++fi)
	if(!(*fi).IsD() && IntersectionSphereTriangle<ScalarType>(ball  ,(*fi), witness , &info))
		closests.push_back(&(*fi));

	res.Clear();
	SubSet(res,closests);
	int i =0;
	while(i<res.fn){
		 bool allIn = ( ball.IsIn(res.face[i].P(0)) && ball.IsIn(res.face[i].P(1))&&ball.IsIn(res.face[i].P(2)));
		if( IntersectionSphereTriangle<ScalarType>(ball  ,res.face[i], witness , &info) && !allIn){
				if(vcg::DoubleArea(res.face[i]) > tol)
				{
				// split the face res.face[i] in four, add the four new faces to the mesh and delete the face res.face[i]
				v0 = vcg::tri::Allocator<TriMeshType>::AddVertices(res,3);	
				fi = vcg::tri::Allocator<TriMeshType>::AddFaces(res,4);	
				
				v1 = v0; ++v1;
				v2 = v1; ++v2;
				(*v0).P() = (res.face[i].P(0) + res.face[i].P(1))*0.5;
				(*v1).P() = (res.face[i].P(1) + res.face[i].P(2))*0.5;
				(*v2).P() = (res.face[i].P(2) + res.face[i].P(0))*0.5;

				(*fi).V(0) = res.face[i].V(0);
				(*fi).V(1) = &(*v0);
				(*fi).V(2) = &(*v2);	
				++fi;

				(*fi).V(0) = res.face[i].V(1);
				(*fi).V(1) = &(*v1);
				(*fi).V(2) = &(*v0);	
				++fi;

				(*fi).V(0) = &(*v0);
				(*fi).V(1) = &(*v1);
				(*fi).V(2) = &(*v2);	
				++fi;

				(*fi).V(0) = &(*v2);
				(*fi).V(1) = &(*v1);
				(*fi).V(2) = res.face[i].V(2) ;	

				vcg::tri::Allocator<TriMeshType>::DeleteFace(res,res.face[i]);
			}
		}// there was no intersection with the boundary

	if(info.first > 0.0) // closest point - radius. If >0 is outside
		vcg::tri::Allocator<TriMeshType>::DeleteFace(res,res.face[i]);
	++i;
	}
}


template < typename  TriMeshType, class ScalarType, class IndexingType>
void IntersectionBallMesh( IndexingType * grid,	 TriMeshType & m, const vcg::Sphere3<ScalarType> &ball, TriMeshType & res,
													float tol = 0){

	typename TriMeshType::VertexIterator v0,v1,v2;
	typename std::vector<typename TriMeshType::FacePointer >::iterator  cfi;
	typename TriMeshType::FaceIterator fi;
	std::vector<typename TriMeshType:: FaceType*> closestsF,closests;
	vcg::Point3<ScalarType>	witness;
	std::vector<vcg::Point3<ScalarType> > witnesses;
	std::vector<ScalarType>	distances;
	std::pair<ScalarType, ScalarType> info;

	if(tol == 0) tol = M_PI * ball.Radius() * ball.Radius() / 100000;

	vcg::tri::GetInSphereFace(m,*grid, ball.Center(), ball.Radius(),closestsF,distances,witnesses);
	for(cfi =closestsF.begin(); cfi != closestsF.end(); ++cfi)
	if(!(**cfi).IsD() && IntersectionSphereTriangle<ScalarType>(ball  ,(**cfi), witness , &info))
		closests.push_back(&(**cfi));

	res.Clear();
	SubSet(res,closests);
	int i =0;
	while(i<res.fn){
		 bool allIn = ( ball.IsIn(res.face[i].P(0)) && ball.IsIn(res.face[i].P(1))&&ball.IsIn(res.face[i].P(2)));
		if( IntersectionSphereTriangle<ScalarType>(ball  ,res.face[i], witness , &info) && !allIn){
				if(vcg::DoubleArea(res.face[i]) > tol)
				{
				// split the face res.face[i] in four, add the four new faces to the mesh and delete the face res.face[i]
				v0 = vcg::tri::Allocator<TriMeshType>::AddVertices(res,3);	
				fi = vcg::tri::Allocator<TriMeshType>::AddFaces(res,4);	
				
				v1 = v0; ++v1;
				v2 = v1; ++v2;
				(*v0).P() = (res.face[i].P(0) + res.face[i].P(1))*0.5;
				(*v1).P() = (res.face[i].P(1) + res.face[i].P(2))*0.5;
				(*v2).P() = (res.face[i].P(2) + res.face[i].P(0))*0.5;

				(*fi).V(0) = res.face[i].V(0);
				(*fi).V(1) = &(*v0);
				(*fi).V(2) = &(*v2);	
				++fi;

				(*fi).V(0) = res.face[i].V(1);
				(*fi).V(1) = &(*v1);
				(*fi).V(2) = &(*v0);	
				++fi;

				(*fi).V(0) = &(*v0);
				(*fi).V(1) = &(*v1);
				(*fi).V(2) = &(*v2);	
				++fi;

				(*fi).V(0) = &(*v2);
				(*fi).V(1) = &(*v1);
				(*fi).V(2) = res.face[i].V(2) ;	

				vcg::tri::Allocator<TriMeshType>::DeleteFace(res,res.face[i]);
			}
		}// there was no intersection with the boundary

	if(info.first > 0.0) // closest point - radius. If >0 is outside
		vcg::tri::Allocator<TriMeshType>::DeleteFace(res,res.face[i]);
	++i;
	}
}

/*@}*/
} // end namespace vcg
#endif
