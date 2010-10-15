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

#ifndef __VCGLIB_POLYGON_SUPPORT
#define __VCGLIB_POLYGON_SUPPORT

#include <vector>
#include <vcg/complex/trimesh/subset.h>
#include <vcg/simplex/face/jumping_pos.h>
#include <vcg/space/planar_polygon_tessellation.h>

namespace vcg
{
	namespace tri{
		/// \ingroup trimesh 

		/// \headerfile polygon_support.h vcg/complex/trimesh/polygon_support.h

		/// \brief This class is used convert between polygonal meshes and triangular meshes

		/**
		This class contains two members that allow to build a triangular mesh from a polygonal mesh
		and viceversa. In a trimesh, the generic polygons with n sides are codified represented by tagging the internal edge of the face 
		with the SetF.
		*/

		template <class TriMeshType,class PolyMeshType >
		struct PolygonSupport{
			
		/**
		Import a  trianglemesh from a polygon mesh
		**/
		static void ImportFromPolyMesh(TriMeshType & tm,  PolyMeshType & pm){
			std::vector<typename PolyMeshType::CoordType> points;
			std::vector<int> faces;

			// the vertices are the same, simply import them
			typename PolyMeshType::VertexIterator vi;
			typename TriMeshType::FaceIterator tfi,tfib ;
			typename TriMeshType ::VertexIterator tvi = Allocator<TriMeshType>::AddVertices(tm,pm.vert.size());
			int cnt = 0;
			for(tvi = tm.vert.begin(),vi = pm.vert.begin(); tvi != tm.vert.end(); ++tvi,++vi,++cnt)
				if(!(*vi).IsD()) (*tvi).ImportData(*vi); else vcg::tri::Allocator<TriMeshType>::DeleteVertex(tm,(*tvi));

			typename PolyMeshType::FaceIterator fi;
			for(fi = pm.face.begin(); fi != pm.face.end(); ++fi)
			if(!((*fi).IsD())){
			points.clear();
				for(int i  = 0; i < (*fi).VN(); ++i) {
					typename	PolyMeshType::VertexType * v = (*fi).V(i);
					points.push_back(v->P());
				}
				faces.clear();
				TessellatePlanarPolygon3(points,faces);
				tfib = tfi  = Allocator<TriMeshType>::AddFaces(tm,faces.size()/3);
				for(int i = 0; tfi !=  tm.face.end();++tfi){
					(*tfi).V(0) = &tm.vert[ (*fi).V( faces[i]  ) - &(*pm.vert.begin())];
					(*tfi).V(1) = &tm.vert[ (*fi).V( faces[i+1]) - &(*pm.vert.begin())];
					(*tfi).V(2) = &tm.vert[ (*fi).V( faces[i+2]) - &(*pm.vert.begin())];
					// set the F flags
					if( (faces[i]+1)%points.size() != faces[i+1]) (*tfi).SetF(0);
					if( (faces[i+1]+1)%points.size() != faces[i+2]) (*tfi).SetF(1);
					if( (faces[i+2]+1)%points.size() != faces[i]) (*tfi).SetF(2);
					i+=3;
				}

			}
		}
		 

		/**
		Import a polygon mesh from a triangle mesh
		**/
		static void ImportFromTriMesh( PolyMeshType & pm,  TriMeshType & tm){

		// the vertices are the same, simply import them
		int cnt = 0;
		typename TriMeshType ::ConstVertexIterator tvi;
		typename PolyMeshType::VertexIterator vi  = vcg::tri::Allocator<PolyMeshType>::AddVertices(pm,tm.vert.size());
		for(tvi = tm.vert.begin(); tvi != tm.vert.end(); ++tvi,++vi,++cnt) 
			if(!(*tvi).IsD())(*vi).ImportData(*tvi); else vcg::tri::Allocator<PolyMeshType> ::DeleteVertex(pm,(*vi));

		// convert the faces
		typename TriMeshType::FaceIterator tfi;
		vcg::face::JumpingPos<typename TriMeshType::FaceType> p;

		for( tfi = tm.face.begin(); tfi != tm.face.end(); ++tfi) if(!(*tfi).IsD() && !(*tfi).IsV())
		{
			std::vector<typename TriMeshType::VertexPointer> vs;// vertices of the polygon
			std::vector<typename TriMeshType::FacePointer> fs;// triangle faces corresponding to the polygon


			// find a non tagged edge
			int se = 0;
			for(;se < 3;++se) if (!(*tfi).IsF(se)) break;

			// initialize a pos on the first non tagged edge
			typename TriMeshType::VertexPointer v0 = (*tfi).V(se);
			p.F() = &(*tfi);
			p.E() = se;
			p.V() = p.F()->V(p.F()->Next(se));
			p.FlipE();

			vs.push_back(p.F()->V(se));
	
			do{
				while(p.F()->IsF(p.E())) { fs.push_back(p.F()); p.FlipF(); p.FlipE(); p.F()->SetV();} 
				vs.push_back(p.F()->V(p.E()));
				p.FlipV();
				p.FlipE();
			} while( p.V() != v0 );

			//now vs  contains all the vertices of the polygon (still in the trimesh)
			typename PolyMeshType::FaceIterator pfi =  vcg::tri::Allocator<PolyMeshType>::AddFaces(pm,1);
			(*pfi).Alloc(vs.size());
			for( int i  = 0 ; i < vs.size(); ++i)
				(*pfi).V(i) = ( typename PolyMeshType::VertexType*)  & pm.vert[vs[i]-&(*tm.vert.begin())];
			// here handle the other compoenents of the face (how the conponents of the n triangle faces goes in the 
			// the property of the polygon (e.g. the normal, the color, the quality and so on)
			// TODO

		}
	} 
	
	static void ExtractPolygon(typename TriMeshType::FacePointer tfi, std::vector<typename TriMeshType::VertexPointer> &vs){
      vs.clear();
			// find a non tagged edge
			int se = -1;
			for(int i=0; i<3; i++) if (!( tfi->IsF(i))) { se = i; break;}
			
			assert(se!=-1); // else, all faux edges!
			
			// initialize a pos on the first non faux edge
			typename TriMeshType::VertexPointer v0 = tfi->V(se);
			
			vcg::face::JumpingPos<typename TriMeshType::FaceType> p;
			
			p.F() = tfi;
			p.E() = se;
			p.V() = p.F()->V(p.F()->Next(se));
			p.FlipE();

			vs.push_back(p.F()->V(se));
	
	    int guard = 0;
			do{
				while(p.F()->IsF(p.E())) {  p.FlipF(); p.FlipE(); p.F()->SetV(); if (guard++>10) break;} 
				if (guard++>10) break;
				vs.push_back(p.F()->V(p.E()));
				p.FlipV();
				p.FlipE();
			} while( p.V() != v0 );
  }

}; // end of struct
}} // end namespace vcg::tri

#endif // __VCGLIB_TRI_CLIP
