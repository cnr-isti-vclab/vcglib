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


#ifndef __VCG_VERTEXMESH_CLOSEST
#define __VCG_VERTEXMESH_CLOSEST
#include <math.h>

#include <vcg/space/point3.h>
#include <vcg/space/box3.h>
#include <vcg/space/point4.h>
#include <vcg/math/base.h>
#include <vcg/simplex/face/distance.h>
#include <vcg/simplex/vertex/distance.h>
#include <vcg/space/intersection3.h>
#include <vcg/space/index/space_iterators.h>

namespace vcg {
	namespace vrt {

		//**MARKER CLASSES**//
		template <class MESH_TYPE,class OBJ_TYPE>
		class Tmark
		{
			MESH_TYPE *m;
		public:
			Tmark(){}
			void UnMarkAll(){m->UnMarkAll();}
			bool IsMarked(OBJ_TYPE* obj){return (m->IsMarked(obj));}
			void Mark(OBJ_TYPE* obj){m->Mark(obj);}
			void SetMesh(MESH_TYPE *_m)
			{m=_m;}
		};


		template <class MESH_TYPE>
		class VertTmark:public Tmark<MESH_TYPE,typename MESH_TYPE::VertexType>
		{};

		//**CLOSEST FUNCTION DEFINITION**//

		/*

		aka MetroCore
		data una mesh m e una ug sulle sue facce trova il punto di m piu' vicino ad
		un punto dato.
		*/

		// input: mesh, punto, griglia (gr), distanza limite (mdist)
		// output: normale (interpolata) alla faccia e punto piu' vicino su di essa, e coord baricentriche del punto trovato

		// Nota che il parametro template GRID non ci dovrebbe essere, visto che deve essere 
		// UGrid<MESH::FaceContainer >, ma non sono riuscito a definirlo implicitamente 

		
		template <class MESH, class GRID>
			typename MESH::VertexType * GetClosestVertex( MESH & mesh,GRID & gr,const typename GRID::CoordType & _p, 
			const typename GRID::ScalarType & _maxDist,typename GRID::ScalarType & _minDist )
		{
			typedef typename GRID::ScalarType ScalarType;
			typedef Point3<ScalarType> Point3x;
			typedef VertTmark<MESH> MarkerVert;
			MarkerVert mv;
			mv.SetMesh(&mesh);
			typedef  PointDistanceFunctor VDistFunct;
			_minDist=_maxDist;
			Point3x _closestPt;
			return (gr.GetClosest(PointDistanceFunctor(),mv,_p,_maxDist,_minDist,_closestPt));
		}

		
		template <class MESH, class GRID, class OBJPTRCONTAINER,class DISTCONTAINER, class POINTCONTAINER>
			unsigned int GetKClosestVertex(MESH & mesh,GRID & gr, const unsigned int _k, 
			const typename GRID::CoordType & _p, const typename GRID::ScalarType & _maxDist,
			OBJPTRCONTAINER & _objectPtrs,DISTCONTAINER & _distances, POINTCONTAINER & _points)
		{
			typedef VertTmark<MESH> MarkerVert;
			MarkerVert mv;
			mv.SetMesh(&mesh);
			typedef vcg::vert::PointDistanceFunctor VDistFunct;
			return (gr.GetKClosest/*<VDistFunct,MarkerVert,OBJPTRCONTAINER,DISTCONTAINER,POINTCONTAINER>*/
				(VDistFunct(),mv,_k,_p,_maxDist,_objectPtrs,_distances,_points));
		}

		template <class MESH, class GRID, class OBJPTRCONTAINER, class DISTCONTAINER, class POINTCONTAINER>
			unsigned int GetInSphereVertex(MESH & mesh,
			GRID & gr,
			const typename GRID::CoordType & _p,
			const typename GRID::ScalarType & _r,
			OBJPTRCONTAINER & _objectPtrs,
			DISTCONTAINER & _distances, 
			POINTCONTAINER & _points)
		{
			typedef VertTmark<MESH> MarkerVert;
			MarkerVert mv;
			mv.SetMesh(&mesh);
			typedef vcg::vert::PointDistanceFunctor VDistFunct;
			return (gr.GetInSphere/*<VDistFunct,MarkerVert,OBJPTRCONTAINER,DISTCONTAINER,POINTCONTAINER>*/
				(VDistFunct(),mv,_p,_r,_objectPtrs,_distances,_points));
		}


		template <class MESH, class GRID, class OBJPTRCONTAINER>
			unsigned int GetInBoxVertex(MESH & mesh,
			GRID & gr,
			const vcg::Box3<typename GRID::ScalarType> _bbox,
			OBJPTRCONTAINER & _objectPtrs) 
		{
			typedef VertTmark<MESH> MarkerVert;
			MarkerVert mv;
			mv.SetMesh(&mesh);
			return(gr.GetInBox(mv,_bbox,_objectPtrs));
		}


		//**ITERATORS DEFINITION**//

		template <class GRID,class MESH>
		class ClosestVertexIterator:public vcg::ClosestIterator<GRID,vcg::vert::PointDistanceFunctor, VertTmark<MESH> >
		{
		public:
			typedef GRID GridType;
			typedef MESH MeshType;
			typedef VertTmark<MESH> MarkerVert;
			typedef vcg::vert::PointDistanceFunctor VDistFunct;
			typedef vcg::ClosestIterator<GRID,VDistFunct, VertTmark<MESH> > ClosestBaseType;

			ClosestVertexIterator(GridType &_Si):ClosestBaseType(_Si,VDistFunct()){}

			void SetMesh(MeshType *m)
			{tm.SetMesh(m);}
		};


	}	 // end namespace vert
}	 // end namespace vcg

#endif
