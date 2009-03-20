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


****************************************************************************/

#ifndef __VCG_EDGEMESH_CLOSEST
#define __VCG_EDGEMESH_CLOSEST
#include <math.h>

#include <vcg/space/point3.h>
#include <vcg/space/box3.h>
#include <vcg/space/point4.h>
#include <vcg/math/base.h>
#include <vcg/simplex/edge/distance.h>
#include <vcg/simplex/vertex/distance.h>
#include <vcg/space/intersection3.h>
#include <vcg/space/index/space_iterators.h>

namespace vcg {
	namespace edgemesh {

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
		class EdgeTmark:public Tmark<MESH_TYPE,typename MESH_TYPE::EdgeType>
		{};

		template <class MESH_TYPE>
		class VertTmark:public Tmark<MESH_TYPE,typename MESH_TYPE::VertexType>
		{};

		//**CLOSEST FUNCTION DEFINITION**//

		template <class MESH, class GRID>
			typename MESH::EdgeType * GetClosestEdge( MESH & mesh,GRID & gr,const typename GRID::CoordType & _p, 
			const typename GRID::ScalarType & _maxDist,typename GRID::ScalarType & _minDist,
			typename GRID::CoordType &_closestPt)
		{
			typedef typename GRID::ScalarType ScalarType;
			typedef Point3<ScalarType> Point3x;
			typedef EdgeTmark<MESH> MarkerEdge;
			MarkerEdge mf;
			mf.SetMesh(&mesh);
			vcg::edge::PointDistanceFunctor<ScalarType> PDistFunct;
			_minDist=_maxDist;
			return (gr.GetClosest(PDistFunct,mf,_p,_maxDist,_minDist,_closestPt));
		}

		template <class MESH, class GRID>
			typename MESH::VertexType * GetClosestVertex( MESH & mesh,GRID & gr,const typename GRID::CoordType & _p, 
			const typename GRID::ScalarType & _maxDist,typename GRID::ScalarType & _minDist )
		{
			typedef typename GRID::ScalarType ScalarType;
			typedef Point3<ScalarType> Point3x;
			typedef VertTmark<MESH> MarkerVert;
			MarkerVert mv;
			mv.SetMesh(&mesh);
			typedef vcg::vertex::PointDistanceFunctor<ScalarType> VDistFunct;
			_minDist=_maxDist;
			Point3x _closestPt;
			return (gr.GetClosest/*<VDistFunct,MarkerVert>*/(VDistFunct(),mv,_p,_maxDist,_minDist,_closestPt));
		}

		template <class MESH, class GRID, class OBJPTRCONTAINER,class DISTCONTAINER, class POINTCONTAINER>
			unsigned int GetKClosestEdge(MESH & mesh,GRID & gr, const unsigned int _k, 
			const typename GRID::CoordType & _p, const typename GRID::ScalarType & _maxDist,
			OBJPTRCONTAINER & _objectPtrs,DISTCONTAINER & _distances, POINTCONTAINER & _points)
		{
			typedef typename GRID::ScalarType ScalarType;
			typedef EdgeTmark<MESH> MarkerEdge;
			MarkerEdge mf;
			mf.SetMesh(&mesh);
		  vcg::face::PointDistanceFunctor<ScalarType> FDistFunct;
			return (gr.GetKClosest /*<vcg::face::PointDistanceFunctor, MarkerFace,OBJPTRCONTAINER,DISTCONTAINER,POINTCONTAINER>*/
				(FDistFunct,mf,_k,_p,_maxDist,_objectPtrs,_distances,_points));
		}

		template <class MESH, class GRID, class OBJPTRCONTAINER,class DISTCONTAINER, class POINTCONTAINER>
			unsigned int GetKClosestVertex(MESH & mesh,GRID & gr, const unsigned int _k, 
			const typename GRID::CoordType & _p, const typename GRID::ScalarType & _maxDist,
			OBJPTRCONTAINER & _objectPtrs,DISTCONTAINER & _distances, POINTCONTAINER & _points)
		{
			typedef typename GRID::ScalarType ScalarType;
			typedef VertTmark<MESH> MarkerVert;
			MarkerVert mv;
			mv.SetMesh(&mesh);
			typedef vcg::vertex::PointDistanceFunctor<ScalarType> VDistFunct;
			return (gr.GetKClosest/* <VDistFunct,MarkerVert,OBJPTRCONTAINER,DISTCONTAINER,POINTCONTAINER>*/
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
			typedef typename GRID::ScalarType ScalarType;
			typedef VertTmark<MESH> MarkerVert;
			MarkerVert mv;
			mv.SetMesh(&mesh);
			typedef vcg::vertex::PointDistanceFunctor<ScalarType> VDistFunct;
			return (gr.GetInSphere/*<VDistFunct,MarkerVert,OBJPTRCONTAINER,DISTCONTAINER,POINTCONTAINER>*/
				(VDistFunct(),mv,_p,_r,_objectPtrs,_distances,_points));
		}

		template <class MESH, class GRID, class OBJPTRCONTAINER>
			unsigned int GetInBoxEdge(MESH & mesh,
			GRID & gr,
			const vcg::Box3<typename GRID::ScalarType> _bbox,
			OBJPTRCONTAINER & _objectPtrs) 
		{
			typedef EdgeTmark<MESH> EdgeTmark;
			EdgeTmark mf;
			mf.SetMesh(&mesh);
			return(gr.GetInBox/*<MarkerFace,OBJPTRCONTAINER>*/(mf,_bbox,_objectPtrs));
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
			return(gr.GetInBox/*<MarkerVert,OBJPTRCONTAINER>*/(mv,_bbox,_objectPtrs));
		}

		
	}	 // end namespace edgemesh
}	 // end namespace vcg

#endif
