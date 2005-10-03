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
Revision 1.14  2005/09/30 13:10:37  pietroni
used functor defined in face/distance.h for distance point-face
used functor defined in intersection3.h for ray-triangle intersection
added GetKClosest and DoRay Functions

Revision 1.13  2005/09/28 08:30:48  cignoni
changed name of include, removed use of an undefined type (scalar instead of Scalar)
removed unused code portions (the old closest code)

Revision 1.12  2005/09/21 09:24:30  pietroni
Added RayIterators.
Added ClosestIterators on Triangles and Vertices.
Added Closest Functions on triangles and Vertices.

Revision 1.11  2005/09/19 13:36:24  pietroni
added ray iterator of faces

Revision 1.10  2005/09/16 11:53:51  cignoni
Small gcc compliling issues

Revision 1.9  2005/09/15 13:16:10  spinelli
fixed bugs

Revision 1.8  2005/09/15 11:15:00  pietroni
minor changes

Revision 1.7  2005/09/14 12:56:47  pietroni
used closest function from grid

Revision 1.6  2005/08/26 09:12:48  cignoni
changed  typedef A2UGridLink  da 'GridStaticPtr<MESH::FaceContainer,double>::Link' a  typedef 'GRID::Link'

Revision 1.5  2005/02/08 17:49:38  pietroni
added  if (!l->Elem()->IsD()) test on each element

Revision 1.4  2005/01/28 12:00:33  cignoni
small gcc compiling issues for namespaces

Revision 1.3  2005/01/24 11:47:23  cignoni
Now used also by the official Metro
Removed using namespace (NEVER IN HEADERS!)
Made  the computation of barycentric coords only when necessary
Renamed Mindistpoint to Closest

Revision 1.2  2005/01/21 17:13:09  pietroni
included distance.h changed Dist to  vcg::face::PointDistance

Revision 1.1  2004/10/04 15:32:16  ganovelli
moved from metro core

Revision 1.6  2004/05/14 00:34:36  ganovelli
header added

****************************************************************************/

#ifndef __VCG_TRIMESH_CLOSEST
#define __VCG_TRIMESH_CLOSEST
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
	namespace trimesh {

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
		class FaceTmark:public Tmark<MESH_TYPE,typename MESH_TYPE::FaceType>
		{};

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
			typename MESH::FaceType * GetClosestFace( MESH & mesh,GRID & gr,const typename GRID::CoordType & _p, 
			const typename GRID::ScalarType & _maxDist,typename GRID::ScalarType & _minDist,
			typename GRID::CoordType _closestPt,typename GRID::CoordType & _normf, 
			typename GRID::CoordType & _ip)
		{
			typedef GRID::ScalarType ScalarType;
			typedef Point3<ScalarType> Point3x;

			typedef FaceTmark<MESH> MarkerFace;
			MarkerFace mf;
			mf.SetMesh(&mesh);
			typedef vcg::face::PointDistanceFunctor FDistFunct;
			_minDist=_maxDist;
			typename MESH::FaceType* bestf= gr.GetClosest<FDistFunct,MarkerFace>(FDistFunct(),mf,_p,_maxDist,_minDist,_closestPt);

			if(_maxDist> ScalarType(fabs(_minDist)))
			{
				// f=bestf;
				typename MESH::ScalarType alfa, beta, gamma;
				//calcolo normale con interpolazione trilineare
				bestf->InterpolationParameters(_closestPt, alfa, beta, gamma);
				_normf =  (bestf->V(0)->cN())*alfa+
					(bestf->V(1)->cN())*beta+
					(bestf->V(2)->cN())*gamma ;
				_ip=Point3x(alfa,beta,gamma);
				//normf.Normalize(); inutile si assume le normali ai vertici benfatte										

				_minDist = fabs(_minDist);
				return(bestf);
			}
			return (0);
		}

		template <class MESH, class GRID>
			typename MESH::FaceType * GetClosestFace( MESH & mesh,GRID & gr,const typename GRID::CoordType & _p, 
			const typename GRID::ScalarType & _maxDist,typename GRID::ScalarType & _minDist,
			typename GRID::CoordType _closestPt,typename GRID::CoordType & _normf)
		{
			Point3<SCALAR> _ip;
			return (GetClosest<MESH,GRID>(mesh,gr,_p,_maxDist,_minDist,_closestPt,_normf,ip));
		}

		template <class MESH, class GRID>
			typename MESH::FaceType * GetClosestFace( MESH & mesh,GRID & gr,const typename GRID::CoordType & _p, 
			const typename GRID::ScalarType & _maxDist,typename GRID::ScalarType & _minDist,
			typename GRID::CoordType _closestPt)
		{
			typedef GRID::ScalarType ScalarType;
			typedef Point3<ScalarType> Point3x;
			typedef FaceTmark<MESH> MarkerFace;
			MarkerFace mf;
			mf.SetMesh(&mesh);
			typedef vcg::face::PointDistanceFunctor PDistFunct;
			_minDist=_maxDist;
			return (gr.GetClosest<FDistFunct,MarkerFace>(PDistFunct(),mf,_p,_maxDist,_minDist,_closestPt));
		}

		template <class MESH, class GRID>
			typename MESH::VertexType * GetClosestVertex( MESH & mesh,GRID & gr,const typename GRID::CoordType & _p, 
			const typename GRID::ScalarType & _maxDist,typename GRID::ScalarType & _minDist,
			typename GRID::CoordType _closestPt,typename GRID::CoordType & _normf)
		{
			typedef GRID::ScalarType ScalarType;
			typedef Point3<ScalarType> Point3x;
			typedef VertTmark<MESH> MarkerVert;
			MarkerVert mv;
			mv.SetMesh(&mesh);
			typedef vcg::vertex::PointDistanceFunctor VDistFunct;
			_minDist=_maxDist;
			return (gr.GetClosest<VDistFunct,MarkerFace>(VDistFunct(),mv,_p,_maxDist,_minDist,_closestPt));
		}

		template <class MESH, class GRID, class OBJPTRCONTAINER,class DISTCONTAINER, class POINTCONTAINER>
			unsigned int GetKClosestFace(MESH & mesh,GRID & gr, const unsigned int _k, 
			const typename GRID::CoordType & _p, const typename GRID::ScalarType & _maxDist,
			OBJPTRCONTAINER & _objectPtrs,DISTCONTAINER & _distances, POINTCONTAINER & _points)
		{
			typedef FaceTmark<MESH> MarkerFace;
			MarkerFace mf;
			mf.SetMesh(&mesh);
			typedef vcg::face::PointDistanceFunctor FDistFunct;
			return (gr.GetKClosest<FDistFunct,MarkerFace,OBJPTRCONTAINER,DISTCONTAINER,POINTCONTAINER>
				(FDistFunct(),mf,_k,_p,_maxDist,_objectPtrs,_distances,_points));
		}

		template <class MESH, class GRID, class OBJPTRCONTAINER,class DISTCONTAINER, class POINTCONTAINER>
			unsigned int GetKClosestVertex(MESH & mesh,GRID & gr, const unsigned int _k, 
			const typename GRID::CoordType & _p, const typename GRID::ScalarType & _maxDist,
			OBJPTRCONTAINER & _objectPtrs,DISTCONTAINER & _distances, POINTCONTAINER & _points)
		{
			typedef VertTmark<MESH> MarkerVert;
			MarkerVert mv;
			mv.SetMesh(&mesh);
			typedef vcg::vertex::PointDistanceFunctor VDistFunct;
			return (gr.GetKClosest<VDistFunct,MarkerFace,OBJPTRCONTAINER,DISTCONTAINER,POINTCONTAINER>
				(VDistFunct(),mv,_k,_p,_maxDist,_objectPtrs,_distances,_points));
		}

		template <class MESH, class GRID, class OBJPTRCONTAINER, class DISTCONTAINER, class POINTCONTAINER>
			unsigned int GetInSphereFace(MESH & mesh,
			GRID & gr,
			const typename GRID::CoordType & _p,
			const typename GRID::ScalarType & _r,
			OBJPTRCONTAINER & _objectPtrs,
			DISTCONTAINER & _distances, 
			POINTCONTAINER & _points)
		{
			typedef FaceTmark<MESH> MarkerFace;
			MarkerFace mf;
			mf.SetMesh(&mesh);
			typedef vcg::face::PointDistanceFunctor FDistFunct;
			return (gr.GetInSphere<FDistFunct,MarkerFace,OBJPTRCONTAINER,DISTCONTAINER,POINTCONTAINER>
				(FDistFunct(),mf,_p,_r,_maxDist,_objectPtrs,_distances,_points));
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
			typedef vcg::vertex::PointDistanceFunctor VDistFunct;
			return (gr.GetInSphere<VDistFunct,MarkerVert,OBJPTRCONTAINER,DISTCONTAINER,POINTCONTAINER>
				(VDistFunct(),mv,_p,_r,_objectPtrs,_distances,_points));
		}

		template <class MESH, class GRID, class OBJPTRCONTAINER>
			unsigned int GetInBoxFace(MESH & mesh,
			GRID & gr,
			const vcg::Box3<typename GRID::ScalarType> _bbox,
			OBJPTRCONTAINER & _objectPtrs) 
		{
			typedef FaceTmark<MESH> MarkerFace;
			MarkerFace mf;
			mf.SetMesh(&mesh);
			return(gr.GetInBox<MarkerFace,OBJPTRCONTAINER>(mf,_bbox,_objectPtrs));
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
			return(gr.GetInBox<MarkerVert,OBJPTRCONTAINER>(mv,_bbox,_objectPtrs));
		}

		template <class MESH, class GRID>
			typename GRID::ObjPtr DoRay(MESH & mesh,GRID & gr, const Ray3<typename GRID::ScalarType> & _ray,
			const typename GRID::ScalarType & _maxDist, typename GRID::ScalarType & _t) 
		{
			typedef MESH::FaceType FaceType;
			typedef MESH::ScalarType ScalarType;
			typedef FaceTmark<MESH> MarkerFace;
			MarkerFace mf;
			mf.SetMesh(&mesh);
			typedef vcg::RayTriangleIntersectionFunctor<true> FintFunct;
			return(gr.DoRay(FintFunct(),mf,_ray,_maxDist,_t));
		}


		//**ITERATORS DEFINITION**//

		template <class GRID,class MESH>
		class ClosestFaceIterator:public vcg::ClosestIterator<GRID,
			vcg::face::PointDistanceFunctor,typename FaceTmark<MESH> >
		{
		public:
			typedef typename GRID GridType;
			typedef typename MESH MeshType;
			typedef typename FaceTmark<MESH> MarkerFace;
			typedef typename vcg::face::PointDistanceFunctor PDistFunct;
			typedef typename vcg::ClosestIterator<GRID,PDistFunct,typename FaceTmark<MESH> > ClosestBaseType;
			typedef typename MESH::FaceType FaceType;
			typedef typename MESH::ScalarType ScalarType;

			//ClosestFaceIterator(GridType &_Si):ClosestBaseType(_Si,PDistFunct<FaceType,ScalarType>()){}
			ClosestFaceIterator(GridType &_Si):ClosestBaseType(_Si,PDistFunct()){}

			void SetMesh(MeshType *m)
			{tm.SetMesh(m);}
		};

		template <class GRID,class MESH>
		class ClosestVertexIterator:public vcg::ClosestIterator<GRID,vcg::vertex::PointDistanceFunctor,typename VertTmark<MESH> >
		{
		public:
			typedef typename GRID GridType;
			typedef typename MESH MeshType;
			typedef typename VertTmark<MESH> MarkerVert;
			typedef typename vcg::vertex::PointDistanceFunctor VDistFunct;
			typedef typename vcg::ClosestIterator<GRID,VDistFunct,typename VertTmark<MESH> > ClosestBaseType;

			ClosestVertexIterator(GridType &_Si):ClosestBaseType(_Si,VDistFunct()){}

			void SetMesh(MeshType *m)
			{tm.SetMesh(m);}
		};

		template <class GRID,class MESH>
		class TriRayIterator:public vcg::RayIterator<GRID,vcg::RayTriangleIntersectionFunctor<true>,typename FaceTmark<MESH> >
		{
		public:
			typedef typename GRID GridType;
			typedef typename MESH MeshType;
			typedef typename FaceTmark<MESH> MarkerFace;
			typedef typename vcg::RayTriangleIntersectionFunctor<true> FintFunct;
			typedef typename vcg::RayIterator<GRID,FintFunct,typename FaceTmark<MESH> > RayBaseType;
			typedef typename MESH::FaceType FaceType;
			typedef typename MESH::ScalarType ScalarType;

			TriRayIterator(GridType &_Si):RayBaseType(_Si,FintFunct()){}

			void SetMesh(MeshType *m)
			{tm.SetMesh(m);}

		};

	}	 // end namespace trimesh
}	 // end namespace vcg

#endif
