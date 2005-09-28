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
//#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/space/index/space_iterators.h>
#include <vcg/space/index/grid_closest.h>

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

//**INTERSECTION FUNCTIONS**//
  ///class of functor used to calculate the radius-triangle intersection
  template <class FACE_TYPE>
  class FaceIntersection {
	  typedef typename FACE_TYPE FaceType;
	  typedef typename FACE_TYPE::ScalarType ScalarType;
	  typedef typename vcg::Line3<ScalarType> RayType;
	  typedef typename FaceType::CoordType CoordType;

  public: 
	  inline bool operator () (const FaceType & f, RayType r,CoordType & Int) {
		  CoordType p0=f.V(0)->P();
		  CoordType p1=f.V(1)->P();
		  CoordType p2=f.V(2)->P();
		  ///2 sides of normal test 
		  if ((vcg::Intersection<ScalarType>(r,p0,p1,p2,Int))
			  ||(vcg::Intersection<ScalarType>(r,p0,p2,p1,Int)))	
			  return ((r.Direction()*(Int-r.Origin()))>0);///control side of intersection
		  return false;
	  }
  };

//**DISTANCE FUNCTIONS**//
  ///class of functor used to calculate the point triangle distance and nearest point
  template <class FACE_TYPE>
  class FaceDistance {
	  typedef typename FACE_TYPE FaceType;
	  typedef typename FACE_TYPE::ScalarType ScalarType;
	  typedef typename FaceType::CoordType CoordType;

  public: 
	  bool operator () (const FaceType & f, const CoordType & pt, ScalarType & mindist, CoordType & result) {
		  return (vcg::face::PointDistance<FaceType>(f,pt, mindist, result));
	  }
  };
	
  ///class of functor used to calculate the point triangle distance and nearest point
  template <class VERTEX_TYPE>
  class VertexDistance {
	  typedef typename VERTEX_TYPE VertexType;
	  typedef typename VERTEX_TYPE::ScalarType ScalarType;
	  typedef typename VERTEX_TYPE::CoordType CoordType;

  public: 
	  bool operator () (const VertexType & v, const CoordType & pt, ScalarType & dist, CoordType & result) {
	     ScalarType d=dist;
		 dist=(v.P()-pt).Norm();
		 result=pt;
		 return(d>dist);
	  }
  };

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

template <class MESH, class GRID, class SCALAR>
void Closest( MESH & mesh, const Point3<SCALAR> & p, GRID & gr, SCALAR & mdist, 
									Point3<SCALAR> & normf, Point3<SCALAR> & bestq, typename MESH::FaceType * &f, Point3<SCALAR> &ip)
{
   typedef SCALAR ScalarType;
   typedef Point3<ScalarType> Point3x;
   typedef Box3<SCALAR> Box3x;
	
	
  ScalarType error = mdist;
  typedef FaceTmark<MESH> MarkerFace;
  MarkerFace t;
  t.SetMesh(&mesh);
  typedef typename FaceDistance<typename MESH::FaceType> FDistFunct;
  typename MESH::FaceType* bestf= vcg::GetClosest<GRID,FDistFunct,MarkerFace>(p,mdist,FDistFunct() ,error,bestq,t,gr);

  if(mdist > ScalarType(fabs(error)))
  {
	  f=bestf;
	  typename MESH::ScalarType alfa, beta, gamma;
	  //calcolo normale con interpolazione trilineare
	  bestf->InterpolationParameters(bestq, alfa, beta, gamma);
	  normf =  (bestf->V(0)->cN())*alfa+
						 (bestf->V(1)->cN())*beta+
						 (bestf->V(2)->cN())*gamma ;
	  ip=Point3x(alfa,beta,gamma);
	  //normf.Normalize(); inutile si assume le normali ai vertici benfatte										
    
    mdist = ScalarType(fabs(error));
  }
}

template <class MESH, class GRID, class SCALAR>
void Closest( MESH & mesh, const Point3<SCALAR> & p, GRID & gr, SCALAR & mdist, 
									Point3<SCALAR> & normf, Point3<SCALAR> & bestq, typename MESH::FaceType * &f)
{
	Point3<SCALAR> ip;
	Closest<MESH,GRID,SCALAR>(mesh,p,gr,mdist,normf,bestq,f,ip);
}

template <class MESH, class GRID, class SCALAR>
void ClosestVertex( MESH & mesh, const Point3<SCALAR> & p, GRID & gr, SCALAR & mdist, 
									Point3<SCALAR> & normf, Point3<SCALAR> & bestq, typename MESH::VertexType * &v)
{
  scalar error = mdist;
  typedef VertTmark<MESH> MarkerVert;
  MarkerVert t;
  t.SetMesh(&mesh);
  typedef typename VertexDistance<typename MESH::VertexType> PDistFunct;
  v= vcg::GetClosest<GRID,PDistFunct,MarkerFace>(p,mdist,PDistFunct() ,error,bestq,t,gr);
}


//**ITERATORS DEFINITION**//

template <class GRID,class MESH>
class ClosestFaceIterator:public vcg::ClosestIterator<GRID,FaceDistance<typename MESH::FaceType>,typename FaceTmark<MESH> >
{
public:
	typedef typename GRID GridType;
	typedef typename MESH MeshType;
	typedef typename FaceTmark<MESH> MarkerFace;
	typedef typename FaceDistance<typename MESH::FaceType> FDistFunct;
	typedef typename vcg::ClosestIterator<GRID,FaceDistance<typename MESH::FaceType>,typename FaceTmark<MESH> > ClosestBaseType;

	ClosestFaceIterator(GridType &_Si):ClosestBaseType(_Si,FDistFunct()){}

	void SetMesh(MeshType *m)
	{tm.SetMesh(m);}
};

template <class GRID,class MESH>
class ClosestVertexIterator:public vcg::ClosestIterator<GRID,VertexDistance<typename MESH::VertexType>,typename VertTmark<MESH> >
{
public:
	typedef typename GRID GridType;
	typedef typename MESH MeshType;
	typedef typename VertTmark<MESH> MarkerVert;
	typedef typename VertexDistance<typename MESH::VertexType> VDistFunct;
	typedef typename vcg::ClosestIterator<GRID,VertexDistance<typename MESH::VertexType>,typename VertTmark<MESH> > ClosestBaseType;

	ClosestVertexIterator(GridType &_Si):ClosestBaseType(_Si,VDistFunct()){}

	void SetMesh(MeshType *m)
	{tm.SetMesh(m);}
};

template <class GRID,class MESH>
class TriRayIterator:public vcg::RayIterator<GRID,FaceIntersection<typename MESH::FaceType>,typename FaceTmark<MESH> >
{
public:
	typedef typename GRID GridType;
	typedef typename MESH MeshType;
	typedef typename FaceTmark<MESH> MarkerFace;
	typedef typename FaceIntersection<typename MESH::FaceType> FintFunct;
	typedef typename vcg::RayIterator<GRID,FaceIntersection<typename MESH::FaceType>,typename FaceTmark<MESH> > RayBaseType;

	TriRayIterator(GridType &_Si):RayBaseType(_Si,FintFunct()){}

	void SetMesh(MeshType *m)
	{tm.SetMesh(m);}

};

}	 // end namespace trimesh
}	 // end namespace vcg

#endif
