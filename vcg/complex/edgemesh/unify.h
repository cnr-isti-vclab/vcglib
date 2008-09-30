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
Revision 1.8  2005/10/03 16:16:54  spinelli
used new version of grid query

Revision 1.7  2005/09/19 13:10:12  spinelli
fixed bugs

Revision 1.6  2005/09/14 14:34:41  spinelli
used new version of Grid_ptr

Revision 1.5  2005/05/30 09:42:05  spinelli
std::std::vector<PVertex> sostituito con std::vector<PVertex>

Revision 1.4  2005/05/30 09:13:08  ganovelli
error in include

Revision 1.3  2005/05/17 21:19:37  ganovelli
some std::and typename  missing  (CRS4)

Revision 1.2  2005/03/08 14:42:22  ganovelli
added vcg header


****************************************************************************/

#ifndef __VCGLIB__UNIFY
#define __VCGLIB__UNIFY

#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/space/box3.h>
#include <vcg/complex/edgemesh/update/bounding.h>
#include <vector>

namespace vcg
	{
/** \addtogroup edgemesh */
/*@{*/
/** 
    Class for computing unification of the vertices or of the edges
*/
		template <class EdgeMeshType>
		struct Unify{

			template <class MESH_TYPE,class OBJ_TYPE>
			class Tmark
			{
				MESH_TYPE m;
			public:
				Tmark(MESH_TYPE &_m):m(_m){}
				void UnMarkAll(){m.UnMarkAll();}
				bool IsMarked(OBJ_TYPE* obj){return (m.IsMarked(obj->v));}
				void Mark(OBJ_TYPE* obj){m.Mark(obj->v);}
			};


			typedef typename EdgeMeshType::VertexPointer VertexPointer;
			typedef typename EdgeMeshType::EdgePointer   EdgePointer;
			typedef typename EdgeMeshType::ScalarType ScalarType;
			typedef typename EdgeMeshType::CoordType CoordType;

			struct PVertex:EdgeMeshType::VertexType
			{
				typedef typename EdgeMeshType::ScalarType ScalarType;
				VertexPointer  v;		// the two Vertex pointer are ordered!
				EdgePointer    e;		  // the edge where this vertex belong
				int      z;				      // index in [0..2] of the edge of the face
				PVertex(EdgePointer  pe, const int nz ):e(pe),z(nz),v(pe->V(nz)){}
				/*bool Dist(Point3<ScalarType> p,ScalarType & d,Point3<ScalarType>& res)
				{
					res = p;
					ScalarType _d =vcg::Distance(p,v->P());
					if(d > _d)
					{
						d = _d;
						return true;
					}
					return false;
				}*/

				void GetBBox(vcg::Box3<ScalarType> & bb){
					bb.Add(v->P());
				}
				bool IsD(){
					return v->IsD();
				}
			};
			typedef GridStaticPtr<	PVertex > GridType;

			static void Join(PVertex pv0,PVertex & pv1){
				pv1.e->V(pv1.z) = pv0.v;
				pv1.e = NULL;
			}
			
			class BackCompDist {
			public:
				inline bool operator () (const PVertex & obj, const CoordType & pt, ScalarType & mindist, CoordType & result) {
					result = pt;
					ScalarType _d =vcg::Distance(result,obj.v->P());
					if(mindist < _d)
					{
						mindist = _d;
						return true;
					}
					return false;
				}
			};

			static GridType & Grid(){static GridType grid; return grid; }

			static void Vertices(EdgeMeshType & em,   ScalarType   epsilon){
				typename EdgeMeshType::EdgeIterator ei;

				typedef Tmark<EdgeMeshType,PVertex> Marker;
				Marker tm=Marker(em);

				bool lastRound ;
				if(em.vn){
					vcg::edg::UpdateBounding<EdgeMeshType>::Box(em);
					//Grid().SetBBox(em.bbox);

					std::vector<PVertex> pv;
					for(ei = em.edges.begin(); ei != em.edges.end();++ei){
						pv.push_back(PVertex(&*ei,0));
						pv.push_back(PVertex(&*ei,1));
					}					
					Grid().Set(pv.begin(), pv.end() );
					typename std::vector<PVertex>::iterator pvi;
					Point3<ScalarType> p;
					PVertex * closest;
					do{
						lastRound = true;
						for(pvi = pv.begin(); pvi != pv.end(); ++pvi)
							if((*pvi).e)
							{
								float eps = epsilon;
								Point3<ScalarType> vpos =(*pvi).v->P() ;
								(*pvi).v->SetD();
								ScalarType max_dist=em.bbox.Diag();
								ScalarType min_dist = 0;
								p =  (vcg::tri::GetClosestVertex<EdgeMeshType, GridType>( em, Grid(), vpos, 
									max_dist, min_dist))->P();

								//closest  = Grid().GetClosest<BackCompDist,Marker>(vpos, max_dist, BackCompDist(), eps, p,tm);
								//closest  = Grid().GetClosest(vpos,eps,p);
								(*pvi).v->ClearD();


								if(closest){

									if(closest->e)
									{
										Join(*pvi,*closest);
										lastRound = false;
									}
								}

							}
					}while(!lastRound);

				}
			} // end Vertices
		}; // end class

	} // end namespace

#endif
