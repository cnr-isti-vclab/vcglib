
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
			typedef typename EdgeMeshType::VertexPointer VertexPointer;
			typedef typename EdgeMeshType::EdgePointer   EdgePointer;
			typedef typename EdgeMeshType::ScalarType ScalarType;

			struct PVertex
			{
				typedef typename EdgeMeshType::ScalarType ScalarType;
				VertexPointer  v;		// the two Vertex pointer are ordered!
				EdgePointer    e;		  // the edge where this vertex belong
				int      z;				      // index in [0..2] of the edge of the face
				PVertex(EdgePointer  pe, const int nz ):e(pe),z(nz),v(pe->V(nz)){}
				bool Dist(Point3<ScalarType> p,ScalarType & d,Point3<ScalarType>& res)
				{
					res = p;
					ScalarType _d =vcg::Distance(p,v->P());
					if(d > _d)
					{
						d = _d;
						return true;
					}
					return false;
				}

				void GetBBox(vcg::Box3<ScalarType> & bb){
					bb.Add(v->P());
				}
				bool IsD(){
					return v->IsD();
				}
			};
			typedef typename GridStaticPtr<	std::vector<PVertex> > GridType;

			static void Join(PVertex pv0,PVertex & pv1){
				pv1.e->V(pv1.z) = pv0.v;
				pv1.e = NULL;
			}

			static GridType & Grid(){static GridType grid; return grid; }

			static void Vertices(EdgeMeshType & em,   ScalarType   epsilon){
				EdgeMeshType::EdgeIterator ei;
				bool lastRound ;
				if(em.vn){
					vcg::edge::UpdateBounding<EdgeMeshType>::Box(em);
					Grid().SetBBox(em.bbox);

					vector<PVertex> pv;
					for(ei = em.edges.begin(); ei != em.edges.end();++ei){
						pv.push_back(PVertex(&*ei,0));
						pv.push_back(PVertex(&*ei,1));
					}					
					Grid().Set(pv);
					vector<PVertex>::iterator pvi;
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
								closest  = Grid().GetClosest(vpos,eps,p);
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