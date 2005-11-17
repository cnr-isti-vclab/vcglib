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
Revision 1.11  2005/11/16 16:33:23  rita_borgo
Changed ComputeSelfintersection

Revision 1.10  2005/11/15 12:16:34  rita_borgo
Changed DegeneratedFaces, sets the D flags for each faces
that is found to be degenerated.
CounEdges and ConnectedComponents check now if a face IsD()
else for degenerated faces many asserts fail.

Revision 1.9  2005/11/14 09:28:18  cignoni
changed access to face functions (border, area)
removed some typecast warnings

Revision 1.8  2005/10/11 16:03:40  rita_borgo
Added new functions belonging to triMeshInfo
Started the Self-Intersection routine

Revision 1.7  2005/10/03 15:57:53  rita_borgo
Alligned with TriMeshInfo Code

Revision 1.6  2005/01/28 11:59:35  cignoni
Add std:: to stl containers

Revision 1.5  2004/09/20 08:37:57  cignoni
Better Doxygen docs

Revision 1.4  2004/08/25 15:15:26  ganovelli
minor changes to comply gcc compiler (typename's and stuff)

Revision 1.3  2004/07/18 06:55:37  cignoni
NewUserBit -> NewBitFlag

Revision 1.2  2004/07/09 15:48:37  tarini
Added an include (<algorithm>)

Revision 1.1  2004/06/24 08:03:59  cignoni
Initial Release


****************************************************************************/

#ifndef __VCGLIB_CLEAN
#define __VCGLIB_CLEAN

#include <map>
#include <algorithm>

#include <vcg/simplex/face/face.h>
#include<vcg/simplex/face/topology.h>
#include <vcg/complex/trimesh/base.h>
#include <vcg/complex/trimesh/closest.h>
#include <vcg/space/index/grid_static_ptr.h>

#include<vcg/complex/trimesh/allocate.h>


namespace vcg {
	namespace tri{
		/// 
		/** \addtogroup trimesh */
		/*@{*/
		/// Class of static functions to clean/correct/restore meshs. 
		template <class CleanMeshType>
		class Clean
		{

		public:
			typedef CleanMeshType MeshType; 
			typedef typename MeshType::VertexType     VertexType;
			typedef typename MeshType::VertexPointer  VertexPointer;
			typedef typename MeshType::VertexIterator VertexIterator;
			typedef	typename MeshType::ScalarType			ScalarType;
			typedef typename MeshType::FaceType       FaceType;
			typedef typename MeshType::FacePointer    FacePointer;
			typedef typename MeshType::FaceIterator   FaceIterator;
			typedef typename MeshType::FaceContainer  FaceContainer;

			typedef GridStaticPtr<FaceType, typename MeshType::ScalarType > TriMeshGrid;
			typedef Point3<typename MeshType::ScalarType> Point3x;

			TriMeshGrid   gM;
			FaceIterator fi;
			FaceIterator gi;
			vcg::face::Pos<FaceType> he;
			vcg::face::Pos<FaceType> hei;

			/* classe di confronto per l'algoritmo di eliminazione vertici duplicati*/
			class RemoveDuplicateVert_Compare{
			public:
				inline bool operator()(VertexPointer &a, VertexPointer &b)
				{
					return (*a).cP() < (*b).cP();
				}
			};

			static int DetectUnreferencedVertex( MeshType& m )   // V1.0
			{
				int count_uv = 0;
				MeshType::VertexIterator v;
				FaceIterator fi;

				for(v=m.vert.begin();v!=m.vert.end();++v)
					(*v).ClearV();

				for(fi=m.face.begin();fi!=m.face.end();++fi)
					for(int j=0;j<3;++j)
						(*fi).V(j)->SetV();

				for(v=m.vert.begin();v!=m.vert.end();++v)
					if( !(*v).IsV() )
						++count_uv;
				return count_uv;

			}

			/** This function removes all duplicate vertices of the mesh by looking only at their spatial positions. 
			Note that it does not update any topology relation that could be affected by this like the VT or TT relation.
			the reason this function is usually performed BEFORE building any topology information.
			*/
			static int RemoveDuplicateVertex( MeshType & m )    // V1.0
			{
				if(m.vert.size()==0 || m.vn==0) return 0;

				std::map<VertexPointer, VertexPointer> mp;
				int i,j;
				VertexIterator vi; 
				int deleted=0;
				int k=0;
				size_t num_vert = m.vert.size();
				std::vector<VertexPointer> perm(num_vert);
				for(vi=m.vert.begin(); vi!=m.vert.end(); ++vi, ++k)
					perm[k] = &(*vi);

				RemoveDuplicateVert_Compare c_obj;

				std::sort(perm.begin(),perm.end(),c_obj);

				j = 0;
				i = j;
				mp[perm[i]] = perm[j];
				++i;
				for(;i!=num_vert;)
				{
					if( (! (*perm[i]).IsD()) && 
						(! (*perm[j]).IsD()) && 
						(*perm[i]).P() == (*perm[j]).cP() )
					{
						VertexPointer t = perm[i];
						mp[perm[i]] = perm[j];
						++i;
						(*t).SetD();
						deleted++;
					}
					else
					{
						j = i;
						++i;
					}
				}
				FaceIterator fi;
				for(fi = m.face.begin(); fi!=m.face.end(); ++fi)
					for(k = 0; k < 3; ++k)
						if( !(*fi).IsD() )
							if( mp.find( (typename MeshType::VertexPointer)(*fi).V(k) ) != mp.end() )
							{
								(*fi).V(k) = &*mp[ (*fi).V(k) ];
							}
							m.vn -= deleted;
							return deleted;
			}


			/** This function removes that are not referenced by any face. The function updates the vn counter.
			@param m The mesh
			@return The number of removed vertices
			*/
			static int RemoveUnreferencedVertex( MeshType& m )   // V1.0
			{
				FaceIterator fi;
				VertexIterator vi;
				int referredBit = VertexType::NewBitFlag();

				int j;
				int deleted = 0;

				for(vi=m.vert.begin();vi!=m.vert.end();++vi)
					(*vi).ClearUserBit(referredBit);

				for(fi=m.face.begin();fi!=m.face.end();++fi)
					if( !(*fi).IsD() )
						for(j=0;j<3;++j)
							(*fi).V(j)->SetUserBit(referredBit);

				for(vi=m.vert.begin();vi!=m.vert.end();++vi)
					if( (!(*vi).IsD()) && (!(*vi).IsUserBit(referredBit)))
					{
						(*vi).SetD();
						++deleted;
					}
					m.vn -= deleted;
					VertexType::DeleteBitFlag(referredBit);
					return deleted;
			}


			static bool IsComplexManifold( MeshType & m ) 
			{
				FaceIterator fi;
				bool Manifold = true;
				for( fi=m.face.begin();fi!=m.face.end();++fi)
				{
					for (int j=0;j<3;j++)
					{
						if(!IsManifold(*fi,j))
						{
							Manifold = false;
							fi= m.face.end();
							--fi;
							j=3;
						}
					}
					if((BorderCount(*fi)>0))
					{
						Manifold = false;
						fi= m.face.end();
						--fi;
					}
				}
				return Manifold;
			}

			static void CountEdges( MeshType & m, int &count_e, int &boundary_e ) 
			{
				FaceIterator fi;
				vcg::face::Pos<FaceType> he;
				vcg::face::Pos<FaceType> hei;
				bool counted =false;
				for(fi=m.face.begin();fi!=m.face.end();fi++)
				{
					if(!((*fi).IsD()))
					{
					(*fi).SetS();
					count_e +=3; //assume that we have to increase the number of edges with three
					for(int j=0; j<3; j++)
					{
            if (face::IsBorder(*fi,j)) //If this edge is a border edge
							boundary_e++; // then increase the number of boundary edges
						else if (IsManifold(*fi,j))//If this edge is manifold
						{
							if((*fi).FFp(j)->IsS()) //If the face on the other side of the edge is already selected
								count_e--; // we counted one edge twice
						}
						else//We have a non-manifold edge
						{
							hei.Set(&(*fi), j , fi->V(j));
							he=hei;
							he.NextF();
							while (he.f!=hei.f)// so we have to iterate all faces that are connected to this edge
							{
								if (he.f->IsS())// if one of the other faces was already visited than this edge was counted already.
								{
									counted=true;
									break;
								}
								else
								{
									he.NextF();
								}
							}
							if (counted)
							{
								count_e--;
								counted=false;
							}
						}
					}
					}
				}
			}


			static int CountHoles( MeshType & m)
			{
				int numholes=0;
				int numholev=0;
				int BEdges=0;
				FaceIterator fi;
				FaceIterator gi;
				vcg::face::Pos<FaceType> he;
				vcg::face::Pos<FaceType> hei;

				vector<vector<Point3x> > holes; //indices of vertices

				for(fi=m.face.begin();fi!=m.face.end();++fi)
					(*fi).ClearS();
				gi=m.face.begin(); fi=gi;



				for(fi=m.face.begin();fi!=m.face.end();fi++)//for all faces do
				{
					for(int j=0;j<3;j++)//for all edges
					{
						if(fi->V(j)->IsS()) continue;

            if(face::IsBorder(*fi,j))//found an unvisited border edge
						{
							he.Set(&(*fi),j,fi->V(j)); //set the face-face iterator to the current face, edge and vertex
							vector<Point3x> hole; //start of a new hole
							hole.push_back(fi->P(j)); // including the first vertex
							numholev++;
							he.v->SetS(); //set the current vertex as selected
							he.NextB(); //go to the next boundary edge


							while(fi->V(j) != he.v)//will we do not encounter the first boundary edge.
							{
								Point3x newpoint = he.v->P(); //select its vertex.
								if(he.v->IsS())//check if this vertex was selected already, because then we have an additional hole.
								{
									//cut and paste the additional hole.
									vector<Point3x> hole2;
									int index = find(hole.begin(),hole.end(),newpoint) - hole.begin();
									for(unsigned int i=index; i<hole.size(); i++)
										hole2.push_back(hole[i]);

									hole.resize(index);
									if(hole2.size()!=0) //annoying in degenerate cases
										holes.push_back(hole2);
								}
								hole.push_back(newpoint);
								numholev++;
								he.v->SetS(); //set the current vertex as selected
								he.NextB(); //go to the next boundary edge
							}
							holes.push_back(hole);
						}
					}
				}
				return holes.size();
			}

			static int BorderEdges( MeshType & m, int numholes)
			{
				int BEdges = 0;
				for(int i=0; i<numholes; i++)
				{
					if(i==numholes-1){ printf("%i)\n",numholes); BEdges++;}
					else{ printf("%i, ",numholes); BEdges++;}
				}
				return BEdges;

			}

			static int ConnectedComponents(MeshType &m)
			{
				FaceIterator fi;
				FaceIterator gi;
				vcg::face::Pos<FaceType> he;
				vcg::face::Pos<FaceType> hei;

				vector<int> nrfaces;
				nrfaces.reserve(1);

				for(fi=m.face.begin();fi!=m.face.end();++fi)
					(*fi).ClearS();
				gi=m.face.begin(); fi=gi;
				int Compindex=0;
				stack<MeshType::FaceIterator> sf;
				MeshType::FaceType *l;
				for(fi=m.face.begin();fi!=m.face.end();++fi)
				{
					if(!((*fi).IsD()))
					{
					if (!(*fi).IsS())
					{
						(*fi).SetS();
						//(*fi).Q()=Compindex;
						nrfaces.push_back(1);
						sf.push(fi);
						while (!sf.empty())
						{
							gi=sf.top();
							he.Set(&(*gi),0,gi->V(0));
							sf.pop();
							for(int j=0;j<3;++j)
							{
								if( !face::IsBorder(*gi,j) )
								{
									l=he.f->FFp(j);
									if( !(*l).IsS() )
									{
										(*l).SetS();
										sf.push(l);
									}
								}
							}
						}
						Compindex++;
					}
					}
				}
				return Compindex;
			}

			static int DegeneratedFaces(MeshType& m)
			{
				FaceIterator fi;
				int count_fd = 0;


				for(fi=m.face.begin(); fi!=m.face.end();++fi)
					if(Area<FaceType>(*fi) == 0)
					{
						count_fd++;
						fi->SetD();
						m.fn--;
					}
				return count_fd;
			}
      /**
      GENUS: A topologically invariant property of a surface defined as:
      the largest number of nonintersecting simple closed curves that can be drawn on the surface without separating it. 
      Roughly speaking, it is the number of holes in a surface. 
      The genus g of a surface, also called the geometric genus, is related to the Euler characteristic $chi$ by $chi==2-2g$.
      
      The genus of a connected, orientable surface is an integer representing the maximum number of cuttings along closed 
      simple curves without rendering the resultant manifold disconnected. It is equal to the number of handles on it.

    */  
			static int MeshGenus(MeshType &m, int count_uv, int numholes, int numcomponents, int count_e)
			{
				int eulernumber = (m.vn-count_uv) + m.fn - count_e;
				return int(-( 0.5 * (eulernumber - numholes) - numcomponents ));
			}
/*
Let a closed surface have genus g. Then the polyhedral formula generalizes to the Poincaré formula
chi=V-E+F,	(1)

where
chi(g)==2-2g	(2)

is the Euler characteristic, sometimes also known as the Euler-Poincaré characteristic. 
The polyhedral formula corresponds to the special case g==0.
*/

      static int EulerCharacteristic()
      {

      }


			static void IsRegularMesh(MeshType &m, bool Regular, bool Semiregular)
			{
				int inc=0;
				VertexIterator v;
				FaceIterator fi;
				vcg::face::Pos<FaceType> he;
				vcg::face::Pos<FaceType> hei;
				for(v=m.vert.begin();v!=m.vert.end();++v)
					(*v).ClearS();
				for(fi=m.face.begin();fi!=m.face.end();++fi)
				{
					for (int j=0; j<3; j++)
					{
						he.Set(&(*fi),j,fi->V(j));
						if (!face::IsBorder(*fi,j) && !face::IsBorder(*fi,(j+2)%3) && !fi->V(j)->IsS())
						{
							hei=he;
							inc=1;
							he.FlipE();
							he.NextF();
							while (he.f!=hei.f)
							{
								he.FlipE();
								if (he.IsBorder())
								{
									inc=6;
									break;
								}
								he.NextF();
								inc++;
							}
							if (inc!=6)
								Regular=false;
							if (inc!=6 && inc!=5)
								Semiregular=false;
							fi->V(j)->SetS();

						}
						else
							fi->V(j)->SetS();
					}
					if (Semiregular==false)
						break;

				}

			}

			static void IsOrientedMesh(MeshType &m, bool Oriented, bool Orientable)
			{
				FaceIterator fi;
				FaceIterator gi;
				vcg::face::Pos<FaceType> he;
				vcg::face::Pos<FaceType> hei;
				stack<MeshType::FaceIterator> sf;	
				MeshType::FacePointer l;

				for(fi=m.face.begin();fi!=m.face.end();++fi)
				{
					(*fi).ClearS();
					(*fi).ClearUserBit(0);
				}
				gi=m.face.begin(); fi=gi;
				for(fi=m.face.begin();fi!=m.face.end();++fi)
				{
					if (!(*fi).IsS())
					{
						(*fi).SetS();
						sf.push(fi);

						while (!sf.empty())
						{
							gi=sf.top();
							sf.pop();
							for(int j=0;j<3;++j)
							{
								if( !face::IsBorder(*gi,j) )
								{
									he.Set(&(*gi),0,gi->V(0));
									l=he.f->FFp(j);
									he.Set(&(*gi),j,gi->V(j));
									hei.Set(he.f->FFp(j),he.f->FFi(j), (he.f->FFp(j))->V(he.f->FFi(j)));
									if( !(*gi).IsUserBit(0) )
									{
										if (he.v!=hei.v)// bene
										{
											if ((*l).IsS() && (*l).IsUserBit(0))
											{
												Orientable=false;
												break;
											}
											else if (!(*l).IsS())
											{
												(*l).SetS();
												sf.push(l);
											}
										}
										else if (!(*l).IsS())
										{
											Oriented=false;
											(*l).SetS();
											(*l).SetUserBit(0);
											sf.push(l);
										}
										else if ((*l).IsS() && !(*l).IsUserBit(0))
										{
											Orientable=false;
											break;
										}
									}
									else if (he.v==hei.v)// bene
									{
										if ((*l).IsS() && (*l).IsUserBit(0))
										{
											Orientable=false;
											break;
										}
										else if (!(*l).IsS())
										{
											(*l).SetS();
											sf.push(l);
										}
									}
									else if (!(*l).IsS())
									{
										Oriented=false;
										(*l).SetS();
										(*l).SetUserBit(0);
										sf.push(l);
									}
									else if ((*l).IsS() && !(*l).IsUserBit(0))
									{
										Orientable=false;
										break;
									}
								}
							}
						}
					}
					if (!Orientable)
						break;
				}
			}

			static bool SelfIntersections(MeshType &m)
			{
				
				Box3< ScalarType> bbox;
				TriMeshGrid   gM;

				int nelem;
				
				double bdiag = m.bbox.Diag();

				bbox.SetNull();
				FaceType   *f=0;
				Point3x             normf, bestq, ip,p;
				FaceIterator fi;

				

				std::vector<FaceType*> ret;
				std::vector<FaceType*> inBox;
				gM.Set<vector<FaceType>::iterator>(m.face.begin(),m.face.end());

				for(fi=m.face.begin();fi!=m.face.end();++fi)
				{
					
			//		for(int i =0; i<3; i++)
				//		bbox.Add((*fi).V(i)->P());
					nelem = vcg::trimesh::GetInBoxFace(m, gM, bbox,inBox);
				// fill the cell
/*....*/

					nelem = inBox.size();
					if (nelem>=2)// in a cell
					{
						//test combinations of elements
						for (int i=0;i<nelem-1;i++)
							for (int j=i+1;j<nelem;j++)
							if ((!inBox[i]->IsD())&&(!inBox[j]->IsD())&&(TestIntersection(inBox[i],inBox[j])))
							{
								ret.push_back(inBox[i]);
								ret.push_back(inBox[j]);
							}
					}
				}	
				return false;
			}

				
	//test real intersection between faces
static	bool TestIntersection(FaceType *f0,FaceType *f1)
	{
		assert((!f0->IsD())&&(!f1->IsD()));
		//no adiacent faces
		if ((f0!=f1)&& (!ShareEdge(f0,f1))
			&& (!ShareVertex(f0,f1)))
			return (vcg::Intersection<FaceType>((*f0),(*f1)));
		return false;
	}

			//control if two faces share an edge
static	bool ShareEdge(FaceType *f0,FaceType *f1)
	{
		assert((!f0->IsD())&&(!f1->IsD()));
		for (int i=0;i<3;i++)
				if (f0->FFp(i)==f1)
					return (true);

		return(false);
	}

	//control if two faces share a vertex
static	bool ShareVertex(FaceType *f0,FaceType *f1)
	{
		assert((!f0->IsD())&&(!f1->IsD()));
		for (int i=0;i<3;i++)
			for (int j=0;j<3;j++)
				if (f0->V(i)==f1->V(j))
					return (true);

		return(false);
	}


		}; // end class
		/*@}*/
	
	} //End Namespace Tri
} // End Namespace vcg
#endif
