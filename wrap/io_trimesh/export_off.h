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

/**
@name Save in OFF format
*/
//@{

#ifndef __VCGLIB_EXPORT_OFF
#define __VCGLIB_EXPORT_OFF

#include<wrap/ply/io_mask.h>

#include <stdio.h>

namespace vcg {
	namespace tri {
		namespace io {
			template <class SaveMeshType>
			class ExporterOFF
			{

			public:
				typedef typename SaveMeshType::VertexPointer VertexPointer;
				typedef typename SaveMeshType::ScalarType ScalarType;
				typedef typename SaveMeshType::VertexType VertexType;
				typedef typename SaveMeshType::FaceType FaceType;
				typedef typename SaveMeshType::FacePointer FacePointer;
				typedef typename SaveMeshType::VertexIterator VertexIterator;
				typedef typename SaveMeshType::FaceIterator FaceIterator;

				static bool Save(SaveMeshType &m, const char * filename )
				{
					vcg::face::Pos<MyMesh::FaceType> he;
					vcg::face::Pos<MyMesh::FaceType> hei;
					FILE * fpout = fopen(filename,"w");
					if(fpout==NULL)	return false;



					if( m.HasPerVertexNormal())
						fprintf(fpout,"N");
					if( m.HasPerVertexColor())
						fprintf(fpout,"C");
					if( m.HasPerVertexTexture())
						fprintf(fpout,"ST");
					fprintf(fpout,"OFF\n");
					fprintf(fpout,"%d %d ", m.vn, m.fn);

					MyMesh::FaceIterator fi;
					int count_e = 0;
					int boundary_e = 0;
					bool counted=false;
					for(fi=m.face.begin();fi!=m.face.end();++fi)
						(*fi).ClearS();



					for(fi=m.face.begin();fi!=m.face.end();fi++)
					{
						(*fi).SetS();
						count_e +=3;								//assume that we have to increase the number of edges with three
						for(int j=0; j<3; j++)
						{
							if (fi->IsBorder(j))			//If this edge is a border edge
								boundary_e++;						//  then increase the number of boundary edges
							else if (IsManifold(*fi,j))		//If this edge is manifold
							{
								if((*fi).FFp(j)->IsS()) //If the face on the other side of the edge is already selected
									count_e--;						//  we counted one edge twice
							}
							else											//We have a non-manifold edge
							{
								hei.Set(&(*fi), j , fi->V(j));
								he=hei;
								he.NextF();
								while (he.f!=hei.f)			//	so we have to iterated all faces that are connected to this edge
								{
									if (he.f->IsS())			//  if one of the other faces was already visited than this edge was counted already.
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
					fprintf(fpout,"%d\n", count_e);

					//vertices
					int j;
					std::vector<int> FlagV; 
					VertexPointer  vp;
					VertexIterator vi;
					for(j=0,vi=m.vert.begin();vi!=m.vert.end();++vi)
					{
						vp=&(*vi);
						FlagV.push_back(vp->UberFlags()); // Salva in ogni caso flag del vertice
						if( ! vp->IsD() )
						{	// ***** ASCII *****

							fprintf(fpout,"Vertex: %g %g %g\n" ,vp->P()[0],vp->P()[1],vp->P()[2]);
							if( m.HasPerVertexColor() )
								fprintf(fpout,"%d %d %d %d\n",vp->C()[0],vp->C()[1],vp->C()[2],vp->C()[3] );

							if( m.HasPerVertexNormal())
								fprintf(fpout,"%g %g %g\n", vp->N()[0],vp->N()[1],vp->N()[2]);

							if( m.HasPerVertexTexture())
								fprintf(fpout,"%g %g\n",vp->T().u(),vp->T().v());
						}

						vp->UberFlags()=j; // Trucco! Nascondi nei flags l'indice del vertice non deletato!
						j++;
					}

					FacePointer fp;
					int vv[3];

					int fcnt=0;
					for(j=0,fi=m.face.begin();fi!=m.face.end();++fi)
					{
						fp=&(*fi);
						if( ! fp->IsD() )
						{ fcnt++;


						fprintf(fpout,"3 %d %d %d\n",
							fp->cV(0)->UberFlags(),	fp->cV(1)->UberFlags(), fp->cV(2)->UberFlags() );
						}
					}


					fclose(fpout);
					// Recupera i flag originali
					for(j=0,vi=m.vert.begin();vi!=m.vert.end();++vi)
						(*vi).UberFlags()=FlagV[j++]; 

					return 0;
				}
			}; // end class
		} // end namespace tri
	} // end namespace io
} // end namespace vcg
//@}
#endif