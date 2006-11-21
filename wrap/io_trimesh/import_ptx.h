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
****************************************************************************/

#ifndef __VCGLIB_IMPORT_PTX
#define __VCGLIB_IMPORT_PTX

#include <io.h>
#include <stdio.h>
#include <wrap/callback.h>
#include <vcg/complex/trimesh/allocate.h>
#include <vcg/complex/trimesh/clean.h>

namespace vcg {
	namespace tri {
		namespace io {
			/** 
			This class encapsulate a filter for importing ptx meshes.
			*/
			template <class OpenMeshType>
			class ImporterPTX
			{
			public:
				enum PTX_OPEN_MASK_ENUM
				{
					PTX_ONLY_POINTS   		= 0x08000000,  //BIT_27 no add faces (PTX_FLIPFACES and PTX_SWITCHSIDE are ignored!)
					PTX_COLOR         		= 0x10000000,  //BIT_28 must be VertexType::HasColor();
					PTX_COMPUTE_AABBOX 		= 0x20000000,  //BIT_29 compute axis aligned bbox
					PTX_FLIPFACES			= 0x40000000,  //BIT_30 flip all faces ( PTX_ONLY_POINTS must be false )
					PTX_SWITCHSIDE    		= 0x80000000   //BIT_31 inverse triangulation order (swaping row->cols) ( PTX_ONLY_POINTS must be false )
				};
				typedef typename OpenMeshType::VertexPointer VertexPointer;
				typedef typename OpenMeshType::ScalarType ScalarType;
				typedef typename OpenMeshType::VertexType VertexType;
				typedef typename OpenMeshType::FaceType FaceType;
				typedef typename OpenMeshType::VertexIterator VertexIterator;
				typedef typename OpenMeshType::FaceIterator FaceIterator;
				
				struct RANGEMAP_INFO
				{
					fpos_t pos;
					int    vn;
					int    fn;
				};

				typedef typename std::vector< RANGEMAP_INFO > RANGEMAP_INFO_TABLE;

				struct PTX_HEAD_INFO
				{
					int    vn;
					int    fn;
					RANGEMAP_INFO_TABLE rmapInfo; 
				};
				
				/// Standard call for knowing the meaning of an error code
				static const char *ErrorMsg(int error)
				{
					static const char * ptx_error_msg[] =
					{
						"No errors",
						"Can't open file",
						"Header not found",
						"Eof in header",
						"Format not found",
						"Syntax error on header",
					};
					if(error>6 || error<0) return "Unknown error";
					else return ptx_error_msg[error];
				};

				static bool skipmesh(FILE* fp, CallBackPos *cb=NULL)
				{
					PTX_HEAD_INFO tab;
					return skipmesh(fp, cb);
				}
				static bool skipmesh(FILE* fp, PTX_HEAD_INFO & tab, CallBackPos *cb=NULL)
				{
					int colnum;
					int rownum;
					int skiplines;
					char linebuf;

					if(feof(fp))	return false;
					RANGEMAP_INFO ptxInfo;
					fgetpos(fp, &ptxInfo.pos );

					// getting mesh size;
					fscanf(fp,"%i\n",&colnum);
					fscanf(fp,"%i\n",&rownum);

					ptxInfo.vn = rownum*colnum;
					ptxInfo.fn = (rownum-1) * (colnum-1) * 2;

					char tmp[255];
					sprintf(tmp, "PTX Mesh analysis... mesh %i  vert %i face %i", (int)tab.rmapInfo.size(), ptxInfo.vn, ptxInfo.fn);

					if ( ( colnum <=0 ) || ( rownum <=0 ) ) return false;

					if(feof(fp))	return false;
					if(cb) cb( rand()%100, tmp);	
					skiplines = (colnum * rownum) + 8; // have to skip (col * row) lines plus 8 lines for the header
					for(int ii=0; ii<skiplines; ii++)
					{
						fread(&linebuf,1,1,fp);
						while(linebuf != '\n')  fread(&linebuf,1,1,fp);
					} 

					if(cb) cb( 100, tmp);
					tab.vn += ptxInfo.vn;
					tab.fn += ptxInfo.fn;
					tab.rmapInfo.push_back( ptxInfo );
					return true;
				}


				static bool Analysis(const char * filename, PTX_HEAD_INFO &info, CallBackPos *cb=NULL)
				{
					info.fn = 0;
					info.vn = 0;
					info.rmapInfo.clear();
					FILE *fp;
					fp = fopen(filename, "rb");
					if(fp == NULL) return false;
					while ( skipmesh( fp, info, cb ) ) {}
					return true;
				};
				static int Open( OpenMeshType &m, const char * filename, int meshNumber, int mask = PTX_ONLY_POINTS, CallBackPos *cb=NULL)
				{
					FILE *fp;
					fp = fopen(filename, "rb");
					if(fp == NULL) return false;
					m.Clear();
					m.vn=0;
					m.fn=0;

					PTX_HEAD_INFO ptxHead;

					if ( meshNumber>0 ) for (int i=0; i!=meshNumber; ++i)  skipmesh(fp, ptxHead, cb);
					if (!readPTX( m, fp, mask, meshNumber, cb))
					{
						m.Clear();
						return 1;
					}
					clearBadVertex(m, mask, cb);
					return 0;
				}

				static void clearBadVertex(OpenMeshType &m, int mask, CallBackPos *cb=NULL)
				{
					if(cb) cb(40,"PTX Mesh Loading - remove bad vertex!");	
					for(VertexIterator vi = m.vert.begin(); vi != m.vert.end(); vi++)
					{
						if((*vi).P() == Point3f(0.0, 0.0, 0.0))
						{
							(*vi).SetD();
							m.vn--;
						}
					}

					if(cb) cb(60,"PTX Mesh Loading - remove bad face!");	
					bool onlypoints  =  ((mask & PTX_ONLY_POINTS) != 0);
					if(! onlypoints)
					{
						for(OpenMeshType::FaceIterator fi = m.face.begin(); fi != m.face.end(); fi++)
						{
							if( ((*fi).V(0)->IsD()) || ((*fi).V(1)->IsD()) || ((*fi).V(2)->IsD()) )
							{
								(*fi).SetD();
								m.fn--;
							}
						}

						// eliminate high angle triangles
						int angle = 88;

						//printf(" culling by angle \n");
						float limit = cos( angle*3.14159265358979323846/180.0 );
						Point3f raggio;

						if(cb) cb(85,"PTX Mesh Loading - remove bad face!");	
						vcg::tri::UpdateNormals<OpenMeshType>::PerFaceNormalized(m);
						for(OpenMeshType::FaceIterator fi = m.face.begin(); fi != m.face.end(); fi++)
							if(!(*fi).IsD())
							{
								raggio = -((*fi).V(0)->P() + (*fi).V(1)->P() + (*fi).V(2)->P()) / 3.0;
								raggio.Normalize();
								if((raggio * (*fi).N()) < limit)
								{
									(*fi).SetD();
									m.fn--;
								}
							}

					}
					/*if(cb) cb(60,"PTX Mesh Loading RemoveDuplicateVertex");	
					tri::Clean<OpenMeshType>::RemoveDuplicateVertex(m);

					if (!onlypoints) 
					{ 
					if(cb) cb(60,"PTX Mesh Loading RemoveUnreferencedVertex");	
					tri::Clean<OpenMeshType>::RemoveUnreferencedVertex(m);
					}*/
					if(cb) cb(100,"PTX Mesh Loading finish!");	


				}

				//if numMesh == -1 load all mesh
				static int Open_( OpenMeshType &m, const char * filename,  int mask = PTX_ONLY_POINTS, CallBackPos *cb=NULL)
				{
					FILE *fp;
					fp = fopen(filename, "rb");
					if(fp == NULL) return 1;

					m.Clear();
					m.vn=0;
					m.fn=0;
					int vn=0;
					int fn=0;
					//PTX_HEAD_INFO tab;
					//tab.clear();
					//while ( skipmesh( fp, tab, cb ) ) {}
					/*if ( (vn<=0) && (fn<=0) ) return false;
					//VertexIterator vi = Allocator<OpenMeshType>::AddVertices(m,vn); 
					//OpenMeshType::FaceIterator fi= Allocator<OpenMeshType>::AddFaces(m,fn);

					VertexIterator vi = Allocator<OpenMeshType>::AddVertices(m, tab[20].vn); 
					FaceIterator fi   = Allocator<OpenMeshType>::AddFaces(m, tab[20].fn);
					readPTX( m, fp, vi, fi, tab[20], mask, 20, cb);
					fclose(fp);   
					/* return true;

					if ( numMesh>0 )
					for (int i=0; i!=numMesh; ++i)  if (!skipmesh(fp, vn, fn, tab)) return false;

					int mn=0;
					if ( numMesh == -1 )
					{
					bool next = true;
					while ( next )
					{
					bool r = readPTX(m, fp, mask, mn, cb); 
					mn++;
					if ((r==false) && (m.vn==0) ) { fclose(fp); return false; }
					else if (r==false) next = false;
					}
					} else 
					{   
					bool r = readPTX(m, fp, mask, numMesh, cb); 
					if ((r==false) && (m.vn==0) ) { fclose(fp); return false; }
					}

					fclose(fp);
					*/
					// now i delete all points in (0,0,0) that are unsampled points
					for(VertexIterator vi = m.vert.begin(); vi != m.vert.end(); vi++)
					{
						if((*vi).P() == Point3f(0.0, 0.0, 0.0))
						{
							(*vi).SetD();
							m.vn--;
						}
					}

					bool onlypoints  =  ((mask & PTX_ONLY_POINTS) != 0);
					if(! onlypoints)
					{


						for(OpenMeshType::FaceIterator fi = m.face.begin(); fi != m.face.end(); fi++)
						{
							if( ((*fi).V(0)->IsD()) || ((*fi).V(1)->IsD()) || ((*fi).V(2)->IsD()) )
							{
								(*fi).SetD();
								m.fn--;
							}
						}

						// eliminate high angle triangles
						int angle = 88;

						printf(" culling by angle \n");
						float limit = cos( angle*3.14159265358979323846/180.0 );
						Point3f raggio;

						vcg::tri::UpdateNormals<OpenMeshType>::PerFaceNormalized(m);
						for(OpenMeshType::FaceIterator fi = m.face.begin(); fi != m.face.end(); fi++)
							if(!(*fi).IsD())
							{
								raggio = -((*fi).V(0)->P() + (*fi).V(1)->P() + (*fi).V(2)->P()) / 3.0;
								raggio.Normalize();
								if((raggio * (*fi).N()) < limit)
								{
									(*fi).SetD();
									m.fn--;
								}
							}

					}
					/*if(cb) cb(60,"PTX Mesh Loading RemoveDuplicateVertex");	
					tri::Clean<OpenMeshType>::RemoveDuplicateVertex(m);

					if (!onlypoints) 
					{ 
					if(cb) cb(60,"PTX Mesh Loading RemoveUnreferencedVertex");	
					tri::Clean<OpenMeshType>::RemoveUnreferencedVertex(m);
					}*/
					if(cb) cb(100,"PTX Mesh Loading finish!");	
					return 0;
				}
				static bool readPTX( OpenMeshType &m, FILE *fp, VertexIterator &vi, FaceIterator &fi, const RANGEMAP_INFO &ptxInfo, int mask, int mn, CallBackPos *cb=NULL)
				{
					int colnum;
					int rownum;
					int numtokens;
					char linebuf[256];
					int ii;
					float xx,yy,zz;	 // position
					float rr,gg,bb;	 // color
					float rf;		     // reflectance
					Matrix44f		currtrasf;

					bool hascolor;

					bool savecolor   =  ((mask & PTX_COLOR) != 0) &&  VertexType::HasColor();
					bool computeBbox =  ((mask & PTX_COMPUTE_AABBOX) != 0);
					bool onlypoints  =  ((mask & PTX_ONLY_POINTS) != 0);
					bool switchside  =  ((mask & PTX_SWITCHSIDE) != 0);
					bool flipfaces   =  ((mask & PTX_FLIPFACES) != 0);
					int total = 50;
					if ( onlypoints ) total = 100;

					if (fsetpos(fp, &ptxInfo.pos)!=0) return false;


					// getting mesh size;
					fscanf(fp,"%i\n",&colnum);
					fscanf(fp,"%i\n",&rownum);

					if ( ( colnum <=0 ) || ( rownum <=0 ) ) return false;

					// initial 4 lines [still don't know what is this :) :)]
					if ( !fscanf(fp,"%f %f %f\n", &xx, &yy, &zz) ) return false;
					if ( !fscanf(fp,"%f %f %f\n", &xx, &yy, &zz) ) return false;
					if ( !fscanf(fp,"%f %f %f\n", &xx, &yy, &zz) ) return false;
					if ( !fscanf(fp,"%f %f %f\n", &xx, &yy, &zz) ) return false;

					// now the transformation matrix
					if ( !fscanf(fp,"%f %f %f %f\n", &(currtrasf.ElementAt(0,0)), &(currtrasf.ElementAt(0,1)), &(currtrasf.ElementAt(0,2)), &(currtrasf.ElementAt(0,3))) )return false;
					if ( !fscanf(fp,"%f %f %f %f\n", &(currtrasf.ElementAt(1,0)), &(currtrasf.ElementAt(1,1)), &(currtrasf.ElementAt(1,2)), &(currtrasf.ElementAt(1,3))) )return false;
					if ( !fscanf(fp,"%f %f %f %f\n", &(currtrasf.ElementAt(2,0)), &(currtrasf.ElementAt(2,1)), &(currtrasf.ElementAt(2,2)), &(currtrasf.ElementAt(2,3))) )return false;
					if ( !fscanf(fp,"%f %f %f %f\n", &(currtrasf.ElementAt(3,0)), &(currtrasf.ElementAt(3,1)), &(currtrasf.ElementAt(3,2)), &(currtrasf.ElementAt(3,3))) )return false;

					// now the real data begins

					// first line, we should know if the format is
					// XX YY ZZ RF
					// or it is
					// XX YY ZZ RF RR GG BB

					// read the entire first line and then count the spaces. it's rude but it works :)
					ii=0;
					fread(&(linebuf[ii++]),1,1,fp);
					while(linebuf[ii-1] != '\n')  if ( fread(&(linebuf[ii++]),1,1,fp)==0 ) return false;
					linebuf[ii-1] = '\0'; // terminate the string
					numtokens=1;
					for(ii=0; ii<(int)strlen(linebuf); ii++) if(linebuf[ii] == ' ') numtokens++;
					if(numtokens == 4)  hascolor = false;
					else if(numtokens == 7)  hascolor = true;
					else  return false;

					Transpose(currtrasf);
					int vn = rownum*colnum;

					//VertexIterator vi = Allocator<OpenMeshType>::AddVertices(m,vn); 
					//m.vn += vn;
					// parse the first line....
					if(hascolor)
					{
						printf("\n hascolor ");
						sscanf(linebuf,"%f %f %f %f %f %f %f", &xx, &yy, &zz, &rf, &rr, &gg, &bb);
					}
					else
					{
						printf("\n no color ");
						sscanf(linebuf,"%f %f %f %f", &xx, &yy, &zz, &rf);
					}

					//if (computeBbox) m.bbox.SetNull();


					//addthefirstpoint
					(*vi).P()[0]=xx;
					(*vi).P()[1]=yy;
					(*vi).P()[2]=zz;
					(*vi).P() = currtrasf * (*vi).P();
					if (computeBbox) m.bbox.Add( (*vi).P() );
					if(hascolor && savecolor)
					{
						(*vi).C()[0]=rr;
						(*vi).C()[1]=gg;
						(*vi).C()[2]=bb;
					}
					vi++;



					// now for each line until end of mesh (row*col)-1
					for(ii=0; ii<((rownum*colnum)-1); ii++)
					{

						char tmp[255];
						sprintf(tmp, "PTX Mesh Loading... mesh %i", mn);
						if(cb) cb((ii*total)/vn, tmp);	

						// read the stream
						if(hascolor)   fscanf(fp,"%f %f %f %f %f %f %f", &xx, &yy, &zz, &rf, &rr, &gg, &bb);
						else  fscanf(fp,"%f %f %f %f", &xx, &yy, &zz, &rf);


						// add the point
						(*vi).P()[0]=xx;
						(*vi).P()[1]=yy;
						(*vi).P()[2]=zz;
						(*vi).P() = currtrasf * (*vi).P();
						if (computeBbox) m.bbox.Add( (*vi).P() );
						if(hascolor && savecolor)
						{
							(*vi).C()[0]=rr;
							(*vi).C()[1]=gg;
							(*vi).C()[2]=bb;
						}
						vi++;



					}

					if(! onlypoints)
					{

						// now i can triangulate
						int trinum = (rownum-1) * (colnum-1) * 2;


						//OpenMeshType::FaceIterator fi= Allocator<OpenMeshType>::AddFaces(m,trinum);

						// m.fn += trinum;


						int v0i,v1i,v2i, t;
						for(int rit=0; rit<rownum-1; rit++)
							for(int cit=0; cit<colnum-1; cit++)
							{

								if(cb) cb(50 + (t*50)/(rownum*colnum),"PTX Mesh Loading");	

								if(!switchside)
								{
									v0i = (rit  ) + ((cit  ) * rownum);
									v1i = (rit+1) + ((cit  ) * rownum);
									v2i = (rit  ) + ((cit+1) * rownum);
								}
								else
								{
									v0i = (cit  ) + ((rit  ) * colnum);
									v1i = (cit+1) + ((rit  ) * colnum);
									v2i = (cit  ) + ((rit+1) * colnum);
								}


								// upper tri
								(*fi).V(2) = &(m.vert[v0i]);
								(*fi).V(1) = &(m.vert[v1i]);
								(*fi).V(0) = &(m.vert[v2i]);

								if(flipfaces)
								{
									(*fi).V(2) = &(m.vert[v1i]);
									(*fi).V(1) = &(m.vert[v0i]);
								}

								//m.fn++;
								fi++;

								if(!switchside)
								{
									v0i = (rit+1) + ((cit  ) * rownum);
									v1i = (rit+1) + ((cit+1) * rownum);
									v2i = (rit  ) + ((cit+1) * rownum);
								}
								else
								{
									v0i = (cit+1) + ((rit  ) * colnum);
									v1i = (cit+1) + ((rit+1) * colnum);
									v2i = (cit  ) + ((rit+1) * colnum);
								}

								// lower tri
								(*fi).V(2) = &(m.vert[v0i]);
								(*fi).V(1) = &(m.vert[v1i]);
								(*fi).V(0) = &(m.vert[v2i]);

								if(flipfaces)
								{
									(*fi).V(2) = &(m.vert[v1i]);
									(*fi).V(1) = &(m.vert[v0i]);
								}

								// m.fn++;
								fi++;
							}



					}

					return true;




				}
				///Call that load a mesh
				static bool readPTX( OpenMeshType &m, FILE *fp, int mask, CallBackPos *cb=NULL)
				{
					int numtokens;
					int colnum;
					int rownum;
					float xx,yy,zz;	 // position
					float rr,gg,bb;	 // color
					float rf;		 // reflectance
					Matrix44f currtrasf;

					bool hascolor;
					bool savecolor   =  ((mask & PTX_COLOR) != 0) &&  VertexType::HasColor();
					bool computeBbox =  ((mask & PTX_COMPUTE_AABBOX) != 0);
					bool onlypoints  =  ((mask & PTX_ONLY_POINTS) != 0);
					bool switchside  =  ((mask & PTX_SWITCHSIDE) != 0);
					bool flipfaces   =  ((mask & PTX_FLIPFACES) != 0);
					int total = 50;
					if ( onlypoints ) total = 100;  
					char linebuf[256];
					
					fscanf(fp,"%i\n",&colnum);					
					fscanf(fp,"%i\n",&rownum);
					
					if ( ( colnum <=0 ) || ( rownum <=0 ) ) return false;
					// initial 4 lines [still don't know what is this :) :)]
     					if ( !fscanf(fp,"%f %f %f\n", &xx, &yy, &zz) ) return false;
					if ( !fscanf(fp,"%f %f %f\n", &xx, &yy, &zz) ) return false;
					if ( !fscanf(fp,"%f %f %f\n", &xx, &yy, &zz) ) return false;
					if ( !fscanf(fp,"%f %f %f\n", &xx, &yy, &zz) ) return false;
											// now the transformation matrix
					if ( !fscanf(fp,"%f %f %f %f\n", &(currtrasf.ElementAt(0,0)), &(currtrasf.ElementAt(0,1)), &(currtrasf.ElementAt(0,2)), &(currtrasf.ElementAt(0,3))) )return false;
					if ( !fscanf(fp,"%f %f %f %f\n", &(currtrasf.ElementAt(1,0)), &(currtrasf.ElementAt(1,1)), &(currtrasf.ElementAt(1,2)), &(currtrasf.ElementAt(1,3))) )return false;
					if ( !fscanf(fp,"%f %f %f %f\n", &(currtrasf.ElementAt(2,0)), &(currtrasf.ElementAt(2,1)), &(currtrasf.ElementAt(2,2)), &(currtrasf.ElementAt(2,3))) )return false;
					if ( !fscanf(fp,"%f %f %f %f\n", &(currtrasf.ElementAt(3,0)), &(currtrasf.ElementAt(3,1)), &(currtrasf.ElementAt(3,2)), &(currtrasf.ElementAt(3,3))) )return false;
   				    	
					//now the real data begins
					// first line, we should know if the format is
					// XX YY ZZ RF
					// or it is
					// XX YY ZZ RF RR GG BB
					// read the entire first line and then count the spaces. it's rude but it works :)
					int ii=0;
					fread(&(linebuf[ii++]),1,1,fp);
					while(linebuf[ii-1] != '\n')  if ( fread(&(linebuf[ii++]),1,1,fp)==0 ) return false;
					linebuf[ii-1] = '\0'; // terminate the string
					numtokens=1;
					for(ii=0; ii<(int)strlen(linebuf); ii++) if(linebuf[ii] == ' ') numtokens++;
					if(numtokens == 4)  hascolor = false;
					else if(numtokens == 7)  hascolor = true;
					else  return false;
					Transpose(currtrasf);
					int vn = rownum*colnum;
					VertexIterator vi = Allocator<OpenMeshType>::AddVertices(m,vn); 
					m.vn = vn;
					// parse the first line....
					if(hascolor)
					{
						printf("\n hascolor ");
						sscanf(linebuf,"%f %f %f %f %f %f %f", &xx, &yy, &zz, &rf, &rr, &gg, &bb);
					}
					else
					{
						printf("\n no color ");
						sscanf(linebuf,"%f %f %f %f", &xx, &yy, &zz, &rf);
					}
					if (computeBbox) m.bbox.SetNull();
					
					//addthefirstpoint
					(*vi).P()[0]=xx;
					(*vi).P()[1]=yy;
					(*vi).P()[2]=zz;
					
					if (computeBbox) m.bbox.Add( (*vi).P() );
					if(hascolor && savecolor)
					{
						(*vi).C()[0]=rr;
						(*vi).C()[1]=gg;
						(*vi).C()[2]=bb;
					}
					vi++;
					// now for each line until end of mesh (row*col)-1
					for(ii=0; ii<((rownum*colnum)-1); ii++)
					{
						char tmp[255];
						sprintf(tmp, "PTX Mesh Loading...");
						if(cb) cb((ii*total)/vn, tmp);	
						// read the stream
						if(hascolor)   fscanf(fp,"%f %f %f %f %f %f %f", &xx, &yy, &zz, &rf, &rr, &gg, &bb);
						else  fscanf(fp,"%f %f %f %f", &xx, &yy, &zz, &rf);
						// add the point
						(*vi).P()[0]=xx;
						(*vi).P()[1]=yy;
						(*vi).P()[2]=zz;
							
						if (computeBbox) m.bbox.Add( (*vi).P() );
						if(hascolor && savecolor)
						{
							(*vi).C()[0]=rr;
							(*vi).C()[1]=gg;
							(*vi).C()[2]=bb;
						}
						vi++;
					}
					if(! onlypoints)
					{
						// now i can triangulate
						int trinum = (rownum-1) * (colnum-1) * 2;
						OpenMeshType::FaceIterator fi= Allocator<OpenMeshType>::AddFaces(m,trinum);
						m.fn = trinum;
						int v0i,v1i,v2i, t;
						t=0;
						for(int rit=0; rit<rownum-1; rit++)
							for(int cit=0; cit<colnum-1; cit++)
							{
								t++;
								if(cb) cb(50 + (t*50)/(rownum*colnum),"PTX Mesh Loading");	

								if(!switchside)
								{
									v0i = (rit  ) + ((cit  ) * rownum);
									v1i = (rit+1) + ((cit  ) * rownum);
									v2i = (rit  ) + ((cit+1) * rownum);
								}
								else
								{
									v0i = (cit  ) + ((rit  ) * colnum);
									v1i = (cit+1) + ((rit  ) * colnum);
									v2i = (cit  ) + ((rit+1) * colnum);
								}
								
								// upper tri
								(*fi).V(2) = &(m.vert[v0i]);
								(*fi).V(1) = &(m.vert[v1i]);
								(*fi).V(0) = &(m.vert[v2i]);
								if(flipfaces)
								{
									(*fi).V(2) = &(m.vert[v1i]);
									(*fi).V(1) = &(m.vert[v0i]);
								}

								fi++;
								if(!switchside)
								{
									v0i = (rit+1) + ((cit  ) * rownum);
									v1i = (rit+1) + ((cit+1) * rownum);
									v2i = (rit  ) + ((cit+1) * rownum);
								}
								else
								{
									v0i = (cit+1) + ((rit  ) * colnum);
									v1i = (cit+1) + ((rit+1) * colnum);
									v2i = (cit  ) + ((rit+1) * colnum);
								}
								// lower tri
								(*fi).V(2) = &(m.vert[v0i]);
								(*fi).V(1) = &(m.vert[v1i]);
								(*fi).V(0) = &(m.vert[v2i]);
								if(flipfaces)
								{
									(*fi).V(2) = &(m.vert[v1i]);
									(*fi).V(1) = &(m.vert[v0i]);
								}

								fi++;
							}
					}	
					if(cb) cb(40,"PTX Mesh Loading - remove bad vertex!");	
					for(OpenMeshType::VertexIterator vi = m.vert.begin(); vi != m.vert.end(); vi++)
					{
						if((*vi).P() == Point3f(0.0, 0.0, 0.0))
						{
							(*vi).SetD();
							m.vn--;
						}
					}
					if(cb) cb(60,"PTX Mesh Loading - remove bad face!");
					onlypoints  =  ((mask & PTX_ONLY_POINTS) != 0);
					if(! onlypoints)
					{
						for(OpenMeshType::FaceIterator fi = m.face.begin(); fi != m.face.end(); fi++)
						{
							if( ((*fi).V(0)->IsD()) || ((*fi).V(1)->IsD()) || ((*fi).V(2)->IsD()) )
							{
								(*fi).SetD();
								m.fn--;
							}
						}
						// eliminate high angle triangles
						int angle = 88;
						printf(" culling by angle \n");
						float limit = cos( angle*3.14159265358979323846/180.0 );
						Point3f raggio;

						vcg::tri::UpdateNormals<OpenMeshType>::PerFaceNormalized(m);
						for(OpenMeshType::FaceIterator fi = m.face.begin(); fi != m.face.end(); fi++)
							if(!(*fi).IsD())
							{
								raggio = -((*fi).V(0)->P() + (*fi).V(1)->P() + (*fi).V(2)->P()) / 3.0;
								raggio.Normalize();
								if((raggio * (*fi).N()) < limit)
								{
									(*fi).SetD();
									m.fn--;
								}
							}
					}
					for(OpenMeshType::VertexIterator vi = m.vert.begin(); vi != m.vert.end(); vi++)
					{
						if(!(*vi).IsD())
							(*vi).P() = currtrasf * (*vi).P();
					}
					vcg::tri::UpdateNormals<OpenMeshType>::PerFaceNormalized(m);
					vcg::tri::UpdateBounding<CMeshO>::Box(m);
					if(cb) cb(100,"PTX Mesh Loading finish!");
					return true;
				}

				///Standard call that reading a mesh
				static int Open( OpenMeshType &m, const char * filename, int mask = PTX_ONLY_POINTS, CallBackPos *cb=NULL)
				{
					FILE *fp;
					fp = fopen(filename, "rb");
					if(fp == NULL) return 1;
					m.Clear();
					m.vn=0;
					m.fn=0;
					if (!readPTX( m, fp, mask,cb))
					{
						m.Clear();
						return 1;
					}
					int endfile,end = 0;
					fscanf(fp,"%i%i",&endfile,&end);
					if(end != 0) 
						return 2;
					return 0;
				}
			}; // end class

		} // end Namespace tri
	} // end Namespace io
} // end Namespace vcg
#endif