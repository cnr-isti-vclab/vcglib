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
Revision 1.22  2008/03/13 08:48:10  granzuglia
added two missing include files:
1) #include <wrap/callback.h>
2) #include <wrap/io_trimesh/io_mask.h>

Revision 1.21  2007/12/13 17:57:33  cignoni
removed harmless gcc warnings

Revision 1.20  2007/03/12 16:40:17  tarini
Texture coord name change!  "TCoord" and "Texture" are BAD. "TexCoord" is GOOD.

Revision 1.19  2006/03/29 08:50:10  corsini
Fix bug in texture coordinates reading

Revision 1.18  2006/03/29 08:15:46  corsini
Fix several bugs
Add LoadMask
Improve parsing capabilities (account for unexpected newline)

Revision 1.17  2006/03/01 08:25:30  cignoni
Corrected bug in wrong counting the parsed tokens during the reading of color components

Revision 1.16  2006/02/28 15:18:10  corsini
Fix loading mask update

Revision 1.15  2006/02/09 15:56:34  corsini
Update load mask

Revision 1.14  2006/02/09 15:18:32  corsini
*** empty log message ***

Revision 1.12  2006/02/06 13:11:01  corsini
Renamed UnexpectedEOF as InvalidFile and
added UnsupportedFormat and ErrorNotTriangularFace (by Laurent Saboret)

Revision 1.11  2006/01/30 15:02:50  cignoni
Added mask filling in open

Revision 1.10  2006/01/10 13:20:42  cignoni
Changed ply::PlyMask to io::Mask

Revision 1.9  2005/12/01 00:58:56  cignoni
Added and removed typenames for gcc compiling...

Revision 1.8  2005/11/12 18:12:16  cignoni
Added casts and changed integral types to remove warnings

Revision 1.7  2005/09/28 10:30:14  rita_borgo
*** empty log message ***

Revision 1.6  2005/01/26 22:44:51  cignoni
Resolved scoping of constant of OFF codes

Revision 1.5  2005/01/18 12:35:18  rita_borgo
Added #include<vcg/complex/trimesh/allocate.h>
it was giving problems with Allocator::

Revision 1.4  2005/01/03 11:18:24  cignoni
changed a .. rfind('OFF') .. in rfind("OFF") and added some casts

Revision 1.3  2004/11/23 11:56:50  cignoni
Corrected small bug in the tokenizer (it would add a fake token for lines ending with a space before /n)


****************************************************************************/

#ifndef __VCGLIB_IMPORT_OFF
#define __VCGLIB_IMPORT_OFF

#include <fstream>
#include <string>
#include <vector>
#include <assert.h>
#include <vcg/space/color4.h>
#include<vcg/complex/trimesh/allocate.h>
#include <wrap/callback.h>
#include <wrap/io_trimesh/io_mask.h>

namespace vcg
{
	namespace tri
	{
		namespace io
		{
			// /** \addtogroup  */
			// /* @{ */
			/** 
			This class encapsulate a filter for importing OFF meshes.
			A basic description of the OFF file format can be found at http://www.geomview.org/docs/html/geomview_41.html
			*/
			template<class MESH_TYPE>
			class ImporterOFF
			{
			public:

				typedef typename MESH_TYPE::VertexType			VertexType;
				typedef typename MESH_TYPE::VertexIterator	VertexIterator;
				typedef typename MESH_TYPE::VertexPointer		VertexPointer;
				typedef typename MESH_TYPE::FaceType				FaceType;
				typedef typename MESH_TYPE::FaceIterator		FaceIterator;
				typedef typename MESH_TYPE::FacePointer			FacePointer;
				typedef typename MESH_TYPE::CoordType				CoordType;
				typedef typename MESH_TYPE::ScalarType			ScalarType;

				// OFF codes
				enum OFFCodes {NoError=0, CantOpen, InvalidFile,
					InvalidFile_MissingOFF,
          UnsupportedFormat, ErrorNotTriangularFace,ErrorHighDimension,ErrorDegenerateFace};

				/*!
				*	Standard call for knowing the meaning of an error code
				* \param message_code	The code returned by <CODE>Open</CODE>
				*	\return							The string describing the error code
				*/
				static const char* ErrorMsg(int message_code)
				{
					static const char* error_msg[] =
					{
						"No errors", "Can't open file", "Invalid file",
						"Invalid file: OFF file should have in the first line the OFF keyword as a first token",
              "Unsupported format", "Face with more than 3 vertices","File with high dimensional vertexes are not supported", "Error Degenerate Face with less than 3 vertices"				};

          if(message_code>6 || message_code<0)
						return "Unknown error";
					else
						return error_msg[message_code];
				};

				/**
				 * Load only the properties of the 3D objects.
				 *
				 * \param filename    the name of the file to read from
				 * \param loadmask    the mask which encodes the properties
				 * \return            the operation result
				 */
				static bool LoadMask(const char *filename, int &loadmask)
				{
					// To obtain the loading mask all the file must be parsed
					// to distinguish between per-vertex and per-face color attribute.
          loadmask=0;
					MESH_TYPE dummyMesh;
					return (Open(dummyMesh, filename, loadmask,0,true)==NoError);
				}

				static int Open(MESH_TYPE &mesh, const char *filename,CallBackPos *cb=0)
				{
					int loadmask;
					return Open(mesh,filename,loadmask,cb);
				}

				/*!
				 *  Standard call for reading a mesh.
				 *
				 *  \param mesh         the destination mesh
				 *  \param filename     the name of the file to read from
				 *  \return             the operation result
				 */
				static int Open(MESH_TYPE &mesh, const char *filename, int &loadmask,
					CallBackPos *cb=0, bool onlyMaskFlag=false )
				{
					std::ifstream stream(filename);
					if (stream.fail())
						return CantOpen;

					std::vector< std::string > tokens;
					TokenizeNextLine(stream, tokens);

					bool isNormalDefined   = false;
					bool isColorDefined    = false;
					bool isTexCoordDefined = false;
					int dimension = 3;
					bool homogeneousComponents = false;


					/*
					 [ST][C][N][4][n]OFF	# Header keyword
					 [Ndim]		# Space dimension of vertices, present only if nOFF
					 NVertices  NFaces  NEdges   # NEdges not used or checked
					 
					 x[0]  y[0]  z[0]	# Vertices, possibly with normals, colors, and/or texture coordinates, in that order,  if the prefixes N, C, ST are present.
					 # If 4OFF, each vertex has 4 components including a final homogeneous component.
					 # If nOFF, each vertex has Ndim components.
					 # If 4nOFF, each vertex has Ndim+1 components.
					 ...
					 x[NVertices-1]  y[NVertices-1]  z[NVertices-1]
					 
					 # Faces
					 # Nv = # vertices on this face
					 # v[0] ... v[Nv-1]: vertex indices
					 #		in range 0..NVertices-1
					 Nv  v[0] v[1] ... v[Nv-1]  colorspec
					 ...
					 # colorspec continues past v[Nv-1] to end-of-line; may be 0 to 4 numbers
					 # nothing: default
					 # integer: colormap index
					 # 3 or 4 integers: RGB[A] values 0..255
					 # 3 or 4 floats: RGB[A] values 0..1				 
					 */
					std::string header = tokens[0];
					if (header.rfind("OFF") != std::basic_string<char>::npos)
					{ // the OFF string is in the header go on parsing it.
						for (int u = static_cast<int>(header.rfind("OFF")-1); u>=0; u--)
						{
							if      (header[u] == 'C')				isColorDefined = true;
							else if (header[u] == 'N')				isNormalDefined = true;
							else if (u>0 && header[u-1] == 'S' && header[u] == 'T') isTexCoordDefined = true;
							else if (header[u] == '4')				homogeneousComponents = true;
							else if (header[u] == 'n') return ErrorHighDimension;
							}
						}
				else return InvalidFile_MissingOFF;

        // If the file is slightly malformed and it has nvert and nface AFTER the OFF string instead of in the next line
				// we manage it here...
				if(tokens.size()==1) TokenizeNextLine(stream, tokens);
				else tokens.erase(tokens.begin(),tokens.begin()+1);
					
					// Update loading mask
					///////////////////////////////////////

					loadmask = Mask::IOM_VERTCOORD | Mask::IOM_FACEINDEX;

					if (isNormalDefined)		loadmask |= Mask::IOM_VERTNORMAL;
					if (isTexCoordDefined)	loadmask |= Mask::IOM_VERTTEXCOORD;
					if (isColorDefined)			{ loadmask |= Mask::IOM_VERTCOLOR;loadmask |= Mask::IOM_FACECOLOR;}


          if(onlyMaskFlag) return NoError;
					
					
					mesh.Clear();
					
					// check on next 2 lines to detect corrupted files
					if(tokens.size() < 3)
						return InvalidFile;

					unsigned int nVertices, nFaces, nEdges;
					nVertices = atoi(tokens[0].c_str());
					nFaces    = atoi(tokens[1].c_str());
					nEdges    = atoi(tokens[2].c_str());

					// dimension is the space dimension of vertices => it must be three(!)
					if (dimension != 3)
						return UnsupportedFormat;

					if (homogeneousComponents)
						return UnsupportedFormat;

					// READ VERTICES
					//////////////////////////////////////////////////////

					VertexIterator v_iter = Allocator<MESH_TYPE>::AddVertices(mesh, nVertices);
					TokenizeNextLine(stream, tokens);
					size_t k = 0; // next token to read

					for (unsigned int i=0; i<nVertices; i++, v_iter++)
					{
						if (cb && (i%1000)==0) 
							cb(i*50/nVertices, "Vertex Loading");

						// Read 3 vertex coordinates
						for (unsigned int j=0; j<3; j++)
						{
							// Go to next line when needed
							if (k == tokens.size())   // if EOL
							{
								TokenizeNextLine(stream, tokens);
								if (tokens.size() == 0) // if EOF
									return InvalidFile;
								k = 0;
							}

							// Read vertex coordinate
							(*v_iter).P()[j] = (ScalarType) atof(tokens[k].c_str());
							k++;
						}

						if (isNormalDefined)
						{
							// Read 3 normal coordinates
							for (unsigned int j=0; j<3; j++)
							{
								// Go to next line when needed
								if (k == tokens.size())   // if EOL
								{
									TokenizeNextLine(stream, tokens);
									if (tokens.size() == 0) // if EOF
										return InvalidFile;
									k = 0;
								}

								// Read normal coordinate
								(*v_iter).N()[j] = (ScalarType) atof(tokens[k].c_str());
								k++;
							}
						}

						// NOTE: It is assumed that colored vertex takes exactly one text line
						//       (otherwise it is impossible to parse color information since
						//        color components can vary)
						if (isColorDefined)
						{
							// The number of color components varies from 0 to 4.
							// The OFF format guaranties that there is 1 vertex per line.
							int nb_color_components = static_cast<int>(tokens.size())
								- static_cast<int>(k) /* tokens already parsed */
								- 2 * (isTexCoordDefined ? 1 : 0);

							if (nb_color_components < 0 || nb_color_components > 4)
								return InvalidFile;

							// set per-vertex color attribute
							if (nb_color_components > 0)
								loadmask |= Mask::IOM_VERTCOLOR;

							// Store color components
							if (tri::HasPerVertexColor(mesh))
							{
								// Read color components

								if (nb_color_components == 1)
								{
									// read color index
									(*v_iter).C().Import(ColorMap(atoi(tokens[k].c_str())));
								}
								else if (nb_color_components == 3)
								{
									// read RGB color
									if (tokens[k].find(".") == size_t(-1))// if it is a float there is a dot
									{ 
										// integers
										unsigned char r = 
											static_cast<unsigned char>(atoi(tokens[k].c_str()));
										unsigned char g = 
											static_cast<unsigned char>(atoi(tokens[k+1].c_str()));
										unsigned char b = 
											static_cast<unsigned char>(atoi(tokens[k+2].c_str()));

										vcg::Color4b color(r, g, b, 255);
										(*v_iter).C().Import(color);
									}
									else
									{ 
										// floats
										float r = static_cast<float>(atof(tokens[k].c_str()));
										float g = static_cast<float>(atof(tokens[k+1].c_str()));
										float b = static_cast<float>(atof(tokens[k+2].c_str()));

										vcg::Color4f color(r, g, b, 1.0);
										(*v_iter).C().Import(color);
									}
								}
								else if (nb_color_components == 4)
								{
									// read RGBA color
									if (tokens[k].find(".") == size_t(-1))
									{ 
										// integers
										unsigned char r = 
											static_cast<unsigned char>(atoi(tokens[k].c_str()));
										unsigned char g = 
											static_cast<unsigned char>(atoi(tokens[k+1].c_str()));
										unsigned char b = 
											static_cast<unsigned char>(atoi(tokens[k+2].c_str()));
										unsigned char a =
											static_cast<unsigned char>(atoi(tokens[k+3].c_str()));

										Color4b color(r, g, b, a);
										(*v_iter).C().Import(color);
									}
									else
									{ 
										// floats
										float r = static_cast<float>(atof(tokens[k].c_str()));
										float g = static_cast<float>(atof(tokens[k+1].c_str()));
										float b = static_cast<float>(atof(tokens[k+2].c_str()));
										float a = static_cast<float>(atof(tokens[k+3].c_str()));

										vcg::Color4f color(r, g, b, a);
										(*v_iter).C().Import(color);
									}
								}
							}

							k += nb_color_components;
						}

						if (isTexCoordDefined)
						{
							for (unsigned int j=0; j<2; j++)
							{
								// Go to next line when needed
								if (k == tokens.size())   // if EOL
								{
									TokenizeNextLine(stream, tokens);
									if (tokens.size() == 0) // if EOF
										return InvalidFile;
									k = 0;
								}

								std::string str = tokens[k];
								k++;

								// Store texture coordinates
								if (tri::HasPerWedgeTexCoord(mesh))
								{
									//...TODO...
								}
							}
						}
					} // for i=...

					// READ FACES
					//////////////////////////////////////////////////////

					Allocator<MESH_TYPE>::AddFaces(mesh, nFaces);
					unsigned int f0=0;
					for (unsigned int f=0; f < nFaces; f++)
					{
						f0 = f;
						if (stream.fail())
							return InvalidFile;

						if(cb && (f%1000)==0)
							cb(50+f*50/nFaces,"Face Loading");

						TokenizeNextLine(stream, tokens);
						int vert_per_face = atoi(tokens[0].c_str());
            if(vert_per_face < 3)
              return ErrorDegenerateFace;
						k = 1;
						if (vert_per_face == 3)
						{
							for (int j = 0; j < 3; j++)
							{
                                if (k == tokens.size())   // if EOL 		// Go to next line when needed
								{
									TokenizeNextLine(stream, tokens);
                                    if (tokens.size() == 0) return InvalidFile; // if EOF
									k = 0;
								}
							
								mesh.face[f].V(j) = &(mesh.vert[ atoi(tokens[k].c_str()) ]);
								k++;
							}
						}
						else
						{
							// The face must be triangulate
							unsigned int trigs = vert_per_face-3; // number of extra faces to add
							nFaces += trigs;
							Allocator<MESH_TYPE>::AddFaces(mesh, trigs);
                            std::vector<int> vertIndices(vert_per_face);
							for (int j=0; j < vert_per_face; j++)
							{
                                if (k == tokens.size())   // if EOL // Go to next line when needed
								{
									TokenizeNextLine(stream, tokens);
                                    if (tokens.size() == 0) return InvalidFile; // if EOF
                                    k = 0;
								}
								vertIndices[j] = atoi(tokens[k].c_str());
								k++;
							}
                            if(vert_per_face==4)
                            {   // To well triangulate the quad:
                                // if the quad is flat and convex we use the shortest diag ,
                                // else we should use the diag that makes the smallest                                
                                
                                const CoordType &P0=mesh.vert[vertIndices[0]].cP();
                                const CoordType &P1=mesh.vert[vertIndices[1]].cP();
                                const CoordType &P2=mesh.vert[vertIndices[2]].cP();
                                const CoordType &P3=mesh.vert[vertIndices[3]].cP();

                                CoordType N00 = Normal(P0,P1,P2);
                                CoordType N01 = Normal(P0,P2,P3);
                                CoordType N10 = Normal(P1,P2,P3);
                                CoordType N11 = Normal(P1,P3,P0);

                                ScalarType Angle0Rad=Angle(N00,N01);
                                ScalarType Angle1Rad=Angle(N10,N11);

                                // QualityRadii is inradius/circumradius; bad when close to zero. 
                                // swap diagonal if the worst triangle improve. 
                                bool qualityImprove = std::min(QualityRadii(P0,P1,P2),QualityRadii(P0,P2,P3)) < std::min(QualityRadii(P1,P2,P3),QualityRadii(P1,P3,P0));
                                bool swapCauseFlip = (Angle1Rad > M_PI/2.0) && (Angle0Rad <M_PI/2.0);
                                if(qualityImprove && ! swapCauseFlip)
                                        std::rotate(vertIndices.begin(), vertIndices.begin()+1, vertIndices.end());

                            }
                            // standard fan triangulation (we hope the polygon is convex...)
                            for (int j=0; j<=vert_per_face-3; j++)
                            {
                                mesh.face[f+j].V(0) = &(mesh.vert[ vertIndices[0  ] ]);
                                mesh.face[f+j].V(1) = &(mesh.vert[ vertIndices[1+j] ]);
                                mesh.face[f+j].V(2) = &(mesh.vert[ vertIndices[2+j] ]);
                                if (tri::HasPerFaceFlags(mesh)) {
                                // tag internal polygonal edges as "faux"
                                  if (j>0) mesh.face[f+j].SetF(0);
                                  if (j<vert_per_face-3) mesh.face[f+j].SetF(2);
                                  loadmask |= Mask::IOM_BITPOLYGONAL;
                              }
                            }
                            f+=trigs;
						}

						// NOTE: It is assumed that colored face takes exactly one text line
						//       (otherwise it is impossible to parse color information since
						//        color components can vary)
						if (isColorDefined && tri::HasPerFaceColor(mesh)) 
						{
							size_t color_elements = tokens.size() - vert_per_face-1;

							// set per-face color attribute
							if (color_elements > 0)
								loadmask |= Mask::IOM_FACECOLOR;

							switch (color_elements)
							{
							case 0:
								{
									for ( ; f0<=f; f0++)
										mesh.face[f0].C().Import(vcg::Color4f(.666f, .666f, .666f, .666f));
									break;
								}
							case 1:
								{
									for ( ; f0<=f; f0++)
										mesh.face[f0].C().Import( ColorMap( atoi(tokens[vert_per_face+1].c_str()) ) );
									break;
								}
							case 3:
								{
									if (tokens[vert_per_face+1].find('.')==size_t(-1))
									{
										int rgb[3];
										rgb[0] = atoi( tokens[vert_per_face+1].c_str() );
										rgb[1] = atoi( tokens[vert_per_face+2].c_str() );
										rgb[2] = atoi( tokens[vert_per_face+3].c_str() );
										for ( ; f0<=f; f0++)
											mesh.face[f0].C().SetRGB(rgb[0], rgb[1], rgb[2]);
									}
									else
									{
										float color[3];
										color[0] = (float) atof( tokens[vert_per_face+1].c_str() );
										color[1] = (float) atof( tokens[vert_per_face+2].c_str() );
										color[2] = (float) atof( tokens[vert_per_face+3].c_str() );
										for ( ; f0<=f; f0++)
											mesh.face[f0].C().Import(vcg::Color4f(color[0], color[1], color[2], 1.0f));
									}
									break;
								}
							case 4:
								{
									if (tokens[vert_per_face+1].find('.')==0) // if it is a float there is a dot
									{
										unsigned char color[4];
										color[0] = (unsigned char) atoi(tokens[vert_per_face+1].c_str());
										color[1] = (unsigned char) atoi(tokens[vert_per_face+2].c_str());
										color[2] = (unsigned char) atoi(tokens[vert_per_face+3].c_str());
										color[3] = (unsigned char) atoi(tokens[vert_per_face+4].c_str());
										for ( ; f0<=f; f0++)
											mesh.face[f0].C().Import(vcg::Color4f(color[0], color[1], color[2], color[3]));
									}
									else
									{
										float color[4];
										color[0] = float( atof(tokens[vert_per_face+1].c_str()) );
										color[1] = float( atof(tokens[vert_per_face+2].c_str()) );
										color[2] = float( atof(tokens[vert_per_face+3].c_str()) );
										color[3] = float( atof(tokens[vert_per_face+4].c_str()) );
										for ( ; f0<=f; f0++)
											mesh.face[f0].C().Import(vcg::Color4f(color[0], color[1], color[2], color[3]));
									}
									break;
								}
							} //end switch
						} // end if (isColorDefined)
					} // end of for f=...

					return NoError;

				} // end Open

    protected:
			
				/*!
				* Read the next valid line and parses it into "tokens", allowing the tokens to be read one at a time.
				* \param stream	The object providing the input stream
				*	\param tokens	The "tokens" in the next line
				*/
				inline static const void TokenizeNextLine(std::ifstream &stream, std::vector< std::string > &tokens)
				{
					std::string line;
					do
						std::getline(stream, line, '\n');
					while (line[0] == '#' || line.length()==0);

					size_t from = 0; 
					size_t to = 0;
					size_t length = line.size();
					tokens.clear();
					do
					{
						while ( (line[from]==' ' || line[from] == '\t'  || line[from] == '\r') && from!=length)
							from++;
						if(from!=length)
						{
							to = from+1;
							while ( (line[to]!=' ' && line[to] != '\t'  || line[to] == '\r') && to!=length)
								to++;
							tokens.push_back(line.substr(from, to-from).c_str());
							from = to;
						}
					}
					while (from<length);
				} // end Tokenize

				/*!
				*	Provide the int->color mapping, according to the Geomview's `cmap.fmap' file.
				*	\param		i	the color index
				*	\return			the corresponding <CODE>vcg::Color4f</CODE> color
				*/
				static const vcg::Color4f ColorMap(int i) 
				{
					static const float colorMap[148][4] = 
					{
						{ 1.0f,	 1.0f,	 1.0f,	 1.0f	 },
						{ 1.0f,	 1.0f,	 1.0f,	 1.0f	 },
						{ 1.0f,	 1.0f,	 1.0f,	 1.0f	 },
						{ 1.0f,	 1.0f,	 1.0f,	 1.0f	 },
						{ 1.0f,	 1.0f,	 1.0f,	 1.0f	 },
						{ 1.0f,	 1.0f,	 1.0f,	 1.0f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.2f,	 0.2f,	 0.2f,	 0.2f	 },
						{ 0.9f,	 0.9f,	 0.9f,	 0.9f	 },
						{ 0.1f,	 0.1f,	 0.1f,	 0.1f	 },
						{ 0.1f,	 0.1f,	 0.1f,	 0.1f	 },
						{ 0.8f,	 0.8f,	 0.8f,	 0.8f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.9f,	 0.9f,	 0.9f,	 0.9f	 },
						{ 0.2f,	 0.2f,	 0.2f,	 0.2f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.75f, 0.75f,	 0.75f,	 0.75f },
						{ 0.8f,	 0.8f,	 0.8f,	 0.8f	 },
						{ 0.8f,	 0.8f,	 0.8f,	 0.8f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.4f,	 0.4f,	 0.4f,	 0.4f	 },
						{ 0.4f,	 0.4f,	 0.4f,	 0.4f	 },
						{ 0.8f,	 0.8f,	 0.8f,	 0.8f	 },
						{ 0.8f,	 0.8f,	 0.8f,	 0.8f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.9f,	 0.9f,	 0.9f,	 0.9f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.75f, 0.75f,	 0.75f,	 0.75f },
						{ 0.8f,	 0.8f,	 0.8f,	 0.8f	 },
						{ 0.4f,	 0.4f,	 0.4f,	 0.4f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.4f,	 0.4f,	 0.4f,	 0.4f	 },
						{ 0.8f,	 0.8f,	 0.8f,	 0.8f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.9f,	 0.9f,	 0.9f,	 0.9f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.75f, 0.75f,	 0.75f,	 0.75f },
						{ 0.8f,	 0.8f,	 0.8f,	 0.8f	 },
						{ 0.4f,	 0.4f,	 0.4f,	 0.4f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.4f,	 0.4f,	 0.4f,	 0.4f	 },
						{ 0.8f,	 0.8f,	 0.8f,	 0.8f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.9f,	 0.9f,	 0.9f,	 0.9f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.75f, 0.75f,	 0.75f,	 0.75f },
						{ 0.8f,	 0.8f,	 0.8f,	 0.8f	 },
						{ 0.4f,	 0.4f,	 0.4f,	 0.4f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.4f,	 0.4f,	 0.4f,	 0.4f	 },
						{ 0.8f,	 0.8f,	 0.8f,	 0.8f	 },
						{ 1.0f,	 1.0f,	 1.0f,	 1.0f	 },
						{ 1.0f,	 1.0f,	 1.0f,	 1.0f	 },
						{ 1.0f,	 1.0f,	 1.0f,	 1.0f	 },
						{ 1.0f,	 1.0f,	 1.0f,	 1.0f	 },
						{ 1.0f,	 1.0f,	 1.0f,	 1.0f	 },
						{ 1.0f,	 1.0f,	 1.0f,	 1.0f	 },
						{ 0.05f, 0.05f,	 0.05f,	 0.05f },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.2f,	 0.2f,	 0.2f,	 0.2f	 },
						{ 0.9f,	 0.9f,	 0.9f,	 0.9f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.1f,	 0.1f,	 0.1f,	 0.1f	 },
						{ 0.8f,	 0.8f,	 0.8f,	 0.8f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.9f,	 0.9f,	 0.9f,	 0.9f	 },
						{ 0.9f,	 0.9f,	 0.9f,	 0.9f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.75f, 0.75f,	 0.75f,	 0.75f },
						{ 0.75f, 0.75f,	 0.75f,	 0.75f },
						{ 0.8f,	 0.8f,	 0.8f,	 0.8f	 },
						{ 0.8f,	 0.8f,	 0.8f,	 0.8f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.4f,	 0.4f,	 0.4f,	 0.4f	 },
						{ 0.4f,	 0.4f,	 0.4f,	 0.4f	 },
						{ 0.8f,	 0.8f,	 0.8f,	 0.8f	 },
						{ 0.8f,	 0.8f,	 0.8f,	 0.8f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.9f,	 0.9f,	 0.9f,	 0.9f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.75f, 0.75f,	 0.75f,	 0.75f },
						{ 0.8f,	 0.8f,	 0.8f,	 0.8f	 },
						{ 0.4f,	 0.4f,	 0.4f,	 0.4f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.4f,	 0.4f,	 0.4f,	 0.4f	 },
						{ 0.8f,	 0.8f,	 0.8f,	 0.8f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.9f,	 0.9f,	 0.9f,	 0.9f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.75f, 0.75f,	 0.75f,	 0.75f },
						{ 0.8f,	 0.8f,	 0.8f,	 0.8f	 },
						{ 0.4f,	 0.4f,	 0.4f,	 0.4f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.4f,	 0.4f,	 0.4f,	 0.4f	 },
						{ 0.8f,	 0.8f,	 0.8f,	 0.8f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.7f,	 0.7f,	 0.7f,	 0.7f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.9f,	 0.9f,	 0.9f,	 0.9f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.75f, 0.75f,	 0.75f,	 0.75f },
						{ 0.8f,	 0.8f,	 0.8f,	 0.8f	 },
						{ 0.4f,	 0.4f,	 0.4f,	 0.4f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.0f,	 0.0f,	 0.0f,	 0.0f	 },
						{ 0.4f,	 0.4f,	 0.4f,	 0.4f	 },
						{ 0.8f,	 0.8f,	 0.8f,	 0.8f	 }
					};
					return Color4f(colorMap[i][0], colorMap[i][1], colorMap[i][2], colorMap[i][3]);
				}
			};
			// /*! @} */
		} //namespace io
	}//namespace tri
} // namespace vcg

#endif //__VCGLIB_IMPORT_OFF
