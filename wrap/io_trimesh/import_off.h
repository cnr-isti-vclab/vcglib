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

#ifndef __VCGLIB_IMPORT_OFF
#define __VCGLIB_IMPORT_OFF

#include <fstream>
#include <string>
#include <vector>
#include <assert.h>
#include <vcg/space/color4.h>

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

				/*!
				*	Standard call for knowing the meaning of an error code
				* \param message_code	The code returned by <CODE>Open</CODE>
				*	\return							The string describing the error code 
				*/
				static const char* ErrorMsg(int message_code)
				{
					static const char* error_msg[3] =
					{
						"No errors", "Can't open file", "Premature End of file",
					};

					if(message_code>2 || message_code<0) 
						return "Unknown error";
					else 
						return error_msg[message_code];
				};

				/*!
				*	Standard call for reading a mesh
				*	\param mesh				the destination mesh
				*	\param filename		the name of the file to read from
				*	\return						the operation result
				*/
				static int Open(MESH_TYPE &mesh, const char *filename)
				{
					mesh.Clear();

					bool isNormalDefined	 = false;
					bool isColorDefined		 = false;
					bool isTexCoordDefined = false;
					int	 dimension = 3;

					std::ifstream stream(filename);
					if (stream.fail())
						return OFFCodes::CantOpen;

					std::vector< std::string > tokens;
					TokenizeNextLine(stream, tokens);
					if (tokens[tokens.size()-1].rfind('OFF')!= std::basic_string<char>::npos)
					{
						for (int u=tokens.size()-2; u>=0; u--)
						{
							std::string header = tokens[u];
							if (header.compare("C")==0)
							{
								isColorDefined = true;
								continue;
							}
							if (header.compare("N")==0)
							{
								isNormalDefined = true;
								continue;
							}
							if (header.compare("ST")==0)
							{
								isTexCoordDefined = true;
								continue;
							}
						}
						if (tokens[tokens.size()-1].compare("4OFF")==0)
							dimension = 4;
						else if (tokens[tokens.size()-1].compare("nOFF")==0)
						{
							TokenizeNextLine(stream, tokens);
							dimension = atoi(tokens[0].c_str());
						}
						else
							dimension = 3;

						TokenizeNextLine(stream, tokens);
					}

					unsigned int nVertices, nFaces, nEdges;
					nVertices = atoi(tokens[0].c_str());
					nFaces		= atoi(tokens[1].c_str());
					nEdges		= atoi(tokens[2].c_str());

					assert(dimension = 3);
					VertexIterator v_iter = Allocator<MESH_TYPE>::AddVertices(mesh, nVertices);
					for (unsigned int i=0; i<nVertices; i++, v_iter++)
					{
						if (stream.fail())
							return OFFCodes::UnexpectedEOF;

						TokenizeNextLine(stream, tokens);
						for (unsigned int j=0; j<3; j++)
							(*v_iter).P()[j] = (ScalarType) atof(tokens[j].c_str());

						if (isNormalDefined)
							for (unsigned int j=3; j<6; j++)
								(*v_iter).N()[j] = (ScalarType)  atof(tokens[j].c_str());

						if (isColorDefined) {}

						if (isTexCoordDefined) {}
					}

					Allocator<MESH_TYPE>::AddFaces(mesh, nFaces);
					unsigned int f0=0;
					for (unsigned int f=0; f<nFaces; f++)
					{
						f0 = f;
						if (stream.fail())
							return OFFCodes::UnexpectedEOF;

						
						TokenizeNextLine(stream, tokens);
						int vert_per_face = atoi(tokens[0].c_str());
						if (vert_per_face == 3)
						{
							mesh.face[f].V(0) = &(mesh.vert[ atoi(tokens[1].c_str()) ]);
							mesh.face[f].V(1) = &(mesh.vert[ atoi(tokens[2].c_str()) ]);
							mesh.face[f].V(2) = &(mesh.vert[ atoi(tokens[3].c_str()) ]);
						}
						else
						{
							unsigned int trigs = vert_per_face-3;
							nFaces += trigs;
							Allocator<MESH_TYPE>::AddFaces(mesh, trigs);
							int *vertIndices = new int[vert_per_face];

							for (int k=0; k<vert_per_face; k++)
								vertIndices[k] = atoi(tokens[1+k].c_str());

							for (int k=0; k<=vert_per_face-3; k++)
							{
								mesh.face[f+k].V(0) = &(mesh.vert[ vertIndices[0  ] ]);
								mesh.face[f+k].V(1) = &(mesh.vert[ vertIndices[1+k] ]);
								mesh.face[f+k].V(2) = &(mesh.vert[ vertIndices[2+k] ]);
							}
							f+=trigs;
							delete []vertIndices;
						}

						if (isColorDefined)
						{
							size_t color_elements = tokens.size()-vert_per_face-1;

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
									if (tokens[vert_per_face+1].find('.')==-1)
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
									if (tokens[vert_per_face+1].find('.')==-1)
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
										color[0] = (ScalarType) atof(tokens[vert_per_face+1].c_str());
										color[1] = (ScalarType) atof(tokens[vert_per_face+2].c_str());
										color[2] = (ScalarType) atof(tokens[vert_per_face+3].c_str());
										color[3] = (ScalarType) atof(tokens[vert_per_face+4].c_str());
										for ( ; f0<=f; f0++)
											mesh.face[f0].C().Import(vcg::Color4f(color[0], color[1], color[2], color[3]));
									}
									break;
								}
							} //end switch
						} // end if (isColorDefined)
					}

					return OFFCodes::NoError;
				} // end Open


			protected:
				enum OFFCodes {NoError, CantOpen, UnexpectedEOF};

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

					size_t from		= 0; 
					size_t to			= 0;
					size_t length = line.size();
					tokens.clear();
					do
					{
						while (line[from]==' ' && from!=length)
							from++;
            if(from!=length)
            {
						  to = from+1;
						  while (line[to]!=' ' && to!=length)
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
		}; //namespace io
	};//namespace tri
}; // namespace vcg

#endif //__VCGLIB_IMPORT_OFF