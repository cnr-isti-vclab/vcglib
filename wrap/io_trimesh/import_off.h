/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
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
#ifndef __VCGLIB_IMPORT_OFF
#define __VCGLIB_IMPORT_OFF

#include <cstdlib>
#include <fstream>
#include <limits>
#include<vcg/complex/algorithms/bitquad_support.h>
#include <vcg/space/planar_polygon_tessellation.h>
#include <wrap/io_trimesh/io_mask.h>
#include <wrap/io_trimesh/c_numeric_locale.h>

namespace vcg {
	namespace tri {
		namespace io {

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
				enum OFFCodes {
					NoError = 0, CantOpen, InvalidFile,
					InvalidFile_MissingOFF,
					UnsupportedFormat, ErrorNotTriangularFace, ErrorHighDimension, ErrorDegenerateFace
				};

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
					  "Unsupported format", "Face with more than 3 vertices","File with high dimensional vertexes are not supported", "Error Degenerate Face with less than 3 vertices" };

					if (message_code > 7 || message_code < 0)
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
					loadmask = 0;
					MESH_TYPE dummyMesh;
					return (Open(dummyMesh, filename, loadmask) == NoError);
				}

				static int Open(MESH_TYPE &mesh, const char *filename, CallBackPos *cb = 0)
				{
					int loadmask;
					return Open(mesh, filename, loadmask, cb);
				}

				static int OpenMem(MESH_TYPE &mesh, const char *mem, size_t sz, int &loadmask,
					CallBackPos *cb = 0)
				{
					std::string str;
					str.append(mem, sz);
					std::istringstream strm(str);
					return OpenStream(mesh, strm, loadmask, cb);
				}

				/*!
							   *  Standard call for reading a mesh.
							   *
							   *  \param mesh         the destination mesh
							   *  \param filename     the name of the file to read from
							   *  \return             the operation result
							   */
				static int Open(MESH_TYPE &mesh, const char *filename, int &loadmask,
					CallBackPos *cb = 0)
				{
					std::ifstream stream(filename);
					if (stream.fail())
						return CantOpen;
					return OpenStream(mesh, stream, loadmask, cb);
				}

				static int OpenStream(MESH_TYPE &mesh, std::istream &stream, int &loadmask,
					CallBackPos *cb = 0)
				{
					CNumericLocale cLocale; // numbers in files always use a dot, whatever the host's locale
					std::vector< std::string > tokens;
					TokenizeNextLine(stream, tokens);
					if (tokens.empty()) return InvalidFile_MissingOFF;

					bool isNormalDefined = false;
					bool isColorDefined = false;
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
						for (int u = static_cast<int>(header.rfind("OFF") - 1); u >= 0; u--)
						{
							if (header[u] == 'C')
								isColorDefined = true;
							else if (header[u] == 'N')
								isNormalDefined = true;
							else if (u > 0 && header[u - 1] == 'S' && header[u] == 'T')
								isTexCoordDefined = true;
							else if (header[u] == '4')
								homogeneousComponents = true;
							else if (header[u] == 'n')
								return ErrorHighDimension;
						}
					}
					else return InvalidFile_MissingOFF;

					// If the file is slightly malformed and it has nvert and nface AFTER the OFF string instead of in the next line
					// we manage it here...
					if (tokens.size() == 1) TokenizeNextLine(stream, tokens);
					else tokens.erase(tokens.begin(), tokens.begin() + 1);

					// Update loading mask
					///////////////////////////////////////

					loadmask = Mask::IOM_VERTCOORD | Mask::IOM_FACEINDEX;

					if (isNormalDefined)		loadmask |= Mask::IOM_VERTNORMAL;
					if (isTexCoordDefined)	loadmask |= Mask::IOM_VERTTEXCOORD;
					//if (isColorDefined)			{ loadmask |= Mask::IOM_VERTCOLOR;loadmask |= Mask::IOM_FACECOLOR;}


					//if(onlyMaskFlag) return NoError;


					mesh.Clear();

					// check on next 2 lines to detect corrupted files
					if (tokens.size() < 3)
						return InvalidFile;

					// The counts come from the file and are not trusted: a negative one is refused, and so is
					// a vertex count the rest of the file is too short to hold (a vertex takes at least 6
					// characters, "0 0 0" and a separator), since the vertices are allocated before any is read.
					const long nVerticesRead = strtol(tokens[0].c_str(), nullptr, 10);
					const long nFacesRead = strtol(tokens[1].c_str(), nullptr, 10);
					if (nVerticesRead < 0 || nFacesRead < 0 ||
						nVerticesRead > std::numeric_limits<int>::max() || nFacesRead > std::numeric_limits<int>::max())
						return InvalidFile;
					{
						const std::streampos here = stream.tellg();
						if (here != std::streampos(-1))
						{
							stream.seekg(0, std::ios::end);
							const std::streamoff left = stream.tellg() - here;
							stream.seekg(here);
							if (nVerticesRead > (left + 1) / 6)
								return InvalidFile;
						}
					}
					const unsigned int nVertices = static_cast<unsigned int>(nVerticesRead);
					const unsigned int nFaces = static_cast<unsigned int>(nFacesRead);

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

					for (unsigned int i = 0; i < nVertices; i++, v_iter++)
					{
						if (cb && (i % 1000) == 0)
							cb(i * 50 / nVertices, "Vertex Loading");

						// Read 3 vertex coordinates
						for (unsigned int j = 0; j < 3; j++)
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
							(*v_iter).P()[j] = (ScalarType)atof(tokens[k].c_str());
							k++;
						}

						if (isNormalDefined)
						{
							// Read 3 normal coordinates
							for (unsigned int j = 0; j < 3; j++)
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
								(*v_iter).N()[j] = (ScalarType)atof(tokens[k].c_str());
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
											static_cast<unsigned char>(atoi(tokens[k + 1].c_str()));
										unsigned char b =
											static_cast<unsigned char>(atoi(tokens[k + 2].c_str()));

										vcg::Color4b color(r, g, b, 255);
										(*v_iter).C().Import(color);
									}
									else
									{
										// floats
										float r = static_cast<float>(atof(tokens[k].c_str()));
										float g = static_cast<float>(atof(tokens[k + 1].c_str()));
										float b = static_cast<float>(atof(tokens[k + 2].c_str()));

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
											static_cast<unsigned char>(atoi(tokens[k + 1].c_str()));
										unsigned char b =
											static_cast<unsigned char>(atoi(tokens[k + 2].c_str()));
										unsigned char a =
											static_cast<unsigned char>(atoi(tokens[k + 3].c_str()));

										Color4b color(r, g, b, a);
										(*v_iter).C().Import(color);
									}
									else
									{
										// floats
										float r = static_cast<float>(atof(tokens[k].c_str()));
										float g = static_cast<float>(atof(tokens[k + 1].c_str()));
										float b = static_cast<float>(atof(tokens[k + 2].c_str()));
										float a = static_cast<float>(atof(tokens[k + 3].c_str()));

										vcg::Color4f color(r, g, b, a);
										(*v_iter).C().Import(color);
									}
								}
							}

							k += nb_color_components;
						}

						if (isTexCoordDefined)
						{
							for (unsigned int j = 0; j < 2; j++)
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
					// Every count and index a face line gives is checked before it is used: a face the
					// header counts that the file does not hold, a line with fewer indices than it declares,
					// and an index outside the vertices all refuse the file, as the PLY importer does
					// (Allocator::AddFace checks the range only by assertion). A triangle naming a vertex
					// twice is skipped, as the OBJ importer does, since AddFace asserts that its vertices
					// differ; addTriangle says whether it added one.
					auto validIndex = [&](int i) { return i >= 0 && static_cast<unsigned int>(i) < nVertices; };
					auto addTriangle = [&](int i0, int i1, int i2) {
						if (i0 == i1 || i1 == i2 || i0 == i2)
							return false;
						Allocator<MESH_TYPE>::AddFace(mesh, size_t(i0), size_t(i1), size_t(i2));
						return true;
					};
					if (FaceType::HasPolyInfo())
					{
						for (unsigned int f = 0; f < nFaces; f++)
						{
							if (cb && (f % 1000) == 0) cb(50 + f * 50 / nFaces, "Face Loading");
							TokenizeNextLine(stream, tokens);
							if (tokens.empty())
								return InvalidFile; // the header counts more faces than the file holds
							const int vert_per_face = atoi(tokens[0].c_str());
							if (vert_per_face < 0)
								return InvalidFile;
							std::vector<int> vInd;
							k = 1;
							for (int j = 0; j < vert_per_face; j++)
							{
								if (k == tokens.size())   // if EOL // Go to next line when needed
								{
									TokenizeNextLine(stream, tokens);
									if (tokens.size() == 0) return InvalidFile; // if EOF
									k = 0;
								}
								vInd.push_back(atoi(tokens[k].c_str()));
								if (!validIndex(vInd.back()))
									return InvalidFile;
								k++;
							}
							// A polygon of any number of corners, as the file gives it: polygons of more than
							// four corners used to be dropped.
							if (vert_per_face >= 3)
							{
								FaceIterator fi = Allocator<MESH_TYPE>::AddFaces(mesh, 1);
								(*fi).Alloc(vert_per_face);
								for (int j = 0; j < vert_per_face; ++j)
									(*fi).V(j) = &mesh.vert[vInd[j]];
							}

						}
					}
					else // Standard Triangular Mesh Loading
					{
						// Initial call to the QuadTriangulate with an empty vector to just reset the static set of existing diagonals
						std::vector<VertexPointer> qtmp;
						BitQuad<MESH_TYPE>::QuadTriangulate(qtmp);

						for (unsigned int ff = 0; ff < nFaces; ff++)
						{
							if (stream.fail())
								return InvalidFile;

							if (cb && (ff % 1000) == 0)
								cb(50 + ff * 50 / nFaces, "Face Loading");

							TokenizeNextLine(stream, tokens);
							if (tokens.empty())
								return InvalidFile; // the header counts more faces than the file holds
							const int vert_per_face = atoi(tokens[0].c_str());
							if (vert_per_face < 2)
								return ErrorDegenerateFace;
							k = 1;
							if(tokens.size()==1) // to handle the case when you have the number of vert is on a different line
								TokenizeNextLine(stream, tokens, true);
							if (tokens.size() < static_cast<size_t>(vert_per_face) + 1)
								return InvalidFile; // fewer indices than the face declares
							std::vector<int> vertIndices(vert_per_face);
							for (int j = 0; j < vert_per_face; j++)
							{
								vertIndices[j] = atoi(tokens[j + 1].c_str());
								if (!validIndex(vertIndices[j]))
									return InvalidFile;
							}
							// The faces this line adds, which its colour, if any, goes to.
							const size_t firstFace = mesh.face.size();
							// A two-corner face is an edge, kept only by a mesh whose edges have vertices (in any
							// other, setting them wrote through the empty component's dummy reference), and only
							// if its ends differ, which AddEdge asserts.
							if (vert_per_face == 2 && tri::HasEVAdjacency(mesh) && vertIndices[0] != vertIndices[1])
								Allocator<MESH_TYPE>::AddEdge(mesh, size_t(vertIndices[0]), size_t(vertIndices[1]));
							if (vert_per_face == 3)
								addTriangle(vertIndices[0], vertIndices[1], vertIndices[2]);
							if(vert_per_face > 3) // The face must be triangulated
							{
								std::vector<CoordType> polygonPoints(vert_per_face);
								for (int j = 0; j < vert_per_face; j++)
									polygonPoints[j] = mesh.vert[vertIndices[j]].P();
								// Polygonal faces are represented by triangles whose introduced
								// diagonals carry the per-face faux-edge bit. This lets algorithms
								// and exporters recover the original polygon boundary.
								if (tri::HasPerFaceFlags(mesh)) loadmask |= Mask::IOM_BITPOLYGONAL;
								if (vert_per_face == 4)
								{   // To well triangulate use the bitquad support function that reorders vertex for a simple fan
									std::vector<VertexPointer> q(4);
									for (int qqi = 0; qqi < 4; ++qqi)
										q[qqi] = &mesh.vert[vertIndices[qqi]];
									BitQuad<MESH_TYPE>::QuadTriangulate(q);
									for (int qqi = 0; qqi < 4; ++qqi)
										vertIndices[qqi] = static_cast<int>(q[qqi] - &mesh.vert[0]);
									// build a two face fan
									if (addTriangle(vertIndices[0],vertIndices[1],vertIndices[2]) && tri::HasPerFaceFlags(mesh))
										mesh.face.back().SetF(2);
									if (addTriangle(vertIndices[0],vertIndices[2],vertIndices[3]) && tri::HasPerFaceFlags(mesh))
										mesh.face.back().SetF(0);
								}
								else
								{
									// OFF indices are mesh-global, whereas the planar tessellator
									// deliberately returns indices local to polygonPoints. Keep the
									// latter for boundary/faux classification and map through
									// vertIndices only when creating mesh faces.
									std::vector<int> localTriangles;
									if (!TessellatePlanarPolygon3(polygonPoints, localTriangles))
										return InvalidFile;
									for (size_t j = 0; j < localTriangles.size(); j += 3)
									{
										if (addTriangle(
											vertIndices[localTriangles[j]],
											vertIndices[localTriangles[j + 1]],
											vertIndices[localTriangles[j + 2]]) && tri::HasPerFaceFlags(mesh))
										{
											for (int edge = 0; edge < 3; ++edge)
											{
												const int from = localTriangles[j + edge];
												const int to = localTriangles[j + (edge + 1) % 3];
												if ((from + 1) % vert_per_face == to)
													mesh.face.back().ClearF(edge); // original polygon boundary
												else
													mesh.face.back().SetF(edge);   // introduced diagonal
											}
										}
									}
								}
							}

							// NOTE: It is assumed that colored face takes exactly one text line
							//       (otherwise it is impossible to parse color information since
							//        color components can vary)
							const size_t color_elements = tokens.size() - vert_per_face - 1;
							if (color_elements > 0)
							{
								loadmask |= Mask::IOM_FACECOLOR;

								if (tri::HasPerFaceColor(mesh))
								{
									// The colour goes to the faces this line added: as many as the
									// tessellation gave, which is not always vert_per_face - 2.
									auto setColor = [&](const Color4f &c) {
										for (size_t i = firstFace; i < mesh.face.size(); ++i)
											mesh.face[i].C().Import(c);
									};
									const std::string *colorTokens = &tokens[vert_per_face + 1];
									if (color_elements == 1)
										setColor(ColorMap(atoi(colorTokens[0].c_str())));
									else if (color_elements == 3 || color_elements == 4)
									{
										if (colorTokens[0].find('.') == std::string::npos) // if it is a float there is a dot
										{
											Color4b cc(Color4b::White);
											for (size_t c = 0; c < color_elements; ++c)
												cc[c] = (unsigned char)atoi(colorTokens[c].c_str());
											for (size_t i = firstFace; i < mesh.face.size(); ++i)
												mesh.face[i].C() = cc;
										}
										else
										{
											float cf[4] = { 0, 0, 0, 1.0f };
											for (size_t c = 0; c < color_elements; ++c)
												cf[c] = float(atof(colorTokens[c].c_str()));
											setColor(Color4f(cf[0], cf[1], cf[2], cf[3]));
										}
									}
								}
							} // end if (isColorDefined)
						} // end of for f=...
					}
					return NoError;

				} // end Open

			protected:

				/*!
							  * Read the next valid line and parses it into "tokens", allowing the tokens to be read one at a time.
							  * \param stream	The object providing the input stream
							  *	\param tokens	the vector that will contain the tokens
							  * \param appendTokens if true the tokens are appended to the token vector. Default the vector is cleared
							  */
				inline static void TokenizeNextLine(std::istream &stream, std::vector< std::string > &tokens, bool appendTokens=false)
				{
					std::string line;
					do
						std::getline(stream, line, '\n');
					while ((line[0] == '#' || line.length() == 0 || line[0] == '\r') && (!stream.eof()));

					size_t from = 0;
					size_t to = 0;
					size_t length = line.size();
					if(!appendTokens) tokens.clear();
					do
					{
						while (from != length && (line[from] == ' ' || line[from] == '\t' || line[from] == '\r'))
							from++;
						if (from != length)
						{
							to = from + 1;
							while (to != length && (((line[to] != ' ') && (line[to] != '\t')) || (line[to] == '\r')))
								to++;
							tokens.push_back(line.substr(from, to - from).c_str());
							from = to;
						}
					} while (from < length);
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
					// The index comes from the file: one outside the map reads as its first entry, white.
					if (i < 0 || i >= 148)
						i = 0;
					return Color4f(colorMap[i][0], colorMap[i][1], colorMap[i][2], colorMap[i][3]);
				}
			};
			// /*! @} */
		} //namespace io
	}//namespace tri
} // namespace vcg

#endif //__VCGLIB_IMPORT_OFF
