/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2009                                           \/)\/    *
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

#include<vcg/simplex/vertex/base.h>
#include<vcg/simplex/vertex/component.h>

#include <vcg/complex/used_types.h>

#include<vcg/simplex/face/base.h>
#include<vcg/simplex/face/component.h>

#include<vcg/simplex/face/topology.h>
#include<vcg/complex/trimesh/base.h>

// input output
#include<wrap/io_trimesh/import.h>
#include<wrap/io_trimesh/export.h>

// topology computation
#include<vcg/complex/trimesh/update/topology.h>

// normals
#include<vcg/complex/trimesh/update/normal.h> //class UpdateNormals

using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>		::AsVertexType,
																				Use<MyEdge>			::AsEdgeType,
																				Use<MyFace>			::AsFaceType>{};

class MyVertex  : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags  >{};
class MyFace    : public Face< MyUsedTypes, face::FFAdj,  face::VertexRef, face::BitFlags > {};
class MyEdge		: public Edge<MyUsedTypes>{};
class MyMesh    : public vcg::tri::TriMesh< vector<MyVertex>, vector<MyFace> , vector<MyEdge>  > {};

int main( int argc, char **argv )
{
		if(argc<2)
		{
				printf("Usage trimesh_base <meshfilename.ply>\n");
				return -1;
		}

		MyMesh m;

		if(vcg::tri::io::ImporterPLY<MyMesh>::Open(m,argv[1])!=0)
		{
				printf("Error reading file  %s\n",argv[1]);
				exit(0);
		}

		vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
		vcg::tri::UpdateFlags<MyMesh>::FaceBorderFromFF(m);
		vcg::tri::UpdateNormals<MyMesh>::PerVertexNormalized(m);
		printf("Input mesh  vn:%i fn:%i\n",m.vn,m.fn);
		printf( "Mesh has %i vert and %i faces\n", m.vn, m.fn );

    return 0;
}
