/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2012                                           \/)\/    *
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
/*! \file trimesh_normal.cpp
\ingroup code_sample

\brief An example of all the methods for computing normals over a mesh.

*/
#include<vcg/complex/complex.h>

#include<wrap/io_trimesh/import.h>
#include<wrap/io_trimesh/export.h>

#include<vcg/complex/algorithms/clean.h>
#include<vcg/complex/algorithms/isotropic_remeshing.h>

using namespace vcg;
using namespace std;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>   ::AsVertexType,
        Use<MyEdge>     ::AsEdgeType,
        Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::VFAdj, vertex::Qualityf, vertex::BitFlags,  vertex::Mark>{};
class MyFace    : public Face< MyUsedTypes, face::Mark,  face::VertexRef, face::VFAdj, face::FFAdj, face::Normal3f, face::BitFlags > {};
class MyEdge    : public Edge<MyUsedTypes>{};
class MyMesh    : public tri::TriMesh< vector<MyVertex>, vector<MyFace> , vector<MyEdge>  > {};

int main( int argc, char **argv )
{
	MyMesh original,toremesh;
	if(argc<2)
	{
		printf("Usage: trimesh_remesh <filename> [[targetLen (bbox perc)] [iterNum] [creaseAngle] [maxSurfDist (bbox perc)]]");
		exit(0);
	}

	if(tri::io::Importer<MyMesh>::Open(original,argv[1])!=0)
	{
		printf("Error reading file  %s\n",argv[1]);
		exit(0);
	}
	float targetLenPerc=.2f;
	int iterNum=20;
	float creaseAngle = 30.f;
	float maxSurfDistPerc = 0.001f;
	if(argc>=3) targetLenPerc = atof(argv[2]);
	if(argc>=4) iterNum = atoi(argv[3]);
	if(argc>=5) creaseAngle = atof(argv[4]);
	if(argc>=6) maxSurfDistPerc = atof(argv[5]);


	// Mesh cleaning
	tri::Clean<MyMesh>::RemoveUnreferencedVertex(original);
	Allocator<MyMesh>::CompactEveryVector(original);


	tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFaceNormalized(original);
	tri::UpdateBounding<MyMesh>::Box(original);

	vcg::tri::Append<MyMesh,MyMesh>::MeshCopy(toremesh,original);
	tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFaceNormalized(toremesh);
	tri::UpdateBounding<MyMesh>::Box(toremesh);

	tri::UpdateTopology<MyMesh>::FaceFace(toremesh);
	float lengthThr = targetLenPerc*(original.bbox.Diag()/100.f);
	float maxSurfDist = maxSurfDistPerc*(original.bbox.Diag()/100.f);
	printf("Length Thr: %8.3f ~ %4.2f %% on %5.3f\n",lengthThr,targetLenPerc,original.bbox.Diag());

	IsotropicRemeshing<MyMesh>::Params params;
	params.SetTargetLen(lengthThr);
	params.SetFeatureAngleDeg(creaseAngle);
	params.iter=iterNum;

	if (maxSurfDistPerc != 0)
	{
		params.surfDistCheck = true;
		params.maxSurfDist = maxSurfDist;
	}
	else
	{
		params.surfDistCheck = false;
	}

	params.cleanFlag = true;
	params.userSelectedCreases = false;



	printf(" Input mesh %8i v %8i f\n",toremesh.VN(),toremesh.FN());
	IsotropicRemeshing<MyMesh>::Do(toremesh, original, params);
	vcg::tri::io::ExporterPLY<MyMesh>::Save(toremesh, "remesh.ply");
	printf("Output mesh %8i v %8i f\n",toremesh.VN(),toremesh.FN());

	return 0;
}
