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
/*! \file trimesh_align_pair.cpp
\ingroup code_sample

\brief the minimal example for aligning two meshes

This file contain a minimal example for aligning two meshs.

Example call:
./trimesh_align_pair mesh1.ply mesh2.ply output.ply

output.ply will contain mesh2.ply rotated in order to be aligned to mesh1.ply

*/

#include <vcg/complex/complex.h>

#include <vcg/complex/algorithms/align_pair.h>

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>

class MyFace;
class MyVertex;

struct MyUsedTypes :
		public vcg::UsedTypes<	vcg::Use<MyVertex>::AsVertexType,
		vcg::Use<MyFace>::AsFaceType>
{};

class MyVertex :
		public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::Color4b, vcg::vertex::BitFlags>
{};

class MyFace :
		public vcg::Face < MyUsedTypes, vcg::face::VertexRef, vcg::face::Normal3d, vcg::face::FFAdj, vcg::face::Mark, vcg::face::BitFlags >
{};

class MyMesh :
		public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> >
{};

using namespace vcg;
using namespace std;

std::vector<vcg::Point3d>* vcg::PointMatchingScale::fix;
std::vector<vcg::Point3d>* vcg::PointMatchingScale::mov;
vcg::Box3d vcg::PointMatchingScale::b;

int main(int argc,char ** argv)
{
	if(argc<3) {
		printf("Usage: trimesh_smooth <filename> <filename>\n");
		return 0;
	}
	std::string outputname = "output.ply";
	if (argc == 4){
		outputname = std::string(argv[3]);
	}

	MyMesh m1, m2;

	//open first mesh
	int err = tri::io::Importer<MyMesh>::Open(m1,argv[1]);
	if(err) { // all the importers return 0 in case of success
		printf("Error in reading %s: '%s'\n", argv[1], tri::io::Importer<MyMesh>::ErrorMsg(err));
		exit(-1);
	}

	//open second mesh
	err = tri::io::Importer<MyMesh>::Open(m2,argv[2]);
	if(err) { // all the importers return 0 in case of success
		printf("Error in reading %s: '%s'\n", argv[2], tri::io::Importer<MyMesh>::ErrorMsg(err));
		exit(-1);
	}

	//update normals
	vcg::tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFace(m1);
	vcg::tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFace(m2);

	////PARAMS
	vcg::AlignPair::Result result;
	vcg::AlignPair::Param ap;
	vcg::AlignPair::A2Mesh fix;
	vcg::AlignPair aa;

	// 1) Convert fixed mesh and put it into the grid.
	aa.convertMesh<MyMesh>(m1,fix);

	vcg::AlignPair::A2Grid UG;
	vcg::AlignPair::A2GridVert VG;

	if(m1.fn==0 || ap.UseVertexOnly) {
		fix.initVert(vcg::Matrix44d::Identity());
		vcg::AlignPair::InitFixVert(&fix,ap,VG);
	}
	else {
		fix.init(vcg::Matrix44d::Identity());
		vcg::AlignPair::initFix(&fix, ap, UG);
	}


	// 2) Convert the second mesh and sample a <ap.SampleNum> points on it.
	std::vector<vcg::AlignPair::A2Vertex> tmpmv;
	aa.convertVertex(m2.vert,tmpmv);
	aa.sampleMovVert(tmpmv, ap.SampleNum, ap.SampleMode);

	aa.mov=&tmpmv;
	aa.fix=&fix;
	aa.ap = ap;

	//use identity as first matrix
	vcg::Matrix44d In;
	In.SetIdentity();

	// Perform the ICP algorithm
	aa.align(In,UG,VG,result);

	//rotate m2 using the resulting transformation
	tri::UpdatePosition<MyMesh>::Matrix(m2, result.Tr, true);
	tri::UpdateBounding<MyMesh>::Box(m2);

	//saves the rotated mesh
	tri::io::ExporterPLY<MyMesh>::Save(m2 ,outputname.c_str());

	return 0;
}

