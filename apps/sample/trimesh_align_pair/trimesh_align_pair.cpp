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

\brief the minimal example of using the lib

This file contain a minimal example of the library

*/

#include <vcg/complex/complex.h>

#include <vcg/complex/algorithms/align_pair.h>

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>

class MyFace;
class MyVertex;

struct MyUsedTypes : public vcg::UsedTypes<	vcg::Use<MyVertex>::AsVertexType,    vcg::Use<MyFace>::AsFaceType>{};

class MyVertex  : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3d, vcg::vertex::Normal3d, vcg::vertex::Color4b, vcg::vertex::BitFlags  >{};
class MyFace    : public vcg::Face < MyUsedTypes, vcg::face::VertexRef, vcg::face::Normal3d, vcg::face::FFAdj, vcg::face::Mark, vcg::face::BitFlags > {};
class MyMesh    : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> > {};

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

	vcg::tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFace(m1);
	vcg::tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFace(m2);

	////PARAMS
	/////TODO
	vcg::AlignPair::Result result;
	vcg::AlignPair::Param ap;

	//MovM
	//vcg::Matrix44d FixM=vcg::Matrix44d::Construct(m1.Tr);
	//MovM=vcg::Matrix44d::Construct(m2.Tr);
	//MovM = Inverse(FixM) * MovM;

	vcg::AlignPair::A2Mesh fix;
	vcg::AlignPair aa;

	// 1) Convert fixed mesh and put it into the grid.
	//m1.face.EnableMark();
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
	//MM(movId)->updateDataMask(MeshModel::MM_FACEMARK);
	//m2.face.EnableMark();
	std::vector<vcg::AlignPair::A2Vertex> tmpmv;
	aa.convertVertex(m2.vert,tmpmv);
	aa.sampleMovVert(tmpmv, ap.SampleNum, ap.SampleMode);

	aa.mov=&tmpmv;
	aa.fix=&fix;
	aa.ap = ap;

	vcg::Matrix44d In;
	In.SetIdentity();
	// Perform the ICP algorithm
	aa.align(In,UG,VG,result);

	tri::UpdatePosition<MyMesh>::Matrix(m2, result.Tr, true);
	tri::UpdateBounding<MyMesh>::Box(m2);

	//result.FixName=fixId;
	//result.MovName=movId;
	//result.as.Dump(stdout);

	//saves the rotated mesh
	tri::io::ExporterPLY<MyMesh>::Save(m2 ,"out.ply");

	return 0;
}

