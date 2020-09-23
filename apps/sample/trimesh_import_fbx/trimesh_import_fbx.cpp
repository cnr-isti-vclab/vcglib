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
/*! \file trimesh_create.cpp
\ingroup code_sample

\brief A very simple example that open a fbx file and save it as an obj. 
Just as an example of what you should add to your project to compile it:
  ../../../wrap/openfbx/src/ofbx.cpp 
  ../../../wrap/openfbx/src/miniz.c


*/
#include <vcg/complex/complex.h>
#include <wrap/io_trimesh/import_fbx.h>
#include <wrap/io_trimesh/export_obj.h>
#include <stdio.h>

using namespace vcg;

class MyFace;
class MyEdge;
class MyVertex;

struct MyUsedTypes : public UsedTypes<Use<MyVertex>::AsVertexType, Use<MyFace>::AsFaceType>{};

class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::Qualityf, vertex::Color4b,vertex::BitFlags  >{};
class MyFace    : public Face  < MyUsedTypes, face::VertexRef,  face::Normal3f, face::Qualityf, face::WedgeTexCoord2f, face::Color4b, face::BitFlags > {};
class MyMesh    : public tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> >{};

int  main()
{
  MyMesh openMesh;
  
//    tri::io::ImporterFBX<MyMesh>::Open(openMesh,"Arabic_Censer/Arabic_Censer.FBX");
//  tri::io::ImporterFBX<MyMesh>::Open(openMesh,"elephant/ELEPHANT_M.fbx");
  tri::io::ImporterFBX<MyMesh>::Open(openMesh,"liontemple/sketchfabTemp.obj.fbx");
 
  tri::io::ExporterOBJ<MyMesh>::Save(openMesh,"sphere.obj",tri::io::Mask::IOM_WEDGTEXCOORD);
  return 0;
}
