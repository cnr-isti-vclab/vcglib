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
#include <vcg/complex/complex.h>
#include <vcg/math/perlin_noise.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/create/resampler.h>
#include <wrap/io_trimesh/export_ply.h>

using namespace std;
using namespace vcg;

typedef float ScalarType;

class MyFace;
class MyVertex;

struct MyUsedTypes : public UsedTypes<	Use<MyVertex>::AsVertexType,
                                        Use<MyFace>  ::AsFaceType>{};

class MyVertex     : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags>{};
class MyFace       : public Face< MyUsedTypes, face::VertexRef, face::Normal3f, face::BitFlags> {};

class MyMesh       : public tri::TriMesh< std::vector< MyVertex>, std::vector< MyFace > > {};

int main(int /*argc*/ , char **/*argv*/)
{
  MyMesh base_mesh,resampled_mesh;
  tri::Torus(base_mesh,10,3);
  vcg::tri::UpdateBounding<MyMesh>::Box(base_mesh);
  Box3f bb = base_mesh.bbox; 
  float cell_side = bb.Diag()/30.0;
  bb.Offset(cell_side);
  Point3i box_size(bb.DimX()/cell_side,bb.DimY()/cell_side,bb.DimZ()/cell_side);
  
  tri::Resampler<MyMesh,MyMesh>::Resample(base_mesh,resampled_mesh,bb,box_size,cell_side*5);
  
  vcg::tri::io::ExporterPLY<MyMesh>::Save( resampled_mesh, "resampled_torus.ply");
  
  printf("OK!\n");
};
