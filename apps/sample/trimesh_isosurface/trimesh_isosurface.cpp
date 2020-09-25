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
#include <vcg/complex/algorithms/create/marching_cubes.h>
#include <vcg/complex/algorithms/create/mc_trivial_walker.h>
#include <wrap/io_trimesh/export_ply.h>
#include <wrap/io_trimesh/export_stl.h>
#include <wrap/io_trimesh/export_obj.h>
#include <wrap/io_trimesh/export_off.h>

using namespace std;
using namespace vcg;

typedef float ScalarType;

class MyFace;
class MyVertex;

struct MyUsedTypes : public UsedTypes<	Use<MyVertex>::AsVertexType,Use<MyFace>::AsFaceType>{};

class MyVertex     : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::BitFlags>{};
class MyFace       : public Face< MyUsedTypes, face::VertexRef, face::BitFlags> {};
class MyMesh	   : public vcg::tri::TriMesh< std::vector< MyVertex>, std::vector< MyFace > > {};

typedef SimpleVolume<SimpleVoxel<float> > MyVolume;

#define SZ 256
#define NOISE_COEFF 2
#define LEN 45

int main(int /*argc*/ , char **/*argv*/)
{
  MyVolume	volume;

  typedef vcg::tri::TrivialWalker<MyMesh,MyVolume>	MyWalker;
  typedef vcg::tri::MarchingCubes<MyMesh, MyWalker>	MyMarchingCubes;
  MyWalker walker;

  // Simple initialization of the volume with some perlin noise
  vcg::Box3f bb(vcg::Point3f(-1,-1,-1),vcg::Point3f(1,1,1));

  volume.Init(Point3i(SZ,SZ,SZ),bb);
  for(short i=0;i<SZ;i++)
    for(short j=0;j<SZ;j++)
      for(short k=0;k<SZ;k++)
        volume.Val(i,j,k)=(j-SZ/2)*(j-SZ/2)+(k-SZ/2)*(k-SZ/2)  + i*NOISE_COEFF*(float)math::Perlin::Noise(i*.2,j*.2,k*.2);

	// MARCHING CUBES
	MyMesh		mc_mesh;
        short mask = 0;

	printf("[MARCHING CUBES] Building mesh...");
	MyMarchingCubes	mc(mc_mesh, walker);
        walker.BuildMesh<MyMarchingCubes>(mc_mesh, volume, mc, LEN*LEN);

    vcg::tri::io::ExporterPLY<MyMesh>::Save( mc_mesh, "marching_cubes.ply");
    vcg::tri::io::ExporterSTL<MyMesh>::Save( mc_mesh, "marching_cubes.stl");
    vcg::tri::io::ExporterOBJ<MyMesh>::Save( mc_mesh, "marching_cubes.obj",mask);
    vcg::tri::io::ExporterOFF<MyMesh>::Save( mc_mesh, "marching_cubes.off");
	
    printf("OK!\n");
};
