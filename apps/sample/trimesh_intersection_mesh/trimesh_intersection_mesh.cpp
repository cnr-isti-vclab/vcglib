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
/*! \file trimesh_intersection_mesh.cpp
\ingroup code_sample

\brief An example of adaptively remeshing the intersection between a cube and a sphere

*/

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/intersection.h>
#include <vcg/complex/algorithms/refine.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include<vcg/complex/algorithms/isotropic_remeshing.h>

#include <wrap/io_trimesh/export_off.h>

using namespace std;
using namespace vcg;

class MyFace;
class MyVertex;

struct MyUsedTypes : public UsedTypes<	
    Use<MyVertex>::AsVertexType,
    Use<MyFace>  ::AsFaceType>{};


class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::VFAdj, vertex::BitFlags, vertex::Normal3f, vertex::Qualityf, vertex::Mark>{};
class MyFace    : public Face  < MyUsedTypes, face::VertexRef, face::VFAdj, face::FFAdj, face::Color4b, face::BitFlags, face::Mark, face::Normal3f> {};

class MyMesh:     public tri::TriMesh< vector<MyVertex>, vector<MyFace > >{};

int main(int ,char **)
{  
  MyMesh m1, m2;
  Sphere(m1);
  Hexahedron(m2);
  printf("Created a sphere mesh of %i vertices %i edges\n",m1.VN(),m1.FN());
  printf("Created a cube mesh of %i vertices %i edges\n",m2.VN(),m2.FN());
  
  math::MarsenneTwisterRNG rnd(clock());
  Point3f direction = vcg::math::GeneratePointOnUnitSphereUniform<float,math::MarsenneTwisterRNG>(rnd);
  tri::UpdatePosition<MyMesh>::Translate(m1,direction);
  for(int i=0;i<6;++i)
  {
    tri::Clean<MyMesh>::SelectIntersectingFaces(m1,m2);
    tri::UpdateSelection<MyMesh>::FaceDilate(m1);
    tri::Clean<MyMesh>::SelectIntersectingFaces(m2,m1);
    tri::UpdateSelection<MyMesh>::FaceDilate(m2);
    IsotropicRemeshing<MyMesh>::Params params;
    
    float len = (tri::Stat<MyMesh>::ComputeFaceEdgeLengthAverage(m1,true) + tri::Stat<MyMesh>::ComputeFaceEdgeLengthAverage(m1,true));
    params.SetTargetLen(len*0.8f);
    params.SetFeatureAngleDeg(10);
    params.iter=1; // just one iteration to avoid overtessellating.
    params.selectedOnly=true;
    printf(" Input mesh %8i v %8i f\n",m1.VN(),m1.FN());
    IsotropicRemeshing<MyMesh>::Do(m1, params);
    IsotropicRemeshing<MyMesh>::Do(m2, params);
    printf(" Input mesh %8i v %8i f\n",m1.VN(),m1.FN());
  }
  tri::Clean<MyMesh>::SelectIntersectingFaces(m1,m2);
  tri::Clean<MyMesh>::SelectIntersectingFaces(m2,m1);  
  int selCnt= tri::UpdateColor<MyMesh>::PerFaceConstant(m1,Color4b::Red,true);
  printf("Intersected %i faces on sphere\n",selCnt);
  selCnt= tri::UpdateColor<MyMesh>::PerFaceConstant(m2,Color4b::Red,true);
  printf("Intersected %i faces on cube\n",selCnt);
  
  tri::io::ExporterOFF<MyMesh>::Save(m1,"sphere.off",tri::io::Mask::IOM_FACECOLOR);
  tri::io::ExporterOFF<MyMesh>::Save(m2,"cube.off",tri::io::Mask::IOM_FACECOLOR);
  
  return 0;
}

