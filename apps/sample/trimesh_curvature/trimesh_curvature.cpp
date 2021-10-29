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
/*! \file trimesh_curvature.cpp
\ingroup code_sample

\brief an example showing various techniques for computing curvatures

*/
#include <vcg/complex/complex.h>

#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/curvature.h>
#include <vcg/complex/algorithms/update/curvature_fitting.h>

#include<wrap/io_trimesh/export_ply.h>

class MyFace;
class MyVertex;
struct MyUsedTypes : public vcg::UsedTypes<	vcg::Use<MyVertex>   ::AsVertexType,
                                            vcg::Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public vcg::Vertex<MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::VFAdj, vcg::vertex::Qualityf, vcg::vertex::CurvatureDirf, vcg::vertex::BitFlags  >{};
class MyFace    : public vcg::Face< MyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::Normal3f, vcg::face::VertexRef, vcg::face::BitFlags > {};
class MyMesh    : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> > {};

int main( int /*argc*/, char **/*argv*/ )
{
  MyMesh m;
  // in a torus with radii 1 and 4
  // on the outside the principal curvature should be 1 and  1/5 
  // and internally the principal curvature  should be 1 and -1/3
  // Gaussian range -0.333 ..  0.200
  // mean     range  0.333 ..  0.600
  //vcg::tri::Torus(m,4,1,32,16);
  
  
  // in a sphere of radius 2 the curvature is everywhere 0.5
  // Gaussian 0.25
  // Mean     0.5
  vcg::tri::Sphere(m,5);
  vcg::tri::UpdatePosition<MyMesh>::Scale(m, 2.0);
  
  vcg::tri::UpdateTopology<MyMesh>::FaceFace(m);
  vcg::tri::UpdateTopology<MyMesh>::VertexFace(m);
  vcg::tri::UpdateNormal<MyMesh>::PerVertexNormalizedPerFaceNormalized(m); 
  vcg::Distribution<float> distr;
  printf("Starting mesh  vn:%i fn:%i\n",m.VN(),m.FN());
  
  // Method 1 (discrete - Optimizing 3d triangulations using discrete curvature analysis 2003) 
  vcg::tri::UpdateCurvature<MyMesh>::PerVertexAbsoluteMeanAndGaussian(m);
  vcg::tri::UpdateQuality<MyMesh>::VertexFromAttributeName(m,"KG");
  vcg::tri::Stat<MyMesh>::ComputePerVertexQualityDistribution(m,distr);
  printf("Gaussian Curvature method 1 Min %f Max %f\n",distr.Min(),distr.Max());
  vcg::tri::UpdateQuality<MyMesh>::VertexFromAttributeName(m,"KH");
  vcg::tri::Stat<MyMesh>::ComputePerVertexQualityDistribution(m,distr);
  printf("Mean     Curvature method 1 Min %f Max %f\n",distr.Min(),distr.Max());
  
  // Method 2 (discrete - Discrete Differential-Geometry Operators for Triangulated 2-Manifolds 2002)
  vcg::tri::UpdateCurvature<MyMesh>::MeanAndGaussian(m);
  vcg::tri::UpdateQuality<MyMesh>::VertexFromAttributeName(m,"KG");
  vcg::tri::Stat<MyMesh>::ComputePerVertexQualityDistribution(m,distr);
  printf("Gaussian Curvature method 2 Min %f Max %f\n",distr.Min(),distr.Max());
  vcg::tri::UpdateQuality<MyMesh>::VertexFromAttributeName(m,"KH");
  vcg::tri::Stat<MyMesh>::ComputePerVertexQualityDistribution(m,distr);
  printf("Mean     Curvature method 2 Min %f Max %f\n",distr.Min(),distr.Max());
  vcg::tri::io::ExporterPLY<MyMesh>::Save(m,"Torus_Discrete_Mean2.ply",vcg::tri::io::Mask::IOM_VERTQUALITY);
  
  // Method 3 (directions - Estimating the Tensor of Curvature of a Surface from a Polyhedral Approximation - 1995)
  vcg::tri::UpdateCurvature<MyMesh>::PrincipalDirections(m);
  vcg::tri::UpdateQuality<MyMesh>::VertexGaussianFromCurvatureDir(m);
  vcg::tri::Stat<MyMesh>::ComputePerVertexQualityDistribution(m,distr);
  printf("Gaussian Curvature method 3 Min %f Max %f\n",distr.Min(),distr.Max());
  vcg::tri::UpdateQuality<MyMesh>::VertexMeanFromCurvatureDir(m);
  vcg::tri::Stat<MyMesh>::ComputePerVertexQualityDistribution(m,distr);
  printf("Mean     Curvature method 3 Min %f Max %f\n",distr.Min(),distr.Max());
  
  // Method 4 (directions - Restricted delaunay triangulations and normal cycle )
  vcg::tri::UpdateCurvature<MyMesh>::PrincipalDirectionsNormalCycle(m);
  vcg::tri::UpdateQuality<MyMesh>::VertexGaussianFromCurvatureDir(m);
  vcg::tri::Stat<MyMesh>::ComputePerVertexQualityDistribution(m,distr);
  printf("Gaussian Curvature method 4 Min %f Max %f\n",distr.Min(),distr.Max());
  vcg::tri::UpdateQuality<MyMesh>::VertexMeanFromCurvatureDir(m);
  vcg::tri::Stat<MyMesh>::ComputePerVertexQualityDistribution(m,distr);
  printf("Mean     Curvature method 4 Min %f Max %f\n",distr.Min(),distr.Max());
  
  
  return 0;
}
