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
/*! \file trimesh_intersection_plane.cpp
\ingroup code_sample

\brief An example of computing the intersection of a mesh with a plane, 
saving this polyline as a well projected 2D SVG 
and splicing the mesh in the two components below and over the plane

*/

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/intersection.h>
#include <vcg/complex/algorithms/refine.h>
#include <vcg/complex/algorithms/create/platonic.h>

#include <wrap/io_edgemesh/export_svg.h>
#include <wrap/io_edgemesh/export_dxf.h>
#include <wrap/io_trimesh/export_off.h>

using namespace std;
using namespace vcg;

class MyFace;
class MyEdge;
class MyVertex;

struct MyUsedTypes : public UsedTypes<	
    Use<MyVertex>::AsVertexType,
    Use<MyEdge>  ::AsEdgeType,
    Use<MyFace>  ::AsFaceType>{};


class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::BitFlags, vertex::Normal3f, vertex::Qualityf, vertex::Mark>{};
class MyEdge    : public Edge  < MyUsedTypes, edge::VertexRef, edge::EVAdj> {};
class MyFace    : public Face  < MyUsedTypes, face::VertexRef, face::FFAdj, face::BitFlags, face::Normal3f> {};

class MyEdgeMesh: public tri::TriMesh< vector<MyVertex>, vector<MyEdge> > {};
class MyMesh:     public tri::TriMesh< vector<MyVertex>, vector<MyFace > >{};

int main(int ,char **)
{  
  MyMesh m, m_over, m_under;;
  Sphere(m);
  printf("Created a sphere mesh of %i vertices %i edges\n",m.VN(),m.FN());
  
  // Compute intersection with a random plane
  math::MarsenneTwisterRNG rnd;
  Point3f direction = vcg::math::GeneratePointOnUnitSphereUniform<float,math::MarsenneTwisterRNG>(rnd);
  float distance = rnd.generate01();
  vcg::Plane3<MyMesh::ScalarType> plane(distance, direction);
  
  printf("Intersecting a sphere with a plane %4.2f %4.2f %4.2f off: %4.2f\n",direction[0],direction[1],direction[2],distance);
  MyEdgeMesh edge_mesh;  // the cross-section
  vcg::IntersectionPlaneMesh<MyMesh, MyEdgeMesh, MyMesh::ScalarType>(m, plane, edge_mesh);
  
  // Compute bounding box
  vcg::tri::UpdateBounding<MyEdgeMesh>::Box(edge_mesh);
  printf("Created a edge mesh of %i vertices %i edges\n",edge_mesh.VN(),edge_mesh.EN());
  
  // export the cross-section (projected along the plane normal direction) 
  tri::io::SVGProperties pro;
  pro.projDir = plane.Direction(); 
  tri::io::ExporterSVG<MyEdgeMesh>::Save(edge_mesh, "trimesh_intersection.svg",pro);
  
  // Now actually cut the mesh using the refine framework 
  tri::UpdateQuality<MyMesh>::VertexFromPlane(m, plane);
  tri::QualityMidPointFunctor<MyMesh> slicingfunc(0.0);
  tri::QualityEdgePredicate<MyMesh> slicingpred(0.0,0.0);
  tri::UpdateTopology<MyMesh>::FaceFace(m);
  tri::RefineE<MyMesh, tri::QualityMidPointFunctor<MyMesh>, tri::QualityEdgePredicate<MyMesh> > (m, slicingfunc, slicingpred, false);
  
  tri::UpdateSelection<MyMesh>::VertexFromQualityRange(m,0,std::numeric_limits<float>::max());
  tri::UpdateSelection<MyMesh>::FaceFromVertexStrict(m);
  tri::Append<MyMesh,MyMesh>::Mesh(m_over,m,true);
  tri::UpdateSelection<MyMesh>::FaceInvert(m);
  tri::Append<MyMesh,MyMesh>::Mesh(m_under,m,true);
  printf("Created a sphere mesh of %i vertices %i edges\n",m_over.VN(),m_over.FN());
  tri::io::ExporterOFF<MyMesh>::Save(m_over,"trimesh_intersection_over.off");
  printf("Created a sphere mesh of %i vertices %i edges\n",m_under.VN(),m_under.FN());
  tri::io::ExporterOFF<MyMesh>::Save(m_under,"trimesh_intersection_under.off");
  
  return 0;
}

