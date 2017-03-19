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

#include<vcg/complex/complex.h>
#include<wrap/io_trimesh/import.h>
#include<wrap/io_trimesh/export.h>


#include<vcg/complex/algorithms/curve_on_manifold.h>
#include<vcg/complex/algorithms/crease_cut.h>

using namespace vcg;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<Use<MyVertex>::AsVertexType, Use<MyEdge>::AsEdgeType, Use<MyFace>::AsFaceType>{};

class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::Qualityf, vertex::Color4b, vertex::VEAdj, vertex::VFAdj,vertex::BitFlags  >{};
class MyEdge    : public Edge<   MyUsedTypes, edge::VertexRef, edge::VEAdj,     edge::EEAdj, edge::BitFlags> {};
class MyFace    : public Face  < MyUsedTypes, face::VertexRef,  face::Normal3f, face::VFAdj, face::FFAdj, face::Mark, face::Color4b, face::BitFlags > {};
class MyMesh : public tri::TriMesh< std::vector<MyVertex>, std::vector<MyEdge>, std::vector<MyFace> >{};


int main(int argc,char ** argv )
{
  MyMesh base,basecopy, poly;
  int ret0 = tri::io::Importer<MyMesh>::Open(base,argv[1]);
  tri::UpdateBounding<MyMesh>::Box(base);
  tri::Append<MyMesh,MyMesh>::MeshCopy(basecopy,base);
  printf( "Mesh %s has %i vert and %i faces\n", argv[1], basecopy.VN(), basecopy.FN() );
  if(ret0 != 0 ) 
  {
    printf("Failed Loading\n");
    exit(-1);
  }
  tri::CoM<MyMesh> cc(base);
  cc.Init();
  cc.BuildVisitTree(poly);
  tri::UpdateBounding<MyMesh>::Box(poly);  
  tri::io::ExporterPLY<MyMesh>::Save(poly,"tree.ply",tri::io::Mask::IOM_EDGEINDEX);  
  while(cc.OptimizeTree(poly));
  tri::io::ExporterPLY<MyMesh>::Save(poly,"tree1.ply",tri::io::Mask::IOM_EDGEINDEX);  

  cc.MarkFauxEdgeWithPolyLine(basecopy,poly);
  tri::UpdateTopology<MyMesh>::FaceFace(basecopy);
  tri::CutMeshAlongNonFauxEdges<MyMesh>(basecopy);
  tri::io::ExporterPLY<MyMesh>::Save(basecopy,"basecut.ply");  
  cc.par.surfDistThr = base.bbox.Diag()/100.0;
  cc.par.maxSimpEdgeLen = base.bbox.Diag()/60.0;
  cc.SmoothProject(poly,5,0.8, 0.2);
  Distribution<float> dist;  
  cc.EvaluateHausdorffDistance(poly, dist );
  tri::io::ExporterPLY<MyMesh>::Save(poly,"poly_adapted.ply",tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);
  cc.par.surfDistThr = base.bbox.Diag()/1000.0;
  cc.par.maxSimpEdgeLen = base.bbox.Diag()/500.0;
  cc.SmoothProject(poly,5,0.3, 0.7);
  tri::io::ExporterPLY<MyMesh>::Save(poly,"poly_adapted2.ply",tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);
  cc.par.minRefEdgeLen = base.bbox.Diag()/200.0;  
  cc.Refine(poly,true);
  cc.Refine(poly,true);
  cc.SnapPolyline(poly);
  tri::io::ExporterPLY<MyMesh>::Save(poly,"poly_adapted2_snap.ply",tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);
  cc.par.maxSimpEdgeLen = base.bbox.Diag()/200.0;
  cc.SmoothProject(poly,10 ,0.8, 0.2);
  cc.RefineBaseMesh(poly);
  tri::io::ExporterPLY<MyMesh>::Save(poly,"poly_adapted2_snapSmooth.ply",tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);
  
  return 0;
}

