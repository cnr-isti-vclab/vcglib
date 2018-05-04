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


#include<vcg/complex/algorithms/cut_tree.h>
#include<vcg/complex/algorithms/curve_on_manifold.h>
#include<vcg/complex/algorithms/crease_cut.h>

using namespace vcg;

class MyEdge;
class MyFace;
class MyVertex;
struct MyUsedTypes : public UsedTypes<Use<MyVertex>::AsVertexType, Use<MyEdge>::AsEdgeType, Use<MyFace>::AsFaceType>{};

class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::Qualityf, vertex::Color4b, vertex::VEAdj, vertex::VFAdj,vertex::BitFlags  >{};
class MyEdge    : public Edge<   MyUsedTypes, edge::VertexRef, edge::VEAdj,     edge::EEAdj, edge::BitFlags> {};
class MyFace    : public Face  < MyUsedTypes, face::VertexRef,  face::Normal3f, face::Qualityf, face::Color4b, face::VFAdj, face::FFAdj, face::Mark, face::Color4b, face::BitFlags > {};
class MyMesh : public tri::TriMesh< std::vector<MyVertex>, std::vector<MyEdge>, std::vector<MyFace> >{};

/**
 * In this sample we take a torus we compute a poly line on it that open it into a disk and we open it. 
 * Then using the COM (Curve On Manifold) framework we smooth this polyline keeping 
 * it on the surface of the torus and then first we refine the torus surface with this 
 * smooth polyline and then we open it along these new edges. 
 * 
 * Optionally you can use your own mesh and polyline by passing them as parameters. 
 */
int main(int argc,char ** argv )
{
  MyMesh base, basecopy, poly;
  int ret0=0, ret1=0;
  if(argc>1) ret0 = tri::io::Importer<MyMesh>::Open(base,argv[1]);
  
  if(base.FN() == 0)   Torus(base,10,4,48,24);

  if(argc>2) ret1 = tri::io::Importer<MyMesh>::Open(poly,argv[2]); 
  tri::UpdateBounding<MyMesh>::Box(base);
  printf( "Mesh %s has %i vert and %i faces\n", argv[1], base.VN(), base.FN() );
  printf( "Poly %s has %i vert and %i edges\n", argv[2], poly.VN(), poly.EN() );
  if(poly.EN() == 0) {
    srand(time(nullptr));
    tri::CutTree<MyMesh> ct(base);
    ct.Build(poly,rand()%base.fn);
  }
  tri::io::ExporterPLY<MyMesh>::Save(poly,"0_cut_tree.ply",tri::io::Mask::IOM_EDGEINDEX);  

  tri::CoM<MyMesh> cc(base);
  cc.Init();
  bool ret = cc.TagFaceEdgeSelWithPolyLine(poly);
  if(ret)
  {
    tri::Append<MyMesh,MyMesh>::MeshCopy(basecopy,base);  
    tri::UpdateTopology<MyMesh>::FaceFace(basecopy);
    tri::CutMeshAlongSelectedFaceEdges<MyMesh>(basecopy);
    tri::io::ExporterPLY<MyMesh>::Save(basecopy,"base_cut_with_tree.ply");  
  }
  // Selected vertices are 'locked' during the smoothing. 
  cc.SelectBoundaryVertex(poly);
//  cc.SelectUniformlyDistributed(poly,10); // lock some vertices uniformly just for fun
  
  // Two smoothing runs,
  // the first that allows fast movement over the surface (long edges that can skim surface details)
  cc.par.surfDistThr = base.bbox.Diag()/100.0f;
  cc.par.maxSimpEdgeLen = base.bbox.Diag()/50.0f;
  cc.par.minRefEdgeLen = base.bbox.Diag()/100.0f;
  cc.SmoothProject(poly,30,0.7f,.3f);
  tri::io::ExporterPLY<MyMesh>::Save(poly,"1_poly_smooth.ply",tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);

  // The second smooting run more accurate to adapt to the surface
  cc.par.surfDistThr = base.bbox.Diag()/1000.0f;
  cc.par.maxSimpEdgeLen = base.bbox.Diag()/1000.0f;
  cc.par.minRefEdgeLen = base.bbox.Diag()/2000.0f;
  cc.SmoothProject(poly,10,0.01f,.99f);
  tri::io::ExporterPLY<MyMesh>::Save(poly,"2_poly_smooth.ply",tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);

  Distribution<float> dist;  
  cc.EvaluateHausdorffDistance(poly, dist );
  
  // Adapt the polyline to the mesh (in the end it will have vertices only on edges and vertices of the base mesh)
  cc.RefineCurveByBaseMesh(poly);
  tri::io::ExporterPLY<MyMesh>::Save(poly,"3_poly_refined.ply",tri::io::Mask::IOM_EDGEINDEX+tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);   
  // Safely split the mesh with this refined polyline
  cc.SplitMeshWithPolyline(poly);
  tri::io::ExporterPLY<MyMesh>::Save(base,"3_mesh_refined.ply",tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);
  // Now the two meshes should have coincident edges
  cc.TagFaceEdgeSelWithPolyLine(poly);
  CutMeshAlongSelectedFaceEdges(base);
  tri::io::ExporterPLY<MyMesh>::Save(base,"4_mesh_cut.ply",tri::io::Mask::IOM_VERTCOLOR+tri::io::Mask::IOM_VERTQUALITY);
  
  return 0;
}

