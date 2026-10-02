/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2026                                           \/)\/    *
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

/*! \file trimesh_handle_tunnel.cpp
\ingroup code_sample

\brief Reeb graph, handle and tunnel loops of a surface

    trimesh_handle_tunnel [mesh.ply]

Without arguments it uses a torus: its handle loop is a meridian, around the tube
(length 2*pi*1 for the 32-gon here: 6.273), its tunnel loop is the inner equator, around
the hole (2*pi*2 for the 64-gon: 12.561). The Reeb graph and the loops are saved as edge
meshes, in reeb_graph.ply, handles.ply and tunnels.ply.
*/

#include <cstdio>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <vcg/complex/algorithms/handle_tunnel_loops.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>

using namespace vcg;

class MyVertex;
class MyEdge;
class MyFace;
struct MyUsedTypes : public UsedTypes<Use<MyVertex>::AsVertexType,
                                      Use<MyEdge>::AsEdgeType,
                                      Use<MyFace>::AsFaceType> {};

class MyVertex : public Vertex<MyUsedTypes, vertex::Coord3f, vertex::Qualityf, vertex::BitFlags> {};
class MyEdge   : public Edge<MyUsedTypes, edge::VertexRef, edge::BitFlags> {};
class MyFace   : public Face<MyUsedTypes, face::VertexRef, face::FFAdj, face::BitFlags> {};
class MyMesh   : public tri::TriMesh<std::vector<MyVertex>, std::vector<MyEdge>,
                                     std::vector<MyFace>> {};

static float Length(const MyMesh &m, const tri::HandleTunnelLoops<MyMesh>::Loop &l)
{
  float len = 0;
  for (size_t i = 0; i < l.size(); ++i) len += Distance(m.vert[l[i]].cP(), m.vert[l[(i + 1) % l.size()]].cP());
  return len;
}

int main(int argc, char **argv)
{
  MyMesh m;
  if (argc < 2) tri::Torus(m, 3, 1, 64, 32);
  else if (tri::io::Importer<MyMesh>::Open(m, argv[1]) != 0) return 1;
  tri::UpdateTopology<MyMesh>::FaceFace(m);

  // The Reeb graph of a height function has one independent cycle per handle.
  tri::UpdateQuality<MyMesh>::VertexFromPlane(m, Plane3f(0, Point3f(0.3f, 0.5f, 0.8f)));
  MyMesh graph;
  tri::ReebGraph<MyMesh> rg;
  rg.Compute(m, graph);
  std::printf("Reeb graph: %d nodes, %d arcs\n", graph.VN(), graph.EN());
  tri::io::ExporterPLY<MyMesh>::Save(graph, "reeb_graph.ply", tri::io::Mask::IOM_EDGEINDEX);

  tri::HandleTunnelLoops<MyMesh> ht;
  ht.Compute(m);
  std::printf("Genus %d\n", ht.genus);
  for (const auto &l : ht.handles) std::printf("  handle of length %f\n", Length(m, l));
  for (const auto &l : ht.tunnels) std::printf("  tunnel of length %f\n", Length(m, l));

  MyMesh hm, tm;
  tri::HandleTunnelLoops<MyMesh>::LoopsToEdgeMesh(m, ht.handles, hm);
  tri::HandleTunnelLoops<MyMesh>::LoopsToEdgeMesh(m, ht.tunnels, tm);
  tri::io::ExporterPLY<MyMesh>::Save(hm, "handles.ply", tri::io::Mask::IOM_EDGEINDEX);
  tri::io::ExporterPLY<MyMesh>::Save(tm, "tunnels.ply", tri::io::Mask::IOM_EDGEINDEX);
  return 0;
}
