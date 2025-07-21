/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2024                                           \/)\/    *
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

#pragma once

namespace vcg {
namespace tri {
/// \ingroup trimesh

/// \headerfile refine_doosabin.h vcg/complex/algorithms/refine_doosabin.h

/// \brief This class is used convert between polygonal meshes and triangular meshes

/**
    This class contains two members that allow to build a triangular mesh from a polygonal mesh
    and viceversa. In a trimesh, the generic polygons with n sides are codified represented by
    tagging the internal edge of the face as 'faux' with the SetF.
    */

template <class PolyMeshType>
class CatmullClark {
    typedef typename PolyMeshType::FaceType FaceType;
    typedef typename PolyMeshType::FacePointer FacePointer;
    typedef typename PolyMeshType::FaceIterator FaceIterator;
    typedef typename PolyMeshType::VertexIterator VertexIterator;
    
public:
static void Refine(PolyMeshType &baseIn, PolyMeshType &refinedOut, int iterationNum=1)
{
    PolyMeshType baseTmp,outTmp;

    tri::Append<PolyMeshType,PolyMeshType>::MeshCopy(baseTmp,baseIn);

    for(int i=0;i<iterationNum;++i)
    {
        RefineSingleStep(baseTmp, outTmp);
        tri::Append<PolyMeshType,PolyMeshType>::MeshCopy(baseTmp,outTmp);
    }
    tri::Append<PolyMeshType,PolyMeshType>::MeshCopy(refinedOut,outTmp);   
}


static void RefineSingleStep(PolyMeshType &baseIn, PolyMeshType &refinedOut)
{
	tri::RequirePolygonalMesh(baseIn);
	tri::RequirePolygonalMesh(refinedOut);
    refinedOut.Clear();
	//tri::RequireFFAdjacency(baseIn);
    // for the output mesh we keep a counter for each vertex as an additional attribute to computing averages easily
    auto v_cnH = tri::Allocator<PolyMeshType>::template AddPerVertexAttribute<int>(refinedOut,"counter");

    // for the input mesh we keep an attribute with the valence of the vertices
    auto v_valH = tri::Allocator<PolyMeshType>:: template AddPerVertexAttribute<int>(baseIn,"valence");

    // step -1 init the vertex valence counter to zero
    for(auto &v : baseIn.vert)
        v_valH[&v]=0;

  // Compute Valence for each vertex by iterating on all the faces
  for(auto &f : baseIn.face)
    for(int i=0;i<f.VN();i++)
      v_valH[f.V(i)]++;

  // step 0 copy all the vertexes in the output mesh
  for(auto &v : baseIn.vert)
    tri::Allocator<PolyMeshType>::AddVertex(refinedOut,Point3f(0,0,0));
    printf( "Mesh has %i vert and %i faces\n", refinedOut.VN(), refinedOut.FN() );

  // We keep two maps to store the new vertices created for the edges and the faces 
  // and retrieve them when actually building the connectivity of the new mesh
  std::map<std::pair<int,int>, int> edgeMap;
  std::map<int, int> faceMap;

  // First Step Create the mid face vertices
  for(auto &f : baseIn.face)
  {
    auto vp = tri::Allocator<PolyMeshType>::AddVertex(refinedOut,PolyBarycenter(f));
    faceMap[tri::Index(baseIn,f)] = tri::Index(refinedOut,*vp);
  }

  // Second step. Create a vertex for each edge.     
  // looping over all faces
  for(auto &f : baseIn.face)
  {
    Point3f center = refinedOut.vert[faceMap[tri::Index(baseIn,f)]].P();
    // loop over all edges
    for(int i=0;i<f.VN();i++)
    {
      int v0 = tri::Index(baseIn,f.V0(i));
      int v1 = tri::Index(baseIn,f.V1(i));
      if(v0>v1) std::swap(v0,v1);
       
      // check if the edge is already in the map
      auto it = edgeMap.find(std::make_pair(v0,v1));

      int edgeVertIndex;
      if(it == edgeMap.end())
      {
        // if not, create a new vertex in the middle of the edge
        auto vp = tri::Allocator<PolyMeshType>::AddVertex(refinedOut,Point3f(0,0,0));        
        edgeMap[std::make_pair(v0,v1)] = tri::Index(refinedOut,*vp);
        edgeVertIndex = tri::Index(refinedOut,*vp);
      }
      else
      { 
        edgeVertIndex = it->second;
      }

        refinedOut.vert[edgeVertIndex].P() += center*0.5;
        refinedOut.vert[edgeVertIndex].P() += f.V0(i)->P()*0.25;
        refinedOut.vert[edgeVertIndex].P() += f.V1(i)->P()*0.25;
        
        v_cnH[edgeVertIndex]+=1; // increment the counter of the vertex for the face
    }
  }
  // the only normalization step we do is for edge vertices to handle the fact 
  // that on boundary faces we do not have two contributions so the weight should be doubled.
  for(auto &v : refinedOut.vert)
    if(v_cnH[v]>0) v.P() /= v_cnH[&v]; 

  // Third step compute new the position of the original vertices
  for(auto &f : baseIn.face)
  {
    Point3f center = refinedOut.vert[faceMap[tri::Index(baseIn,f)]].P();

    // for each vertex of the face we add to its position 
    // the coords of the baricenter and of the two vertices adjacent to it
    for(int i=0;i<f.VN();i++)
    {
      float val = v_valH[f.V(i)];
      float val1 = (val -2.0)/(val*val);
      float val2 = 1.0/(val*val);

      int v0 = tri::Index(baseIn,f.V0(i));
      int v1 = tri::Index(baseIn,f.V1(i));
      if(v0>v1) std::swap(v0,v1);
      Point3f edgep = refinedOut.vert[edgeMap[std::make_pair(v0,v1)]].P();
      
      refinedOut.vert[tri::Index(baseIn,f.V(i))].P() += f.V(i)->P() * val1;
      refinedOut.vert[tri::Index(baseIn,f.V(i))].P() += center * val2;
      refinedOut.vert[tri::Index(baseIn,f.V(i))].P() += edgep * val2;
    }
  }
  
  printf( "Mesh has %i vert and %i faces\n", refinedOut.VN(), refinedOut.FN() );

  // Final step. Create Connectivity: a new face for each wedge of each face of the original mesh
  for(auto &f : baseIn.face)
  {
    // loop over all wedge of the face 
    for(int i=0;i<f.VN();i++)
    {
      int e00 = tri::Index(baseIn,f.V0(i));
      int e01 = tri::Index(baseIn,f.V1(i));
      int e10 = tri::Index(baseIn,f.V0(i));
      int e11 = tri::Index(baseIn,f.V((i+f.VN()-1)%f.VN()));
      if(e00>e01) std::swap(e00,e01);
      if(e10>e11) std::swap(e10,e11);

      // retrieve the edge vertices 
      auto e1it = edgeMap.find(std::make_pair(e00,e01));
      auto e2it = edgeMap.find(std::make_pair(e10,e11));

      // retrieve the face vertex
      auto fit = faceMap.find(tri::Index(baseIn,f));
      
      auto v0p = &refinedOut.vert[tri::Index(baseIn,f.V(i))];
      auto v1p = &refinedOut.vert[e1it->second];
      auto v2p = &refinedOut.vert[fit->second];
      auto v3p = &refinedOut.vert[e2it->second];

      tri::Allocator<PolyMeshType>::AddQuadFace(refinedOut, v0p, v1p, v2p, v3p);
    }
  }//	fprintf(stdout,"Refining starting \n");fflush(stdout);
	tri::Allocator<PolyMeshType>::template DeletePerVertexAttribute<int>(refinedOut,v_cnH);

    // for the input mesh we keep an attribute with the valence of the vertices
    tri::Allocator<PolyMeshType>:: template DeletePerVertexAttribute<int>(baseIn,v_valH);

} // end refine Function

}; // end CatmullClark class 
} // end namespace tri
} // end namespace vcg

