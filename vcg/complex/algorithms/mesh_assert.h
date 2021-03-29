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
#ifndef __VCGLIB_MESH_ASSERT
#define __VCGLIB_MESH_ASSERT

#include <vcg/complex/exception.h>
#include <vcg/simplex/edge/pos.h>
#include <vcg/complex/algorithms/update/flag.h>

namespace vcg {
namespace tri {
/**
 * \brief For checking the adequacy of a mesh to a given algorithm.
 *
 * While the many RequireXXX functions allow to check the static correctness of a mesh and
 * have a O(1) complexity, in many cases we need to run more complex checks to be sure that
 * the subsequent algorithm can run without issues.
 * Typical cases are the fact that there are no unreferenced vertices (NoUnreferencedVertex)
 * or a given adjacency is correctly initialized (and not only statically present as a type component).
 * 
 * Naming Notes: 
 * - First two letters of the function name code the required adjacency involved 
 * - The exception text completes the sentence "This exception is thronw when/because...
 * 
 */
template <class MeshType>
class MeshAssert
{
public:
  typedef typename MeshType::VertexType VertexType;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::FaceType   FaceType;
  typedef typename MeshType::FaceIterator   FaceIterator;
  typedef typename MeshType::CoordType  CoordType;
  typedef typename MeshType::ScalarType ScalarType;

  /// \brief Throw vcg::MissingPreconditionException if FF adjacency is not initialized
  static void FFAdjacencyIsInitialized(MeshType &m)
  {
    for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
    {
      if(!fi->IsD())
        for(int i=0;i<fi->VN();++i)
        {
          if(fi->FFp(i)==0)
            throw vcg::MissingPreconditionException("FF adjacency is not initialized");
        }
    }
  }

  /// \brief Throw vcg::MissingPreconditionException if VF adjacency is not initialized
  static void VFAdjacencyIsInitialized(MeshType &m)
  {
    for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi) if(!vi->IsD())
      {
        if(vi->VFp() == nullptr)
          throw vcg::MissingPreconditionException("VF adjacency is not initialized");
      }
  }

  /// \brief Throw vcg::MissingPreconditionException if EE adjacency is not initialized
  static void EEAdjacencyIsInitialized(MeshType &m)
  {
    for(auto ei=m.edge.begin();ei!=m.edge.end();++ei) if(!ei->IsD())
        {
          if(ei->EEp(0)==0)
            throw vcg::MissingPreconditionException("EE adjacency is not initialized");
        }
  }
  
  /// \brief Throw vcg::MissingPreconditionException if According to EE adjacency, the edge mesh is not 1-manifold (e.g there are more than 2 edges on a vertex)
  static void EEOneManifold(MeshType &m)
  {
    EEAdjacencyIsInitialized(m);
    for(auto ei=m.edge.begin();ei!=m.edge.end();++ei) if(!ei->IsD())
    {
     if(! edge::IsEdgeManifold(*ei,0) )
       throw vcg::MissingPreconditionException("According to EE adjacency, the edge mesh is not 1-manifold (e.g there are more than 2 edges on a vertex)");
    }
  }

  /// \brief Throw vcg::MissingPreconditionException if Vertex Normals are not unit lenght
  static void VertexNormalNormalized(MeshType &m)
  {
    for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi) if(!vi->IsD())
      {
        if(fabs(vi->cN().Norm()-1.0)>0.000001)
          throw vcg::MissingPreconditionException("Vertex Normal are not normalized");
      }
  }

  
  
  /// \brief Throw vcg::MissingPreconditionException if There are unreferenced vertices
  static void NoUnreferencedVertex(MeshType &m)
  {
    UpdateFlags<MeshType>::VertexClearV(m);
    for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi) if(!fi->IsD())
    {
      for(int i=0;i<fi->VN();++i) fi->V(i)->SetV();
    }

    for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
      if(!vi->IsD())
      {
        if(!vi->IsV())
          throw vcg::MissingPreconditionException("There are unreferenced vertices");
      }
  }

  /// \brief Throw vcg::MissingPreconditionException if There are faces with more than three vertices
  static void OnlyTriFace(MeshType &m)
  {
    for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi) if(!fi->IsD())
    {
      if(fi->VN()!=3)
        throw vcg::MissingPreconditionException("There are faces with more than three vertices");
    }
  }

  /// \brief Throw vcg::MissingPreconditionException if There are non quadrilateral faces
  static void OnlyQuadFace(MeshType &m)
  {
    for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi) if(!fi->IsD())
    {
      if(fi->VN()!=4)
        throw vcg::MissingPreconditionException("There are non quadrilateral faces");
    }
  }

  /// \brief Throw vcg::MissingPreconditionException if The mesh is not composed only by edges (no faces needed or allowed)
  static void OnlyEdgeMesh(MeshType &m)
  {
      if(m.FN()>0)
        throw vcg::MissingPreconditionException("The mesh is not composed only by edges (no faces needed or allowed)");
  }
  
  /// \brief Throw vcg::MissingPreconditionException if According to FF adjacency, the mesh is not two manifold (e.g. there are more than two faces on an edge)
  static void FFTwoManifoldEdge(MeshType & m)
  {
    for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi) if(!fi->IsD())
    {
      for(int i=0;i<fi->VN();++i){
        if(!face::IsManifold(*fi,i))
          throw vcg::MissingPreconditionException("According to FF adjacency, the mesh is not two manifold (e.g. there are more than two faces on an edge)");
      }
    }    
}
  
};

} // end namespace tri
} // end namespace vcg
#endif // __VCGLIB_MESH_ASSERT
