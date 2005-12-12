/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    *
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
/****************************************************************************
  History

$Log: not supported by cvs2svn $
Revision 1.8  2004/10/28 00:54:34  cignoni
Better Doxygen documentation

Revision 1.7  2004/09/01 12:17:25  pietroni
minor changes to comply gcc compiler (typename's )

Revision 1.6  2004/06/25 11:25:07  pietroni
changrd update for a vertex normal functions void PerVertex(VertexType *v)

Revision 1.5  2004/06/15 16:01:26  pietroni
added functions to update local vertices normal

Revision 1.4  2004/05/26 11:59:09  turini
Changed : Use Of vector In Use Of std::vector.

Revision 1.3  2004/05/26 11:53:17  pietroni
modified face declaration

Revision 1.2  2004/05/26 11:48:53  turini
Changed : include Path.

Revision 1.1  2004/05/12 13:52:11  pietroni
*** empty log message ***

Revision 1.2  2004/03/12 15:22:19  pietroni
Written some documentation and added to the trimes doxygen module


****************************************************************************/
#ifndef __VCG_TETRA_UPDATE_NORMALS
#define __VCG_TETRA_UPDATE_NORMALS

#include<vcg\complex\tetramesh\update\triconvert.h>
#include<vcg\simplex\face\face.h>
#include<vcg\complex\trimesh\base.h>
#include<vcg\complex\trimesh\update\normal.h>
#include<vcg\simplex\tetrahedron\pos.h>
#include<vector>

namespace vcg {
namespace tetra {

/** \addtogroup tetramesh */
/*@{*/

/// Management, updating and computation of per-vertex and per-face normals.
/// This class is used to compute or update the normals that can be stored in the vertex or face component of a mesh.
template <class ComputeMeshType>
class UpdateNormals
{

public:
typedef ComputeMeshType TetraMeshType; 
typedef typename TetraMeshType::VertexType     VertexType;
typedef typename TetraMeshType::VertexPointer  VertexPointer;
typedef typename TetraMeshType::VertexIterator VertexIterator;
typedef typename TetraMeshType::TetraType      TetraType;
typedef typename TetraMeshType::TetraPointer   TetraPointer;
typedef typename TetraMeshType::TetraIterator  TetraIterator;
typedef typename VertexType::NormalType		   NormalType;

typedef vcg::Face<VertexType,vcg::DUMMYEDGETYPE,vcg::DUMMYFACETYPE> FaceTemp;
typedef vcg::tri::TriMesh< std::vector<VertexType>,std::vector<FaceTemp> > TriMeshTemp;

/// Calculates the vertex normal (if stored in the current face type)
static void PerTetraFace(TetraMeshType &m)
{
	if( !m.HasPerTetraNormal()) return;
	TetraIterator t;
	for(t=m.tetra.begin();t!=m.tetra.end();++t)
			if( !(*t).IsD() )	(*t).ComputeNormal();
}

/// Calculates the vertex normal of a vertex.
static void PerVertex(VertexType *v)
{
 
 if( !VertexType::HasNormal()) return;

 VTIterator<TetraType> VTi=VTIterator<TetraType>(v->VTb(),v->VTi());
 NormalType  Norm=NormalType(0,0,0);
 int iter=0;
 while (!VTi.End())
 {
   //take the tree faces on the vertex's tetrahedron
   for (int j=0;j<3;j++)
   {
    int f=Tetra::FofV(VTi.Vi(),j);
    if (VTi.Vt()->IsBorderF(f))
      {
        iter++; 
        Norm+=VTi.Vt()->N(f);
      }
   }
  ++VTi;
 }
 Norm/=(float)iter;
 v->N()=Norm.Normalize();
}

/// Calculates the vertex normal. Without exploiting or touching face normals
/// The normal of a vertex v is the weigthed average of the normals of the faces incident on v.
static void PerVertex(TetraMeshType &m)
{
 
 if( !m.HasPerVertexNormal()) return;
 _ClearNormal(m);
 TriMeshTemp tri_mesh=TriMeshTemp();
 TriConverter <TetraMeshType,TriMeshTemp>tric=TriConverter<TetraMeshType,TriMeshTemp>();
 tric.Convert(m.tetra,tri_mesh);
 vcg::tri::UpdateNormals<TriMeshTemp> UNT=vcg::tri::UpdateNormals<TriMeshTemp>();
 UNT.PerVertexNormalized(tri_mesh);
}
private:
static void _ClearNormal(TetraMeshType &m)
{
 if( !m.HasPerVertexNormal()) return;
 VertexIterator vi;
 for(vi=m.vert.begin();vi!=m.vert.end();++vi)
   if( !(*vi).IsD() && (*vi).IsRW() )
     (*vi).N() = VertexType::NormalType(0,0,0);
}

///// Calculates both vertex and face normals.
///// The normal of a vertex v is the weigthed average of the normals of the faces incident on v.
//static void PerVertexPerFace(ComputeTetraMeshType &m)
//{
// if( !m.HasPerVertexNormal() || !m.HasPerFaceNormal()) return;
// 
// 
//}
//
//
//static void PerFaceNormalized(ComputeTetraMeshType &m)
//{
//	
//}
//
//
///// Calculates the vertex normal
//static void PerVertexNormalized(ComputeTetraMeshType &m)
//{
//  
//}


}; // end class

}	// End namespace
}	// End namespace


#endif
