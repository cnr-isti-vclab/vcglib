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


****************************************************************************/

#ifndef __VCGLIB_VERTEXALLOCATOR
#define __VCGLIB_VERTEXALLOCATOR

namespace vcg {
namespace vrt {
/** \addtogroup vertexmesh */
/*@{*/
/// Class to safely add vertexes and faces to a mesh updating all the involved pointers.
/// It provides static memeber to add either vertex or faces to a edgemesh.
template <class AllocateMeshType>
class Allocator
{
 
public:
typedef AllocateMeshType MeshType; 
typedef typename MeshType::VertexType     VertexType;
typedef typename MeshType::VertexPointer  VertexPointer;
typedef typename MeshType::VertexIterator VertexIterator;

/** This class is used when allocating new vertexes and faces to update 
   the pointers that can be changed when resizing the involved vectors of vertex or faces.
   It can also be used to prevent any update of the various mesh fields
   (e.g. in case you are building all the connections by hand as in a importer);
*/ 
template<class SimplexPointerType>
class PointerUpdater
{
public:
  void Clear(){newBase=oldBase=newEnd=oldEnd=0;preventUpdateFlag=false;};
  void Update(SimplexPointerType &vp)
  {
    vp=newBase+(vp-oldBase);
  }
  bool NeedUpdate() {if(newBase!=oldBase && !preventUpdateFlag) return true; else return false;}
  
  SimplexPointerType oldBase;
  SimplexPointerType newBase;
  SimplexPointerType newEnd;
  SimplexPointerType oldEnd;
  bool preventUpdateFlag; /// when true no update is considered necessary.
};


/** Function to safely add n vertices to a mesh. 

	@param m The mesh to be expanded
	@param n the number of vertexes to be added
	@param pu A PointerUpdater that stores the relocation that can be happened.
*/
static VertexIterator AddVertices(MeshType &m,int n, PointerUpdater<VertexPointer> &pu)
{
  VertexIterator last=m.vert.end();
  pu.Clear();
	if(m.vert.empty()) pu.oldBase=0;  // if the vector is empty we cannot find the last valid element
	else               pu.oldBase=&*m.vert.begin(); 
    
	for(int i=0; i<n; ++i)
	{
    m.vert.push_back(MeshType::VertexType());
		m.vert.back().ClearFlags();
	}

	m.vn+=n;
	return last;// iterator to the first added vertex
}

static VertexIterator AddVertices(MeshType &m, int n)
{
    PointerUpdater<VertexPointer> pu;
    return AddVertices(m, n,pu);
}
}; // end class
/*@}*/
} // End Namespace TriMesh
} // End Namespace vcg

#endif
