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
Revision 1.2  2004/03/03 15:35:52  cignoni
Yet another cr lf mismatch

Revision 1.1  2004/02/24 21:36:42  cignoni
grouped documentation, changed typenames and reflection mechanism

Revision 1.1  2004/02/19 13:11:06  cignoni
Initial commit


****************************************************************************/

namespace vcg {
namespace tri {
/** \addtogroup trimesh */
/*@{*/
/// Class to safely add vertexes and faces to a mesh updating all the involved pointers.
/// It provides static memeber to add either vertex or faces to a trimesh.
template <class AllocateMeshType>
class Allocator
{
 
public:
typedef AllocateMeshType MeshType; 
typedef typename MeshType::VertexType     VertexType;
typedef typename MeshType::VertexPointer  VertexPointer;
typedef typename MeshType::VertexIterator VertexIterator;
typedef typename MeshType::FaceType       FaceType;
typedef typename MeshType::FacePointer    FacePointer;
typedef typename MeshType::FaceIterator   FaceIterator;

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


/** Function to add n vertices to the mesh. The second parameter hold a vector of 
	pointers to pointer to elements of the mesh that should be updated after a 
	possible vector realloc. 
	@param n Il numero di vertici che si vuole aggiungere alla mesh.
	@param local_var Vettore di variabili locali che rappresentano puntatori a vertici. 
	restituisce l'iteratore al primo elemento aggiunto.
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

	pu.newBase = &*m.vert.begin();
	if(pu.NeedUpdate())
		{
			FaceIterator fi;
			for (fi=m.face.begin(); fi!=m.face.end(); ++fi)
        if(!(*fi).IsD())
        {
          pu.Update((*fi).V(0));
          pu.Update((*fi).V(1));
          pu.Update((*fi).V(2));
        }
		
		// e poiche' lo spazio e' cambiato si ricalcola anche last da zero  
			unsigned int siz=m.vert.size()-n;	
		  if(last!=0)  
			{ 
				last = m.vert.begin(); 
				advance(last,siz);
			}
			else last=m.vert.begin(); 
		}
 
	return last;// deve restituire l'iteratore alla prima faccia aggiunta;
}

static VertexIterator AddVertices(MeshType &m, int n)
{
    PointerUpdater<VertexPointer> pu;
    return AddVertices(m, n,pu);
}

/** Function to add n faces to the mesh.
	@param n Il numero di facce che si vuole aggiungere alla mesh
*/
static FaceIterator AddFaces(MeshType &m, int n)
{
	PointerUpdater<FacePointer> pu;
  return AddFaces(m,n,pu);
}
/** Function to add n faces to the mesh. 
  NOTA: Aggiorna fn;
*/
static FaceIterator AddFaces(MeshType &m, int n, PointerUpdater<FacePointer> &pu)
{
  FaceIterator last=0;
  pu.Clear();
  if(m.face.empty()) {
    pu.oldBase=0;  // if the vector is empty we cannot find the last valid element
    last=0;
  }	else  {
    pu.oldBase=&*m.face.begin(); 
    last=m.face.end();
  } 
	for(int i=0; i<n; ++i)
	{
    m.face.push_back(MeshType::FaceType());
		m.face.back().ClearFlags();
	}

	m.fn+=n;
	
  pu.newBase = &*m.face.begin();

  if(pu.NeedUpdate())
		{
      FaceIterator fi;
			for (fi=m.face.begin(); fi!=m.face.end(); ++fi)
        if(!(*fi).IsD())
        {
          if(FaceType::HasFFAdjacency())
          {
            pu.Update((*fi).F(0));
            pu.Update((*fi).F(1));
            pu.Update((*fi).F(2));
          }
        }
      VertexIterator vi;
			for (vi=m.vert.begin(); vi!=m.vert.end(); ++vi)
        if(!(*vi).IsD())
        {
          if(VertexType::HasVFAdjacency())
          {
            pu.Update((*fi).F(0));
            pu.Update((*fi).F(1));
            pu.Update((*fi).F(2));
          }
        }
        		// e poiche' lo spazio e' cambiato si ricalcola anche last da zero  
		unsigned int siz=m.face.size()-n;	
		if(last!=0)  
			{ 
				last = m.face.begin(); 
				advance(last,siz);
			}
 		else last=m.face.begin(); 
		}

	return last;
}

}; // end class
/*@}*/
} // End Namespace TriMesh
} // End Namespace vcg
