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
Revision 1.5  2004/09/20 08:37:57  cignoni
Better Doxygen docs

Revision 1.4  2004/08/25 15:15:26  ganovelli
minor changes to comply gcc compiler (typename's and stuff)

Revision 1.3  2004/07/18 06:55:37  cignoni
NewUserBit -> NewBitFlag

Revision 1.2  2004/07/09 15:48:37  tarini
Added an include (<algorithm>)

Revision 1.1  2004/06/24 08:03:59  cignoni
Initial Release


****************************************************************************/

#ifndef __VCGLIB_CLEAN
#define __VCGLIB_CLEAN

#include <map>
#include <algorithm>

#include<vcg/complex/trimesh/allocate.h>
namespace vcg {
namespace tri {
/// 
/** \addtogroup trimesh */
/*@{*/
/// Class of static functions to clean/correct/restore meshs. 
template <class CleanMeshType>
class Clean
{
  public:
  typedef CleanMeshType MeshType; 
  typedef typename MeshType::VertexType     VertexType;
  typedef typename MeshType::VertexPointer  VertexPointer;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::FaceType       FaceType;
  typedef typename MeshType::FacePointer    FacePointer;
  typedef typename MeshType::FaceIterator   FaceIterator;
/* classe di confronto per l'algoritmo di eliminazione vertici duplicati*/
template <class VertexIterator>
class RemoveDuplicateVert_Compare{
public:
	inline bool operator() (VertexIterator a, VertexIterator b)
		{
			return *a < *b;
		}
};

/** This function removes all duplicate vertices of the mesh by looking only at their spatial positions. 
 Note that it does not update any topology relation that could be affected by this like the VT or TT relation.
 the reason this function is usually performed BEFORE building any topology information.
*/
static int RemoveDuplicateVertex( MeshType & m )    // V1.0
{
	if(m.vert.size()==0 || m.vn==0) return 0;

	std::map<VertexPointer, VertexPointer> mp;
	int i,j;
	VertexIterator vi; 
	int deleted=0;
	int k=0;
	int num_vert = m.vert.size();
	std::vector<VertexPointer> perm(num_vert);
	for(vi=m.vert.begin(); vi!=m.vert.end(); ++vi, ++k)
		perm[k] = &(*vi);

	RemoveDuplicateVert_Compare<VertexPointer> c_obj;

	std::sort(perm.begin(),perm.end(),c_obj);

  j = 0;
  i = j;
  mp[perm[i]] = perm[j];
  ++i;
  for(;i!=num_vert;)
	{
		if( (! (*perm[i]).IsD()) && 
        (! (*perm[j]).IsD()) && 
				(*perm[i]).P() == (*perm[j]).cP() )
		{
			VertexPointer t = perm[i];
	    mp[perm[i]] = perm[j];
	    ++i;
			(*t).SetD();
			deleted++;
		}
		else
		{
			j = i;
	    ++i;
		}
	}
  FaceIterator fi;
  for(fi = m.face.begin(); fi!=m.face.end(); ++fi)
		for(k = 0; k < 3; ++k)
			if( !(*fi).IsD() )
				if( mp.find( (typename MeshType::VertexPointer)(*fi).V(k) ) != mp.end() )
				{
					(*fi).V(k) = &*mp[ (*fi).V(k) ];
				}
	m.vn -= deleted;
	return deleted;
}


/** This function removes that are not referenced by any face. The function updates the vn counter.
		@param m The mesh
		@return The number of removed vertices
*/
static int RemoveUnreferencedVertex( CleanMeshType& m )   // V1.0
{
	FaceIterator fi;
	VertexIterator vi;
	int referredBit = VertexType::NewBitFlag();
		
	int j;
	int deleted = 0;
	
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
		(*vi).ClearUserBit(referredBit);

	for(fi=m.face.begin();fi!=m.face.end();++fi)
		if( !(*fi).IsD() )
			for(j=0;j<3;++j)
					(*fi).V(j)->SetUserBit(referredBit);

	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
		if( (!(*vi).IsD()) && (!(*vi).IsUserBit(referredBit)))
		{
			(*vi).SetD();
			++deleted;
		}
	m.vn -= deleted;
  VertexType::DeleteBitFlag(referredBit);
	return deleted;
}


}; // end class
/*@}*/
} // End Namespace TriMesh
} // End Namespace vcg
#endif
