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
Revision 1.8  2005/11/21 21:44:43  cignoni
Moved ComputeNormal and ComputeNormalizedNormal out of the face class (no more a member function!)

Revision 1.7  2005/10/13 08:38:00  cignoni
removed the access to the face member function normal and substituted with vcg::normal(*f);

Revision 1.6  2005/06/17 00:46:09  cignoni
Added a PerVertexNormalizedPerFace (vertex are face/area weighted AND normalized)

Revision 1.5  2005/04/01 13:04:55  fiorin
Minor changes

Revision 1.4  2004/09/09 14:35:14  ponchio
Typename changes for linux

Revision 1.3  2004/08/31 15:18:54  pietroni
minor changes to comply gcc compiler (typename's )

Revision 1.2  2004/03/12 15:22:19  cignoni
Written some documentation and added to the trimes doxygen module

Revision 1.1  2004/03/05 10:59:24  cignoni
Changed name from plural to singular (normals->normal)

Revision 1.1  2004/03/04 00:05:50  cignoni
First working version!

Revision 1.1  2004/02/19 13:11:06  cignoni
Initial commit


****************************************************************************/
#ifndef __VCG_TRI_UPDATE_NORMALS
#define __VCG_TRI_UPDATE_NORMALS

namespace vcg {
namespace tri {

/** \addtogroup trimesh */
/*@{*/

/// Management, updating and computation of per-vertex and per-face normals.
/// This class is used to compute or update the normals that can be stored in the vertex or face component of a mesh.
template <class ComputeMeshType>
class UpdateNormals
{

public:
typedef ComputeMeshType MeshType; 	
typedef typename MeshType::VertexType     VertexType;
typedef typename VertexType::NormalType     NormalType;
typedef typename VertexType::ScalarType ScalarType;
typedef typename MeshType::VertexPointer  VertexPointer;
typedef typename MeshType::VertexIterator VertexIterator;
typedef typename MeshType::FaceType       FaceType;
typedef typename MeshType::FacePointer    FacePointer;
typedef typename MeshType::FaceIterator   FaceIterator;

/// Calculates the vertex normal (if stored in the current face type)
static void PerFace(ComputeMeshType &m)
{
	if( !m.HasPerFaceNormal()) return;
	FaceIterator f;
	for(f=m.face.begin();f!=m.face.end();++f)
    if( !(*f).IsD() )	face::ComputeNormal(*f);
}


/// Calculates the vertex normal. Without exploiting or touching face normals
/// The normal of a vertex v is the weigthed average of the normals of the faces incident on v.
static void PerVertex(ComputeMeshType &m)
{
 if( !m.HasPerVertexNormal()) return;
 
 VertexIterator vi;
 for(vi=m.vert.begin();vi!=m.vert.end();++vi)
   if( !(*vi).IsD() && (*vi).IsRW() )
     (*vi).N() = NormalType((ScalarType)0,(ScalarType)0,(ScalarType)0);

 FaceIterator f;

 for(f=m.face.begin();f!=m.face.end();++f)
   if( !(*f).IsD() && (*f).IsR() )
   {
    //typename FaceType::NormalType t = (*f).Normal();
    typename FaceType::NormalType t = vcg::Normal(*f);
 
    for(int j=0; j<3; ++j)
     if( !(*f).V(j)->IsD() && (*f).V(j)->IsRW() )  
      (*f).V(j)->N() += t;
   }
}

/// Calculates both vertex and face normals.
/// The normal of a vertex v is the weigthed average of the normals of the faces incident on v.
static void PerVertexPerFace(ComputeMeshType &m)
{
 if( !m.HasPerVertexNormal() || !m.HasPerFaceNormal()) return;
 
 PerFace(m);
 VertexIterator vi;
 for(vi=m.vert.begin();vi!=m.vert.end();++vi)
   if( !(*vi).IsD() && (*vi).IsRW() )
     (*vi).N() = NormalType((ScalarType)0,(ScalarType)0,(ScalarType)0);

 FaceIterator f;

 for(f=m.face.begin();f!=m.face.end();++f)
   if( !(*f).IsD() && (*f).IsR() )
   {
     for(int j=0; j<3; ++j)
     if( !(*f).V(j)->IsD() && (*f).V(j)->IsRW() )  
      (*f).V(j)->N() += (*f).cN();
   }
}

/// Calculates both vertex and face normals.
/// The normal of a vertex v is the weigthed average of the normals of the faces incident on v.
static void PerVertexNormalizedPerFace(ComputeMeshType &m)
{
PerVertexPerFace(m);
 VertexIterator vi;
 for(vi=m.vert.begin();vi!=m.vert.end();++vi)
   if( !(*vi).IsD() && (*vi).IsRW() ) 
     (*vi).N().Normalize();
}

static void PerFaceRW(ComputeMeshType &m, bool normalize=false)
{
	if( !m.HasPerFaceNormal()) return;

	FaceIterator f;
	bool cn = true;

	if(normalize)
	{
		for(f=m.m.face.begin();f!=m.m.face.end();++f)
		if( !(*f).IsD() && (*f).IsRW() )
		{
			for(int j=0; j<3; ++j)
				if( !(*f).V(j)->IsR()) 	cn = false;
      if( cn ) face::ComputeNormalizedNormal(*f);
			cn = true;
		}
	}
	else
	{
		for(f=m.m.face.begin();f!=m.m.face.end();++f)
			if( !(*f).IsD() && (*f).IsRW() )
			{
				for(int j=0; j<3; ++j)
					if( !(*f).V(j)->IsR()) 	cn = false;

				if( cn )
					(*f).ComputeNormal();
				cn = true;
			}
	}
}


static void PerFaceNormalized(ComputeMeshType &m)
{
	if( !m.HasPerFaceNormal()) return;
	FaceIterator f;
		for(f=m.face.begin();f!=m.face.end();++f)
      if( !(*f).IsD() )	face::ComputeNormalizedNormal(*f);
}


/// Calculates the vertex normal
static void PerVertexNormalized(ComputeMeshType &m)
{
  if( !m.HasPerVertexNormal()) return;
  PerVertex(m);
  for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
   if( !(*vi).IsD() && (*vi).IsRW() )
     (*vi).N().Normalize();
}

}; // end class

}	// End namespace
}	// End namespace


#endif
