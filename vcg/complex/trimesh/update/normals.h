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
Revision 1.1  2004/02/19 13:11:06  cignoni
Initial commit


****************************************************************************/
#ifndef __VCG_TRI_UPDATE_NORMALS
#define __VCG_TRI_UPDATE_NORMALS

namespace vcg {
namespace tri {

template <class ComputeMeshType>
class UpdateNormals
{

public:
typedef ComputeMeshType MeshType; 
typedef typename MeshType::VertexType     VertexType;
typedef typename MeshType::VertexPointer  VertexPointer;
typedef typename MeshType::VertexIterator VertexIterator;
typedef typename MeshType::FaceType       FaceType;
typedef typename MeshType::FacePointer    FacePointer;
typedef typename MeshType::FaceIterator   FaceIterator;

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
			for(int j=0; j<3); ++j)
				if( !(*f).V(j)->IsR()) 	cn = false;
			if( cn ) (*f).ComputeNormalizedNormal();
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
			if( !(*f).IsD() )	(*f).ComputeNormalizedNormal();
}

static void PerFace(ComputeMeshType &m)
{
	if( !m.HasPerFaceNormal()) return;
	FaceIterator f;
	for(f=m.face.begin();f!=m.face.end();++f)
			if( !(*f).IsD() )	(*f).ComputeNormal();
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

/// Calculates the vertex normal
static void PerVertex(ComputeMeshType &m)
{
 if( !m.HasPerVertexNormal()) return;
 
 VertexIterator vi;
 for(vi=m.vert.begin();vi!=m.vert.end();++vi)
   if( !(*vi).IsD() && (*vi).IsRW() )
     (*vi).N() = VertexType::NormalType(0,0,0);

 FaceIterator f;

 for(f=m.face.begin();f!=m.face.end();++f)
   if( !(*f).IsD() && (*f).IsR() )
   {
    FaceType::NormalType t = (*f).Normal();

    for(int j=0; j<3; ++j)
     if( !(*f).V(j)->IsD() && (*f).V(j)->IsRW() )  
      (*f).V(j)->N() += t;
   }
}



void ComputeE()
{
	FaceIterator f;
 
	for(f = m.face.begin(); f!=m.face.end(); ++f)
		(*f).ComputeE();
}

}; // end class

}	// End namespace
}	// End namespace


#endif
