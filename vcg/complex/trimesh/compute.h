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

template <class ComputeMeshType>
class TriMeshCompute
{
public:
static void FaceNormalRW(ComputeMeshType &m, bool normalize=false)
{
	if( !m.HasPerFaceNormal()) return;

	face_iterator f;
	bool cn = true;

	if(normalize)
	{
		for(f=m.face.begin();f!=m.face.end();++f)
		if( !(*f).IsD() && (*f).IsRW() )
		{
			for(int j=0; j<(*f).size(); ++j)
				if( !(*f).V(j)->IsR()) 	cn = false;
			if( cn ) (*f).ComputeNormalizedNormal();
			cn = true;
		}
	}
	else
	{
		for(f=m.face.begin();f!=m.face.end();++f)
			if( !(*f).IsD() && (*f).IsRW() )
			{
				for(int j=0; j<(*f).size(); ++j)
					if( !(*f).V(j)->IsR()) 	cn = false;

				if( cn )
					(*f).ComputeNormal();
				cn = true;
			}
	}
}


static void NormalizedFaceNormal(ComputeMeshType &m)
{
   ComputeFaceNormal(true);
}


/// Calculates the vertex normal
static void NormalizedVertexNormal(ComputeMeshType &m)
{
 VertexNormal(true);
}
/// Calculates the vertex normal
static void VertexNormal(ComputeMeshType &m, bool normalize=false)
{
 if( m.HasPerVertexNormal())
 {
  vertex_iterator vi;

  for(vi=vert.begin();vi!=vert.end();++vi)
   if( !(*vi).IsDeleted() && (*vi).IsRW() )
    (*vi).Normal() = vectorial_type(0,0,0);

  face_iterator f;

  for(f=face.begin();f!=face.end();++f)
   if( !(*f).IsDeleted() && (*f).IsR() )
   {
    vectorial_type t = (*f).Normal();

    for(int j=0; j<(*f).size(); ++j)
     if( !(*f).V(j)->IsD() && (*f).V(j)->IsR() && (*f).V(j)->IsW() )  
      (*f).V(j)->Normal() += t;
   }
 if(normalize)
  for(vi=vert.begin();vi!=vert.end();++vi)
   if( !(*vi).IsDeleted() && (*vi).IsRW() )
     (*vi).Normal().Normalize();
 }
}
void ComputeE()
{
	face_iterator f;
 
	for(f = face.begin(); f!=face.end(); ++f)
		(*f).ComputeE();
}

}; // end class
