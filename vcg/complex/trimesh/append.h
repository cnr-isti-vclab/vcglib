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
Revision 1.7  2008/01/28 08:39:56  cignoni
added management of normals

Revision 1.6  2007/03/12 15:38:03  tarini
Texture coord name change!  "TCoord" and "Texture" are BAD. "TexCoord" is GOOD.

Revision 1.5  2006/05/25 04:40:57  cignoni
Updated HasPerFaceColor/Quality to the new style with mesh param.

Revision 1.4  2006/04/11 13:51:21  zifnab1974
commented out one function which does not compile on linux with gcc 3.4.5

Revision 1.3  2006/01/30 09:00:40  cignoni
Corrected use of HasPerWedgeTexture

Revision 1.2  2006/01/22 17:08:50  cignoni
Bug due to wrong compuation of size of auxiliary vector (vn instead of vert.size() )

Revision 1.1  2006/01/11 15:45:21  cignoni
Initial Release

****************************************************************************/


#ifndef __VCGLIB_APPEND
#define __VCGLIB_APPEND

#include <vcg/complex/trimesh/allocate.h>
#include <vcg/complex/trimesh/update/flag.h>

namespace vcg {
namespace tri {

template<class MeshLeft, class MeshRight> 
class Append
{
public:
 typedef typename MeshLeft::ScalarType  ScalarLeft;
 typedef typename MeshLeft::CoordType CoordLeft;
 typedef typename MeshLeft::VertexType  VertexLeft;
 typedef typename MeshLeft::FaceType  FaceLeft;
 typedef typename MeshLeft::VertexPointer  VertexPointerLeft;
 typedef typename MeshLeft::VertexIterator VertexIteratorLeft;
 typedef typename MeshLeft::FaceIterator   FaceIteratorLeft;

 typedef typename MeshRight::ScalarType ScalarRight;
 typedef typename MeshRight::CoordType CoordRight;
 typedef typename MeshRight::VertexType  VertexRight;
 typedef typename MeshRight::FaceType  FaceRight;
 typedef typename MeshRight::VertexPointer  VertexPointerRight;
 typedef typename MeshRight::VertexIterator VertexIteratorRight;
 typedef typename MeshRight::FaceIterator   FaceIteratorRight;
 typedef typename MeshRight::FacePointer   FacePointerRight;

static void ImportFace(MeshLeft &ml, MeshRight &mr, FaceLeft &fl, const FaceRight &fr, std::vector<int> &remap)
{
	fl.template ImportLocal<FaceRight>(fr);
  fl.V(0)=&ml.vert[remap[ Index(mr,fr.V(0))]];
  fl.V(1)=&ml.vert[remap[ Index(mr,fr.V(1))]];
  fl.V(2)=&ml.vert[remap[ Index(mr,fr.V(2))]];

	if(HasPerWedgeTexCoord(mr) && HasPerWedgeTexCoord(ml)) 
		for(int i=0;i<3;++i){
				fl.WT(i).P()=fr.cWT(i).P();
				fl.WT(i).N()=fr.cWT(i).N()+ml.textures.size();
		}
}

// Append Right Mesh to the Left Mesh
// Append::Mesh(ml, mr) is equivalent to ml += mr. 
// Note MeshRigth could be costant...
static void Mesh(MeshLeft& ml, MeshRight& mr, const bool selected = false, const bool copyUnrefFlag=false)
{
 // remap[i] keep where the position of where the i-th vertex of meshright has landed in meshleft
  std::vector<int> remap(mr.vert.size(),-1); 
 
 if(copyUnrefFlag) // copy ALL the vertices of MR onto ML
 {
	VertexIteratorRight vi;
	for(vi=mr.vert.begin();vi!=mr.vert.end();++vi)
	{
	 int vind=Index(mr,*vi);
		if(remap[vind]==-1)
			{
				VertexIteratorLeft vp;
				vp=Allocator<MeshLeft>::AddVertices(ml,1);
				(*vp).ImportLocal(*(vi));
				remap[vind]=Index(ml,*vp);
			}
	}
 } 
 // first loop to find the referenced vertices and copy them preparing the remap vector
 FaceIteratorRight fi;
 int FaceToAdd=0;
 for(fi=mr.face.begin();fi!=mr.face.end();++fi)
    if((!(*fi).IsD()) && (!selected || (*fi).IsS() ))
        {
          ++FaceToAdd;
          for(int i=0;i<3;++i)
          {
            int vind=Index(mr, *(*fi).V(i));
            if(remap[vind]==-1)
            {
              VertexIteratorLeft vp;
              vp=Allocator<MeshLeft>::AddVertices(ml,1);
              (*vp).ImportLocal(*(*fi).V(i));
              remap[vind]=Index(ml,*vp);
            }
          }
        }
 // second loop copy the faces updating the vertex references
 FaceIteratorLeft fp=Allocator<MeshLeft>::AddFaces(ml,FaceToAdd); 
 for(fi=mr.face.begin();fi!=mr.face.end();++fi)
 if(!(*fi).IsD() && (!selected || (*fi).IsS() ))
    {
        ImportFace(ml,mr,(*fp),(*fi),remap);
        ++fp;
    }
	// At the end concatenate the vector with texture names.
	 ml.textures.insert(ml.textures.end(),mr.textures.begin(),mr.textures.end());
}


// static void Subset(MeshLeft& ml, std::vector<FacePointerRight> & vfpr)
// {
//  // remap[i] keep where the position of where the i-th vertex of meshright has landed in meshleft
//   std::vector<int> remap(mr.vert.size(),-1); 
 
//  // first loop to find the referenced vertices and copy them preparing the remap vector
// typename  std::vector<FacePointerRight>::iterator  fi;
//  int FaceToAdd=0;
//  for(fi=vfpr.begin();fi!=vfpr.end();++fi)
//     if(!(*fi)->IsD())
//         {
//           FaceToAdd++;
//           for(int i=0;i<3;++i)
//           {
//             int vind=Index(mr, *(**fi).V(i));
//             if(remap[vind]==-1)
//             {
//               VertexIteratorLeft vp;
//               vp=Allocator<MeshLeft>::AddVertices(ml,1);
//               ImportVertex((*vp),*(**fi).V(i));
//               remap[vind]=Index(ml,*vp);
//             }
//           }
//         }
//  // second loop copy the faces updating the vertex references
//  FaceIteratorLeft fp=Allocator<MeshLeft>::AddFaces(ml,FaceToAdd); 
//  for(fi=vfpr.begin();fi!=vfpr.end();++fi)
//     if(!(*fi).IsD())
//     {
//         ImportFace(ml,mr,(*fp),(*fi),remap);
//         ++fp;
//     }
// }


static void Selected(MeshLeft& ml, MeshRight& mr)
{
  Mesh(ml,mr,true);
}

}; // end of class Append





} // End Namespace TriMesh
} // End Namespace vcg


#endif


