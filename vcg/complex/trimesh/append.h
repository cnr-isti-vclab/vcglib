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


#ifndef __VCGLIB_APPEND
#define __VCGLIB_APPEND

#include <vcg/complex/trimesh/update/flag.h>

namespace vcg {
namespace tri {

template<class MeshType>
int Index(MeshType &m, typename MeshType::VertexType &v) {return &v-&*m.vert.begin();}
template<class MeshType>
int Index(MeshType &m, typename MeshType::FaceType &f) {return &f-&*m.face.begin();}

template<class MeshType>
int Index(MeshType &m, typename MeshType::VertexPointer &vp) {return vp-&*m.vert.begin();}
template<class MeshType>
int Index(MeshType &m, typename MeshType::FacePointer &fp) {return fp-&*m.face.begin();}

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

static void ImportVertex(VertexLeft &vl, VertexRight &vr)
{
  vl.P().Import(vr.P());
  if(vl.HasColor()   && vl.HasColor()) vl.C()=vr.C();
  if(vl.HasQuality() && vl.HasQuality()) vl.Q()=vr.Q();
  if(vl.HasTexture() && vl.HasTexture()) vl.T()=vr.T();
}

static void ImportFace(MeshLeft &ml, MeshRight &mr, FaceLeft &fl, FaceRight &fr, std::vector<int> &remap)
{
  fl.V(0)=&ml.vert[remap[ Index(mr,fr.V(0))]];
  fl.V(1)=&ml.vert[remap[ Index(mr,fr.V(1))]];
  fl.V(2)=&ml.vert[remap[ Index(mr,fr.V(2))]];
  if(fl.HasFaceColor()   && fl.HasFaceColor()) fl.C()=fr.C();
  if(fl.HasFaceQuality() && fl.HasFaceQuality()) fl.Q()=fr.Q();
  if(fl.HasWedgeTexture() && fl.HasWedgeTexture()) 
  {
    fl.WT(0)=fr.WT(0);
    fl.WT(1)=fr.WT(1);
    fl.WT(2)=fr.WT(2);
  };
}

static void Mesh(MeshLeft& ml, MeshRight& mr, const bool selected = false)
{
 // remap[i] keep where the position of where the i-th vertex of meshright has landed in meshleft
  std::vector<int> remap(mr.vn,-1); 
 
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
              ImportVertex((*vp),*(*fi).V(i));
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
}


static void Subset(MeshLeft& ml, std::vector<FacePointerRight> & vfpr)
{
 // remap[i] keep where the position of where the i-th vertex of meshright has landed in meshleft
  std::vector<int> remap(mr.vn,-1);
 
 // first loop to find the referenced vertices and copy them preparing the remap vector
 vector<FacePointerRight>::iterator  fi;
 int FaceToAdd=0;
 for(fi=vfpr.begin();fi!=vfpr.end();++fi)
    if(!(*fi)->IsD())
        {
          FaceToAdd++;
          for(int i=0;i<3;++i)
          {
            int vind=Index(mr, *(**fi).V(i));
            if(remap[vind]==-1)
            {
              VertexIteratorLeft vp;
              vp=Allocator<MeshLeft>::AddVertices(ml,1);
              ImportVertex((*vp),*(**fi).V(i));
              remap[vind]=Index(ml,*vp);
            }
          }
        }
 // second loop copy the faces updating the vertex references
 FaceIteratorLeft fp=Allocator<MeshLeft>::AddFaces(ml,FaceToAdd); 
 for(fi=vfpr.begin();fi!=vfpr.end();++fi)
    if(!(*fi).IsD())
    {
        ImportFace(ml,mr,(*fp),(*fi),remap);
        ++fp;
    }
}


static void Selected(MeshLeft& ml, MeshRight& mr)
{
  Mesh(mr,mr,true);
}

}; // end of class Append





} // End Namespace TriMesh
} // End Namespace vcg


#endif


