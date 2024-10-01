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
#ifndef __VCGLIB_IOTRIMESH_IO_PLY
#define __VCGLIB_IOTRIMESH_IO_PLY

/**
@name Ply I/O
It contains the definition of the accessory class PlyInfo
It is directly included by the importer and exporter ply classes
*/
//@{

#include "vcg/complex/allocate.h"
#include "vcg/complex/base.h"
#include "vcg/space/point3.h"
#include "wrap/callback.h"
#include "wrap/ply/plylib.h"

namespace vcg {
namespace tri {
namespace io {

/** Additional data needed or useful for parsing a ply mesh.
This class can be passed to the ImporterPLY::Open() function for
- retrieving additional per-vertex per-face data
- specifying a callback for long ply parsing
- knowing what data is  contained in a ply file
*/
class PlyInfo
{
public:
	typedef ::vcg::ply::PropDescriptor PropDescriptor ;
    typedef int PlyElemType;
    const static PlyElemType PlyElemVertex = 0;
    const static PlyElemType PlyElemFace = 1;
    
  /// Store the error codes enconutered when parsing a ply
  int status=0;
  /// It returns a bit mask describing the field present in the ply file
  int mask=0;

  /// a Simple callback that can be used for long ply parsing.
  // it returns the current position, and formats a string with a description of what th efunction is doing (loading vertexes, faces...)
  CallBackPos *cb=0;

  /// The additional vertex descriptor that a user can specify to load additional per-vertex non-standard data stored in a ply
  std::vector<PropDescriptor> VertDescriptorVec;
  
  /// VertAttrNameVec is a vector containing the names of the attributes to be saved (loaded).
  /// We assume that VertAttrNameVec, if not empty, is exactly of the same size of VertDescriptorVec
  /// If VertAttrNameVec[i] is not empty we use it to retrieve/store the info instead of the offsetted space in the current vertex
  std::vector<std::string> VertAttrNameVec;

  std::vector<PropDescriptor> FaceDescriptorVec;
  std::vector<std::string> FaceAttrNameVec;

  /// a string containing the current ply header. Useful for showing it to the user.
  std::string header;

enum Error
{
	// Funzioni superiori
	E_NO_VERTEX       = ply::E_MAXPLYERRORS+1,       // 15
	E_NO_FACE         = ply::E_MAXPLYERRORS+2,       // 16
	E_SHORTFILE       = ply::E_MAXPLYERRORS+3,       // 17
	E_NO_3VERTINFACE  = ply::E_MAXPLYERRORS+4,       // 18
	E_BAD_VERT_INDEX  = ply::E_MAXPLYERRORS+5,       // 19
	E_NO_6TCOORD      = ply::E_MAXPLYERRORS+6,		 // 20
	E_DIFFER_COLORS   = ply::E_MAXPLYERRORS+7,	     // 21
	E_BAD_VERT_INDEX_EDGE  = ply::E_MAXPLYERRORS+8,  // 22
	E_MAXPLYINFOERRORS= ply::E_MAXPLYERRORS+9        // 23
};
  
  /**
	 * @brief addPerElemScalarAttribute
	 * @param elemType  Vertex or Face
	 * @param propertyType the type of the property
	 * @param attrName the name of the attribute (both in the ply file and in the mesh)
	 * @param propName, optional only in case we want to load a certain ply property in an attribute with a different name
	 * 
	 * This function is used to specify additional per-vertex or per-face data 
	 * that is stored in the ply file and that should loaded or saved using attributes
	 * 
	 */ 
  void addPerElemScalarAttribute(PlyElemType elemType, vcg::ply::PlyTypes propertyType, const std::string& attrName, std::string propName="")
  {
      if(propName=="")
          propName=attrName;
      
      PropDescriptor p;
      p.propname=propName;
      p.islist = false;
      p.stotype1 = propertyType;
      p.memtype1 = propertyType;
      
      if (elemType == PlyElemVertex){
          p.elemname="vertex";
          VertDescriptorVec.push_back(p);
          VertAttrNameVec.resize(VertDescriptorVec.size());
          VertAttrNameVec.back() = attrName;          
      }
      else if (elemType == PlyElemFace){
          p.elemname="face";
          FaceDescriptorVec.push_back(p);
          FaceAttrNameVec.resize(FaceDescriptorVec.size());
          FaceAttrNameVec.back() = attrName;
      }
  }
  
  void AddPerElemFloatAttribute(PlyElemType elemType, const std::string& attrName, std::string propName="")
  {
      addPerElemScalarAttribute(elemType, vcg::ply::T_FLOAT, attrName, propName);
  }
  
  void AddPerElemDoubleAttribute(PlyElemType elemType, const std::string& attrName, std::string propName="")
  {
      addPerElemScalarAttribute(elemType, vcg::ply::T_DOUBLE, attrName, propName);
  }
  
  void addPerVertexScalarAttribute(const std::string& attrName, vcg::ply::PlyTypes attrType, std::string propName="")
  {
      addPerElemScalarAttribute(PlyInfo::PlyElemVertex, attrType, attrName,propName);
  }
  
  void AddPerVertexFloatAttribute(const std::string& attrName, std::string propName="")
  {
      AddPerElemFloatAttribute(PlyInfo::PlyElemVertex,attrName,propName);
  }
  
  void addPerFaceScalarAttribute(const std::string& attrName,  vcg::ply::PlyTypes attrType, std::string propName="")
  {
      addPerElemScalarAttribute(PlyInfo::PlyElemFace,attrType, attrName,propName);
  }
  
  void AddPerFaceFloatAttribute(const std::string& attrName, std::string propName="")
  {
      AddPerElemFloatAttribute(PlyInfo::PlyElemFace,attrName,propName);
  }
  
  void addPerElemPointAttribute(int elemType, vcg::ply::PlyTypes propertyType, const std::string& attrName, std::string propName="")
  {
      if(propName=="")
          propName=attrName;
      
      PropDescriptor p;
      p.propname=propName;
      p.stotype1 = propertyType;
      p.memtype1 = propertyType;
      p.islist = true;
      p.memtype2 = vcg::ply::PlyTypes::T_UCHAR;
      p.stotype2 = vcg::ply::PlyTypes::T_UCHAR;
      
      if (elemType == 0){ //vertex
          VertAttrNameVec.push_back(attrName);
          p.elemname="vertex";
          VertDescriptorVec.push_back(p);
      }
      else if (elemType == 1){ //face
          FaceAttrNameVec.push_back(attrName);
          p.elemname="face";
          FaceDescriptorVec.push_back(p);
      }
  }
  
  void addPerVertexPoint3mAttribute(const std::string& attrName, vcg::ply::PlyTypes attrType, std::string propName="")
  {
      addPerElemPointAttribute(PlyInfo::PlyElemVertex,attrType, attrName,propName);
  }
  
  void addPerVertexPoint3fAttribute(const std::string& attrName, std::string propName="")
  {
      addPerElemPointAttribute(PlyInfo::PlyElemVertex,vcg::ply::PlyTypes::T_FLOAT, attrName,propName);
  }
  
  void addPerVertexPoint3dAttribute(const std::string& attrName, std::string propName="")
  {
      addPerElemPointAttribute(PlyInfo::PlyElemVertex,vcg::ply::PlyTypes::T_DOUBLE, attrName,propName);
  }
  
  void addPerFacePoint3mAttribute(const std::string& attrName, vcg::ply::PlyTypes attrType, std::string propName="")
  {
      addPerElemPointAttribute(PlyInfo::PlyElemFace,attrType, attrName,propName);
  }
  
  void addPerFacePoint3fAttribute(const std::string& attrName, std::string propName="")
  {
      addPerElemPointAttribute(PlyInfo::PlyElemFace,vcg::ply::PlyTypes::T_FLOAT, attrName,propName);
  }
  
  void addPerFacePoint3dAttribute(const std::string& attrName, std::string propName="")
  {
      addPerElemPointAttribute(PlyInfo::PlyElemFace,vcg::ply::PlyTypes::T_DOUBLE, attrName,propName);
  }
  
  /* Note that saving a per vertex point3 attribute is a mess.
   * Actually it requires to allocate 3 float attribute and save them. And they are never deallocated... */
  template<class MeshType>
  void AddPerVertexPoint3fAttribute(MeshType &m, const std::string& attrName, std::string propName="")
  {
      if(propName=="")
          propName=attrName;
      
      const char *attrxyz[3] = {
          strdup((std::string(attrName)+std::string("_x")).c_str()),
          strdup((std::string(attrName)+std::string("_y")).c_str()),
          strdup((std::string(attrName)+std::string("_z")).c_str()),
      };
      typename MeshType::template PerVertexAttributeHandle <vcg::Point3f>
          ht = vcg::tri::Allocator<MeshType>:: template GetPerVertexAttribute <vcg::Point3f> (m,attrName);
      
      typename MeshType::template PerVertexAttributeHandle <float> htt[3];
      
      for(int i=0;i<3;++i)
      {
          htt[i] = vcg::tri::Allocator<MeshType>:: template GetPerVertexAttribute<float> (m,std::string(attrxyz[i]));
          //      ForEachVertex (m, [&](typename MeshType::VertexType &v) {
          //        htt[i][v] = ht[v][i];
          //      });
          for(auto vi=m.vert.begin();vi!=m.vert.end();++vi)
              if(!vi->IsD())
                  htt[i][vi] = ht[vi][i];
          AddPerVertexFloatAttribute(attrxyz[i]);
      }
  }
  
}; // end PlyInfo class

// This class holds the vectors of handles to the attributes of a given mesh type
// it is intialized with a PlyInfo object and a mesh and try to retrieve and initialize
// a set of handles for all the attributes that we want to load/save into ply elements
// It is used internally by the importer and exporter ply classes
// Note: we have two sets of handles, const and not const ones to handle the fact
// that for loading we need non const handles while for saving we need const handles
// (saved meshes are const) 

template <class MeshType> class PlyAttributeHelper
{    
   public:
    std::vector<typename MeshType:: template ConstPerVertexAttributeHandle<float > >        tchfv;
    std::vector<typename MeshType:: template ConstPerVertexAttributeHandle<double> >        tchdv;
    std::vector<typename MeshType:: template ConstPerVertexAttributeHandle<int   > >        tchiv;
    std::vector<typename MeshType:: template ConstPerVertexAttributeHandle<short > >        tchsv;
    std::vector<typename MeshType:: template ConstPerVertexAttributeHandle<char  > >        tchcv;
    std::vector<typename MeshType:: template ConstPerVertexAttributeHandle<unsigned char> > tchuv;
    std::vector<typename MeshType:: template ConstPerVertexAttributeHandle<vcg::Point3f> >  tchp3fv;
    std::vector<typename MeshType:: template ConstPerVertexAttributeHandle<vcg::Point3d> >  tchp3dv;
    
    std::vector<typename MeshType:: template ConstPerFaceAttributeHandle<float > >          tchff;
    std::vector<typename MeshType:: template ConstPerFaceAttributeHandle<double> >          tchdf;
    std::vector<typename MeshType:: template ConstPerFaceAttributeHandle<int   > >          tchif;
    std::vector<typename MeshType:: template ConstPerFaceAttributeHandle<short > >          tchsf;
    std::vector<typename MeshType:: template ConstPerFaceAttributeHandle<char  > >          tchcf;
    std::vector<typename MeshType:: template ConstPerFaceAttributeHandle<unsigned char> >   tchuf;
    std::vector<typename MeshType:: template ConstPerFaceAttributeHandle<vcg::Point3f> >    tchp3ff;
    std::vector<typename MeshType:: template ConstPerFaceAttributeHandle<vcg::Point3d> >    tchp3df;
    
    std::vector<typename MeshType:: template      PerVertexAttributeHandle<float > >        thfv;
    std::vector<typename MeshType:: template      PerVertexAttributeHandle<double> >        thdv;
    std::vector<typename MeshType:: template      PerVertexAttributeHandle<int   > >        thiv;
    std::vector<typename MeshType:: template      PerVertexAttributeHandle<short > >        thsv;
    std::vector<typename MeshType:: template      PerVertexAttributeHandle<char  > >        thcv;
    std::vector<typename MeshType:: template      PerVertexAttributeHandle<unsigned char> > thuv;
    std::vector<typename MeshType:: template      PerVertexAttributeHandle<vcg::Point3f> >  thp3fv;
    std::vector<typename MeshType:: template      PerVertexAttributeHandle<vcg::Point3d> >  thp3dv;
    
    std::vector<typename MeshType:: template      PerFaceAttributeHandle<float > >          thff;
    std::vector<typename MeshType:: template      PerFaceAttributeHandle<double> >          thdf;
    std::vector<typename MeshType:: template      PerFaceAttributeHandle<int   > >          thif;
    std::vector<typename MeshType:: template      PerFaceAttributeHandle<short > >          thsf;
    std::vector<typename MeshType:: template      PerFaceAttributeHandle<char  > >          thcf;
    std::vector<typename MeshType:: template      PerFaceAttributeHandle<unsigned char> >   thuf;
    std::vector<typename MeshType:: template      PerFaceAttributeHandle<vcg::Point3f> >    thp3ff;
    std::vector<typename MeshType:: template      PerFaceAttributeHandle<vcg::Point3d> >    thp3df;
   
    /* The constructor take care of initializing the handles to the attributes of a mesh
     * given a PlyInfo object and a mesh
     * Note that there are two version of the constructor one for const handler (used by save) and one for non-const handlers (used by load)
     */
    
    PlyAttributeHelper(const PlyInfo &pi, const MeshType &m)
    {
        tchfv.resize(pi.VertDescriptorVec.size());
        tchdv.resize(pi.VertDescriptorVec.size());
        tchiv.resize(pi.VertDescriptorVec.size());
        tchsv.resize(pi.VertDescriptorVec.size());
        tchcv.resize(pi.VertDescriptorVec.size());
        tchuv.resize(pi.VertDescriptorVec.size());
        tchp3fv.resize(pi.VertDescriptorVec.size());
        tchp3dv.resize(pi.VertDescriptorVec.size());
        
        for(size_t i=0;i<pi.VertDescriptorVec.size();i++)
        {
            if(!pi.VertAttrNameVec.empty() && !pi.VertAttrNameVec[i].empty())
            { // trying to use named attribute to retrieve the value to store
                assert(vcg::tri::HasPerVertexAttribute(m,pi.VertAttrNameVec[i]));
                if (!pi.VertDescriptorVec[i].islist){
                    switch (pi.VertDescriptorVec[i].stotype1)
                    {
                    case ply::T_FLOAT  : tchfv[i] = vcg::tri::Allocator<MeshType>::template FindPerVertexAttribute<float>(m,pi.VertAttrNameVec[i]); break;
                    case ply::T_DOUBLE : tchdv[i] = vcg::tri::Allocator<MeshType>::template FindPerVertexAttribute<double>(m,pi.VertAttrNameVec[i]); break;
                    case ply::T_INT    : tchiv[i] = vcg::tri::Allocator<MeshType>::template FindPerVertexAttribute<int   >(m,pi.VertAttrNameVec[i]); break;
                    case ply::T_SHORT  : tchsv[i] = vcg::tri::Allocator<MeshType>::template FindPerVertexAttribute<short >(m,pi.VertAttrNameVec[i]); break;
                    case ply::T_CHAR   : tchcv[i] = vcg::tri::Allocator<MeshType>::template FindPerVertexAttribute<char>(m,pi.VertAttrNameVec[i]); break;
                    case ply::T_UCHAR  : tchuv[i] = vcg::tri::Allocator<MeshType>::template FindPerVertexAttribute<unsigned char>(m,pi.VertAttrNameVec[i]); break;
                    default : assert(0);
                    }
                }
                else {
                    switch (pi.VertDescriptorVec[i].stotype1)
                    {
                    case ply::T_FLOAT  : tchp3fv[i] = vcg::tri::Allocator<MeshType>::template FindPerVertexAttribute<vcg::Point3f>(m,pi.VertAttrNameVec[i]); break;
                    case ply::T_DOUBLE : tchp3dv[i] = vcg::tri::Allocator<MeshType>::template FindPerVertexAttribute<vcg::Point3d>(m,pi.VertAttrNameVec[i]); break;
                    default : assert(0);
                    }
                }
            }
        }
        
        
        tchff.resize(pi.FaceDescriptorVec.size());
        tchdf.resize(pi.FaceDescriptorVec.size());
        tchif.resize(pi.FaceDescriptorVec.size());
        tchsf.resize(pi.FaceDescriptorVec.size());
        tchcf.resize(pi.FaceDescriptorVec.size());
        tchuf.resize(pi.FaceDescriptorVec.size());
        tchp3ff.resize(pi.FaceDescriptorVec.size());
        tchp3df.resize(pi.FaceDescriptorVec.size());
        
        for(size_t i=0;i<pi.FaceDescriptorVec.size();i++)
        {
            if(!pi.FaceAttrNameVec.empty() && !pi.FaceAttrNameVec[i].empty())
            { // trying to use named attribute to retrieve the value to store
                assert(vcg::tri::HasPerFaceAttribute(m,pi.FaceAttrNameVec[i]));
                if (!pi.FaceDescriptorVec[i].islist) {
                    switch (pi.FaceDescriptorVec[i].stotype1)
                    {
                    case ply::T_FLOAT  : tchff[i] = vcg::tri::Allocator<MeshType>::template FindPerFaceAttribute<float>(m,pi.FaceAttrNameVec[i]); break;
                    case ply::T_DOUBLE : tchdf[i] = vcg::tri::Allocator<MeshType>::template FindPerFaceAttribute<double>(m,pi.FaceAttrNameVec[i]); break;
                    case ply::T_INT    : tchif[i] = vcg::tri::Allocator<MeshType>::template FindPerFaceAttribute<int   >(m,pi.FaceAttrNameVec[i]); break;
                    case ply::T_SHORT  : tchsf[i] = vcg::tri::Allocator<MeshType>::template FindPerFaceAttribute<short >(m,pi.FaceAttrNameVec[i]); break;
                    case ply::T_CHAR   : tchcf[i] = vcg::tri::Allocator<MeshType>::template FindPerFaceAttribute<char>(m,pi.FaceAttrNameVec[i]); break;
                    case ply::T_UCHAR  : tchuf[i] = vcg::tri::Allocator<MeshType>::template FindPerFaceAttribute<unsigned char>(m,pi.FaceAttrNameVec[i]); break;
                    default : assert(0);
                    }
                }
                else {
                    switch (pi.FaceDescriptorVec[i].stotype1)
                    {
                    case ply::T_FLOAT  : tchp3ff[i] = vcg::tri::Allocator<MeshType>::template FindPerFaceAttribute<vcg::Point3f>(m,pi.FaceAttrNameVec[i]); break;
                    case ply::T_DOUBLE : tchp3df[i] = vcg::tri::Allocator<MeshType>::template FindPerFaceAttribute<vcg::Point3d>(m,pi.FaceAttrNameVec[i]); break;
                    default : assert(0);
                    }
                }
            }
        }
    }
        
    // Second version of the constructor
    // This one is used when we want to load attributes from a ply file and we need non const handles to the attributes
    PlyAttributeHelper(const PlyInfo &pi, MeshType &m)
    {
        thfv.resize(pi.VertDescriptorVec.size());
        thdv.resize(pi.VertDescriptorVec.size());
        thiv.resize(pi.VertDescriptorVec.size());
        thsv.resize(pi.VertDescriptorVec.size());
        thcv.resize(pi.VertDescriptorVec.size());
        thuv.resize(pi.VertDescriptorVec.size());
        thp3fv.resize(pi.VertDescriptorVec.size());
        thp3dv.resize(pi.VertDescriptorVec.size());
        
        for(size_t i=0;i<pi.VertDescriptorVec.size();i++)
        {
            if(!pi.VertAttrNameVec.empty() && !pi.VertAttrNameVec[i].empty())
            { // trying to use named attribute to retrieve the value to store
                assert(vcg::tri::HasPerVertexAttribute(m,pi.VertAttrNameVec[i]));
                if (!pi.VertDescriptorVec[i].islist){
                    switch (pi.VertDescriptorVec[i].stotype1)
                    {
                    case ply::T_FLOAT  : thfv[i] = vcg::tri::Allocator<MeshType>::template FindPerVertexAttribute<float>(m,pi.VertAttrNameVec[i]); break;
                    case ply::T_DOUBLE : thdv[i] = vcg::tri::Allocator<MeshType>::template FindPerVertexAttribute<double>(m,pi.VertAttrNameVec[i]); break;
                    case ply::T_INT    : thiv[i] = vcg::tri::Allocator<MeshType>::template FindPerVertexAttribute<int   >(m,pi.VertAttrNameVec[i]); break;
                    case ply::T_SHORT  : thsv[i] = vcg::tri::Allocator<MeshType>::template FindPerVertexAttribute<short >(m,pi.VertAttrNameVec[i]); break;
                    case ply::T_CHAR   : thcv[i] = vcg::tri::Allocator<MeshType>::template FindPerVertexAttribute<char>(m,pi.VertAttrNameVec[i]); break;
                    case ply::T_UCHAR  : thuv[i] = vcg::tri::Allocator<MeshType>::template FindPerVertexAttribute<unsigned char>(m,pi.VertAttrNameVec[i]); break;
                    default : assert(0);
                    }
                }
                else {
                    switch (pi.VertDescriptorVec[i].stotype1)
                    {
                    case ply::T_FLOAT  : thp3fv[i] = vcg::tri::Allocator<MeshType>::template FindPerVertexAttribute<vcg::Point3f>(m,pi.VertAttrNameVec[i]); break;
                    case ply::T_DOUBLE : thp3dv[i] = vcg::tri::Allocator<MeshType>::template FindPerVertexAttribute<vcg::Point3d>(m,pi.VertAttrNameVec[i]); break;
                    default : assert(0);
                    }
                }
            }
        }
        
        
        thff.resize(pi.FaceDescriptorVec.size());
        thdf.resize(pi.FaceDescriptorVec.size());
        thif.resize(pi.FaceDescriptorVec.size());
        thsf.resize(pi.FaceDescriptorVec.size());
        thcf.resize(pi.FaceDescriptorVec.size());
        thuf.resize(pi.FaceDescriptorVec.size());
        thp3ff.resize(pi.FaceDescriptorVec.size());
        thp3df.resize(pi.FaceDescriptorVec.size());
        
        for(size_t i=0;i<pi.FaceDescriptorVec.size();i++)
        {
            if(!pi.FaceAttrNameVec.empty() && !pi.FaceAttrNameVec[i].empty())
            { // trying to use named attribute to retrieve the value to store
                assert(vcg::tri::HasPerFaceAttribute(m,pi.FaceAttrNameVec[i]));
                if (!pi.FaceDescriptorVec[i].islist) {
                    switch (pi.FaceDescriptorVec[i].stotype1)
                    {
                    case ply::T_FLOAT  : thff[i] = vcg::tri::Allocator<MeshType>::template FindPerFaceAttribute<float>(m,pi.FaceAttrNameVec[i]); break;
                    case ply::T_DOUBLE : thdf[i] = vcg::tri::Allocator<MeshType>::template FindPerFaceAttribute<double>(m,pi.FaceAttrNameVec[i]); break;
                    case ply::T_INT    : thif[i] = vcg::tri::Allocator<MeshType>::template FindPerFaceAttribute<int   >(m,pi.FaceAttrNameVec[i]); break;
                    case ply::T_SHORT  : thsf[i] = vcg::tri::Allocator<MeshType>::template FindPerFaceAttribute<short >(m,pi.FaceAttrNameVec[i]); break;
                    case ply::T_CHAR   : thcf[i] = vcg::tri::Allocator<MeshType>::template FindPerFaceAttribute<char>(m,pi.FaceAttrNameVec[i]); break;
                    case ply::T_UCHAR  : thuf[i] = vcg::tri::Allocator<MeshType>::template FindPerFaceAttribute<unsigned char>(m,pi.FaceAttrNameVec[i]); break;
                    default : assert(0);
                    }
                }
                else {
                    switch (pi.FaceDescriptorVec[i].stotype1)
                    {
                    case ply::T_FLOAT  : thp3ff[i] = vcg::tri::Allocator<MeshType>::template FindPerFaceAttribute<vcg::Point3f>(m,pi.FaceAttrNameVec[i]); break;
                    case ply::T_DOUBLE : thp3df[i] = vcg::tri::Allocator<MeshType>::template FindPerFaceAttribute<vcg::Point3d>(m,pi.FaceAttrNameVec[i]); break;
                    default : assert(0);
                    }
                }
            }
        }
    }

}; // end PlyAttributeHelper class


} // end namespace tri
} // end namespace io
} // end namespace vcg
#endif
