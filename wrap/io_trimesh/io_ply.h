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
/****************************************************************************
  History

$Log: not supported by cvs2svn $
Revision 1.3  2004/05/12 10:19:30  ganovelli
new line added at the end of file

Revision 1.2  2004/03/09 21:26:47  cignoni
cr lf mismatch

Revision 1.1  2004/03/08 09:21:34  cignoni
Initial commit

Revision 1.1  2004/03/03 15:00:51  cignoni
Initial commit

****************************************************************************/
#ifndef __VCGLIB_IOTRIMESH_IO_PLY
#define __VCGLIB_IOTRIMESH_IO_PLY


#include <vcg/complex/allocate.h>

/**
@name Load and Save in Ply format
*/
//@{
#include<wrap/callback.h>
#include<wrap/ply/plylib.h>

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

	void addPerElemScalarAttribute(int elemType, vcg::ply::PlyTypes propertyType, const std::string& attrName, std::string propName="")
	{
		if(propName=="")
			propName=attrName;

		PropDescriptor p;
		p.propname=propName;
		p.islist = false;
		p.stotype1 = propertyType;
		p.memtype1 = propertyType;

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

	void AddPerElemFloatAttribute(int elemType, const std::string& attrName, std::string propName="")
	{
		addPerElemScalarAttribute(elemType, vcg::ply::T_FLOAT, attrName, propName);
	}

	void AddPerElemDoubleAttribute(int elemType, const std::string& attrName, std::string propName="")
	{
		addPerElemScalarAttribute(elemType, vcg::ply::T_DOUBLE, attrName, propName);
	}

	void addPerVertexScalarAttribute(const std::string& attrName, vcg::ply::PlyTypes attrType, std::string propName="")
	{
		addPerElemScalarAttribute(0,attrType, attrName,propName);
	}

	void AddPerVertexFloatAttribute(const std::string& attrName, std::string propName="")
	{
		AddPerElemFloatAttribute(0,attrName,propName);
	}

	void addPerFaceScalarAttribute(const std::string& attrName,  vcg::ply::PlyTypes attrType, std::string propName="")
	{
		addPerElemScalarAttribute(1,attrType, attrName,propName);
	}

	void AddPerFaceFloatAttribute(const std::string& attrName, std::string propName="")
	{
		AddPerElemFloatAttribute(1,attrName,propName);
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
		addPerElemPointAttribute(0,attrType, attrName,propName);
	}

	void addPerVertexPoint3fAttribute(const std::string& attrName, std::string propName="")
	{
		addPerElemPointAttribute(0,vcg::ply::PlyTypes::T_FLOAT, attrName,propName);
	}

	void addPerVertexPoint3dAttribute(const std::string& attrName, std::string propName="")
	{
		addPerElemPointAttribute(0,vcg::ply::PlyTypes::T_DOUBLE, attrName,propName);
	}

	void addPerFacePoint3mAttribute(const std::string& attrName, vcg::ply::PlyTypes attrType, std::string propName="")
	{
		addPerElemPointAttribute(1,attrType, attrName,propName);
	}

	void addPerFacePoint3fAttribute(const std::string& attrName, std::string propName="")
	{
		addPerElemPointAttribute(1,vcg::ply::PlyTypes::T_FLOAT, attrName,propName);
	}

	void addPerFacePoint3dAttribute(const std::string& attrName, std::string propName="")
	{
		addPerElemPointAttribute(1,vcg::ply::PlyTypes::T_DOUBLE, attrName,propName);
	}

  /* Note that saving a per vertex point3 attribute is a mess.
   * Actually require to allocate 3 float attribute and save them. And they are never deallocated... */
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


  PlyInfo()
  {
    status=0;
    mask=0;
    cb=0;
  }
  /// Store the error codes enconutered when parsing a ply
  int status;
  /// It returns a bit mask describing the field present in the ply file
  int mask;

  /// a Simple callback that can be used for long ply parsing.
  // it returns the current position, and formats a string with a description of what th efunction is doing (loading vertexes, faces...)
  CallBackPos *cb;

  /// The additional vertex descriptor that a user can specify to load additional per-vertex non-standard data stored in a ply
  std::vector<PropDescriptor> VertDescriptorVec;
  /// AttributeName is an array, externally allocated, containing the names of the attributes to be saved (loaded).
  /// We assume that AttributeName[], if not empty, is exactly of the same size of VertexdData[]
  /// If AttributeName[i] is not empty we use it to retrieve/store the info instead of the offsetted space in the current vertex
  std::vector<std::string> VertAttrNameVec;

  /// The additional vertex descriptor that a user can specify to load additional per-face non-standard data stored in a ply
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

}; // end class
} // end namespace tri
} // end namespace io
} // end namespace vcg
#endif
