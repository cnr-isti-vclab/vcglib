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



#ifndef __VCGLIB_EXPORT_OBJ
#define __VCGLIB_EXPORT_OBJ

#include <vcg/complex/algorithms/polygon_support.h>
#include <wrap/system/utf8_file.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <wrap/callback.h>
#include <wrap/io_trimesh/io_mask.h>
#include <wrap/io_trimesh/io_material.h>
#include <wrap/io_trimesh/precision.h>
#include <iostream>
#include <fstream>
#include <map>

namespace vcg {
namespace tri {
namespace io {

template <class SaveMeshType>
class ExporterOBJ
{
public:
  typedef typename SaveMeshType::FaceIterator FaceIterator;
  typedef typename SaveMeshType::ConstFaceIterator ConstFaceIterator;
  typedef typename SaveMeshType::EdgeIterator EdgeIterator;
  typedef typename SaveMeshType::ConstEdgeIterator ConstEdgeIterator;
  typedef typename SaveMeshType::VertexIterator VertexIterator;
  typedef typename SaveMeshType::VertexType VertexType;
  typedef typename SaveMeshType::VertexPointer VertexPointer;
  typedef typename SaveMeshType::FacePointer FacePointer;
  typedef typename SaveMeshType::ScalarType ScalarType;
  typedef typename SaveMeshType::CoordType CoordType;
  typedef typename SaveMeshType::FaceType::TexCoordType TexCoordType;

  /*
            enum of all the types of error
        */
  enum SaveError
  {
    E_NOERROR,            // 0
    E_CANTOPENFILE,       // 1
    E_CANTCLOSEFILE,      // 2
    E_UNESPECTEDEOF,      // 3
    E_ABORTED,            // 4
    E_NOTDEFINITION,      // 5
    E_NO_VERTICES,        // 6
    E_NOTFACESVALID,      // 7
    E_NO_VALID_MATERIAL,  // 8
    E_STREAMERROR         // 9
  };

  /*
            this function takes an index and the relative error message gets back
        */
  static const char* ErrorMsg(int error)
  {
    static const char* obj_error_msg[] =
    {
      "No errors",  // 0
      "Can't open file",  // 1
      "can't close file",  // 2
      "Premature End of file",  // 3
      "File saving aborted",  // 4
      "Function not defined",  // 5
      "Vertices not valid",  // 6
      "Faces not valid",  // 7
      "The mesh has not a attribute containing the vector of materials",  // 8
      "Output Stream Error" //9
    };

    if(error>9 || error<0) return "Unknown error";
    else return obj_error_msg[error];
  };

  /*
            returns mask of capability one define with what are the saveable information of the format.
        */
  static int GetExportMaskCapability()
  {
    int capability = 0;

    //vert
    capability |= vcg::tri::io::Mask::IOM_VERTCOORD;
    capability |= vcg::tri::io::Mask::IOM_VERTNORMAL;
    capability |= vcg::tri::io::Mask::IOM_VERTTEXCOORD;
    capability |= vcg::tri::io::Mask::IOM_VERTCOLOR;

    //face
    capability |= vcg::tri::io::Mask::IOM_FACEINDEX;
    capability |= vcg::tri::io::Mask::IOM_FACECOLOR;

    capability |= vcg::tri::io::Mask::IOM_EDGEINDEX;
    //wedg
    capability |= vcg::tri::io::Mask::IOM_WEDGTEXCOORD;
    capability |= vcg::tri::io::Mask::IOM_WEDGNORMAL;
    capability |= vcg::tri::io::Mask::IOM_BITPOLYGONAL;

    return capability;
  }
  
  static int Save(SaveMeshType &m, const char * filename, int mask, CallBackPos *cb=0)
  {
    return Save(m,filename,mask,false,cb);
  }
  
  /** 
   * main function to export a mesh in OBJ file format
   * 
   * if you enable the useMaterialAttribute flag, the exporter will assume that the mesh has a consistent per mesh attribute and a per face attribute containing the index of the material
   */
  static int Save(SaveMeshType &m, const char * filename, int mask, bool useMaterialAttribute ,CallBackPos *cb=0)
  {
    const int DGT = vcg::tri::io::Precision<ScalarType>::digits();
    const bool savePolygons = (mask & Mask::IOM_BITPOLYGONAL) != 0;
    if (savePolygons) {
      tri::RequireFFAdjacency(m);
      tri::UpdateTopology<SaveMeshType>::FaceFace(m);
    }
    // OBJ has one index stream for each attribute, so prefer per-wedge data
    // whenever both per-vertex and per-wedge variants were requested.
    if (mask & vcg::tri::io::Mask::IOM_WEDGTEXCOORD &&
        mask & vcg::tri::io::Mask::IOM_VERTTEXCOORD ) {
      mask &= ~vcg::tri::io::Mask::IOM_VERTTEXCOORD;
    }
    if (mask & vcg::tri::io::Mask::IOM_WEDGNORMAL &&
        mask & vcg::tri::io::Mask::IOM_VERTNORMAL ) {
      mask &= ~vcg::tri::io::Mask::IOM_VERTNORMAL;
    }
    if (mask & vcg::tri::io::Mask::IOM_WEDGCOLOR &&
        mask & vcg::tri::io::Mask::IOM_VERTCOLOR ) {
      mask &= ~vcg::tri::io::Mask::IOM_VERTCOLOR;
    }
    if(m.vn == 0)	return E_NO_VERTICES;
    
    typename SaveMeshType::template PerMeshAttributeHandle<std::vector<Material> > materialVecHandle =
        vcg::tri::Allocator<SaveMeshType>::template FindPerMeshAttribute<std::vector<Material> >(m, "materialVector");
    typename SaveMeshType::template PerFaceAttributeHandle<int> materialIndexHandle =
            vcg::tri::Allocator<SaveMeshType>::template FindPerFaceAttribute<int>(m, "materialIndex");
    
    if(useMaterialAttribute && (!Allocator<SaveMeshType>::IsValidHandle(m,materialVecHandle)) && 
       (!Allocator<SaveMeshType>::IsValidHandle(m,materialIndexHandle)) ) 
        return E_NO_VALID_MATERIAL;

    FILE *fp = vcg::utf8::FOpen(filename,"w");
    if(fp == NULL) return E_CANTOPENFILE;
    std::string shortFilename(filename);
    int LastSlash=shortFilename.size()-1;
    while(LastSlash>=0 && shortFilename[LastSlash]!='/')
      --LastSlash;
    shortFilename = shortFilename.substr(LastSlash+1);
    
    fprintf(fp,"####\n#\n# OBJ File Generated by Meshlab\n#\n####\n");
    fprintf(fp,"# Object %s\n#\n# Vertices: %d\n# Faces: ",shortFilename.c_str(),m.vn);
    // The polygon count is only known after faux-edge groups have been
    // traversed. Reserve a fixed-width field and update it after writing.
    const long faceCountPosition = ftell(fp);
    fprintf(fp,"%-20d\n#\n####\n",m.fn);
    
    //library materialVec
    if( (mask & vcg::tri::io::Mask::IOM_FACECOLOR)  || (mask & Mask::IOM_WEDGTEXCOORD) || (mask & Mask::IOM_VERTTEXCOORD) )
      fprintf(fp,"mtllib ./%s.mtl\n\n",shortFilename.c_str());

    std::map<CoordType,int> NormalVertex;
    std::vector<int> VertexId(m.vert.size());
    int numvert = 0;
    int curNormalIndex = 1;
    int current = 0;
    const int totalPrimitives = m.vn+m.fn;    
    /*********************************** VERTICES *********************************/
    for(auto vi=m.vert.begin(); vi!=m.vert.end(); ++vi) if( !(*vi).IsD() )
    {
      VertexId[vi-m.vert.begin()]=numvert;
      if (mask & Mask::IOM_VERTNORMAL ) {
        fprintf(fp,"vn %.*f %.*f %.*f\n", DGT, (*vi).N()[0], DGT, (*vi).N()[1], DGT, (*vi).N()[2]);
      }
      
      if (mask & Mask::IOM_VERTTEXCOORD ) {
        fprintf(fp,"vt %.*f %.*f\n", DGT, (*vi).T().P()[0], DGT, (*vi).T().P()[1]);
      }

      fprintf(fp,"v %.*f %.*f %.*f", DGT, (*vi).P()[0], DGT, (*vi).P()[1], DGT, (*vi).P()[2]);
      
      if(mask & Mask::IOM_VERTCOLOR) // the socially accepted extension to the obj format. 
        fprintf(fp," %f %f %f", double((*vi).C()[0])/255., double((*vi).C()[1])/255., double((*vi).C()[2])/255.);
      fprintf(fp,"\n");

      if (cb !=NULL)
      {
        if(!(*cb)((100*++current)/totalPrimitives, "writing vertices "))
        {
          fclose(fp);
          return E_ABORTED;
        }
      }
      numvert++;
    }
    assert(numvert == m.vn);
    fprintf(fp,"# %d vertices, %d vertices normals\n\n",m.vn,int(NormalVertex.size()));

    /********************* FACES ************************/      
    //faces + texture coords
    std::map<TexCoordType,int> CoordIndexTexture;
    int curTexCoordIndex = 1;
    int curMatIndex = -1;
    int writtenFaces = 0;
    std::vector<Material> materialVec; //used if we do not have material attributes 

    typedef std::pair<FacePointer, int> BoundaryCorner;
    std::vector<VertexPointer> polygonVertices;
    std::vector<FacePointer> polygonFaces;
    std::vector<BoundaryCorner> boundaryCorners;
    boundaryCorners.reserve(4);
    if(savePolygons)
      tri::UpdateFlags<SaveMeshType>::FaceClearV(m);

    for(FaceIterator fi=m.face.begin(); fi!=m.face.end(); ++fi)
    {
      if((*fi).IsD() || (savePolygons && (*fi).IsV()))
        continue;

      FacePointer representativeFace = &*fi;
      boundaryCorners.clear();
      polygonFaces.clear();
      if(savePolygons)
      {
        vcg::tri::PolygonSupport<SaveMeshType,SaveMeshType>::ExtractPolygon(
          representativeFace, polygonVertices, polygonFaces, boundaryCorners);
        // ExtractPolygon follows the boundary clockwise for multi-triangle
        // polygons; restore the source face winding expected by OBJ.
        if(boundaryCorners.size()>3)
          std::reverse(boundaryCorners.begin(), boundaryCorners.end());
        if(boundaryCorners.empty())
        {
          tri::UpdateFlags<SaveMeshType>::FaceClearV(m);
          fclose(fp);
          return E_NOTFACESVALID;
        }
      }
      else
      {
        polygonFaces.push_back(representativeFace);
        for(int k=0;k<representativeFace->VN();++k)
          boundaryCorners.emplace_back(representativeFace,k);
      }

      if((mask & Mask::IOM_FACECOLOR) || (mask & Mask::IOM_WEDGTEXCOORD) || (mask & Mask::IOM_VERTTEXCOORD))
      {
        int index=-1;
        if(useMaterialAttribute) index = materialIndexHandle[representativeFace];
        else                     index = Materials<SaveMeshType>::CreateNewMaterial(m,materialVec,*representativeFace);
                  
        if(index != curMatIndex) {
          fprintf(fp,"\nusemtl material_%d\n", index);
          curMatIndex = index;
        }
      }

      //saves texture coord x wedge
      if(HasPerWedgeTexCoord(m) && (mask & Mask::IOM_WEDGTEXCOORD))
        for(const BoundaryCorner &corner : boundaryCorners)
        {
            const TexCoordType &texCoord = corner.first->cWT(corner.second);
            if(AddNewTextureCoord(CoordIndexTexture,texCoord,curTexCoordIndex))
            {
              fprintf(fp,"vt %.*f %.*f\n", DGT, texCoord.u(), DGT, texCoord.v());
              curTexCoordIndex++; //ncreases the value number to be associated to the Texture
            }
        }

      if(HasPerWedgeNormal(m) && (mask & Mask::IOM_WEDGNORMAL))
        for(const BoundaryCorner &corner : boundaryCorners)
        {
          const CoordType normal(corner.first->cWN(corner.second));
          if(AddNewNormalVertex(NormalVertex,normal,curNormalIndex))
          {
            fprintf(fp,"vn %.*f %.*f %.*f\n", DGT, normal[0], DGT, normal[1], DGT, normal[2]);
            curNormalIndex++;
          }
        }

      fprintf(fp,"f ");
      for(size_t k=0;k<boundaryCorners.size();k++)
      {
        if(k!=0) fprintf(fp," ");
        const BoundaryCorner &corner = boundaryCorners[k];
        int vInd = -1;
        // +1 because Obj file format begins from index = 1 but not from index = 0.
        vInd = VertexId[tri::Index(m, corner.first->V(corner.second))] + 1;//index of vertex per face

        int vt = -1;
        if(mask & Mask::IOM_WEDGTEXCOORD)
          vt = GetIndexVertexTexture(CoordIndexTexture,corner.first->cWT(corner.second));//index of vertex texture per face
        if (mask & Mask::IOM_VERTTEXCOORD)
          vt = vInd;

        int vn = -1;
        if(HasPerWedgeNormal(m) && (mask & Mask::IOM_WEDGNORMAL))
          vn = GetIndexVertexNormal(m, NormalVertex, corner.first->cWN(corner.second));//index of wedge normal
        if (mask & Mask::IOM_VERTNORMAL)
          vn = vInd;

        //writes elements on file obj
        WriteFacesElement(fp,vInd,vt,vn);
      }
      fprintf(fp,"\n");
      ++writtenFaces;

      if (cb !=NULL) {
        current += int(polygonFaces.size());
        if(!(*cb)((100*current)/totalPrimitives, "writing faces "))
        {
          if(savePolygons)
            tri::UpdateFlags<SaveMeshType>::FaceClearV(m);
          fclose(fp);
          return E_ABORTED;
        }
      }

    } // end for faces

    if(savePolygons)
      tri::UpdateFlags<SaveMeshType>::FaceClearV(m);

    for(ConstEdgeIterator ei=m.edge.begin(); ei!=m.edge.end(); ++ei) if( !(*ei).IsD() )
    {
      fprintf(fp,"l %i %i\n",
              VertexId[tri::Index(m, (*ei).V(0))] + 1,
              VertexId[tri::Index(m, (*ei).V(1))] + 1);
    }

    fprintf(fp,"# %d faces, %d coords texture\n\n",writtenFaces,int(CoordIndexTexture.size()));

    const long endPosition = ftell(fp);
    if(faceCountPosition >= 0 && endPosition >= 0 && fseek(fp,faceCountPosition,SEEK_SET) == 0)
    {
      fprintf(fp,"%-20d",writtenFaces);
      fseek(fp,endPosition,SEEK_SET);
    }

    fprintf(fp,"# End of File\n");

    int errCode = E_NOERROR;
    if((mask & Mask::IOM_WEDGTEXCOORD) || (mask & Mask::IOM_FACECOLOR) || (mask & Mask::IOM_VERTTEXCOORD) )
    {
      if(useMaterialAttribute) errCode = WriteMaterials(materialVecHandle(), filename,cb);
      else                     errCode = WriteMaterials(materialVec, filename,cb);
    }

    int result = E_NOERROR;
    if (errCode != E_NOERROR)
        result = errCode;
    else if (ferror(fp))
        result = E_STREAMERROR;
    fclose(fp);
    return result;
  }

   /*
            returns index of the texture coord
        */
  inline static int GetIndexVertexTexture(typename std::map<TexCoordType,int> &mapTexToInt, const TexCoordType &wt)
  {
    typename std::map<TexCoordType,int>::iterator iter= mapTexToInt.find(wt);
    if(iter != mapTexToInt.end()) return (*iter).second;
    else 		return -1;
    // Old wrong version.
    // int index = mapTexToInt[wt];
    // if(index!=0){return index;}
  }

  /*
            returns index of the vertex normal
        */
  inline static int GetIndexVertexNormal(const SaveMeshType &/*m*/, std::map<CoordType,int> &mapNormToInt, const CoordType &norm )
  {
    typename std::map<CoordType,int>::iterator iter= mapNormToInt.find(norm);
    if(iter != mapNormToInt.end()) return (*iter).second;
    else 		return -1;
  }

  /*
            write elements on file

                f v/tc/vn v/tc/vn v/tc/vn ...
                f v/tc v/tc v/tc ...
                f v//vn v//vn v//vn ...
                f v v v ...

        */
  inline static void WriteFacesElement(FILE *fp,int v,int vt, int vn)
  {
    fprintf(fp,"%d",v);
    if(vt!=-1)
    {
      fprintf(fp,"/%d",vt);
      if(vn!=-1)
        fprintf(fp,"/%d",vn);
    }
    else if(vn!=-1)
      fprintf(fp,"//%d",vn);
  }

  /*
            adds a new index to the coordinate of Texture if it is the first time
            which is otherwise met not execute anything
        */
  template <class TexScalarType>
  inline static bool AddNewTextureCoord(std::map<typename vcg::TexCoord2<TexScalarType>,int> &m,
                                        const typename vcg::TexCoord2<TexScalarType> &wt,int value)
  {
    int index = m[wt];
    if(index==0){m[wt]=value;return true;}
    return false;
  }

  /*
            adds a new index to the normal per vertex if it is the first time
            which is otherwise met does not execute anything
        */
  inline static bool AddNewNormalVertex(typename std::map<CoordType,int> &m, const CoordType &n ,int value)
  {
    int index = m[n];
    if(index==0){m[n]=value;return true;}
    return false;
  }

  /*
            writes material into file
        */
  inline static int WriteMaterials(std::vector<Material> &materialVec, const char * filename, CallBackPos *cb=0)
  {
    std::string fileName = std::string(filename);
    fileName+=".mtl";

    if(materialVec.size() > 0)
    {
      FILE *fp;
      fp = vcg::utf8::FOpen(fileName.c_str(),"w");
      if(fp==NULL)return E_ABORTED;

      fprintf(fp,"#\n# Wavefront material file\n# Converted by Meshlab Group\n#\n\n");

      int current = 0;

      for(unsigned int i=0;i<materialVec.size();i++)
      {
        if (cb !=NULL)
          (*cb)((100 * ++current)/materialVec.size(), "saving material file ");
        else
        { /* fclose(fp); return E_ABORTED; */ }

        fprintf(fp,"newmtl material_%d\n",i);
        fprintf(fp,"Ka %f %f %f\n",materialVec[i].Ka[0],materialVec[i].Ka[1],materialVec[i].Ka[2]);
        fprintf(fp,"Kd %f %f %f\n",materialVec[i].Kd[0],materialVec[i].Kd[1],materialVec[i].Kd[2]);
        fprintf(fp,"Ks %f %f %f\n",materialVec[i].Ks[0],materialVec[i].Ks[1],materialVec[i].Ks[2]);
        fprintf(fp,"Tr %f\n",materialVec[i].Tr);
        fprintf(fp,"illum %d\n",materialVec[i].illum);
        fprintf(fp,"Ns %f\n",materialVec[i].Ns);

        if(materialVec[i].map_Kd.size()>0)
          fprintf(fp,"map_Kd %s\n",materialVec[i].map_Kd.c_str());
        fprintf(fp,"\n");
      }
      fclose(fp);
    }
    return E_NOERROR;
  }

}; // end class
} // end Namespace tri
} // end Namespace io
} // end Namespace vcg

#endif
