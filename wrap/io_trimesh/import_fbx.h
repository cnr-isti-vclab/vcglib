/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2018                                           \/)\/    *
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

/**
@name Open FBX format
*/
//@{

#ifndef VCGLIB_IMPORT_FBX
#define VCGLIB_IMPORT_FBX
#include <wrap/openfbx/src/ofbx.h>

#include <fstream>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/clean.h>
#include <wrap/io_trimesh/io_mask.h>

/**
\brief ImporterFBX is the class devoted to import the info contained in a FBX file inside a mesh defined following the vcg standards.

To use it you need to add the following two files to your project:

wrap/openfbx/src/ofbx.cpp  
wrap/openfbx/src/miniz.c  


*/
namespace vcg {
namespace tri {
namespace io {

template <class OpenMeshType>
class ImporterFBX
{
public:
  
  typedef typename OpenMeshType::VertexPointer VertexPointer;
  typedef typename OpenMeshType::ScalarType ScalarType;
  typedef typename OpenMeshType::VertexType VertexType;
  typedef typename OpenMeshType::FaceType FaceType;
  typedef typename OpenMeshType::VertexIterator VertexIterator;
  typedef typename OpenMeshType::FaceIterator FaceIterator;
  typedef typename OpenMeshType::CoordType CoordType;
  typedef typename OpenMeshType::VertexType::ColorType ColorType;
  typedef typename OpenMeshType::VertexType::ColorType::ScalarType CSType;
  typedef typename vcg::Matrix44<ScalarType> Matrix44Type;
  
public:
  /*****
        Enum containing possible error codes you can get opening an FBX file.
    */
  enum FBXError 
  {
    // Successfull opening
    E_NOERROR = 0,													
    
    // Critical Opening Errors (only even numbers)
  };
  
  /****
        Function devoted to check if a given error is critical or the elaboration can go on.
    */
  static bool ErrorCritical(int err)
  { 
    if(err == E_NOERROR) 
      return false;
    return true;
  }
  
  /****
        The function return the sting info associated with an error code.
    */
  static const char* ErrorMsg(int error)
  {
    static const char* fbx_error_msg[] =
    {
      "No errors",
    };
    
    if(error>1 || error<0) 
      return "Unknown error";
    else 
      return fbx_error_msg[error];
  }
  
  /****
        The open function import info contained in a FBX file into a mesh defined using the VCG standards.
        Parameters:
            - m is the mesh in which you want to load the file info
            - filename is the filename path of a valid FBX file. It could be a relative path or an absolute one.
            - oi is a small structure devoted to contain some feedbacks about the loaded file (like vertex number or face number). Typically you can ignore it. 
            
        Semantics:
            If the file has been correctly loaded (i.e. Open function returned a E_NOERROR code, you will find in the mesh m at least:
                - a vertex position for each vertex
                - the triangles componing the mesh surface
                
            Depending on which type mesh you passed to the function you could also have:
                - a vertex normal for each vertex
                - a vertex color (ONLY IF THE FBX FILE HAS ALL THE MATERIALS WITHOUT TEXTURES AT ALL)
                - a single material for each triangle. We are using the one containing a diffuse texture and/or a diffuse color. 
                - a UV parameterization for each vertex componing a triangle
                
            Please note that if the file is composed by a set of meshes we compact all these meshes in a single one.
    */
  
  static int Open( OpenMeshType &m, const char * filename, CallBackPos *cb = nullptr)
  {
    
    FILE* fp = fopen(filename, "rb");
    if (!fp) return false;
    
    fseek(fp, 0, SEEK_END);
    long file_size = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    ofbx::u8* content = new ofbx::u8[file_size];
    fread(content, 1, size_t(file_size), fp);
    fclose(fp);
    
    ofbx::IScene* scene =ofbx::load(content, int(file_size));   
    
    m.Clear();
    int mesh_count = scene->getMeshCount();
    printf("mesh count %i\n",mesh_count);
    int totVert=0;
    for (int i = 0; i < mesh_count; ++i)
      totVert += scene->getMesh(i)->getGeometry()->getVertexCount();
    
    int baseVertexCount;
    printf("Scene has %i meshes\n",mesh_count);
    for (int i = 0; i < mesh_count; ++i)
    {
      baseVertexCount=m.vert.size();
      int baseTextureIndexCount = m.textures.size();

      const ofbx::Mesh& mesh = *scene->getMesh(i);
      const ofbx::Geometry& geom = *mesh.getGeometry();
      const int *materials = geom.getMaterials();

      // Load textured materials. The mesh object has pointers to the materials,
      // the geometry object instead has a vector of FN material indices
      for (int mat_i = 0; mat_i < mesh.getMaterialCount(); ++mat_i) {
        const ofbx::Material *mat = mesh.getMaterial(mat_i);
        assert(mat);
        const ofbx::Texture  *tex = mat->getTexture(ofbx::Texture::DIFFUSE);
        char buf[1024];
        if (tex) {
          tex->getRelativeFileName().toString(buf); buf[1023] = 0;
          if (std::string(buf) == "")  {
            tex->getFileName().toString(buf); buf[1023] = 0;
          }
          printf("Texture %s\n", buf);
          m.textures.push_back(buf);
        }
      }

      Matrix44d globTransfd(mesh.getGlobalTransform().m);
      /*
      Matrix44f globTransf;
      globTransf.SetIdentity();
      globTransf.Import(globTransfd);
      Transpose(globTransf);
      */
      Transpose(globTransfd);
      int vertex_count = geom.getVertexCount();
      const ofbx::Vec3* vertices = geom.getVertices();
      const ofbx::Vec2* texcoords = geom.getUVs();
      printf("  Mesh %i has %i vertices\n",i,vertex_count);
      
      for (int j = 0; j < vertex_count; ++j)
      {
        ofbx::Vec3 vt = vertices[j];
        Point3d vd = globTransfd*Point3d(vt.x, vt.y, vt.z);
        CoordType v(vd.X(), vd.Y(), vd.Z());
        tri::Allocator<OpenMeshType>::AddVertex(m,v);
        if(cb && (m.vert.size()%1000) == 0 ) cb(m.vert.size() * 100 / totVert, "Vertex Loading");
      }
      
//      bool has_normals = geom.getNormals() != nullptr;
      bool has_uvs = texcoords != nullptr;           
      
      const ofbx::Skin* skin = geom.getSkin();
      if(!skin)
      {
        printf("Geom has no skin\n");
        for (int j = 0; j < vertex_count/3; ++j)
        {
          int curTriBase = baseVertexCount + j*3;
          tri::Allocator<OpenMeshType>::AddFace(m,curTriBase,curTriBase+1,curTriBase+2);   
          if(has_uvs) 
          {
            for(int k=0;k<3;++k)
            {
              m.face.back().WT(k).u() = texcoords[j*3+k].x;
              m.face.back().WT(k).v() = texcoords[j*3+k].y;
              if (materials != nullptr)
                m.face.back().WT(k).n() = baseTextureIndexCount + materials[j];
            }
          }
        }
      }
      else
      {
        int cluster_count = skin->getClusterCount();
        
        printf("Mesh %i has %i clusters\n",i,cluster_count);       
        for (int k = 0; k < cluster_count; ++k)
        {        
//          const ofbx::Cluster *cluster=skin->getCluster(k);
        }
      }
    }
    
    tri::Clean<OpenMeshType>::RemoveDuplicateVertex(m);
    tri::Allocator<OpenMeshType>::CompactEveryVector(m);
    return E_NOERROR;
    
  } 
  
  /****
    The LoadMask function gives you a preview of which attributes you will find inside the fbx file. 
    It's useful only if you are planning to use a mesh with dynamic attributes.
    
    Parameters:
        - filename is the filename path of a valid FBX file. It could be a relative path or an absolute one.
        - oi is a small structure devoted to contain some feedbacks about the loaded file (like vertex number or face number). 
        In oi.mask you will find the mask containing the attributes contained in the file. You can query the presence of a specific attribute using the operators on bits.
        ( for example if I want to check if the file contains the vertex normal layer I have to call the LoadMask function and after check if (oi.mask & vcg::tri::io::Mask::IOM_VERTNORMAL) is true).
    */
  
  static int LoadMask(const char * /*filename*/)
  {
    return tri::io::Mask::IOM_WEDGTEXCOORD;
  }
  
}; // end class
} // end namespace tri
} // end namespace io
} // end namespace vcg


#endif 
