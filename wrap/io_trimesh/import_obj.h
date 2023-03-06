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


#ifndef __VCGLIB_IMPORT_OBJ
#define __VCGLIB_IMPORT_OBJ

#include <vcg/complex/allocate.h>

#include <wrap/callback.h>
#include <wrap/io_trimesh/io_mask.h>
#include <wrap/io_trimesh/io_material.h>
#include <wrap/io_trimesh/io_fan_tessellator.h>
#ifdef __gl_h_
#include <wrap/gl/glu_tesselator.h>
#endif
#include <vcg/space/color4.h>


#include <fstream>
#include <string>
#include <vector>


namespace vcg {
namespace tri {
namespace io {

/**
This class encapsulate a filter for importing obj (Alias Wavefront) meshes.
Warning: this code assume little endian (PC) architecture!!!
*/

template <class OpenMeshType>
class ImporterOBJ
{
public:
  static int &MRGBLineCount(){static int _MRGBLineCount=0; return _MRGBLineCount;}
  
  typedef typename OpenMeshType::VertexPointer VertexPointer;
  typedef typename OpenMeshType::ScalarType ScalarType;
  typedef typename OpenMeshType::VertexType VertexType;
  typedef typename OpenMeshType::EdgeType   EdgeType;
  typedef typename OpenMeshType::FaceType FaceType;
  typedef typename OpenMeshType::VertexIterator VertexIterator;
  typedef typename OpenMeshType::FaceIterator FaceIterator;
  typedef typename OpenMeshType::CoordType CoordType;
  
  class Info
  {
  public:
    
    Info() {}
    
    /// It returns a bit mask describing the field present in the obj file
    int mask=0;
    
    /// a Simple callback that can be used for long obj parsing.
    CallBackPos *cb=nullptr;
    
    int numVertices;
    int numEdges;
    int numFaces; // Note that numFaces can be different from the  final number of triangles
    int numTexCoords=0;
    int numNormals;
    
  }; // end class
  
  
  struct ObjIndexedFace
  {
    void set(const int & num){v.resize(num);n.resize(num); t.resize(num);}
    std::vector<int> v;
    std::vector<int> n;
    std::vector<int> t;
    int tInd;
    bool edge[3];// useless if the face is a polygon, no need to have variable length array
    Color4b c;
    int mInd;
  };
  
  struct ObjEdge
  {
    int v0;
    int v1;
  };
  
  struct ObjTexCoord
  {
    float u;
    float v;
  };
  
  enum OBJError {
    // Successfull opening
    E_NOERROR                           = 0*2+0,  //  A*2+B  (A position of correspondig string in the array, B=1 if not critical)
    
    // Non Critical Errors (only odd numbers)
    E_NON_CRITICAL_ERROR                = 0*2+1,
    E_MATERIAL_FILE_NOT_FOUND           = 1*2+1,
    E_MATERIAL_NOT_FOUND                = 2*2+1,
    E_TEXTURE_NOT_FOUND                 = 3*2+1,
    E_VERTICES_WITH_SAME_IDX_IN_FACE    = 4*2+1,
    E_LESS_THAN_3_VERT_IN_FACE          = 5*2+1,
    
    // Critical Opening Errors (only even numbers)
    E_CANTOPEN                          = 6*2+0,
    E_UNEXPECTED_EOF                    = 7*2+0,
    E_ABORTED                           = 8*2+0,
    E_NO_VERTEX                         = 9*2+0,
    E_NO_FACE                           =10*2+0,
    E_BAD_VERTEX_STATEMENT              =11*2+0,
    E_BAD_VERT_TEX_STATEMENT            =12*2+0,
    E_BAD_VERT_NORMAL_STATEMENT         =13*2+0,
    E_BAD_VERT_INDEX                    =14*2+0,
    E_BAD_VERT_TEX_INDEX                =15*2+0,
    E_BAD_VERT_NORMAL_INDEX             =16*2+0,
    E_LESS_THAN_4_VERT_IN_QUAD          =17*2+0
  };
  
  // to check if a given error is critical or not.
  static bool ErrorCritical(int err)
  {
    if (err==0) return false;
    if (err&1) return false;
    return true;
  }
  
  static const char* ErrorMsg(int error)
  {
    const int MAXST = 18;
    static const char* obj_error_msg[MAXST] =
    {
      /*  0 */ "No errors",
      
      /*  1 */ "Material library file wrong or not found, a default white material is used",
      /*  2 */ "Some materials definitions were not found, a default white material is used where no material was available",
      /*  3 */ "Texture file not found",
      /*  4 */ "Identical vertex indices found in the same faces -- faces ignored",
      /*  5 */ "Faces with fewer than 3 vertices  -- faces ignored",
      
      /*  6 */ "Can't open file",
      /*  7 */ "Premature End of File. File truncated?",
      /*  8 */ "Loading aborted by user",
      /*  9 */ "No vertex found",
      /* 10 */ "No face found",
      /* 11 */ "Vertex statement with fewer than 3 coords",
      /* 12 */ "Texture coords statement with fewer than 2 coords",
      /* 13 */ "Vertex normal statement with fewer than 3 coords",
      /* 14 */ "Bad vertex index in face",
      /* 15 */ "Bad texture coords index in face",
      /* 16 */ "Bad vertex normal index in face",
      /* 17 */ "Quad faces with number of corners different from 4"
    };
    
    error >>= 1;
    
    if( (error>=MAXST) || (error<0) ) return "Unknown error";
    else return obj_error_msg[error];
  }
  
  // Helper functions that checks the range of indexes
  // putting them in the correct range if less than zero (as in the obj style)
  
  static bool GoodObjIndex(int &index, const int maxVal)
  {
    if (index > maxVal)	return false;
    if (index < 0)
    {
      index += maxVal+1;
      if (index<0 || index > maxVal)	return false;
    }
    return true;
  }
  
  static int Open(OpenMeshType &mesh, const char *filename, int &loadmask, CallBackPos *cb=0)
  {
    Info oi;
    oi.mask=0;
    oi.cb=cb;
    int ret=Open(mesh,filename,oi);
    loadmask=oi.mask;
    return ret;
  }
  
  /*!
  * Opens an object file (in ascii format) and populates the mesh passed as first
  * accordingly to read data
  * \param m The mesh model to be populated with data stored into the file
  * \param filename The name of the file to be opened
  * \param oi A structure containing infos about the object to be opened
  */
  static int Open( OpenMeshType &m, const char * filename, Info &oi)
  {
    int result = E_NOERROR;
    
    m.Clear();
    CallBackPos *cb = oi.cb;
    
    // if LoadMask has not been called yet, we call it here
    if (oi.mask == 0)
      LoadMask(filename, oi);
    
    const int inputMask = oi.mask;
    Mask::ClampMask<OpenMeshType>(m,oi.mask);
    
    if (oi.numVertices == 0)
      return E_NO_VERTEX;
          
    std::ifstream stream(filename);
    if (stream.fail())
    {
      stream.close();
      return E_CANTOPEN;
    }
    
    typename OpenMeshType::template PerMeshAttributeHandle<std::vector<Material> > materialsHandle =
        vcg::tri::Allocator<OpenMeshType>:: template GetPerMeshAttribute<std::vector<Material> >(m, std::string("materialVector"));
    typename OpenMeshType::template PerFaceAttributeHandle<int> mIndHandle =
        vcg::tri::Allocator<OpenMeshType>:: template GetPerFaceAttribute<int>(m, std::string("materialIndex"));
    std::vector<Material>&	materials = materialsHandle();  // materials vector
    std::vector<ObjTexCoord>	texCoords;  // texture coordinates
    std::vector<CoordType>  normals;		// vertex normals
    std::vector<ObjIndexedFace> indexedFaces;
    std::vector< std::string > tokens;
    std::string line;
    std::string	header;
    
    short currentMaterialIdx = 0;			// index of current material into materials vector
    Color4b currentColor=Color4b::LightGray;	// we declare this outside code block since other
    // triangles of this face will share the same color
    
    Material defaultMaterial;					// default material: white
    defaultMaterial.index=currentMaterialIdx;
    materials.push_back(defaultMaterial);
    
    int numVertices  = 0;  // stores the number of vertices been read till now
    int numEdges     = 0;  // stores the number of edges read till now
    int numTriangles = 0;  // stores the number of faces been read till now
    int numTexCoords = 0;  // stores the number of texture coordinates been read till now
    int numVNormals	 = 0;  // stores the number of vertex normals been read till now
    
    int numVerticesPlusFaces = oi.numVertices + oi.numFaces;
    int extraTriangles=0;
    // vertices and faces allocation
    VertexIterator vi = vcg::tri::Allocator<OpenMeshType>::AddVertices(m,oi.numVertices);
    //FaceIterator   fi = Allocator<OpenMeshType>::AddFaces(m,oi.numFaces);
    // edges found
    std::vector<ObjEdge> ev;
    std::vector<Color4b> vertexColorVector;
    ObjIndexedFace	ff;
    const char *loadingStr = "Loading";
    while (!stream.eof())
    {
      tokens.clear();
      TokenizeNextLine(stream, tokens, line, &vertexColorVector);
      
      unsigned int numTokens = static_cast<unsigned int>(tokens.size());
      if (numTokens > 0)
      {
        header.clear();
        header = tokens[0];
        
        // callback invocation, abort loading process if the call returns false
        if ((cb !=NULL) && (((numTriangles + numVertices)%100)==0) && !(*cb)((100*(numTriangles + numVertices))/numVerticesPlusFaces, loadingStr))
        {
          stream.close();
          return E_ABORTED;
        }
        if (header.compare("v")==0)	// vertex
        {
          loadingStr="Vertex Loading";
          if (numTokens < 4)
          {
            stream.close();
            return E_BAD_VERTEX_STATEMENT;
          }
          (*vi).P()[0] = (ScalarType) atof(tokens[1].c_str());
          (*vi).P()[1] = (ScalarType) atof(tokens[2].c_str());
          (*vi).P()[2] = (ScalarType) atof(tokens[3].c_str());
          ++numVertices;
          
          // assigning vertex color
          // ----------------------
          if (((oi.mask & vcg::tri::io::Mask::IOM_VERTCOLOR) != 0) && (HasPerVertexColor(m)))
          {
            if(numTokens>=7)
            {
              ScalarType rf(atof(tokens[4].c_str())), gf(atof(tokens[5].c_str())), bf(atof(tokens[6].c_str()));
              ScalarType scaling = (rf<=1 && gf<=1 && bf<=1) ? 255. : 1;
              
              unsigned char r			= (unsigned char) ((ScalarType) atof(tokens[4].c_str()) * scaling);
              unsigned char g			= (unsigned char) ((ScalarType) atof(tokens[5].c_str()) * scaling);
              unsigned char b			= (unsigned char) ((ScalarType) atof(tokens[6].c_str()) * scaling);
              unsigned char alpha = (unsigned char) ((numTokens>=8 ? (ScalarType) atof(tokens[7].c_str()) : 1)  * scaling);
              (*vi).C() = Color4b(r, g, b, alpha);
            }
            else
            {
              (*vi).C() = currentColor;
            }
          }
          
          ++vi;  // move to next vertex iterator
        }
        else if (header.compare("vt")==0)	// vertex texture coords
        {
          loadingStr="Vertex Texture Loading";
          
          if (numTokens < 3)
          {
            stream.close();
            return E_BAD_VERT_TEX_STATEMENT;
          }
          ObjTexCoord t;
          t.u = static_cast<float>(atof(tokens[1].c_str()));
          t.v = static_cast<float>(atof(tokens[2].c_str()));
          texCoords.push_back(t);
          
          numTexCoords++;
        }
        else if (header.compare("vn")==0)  // vertex normal
        {
          loadingStr="Vertex Normal Loading";
          
          if (numTokens != 4)
          {
            stream.close();
            return E_BAD_VERT_NORMAL_STATEMENT;
          }
          CoordType n;
          n[0] = (ScalarType) atof(tokens[1].c_str());
          n[1] = (ScalarType) atof(tokens[2].c_str());
          n[2] = (ScalarType) atof(tokens[3].c_str());
          normals.push_back(n);
          
          numVNormals++;
        }
        else if ( header.compare("l")==0 )
        {
          loadingStr = "Edge Loading";
          
          if (numTokens < 3)
          {
            result = E_LESS_THAN_3_VERT_IN_FACE; // TODO add proper/handling error code
            continue;
          }
          
          ObjEdge e = { (atoi(tokens[1].c_str()) - 1),
                        (atoi(tokens[2].c_str()) - 1) };
          ev.push_back(e);
          
          numEdges++;
        }
        else if( (header.compare("f")==0) || (header.compare("q")==0) )  // face
        {
          loadingStr="Face Loading";
          
          int vertexesPerFace = static_cast<int>(tokens.size()-1);
          
          bool QuadFlag = false; // QOBJ format by Silva et al for simply storing quadrangular meshes.
          if(header.compare("q")==0) {
            QuadFlag=true;
            if (vertexesPerFace != 4) {
              stream.close();
              return E_LESS_THAN_4_VERT_IN_QUAD;
            }
          }
          
          
          if (vertexesPerFace < 3) {
            // face with fewer than 3 vertices found: ignore this face
            extraTriangles--;
            result = E_LESS_THAN_3_VERT_IN_FACE;
            continue;
          }
          
          
          if( (vertexesPerFace>3) && OpenMeshType::FaceType::HasPolyInfo() )
          {
            //_BEGIN___ if you are filling a vcg mesh with GENERIC POLYGON 
            ff.set(vertexesPerFace);
            for(int i=0;i<vertexesPerFace;++i) { // remember index starts from 1 instead of 0
              SplitToken(tokens[i+1], ff.v[i], ff.n[i], ff.t[i], inputMask);
              if(QuadFlag) ff.v[i]++; // NOTE THAT THE STUPID QOBJ FORMAT IS ZERO INDEXED!!!!
            }
            if ( oi.mask & vcg::tri::io::Mask::IOM_WEDGTEXCOORD )
            {
              // verifying validity of texture coords indices
              for(int i=0;i<vertexesPerFace;i++)
                if(!GoodObjIndex(ff.t[i],oi.numTexCoords))
                {
                  stream.close();
                  return E_BAD_VERT_TEX_INDEX;
                }
              ff.tInd=materials[currentMaterialIdx].index;
            }
            
            // verifying validity of vertex indices
            std::vector<int> tmp = ff.v;
            std::sort(tmp.begin(),tmp.end());
            std::unique(tmp.begin(),tmp.end());
            if(tmp.size() != ff.v.size()) {
              result = E_VERTICES_WITH_SAME_IDX_IN_FACE;
              extraTriangles--;
              continue;
            }
            
            for(int i=0;i<vertexesPerFace;i++)
              if(!GoodObjIndex(ff.v[i],numVertices))
              {
                stream.close();
                return E_BAD_VERT_INDEX;
              }
            
            if(( oi.mask & vcg::tri::io::Mask::IOM_WEDGNORMAL ) ||
               ( oi.mask & vcg::tri::io::Mask::IOM_VERTNORMAL  ) )
            {
              // verifying validity of vertex normal indices
              for(int i=0;i<vertexesPerFace;i++)
                if(!GoodObjIndex(ff.n[i],numVNormals))
                {
                  stream.close();
                  return E_BAD_VERT_NORMAL_INDEX;
                }
            }
            
            
            if(((oi.mask & vcg::tri::io::Mask::IOM_FACECOLOR) != 0) && HasPerFaceColor(m)) // assigning face color
              ff.c = currentColor;
            
            ++numTriangles;
            indexedFaces.push_back(ff);
            
            //_END  ___ if you are filling a vcg mesh with GENERIC POLYGON 
          }
          else
          {
            //_BEGIN___ if you are filling a vcg mesh with TRIANGLES 
            std::vector<std::vector<vcg::Point3f> > polygonVect(1); // it is a vector of polygon loops
            polygonVect[0].resize(vertexesPerFace);
            std::vector<int> indexVVect(vertexesPerFace);
            std::vector<int> indexNVect(vertexesPerFace);
            std::vector<int> indexTVect(vertexesPerFace);
            std::vector<int> indexTriangulatedVect;
            
            for(int pi=0;pi<vertexesPerFace;++pi)
            {
              SplitToken(tokens[pi+1], indexVVect[pi],indexNVect[pi],indexTVect[pi], inputMask);
              if(QuadFlag) indexVVect[pi]++; // NOTE THAT THE STUPID QOBJ FORMAT IS ZERO INDEXED!!!!
              GoodObjIndex(indexVVect[pi],numVertices);
              GoodObjIndex(indexTVect[pi],oi.numTexCoords);
              polygonVect[0][pi].Import(m.vert[indexVVect[pi]].cP());
            }
            if(vertexesPerFace>3)
              oi.mask |= Mask::IOM_BITPOLYGONAL;
            
            if(vertexesPerFace<5)
              FanTessellator(polygonVect, indexTriangulatedVect);
            else
            {
#ifdef __gl_h_
              //qDebug("OK: using opengl tessellation for a polygon of %i verteces",vertexesPerFace);
              vcg::glu_tesselator::tesselate<vcg::Point3f>(polygonVect, indexTriangulatedVect);
              if(indexTriangulatedVect.size()==0)
                FanTessellator(polygonVect, indexTriangulatedVect);
#else
              //qDebug("Warning: using fan tessellation for a polygon of %i verteces",vertexesPerFace);
              FanTessellator(polygonVect, indexTriangulatedVect);
#endif
            }
            extraTriangles+=((indexTriangulatedVect.size()/3) -1);
#ifdef QT_VERSION
            if( int(indexTriangulatedVect.size()/3) != vertexesPerFace-2)
            {
              qDebug("Warning there is a degenerate poligon of %i verteces that was triangulated into %i triangles",vertexesPerFace,int(indexTriangulatedVect.size()/3));
              for(size_t qq=0;qq<polygonVect[0].size();++qq)
                qDebug("      (%f %f %f)",polygonVect[0][qq][0],polygonVect[0][qq][1],polygonVect[0][qq][2]);
              for(size_t qq=0;qq<tokens.size();++qq) qDebug("<%s>",tokens[qq].c_str());
            }
#endif
            //qDebug("Triangulated a face of %i vertexes into %i triangles",polygonVect[0].size(),indexTriangulatedVect.size());
            
            for(size_t pi=0;pi<indexTriangulatedVect.size();pi+=3)
            {
              ff.set(3);
              int locInd[3];
              for(int iii=0;iii<3;++iii)
              {
                locInd[iii]=indexTriangulatedVect[pi+iii];
                ff.v[iii]=indexVVect[ locInd[iii] ];
                ff.n[iii]=indexNVect[ locInd[iii] ];
                ff.t[iii]=indexTVect[ locInd[iii] ];
              }
              
              // Setting internal edges: only edges formed by consecutive edges are external.
              for(int iii=0;iii<3;++iii)
              {
                if( (locInd[iii]+1)%vertexesPerFace == locInd[(iii+1)%3]) ff.edge[iii]=false;
                else ff.edge[iii]=true;
              }
              
              if ( oi.mask & vcg::tri::io::Mask::IOM_WEDGTEXCOORD )
              { // verifying validity of texture coords indices
                bool invalid = false;
                for(int i=0;i<3;i++)
                  if(!GoodObjIndex(ff.t[i],oi.numTexCoords))
                  {
                    //return E_BAD_VERT_TEX_INDEX;
                    invalid = true;
                    break;
                  }
                if (invalid) continue;
                ff.tInd=materials[currentMaterialIdx].index;
              }
              
              // verifying validity of vertex indices
              if ((ff.v[0] == ff.v[1]) || (ff.v[0] == ff.v[2]) || (ff.v[1] == ff.v[2])) {
                result = E_VERTICES_WITH_SAME_IDX_IN_FACE;
                extraTriangles--;
                continue;
              }
              
              {
                bool invalid = false;
                for(int i=0;i<3;i++)
                  if(!GoodObjIndex(ff.v[i],numVertices))
                  {
                    //return E_BAD_VERT_INDEX;
                    invalid = true;
                    break;
                  }
                if (invalid) continue;
              }
              
              // assigning face normal
              if ( ( oi.mask & vcg::tri::io::Mask::IOM_WEDGNORMAL  ) ||
                   ( oi.mask & vcg::tri::io::Mask::IOM_VERTNORMAL  ) )
              {   // verifying validity of vertex normal indices
                bool invalid = false;
                for(int i=0;i<3;i++)
                  if(!GoodObjIndex(ff.n[i],numVNormals))
                  {
                    //return E_BAD_VERT_NORMAL_INDEX;
                    invalid = true;
                    break;
                  }
                if (invalid) continue;
              }
              
              // assigning face color
              if( ((oi.mask & vcg::tri::io::Mask::IOM_FACECOLOR) != 0) && HasPerFaceColor(m))
			  {
				  ff.c = currentColor;
			  }
              
              ff.mInd = currentMaterialIdx;
              
              ++numTriangles;
              indexedFaces.push_back(ff);
            }
            
          }  //_END  ___ if you are filling a vcg mesh with TRIANGLES
        }
        else if ((header.compare("mtllib")==0) && (tokens.size() > 1))	// material library
        {
          // obtain the name of the file containing materials library
          std::string materialFileName;
          if (tokens.size() == 2)
            materialFileName = tokens[1]; //play it safe
          else
            materialFileName = line.substr(7); //get everything after "mtllib "
          
          if (!LoadMaterials( materialFileName.c_str(), materials, m.textures))
            result = E_MATERIAL_FILE_NOT_FOUND;
        }
        else if ((header.compare("usemtl")==0) && (tokens.size() > 1))	// material usage
        {
          // emergency check. If there are no materials, the material library failed to load or was not specified
          // but there are tools that save the material library with the same name of the file, but do not add the 
          // "mtllib" definition in the header. So, we can try to see if this is the case
          if ((materials.size() == 1)&&(materials[0].materialName == "")){
            std::string materialFileName(filename);
            materialFileName.replace(materialFileName.end()-4, materialFileName.end(), ".mtl");
            LoadMaterials(materialFileName.c_str(), materials, m.textures);
          }
          
          std::string materialName;
          if (tokens.size() == 2)
            materialName = tokens[1]; //play it safe
          else
            materialName = line.substr(7); //get everything after "usemtl "
          
          bool found = false;
          unsigned i = 0;
          while (!found && (i < materials.size()))
          {
            std::string currentMaterialName = materials[i].materialName;
            if (currentMaterialName == materialName)
            {
              currentMaterialIdx = i;
              Material &material = materials[currentMaterialIdx];
              Point3f diffuseColor = material.Kd;
              unsigned char r			= (unsigned char) (diffuseColor[0] * 255.0);
              unsigned char g			= (unsigned char) (diffuseColor[1] * 255.0);
              unsigned char b			= (unsigned char) (diffuseColor[2] * 255.0);
              unsigned char alpha = (unsigned char) (material.Tr  * 255.0);
              currentColor= Color4b(r, g, b, alpha);
              found = true;
            }
            ++i;
          }
          
          if (!found)
          {
            currentMaterialIdx = 0;
            result = E_MATERIAL_NOT_FOUND;
          }
        }
        // we simply ignore other situations
      } // end for each line...
    } // end while stream not eof
    assert((numTriangles +numVertices) == numVerticesPlusFaces+extraTriangles);
    vcg::tri::Allocator<OpenMeshType>::AddFaces(m,numTriangles);
    
    // Add found edges
    if (numEdges > 0)
    {
      vcg::tri::Allocator<OpenMeshType>::AddEdges(m,numEdges);
      
      assert(m.edge.size() == size_t(m.en));
      
      for(int i=0; i<numEdges; ++i)
      {
        ObjEdge &  e    = ev[i];
        assert(e.v0 >= 0 && size_t(e.v0) < m.vert.size() &&
               e.v1 >= 0 && size_t(e.v1) < m.vert.size());
        // TODO add proper handling of bad indices
        m.edge[i].V(0) = &(m.vert[e.v0]);
        m.edge[i].V(1) = &(m.vert[e.v1]);
      }
    }
    //-------------------------------------------------------------------------------
    
    // Now the final passes:
    // First Pass to convert indexes into pointers for face to vert/norm/tex references
    for(int i=0; i<numTriangles; ++i)
    {
      assert(m.face.size() == size_t(m.fn));
      m.face[i].Alloc(indexedFaces[i].v.size()); // it does not do anything if it is a trimesh
      
      for(unsigned int j=0;j<indexedFaces[i].v.size();++j)
      {   
        int vertInd = indexedFaces[i].v[j];
        assert(vertInd >=0 && vertInd < m.vn); (void)vertInd;
        m.face[i].V(j) = &(m.vert[indexedFaces[i].v[j]]);
        
        if (((oi.mask & vcg::tri::io::Mask::IOM_WEDGTEXCOORD) != 0) && HasPerWedgeTexCoord(m))
        {
          ObjTexCoord t = texCoords[indexedFaces[i].t[j]];
          m.face[i].WT(j).u() = t.u;
          m.face[i].WT(j).v() = t.v;
          m.face[i].WT(j).n() = indexedFaces[i].tInd;
        }
        if (((oi.mask & vcg::tri::io::Mask::IOM_VERTTEXCOORD) != 0 ) && HasPerVertexTexCoord(m))
		{
          ObjTexCoord t = texCoords[indexedFaces[i].t[j]];
          m.face[i].V(j)->T().u() = t.u;
          m.face[i].V(j)->T().v() = t.v;
          m.face[i].V(j)->T().n() = indexedFaces[i].tInd;
        }
        if (((oi.mask & vcg::tri::io::Mask::IOM_WEDGNORMAL) != 0) && HasPerWedgeNormal(m))
        {
          m.face[i].WN(j).Import(normals[indexedFaces[i].n[j]]);
        }
        
        if (((oi.mask & vcg::tri::io::Mask::IOM_VERTNORMAL) != 0) && HasPerVertexNormal(m))
        {
          m.face[i].V(j)->N().Import(normals[indexedFaces[i].n[j]]);
        }
        
        // set faux edge flags according to internals faces
        if (indexedFaces[i].edge[j]) 
          m.face[i].SetF(j);
        else 
          m.face[i].ClearF(j);
      }
      
      if (HasPerFaceNormal(m))
      {
        if (((oi.mask & vcg::tri::io::Mask::IOM_FACECOLOR) != 0) && (HasPerFaceColor(m)))
        {
          m.face[i].C() = indexedFaces[i].c;
          mIndHandle[i] = indexedFaces[i].mInd;
        }
        
        if (((oi.mask & vcg::tri::io::Mask::IOM_WEDGNORMAL) != 0) && (HasPerWedgeNormal(m)))
        {
          // face normal is computed as an average of wedge normals
          m.face[i].N().Import(m.face[i].WN(0)+m.face[i].WN(1)+m.face[i].WN(2));
        }
        else
        {
          m.face[i].N().Import(TriangleNormal(m.face[i]).Normalize());
        }
      }
    }
    // final pass to manage the ZBrush PerVertex Color that are managed into comments
    if(vertexColorVector.size()>0)
    {
      //	  if(vertexColorVector.size()!=m.vn){
      //		qDebug("Warning Read %i vertices and %i vertex colors",m.vn,vertexColorVector.size());
      //		qDebug("line count %i x 64 = %i",MRGBLineCount(), MRGBLineCount()*64);
      //	  }
      for(int i=0;i<m.vn;++i)
      {
        m.vert[i].C()=vertexColorVector[i];
      }
    }
    stream.close();
    return result;
  } // end of Open
  
  
  /*!
  * Read the next valid line and parses it into "tokens" (e.g. groups like 234/234/234), allowing
  * the tokens to be read one at a time. It read multiple lines  concatenating them if they end with '\'
  *  \param stream  The object providing the input stream
  *  \param tokens  The "tokens" in the next line
  */
  inline static void TokenizeNextLine(std::ifstream &stream, std::vector< std::string > &tokens, std::string &line, std::vector<Color4b> *colVec)
  {
    if(stream.eof()) return;
    
    do
    {
      std::getline(stream, line);
      // We have to manage backspace terminated lines, 
      // joining them together before parsing them
      if(!line.empty() && line.back()==13) line.pop_back();
      while(!line.empty() && line.back()=='\\') {
        std::string tmpLine;
        std::getline(stream, tmpLine);
        if(tmpLine.back()==13) line.pop_back();
        line.pop_back(); 
        line.append(tmpLine);
      }
      const size_t len = line.length();
      if((len > 0) && colVec && line[0] == '#')
      {
        // The following MRGB block contains ZBrush Vertex Color (Polypaint)
        // and masking output as 4 hexadecimal values per vertex. The vertex color format is MMRRGGBB with up to 64 entries per MRGB line.
        if((len >= 5) && line[1] == 'M' && line[2] == 'R' && line[3] == 'G' && line[4] == 'B')
        { // Parsing the polycolor of ZBrush
          MRGBLineCount()++;
          char buf[3]="00";
          Color4b cc(Color4b::Black);
          for(size_t i=6;(i+7)<len;i+=8)
          {
            for(size_t j=1;j<4;j++)
            {
              buf[0]=line[i+j*2+0];
              buf[1]=line[i+j*2+1];
              buf[2]=0;
              char *p;
              int val=strtoul(buf,&p,16);
              cc[j-1]= val;
            }
            colVec->push_back(cc);
          }
        }
      }
    }
    while (( line.length()==0 || line[0] == '#') && !stream.eof());  // skip comments and empty lines
    
    if ( (line.length() == 0)||(line[0] == '#') )  // can be true only on last line of file
      return;
    
    size_t from		= 0;
    size_t to			= 0;
    size_t length = line.size();
    
    tokens.clear();
    do
    {
      while (from!=length && (line[from]==' ' || line[from]=='\t' || line[from]=='\r') )
        from++;
      if(from!=length)
      {
        to = from+1;
        while (to!=length && line[to]!=' ' && line[to] != '\t' && line[to]!='\r')
          to++;
        tokens.push_back(line.substr(from, to-from).c_str());
        from = to;
      }
    }
    while (from<length);
  } // end TokenizeNextLine
  
  // This function takes a token and, according to the mask, it returns the indexes of the involved vertex, normal and texcoord indexes.
  // Example. if the obj file has vertex texcoord (e.g. lines 'vt 0.444 0.5555')
  // when parsing  a line like
  // f 46/303 619/325 624/326 623/327
  // if in the mask you have specified to read wedge tex coord
  // for the first token it will return inside vId and tId the corresponding indexes 46 and 303 )                
  inline static void SplitToken(const std::string & token, int & vId, int & nId, int & tId, int mask)
  {
    static const char delimiter = '/';
    
    vId = nId = tId = 0;
    if (token.empty()) return;
    
    size_t firstSep  = token.find_first_of(delimiter);
    size_t secondSep = (firstSep == std::string::npos) ? (std::string::npos) : (token.find_first_of(delimiter, firstSep + 1));
    
    const bool hasPosition = true;
    const bool hasTexcoord = (firstSep  != std::string::npos) && ((firstSep + 1) < secondSep);
    const bool hasNormal   = (secondSep != std::string::npos) || (mask & Mask::IOM_WEDGNORMAL) || (mask & Mask::IOM_VERTNORMAL);
    
    if (hasPosition) vId = atoi(token.substr(0, firstSep).c_str()) - 1;
    if (hasTexcoord) tId = atoi(token.substr(firstSep + 1, secondSep - firstSep - 1).c_str()) - 1;
    if (hasNormal)
      nId = atoi(token.substr(secondSep + 1).c_str()) - 1;
  }
  /** returns a Point3f done from (tokens[pos],tokens[pos+1],tokens[pos+2])
    
    */
  static Point3f Point3fFrom3Tokens(std::vector< std::string > &tokens, int pos)
  {
    float r = (float) atof(tokens[pos+0].c_str());
    float g = (float) atof(tokens[pos+1].c_str());
    float b = (float) atof(tokens[pos+2].c_str());
    return Point3f(r, g, b);
  }
  
  /*!
  * Retrieves infos about kind of data stored into the file and fills a mask appropriately
  * \param filename The name of the file to open
  * \param mask     A mask which will be filled according to type of data found in the object
  * \param oi       A structure which will be filled with infos about the object to be opened
  */
  
  static bool LoadMask(const char * filename, Info &oi)
  {
    
    std::ifstream stream(filename);
    if (stream.fail())
    {
      stream.close();
      return false;
    }
    // obtain length of file:
    stream.seekg (0, std::ios::end);
    int length = stream.tellg();
    stream.seekg (0, std::ios::beg);
    
    if (length == 0) return false;
    
    bool bHasPerFaceColor = false;
    bool bHasNormals     = false;
    bool bHasPerVertexColor = false;
    
    oi.numVertices=0;
    oi.numEdges=0;
    oi.numFaces=0;
    oi.numTexCoords=0;
    oi.numNormals=0;
    int lineCount=0;
    int totRead=0;
    bool firstV = true;
    std::string line;
    while (!stream.eof())
    {
      lineCount++;
      std::getline(stream, line);
      totRead+=line.size();
      if(oi.cb && (lineCount%1000)==0)
        (*oi.cb)( (int)(100.0*(float(totRead))/float(length)), "Loading mask...");
      if(line.size()>2)
      {
        if(line[0]=='v')
        {
          if ((line[1] == ' ') || (line[1] == '\t'))
          {
            oi.numVertices++;
            if (firstV)
            {
              int sepN = 0;
              for (size_t lit = 0; lit < line.size(); lit++){
                if ((line[lit] == ' ') || (line[lit] == '\t'))
                  sepN++;
              }
              if (sepN >= 6)
                bHasPerVertexColor = true;
              firstV = false;
            }
          }
          if(line[1]=='t') oi.numTexCoords++;
          if(line[1]=='n') {
            oi.numNormals ++;
            bHasNormals = true;
          }
        }
        else {
          if((line[0]=='f') || (line[0]=='q')) oi.numFaces++;
          else
            if (line[0]=='l') oi.numEdges++;
            else
              if(line[0]=='u' && line[1]=='s') bHasPerFaceColor = true; // there is a usematerial so add per face color
        }
      }
    }
    oi.mask = 0;
    if (oi.numTexCoords)
    {
      if (oi.numTexCoords==oi.numVertices)
        oi.mask |= vcg::tri::io::Mask::IOM_VERTTEXCOORD;
      
      oi.mask |= vcg::tri::io::Mask::IOM_WEDGTEXCOORD;
      // Usually if you have tex coords you also have materials
      oi.mask |= vcg::tri::io::Mask::IOM_FACECOLOR;
    }
    if(bHasPerFaceColor)  oi.mask |= vcg::tri::io::Mask::IOM_FACECOLOR;
    if(bHasPerVertexColor)  oi.mask |= vcg::tri::io::Mask::IOM_VERTCOLOR;
    if (bHasNormals) {
      if (oi.numNormals == oi.numVertices)
        oi.mask |= vcg::tri::io::Mask::IOM_VERTNORMAL;
      else
        oi.mask |= vcg::tri::io::Mask::IOM_WEDGNORMAL;
    }
    if (oi.numEdges)
      oi.mask |= vcg::tri::io::Mask::IOM_EDGEINDEX;
    
    stream.close();
    
    return true;
  }
  
  static bool LoadMask(const char * filename, int &mask)
  {
    Info oi;
    bool ret=LoadMask(filename, oi);
    mask= oi.mask;
    return ret;
  }
  
	static bool LoadMaterials(const char * filename, std::vector<Material> &materials, std::vector<std::string> &textures)
	{
		// assumes we are in the right directory

		std::ifstream stream(filename);
		if (stream.fail())
			return false;

		std::vector< std::string > tokens;
		std::string line;
		std::string	header;

		materials.clear();
		Material currentMaterial;

		// Fill in some default values for the material
		currentMaterial.index = (unsigned int)(-1);
		currentMaterial.Ka = Point3f(0.2, 0.2, 0.2);
		currentMaterial.Kd = Point3f(1, 1, 1);
		currentMaterial.Ks = Point3f(1, 1, 1);
		currentMaterial.Tr = 1;
		currentMaterial.Ns = 0;
		currentMaterial.illum = 2;

		bool first = true;
		while (!stream.eof()) {
			tokens.clear();
			TokenizeNextLine(stream, tokens, line, 0);

			if (tokens.size() > 0) {
				header.clear();
				header = tokens[0];

				if (header.compare("newmtl")==0) {
					if (!first) {
						materials.push_back(currentMaterial);
						currentMaterial = Material();
						currentMaterial.index = (unsigned int)(-1);
					}
					else
						first = false;

					if(tokens.size() < 2)
						return false;
					else if (tokens.size() == 2)
						currentMaterial.materialName = tokens[1]; //play it safe
					else
						currentMaterial.materialName = line.substr(7); //space in the name, get everything after "newmtl "
				}
				else if (header.compare("Ka")==0) {
					if (tokens.size() >= 4) {
						currentMaterial.Ka = Point3fFrom3Tokens(tokens,1);
					}
				}
				else if (header.compare("Kd")==0) {
					if (tokens.size() >= 4) {
						currentMaterial.Kd = Point3fFrom3Tokens(tokens,1);
					}
				}
				else if (header.compare("Ks")==0) {
					if (tokens.size() >= 4) {
						currentMaterial.Ks = Point3fFrom3Tokens(tokens,1);
					}
				}
				else if ((header.compare("d")==0) || (header.compare("Tr")==0)) { // alpha
					if (tokens.size() < 2) {
						currentMaterial.Tr = (float) atof(tokens[1].c_str());
					}
				}
				else if (header.compare("Ns")==0) { // shininess
					if (tokens.size() < 2) {
						currentMaterial.Ns = float(atoi(tokens[1].c_str()));
					}
				}
				else if (header.compare("illum")==0) { // specular illumination on/off
					if (tokens.size() < 2) {
						currentMaterial.illum = atoi(tokens[1].c_str());
					}
				}
				else if(header.compare("map_Kd")==0) { // texture name
					std::string textureName;
					if (tokens.size() == 2) {
						//the tex name is the last one (after any option)
						textureName = tokens[tokens.size()-1];
					}

					currentMaterial.map_Kd=textureName;

					// adding texture name into textures vector (if not already present)
					// avoid adding the same name twice
					auto it = std::find(textures.begin(), textures.end(), textureName);
					if(it==textures.end()) {
						currentMaterial.index = textures.size();
						textures.push_back(textureName);
					} else {
						 currentMaterial.index = std::distance(textures.begin(),it);
					}
				}
				// we simply ignore other situations
			}
		}
		materials.push_back(currentMaterial);  // add last read material

		stream.close();
		// Sometimes some materials have texture and no texture
		// in this case for sake of uniformity we just use the first texture.
		if(!textures.empty()) {
			for(size_t i=0;i<materials.size();++i) {
				if(materials[i].map_Kd.empty()) {
					materials[i].map_Kd=textures[0];
					materials[i].index=0;
				}
			}
		}

		return true;
	}

}; // end class
} // end Namespace tri
} // end Namespace io
} // end Namespace vcg

#endif  // ndef __VCGLIB_IMPORT_OBJ
