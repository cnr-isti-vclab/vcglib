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

#ifndef __VCGLIB_IMPORT_STL
#define __VCGLIB_IMPORT_STL
#include <stdio.h>
#include <algorithm>
#include <wrap/io_trimesh/io_mask.h>

namespace vcg {
namespace tri {
namespace io {

/**
This class encapsulate a filter for importing stl (stereolitograpy) meshes.
The stl format is quite simple and rather un-flexible. It just stores, in ascii or binary the, unindexed, geometry of the faces.
Warning: this code assume little endian (PC) architecture!!!
*/
template <class OpenMeshType>
class ImporterSTL
{
public:

typedef typename OpenMeshType::VertexPointer VertexPointer;
typedef typename OpenMeshType::ScalarType ScalarType;
typedef typename OpenMeshType::VertexType VertexType;
typedef typename OpenMeshType::FaceType FaceType;
typedef typename OpenMeshType::VertexIterator VertexIterator;
typedef typename OpenMeshType::FaceIterator FaceIterator;

// if it is binary there are 80 char of comment, the number fn of faces and then exactly fn*4*3 bytes.

enum {STL_LABEL_SIZE=80};

class STLFacet
{
public:
  Point3f n;
  Point3f v[3];
//  short attr;
};

enum STLError {
    E_NOERROR,       // 0
    E_CANTOPEN,      // 1
    E_UNESPECTEDEOF, // 2
    E_MALFORMED,     // 3 
  E_LAST
};

static const char *ErrorMsg(int error)
{
  static const char * stl_error_msg[] =
  {
    "No errors",
    "Can't open file",
    "Premature end of file",
    "Malformed file",
    };

  if(error>=E_LAST || error<0) return "Unknown error";
  else return stl_error_msg[error];
}

static bool LoadMask(const char * filename, int &mask)
{
  bool magicMode,colored;
  mask = Mask::IOM_VERTCOORD | Mask::IOM_FACEINDEX;
  if(!IsSTLColored(filename, colored, magicMode))
    return false;
   
  if(colored) mask |= Mask::IOM_FACECOLOR;
    return true;
}

/* Try to guess if a stl has color
 *
 * rules:
 * - It has to be binary
 * - The per face attribute should be not zero
 *
 * return false in case of malformed files
 */
static bool IsSTLColored(const char * filename, bool &coloredFlag, bool &magicsMode)
{
  coloredFlag=false;
  magicsMode=false;
  bool binaryFlag;
  if(IsSTLMalformed(filename,binaryFlag)==false)
    return false;
  
  if(binaryFlag==false)
     return true; 
   
   FILE *fp = fopen(filename, "rb");
   char buf[STL_LABEL_SIZE+1];
   fread(buf,sizeof(char),STL_LABEL_SIZE,fp);
   std::string strInput(buf);
   size_t cInd = strInput.rfind("COLOR=");
   size_t mInd = strInput.rfind("MATERIAL=");
   if(cInd!=std::string::npos && mInd!=std::string::npos)
     magicsMode = true;
   else
     magicsMode = false;
   int facenum;
   fread(&facenum, sizeof(int), 1, fp);

   for(int i=0;i<std::min(facenum,1000);++i)
   {
     unsigned short attr;
     Point3f norm;
     Point3f tri[3];
     fread(&norm,sizeof(Point3f),1,fp);
     fread(&tri,sizeof(Point3f),3,fp);
     fread(&attr,sizeof(unsigned short),1,fp);
     if(attr!=0)
     {
      if(Color4b::FromUnsignedR5G5B5(attr) != Color4b(Color4b::White))
        coloredFlag=true;
     }
   }

   fclose(fp);
   return true;
}

/*
 * return false in case of malformed files
 * Try to guess if a stl is in binary format, and sets
 * the binaryFlag accordingly
 */
static bool IsSTLMalformed(const char * filename, bool &binaryFlag)
{
  binaryFlag=false;
  FILE *fp = fopen(filename, "rb");
  /* Find size of file */
  fseek(fp, 0, SEEK_END);
  std::size_t file_size = ftell(fp);
  unsigned int facenum;
  /* Check for binary or ASCII file */
  int ret = fseek(fp, STL_LABEL_SIZE, SEEK_SET);
  if (ret != 0) return false;
  ret = fread(&facenum, sizeof(unsigned int), 1, fp);
  if (ret != 1) return false;
  
  std::size_t expected_file_size=STL_LABEL_SIZE + 4 + (sizeof(short)+sizeof(STLFacet) )*facenum ;
  if(file_size ==  expected_file_size) 
  {
    binaryFlag = true;
    fclose(fp);
    return true;
  }
  
  // second check, sometimes the size is a bit wrong, 
  // lets'make a test to check that we find only ascii stuff before assuming it is ascii
  unsigned char tmpbuf[1000];
  std::size_t byte_to_read = std::min(sizeof(tmpbuf), (size_t)file_size - 80);
  fread(tmpbuf, byte_to_read,1,fp);
  fclose(fp);
  for(std::size_t i = 0; i < byte_to_read; i++)
    {
      if(tmpbuf[i] > 127)
          {
            binaryFlag=true; 
			std::size_t diff = (file_size > expected_file_size) ? file_size-expected_file_size : expected_file_size-file_size;
            if(diff > file_size/20 )
              return false; // 
            break;
          }
    }
  // Now we know if the stl file is ascii or binary.
  return true;
}

static int Open( OpenMeshType &m, const char * filename, int &loadMask, CallBackPos *cb=0)
{
  FILE *fp = fopen(filename, "r");
  if(fp == NULL)
      return E_CANTOPEN;
  fclose(fp);
  loadMask |= Mask::IOM_VERTCOORD | Mask::IOM_FACEINDEX;
  bool binaryFlag;
  if(!IsSTLMalformed(filename,binaryFlag))
    return E_MALFORMED;
  
  if(binaryFlag) return OpenBinary(m,filename,loadMask,cb);
  else return OpenAscii(m,filename,cb);
}

static int OpenBinary( OpenMeshType &m, const char * filename, int &loadMask, CallBackPos *cb=0)
{
  FILE *fp;
  fp = fopen(filename, "rb");
  if(fp == NULL)
  {
    return E_CANTOPEN;
  }

  bool magicsMode,coloredFlag;
  if(!IsSTLColored(filename,coloredFlag, magicsMode))
    return E_MALFORMED;
  if(!coloredFlag) 
    loadMask = loadMask & (~Mask::IOM_FACECOLOR);

  int facenum;
  fseek(fp, STL_LABEL_SIZE, SEEK_SET);
  fread(&facenum, sizeof(int), 1, fp);

  m.Clear();
  FaceIterator fi=Allocator<OpenMeshType>::AddFaces(m,facenum);
  VertexIterator vi=Allocator<OpenMeshType>::AddVertices(m,facenum*3);
  // For each triangle read the normal, the three coords and a short set to zero
    for(int i=0;i<facenum;++i)
    {
      unsigned short attr;
      Point3f norm;
      Point3f tri[3];
      fread(&norm,sizeof(Point3f),1,fp);
      fread(&tri,sizeof(Point3f),3,fp);
      fread(&attr,sizeof(unsigned short),1,fp);
      if(tri::HasPerFaceColor(m) && (loadMask & Mask::IOM_FACECOLOR) )
      {
        if(magicsMode) (*fi).C()= Color4b::FromUnsignedR5G5B5(attr);
                  else (*fi).C()= Color4b::FromUnsignedB5G5R5(attr);
      }
      for(int k=0;k<3;++k)
      {
        (*vi).P().Import(tri[k]);
        (*fi).V(k)=&*vi;
        ++vi;
      }
      ++fi;
      if(cb && (i%1000)==0) cb((i*100)/facenum,"STL Mesh Loading");
    }
    fclose(fp);
    return E_NOERROR;
  }


  static int OpenAscii( OpenMeshType &m, const char * filename, CallBackPos *cb=0)
  {
    FILE *fp;
    fp = fopen(filename, "r");
    if(fp == NULL)
    {
      return E_CANTOPEN;
    }
        long currentPos = ftell(fp);
        fseek(fp,0L,SEEK_END);
        long fileLen = ftell(fp);
        fseek(fp,currentPos,SEEK_SET);

    m.Clear();

    /* Skip the first line of the file */
    while(getc(fp) != '\n') { }

    STLFacet f;
    int cnt=0;
        int lineCnt=0;
        int ret;
    /* Read a single facet from an ASCII .STL file */
    while(!feof(fp))
    {
      if(cb && (++cnt)%1000)   cb( int(double(ftell(fp))*100.0/fileLen), "STL Mesh Loading");
        ret=fscanf(fp, "%*s %*s %f %f %f\n", &f.n.X(), &f.n.Y(), &f.n.Z()); // --> "facet normal 0 0 0"
            if(ret!=3)
            {
                // we could be in the case of a multiple solid object, where after a endfaced instead of another facet we have to skip two lines:
                //     endloop
                //	 endfacet
                //endsolid     <- continue on ret==0 will skip this line
                //solid ascii  <- and this one.
                //   facet normal 0.000000e+000 7.700727e-001 -6.379562e-001
                lineCnt++;
                continue;
            }
      ret=fscanf(fp, "%*s %*s"); // --> "outer loop"
      ret=fscanf(fp, "%*s %f %f %f\n", &f.v[0].X(),  &f.v[0].Y(),  &f.v[0].Z()); // --> "vertex x y z"
            if(ret!=3)
                return E_UNESPECTEDEOF;
      ret=fscanf(fp, "%*s %f %f %f\n", &f.v[1].X(),  &f.v[1].Y(),  &f.v[1].Z()); // --> "vertex x y z"
            if(ret!=3)
                return E_UNESPECTEDEOF;
      ret=fscanf(fp, "%*s %f %f %f\n", &f.v[2].X(),  &f.v[2].Y(),  &f.v[2].Z()); // --> "vertex x y z"
            if(ret!=3)
                return E_UNESPECTEDEOF;
      ret=fscanf(fp, "%*s"); // --> "endloop"
      ret=fscanf(fp, "%*s"); // --> "endfacet"
            lineCnt+=7;
      if(feof(fp)) break;
      FaceIterator fi=Allocator<OpenMeshType>::AddFaces(m,1);
      VertexIterator vi=Allocator<OpenMeshType>::AddVertices(m,3);
      for(int k=0;k<3;++k)
      {
        (*vi).P().Import(f.v[k]);
        (*fi).V(k)=&*vi;
        ++vi;
      }
    }
    fclose(fp);
    return E_NOERROR;
  }
}; // end class
} // end Namespace tri
} // end Namespace io
} // end Namespace vcg

#endif
