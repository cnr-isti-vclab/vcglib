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
Revision 1.8  2004/10/28 00:52:45  cignoni
Better Doxygen documentation

Revision 1.7  2004/08/25 15:27:51  ponchio
Comma at end of enum.

Revision 1.6  2004/08/25 15:15:27  ganovelli
minor changes to comply gcc compiler (typename's and stuff)

Revision 1.5  2004/06/23 15:36:44  cignoni
Restructured management of error, now the standard open for any mesh type return the error code, the default success value is zero
Any import class has a method ErrorMsg that give a verbal description of an error code.

Revision 1.4  2004/03/18 15:30:57  cignoni
Removed float/double warning

Revision 1.3  2004/03/12 21:42:52  cignoni
First working version!

****************************************************************************/

#ifndef __VCGLIB_IMPORT_STL
#define __VCGLIB_IMPORT_STL

#include <stdio.h>
#include <wrap/callback.h>
#include <vcg/complex/trimesh/allocate.h>

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
	E_NOERROR,				// 0
		// Errori di open
	E_CANTOPEN,				// 1
	E_UNESPECTEDEOF       		        // 2
};

static const char *ErrorMsg(int error)
{
  static const char * stl_error_msg[] =
  {
	"No errors",
	"Can't open file",
	"Premature End of file",
	};

  if(error>2 || error<0) return "Unknown error";
  else return stl_error_msg[error];
};

static int Open( OpenMeshType &m, const char * filename, CallBackPos *cb=0)
{
  FILE *fp;
  bool binary=false;
  fp = fopen(filename, "r");
  if(fp == NULL)
    {
      return E_CANTOPEN;
    }

  /* Find size of file */
  fseek(fp, 0, SEEK_END);
  int file_size = ftell(fp);
  int facenum;
  /* Check for binary or ASCII file */
  fseek(fp, STL_LABEL_SIZE, SEEK_SET);
  fread(&facenum, sizeof(int), 1, fp);
  int expected_file_size=STL_LABEL_SIZE + 4 + (sizeof(short)+sizeof(STLFacet) )*facenum ;
  if(file_size ==  expected_file_size) binary = true;
  unsigned char tmpbuf[128];
  fread(tmpbuf,sizeof(tmpbuf),1,fp);
  for(int i = 0; i < sizeof(tmpbuf); i++)
    {
      if(tmpbuf[i] > 127)
 	      {
	        binary=true;
	        break;
	      }
    }
  // Now we know if the stl file is ascii or binary.
  fclose(fp);
  if(binary) return OpenBinary(m,filename,cb);
  else return OpenAscii(m,filename,cb);
}

static int OpenBinary( OpenMeshType &m, const char * filename, CallBackPos *cb=0)
{
  FILE *fp;
  fp = fopen(filename, "rb");
  if(fp == NULL)
  {
    return E_CANTOPEN;
  }
   
  int facenum;
  fseek(fp, STL_LABEL_SIZE, SEEK_SET);
  fread(&facenum, sizeof(int), 1, fp);
  
  m.Clear();
  FaceIterator fi=Allocator<OpenMeshType>::AddFaces(m,facenum);
  VertexIterator vi=Allocator<OpenMeshType>::AddVertices(m,facenum*3);
  // For each triangle read the normal, the three coords and a short set to zero
	for(int i=0;i<facenum;++i)
    {
      short attr;
      Point3f norm;
      Point3f tri[3];
      fread(&norm,sizeof(Point3f),1,fp);
      fread(&tri,sizeof(Point3f),3,fp);
      fread(&attr,sizeof(short),1,fp);
      for(int k=0;k<3;++k)
      {
        (*vi).P().Import(tri[k]); 
        (*fi).V(k)=&*vi; 
        ++vi;
      }
      ++fi;
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

    m.Clear();
  
    /* Skip the first line of the file */
    while(getc(fp) != '\n');

    STLFacet f;
    /* Read a single facet from an ASCII .STL file */
    while(!feof(fp))
    {
      fscanf(fp, "%*s %*s %f %f %f\n", &f.n.X(), &f.n.Y(), &f.n.Z());
      fscanf(fp, "%*s %*s");
      fscanf(fp, "%*s %f %f %f\n", &f.v[0].X(),  &f.v[0].Y(),  &f.v[0].Z());
      fscanf(fp, "%*s %f %f %f\n", &f.v[1].X(),  &f.v[1].Y(),  &f.v[1].Z());
      fscanf(fp, "%*s %f %f %f\n", &f.v[2].X(),  &f.v[2].Y(),  &f.v[2].Z());
      fscanf(fp, "%*s"); // end loop
      fscanf(fp, "%*s"); // end facet
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
