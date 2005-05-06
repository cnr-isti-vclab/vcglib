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

#ifndef __VCGLIB_IMPORT_RAW
#define __VCGLIB_IMPORT_RAW

#include <stdio.h>

namespace vcg {
namespace tri {
namespace io {

/** 
This class encapsulate a filter for importing raw format pointcloud.
there exists many raw formats. each one with a particular sintax even if they only contains 

*/
template <class MESH_TYPE>
class ImporterRAW
{
public:

typedef typename MESH_TYPE::VertexPointer VertexPointer;
typedef typename MESH_TYPE::ScalarType ScalarType;
typedef typename MESH_TYPE::VertexType VertexType;
typedef typename MESH_TYPE::FaceType FaceType;
typedef typename MESH_TYPE::VertexIterator VertexIterator;
typedef typename MESH_TYPE::FaceIterator FaceIterator;

// max token number
#define RAW_MAX_TOKEN_LINE_DESCRIPTOR 32


enum RAWError {
	E_NOERROR,				// 0
		// Errori di open
	E_CANTOPEN,				// 1
	E_UNESPECTEDEOF,        // 2
	    // errore line descriptor
	E_INVALIDLINEDESC       // 3
};

static const char *ErrorMsg(int error)
{
  static const char * raw_error_msg[] =
  {
	"No errors",
	"Can't open file",
	"Premature End of file",
	"Invalid line Descriptor",
	};

  if(error>2 || error<0) return "Unknown error";
  else return stl_error_msg[error];
};

// line format is a string describing which data is stored for every data line
// PX PY PZ   posizione
// NX NY NZ   normale
// CR CG CB   colore
// RF         riflettanza (qualita')
//
// the string is parsed to know how many value are contained in each line
// and which is the order. the result is a number (how many) and a vector 
// describing the order
//
//
// during reading a data structure is used to store intermediate values
// it is basically an array of float 
//
// float linebuffer[]
//[0][1][2][3][4][5][6][7][8][9]
// p  p  p  n  n  n	 c  c  c  r
// x  y  z  x  y  z	 r  g  b  f
//
// given the number of tokens and the order vector it is possible to scan a line using the command
//
// for(...n 0->tokennumber...)
//   fscanf(fp,"%f", &linebuffer[tokenorder[n]])

static int Parselinedescription(const char * linedesc, int &tokennumber, int *order)
{
 int ii;
 char tok[3];
 int index;

 // controllo lunghezza
 // se non e' multiplo di 3 allora e' errato
 int len = strlen(linedesc) + 1;
  if(len%3 != 0)
   return E_INVALIDLINEDESC;

 index=0;
 tok[2] = '\0';
 tokennumber = 0;
 for(ii=0; ii<RAW_MAX_TOKEN_LINE_DESCRIPTOR; ii++)
   order[ii] = -1;

 while(index <= (len-3))
 {
  tok[0] = linedesc[index  ];
  tok[1] = linedesc[index+1];

  if(strcmp(tok,"PX") == 0)						// pos  x
    {order[tokennumber] = 0; tokennumber++;}
  else if (strcmp(tok,"PY") == 0)				// pos  y
    {order[tokennumber] = 1; tokennumber++;}
  else if (strcmp(tok,"PZ") == 0)				// pos  z
    {order[tokennumber] = 2; tokennumber++;}
  else if (strcmp(tok,"NX") == 0)				// norm x
    {order[tokennumber] = 3; tokennumber++;}
  else if (strcmp(tok,"NY") == 0)				// norm y
    {order[tokennumber] = 4; tokennumber++;}
  else if (strcmp(tok,"NZ") == 0)				// norm z
    {order[tokennumber] = 5; tokennumber++;}
  else if (strcmp(tok,"CR") == 0)				// col  r
    {order[tokennumber] = 6; tokennumber++;}
  else if (strcmp(tok,"CG") == 0)				// col  g
    {order[tokennumber] = 7; tokennumber++;}
  else if (strcmp(tok,"CB") == 0)				// col  b
    {order[tokennumber] = 8; tokennumber++;}
  else if (strcmp(tok,"RF") == 0)				// rifl
    {order[tokennumber] = 9; tokennumber++;}
  else											// nessuno dei suddetti... errore
  { return E_INVALIDLINEDESC; }

  index +=3;
 }

 return E_NOERROR;
}


// function to skip a line
static void Skipline(FILE *fp)
{
 char buf;

 fread(&(buf),sizeof(char),1,fp);

 while(buf != '\n')
  fread(&(buf),sizeof(char),1,fp);
}



/*!
*	Standard call for reading a mesh
*	\param m			the destination mesh
*	\param filename		the name of the file to read from
*	\param triangulate	if true, the mesh will be triangulated, otherwise only points will be stored
*	\param lineskip  	number of lines to be skipped at the begin of the file
*	\return				the operation result
*/
static int Open( MESH_TYPE &m, const char * filename, bool triangulate=false, int lineskip = 0, const char * linedesc = "PX PY PZ")
{
	if(triangulate)
		return Opentriang( m, filename, lineskip, linedesc);
	else
		return Openpoints( m, filename, lineskip, linedesc);
}


// in raw files with indication of rows and columns
static int Opentriang( MESH_TYPE &m, const char * filename, int lineskip = 0, const char * linedesc = "PX PY PZ")
{
  return E_NOERROR;
}

// generic raw reader, only points are imported
static int Openpoints( MESH_TYPE &m, const char * filename, int lineskip = 0, const char * linedesc = "PX PY PZ")
{
  int ii;
  int ret;
  FILE *fp;

  // line description
  int   tokennumber;
  int   tokenorder[RAW_MAX_TOKEN_LINE_DESCRIPTOR];

  //line data buffer
  float linebuffer[10];
  // fill buffer with standard values
  linebuffer[0] = linebuffer[1] = linebuffer[2] = 0.0;
  linebuffer[3] = linebuffer[4] = linebuffer[5] = 1.0;
  linebuffer[6] = 0.0; linebuffer[7] = 1.0; linebuffer[8] = 0.0;
  linebuffer[9] = 1.0;


  fp = fopen(filename, "r");
  if(fp == NULL)
  {
   return E_CANTOPEN;
  }

  // skip initial lines
  for(ii=0; ii<lineskip; ii++)
   Skipline(fp);

  // parsing line description
  ret = Parselinedescription(linedesc, tokennumber, tokenorder);
  if(ret)
   return ret;

  m.Clear();
//  FaceIterator fi=Allocator<OpenMeshType>::AddFaces(m,facenum);
//  VertexIterator vi=Allocator<OpenMeshType>::AddVertices(m,facenum*3);

  while(!feof(fp))
  {
   for(ii=0; ii<tokennumber; ii++)
     fscanf(fp,"%f", &(linebuffer[tokenorder[ii]]));

   // new vertex
   VertexType nv;
   //nv.Supervisor_Flags() = 0;

   // store the position
   nv.P()[0] = linebuffer[0];   nv.P()[1] = linebuffer[1];   nv.P()[2] = linebuffer[2];
   // store the normal
   if(m.HasPerVertexNormal())
   {
    nv.N()[0] = linebuffer[3];    nv.N()[1] = linebuffer[4];    nv.N()[2] = linebuffer[5];
   }

   // store the color
   if(m.HasPerVertexColor())
   {
    nv.C()[0] = linebuffer[6];    nv.C()[1] = linebuffer[7];    nv.C()[2] = linebuffer[8];
   }

   // store the reflectance
   if(m.HasPerVertexQuality())
   {
    nv.Q() = linebuffer[9];
   }

   m.vert.push_back(nv);
  }

  // update model point number
  m.vn = m.vert.size();

  fclose(fp);
  return E_NOERROR;
}

}; // end class
} // end Namespace tri
} // end Namespace io
} // end Namespace vcg

#endif
