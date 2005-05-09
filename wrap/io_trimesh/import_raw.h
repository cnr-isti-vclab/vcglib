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
Revision 1.1  2005/05/06 13:58:26  callieri
First working version (callieri)


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
		// Error open
	E_CANTOPEN,				// 1
	E_UNESPECTEDEOF,        // 2
	    // error line descriptor
	E_INVALIDLINEDESC,      // 3
	    // error line parsing
	E_LINEERROR,            // 4
        // wrong number of points
	E_WRONGPOINTNUM			// 5
};

static const char *ErrorMsg(int error)
{
  static const char * raw_error_msg[] =
  {
	"No errors",
	"Can't open file",
	"Premature End of file",
	"Invalid line Descriptor",
	"Error parsing a line",
	"Point number different from expected"
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

// function to parse the line read from the raw file
// all characters besides numerals,dot,minus,plus 
// if e is found between numbers it's ketpt (12e4)
static int Parseline(int tokennumber, int *tokenorder, char *rawline, float *linebuffer)
{
 int linelen;
 int ii;
 bool change;
 int foundtok;

 // length
 linelen = strlen(rawline);

 // cleaning the line
 for(ii=0;ii<(linelen-1);)
 {
   change = true;

   if(isdigit(rawline[ii]))		// is a number
    change = false;
   else if(rawline[ii]==' ')	// is a space
    change = false;
   else if((rawline[ii]=='-') || (rawline[ii]=='+') || (rawline[ii]=='-') || (rawline[ii]=='.'))	// is + - .
    change = false;
   else if(rawline[ii]=='e')	// is e... is perhaps an exponential ?
     {  
	   if((ii>0) || (ii<(linelen-2)))	// if it's at the begin or the end of line, it's not an exponential
	    if(isdigit(rawline[ii-1]))		// a number before it
		  if(isdigit(rawline[ii+1]) || (rawline[ii+1]=='+') || (rawline[ii+1]=='-')) // after it a number a plus or a minus
            change = false;
	 }

   if(change)
      rawline[ii++] = ' ';			// then change it to ' ' 
   else
	  ii++;

 }
 rawline[linelen] = '\0';

 // now parsing the line
 foundtok = 0;
 ii = 0;
 while((foundtok<tokennumber)&&(ii<linelen))
 {
  //find the next token begin
  while(rawline[ii] == ' ')
   ii++;

  foundtok++;
  sscanf(&(rawline[ii]),"%f", &(linebuffer[tokenorder[foundtok-1]]));
 
  // going after the token
  while((rawline[ii] != ' ')&&(rawline[ii] != '\0'))
   ii++;
 }

 if(foundtok<tokennumber)
   return E_LINEERROR;
 else
   return E_NOERROR;
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
  int ii;
  int ret;
  FILE *fp;
  int rownumber;
  int colnumber;

  // line description
  int   tokennumber;
  int   tokenorder[RAW_MAX_TOKEN_LINE_DESCRIPTOR];

  // line read from file, to be parsed
  char rawline[512];

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

  // in raw files with indication of rows and columns it's also possible to triangulate points
  // after the skipped lines there should be the number of row and columns
  if(triangulate)
  {
   fscanf(fp,"%i", &(rownumber));
   fscanf(fp,"%i\n", &(colnumber));
  }

  // parsing line description
  ret = Parselinedescription(linedesc, tokennumber, tokenorder);
  if(ret)
   return ret;

  m.Clear();
//  FaceIterator fi=Allocator<OpenMeshType>::AddFaces(m,facenum);
//  VertexIterator vi=Allocator<OpenMeshType>::AddVertices(m,facenum*3);

  while(!feof(fp))
  {
	  /**/
   //read a new line
   ii=0;
   fread(&(rawline[ii++]),sizeof(char),1,fp);
   while(rawline[ii-1] != '\n')
    fread(&(rawline[ii++]),sizeof(char),1,fp);
   rawline[ii-1] = '\0';

   if(strlen(rawline) >0)  // empty line, just skip
   {

    ret = Parseline(tokennumber, tokenorder, rawline, linebuffer);
    if(ret)
     return ret;

    /*
    // old code reading directly fom file stream
    for(ii=0; ii<tokennumber; ii++)
      fscanf(fp,"%f", &(linebuffer[tokenorder[ii]]));
    */

    // new vertex
    VertexType nv;

    // store the position
    nv.P()[0] = linebuffer[0];    nv.P()[1] = linebuffer[1];    nv.P()[2] = linebuffer[2];
    // store the normal
    if(m.HasPerVertexNormal())
    {
     nv.N()[0] = linebuffer[3];     nv.N()[1] = linebuffer[4];     nv.N()[2] = linebuffer[5];
    }

    // store the color
    if(m.HasPerVertexColor())
    {
     nv.C()[0] = linebuffer[6];     nv.C()[1] = linebuffer[7];     nv.C()[2] = linebuffer[8];
    }

    // store the reflectance
    if(m.HasPerVertexQuality())
    {
     nv.Q() = linebuffer[9];
    }

    m.vert.push_back(nv);
   } // end if zero length
  }

  // update model point number
  m.vn = m.vert.size();

  // now generate the triangles
  if(triangulate)
  {
    int rr,cc;
    if(m.vn != (rownumber * colnumber))
      return E_WRONGPOINTNUM;

    int trinum = (rownumber-1) * (colnumber-1) * 2;

    FaceIterator fi=Allocator<MESH_TYPE>::AddFaces(m,trinum);

	m.fn = 0;
	for(rr=0; rr<rownumber-1; rr++)
	 for(cc=0; cc<colnumber-1; cc++)
	 {
		// upper tri
		(*fi).V(0) = &(m.vert[(rr  ) + ((cc  ) * rownumber)]);
		(*fi).V(1) = &(m.vert[(rr+1) + ((cc  ) * rownumber)]);
		(*fi).V(2) = &(m.vert[(rr  ) + ((cc+1) * rownumber)]);

		 m.fn++;
		 fi++;

		// lower tri
		(*fi).V(0) = &(m.vert[(rr+1) + ((cc  ) * rownumber)]);
		(*fi).V(1) = &(m.vert[(rr+1) + ((cc+1) * rownumber)]);
		(*fi).V(2) = &(m.vert[(rr  ) + ((cc+1) * rownumber)]);
/**/
		 m.fn++;
		 fi++;
	 
	 }
	 
		 
      ++fi;
  }

  fclose(fp);
  return E_NOERROR;
}

}; // end class
} // end Namespace tri
} // end Namespace io
} // end Namespace vcg

#endif
