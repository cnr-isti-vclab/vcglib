/*****************************************************************************
 * VCGLib                                                                    *
 *																																					 *
 * Visual Computing Group                                                o>  *
 * IEI Institute, CNUCE Institute, CNR Pisa                             <|   *
 *                                                                      / \  *
 * Copyright(C) 1999 by Paolo Cignoni, Claudio Rocchini                      *
 * All rights reserved.                                                      *
 *																																					 *
 * Permission  to use, copy, modify, distribute  and sell this  software and *
 * its documentation for any purpose is hereby granted without fee, provided *
 * that  the above copyright notice appear  in all copies and that both that *
 * copyright   notice  and  this  permission  notice  appear  in  supporting *
 * documentation. the author makes  no representations about the suitability *
 * of this software for any purpose. It is provided  "as is" without express *
 * or implied warranty.                                                      *
 *					                                                         *
 *****************************************************************************/
/****************************************************************************
  History
	2003 Dic 17 modifiche per conversione alla versione template	(gano)
****************************************************************************/

// -----------------------------------------------------------------------------------------------
#ifndef _CMESH_H
#define _CMESH_H
// -----------------------------------------------------------------------------------------------

#pragma warning(disable:4786 4804 4666) 
// standard libraries.
#include <math.h>

// VCG library.
//#include <vcg/tools/plylib.h>
#include <vcg/simplex/vertex/with/vcvmvnvq.h>
#include <vcg/simplex/face/with/rtfmfn.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/complex/trimesh/base.h>



using namespace vcg;
using namespace std;

// Vertex, Face, Mesh and Grid definitions.
class MyEdge;
class CFace;
class CVertex   : public VertexCMNQ< double,MyEdge,CFace > {};
class CFace     : public FaceRTFMFN< CVertex,MyEdge,CFace > {};
class CMesh     : public tri::TriMesh< vector<CVertex>, vector<CFace> > {};

// -----------------------------------------------------------------------------------------------
#endif
// -----------------------------------------------------------------------------------------------
