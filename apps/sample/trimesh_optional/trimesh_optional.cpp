/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2012                                           \/)\/    *
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
#include<vcg/complex/complex.h>

#include<vcg/simplex/vertex/component_ocf.h>
#include<vcg/simplex/face/component_ocf.h>

#include<vcg/complex/algorithms/create/platonic.h>
#include<vcg/complex/algorithms/update/topology.h>
#include<vcg/complex/algorithms/update/flag.h>
#include<vcg/complex/algorithms/update/normal.h>
#include<vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/refine.h>

class CFace;
class CFaceOcf;

class CVertex;
class CVertexOcf;

struct MyUsedTypes:		 public	vcg::UsedTypes<vcg::Use<CVertex>::AsVertexType,vcg::Use<CFace>::AsFaceType>{};
struct MyUsedTypesOcf: public vcg::UsedTypes<vcg::Use<CVertexOcf>::AsVertexType,vcg::Use<CFaceOcf>::AsFaceType>{};

// Optional stuff has two suffixes:
// OCF Optional Component Fast
// OCC Optional Component Compact

class CVertex     : public vcg::Vertex<	MyUsedTypes,  vcg::vertex::Coord3f, vcg::vertex::BitFlags,vcg::vertex::Normal3f >{};
class CVertexOcf  : public vcg::Vertex< MyUsedTypesOcf,vcg::vertex::InfoOcf,vcg::vertex::Coord3f,vcg::vertex::QualityfOcf, vcg::vertex::BitFlags,vcg::vertex::Normal3f,vcg::vertex::RadiusfOcf >{};
class CFace       : public vcg::Face< MyUsedTypes,    vcg::face::FFAdj,    vcg::face::VertexRef, vcg::face::BitFlags, vcg::face::Normal3f > {};
class CFaceOcf    : public vcg::Face< MyUsedTypesOcf, vcg::face::InfoOcf, vcg::face::FFAdjOcf, vcg::face::VertexRef, vcg::face::BitFlags, vcg::face::Normal3fOcf > {};

class CMesh       : public vcg::tri::TriMesh<     std::vector<CVertex   >,           std::vector<CFace   > > {};
class CMeshOcf    : public vcg::tri::TriMesh<     vcg::vertex::vector_ocf<CVertexOcf>, vcg::face::vector_ocf<CFaceOcf> > {};


using namespace vcg;
using namespace std;

int main(int , char **)
{
	CMesh cm;
	CMeshOcf cmof;
 
  tri::Tetrahedron(cm);
  tri::Tetrahedron(cmof);
 
  printf("Generated mesh has %i vertices and %i triangular faces\n",cm.vn,cm.fn);
  
  /// Calculates both vertex and face normals.
  /// The normal of a vertex v is the weigthed average of the normals of the faces incident on v.
  /// normals are not normalized
  cmof.face.EnableFFAdjacency();

  tri::UpdateTopology<CMesh   >::FaceFace(cm);
  tri::UpdateTopology<CMeshOcf>::FaceFace(cmof);

  tri::UpdateFlags<CMesh   >::FaceBorderFromFF(cm);
  tri::UpdateFlags<CMeshOcf>::FaceBorderFromFF(cmof);

  tri::UpdateNormal<CMesh   >::PerVertexNormalized(cm);
  tri::UpdateNormal<CMeshOcf>::PerVertexNormalized(cmof);

  printf("Normal of face 0 is %f %f %f\n\n",cm.face[0].N()[0],cm.face[0].N()[1],cm.face[0].N()[2]);
  int t0=0,t1=0;
  while(t1-t0<200)
  {
    t0=clock();
    tri::Refine(cm,tri::MidPointButterfly<CMesh>(cm),0);
    t1=clock();
    tri::Refine(cmof,tri::MidPointButterfly<CMeshOcf>(cmof),0);
  }

	cmof.vert.EnableRadius();
	cmof.vert.EnableQuality();

	unsigned int hh = 0;
	for(CMeshOcf::VertexIterator vi = cmof.vert.begin(); vi!=cmof.vert.end();++vi,++hh){
		if(hh%3==0)
			vcg::tri::Allocator<CMeshOcf>::DeleteVertex(cmof,*vi);
	}

  for(CMeshOcf::VertexIterator vi = cmof.vert.begin(); vi!=cmof.vert.end();++vi)
      {
        if(!(*vi).IsD())
        {
          float q =vi->Q();
          float r =vi->R();
          assert(q==r);
        }
      }

      tri::Allocator<CMeshOcf>::CompactVertexVector(cmof);
      tri::UpdateBounding<CMeshOcf>::Box(cmof);
      for(CMeshOcf::VertexIterator vi = cmof.vert.begin(); vi!=cmof.vert.end();++vi)
      {
        if(!(*vi).IsD())
        {
          float q =vi->Q();
          float r =vi->R();
//          int ii = vcg::tri::Index(cmof, *vi);
          assert(q==r);
        }
      }


  return 0;
}
