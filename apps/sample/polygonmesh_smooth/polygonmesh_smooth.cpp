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

#include <vcg/complex/complex.h>

#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/polygonal_algorithms.h>

#include <wrap/io_trimesh/export_obj.h>
#include <wrap/io_trimesh/import_ply.h>
#include <wrap/io_trimesh/export_off.h>
#include <wrap/io_trimesh/import_off.h>

using namespace vcg;
using namespace std;

// forward declarations
class TFace;
class TVertex;

struct TUsedTypes: public vcg::UsedTypes< vcg::Use<TVertex>::AsVertexType, vcg::Use<TFace>::AsFaceType >{};

/* Definition of a mesh of triangles
*/
class TVertex : public Vertex< TUsedTypes,
    vertex::BitFlags,
    vertex::Coord3f,
    vertex::Normal3f,
    vertex::Mark >{};

class TFace   : public Face<   TUsedTypes,
    face::VertexRef,	// three pointers to vertices
    face::Normal3f,		// normal
    face::BitFlags,		// flags
    face::Mark,     // incremental mark
    face::FFAdj			// three pointers to adjacent faces
> {};

/* the mesh is a container of vertices and a container of faces */
class TMesh   : public vcg::tri::TriMesh< vector<TVertex>, vector<TFace> > {};


/* Definition of a mesh of polygons that also supports half-edges
*/
class PFace;
class PVertex;

struct PUsedTypes: public vcg::UsedTypes<vcg::Use<PVertex>  ::AsVertexType,
                                          vcg::Use<PFace>	::AsFaceType>{};

class PVertex:public vcg::Vertex<	PUsedTypes,
    vcg::vertex::Coord3f,
    vcg::vertex::Normal3f,
    vcg::vertex::Mark,
    vcg::vertex::Qualityf,
    vcg::vertex::Mark,
    vcg::vertex::BitFlags>{} ;

class PFace:public vcg::Face<
     PUsedTypes
    ,vcg::face::PolyInfo // this is necessary  if you use component in vcg/simplex/face/component_polygon.h
                         // It says "this class is a polygon and the memory for its components (e.g. pointer to its vertices
                         // will be allocated dynamically")
    ,vcg::face::PFVAdj	 // Pointer to the vertices (just like FVAdj )
    ,vcg::face::PFVAdj
    ,vcg::face::PFFAdj	 // Pointer to edge-adjacent face (just like FFAdj )
    ,vcg::face::BitFlags // bit flags
    ,vcg::face::Qualityf // quality
    ,vcg::face::Mark     // incremental mark
    ,vcg::face::Normal3f // normal
> {};

class PMesh: public
    vcg::tri::TriMesh<
    std::vector<PVertex>,	// the vector of vertices
    std::vector<PFace >     // the vector of faces
    >{};

using namespace vcg;

int	main(int, char **)
{
  TMesh tm;
  PMesh pm;
  math::MarsenneTwisterRNG RndGen;
  tri::io::ImporterOFF<PMesh>::Open(pm,"../../meshes/fertility_quad.off");
  tri::UpdateBounding<PMesh>::Box(pm);
  tri::UpdateSelection<PMesh>::VertexClear(pm);
  // randomly displace half of the vertices
  float randScale = pm.bbox.Diag()/200.0f;
  ForEachVertex(pm,[&](PVertex &v){
    if(v.P()[0] < pm.bbox.Center()[0]) 
      v.P() += math::GeneratePointInUnitBallUniform<float,math::MarsenneTwisterRNG>(RndGen)*randScale;       
  });
  
  // Select half of the vertices 
  ForEachVertex(pm,[&](PVertex &v){
    if(v.P()[2] < pm.bbox.Center()[2]) 
      v.SetS();  
  });
  printf("Input quad mesh has %i %i (sel %i)\n",pm.VN(), pm.FN(), tri::UpdateSelection<PMesh>::VertexCount(pm) );
  
  tri::io::ImporterPLY<TMesh>::Open(tm,"../../meshes/fertility_tri.ply");
  tri::UpdateNormal<TMesh>::PerVertexNormalizedPerFaceNormalized(tm);
  printf("Input tri mesh has %i %i\n",tm.VN(), tm.FN());
  tri::io::ExporterOFF<PMesh>::Save(pm,"./smooth.off");
  PolygonalAlgorithm<PMesh>::LaplacianReproject(pm,tm,100,0.6,0.3,true);
  tri::io::ExporterOFF<PMesh>::Save(pm,"./smooth2.off");  
}


