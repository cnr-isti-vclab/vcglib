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

#ifndef __VCG_TRI_UPDATE_NORMALS
#define __VCG_TRI_UPDATE_NORMALS

#include <vcg/space/triangle3.h>
#include <vcg/math/matrix33.h>
#include <vcg/complex/trimesh/update/flag.h>


namespace vcg {
namespace tri {

/// \ingroup trimesh

/// \headerfile normal.h vcg/complex/trimesh/update/normal.h

/// \brief Management, updating and computation of per-vertex and per-face normals.
/**
This class is used to compute or update the normals that can be stored in the vertex or face component of a mesh.
*/

template <class ComputeMeshType>
class UpdateNormals
{
public:
typedef ComputeMeshType MeshType; 	
typedef typename MeshType::VertexType     VertexType;
typedef typename MeshType::CoordType     CoordType;
typedef typename VertexType::NormalType     NormalType;
typedef typename VertexType::ScalarType ScalarType;
typedef typename MeshType::VertexPointer  VertexPointer;
typedef typename MeshType::VertexIterator VertexIterator;
typedef typename MeshType::FaceType       FaceType;
typedef typename MeshType::FacePointer    FacePointer;
typedef typename MeshType::FaceIterator   FaceIterator;

/**
 Set to zero all the normals. Usued by all the face averaging algorithms.
 by default it does not clear the normals of unreferenced vertices because they could be still useful
 */
static void PerVertexClear(ComputeMeshType &m, bool ClearAllVertNormal=false)
{
  assert(HasPerVertexNormal(m));
  if(ClearAllVertNormal)
    UpdateFlags<ComputeMeshType>::VertexClearV(m);
  else
  {
    UpdateFlags<ComputeMeshType>::VertexSetV(m);
    for(FaceIterator f=m.face.begin();f!=m.face.end();++f)
     if( !(*f).IsD() )
       for(int i=0;i<3;++i) (*f).V(i)->ClearV();
   }
  VertexIterator vi;
  for(vi=m.vert.begin();vi!=m.vert.end();++vi)
     if( !(*vi).IsD() && (*vi).IsRW() && (!(*vi).IsV()) )
         (*vi).N() = NormalType((ScalarType)0,(ScalarType)0,(ScalarType)0);
}

/// \brief Calculates the face normal (if stored in the current face type)

static void PerFace(ComputeMeshType &m)
{
	if( !m.HasPerFaceNormal()) return;
	FaceIterator f;
	for(f=m.face.begin();f!=m.face.end();++f)
		if( !(*f).IsD() )	face::ComputeNormal(*f);
}

/// \brief Calculates the vertex normal. Exploiting or current face normals.
/**
	The normal of a vertex v is the weigthed average of the normals of the faces incident on v.
*/
static void PerVertexFromCurrentFaceNormal(ComputeMeshType &m)
{
 if( !m.HasPerVertexNormal()) return;
 
 VertexIterator vi;
 for(vi=m.vert.begin();vi!=m.vert.end();++vi)
   if( !(*vi).IsD() && (*vi).IsRW() )
     (*vi).N()=CoordType(0,0,0);

 FaceIterator fi;
 for(fi=m.face.begin();fi!=m.face.end();++fi)
   if( !(*fi).IsD())
   { 
    for(int j=0; j<3; ++j)
			if( !(*fi).V(j)->IsD())  
					(*fi).V(j)->N() += (*fi).cN();
   }
}
/// \brief Calculates the vertex normal. Exploiting or current face normals.
/**
	The normal of a face f is the average of the normals of the vertices of f.
*/
static void PerFaceFromCurrentVertexNormal(ComputeMeshType &m)
{
	for (FaceIterator fi=m.face.begin(); fi!=m.face.end(); ++fi)
   if( !(*fi).IsD())
	 	{
		NormalType n;
		n.SetZero();
		for(int j=0; j<3; ++j)
			n += fi->V(j)->cN();
		n.Normalize();
		fi->N() = n;
	}
}


///  \brief Calculates the vertex normal. Without exploiting or touching face normals.
/**
 The normal of a vertex v computed as a weighted sum f the incident face normals. 
 The weight is simlply the angle of the involved wedge.  Described in:
 
G. Thurmer, C. A. Wuthrich 
"Computing vertex normals from polygonal facets"
Journal of Graphics Tools, 1998
 */
 
static void PerVertexAngleWeighted(ComputeMeshType &m)
{
	assert(HasPerVertexNormal(m));
  PerVertexClear(m);
 FaceIterator f;
 for(f=m.face.begin();f!=m.face.end();++f)
   if( !(*f).IsD() && (*f).IsR() )
   {
    typename FaceType::NormalType t = vcg::NormalizedNormal(*f);
		NormalType e0 = ((*f).V1(0)->cP()-(*f).V0(0)->cP()).Normalize();
		NormalType e1 = ((*f).V1(1)->cP()-(*f).V0(1)->cP()).Normalize();
		NormalType e2 = ((*f).V1(2)->cP()-(*f).V0(2)->cP()).Normalize();
		
		(*f).V(0)->N() += t*AngleN(e0,-e2);
		(*f).V(1)->N() += t*AngleN(-e0,e1);
		(*f).V(2)->N() += t*AngleN(-e1,e2);
   }
}

///  \brief Calculates the vertex normal. Without exploiting or touching face normals.
/**
 The normal of a vertex v is computed according to the formula described by Nelson Max in 
 Max, N., "Weights for Computing Vertex Normals from Facet Normals", Journal of Graphics Tools, 4(2) (1999)
 
 The weight for each wedge is the cross product of the two edge over the product of the square of the two edge lengths. 
 According to the original paper it is perfect only for spherical surface, but it should perform well...
 */
static void PerVertexWeighted(ComputeMeshType &m)
{
 assert(HasPerVertexNormal(m));

 PerVertexClear(m);

 FaceIterator f;
 for(f=m.face.begin();f!=m.face.end();++f)
   if( !(*f).IsD() && (*f).IsR() )
   {
    typename FaceType::NormalType t = vcg::Normal(*f);
		ScalarType e0 = SquaredDistance((*f).V0(0)->cP(),(*f).V1(0)->cP());
		ScalarType e1 = SquaredDistance((*f).V0(1)->cP(),(*f).V1(1)->cP());
		ScalarType e2 = SquaredDistance((*f).V0(2)->cP(),(*f).V1(2)->cP());
		
		(*f).V(0)->N() += t/(e0*e2);
		(*f).V(1)->N() += t/(e0*e1);
		(*f).V(2)->N() += t/(e1*e2);
   }
}

///  \brief Calculates the vertex normal. Without exploiting or touching face normals.
/**
 The normal of a vertex v is the classical area weigthed average of the normals of the faces incident on v.
 */
 
static void PerVertex(ComputeMeshType &m)
{
 assert(HasPerVertexNormal(m));
 
 PerVertexClear(m);

 FaceIterator f;
 for(f=m.face.begin();f!=m.face.end();++f)
   if( !(*f).IsD() && (*f).IsR() )
   {
    //typename FaceType::NormalType t = (*f).Normal();
    typename FaceType::NormalType t = vcg::Normal(*f);
 
    for(int j=0; j<3; ++j)
     if( !(*f).V(j)->IsD() && (*f).V(j)->IsRW() )  
      (*f).V(j)->N() += t;
   }
}


/// \brief Calculates both vertex and face normals.
/**
 The normal of a vertex v is the weigthed average of the normals of the faces incident on v.
*/

static void PerVertexPerFace(ComputeMeshType &m)
{
 if( !m.HasPerVertexNormal() || !m.HasPerFaceNormal()) return;
 
 PerFace(m);
 PerVertexClear(m);

 FaceIterator f;

 for(f=m.face.begin();f!=m.face.end();++f)
   if( !(*f).IsD() && (*f).IsR() )
   {
     for(int j=0; j<3; ++j)
     if( !(*f).V(j)->IsD() && (*f).V(j)->IsRW() )  
      (*f).V(j)->N() += (*f).cN();
   }
}

/// \brief Calculates both vertex and face normals.
/**
 The normal of a vertex v is the weigthed average of the normals of the faces incident on v.
*/

static void PerVertexNormalizedPerFace(ComputeMeshType &m)
{
	PerVertexPerFace(m);
	NormalizeVertex(m);
}

/// \brief Normalize the lenght of the face normals.
static void NormalizeVertex(ComputeMeshType &m)
{
	VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
		if( !(*vi).IsD() && (*vi).IsRW() ) 
			(*vi).N().Normalize();
}

/// \brief Normalize the lenght of the face normals.
static void NormalizeFace(ComputeMeshType &m)
{
	FaceIterator fi;
	for(fi=m.face.begin();fi!=m.face.end();++fi)
      if( !(*fi).IsD() )	(*fi).N().Normalize();
}

static void AreaNormalizeFace(ComputeMeshType &m)
{
	FaceIterator fi;
	for(fi=m.face.begin();fi!=m.face.end();++fi)
      if( !(*fi).IsD() )	
			{
				(*fi).N().Normalize();
				(*fi).N() = (*fi).N() * DoubleArea(*fi);
			}
}

static void PerVertexNormalizedPerFaceNormalized(ComputeMeshType &m)
{
	PerVertexNormalizedPerFace(m);
	NormalizeFace(m);
}

static void PerFaceRW(ComputeMeshType &m, bool normalize=false)
{
	if( !m.HasPerFaceNormal()) return;

	FaceIterator f;
	bool cn = true;

	if(normalize)
	{
		for(f=m.m.face.begin();f!=m.m.face.end();++f)
		if( !(*f).IsD() && (*f).IsRW() )
		{
			for(int j=0; j<3; ++j)
				if( !(*f).V(j)->IsR()) 	cn = false;
      if( cn ) face::ComputeNormalizedNormal(*f);
			cn = true;
		}
	}
	else
	{
		for(f=m.m.face.begin();f!=m.m.face.end();++f)
			if( !(*f).IsD() && (*f).IsRW() )
			{
				for(int j=0; j<3; ++j)
					if( !(*f).V(j)->IsR()) 	cn = false;

				if( cn )
					(*f).ComputeNormal();
				cn = true;
			}
	}
}


static void PerFaceNormalized(ComputeMeshType &m)
{
	if( !m.HasPerFaceNormal()) return;
	FaceIterator f;
		for(f=m.face.begin();f!=m.face.end();++f)
      if( !(*f).IsD() )	face::ComputeNormalizedNormal(*f);
}

static void PerBitQuadFaceNormalized(ComputeMeshType &m)
{
	if( !m.HasPerFaceNormal()) return;
	PerFace(m);

	FaceIterator f;
	for(f=m.face.begin();f!=m.face.end();++f) {
      if( !(*f).IsD() )	{
        for (int k=0; k<3; k++) if (f->IsF(k)) 
        if (&*f < f->FFp(k)) {
          f->N() = f->FFp(k)->N() = (f->FFp(k)->N() + f->N()).Normalize();
        }
      }
  }
}


/// \brief Calculates the vertex normal.
static void PerVertexNormalized(ComputeMeshType &m)
{
  if( !m.HasPerVertexNormal()) return;
  PerVertex(m);
  for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
   if( !(*vi).IsD() && (*vi).IsRW() )
     (*vi).N().Normalize();
}

/// \brief Multiply the vertex normals by the matrix passed. By default, the scale component is removed.
static void PerVertexMatrix(ComputeMeshType &m, const Matrix44<ScalarType> &mat, bool remove_scaling= true){
	float scale;

	Matrix33<ScalarType> mat33(mat,3);
	
	if( !m.HasPerVertexNormal()) return;

	if(remove_scaling){
		scale = pow(mat33.Determinant(),(ScalarType)(1.0/3.0));
		mat33[0][0]/=scale;
		mat33[1][1]/=scale;
		mat33[2][2]/=scale;
	}
	
  for(VertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
   if( !(*vi).IsD() && (*vi).IsRW() )
     (*vi).N()  = mat33*(*vi).N();
}

/// \brief Multiply the face normals by the matrix passed. By default, the scale component is removed.
static void PerFaceMatrix(ComputeMeshType &m, const Matrix44<ScalarType> &mat, bool remove_scaling= true){
	float scale; 

	Matrix33<ScalarType> mat33(mat,3);

	if( !m.HasPerFaceNormal()) return;

	if(remove_scaling){
		scale = pow(mat33.Determinant(),ScalarType(1.0/3.0));
		mat33[0][0]/=scale;
		mat33[1][1]/=scale;
		mat33[2][2]/=scale;
	}
	
  for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
   if( !(*fi).IsD() && (*fi).IsRW() )
     (*fi).N() = mat33* (*fi).N();
}

}; // end class

}	// End namespace
}	// End namespace


#endif
