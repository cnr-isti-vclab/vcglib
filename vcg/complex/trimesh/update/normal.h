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
Revision 1.20  2008/04/18 17:52:08  cignoni
added PerVertexFromCurrentFaceNormal
AreaNormalizeFace NormalizeFace
and shortened PerVertexNormalizedPerFaceNormalized

Revision 1.19  2008/02/15 08:08:59  cignoni
added missing include matrix33

Revision 1.18  2007/12/13 17:57:27  cignoni
removed harmless gcc warnings

Revision 1.17  2007/11/23 17:02:47  cignoni
disambiguated pow call (again)

Revision 1.16  2007/11/23 15:42:11  cignoni
disambiguated pow call

Revision 1.15  2007/11/14 11:56:23  ganovelli
added updating of vertex and face normals

Revision 1.14  2007/07/12 23:11:35  cignoni
added the missing PerVertexNormalizedPerFaceNormalized

Revision 1.13  2007/01/10 17:25:14  matteodelle
*** empty log message ***

Revision 1.12  2006/11/07 15:13:56  zifnab1974
Necessary changes for compilation with gcc 3.4.6. Especially the hash function is a problem

Revision 1.11  2005/12/06 18:22:31  pietroni
changed FaceType::ComputeNormal and FaceType::ComputeNormalizedNormal
with face::ComputeNormal and face::ComputeNormalizedNormal

Revision 1.10  2005/12/06 15:30:45  ponchio
added #include triangle3.h for Normal(...)

added a few FaceType:: instead of face::

Revision 1.9  2005/11/22 15:47:34  cignoni
Moved ComputeNormal and ComputeNormalizedNormal out of the face class (no more a member function!)

Revision 1.8  2005/11/21 21:44:43  cignoni
Moved ComputeNormal and ComputeNormalizedNormal out of the face class (no more a member function!)

Revision 1.7  2005/10/13 08:38:00  cignoni
removed the access to the face member function normal and substituted with vcg::normal(*f);

Revision 1.6  2005/06/17 00:46:09  cignoni
Added a PerVertexNormalizedPerFace (vertex are face/area weighted AND normalized)

Revision 1.5  2005/04/01 13:04:55  fiorin
Minor changes

Revision 1.4  2004/09/09 14:35:14  ponchio
Typename changes for linux

Revision 1.3  2004/08/31 15:18:54  pietroni
minor changes to comply gcc compiler (typename's )

Revision 1.2  2004/03/12 15:22:19  cignoni
Written some documentation and added to the trimes doxygen module

Revision 1.1  2004/03/05 10:59:24  cignoni
Changed name from plural to singular (normals->normal)

Revision 1.1  2004/03/04 00:05:50  cignoni
First working version!

Revision 1.1  2004/02/19 13:11:06  cignoni
Initial commit


****************************************************************************/
#ifndef __VCG_TRI_UPDATE_NORMALS
#define __VCG_TRI_UPDATE_NORMALS

#include <vcg/space/triangle3.h>
#include <vcg/math/matrix33.h>
#include <vcg/simplex/face/component.h>

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
	VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
		 if( !(*vi).IsD() && (*vi).IsRW() )
				 (*vi).N() = NormalType((ScalarType)0,(ScalarType)0,(ScalarType)0);

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
  VertexIterator vi;
 for(vi=m.vert.begin();vi!=m.vert.end();++vi)
   if( !(*vi).IsD() && (*vi).IsRW() )
     (*vi).N() = NormalType((ScalarType)0,(ScalarType)0,(ScalarType)0);

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
 
 VertexIterator vi;
 for(vi=m.vert.begin();vi!=m.vert.end();++vi)
   if( !(*vi).IsD() && (*vi).IsRW() )
     (*vi).N() = NormalType((ScalarType)0,(ScalarType)0,(ScalarType)0);

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
 VertexIterator vi;
 for(vi=m.vert.begin();vi!=m.vert.end();++vi)
   if( !(*vi).IsD() && (*vi).IsRW() )
     (*vi).N() = NormalType((ScalarType)0,(ScalarType)0,(ScalarType)0);

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
