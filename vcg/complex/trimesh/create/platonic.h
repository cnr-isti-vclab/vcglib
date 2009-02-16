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
Revision 1.13  2008/02/15 14:38:32  ganovelli
added Cylinder(..). the filename platonic.h is lesser and lesser significant...

Revision 1.12  2007/10/19 22:29:36  cignoni
Re-Wrote basic build function

Revision 1.11  2007/02/01 06:38:27  cignoni
Added small comment to grid function

Revision 1.10  2007/01/27 13:14:34  marfr960
Removed unuseful CoordType test

Revision 1.9  2006/08/23 15:29:44  marfr960
*** empty log message ***

Revision 1.8  2006/03/27 04:18:35  cignoni
Double->Scalar in dodecahedron

Revision 1.7  2006/01/30 08:09:05  cignoni
Corrected Grid

Revision 1.6  2006/01/22 17:10:15  cignoni
Added Grid function (to build range map meshes...)

Revision 1.5  2005/07/11 13:16:34  cignoni
small gcc-related compiling issues (typenames,ending cr, initialization order)

Revision 1.4  2005/07/01 11:17:06  cignoni
Added option of passing a base mesh to Sphere for spherifying it

Revision 1.3  2005/06/17 00:49:29  cignoni
Added missing Sphere function

Revision 1.2  2005/02/25 11:41:08  pietroni
Fixed bug in Square

Revision 1.1  2005/01/19 15:43:15  fiorin
Moved from vcg/complex/trimesh to vcg/complex/trimesh/create

Revision 1.10  2004/10/28 00:54:34  cignoni
Better Doxygen documentation

Revision 1.9  2004/09/24 10:14:38  fiorin
Corrected bug in cone

Revision 1.8  2004/09/22 15:12:38  fiorin
Corrected bug in hexahedron

Revision 1.7  2004/07/09 15:34:29  tarini
Dodecahedron added! (and doxigened a little bit)

Revision 1.6  2004/05/13 21:08:00  cignoni
Conformed C++ syntax to GCC requirements

Revision 1.5  2004/03/18 15:29:07  cignoni
Completed Octahedron and Icosahedron

Revision 1.2  2004/03/03 16:11:46  cignoni
First working version (tetrahedron!)


****************************************************************************/

#ifndef __VCGLIB_PLATONIC
#define __VCGLIB_PLATONIC

#include<vcg/math/base.h>
#include<vcg/complex/trimesh/allocate.h>
#include<vcg/complex/trimesh/refine.h>
#include<vcg/complex/trimesh/update/flag.h>

namespace vcg {
namespace tri {
/** \addtogroup trimesh */
//@{
    /**
        A set of functions that builds meshes 
        that represent surfaces of platonic solids,
				and other simple shapes.

				 The 1st parameter is the mesh that will
				be filled with the solid.
		*/
template <class TetraMeshType>
void Tetrahedron(TetraMeshType &in)
{
 typedef TetraMeshType MeshType; 
 typedef typename TetraMeshType::CoordType CoordType;
 typedef typename TetraMeshType::VertexPointer  VertexPointer;
 typedef typename TetraMeshType::VertexIterator VertexIterator;
 typedef typename TetraMeshType::FaceIterator   FaceIterator;

 in.Clear();
 Allocator<TetraMeshType>::AddVertices(in,4);
 Allocator<TetraMeshType>::AddFaces(in,4);

 VertexPointer ivp[4];
 VertexIterator vi=in.vert.begin();
 ivp[0]=&*vi;(*vi).P()=CoordType ( 1.0, 1.0, 1.0); ++vi;
 ivp[1]=&*vi;(*vi).P()=CoordType (-1.0, 1.0,-1.0); ++vi;
 ivp[2]=&*vi;(*vi).P()=CoordType (-1.0,-1.0, 1.0); ++vi;
 ivp[3]=&*vi;(*vi).P()=CoordType ( 1.0,-1.0,-1.0); 

 FaceIterator fi=in.face.begin();
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[1]; (*fi).V(2)=ivp[2]; ++fi;
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[2]; (*fi).V(2)=ivp[3]; ++fi;
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[3]; (*fi).V(2)=ivp[1]; ++fi;
 (*fi).V(0)=ivp[3];  (*fi).V(1)=ivp[2]; (*fi).V(2)=ivp[1]; 
}


/// builds a Dodecahedron, 
/// (each pentagon is composed of 5 triangles)
template <class DodMeshType>
void Dodecahedron(DodMeshType & in)
{
 typedef DodMeshType MeshType; 
 typedef typename MeshType::CoordType CoordType;
 typedef typename MeshType::VertexPointer  VertexPointer;
 typedef typename MeshType::VertexIterator VertexIterator;
 typedef typename MeshType::FaceIterator   FaceIterator;
 typedef typename MeshType::ScalarType     ScalarType;
 const int N_penta=12;
 const int N_points=62;

 int penta[N_penta*3*3]=
	{20,11, 18,  18, 11,  8,  8, 11,  4,   
		13,23,  4,  4, 23,  8,  8, 23, 16, 
    13, 4, 30, 30,  4, 28, 28, 4,  11, 
    16,34,  8,  8, 34, 18, 18, 34, 36, 
    11,20, 28, 28, 20, 45, 45, 20, 38, 
    13,30, 23, 23, 30, 41, 41, 30, 47,
    16,23, 34, 34, 23, 50, 50, 23, 41, 
    20,18, 38, 38, 18, 52, 52, 18, 36,  
    30,28, 47, 47, 28, 56, 56, 28, 45,  
    50,60, 34, 34, 60, 36, 36, 60, 52, 
    45,38, 56, 56, 38, 60, 60, 38, 52, 
    50,41, 60, 60, 41, 56, 56, 41, 47 };
   //A B   E                D       C
  const ScalarType p=(1.0 + math::Sqrt(5.0)) / 2.0;
  const ScalarType p2=p*p;
  const ScalarType p3=p*p*p;
	ScalarType vv[N_points*3]=
	{
   0, 0, 2*p2,     p2, 0, p3,      p, p2, p3, 
   0, p, p3,       -p, p2, p3,     -p2, 0, p3, 
   -p, -p2, p3,    0,   -p, p3,    p,  -p2, p3,
   p3,  p, p2,     p2,  p2, p2,    0,   p3, p2, 
   -p2, p2, p2,    -p3, p, p2,     -p3, -p, p2, 
   -p2, -p2, p2,   0, -p3, p2,     p2, -p2, p2, 
   p3, -p, p2,     p3, 0, p,       p2, p3, p, 
   -p2, p3, p,     -p3, 0, p,      -p2, -p3, p, 
   p2, -p3, p,     2*p2, 0, 0,     p3, p2, 0, 
   p, p3, 0,       0, 2*p2, 0,     -p, p3, 0, 
   -p3, p2, 0,     -2*p2, 0, 0,    -p3, -p2, 0, 
   -p, -p3, 0,     0, -2*p2, 0,    p, -p3, 0, 
   p3, -p2, 0,     p3, 0, -p,      p2, p3, -p, 
   -p2, p3, -p,    -p3, 0, -p,     -p2, -p3, -p, 
   p2, -p3, -p,    p3, p, -p2,     p2, p2, -p2, 
   0, p3, -p2,     -p2, p2, -p2,   -p3, p, -p2, 
   -p3, -p, -p2,   -p2, -p2, -p2,  0, -p3, -p2, 
   p2, -p2, -p2,   p3, -p, -p2,    p2, 0, -p3, 
   p, p2, -p3,     0, p, -p3,      -p, p2, -p3, 
   -p2, 0, -p3,    -p, -p2, -p3,   0, -p, -p3, 
   p, -p2, -p3,    0, 0, -2*p2
	};
	in.Clear();
	//in.face.clear();
  Allocator<DodMeshType>::AddVertices(in,20+12); 
  Allocator<DodMeshType>::AddFaces(in, 5*12); // five pentagons, each made by 5 tri

	int h,i,j,k=0,m=0;

	bool used[N_points];
	for (i=0; i<N_points; i++) used[i]=false;

	int reindex[20+12 *10];
	ScalarType xx,yy,zz, sx,sy,sz;

	int order[5]={0,1,8,6,2};
	int added[12];

	VertexIterator vi=in.vert.begin();

	for (i=0; i<12; i++) {
		sx=sy=sz=0;
		for (int j=0; j<5; j++) {
			h= penta[ i*9 + order[j]  ]-1;
		  xx=vv[h*3];yy=vv[h*3+1];zz=vv[h*3+2]; sx+=xx; sy+=yy; sz+=zz;
			if (!used[h]) {
				(*vi).P()=CoordType( xx, yy, zz ); vi++;
				used[h]=true;
				reindex[ h ] = m++;
			}
		};
		(*vi).P()=CoordType( sx/5.0, sy/5.0, sz/5.0 ); 	vi++;
		added[ i ] = m++;
	}

  std::vector<VertexPointer> index(in.vn);
	
	for(j=0,vi=in.vert.begin();j<in.vn;++j,++vi)	index[j] = &(*vi);

  FaceIterator fi=in.face.begin();

	for (i=0; i<12; i++) {
		for (j=0; j<5; j++){
	    (*fi).V(0)=index[reindex[penta[i*9 + order[j      ] ] -1 ] ]; 
	    (*fi).V(1)=index[reindex[penta[i*9 + order[(j+1)%5] ] -1 ] ];  
		  (*fi).V(2)=index[added[i] ];
		  fi++;
		}
	};
};

template <class OctMeshType>
void Octahedron(OctMeshType &in)
{
 typedef OctMeshType MeshType; 
 typedef typename MeshType::CoordType CoordType;
 typedef typename MeshType::VertexPointer  VertexPointer;
 typedef typename MeshType::VertexIterator VertexIterator;
 typedef typename MeshType::FaceIterator   FaceIterator;

 in.Clear();
 Allocator<OctMeshType>::AddVertices(in,6);
 Allocator<OctMeshType>::AddFaces(in,8);

 VertexPointer ivp[6];

 VertexIterator vi=in.vert.begin();
 ivp[0]=&*vi;(*vi).P()=CoordType ( 1, 0, 0); ++vi;
 ivp[1]=&*vi;(*vi).P()=CoordType ( 0, 1, 0); ++vi;
 ivp[2]=&*vi;(*vi).P()=CoordType ( 0, 0, 1); ++vi;
 ivp[3]=&*vi;(*vi).P()=CoordType (-1, 0, 0); ++vi;
 ivp[4]=&*vi;(*vi).P()=CoordType ( 0,-1, 0); ++vi;
 ivp[5]=&*vi;(*vi).P()=CoordType ( 0, 0,-1); 

 FaceIterator fi=in.face.begin();
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[1]; (*fi).V(2)=ivp[2]; ++fi;
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[2]; (*fi).V(2)=ivp[4]; ++fi;
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[4]; (*fi).V(2)=ivp[5]; ++fi;
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[5]; (*fi).V(2)=ivp[1]; ++fi;
 (*fi).V(0)=ivp[3];  (*fi).V(1)=ivp[1]; (*fi).V(2)=ivp[5]; ++fi;
 (*fi).V(0)=ivp[3];  (*fi).V(1)=ivp[5]; (*fi).V(2)=ivp[4]; ++fi;
 (*fi).V(0)=ivp[3];  (*fi).V(1)=ivp[4]; (*fi).V(2)=ivp[2]; ++fi;
 (*fi).V(0)=ivp[3];  (*fi).V(1)=ivp[2]; (*fi).V(2)=ivp[1]; 
}

template <class IcoMeshType>
void Icosahedron(IcoMeshType &in)
{
 typedef IcoMeshType MeshType; 
 typedef typename MeshType::ScalarType ScalarType;
 typedef typename MeshType::CoordType CoordType;
 typedef typename MeshType::VertexPointer  VertexPointer;
 typedef typename MeshType::VertexIterator VertexIterator;
 typedef typename MeshType::FaceIterator   FaceIterator;

  ScalarType L=ScalarType((math::Sqrt(5.0)+1.0)/2.0);
	CoordType vv[12]={
	CoordType ( 0, L, 1),
	CoordType ( 0, L,-1),
	CoordType ( 0,-L, 1),
	CoordType ( 0,-L,-1),

	CoordType ( L, 1, 0),
	CoordType ( L,-1, 0),
	CoordType (-L, 1, 0),
	CoordType (-L,-1, 0),

	CoordType ( 1, 0, L),
	CoordType (-1, 0, L),
	CoordType ( 1, 0,-L),
	CoordType (-1, 0,-L)
	};

	int ff[20][3]={
		{1,0,4},{0,1,6},{2,3,5},{3,2,7},
		{4,5,10},{5,4,8},{6,7,9},{7,6,11},
		{8,9,2},{9,8,0},{10,11,1},{11,10,3},
		{0,8,4},{0,6,9},{1,4,10},{1,11,6},
		{2,5,8},{2,9,7},{3,10,5},{3,7,11}
	};


  in.Clear();
  Allocator<IcoMeshType>::AddVertices(in,12);
  Allocator<IcoMeshType>::AddFaces(in,20);
  VertexPointer ivp[12];

  VertexIterator vi;
  int i;
  for(i=0,vi=in.vert.begin();vi!=in.vert.end();++i,++vi){
    (*vi).P()=vv[i];
	  ivp[i]=&*vi;
  }

 FaceIterator fi;
 for(i=0,fi=in.face.begin();fi!=in.face.end();++i,++fi){	
		(*fi).V(0)=ivp[ff[i][0]];
		(*fi).V(1)=ivp[ff[i][1]];
		(*fi).V(2)=ivp[ff[i][2]];
	}
}

template <class MeshType>
void Hexahedron(MeshType &in)
{
 typedef typename MeshType::ScalarType ScalarType;
 typedef typename MeshType::CoordType CoordType;
 typedef typename MeshType::VertexPointer  VertexPointer;
 typedef typename MeshType::VertexIterator VertexIterator;
 typedef typename MeshType::FaceIterator   FaceIterator;

 in.Clear();
 Allocator<MeshType>::AddVertices(in,8);
 Allocator<MeshType>::AddFaces(in,12);

 VertexPointer ivp[8];

 VertexIterator vi=in.vert.begin();

 ivp[7]=&*vi;(*vi).P()=CoordType (-1,-1,-1); ++vi;
 ivp[6]=&*vi;(*vi).P()=CoordType ( 1,-1,-1); ++vi;
 ivp[5]=&*vi;(*vi).P()=CoordType (-1, 1,-1); ++vi;
 ivp[4]=&*vi;(*vi).P()=CoordType ( 1, 1,-1); ++vi;
 ivp[3]=&*vi;(*vi).P()=CoordType (-1,-1, 1); ++vi;
 ivp[2]=&*vi;(*vi).P()=CoordType ( 1,-1, 1); ++vi;
 ivp[1]=&*vi;(*vi).P()=CoordType (-1, 1, 1); ++vi;
 ivp[0]=&*vi;(*vi).P()=CoordType ( 1, 1, 1); 

 FaceIterator fi=in.face.begin();
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[1]; (*fi).V(2)=ivp[2]; ++fi;
 (*fi).V(0)=ivp[3];  (*fi).V(1)=ivp[2]; (*fi).V(2)=ivp[1]; ++fi;
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[2]; (*fi).V(2)=ivp[4]; ++fi;
 (*fi).V(0)=ivp[6];  (*fi).V(1)=ivp[4]; (*fi).V(2)=ivp[2]; ++fi;
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[4]; (*fi).V(2)=ivp[1]; ++fi;
 (*fi).V(0)=ivp[5];  (*fi).V(1)=ivp[1]; (*fi).V(2)=ivp[4]; ++fi;
 (*fi).V(0)=ivp[7];  (*fi).V(1)=ivp[5]; (*fi).V(2)=ivp[6]; ++fi;
 (*fi).V(0)=ivp[4];  (*fi).V(1)=ivp[6]; (*fi).V(2)=ivp[5]; ++fi;
 (*fi).V(0)=ivp[7];  (*fi).V(1)=ivp[6]; (*fi).V(2)=ivp[3]; ++fi;
 (*fi).V(0)=ivp[2];  (*fi).V(1)=ivp[3]; (*fi).V(2)=ivp[6]; ++fi;
 (*fi).V(0)=ivp[7];  (*fi).V(1)=ivp[3]; (*fi).V(2)=ivp[5]; ++fi;
 (*fi).V(0)=ivp[1];  (*fi).V(1)=ivp[5]; (*fi).V(2)=ivp[3]; 
}

template <class MeshType>
void Square(MeshType &in)
{
 typedef typename MeshType::ScalarType ScalarType;
 typedef typename MeshType::CoordType CoordType;
 typedef typename MeshType::VertexPointer  VertexPointer;
 typedef typename MeshType::VertexIterator VertexIterator;
 typedef typename MeshType::FaceIterator   FaceIterator;

 in.Clear();
 Allocator<MeshType>::AddVertices(in,4);
 Allocator<MeshType>::AddFaces(in,2);

 VertexPointer ivp[4];

 VertexIterator vi=in.vert.begin();
 ivp[0]=&*vi;(*vi).P()=CoordType ( 1, 0, 0); ++vi;
 ivp[1]=&*vi;(*vi).P()=CoordType ( 0, 1, 0); ++vi;
 ivp[2]=&*vi;(*vi).P()=CoordType (-1, 0, 0); ++vi;
 ivp[3]=&*vi;(*vi).P()=CoordType ( 0,-1, 0); 

 FaceIterator fi=in.face.begin();
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[1]; (*fi).V(2)=ivp[2]; ++fi;
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[2]; (*fi).V(2)=ivp[3]; 
}

// this function build a sphere starting from a eventually not empty mesh.
// If the mesh is not empty it is 'spherified' and used as base for the subdivision process.
// otherwise an icosahedron is used.
template <class MeshType>
void Sphere(MeshType &in, const int subdiv = 3 )
{
 typedef typename MeshType::ScalarType ScalarType;
 typedef typename MeshType::CoordType CoordType;
 typedef typename MeshType::VertexPointer  VertexPointer;
 typedef typename MeshType::VertexIterator VertexIterator;
 typedef typename MeshType::FaceIterator   FaceIterator;
	if(in.vn==0 && in.fn==0) Icosahedron(in);

	VertexIterator vi;
	for(vi = in.vert.begin(); vi!=in.vert.end();++vi)
		vi->P().Normalize();

	tri::UpdateFlags<MeshType>::FaceBorderFromNone(in);
	tri::UpdateTopology<MeshType>::FaceFace(in);

	size_t lastsize = 0;
	for(int i = 0 ; i < subdiv; ++i)
	{
		Refine< MeshType, MidPoint<MeshType> >(in, MidPoint<MeshType>(), 0);

		for(vi = in.vert.begin() + lastsize; vi != in.vert.end(); ++vi)
			vi->P().Normalize();
 
		lastsize = in.vert.size();
	}
}


	/// r1 = raggio 1, r2 = raggio2, h = altezza (asse y)
template <class MeshType>
void Cone( MeshType& in,
		  const typename MeshType::ScalarType r1,
		  const typename MeshType::ScalarType r2,
		  const typename MeshType::ScalarType h  )
{
 typedef typename MeshType::ScalarType ScalarType;
 typedef typename MeshType::CoordType CoordType;
 typedef typename MeshType::VertexPointer  VertexPointer;
 typedef typename MeshType::VertexIterator VertexIterator;
 typedef typename MeshType::FaceIterator   FaceIterator;

 const int D = 24;
	int i,b1,b2;
  in.Clear();
  int VN,FN;
	if(r1==0 || r2==0) {
		VN=D+2;
		FN=D*2;
	}	else {
		VN=D*2+2;
		FN=D*4;
	}

  Allocator<MeshType>::AddVertices(in,VN);
  Allocator<MeshType>::AddFaces(in,FN);
	VertexPointer  *ivp = new VertexPointer[VN];

  VertexIterator vi=in.vert.begin();
  ivp[0]=&*vi;(*vi).P()=CoordType ( 0,-h/2,0 ); ++vi;
  ivp[1]=&*vi;(*vi).P()=CoordType ( 0, h/2,0 ); ++vi;
 
	b1 = b2 = 2;
 int cnt=2;
	if(r1!=0)
	{
		for(i=0;i<D;++i)
		{
			double a = i*3.14159265358979323846*2/D;
			double s = sin(a);
			double c = cos(a);
			double x,y,z;
			x = r1*c;
			z = r1*s;
			y = -h/2;
			
      ivp[cnt]=&*vi; (*vi).P()= CoordType( x,y,z ); ++vi;++cnt;
		}
		b2 += D;
	}

	if(r2!=0)
	{
		for(i=0;i<D;++i)
		{
			double a = i*3.14159265358979323846*2/D;
			double s = sin(a);
			double c = cos(a);
			double x,y,z;
			x = r2*c;
			z = r2*s;
			y =  h/2;
			
      ivp[cnt]=&*vi; (*vi).P()= CoordType( x,y,z ); ++vi;++cnt;
  		}
	}
	
  FaceIterator fi=in.face.begin();
 
  if(r1!=0) for(i=0;i<D;++i,++fi)	{
      (*fi).V(0)=ivp[0];
      (*fi).V(1)=ivp[b1+i]; 
      (*fi).V(2)=ivp[b1+(i+1)%D]; 
		}

	if(r2!=0) for(i=0;i<D;++i,++fi) {
      (*fi).V(0)=ivp[1];
      (*fi).V(2)=ivp[b2+i]; 
      (*fi).V(1)=ivp[b2+(i+1)%D]; 
		}

	if(r1==0) for(i=0;i<D;++i,++fi)
		{
      (*fi).V(0)=ivp[0];
      (*fi).V(1)=ivp[b2+i]; 
      (*fi).V(2)=ivp[b2+(i+1)%D]; 
			//in.face.push_back(*fi);
		}
  if(r2==0)	for(i=0;i<D;++i,++fi){
      (*fi).V(0)=ivp[1];
      (*fi).V(2)=ivp[b1+i]; 
      (*fi).V(1)=ivp[b1+(i+1)%D];
		}
	
	if(r1!=0 && r2!=0)for(i=0;i<D;++i)
		{
      (*fi).V(0)=ivp[b1+i];
      (*fi).V(1)=ivp[b2+i]; 
      (*fi).V(2)=ivp[b2+(i+1)%D]; 
      ++fi;
      (*fi).V(0)=ivp[b1+i]; 
      (*fi).V(1)=ivp[b2+(i+1)%D]; 
      (*fi).V(2)=ivp[b1+(i+1)%D]; 
      ++fi;
		}		
}


template <class MeshType >
void Box(MeshType &in, const typename MeshType::BoxType & bb )
{
 typedef typename MeshType::ScalarType ScalarType;
 typedef typename MeshType::CoordType CoordType;
 typedef typename MeshType::VertexPointer  VertexPointer;
 typedef typename MeshType::VertexIterator VertexIterator;
 typedef typename MeshType::FaceIterator   FaceIterator;

 in.Clear();
 Allocator<MeshType>::AddVertices(in,8);
 Allocator<MeshType>::AddFaces(in,12);

 VertexPointer ivp[8];

 VertexIterator vi=in.vert.begin();
 ivp[0]=&*vi;(*vi).P()=CoordType (bb.min[0],bb.min[1],bb.min[2]); ++vi;
 ivp[1]=&*vi;(*vi).P()=CoordType (bb.max[0],bb.min[1],bb.min[2]); ++vi;
 ivp[2]=&*vi;(*vi).P()=CoordType (bb.min[0],bb.max[1],bb.min[2]); ++vi;
 ivp[3]=&*vi;(*vi).P()=CoordType (bb.max[0],bb.max[1],bb.min[2]); ++vi;
 ivp[4]=&*vi;(*vi).P()=CoordType (bb.min[0],bb.min[1],bb.max[2]); ++vi;
 ivp[5]=&*vi;(*vi).P()=CoordType (bb.max[0],bb.min[1],bb.max[2]); ++vi;
 ivp[6]=&*vi;(*vi).P()=CoordType (bb.min[0],bb.max[1],bb.max[2]); ++vi;
 ivp[7]=&*vi;(*vi).P()=CoordType (bb.max[0],bb.max[1],bb.max[2]); 

 FaceIterator fi=in.face.begin();
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[1]; (*fi).V(2)=ivp[2]; ++fi;
 (*fi).V(0)=ivp[3];  (*fi).V(1)=ivp[2]; (*fi).V(2)=ivp[1]; ++fi;
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[2]; (*fi).V(2)=ivp[4]; ++fi;
 (*fi).V(0)=ivp[6];  (*fi).V(1)=ivp[4]; (*fi).V(2)=ivp[2]; ++fi;
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[4]; (*fi).V(2)=ivp[1]; ++fi;
 (*fi).V(0)=ivp[5];  (*fi).V(1)=ivp[1]; (*fi).V(2)=ivp[4]; ++fi;
 (*fi).V(0)=ivp[7];  (*fi).V(1)=ivp[5]; (*fi).V(2)=ivp[6]; ++fi;
 (*fi).V(0)=ivp[4];  (*fi).V(1)=ivp[6]; (*fi).V(2)=ivp[5]; ++fi;
 (*fi).V(0)=ivp[7];  (*fi).V(1)=ivp[6]; (*fi).V(2)=ivp[3]; ++fi;
 (*fi).V(0)=ivp[2];  (*fi).V(1)=ivp[3]; (*fi).V(2)=ivp[6]; ++fi;
 (*fi).V(0)=ivp[7];  (*fi).V(1)=ivp[3]; (*fi).V(2)=ivp[5]; ++fi;
 (*fi).V(0)=ivp[1];  (*fi).V(1)=ivp[5]; (*fi).V(2)=ivp[3]; 
}


// this function build a mesh starting from a vector of generic coords (objects having a triple of float at their beginning)
// and a vector of faces (objects having a triple of ints at theri beginning). 


template <class MeshType,class V, class F >
void Build( MeshType & in, const V & v, const F & f)
{
 typedef typename MeshType::ScalarType ScalarType;
 typedef typename MeshType::CoordType CoordType;
 typedef typename MeshType::VertexPointer  VertexPointer;
 typedef typename MeshType::VertexIterator VertexIterator;
 typedef typename MeshType::FaceIterator   FaceIterator;

 Allocator<MeshType>::AddVertices(in,v.size());
 Allocator<MeshType>::AddFaces(in,f.size());

	typename V::const_iterator vi;

	typename MeshType::VertexType tv;
	
	for(int i=0;i<v.size();++i)
	{
		float *vv=(float *)(&v[i]);
		in.vert[i].P() = CoordType( vv[0],vv[1],vv[2]);
	}

  std::vector<VertexPointer> index(in.vn);
	VertexIterator j;
	int k;
	for(k=0,j=in.vert.begin();j!=in.vert.end();++j,++k)
		index[k] = &*j;
	
	typename F::const_iterator fi;

	typename MeshType::FaceType ft;
	
	for(int i=0;i<f.size();++i)
	{
	int * ff=(int *)(&f[i]);
		assert( ff[0]>=0 );
		assert( ff[1]>=0 );
		assert( ff[2]>=0 );
		assert( ff[0]<in.vn );
		assert( ff[1]<in.vn );
		assert( ff[2]<in.vn );
		in.face[i].V(0) = &in.vert[ ff[0] ];
		in.face[i].V(1) = &in.vert[ ff[0] ];
		in.face[i].V(2) = &in.vert[ ff[0] ];

	}
}

// Build a regular grid mesh as a typical height field mesh
// x y are the position on the grid scaled by wl and hl (at the end x is in the range 0..wl and y is in 0..hl)
// z is taken from the <data> array 
// Once generated the vertex positions it uses the FaceGrid function to generate the faces;

template <class MeshType>
void Grid(MeshType & in, int w, int h, float wl, float hl, float *data)
{
 typedef typename MeshType::CoordType CoordType;
 typedef typename MeshType::VertexPointer  VertexPointer;
 typedef typename MeshType::VertexIterator VertexIterator;
 typedef typename MeshType::FaceIterator   FaceIterator;

  in.Clear();
	Allocator<MeshType>::AddVertices(in,w*h);


  float wld=wl/float(w);
  float hld=hl/float(h);

  for(int i=0;i<h;++i)
    for(int j=0;j<w;++j)
      in.vert[i*w+j].P()=CoordType ( j*wld, i*hld, data[i*w+j]);
 FaceGrid(in,w,h);
}


// Build a regular grid mesh of faces as a typical height field mesh
// Vertexes are assumed to be already be allocated.

template <class MeshType>
void FaceGrid(MeshType & in, int w, int h)
{
	assert(in.vn == (int)in.vert.size()); // require a compact vertex vector
	assert(in.vn >= w*h); // the number of vertices should match the number of expected grid vertices
	
	Allocator<MeshType>::AddFaces(in,(w-1)*(h-1)*2);
	
//   i+0,j+0 -- i+0,j+1
//      |   \     |
//      |    \    |
//      |     \   |
//      |      \  |
//   i+1,j+0 -- i+1,j+1
//
  for(int i=0;i<h-1;++i)
    for(int j=0;j<w-1;++j)
    {
      in.face[2*(i*(w-1)+j)+0].V(0) = &(in.vert[(i+0)*w+j+0]);
      in.face[2*(i*(w-1)+j)+0].V(1) = &(in.vert[(i+1)*w+j+1]);
      in.face[2*(i*(w-1)+j)+0].V(2) = &(in.vert[(i+0)*w+j+1]);
  
      in.face[2*(i*(w-1)+j)+1].V(0) = &(in.vert[(i+0)*w+j+0]);
      in.face[2*(i*(w-1)+j)+1].V(1) = &(in.vert[(i+1)*w+j+0]);
      in.face[2*(i*(w-1)+j)+1].V(2) = &(in.vert[(i+1)*w+j+1]);
    }
}


// Build a regular grid mesh of faces as a typical height field mesh
// Vertexes are assumed to be already be allocated, but not oll the grid vertexes are present.
// For this purpos a grid of indexes is also passed. negative indexes means that there is no vertex.

template <class MeshType>
void FaceGrid(MeshType & in, std::vector<int> &grid, int w, int h)
{
	assert(in.vn == (int)in.vert.size()); // require a compact vertex vector
	assert(in.vn <= w*h); // the number of vertices should match the number of expected grid vertices
	
//	    V0       V1
//   i+0,j+0 -- i+0,j+1
//      |   \     |
//      |    \    |
//      |     \   |
//      |      \  |
//   i+1,j+0 -- i+1,j+1
//     V2        V3


  for(int i=0;i<h-1;++i)
    for(int j=0;j<w-1;++j)
    {
			int V0i= grid[(i+0)*w+j+0];
			int V1i= grid[(i+0)*w+j+1];
			int V2i= grid[(i+1)*w+j+0];
			int V3i= grid[(i+1)*w+j+1];
			
			if(V0i>=0 && V2i>=0 && V3i>=0 )
			{
				typename MeshType::FaceIterator f= Allocator<MeshType>::AddFaces(in,1);
				f->V(0)=&(in.vert[V0i]);
				f->V(1)=&(in.vert[V3i]);
				f->V(2)=&(in.vert[V2i]);
			}
			if(V0i>=0 && V1i>=0 && V3i>=0 )
			{
				typename MeshType::FaceIterator f= Allocator<MeshType>::AddFaces(in,1);
				f->V(0)=&(in.vert[V0i]);
				f->V(1)=&(in.vert[V1i]);
				f->V(2)=&(in.vert[V3i]);
			}
    }
}

template <class MeshType>
void Cylinder(const int &slices, const int &stacks, MeshType & m){
	float rad_step = M_PI / (float)stacks;
	float sli_step = 1.0 / (float)slices;
  float angle = 0,heigth = 0;

	float x,y,h;
	typename MeshType::VertexIterator vi = vcg::tri::Allocator<MeshType>::AddVertices(m,slices*stacks);
	for ( int j = 0; j < slices; ++j)
		for ( int i = 0; i < stacks; ++i){

				x = cos(	2.0 * M_PI / stacks * i);
				y = sin(	2.0 * M_PI / stacks * i);
				h = j / (float)slices;

				(*vi).P() = MeshType::CoordType(x,y,h);
				++vi;
			}

	typename MeshType::FaceIterator fi ;
	for ( int j = 0; j < slices-1; ++j)
		for ( int i = 0; i < stacks; ++i){

			 if(((i+j)%2) == 0){
				fi = vcg::tri::Allocator<MeshType>::AddFaces(m,1);
				(*fi).V(0) = &m.vert[j    *stacks+ i];
				(*fi).V(1) = &m.vert[j    *stacks+(i+1)%stacks];
				(*fi).V(2) = &m.vert[(j+1)*stacks+(i+1)%stacks];

			 	fi = vcg::tri::Allocator<MeshType>::AddFaces(m,1);
			 	(*fi).V(0) = &m.vert[j    *stacks+i];
			 	(*fi).V(1) = &m.vert[(j+1)*stacks+(i+1)%stacks];
			 	(*fi).V(2) = &m.vert[(j+1)*stacks+i];
			 }
			 else{
				fi = vcg::tri::Allocator<MeshType>::AddFaces(m,1);
				(*fi).V(0) = &m.vert[(j+1)    *stacks+i];
				(*fi).V(1) = &m.vert[ j      *stacks+i];
				(*fi).V(2) = &m.vert[ j      *stacks+(i+1)%stacks];

			 	fi = vcg::tri::Allocator<MeshType>::AddFaces(m,1);
				(*fi).V(0) = &m.vert[(j+1)    *stacks+i];
				(*fi).V(1) = &m.vert[ j       *stacks+(i+1)%stacks];
				(*fi).V(2) = &m.vert[ (j+1)   *stacks+(i+1)%stacks];

			 }

}

	}
//@}

} // End Namespace TriMesh
} // End Namespace vcg
#endif
