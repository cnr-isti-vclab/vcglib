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
Revision 1.5  2004/03/18 15:29:07  cignoni
Completed Octahedron and Icosahedron

Revision 1.4  2004/03/11 18:14:19  tarini
prova

Revision 1.2  2004/03/03 16:11:46  cignoni
First working version (tetrahedron!)


****************************************************************************/

#ifndef __VCGLIB_PLATONIC
#define __VCGLIB_PLATONIC

//#include <vcg/Mesh/Refine.h>
#include<vcg/complex/trimesh/allocate.h>
namespace vcg {
namespace tri {

template <class TetraMeshType>
void Tetrahedron(TetraMeshType &in)
{
 typedef TetraMeshType MeshType; 
 typedef typename MeshType::CoordType CoordType;
 typedef typename MeshType::VertexPointer  VertexPointer;
 typedef typename MeshType::VertexIterator VertexIterator;
 typedef typename MeshType::FaceIterator   FaceIterator;

 in.Clear();
 Allocator<TetraMeshType>::AddVertices(in,4);
 Allocator<TetraMeshType>::AddFaces(in,4);

 VertexPointer ivp[4];

 VertexIterator vi=in.vert.begin();
 ivp[0]=&*vi;(*vi).P()=TetraMeshType::CoordType ( 1, 1, 1); ++vi;
 ivp[1]=&*vi;(*vi).P()=TetraMeshType::CoordType (-1, 1,-1); ++vi;
 ivp[2]=&*vi;(*vi).P()=TetraMeshType::CoordType (-1,-1, 1); ++vi;
 ivp[3]=&*vi;(*vi).P()=TetraMeshType::CoordType ( 1,-1,-1); 

 FaceIterator fi=in.face.begin();
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[1]; (*fi).V(2)=ivp[2]; ++fi;
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[2]; (*fi).V(2)=ivp[3]; ++fi;
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[3]; (*fi).V(2)=ivp[1]; ++fi;
 (*fi).V(0)=ivp[3];  (*fi).V(1)=ivp[2]; (*fi).V(2)=ivp[1]; 
}

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
 ivp[0]=&*vi;(*vi).P()=CoordType (-1,-1,-1); ++vi;
 ivp[1]=&*vi;(*vi).P()=CoordType ( 1,-1,-1); ++vi;
 ivp[2]=&*vi;(*vi).P()=CoordType (-1, 1,-1); ++vi;
 ivp[3]=&*vi;(*vi).P()=CoordType ( 1, 1,-1); ++vi;
 ivp[3]=&*vi;(*vi).P()=CoordType (-1,-1, 1); ++vi;
 ivp[3]=&*vi;(*vi).P()=CoordType ( 1,-1, 1); ++vi;
 ivp[4]=&*vi;(*vi).P()=CoordType (-1, 1, 1); ++vi;
 ivp[5]=&*vi;(*vi).P()=CoordType ( 1, 1, 1); 

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
 ivp[5]=&*vi;(*vi).P()=CoordType ( 0,-1, 0); 

 FaceIterator fi=in.face.begin();
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[1]; (*fi).V(2)=ivp[2]; ++fi;
 (*fi).V(0)=ivp[0];  (*fi).V(1)=ivp[2]; (*fi).V(2)=ivp[3]; 
}

//template <class MESH_TYPE>
//void Sphere(MESH_TYPE &in, const int subdiv = 3 )
//{
//	Icosahedron(in);
//	in.ComputeBorderFlag();
//	int lastsize = 0;
//	for(int i=0;i<subdiv;++i)
//	{
//		Refine<MESH_TYPE, MidPoint<MESH_TYPE> >(in,MidPoint<MESH_TYPE>(),0);
//		MESH_TYPE::vertex_iterator vi;
//
//		for(vi = in.vert.begin()+lastsize;vi!=in.vert.end();++vi)
//			vi->P().Normalize();
//
//		lastsize = in.vert.size();
//	}
//}


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
  VertexPointer ivp[VN];

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
      (*fi).V(1)=ivp[b2+i]; 
      (*fi).V(2)=ivp[b2+(i+1)%D]; 
		}

	if(r1==0) for(i=0;i<D;++i,++fi)
		{
      (*fi).V(0)=ivp[0];
      (*fi).V(1)=ivp[b2+i]; 
      (*fi).V(2)=ivp[b2+(i+1)%D]; 
			in.face.push_back(f);
		}
  if(r2==0)	for(i=0;i<D;++i,++fi){
      (*fi).V(0)=ivp[1];
      (*fi).V(1)=ivp[b1+i]; 
      (*fi).V(2)=ivp[b1+(i+1)%D]; 
		}
	
	if(r1!=0 && r2!=0)for(i=0;i<D;++i)
		{
      (*fi).V(0)=ivp[b1+i];
      (*fi).V(1)=ivp[b2+i]; 
      (*fi).V(2)=ivp[b2+(i+1)%D]; 
      ++fi;
      (*fi).V(0)=ivp[b1+1]; 
      (*fi).V(2)=ivp[b2+(i+1)%D]; 
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
 ivp[3]=&*vi;(*vi).P()=CoordType (bb.min[0],bb.min[1],bb.max[2]); ++vi;
 ivp[3]=&*vi;(*vi).P()=CoordType (bb.max[0],bb.min[1],bb.max[2]); ++vi;
 ivp[4]=&*vi;(*vi).P()=CoordType (bb.min[0],bb.max[1],bb.max[2]); ++vi;
 ivp[5]=&*vi;(*vi).P()=CoordType (bb.max[0],bb.max[1],bb.max[2]); 

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


	/// Questa funzione costruisce una mesh a partire da un insieme di coordiante
	/// ed un insieme di terne di indici di vertici

template <class M,class V, class F >
void Build( M & in, const V & v, const F & f)
{
	in.vn = v.size();
	in.fn = f.size();

	in.vert.clear();
	in.face.clear();

	V::const_iterator vi;

	M::VertexType tv;
	tv.Supervisor_Flags()=0;
	
	for(vi=v.begin();vi!=v.end();++vi)
	{
		tv.P() = M::CoordType( 
			(M::ScalarType)(*vi).Ext(0),
			(M::ScalarType)(*vi).Ext(1),
			(M::ScalarType)(*vi).Ext(2)
		);
		in.vert.push_back(tv);
	}

	vector<M::vertex_pointer> index(in.vn);
	M::vertex_iterator j;
	int k;
	for(k=0,j=in.vert.begin();j!=in.vert.end();++j,++k)
		index[k] = &*j;
	
	F::const_iterator fi;

	M::face_type ft;
	ft.Supervisor_Flags()=0;
	
	for(fi=f.begin();fi!=f.end();++fi)
	{
		assert( (*fi)[0]>=0 );
		assert( (*fi)[1]>=0 );
		assert( (*fi)[2]>=0 );
		assert( (*fi)[0]<in.vn );
		assert( (*fi)[1]<in.vn );
		assert( (*fi)[2]<in.vn );
		ft.V(0) = index[ (*fi)[0] ];
		ft.V(1) = index[ (*fi)[1] ];
		ft.V(2) = index[ (*fi)[2] ];
		in.face.push_back(ft);
	}
}

} // End Namespace TriMesh
} // End Namespace vcg
#endif
