/*#***************************************************************************
 * VertexBase.h                                                     o o      *
 *                                                                o     o    *
 * Visual Computing Group                                         _  O  _    *
 * IEI Institute, CNUCE Institute, CNR Pisa                        \/)\/     *
 *                                                                /\/|       *
 * Copyright(C) 1999 by Paolo Cignoni, Claudio Rocchini              |       *
 * All rights reserved.                                              \       *
 *                                                                           *
 * Permission  to use, copy, modify, distribute  and sell this  software and *
 * its documentation for any purpose is hereby granted without fee, provided *
 * that  the above copyright notice appear  in all copies and that both that *
 * copyright   notice  and  this  permission  notice  appear  in  supporting *
 * documentation. the author makes  no representations about the suitability *
 * of this software for any purpose. It is provided  "as is" without express *
 * or implied warranty.                                                      *
 *                                                                           *
 *****************************************************************************/

/****************************************************************************
  History

 2000 Jul 23 First Working release (only the tetrahedron)
          30 Added Octahedron
			Aug	 4 Added Icosahedron
			Nov 22 Added Hexahedron (cube)
 2001 Apr 19 Added HalfOctahedron
          20 Added Square
      May 16 Added Sphere Function	(CR) (PC)
              "     Cone    "         "    "
              "     Box             (CR)
	  Jul 09 Aggiunta Build (CR)
****************************************************************************/


#ifndef __VCGLIB_PLATONIC
#define __VCGLIB_PLATONIC

#include <vcg/Mesh/Refine.h>

template <class MESH_TYPE>
void Tetrahedron(MESH_TYPE &in)
{
	in.vn=4;
	in.fn=4;
	in.vert.clear();
	in.face.clear();
	MESH_TYPE::vertex_type tv;tv.Supervisor_Flags()=0;
	MESH_TYPE::vectorial_type tp;
	tp=MESH_TYPE::vectorial_type ( 1, 1, 1); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type (-1, 1,-1); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type (-1,-1, 1); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type ( 1,-1,-1); tv.P()=tp; in.vert.push_back(tv);

	vector<MESH_TYPE::vertex_pointer> index(in.vn);
	
	MESH_TYPE::face_type f;f.Supervisor_Flags()=0;
	
	MESH_TYPE::vertex_iterator vi;
	int j;
	for(j=0,vi=in.vert.begin();j<in.vn;++j,++vi)	index[j] = &*vi;
	f.V(0)=index[0];  f.V(1)=index[1];f.V(2)=index[2];  in.face.push_back(f);
	f.V(0)=index[0];  f.V(1)=index[2];f.V(2)=index[3];  in.face.push_back(f);
	f.V(0)=index[0];  f.V(1)=index[3];f.V(2)=index[1];  in.face.push_back(f);
	f.V(0)=index[3];  f.V(1)=index[2];f.V(2)=index[1];  in.face.push_back(f);
}

template <class MESH_TYPE>
void Octahedron(MESH_TYPE &in)
{
	in.vn=6;
	in.fn=8;
	in.vert.clear();
	in.face.clear();
	MESH_TYPE::vertex_type tv;tv.Supervisor_Flags()=0;
	MESH_TYPE::vectorial_type tp;
	tp=MESH_TYPE::vectorial_type ( 1, 0, 0); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type ( 0, 1, 0); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type ( 0, 0, 1); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type (-1, 0, 0); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type ( 0,-1, 0); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type ( 0, 0,-1); tv.P()=tp; in.vert.push_back(tv);

	vector<MESH_TYPE::vertex_pointer> index(in.vn);
	
	MESH_TYPE::face_type f;f.Supervisor_Flags()=0;
	
	MESH_TYPE::vertex_iterator vi;
	int j;
	for(j=0,vi=in.vert.begin();j<in.vn;++j,++vi)	index[j] = &*vi;

	f.V(0)=index[0];  f.V(1)=index[1];f.V(2)=index[2];  in.face.push_back(f);
	f.V(0)=index[0];  f.V(1)=index[2];f.V(2)=index[4];  in.face.push_back(f);
	f.V(0)=index[0];  f.V(1)=index[4];f.V(2)=index[5];  in.face.push_back(f);
	f.V(0)=index[0];  f.V(1)=index[5];f.V(2)=index[1];  in.face.push_back(f);
	f.V(0)=index[3];  f.V(1)=index[1];f.V(2)=index[5];  in.face.push_back(f);
	f.V(0)=index[3];  f.V(1)=index[5];f.V(2)=index[4];  in.face.push_back(f);
	f.V(0)=index[3];  f.V(1)=index[4];f.V(2)=index[2];  in.face.push_back(f);
	f.V(0)=index[3];  f.V(1)=index[2];f.V(2)=index[1];  in.face.push_back(f);
}

template <class MESH_TYPE>
void Icosahedron(MESH_TYPE &in)
{
	MESH_TYPE::scalar_type L=(Sqrt(5.0)+1.0)/2.0;
	MESH_TYPE::vectorial_type vv[12]={
	MESH_TYPE::vectorial_type ( 0, L, 1),
	MESH_TYPE::vectorial_type ( 0, L,-1),
	MESH_TYPE::vectorial_type ( 0,-L, 1),
	MESH_TYPE::vectorial_type ( 0,-L,-1),

	MESH_TYPE::vectorial_type ( L, 1, 0),
	MESH_TYPE::vectorial_type ( L,-1, 0),
	MESH_TYPE::vectorial_type (-L, 1, 0),
	MESH_TYPE::vectorial_type (-L,-1, 0),

	MESH_TYPE::vectorial_type ( 1, 0, L),
	MESH_TYPE::vectorial_type (-1, 0, L),
	MESH_TYPE::vectorial_type ( 1, 0,-L),
	MESH_TYPE::vectorial_type (-1, 0,-L)
	};

	int ff[20][3]={
		{1,0,4},{0,1,6},{2,3,5},{3,2,7},
		{4,5,10},{5,4,8},{6,7,9},{7,6,11},
		{8,9,2},{9,8,0},{10,11,1},{11,10,3},
		{0,8,4},{0,6,9},{1,4,10},{1,11,6},
		{2,5,8},{2,9,7},{3,10,5},{3,7,11}
	};

	in.vn=12;
	in.fn=20;
	in.vert.clear();
	in.face.clear();
	MESH_TYPE::vertex_type tv;tv.Supervisor_Flags()=0;
	MESH_TYPE::vectorial_type tp;
	for(int i=0;i<in.vn;i++)
		{
			tv.P()=vv[i];
			in.vert.push_back(tv);
		}
	
	vector<MESH_TYPE::vertex_pointer> index(in.vn);
	
	MESH_TYPE::face_type f;f.Supervisor_Flags()=0;
	
	MESH_TYPE::vertex_iterator vi;
	int j;
	for(j=0,vi=in.vert.begin();j<in.vn;++j,++vi)	index[j] = &*vi;
	for(j=0;j<in.fn;++j)
	{	
		f.V(0)=index[ff[j][0]];
		f.V(1)=index[ff[j][1]];
		f.V(2)=index[ff[j][2]];
		in.face.push_back(f);
	}
}

template <class MESH_TYPE>
void Hexahedron(MESH_TYPE &in)
{
	in.vn=8;
	in.fn=12;
	in.vert.clear();
	in.face.clear();
	MESH_TYPE::vertex_type tv;tv.Supervisor_Flags()=0;
	MESH_TYPE::vectorial_type tp;
	tp=MESH_TYPE::vectorial_type (-1,-1,-1); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type ( 1,-1,-1); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type (-1, 1,-1); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type ( 1, 1,-1); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type (-1,-1, 1); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type ( 1,-1, 1); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type (-1, 1, 1); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type ( 1, 1, 1); tv.P()=tp; in.vert.push_back(tv);
	
	vector<MESH_TYPE::vertex_pointer> index(in.vn);
	
	MESH_TYPE::face_type f;f.Supervisor_Flags()=0;
	
	MESH_TYPE::vertex_iterator vi;
	int j;
	for(j=0,vi=in.vert.begin();j<in.vn;++j,++vi)	index[j] = &*vi;
	f.V(0)=index[0];  f.V(1)=index[1];f.V(2)=index[2];  in.face.push_back(f);
	f.V(0)=index[3];  f.V(1)=index[2];f.V(2)=index[1];  in.face.push_back(f);
	
	f.V(0)=index[0];  f.V(1)=index[2];f.V(2)=index[4];  in.face.push_back(f);
	f.V(0)=index[6];  f.V(1)=index[4];f.V(2)=index[2];  in.face.push_back(f);

	f.V(0)=index[0];  f.V(1)=index[4];f.V(2)=index[1];  in.face.push_back(f);
	f.V(0)=index[5];  f.V(1)=index[1];f.V(2)=index[4];  in.face.push_back(f);

	f.V(0)=index[7];  f.V(1)=index[5];f.V(2)=index[6];  in.face.push_back(f);
	f.V(0)=index[4];  f.V(1)=index[6];f.V(2)=index[5];  in.face.push_back(f);

	f.V(0)=index[7];  f.V(1)=index[6];f.V(2)=index[3];  in.face.push_back(f);
	f.V(0)=index[2];  f.V(1)=index[3];f.V(2)=index[6];  in.face.push_back(f);

	f.V(0)=index[7];  f.V(1)=index[3];f.V(2)=index[5];  in.face.push_back(f);
	f.V(0)=index[1];  f.V(1)=index[5];f.V(2)=index[3];  in.face.push_back(f);
}

template <class MESH_TYPE>
void HalfOctahedron(MESH_TYPE &in)
{
	in.vn=5;
	in.fn=4;
	in.vert.clear();
	in.face.clear();
	MESH_TYPE::vertex_type tv;tv.Supervisor_Flags()=0;
	MESH_TYPE::vectorial_type tp;
	tp=MESH_TYPE::vectorial_type ( 1, 0, 0); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type ( 0, 1, 0); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type ( 0, 0, 1); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type (-1, 0, 0); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type ( 0,-1, 0); tv.P()=tp; in.vert.push_back(tv);

	vector<MESH_TYPE::vertex_pointer> index(in.vn);
	
	MESH_TYPE::face_type f;f.Supervisor_Flags()=0;
	
	MESH_TYPE::vertex_iterator vi;
	int j;
	for(j=0,vi=in.vert.begin();j<in.vn;++j,++vi)	index[j] = &*vi;

	f.V(0)=index[0];  f.V(1)=index[1];f.V(2)=index[2];  in.face.push_back(f);
	f.V(0)=index[0];  f.V(1)=index[2];f.V(2)=index[4];  in.face.push_back(f);
	//f.V(0)=index[0];  f.V(1)=index[4];f.V(2)=index[5];  in.face.push_back(f);
	//f.V(0)=index[0];  f.V(1)=index[5];f.V(2)=index[1];  in.face.push_back(f);
	//f.V(0)=index[3];  f.V(1)=index[1];f.V(2)=index[5];  in.face.push_back(f);
	//f.V(0)=index[3];  f.V(1)=index[5];f.V(2)=index[4];  in.face.push_back(f);
	f.V(0)=index[3];  f.V(1)=index[4];f.V(2)=index[2];  in.face.push_back(f);
	f.V(0)=index[3];  f.V(1)=index[2];f.V(2)=index[1];  in.face.push_back(f);
}

template <class MESH_TYPE>
void Square(MESH_TYPE &in)
{
	in.vn=4;
	in.fn=2;
	in.vert.clear();
	in.face.clear();
	MESH_TYPE::vertex_type tv;tv.Supervisor_Flags()=0;
	MESH_TYPE::vectorial_type tp;
	tp=MESH_TYPE::vectorial_type ( 1, 0, 0); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type ( 0, 1, 0); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type (-1, 0, 0); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type ( 0,-1, 0); tv.P()=tp; in.vert.push_back(tv);

	vector<MESH_TYPE::vertex_pointer> index(in.vn);
	
	MESH_TYPE::face_type f;f.Supervisor_Flags()=0;
	
	MESH_TYPE::vertex_iterator vi;
	int j;
	for(j=0,vi=in.vert.begin();j<in.vn;++j,++vi)	index[j] = &*vi;

	f.V(0)=index[0];  f.V(1)=index[1];f.V(2)=index[2];  in.face.push_back(f);
	f.V(0)=index[0];  f.V(1)=index[2];f.V(2)=index[3];  in.face.push_back(f);
}

template <class MESH_TYPE>
void Sphere(MESH_TYPE &in, const int subdiv = 3 )
{
	Icosahedron(in);
	in.ComputeBorderFlag();
	int lastsize = 0;
	for(int i=0;i<subdiv;++i)
	{
		Refine<MESH_TYPE, MidPoint<MESH_TYPE> >(in,MidPoint<MESH_TYPE>(),0);
		MESH_TYPE::vertex_iterator vi;

		for(vi = in.vert.begin()+lastsize;vi!=in.vert.end();++vi)
			vi->P().Normalize();

		lastsize = in.vert.size();
	}
}


	/// r1 = raggio 1, r2 = raggio2, h = altezza (asse y)
template <class MESH_TYPE>
void Cone( MESH_TYPE & in,
		  const MESH_TYPE::scalar_type r1,
		  const MESH_TYPE::scalar_type r2,
		  const MESH_TYPE::scalar_type h  )
{
	const int D = 24;
	int i,b1,b2;

	if(r1==0 || r2==0)
	{
		in.vn=D+2;
		in.fn=D*2;
	}
	else
	{
		in.vn=D*2+2;
		in.fn=D*4;
	}

	in.vert.clear();
	in.face.clear();

	MESH_TYPE::vertex_type tv;tv.Supervisor_Flags()=0;
	MESH_TYPE::vectorial_type tp;
	
	tp=MESH_TYPE::vectorial_type ( 0,-h/2,0 );
	tv.P()=tp;
	in.vert.push_back(tv);

	tp=MESH_TYPE::vectorial_type ( 0, h/2,0 );
	tv.P()=tp;
	in.vert.push_back(tv);

	b1 = b2 = 2;

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
			tp=MESH_TYPE::vectorial_type ( x,y,z );
			tv.P()=tp;
			in.vert.push_back(tv);
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
			tp=MESH_TYPE::vectorial_type ( x,y,z );
			tv.P()=tp;
			in.vert.push_back(tv);
		}
	}

	vector<MESH_TYPE::vertex_pointer> index(in.vn);
	
	MESH_TYPE::face_type f;
	f.Supervisor_Flags()=0;
	
	MESH_TYPE::vertex_iterator vi;
	int j;
	for(j=0,vi=in.vert.begin();j<in.vn;++j,++vi)	index[j] = &*vi;

	if(r1!=0)
	{
		for(i=0;i<D;++i)
		{
			f.V(0)=index[0];
			f.V(1)=index[b1+i];
			f.V(2)=index[b1+(i+1)%D];
			in.face.push_back(f);
		}
	}

	if(r2!=0)
	{
		for(i=0;i<D;++i)
		{
			f.V(0)=index[1];
			f.V(1)=index[b2+(i+1)%D];
			f.V(2)=index[b2+i];
			in.face.push_back(f);
		}
	}

	if(r1==0)
	{
		for(i=0;i<D;++i)
		{
			f.V(0)=index[0];
			f.V(1)=index[b2+i];
			f.V(2)=index[b2+(i+1)%D];
			in.face.push_back(f);
		}
	}
	else if(r2==0)
	{
		for(i=0;i<D;++i)
		{
			f.V(0)=index[1];
			f.V(2)=index[b1+i];
			f.V(1)=index[b1+(i+1)%D];
			in.face.push_back(f);
		}
	}
	else
	{
		for(i=0;i<D;++i)
		{
			f.V(0)=index[b1+i];
			f.V(1)=index[b2+i];
			f.V(2)=index[b2+(i+1)%D];
			in.face.push_back(f);
			f.V(0)=index[b1+i];
			f.V(1)=index[b2+(i+1)%D];
			f.V(2)=index[b1+(i+1)%D];
			in.face.push_back(f);
		}
	}
}


template <class MESH_TYPE>
void Box(MESH_TYPE &in, const MESH_TYPE::BOX_TYPE & bb )
{
	in.vn=8;
	in.fn=12;
	in.vert.clear();
	in.face.clear();
	MESH_TYPE::vertex_type tv;tv.Supervisor_Flags()=0;
	MESH_TYPE::vectorial_type tp;
	tp=MESH_TYPE::vectorial_type (bb.min[0],bb.min[1],bb.min[2]); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type (bb.max[0],bb.min[1],bb.min[2]); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type (bb.min[0],bb.max[1],bb.min[2]); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type (bb.max[0],bb.max[1],bb.min[2]); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type (bb.min[0],bb.min[1],bb.max[2]); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type (bb.max[0],bb.min[1],bb.max[2]); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type (bb.min[0],bb.max[1],bb.max[2]); tv.P()=tp; in.vert.push_back(tv);
	tp=MESH_TYPE::vectorial_type (bb.max[0],bb.max[1],bb.max[2]); tv.P()=tp; in.vert.push_back(tv);
	
	vector<MESH_TYPE::vertex_pointer> index(in.vn);
	
	MESH_TYPE::face_type f;f.Supervisor_Flags()=0;
	
	MESH_TYPE::vertex_iterator vi;
	int j;
	for(j=0,vi=in.vert.begin();j<in.vn;++j,++vi)	index[j] = &*vi;
	f.V(0)=index[0];  f.V(1)=index[1];f.V(2)=index[2];  in.face.push_back(f);
	f.V(0)=index[3];  f.V(1)=index[2];f.V(2)=index[1];  in.face.push_back(f);
	
	f.V(0)=index[0];  f.V(1)=index[2];f.V(2)=index[4];  in.face.push_back(f);
	f.V(0)=index[6];  f.V(1)=index[4];f.V(2)=index[2];  in.face.push_back(f);

	f.V(0)=index[0];  f.V(1)=index[4];f.V(2)=index[1];  in.face.push_back(f);
	f.V(0)=index[5];  f.V(1)=index[1];f.V(2)=index[4];  in.face.push_back(f);

	f.V(0)=index[7];  f.V(1)=index[5];f.V(2)=index[6];  in.face.push_back(f);
	f.V(0)=index[4];  f.V(1)=index[6];f.V(2)=index[5];  in.face.push_back(f);

	f.V(0)=index[7];  f.V(1)=index[6];f.V(2)=index[3];  in.face.push_back(f);
	f.V(0)=index[2];  f.V(1)=index[3];f.V(2)=index[6];  in.face.push_back(f);

	f.V(0)=index[7];  f.V(1)=index[3];f.V(2)=index[5];  in.face.push_back(f);
	f.V(0)=index[1];  f.V(1)=index[5];f.V(2)=index[3];  in.face.push_back(f);
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

	M::vertex_type tv;
	tv.Supervisor_Flags()=0;
	
	for(vi=v.begin();vi!=v.end();++vi)
	{
		tv.P() = M::vectorial_type( 
			(M::scalar_type)(*vi).Ext(0),
			(M::scalar_type)(*vi).Ext(1),
			(M::scalar_type)(*vi).Ext(2)
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

#endif
