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
Revision 1.2  2004/03/09 21:26:47  cignoni
cr lf mismatch

Revision 1.1  2004/03/03 15:00:51  cignoni
Initial commit

****************************************************************************/

#include<wrap/callback.h>
#include<wrap/ply/plylib.h>
#include<wrap/io_trimesh/io_mask.h>
#include<wrap/io_trimesh/io_ply.h>
#include<vcg/complex/trimesh/allocate.h>



namespace vcg {
namespace tri {
namespace io {

template <class TYPE>
int PlyType ()  { return 0;}

template <> int PlyType <float >()  { return ply::T_FLOAT; }
template <> int PlyType <double>()  { return ply::T_DOUBLE; }
template <> int PlyType <int   >()  { return ply::T_INT; } 
template <> int PlyType <short >()  { return ply::T_SHORT; }
template <> int PlyType <unsigned char >()  { return ply::T_UCHAR; }

/** 
This class encapsulate a filter for opening ply meshes.
The ply file format is quite extensible...
*/
template <class OpenMeshType>
class ImporterPLY
{
public:

typedef ::vcg::ply::PropDescriptor PropDescriptor ;
typedef typename OpenMeshType::VertexPointer VertexPointer;
typedef typename OpenMeshType::ScalarType ScalarType;
typedef typename OpenMeshType::VertexType VertexType;
typedef typename OpenMeshType::FaceType FaceType;
typedef typename OpenMeshType::VertexIterator VertexIterator;
typedef typename OpenMeshType::FaceIterator FaceIterator;

//template <class T> int PlyType () {	assert(0);  return 0;}

#define MAX_USER_DATA 256
// Struttura ausiliaria per la lettura del file ply
struct LoadPly_FaceAux
{
	unsigned char size;
	int v[512];
	int flags;
	float q;
	float tcoord[32];
	unsigned char ntcoord;
	int tcoordind;
	float colors[32];
	unsigned char ncolors;
	
	unsigned char r;
	unsigned char g;
	unsigned char b;

	unsigned char data[MAX_USER_DATA];  
};

struct LoadPly_TristripAux
{
	int size;
	int *v;
	unsigned char data[MAX_USER_DATA];  
};

// Struttura ausiliaria per la lettura del file ply
template<class S>
struct LoadPly_VertAux
{
	S p[3];
	int flags;
	float q;
	unsigned char r;
	unsigned char g;
	unsigned char b;
	unsigned char data[MAX_USER_DATA];  
};

// Struttura ausiliaria caricamento camera
struct LoadPly_Camera
{
	float view_px;
	float view_py;
	float view_pz;
	float x_axisx;
	float x_axisy;
	float x_axisz;
	float y_axisx;
	float y_axisy;
	float y_axisz;
	float z_axisx;
	float z_axisy;
	float z_axisz;
	float focal;
	float scalex;
	float scaley;
	float centerx;
	float centery;
	int   viewportx;
	int   viewporty;
	float k1;
	float k2;
	float k3;
	float k4;
};

static const  PropDescriptor &VertDesc(int i)  
{
	const static PropDescriptor pv[9]={
		{"vertex", "x",         ply::T_FLOAT, PlyType<ScalarType>(),offsetof(LoadPly_VertAux<ScalarType>,p[0]),0,0,0,0,0},
		{"vertex", "y",         ply::T_FLOAT, PlyType<ScalarType>(),offsetof(LoadPly_VertAux<ScalarType>,p[1]),0,0,0,0,0},
		{"vertex", "z",         ply::T_FLOAT, PlyType<ScalarType>(),offsetof(LoadPly_VertAux<ScalarType>,p[2]),0,0,0,0,0},
		{"vertex", "flags",     ply::T_INT,   ply::T_INT,           offsetof(LoadPly_VertAux<ScalarType>,flags),0,0,0,0,0},
		{"vertex", "quality",   ply::T_FLOAT, ply::T_FLOAT,         offsetof(LoadPly_VertAux<ScalarType>,q),0,0,0,0,0},
		{"vertex", "red"  ,     ply::T_UCHAR, ply::T_UCHAR,         offsetof(LoadPly_VertAux<ScalarType>,r),0,0,0,0,0},
		{"vertex", "green",     ply::T_UCHAR, ply::T_UCHAR,         offsetof(LoadPly_VertAux<ScalarType>,g),0,0,0,0,0},
		{"vertex", "blue" ,     ply::T_UCHAR, ply::T_UCHAR,         offsetof(LoadPly_VertAux<ScalarType>,b),0,0,0,0,0},
		{"vertex", "confidence",ply::T_FLOAT, ply::T_FLOAT,         offsetof(LoadPly_VertAux<ScalarType>,q),0,0,0,0,0},
	};
	return pv[i];
}


static const  PropDescriptor &FaceDesc(int i)  
{		
	const static 	PropDescriptor qf[10]=
	{
		{"face", "vertex_indices", ply::T_INT,   ply::T_INT,   offsetof(LoadPly_FaceAux,v),		     1,0,ply::T_UCHAR,ply::T_UCHAR,offsetof(LoadPly_FaceAux,size) },
		{"face", "flags",          ply::T_INT,   ply::T_INT,   offsetof(LoadPly_FaceAux,flags),     0,0,0,0,0},
		{"face", "quality",        ply::T_FLOAT, ply::T_FLOAT, offsetof(LoadPly_FaceAux,q),         0,0,0,0,0},
		{"face", "texcoord",       ply::T_FLOAT, ply::T_FLOAT, offsetof(LoadPly_FaceAux,tcoord),    1,0,ply::T_UCHAR,ply::T_UCHAR,offsetof(LoadPly_FaceAux,ntcoord) },
		{"face", "color",          ply::T_FLOAT, ply::T_FLOAT, offsetof(LoadPly_FaceAux,colors),    1,0,ply::T_UCHAR,ply::T_UCHAR,offsetof(LoadPly_FaceAux,ncolors) },
		{"face", "texnumber",      ply::T_INT,   ply::T_INT,   offsetof(LoadPly_FaceAux,tcoordind), 0,0,0,0,0},
		{"face", "red"  ,          ply::T_UCHAR, ply::T_UCHAR, offsetof(LoadPly_FaceAux,r),         0,0,0,0,0},
		{"face", "green",          ply::T_UCHAR, ply::T_UCHAR, offsetof(LoadPly_FaceAux,g),         0,0,0,0,0},
		{"face", "blue" ,          ply::T_UCHAR, ply::T_UCHAR, offsetof(LoadPly_FaceAux,b),         0,0,0,0,0},
		{"face", "vertex_index",   ply::T_INT,   ply::T_INT,   offsetof(LoadPly_FaceAux,v),		     1,0,ply::T_UCHAR,ply::T_CHAR,offsetof(LoadPly_FaceAux,size) },
	};
	return qf[i];
}
static const  PropDescriptor &TristripDesc(int i)  
{		
	const static 	PropDescriptor qf[1]=
	{
		{"tristrips","vertex_indices", ply::T_INT,  ply::T_INT,  offsetof(LoadPly_TristripAux,v),		  1,1,ply::T_INT,ply::T_INT,offsetof(LoadPly_TristripAux,size) },				
	};
	return qf[i];
}


static const  PropDescriptor &CameraDesc(int i)  
{		
	const static PropDescriptor cad[23] =
	{
		{"camera","view_px",ply::T_FLOAT,ply::T_FLOAT,offsetof(LoadPly_Camera,view_px),0,0,0,0,0},
		{"camera","view_py",ply::T_FLOAT,ply::T_FLOAT,offsetof(LoadPly_Camera,view_py),0,0,0,0,0},
		{"camera","view_pz",ply::T_FLOAT,ply::T_FLOAT,offsetof(LoadPly_Camera,view_pz),0,0,0,0,0},
		{"camera","x_axisx",ply::T_FLOAT,ply::T_FLOAT,offsetof(LoadPly_Camera,x_axisx),0,0,0,0,0},
		{"camera","x_axisy",ply::T_FLOAT,ply::T_FLOAT,offsetof(LoadPly_Camera,x_axisy),0,0,0,0,0},
		{"camera","x_axisz",ply::T_FLOAT,ply::T_FLOAT,offsetof(LoadPly_Camera,x_axisz),0,0,0,0,0},
		{"camera","y_axisx",ply::T_FLOAT,ply::T_FLOAT,offsetof(LoadPly_Camera,y_axisx),0,0,0,0,0},
		{"camera","y_axisy",ply::T_FLOAT,ply::T_FLOAT,offsetof(LoadPly_Camera,y_axisy),0,0,0,0,0},
		{"camera","y_axisz",ply::T_FLOAT,ply::T_FLOAT,offsetof(LoadPly_Camera,y_axisz),0,0,0,0,0},
		{"camera","z_axisx",ply::T_FLOAT,ply::T_FLOAT,offsetof(LoadPly_Camera,z_axisx),0,0,0,0,0},
		{"camera","z_axisy",ply::T_FLOAT,ply::T_FLOAT,offsetof(LoadPly_Camera,z_axisy),0,0,0,0,0},
		{"camera","z_axisz",ply::T_FLOAT,ply::T_FLOAT,offsetof(LoadPly_Camera,z_axisz),0,0,0,0,0},
		{"camera","focal"  ,ply::T_FLOAT,ply::T_FLOAT,offsetof(LoadPly_Camera,focal  ),0,0,0,0,0},
		{"camera","scalex" ,ply::T_FLOAT,ply::T_FLOAT,offsetof(LoadPly_Camera,scalex ),0,0,0,0,0},
		{"camera","scaley" ,ply::T_FLOAT,ply::T_FLOAT,offsetof(LoadPly_Camera,scaley ),0,0,0,0,0},
		{"camera","centerx",ply::T_FLOAT,ply::T_FLOAT,offsetof(LoadPly_Camera,centerx),0,0,0,0,0},
		{"camera","centery",ply::T_FLOAT,ply::T_FLOAT,offsetof(LoadPly_Camera,centery),0,0,0,0,0},
		{"camera","viewportx",ply::T_INT,ply::T_INT  ,offsetof(LoadPly_Camera,viewportx),0,0,0,0,0},
		{"camera","viewporty",ply::T_INT,ply::T_INT  ,offsetof(LoadPly_Camera,viewporty),0,0,0,0,0},
		{"camera","k1"     ,ply::T_FLOAT,ply::T_FLOAT,offsetof(LoadPly_Camera,k1 ),0,0,0,0,0},
		{"camera","k2"     ,ply::T_FLOAT,ply::T_FLOAT,offsetof(LoadPly_Camera,k2 ),0,0,0,0,0},
		{"camera","k3"     ,ply::T_FLOAT,ply::T_FLOAT,offsetof(LoadPly_Camera,k3 ),0,0,0,0,0},
		{"camera","k4"     ,ply::T_FLOAT,ply::T_FLOAT,offsetof(LoadPly_Camera,k4 ),0,0,0,0,0}
	};
	return cad[i];
}


/// Standard call for reading a mesh
static int Open( OpenMeshType &m, const char * filename, CallBackPos *cb=0)
{
	PlyInfo pi;
  pi.cb=cb; 
  return Open(m, filename, pi);
}

/// Read a mesh and store in loadmask the loaded field
static int Open( OpenMeshType &m, const char * filename, int & loadmask, CallBackPos *cb =0)
{
  PlyInfo pi;
  pi.mask=loadmask;
	return Open(m, filename,pi);
  loadmask=pi.mask;
}


/// read a mesh with all the possible option specified in the PlyInfo obj.
static int Open( OpenMeshType &m, const char * filename, PlyInfo &pi )
{
  assert(filename!=0);
	vector<VertexPointer> index;
	LoadPly_FaceAux fa;
	LoadPly_TristripAux tsa;
	LoadPly_VertAux<ScalarType> va;
  
	pi.mask = 0;
	bool multit = false; // true if texture has a per face int spec the texture index

	va.flags = 42;

  pi.status = ::vcg::ply::E_NOERROR;

  // init defaults
	VertexType tv;
	tv.UberFlags() = 0;
	if( VertexType::HasQuality() ) tv.Q()=1.0;
	if( VertexType::HasColor() )     tv.C()=Color4b(Color4b::White);
	
	FaceType tf;
	tf.UberFlags() = 0;
	if( FaceType::HasFaceQuality() )  tf.Q()=1.0;
	if( FaceType::HasWedgeColor() )   tf.WC(0)=tf.WC(1)=tf.WC(2)=Color4b(Color4b::White);
	if( FaceType::HasFaceColor() )    tf.C()=Color4b(Color4b::White);			
	// Descrittori delle strutture

  //bool isvflags = false;	// Il file contiene i flags 


	// The main descriptor of the ply file 
  vcg::ply::PlyFile pf;
  
  // Open the file and parse the header
  if( pf.Open(filename,vcg::ply::PlyFile::MODE_READ)==-1 )
	{
		pi.status = pf.GetError();
		return -1;
	}
  pi.header = pf.GetHeader();

	// Descrittori della camera
	{  // Check that all the camera properties are present.
		bool found = true;
		for(int i=0;i<23;++i)
		{
			if( pf.AddToRead(CameraDesc(i))==-1 ) {
				found = false;
				break;
			}
		}
		if(found) pi.mask |= PLYMask::PM_CAMERA;
	}

  // Descrittori dati standard (vertex coord e faces)
  if( pf.AddToRead(VertDesc(0))==-1 ) { pi.status = PlyInfo::E_NO_VERTEX; return -1; }
	if( pf.AddToRead(VertDesc(1))==-1 ) { pi.status = PlyInfo::E_NO_VERTEX; return -1; }
	if( pf.AddToRead(VertDesc(2))==-1 ) { pi.status = PlyInfo::E_NO_VERTEX; return -1; }
	if( pf.AddToRead(FaceDesc(0))==-1 ) // Se fallisce si prova anche la sintassi di rapidform con index al posto di indices
		if( pf.AddToRead(FaceDesc(9))==-1 ) 
			if(pf.AddToRead(TristripDesc(0))==-1) // Se fallisce tutto si prova a vedere se ci sono tristrip alla levoy.
						{ pi.status = PlyInfo::E_NO_FACE;   return -1; }

		// Descrittori facoltativi dei flags
	if( pf.AddToRead(VertDesc(3))!=-1 )
		pi.mask |= PLYMask::PM_VERTFLAGS;

  if( VertexType::HasQuality() )
	{
		if( pf.AddToRead(VertDesc(4))!=-1 ||
		    pf.AddToRead(VertDesc(8))!=-1 )
			pi.mask |= PLYMask::PM_VERTQUALITY;
	}

	if( VertexType::HasColor() )
	{
		if( pf.AddToRead(VertDesc(5))!=-1 )
		{
			pf.AddToRead(VertDesc(6));
			pf.AddToRead(VertDesc(7));
			pi.mask |= PLYMask::PM_VERTCOLOR;
		}
	}

			// se ci sono i flag per vertice ci devono essere anche i flag per faccia
	if( pf.AddToRead(FaceDesc(1))!=-1 )
		pi.mask |= PLYMask::PM_FACEFLAGS;

  if( FaceType::HasFaceQuality())
	{
		if( pf.AddToRead(FaceDesc(2))!=-1 )
			pi.mask |= PLYMask::PM_FACEQUALITY;
	}

	if( FaceType::HasFaceColor() )
	{
		if( pf.AddToRead(FaceDesc(6))!=-1 )
		{
			pf.AddToRead(FaceDesc(7));
			pf.AddToRead(FaceDesc(8));
			pi.mask |= PLYMask::PM_FACECOLOR;
		}
	}


	if( FaceType::HasWedgeTexture() )
	{
		if( pf.AddToRead(FaceDesc(3))!=-1 )
		{
			if(pf.AddToRead(FaceDesc(5))==0) {
				multit=true; // try to read also the multi texture indicies
				pi.mask |= PLYMask::PM_WEDGTEXMULTI;
			}
			pi.mask |= PLYMask::PM_WEDGTEXCOORD;
		}
	}

  if( FaceType::HasWedgeColor() || FaceType::HasFaceColor() || VertexType::HasColor())
	{
		if( pf.AddToRead(FaceDesc(4))!=-1 )
		{
			pi.mask |= PLYMask::PM_WEDGCOLOR;
		}
	}

	// Descrittori definiti dall'utente, 
	vector<PropDescriptor> VPV(pi.vdn); // property descriptor relative al tipo LoadPly_VertexAux
  vector<PropDescriptor> FPV(pi.fdn); // property descriptor relative al tipo LoadPly_FaceAux
	if(pi.vdn>0){
	// Compute the total size needed to load additional per vertex data.
		size_t totsz=0;
		for(int i=0;i<pi.vdn;i++){
			VPV[i] = pi.VertexData[i];
			VPV[i].offset1=offsetof(LoadPly_VertAux<ScalarType>,data)+totsz;
			totsz+=pi.VertexData[i].memtypesize();
			if( pf.AddToRead(VPV[i])==-1 ) { pi.status = pf.GetError(); return -1; }
		}
		if(totsz > MAX_USER_DATA) 
		{ 
			pi.status = vcg::ply::E_BADTYPE; 
			return -1; 
		}
	}
	if(pi.fdn>0){
		size_t totsz=0;
		for(int i=0;i<pi.fdn;i++){
			FPV[i] = pi.FaceData[i];
			FPV[i].offset1=offsetof(LoadPly_FaceAux,data)+totsz;
			totsz+=pi.FaceData[i].memtypesize();
			if( pf.AddToRead(FPV[i])==-1 ) { pi.status = pf.GetError(); return -1; }
		}
		if(totsz > MAX_USER_DATA) 
		{ 
      pi.status = vcg::ply::E_BADTYPE; 
			return -1; 
		}
	}

  /**************************************************************/
  /* Main Reading Loop */
  /**************************************************************/
  m.Clear();
  for(int i=0;i<int(pf.elements.size());i++)
	{
		int n = pf.ElemNumber(i);

		if( !strcmp( pf.ElemName(i),"camera" ) )
		{
			pf.SetCurElement(i);

			LoadPly_Camera ca;

			for(int j=0;j<n;++j)
			{
				if( pf.Read( (void *)&(ca) )==-1 )
				{
					pi.status = PlyInfo::E_SHORTFILE;
					return -1;
				}	
				//camera.valid     = true;
				//camera.view_p[0] = ca.view_px;
				//camera.view_p[1] = ca.view_py;
				//camera.view_p[2] = ca.view_pz;
				//camera.x_axis[0] = ca.x_axisx;
				//camera.x_axis[1] = ca.x_axisy;
				//camera.x_axis[2] = ca.x_axisz;
				//camera.y_axis[0] = ca.y_axisx;
				//camera.y_axis[1] = ca.y_axisy;
				//camera.y_axis[2] = ca.y_axisz;
				//camera.z_axis[0] = ca.z_axisx;
				//camera.z_axis[1] = ca.z_axisy;
				//camera.z_axis[2] = ca.z_axisz;
				//camera.f         = ca.focal;
				//camera.s[0]      = ca.scalex;
				//camera.s[1]      = ca.scaley;
				//camera.c[0]      = ca.centerx;
				//camera.c[1]      = ca.centery;
				//camera.viewport[0] = ca.viewportx;
				//camera.viewport[1] = ca.viewporty;
				//camera.k[0]      = ca.k1;
				//camera.k[1]      = ca.k2;
				//camera.k[2]      = ca.k3;
				//camera.k[3]      = ca.k4;
			}
		}
		else if( !strcmp( pf.ElemName(i),"vertex" ) )
		{
			int j;

			pf.SetCurElement(i);
      VertexIterator vi=Allocator<OpenMeshType>::AddVertices(m,n);
			
      for(j=0;j<n;++j)
			{
				if(pi.cb && (j%1000)==0) pi.cb(j*50/n,"Vertex Loading");
				(*vi).UberFlags()=0;
			  if( pf.Read( (void *)&(va) )==-1 )
				{
					pi.status = PlyInfo::E_SHORTFILE;
					return -1;
				}
				
				(*vi).P()[0] = va.p[0];
				(*vi).P()[1] = va.p[1];
				(*vi).P()[2] = va.p[2];

				if( pi.mask & PLYMask::PM_VERTFLAGS )
					(*vi).UberFlags() = va.flags;

				if( pi.mask & PLYMask::PM_VERTQUALITY )
					(*vi).Q() = va.q;

				if( pi.mask & PLYMask::PM_VERTCOLOR )
				{
					(*vi).C()[0] = va.r;
					(*vi).C()[1] = va.g;
					(*vi).C()[2] = va.b;
					(*vi).C()[3] = 255;
				}

				
				for(int k=0;k<pi.vdn;k++)	
					memcpy((char *)(&*vi) + pi.VertexData[k].offset1,
									(char *)(&va) + VPV[k].offset1, 
									VPV[k].memtypesize());
			++vi;
      }

			index.resize(n);
			for(j=0,vi=m.vert.begin();j<n;++j,++vi)
				index[j] = &*vi;
		}
		else if( !strcmp( pf.ElemName(i),"face") )/************************************************************/
		{
			int j;
			
      FaceIterator fi=Allocator<OpenMeshType>::AddFaces(m,n);
			pf.SetCurElement(i);

			for(j=0;j<n;++j)
			{
				int k;

				if(pi.cb && (j%1000)==0) pi.cb(50+j*50/n,"Face Loading");
				if( pf.Read(&fa)==-1 )
				{
					pi.status = PlyInfo::E_SHORTFILE;
					return -1;
				}
				if(fa.size!=3)
				{
					pi.status = PlyInfo::E_NO_3VERTINFACE;
					return -1;
				}

				for(k=0;k<3;++k)
				{
					if( fa.v[k]<0 || fa.v[k]>=m.vn )
					{
						pi.status = PlyInfo::E_BAD_VERT_INDEX;
						return -1;
					}
					(*fi).V(k) = index[ fa.v[k] ];
				}

				if( pi.mask & PLYMask::PM_FACEFLAGS )
				{
					(*fi).UberFlags() = fa.flags;
				}

				if( pi.mask & PLYMask::PM_FACEQUALITY )
				{
					(*fi).Q() = fa.q;
				}

				if( pi.mask & PLYMask::PM_FACECOLOR )
				{
					(*fi).C()[0] = fa.r;
					(*fi).C()[1] = fa.g;
					(*fi).C()[2] = fa.b;
					(*fi).C()[3] = 255;
				}

				if( pi.mask & PLYMask::PM_WEDGTEXCOORD )
				{
					for(int k=0;k<3;++k)
					{
						(*fi).WT(k).u() = fa.tcoord[k*2+0];
						(*fi).WT(k).v() = fa.tcoord[k*2+1];
						if(multit) (*fi).WT(k).n() = fa.tcoordind;
					}
				}

				if( pi.mask & PLYMask::PM_WEDGCOLOR )
				{
					if(FaceType::HasWedgeColor()){
						for(int k=0;k<3;++k)
						{
							(*fi).WC(k)[0] = (unsigned char)(fa.colors[k*3+0]*255);
							(*fi).WC(k)[1] = (unsigned char)(fa.colors[k*3+1]*255);
							(*fi).WC(k)[2] = (unsigned char)(fa.colors[k*3+2]*255);
						}
					}
					if(FaceType::HasFaceColor()){
						{
							(*fi).C()[0] = (unsigned char)((fa.colors[0*3+0]*255+fa.colors[1*3+0]*255+fa.colors[2*3+0]*255)/3.0f);
							(*fi).C()[1] = (unsigned char)((fa.colors[0*3+1]*255+fa.colors[1*3+1]*255+fa.colors[2*3+1]*255)/3.0f);
							(*fi).C()[2] = (unsigned char)((fa.colors[0*3+2]*255+fa.colors[1*3+2]*255+fa.colors[2*3+2]*255)/3.0f);
						}
					}
				}

				for(k=0;k<pi.fdn;k++)	
					memcpy((char *)(&(*fi)) + pi.FaceData[k].offset1,
								 (char *)(&fa) + FPV[k].offset1, 
									FPV[k].memtypesize());
      ++fi;
      }
		}else if( !strcmp( pf.ElemName(i),"tristrips") )//////////////////// LETTURA TRISTRIP DI STANFORD
		{
			int j;
			pf.SetCurElement(i);
			int numvert_tmp = m.vert.size();
			for(j=0;j<n;++j)
			{
				int k;
				if(pi.cb && (j%1000)==0) pi.cb(50+j*50/n,"Tristrip Face Loading");
				if( pf.Read(&tsa)==-1 )
				{
					pi.status = PlyInfo::E_SHORTFILE;
					return -1;
				}
				int remainder=0;
				//int startface=m.face.size();
				for(k=0;k<tsa.size-2;++k)
				{
					if(pi.cb && (k%1000)==0) pi.cb(50+k*50/tsa.size,"Tristrip Face Loading");				
          if(tsa.v[k]<0 || tsa.v[k]>=numvert_tmp )	{
						pi.status = PlyInfo::E_BAD_VERT_INDEX;
						return -1;
					}
				  if(tsa.v[k+2]==-1)
					{
						k+=2;
						if(k%2) remainder=0;
							 else remainder=1;
						continue;
					}
					tf.V(0) = index[ tsa.v[k+0] ];
					tf.V(1) = index[ tsa.v[k+1] ];
					tf.V(2) = index[ tsa.v[k+2] ];
					if((k+remainder)%2) swap (tf.V(0), tf.V(1) );
					m.face.push_back( tf );
				}
			}
		}
		else
		{
				// Skippaggio elementi non gestiti
			int n = pf.ElemNumber(i);
			pf.SetCurElement(i);

			for(int j=0;j<n;j++)
			{
				if( pf.Read(0)==-1)
				{
					pi.status = PlyInfo::E_SHORTFILE;
					return -1;
				}
			}
		}
	}

 // // Parsing texture names
	//textures.clear();
	//normalmaps.clear();

	//for(int co=0;co<int(pf.comments.size());++co)
	//{
	//	const char * TFILE = "TextureFile";
	//	const char * NFILE = "TextureNormalFile";
	//	const char * c = pf.comments[co];
	//	char buf[256];
	//	int i,j,n;

	//	if( !strncmp(c,TFILE,strlen(TFILE)) )
	//	{
	//		strcpy(buf,c+strlen(TFILE)+1);
	//		n = strlen(buf);
	//		for(i=j=0;i<n;i++)
	//			if( buf[i]!=' ' && buf[i]!='\t' && buf[i]>32 && buf[i]<125 )	buf[j++] = buf[i];
	//		
	//		buf[j] = 0;
	//		char buf2[255];
	//		__interpret_texture_name( buf,filename,buf2 );
	//		textures.push_back( xstring(buf2) );
	//	}
	//	if( !strncmp(c,NFILE,strlen(NFILE)) )
	//	{
	//		strcpy(buf,c+strlen(NFILE)+1);
	//		n = strlen(buf);
	//		for(i=j=0;i<n;i++)
	//			if( buf[i]!=' ' && buf[i]!='\t' && buf[i]>32 && buf[i]<125 )	buf[j++] = buf[i];
	//		
	//		buf[j] = 0;
	//		char buf2[255];
	//		__interpret_texture_name( buf,filename,buf2 );
	//		normalmaps.push_back( xstring(buf2) );
	//	}
	//}

  // vn and fn should be correct but if someone wrongly saved some deleted elements they can be wrong.
	m.vn = 0;
	VertexIterator vi;
	for(vi=m.vert.begin();vi!=m.vert.end();++vi)
		if( ! (*vi).IsD() )
			++m.vn;

	m.fn = 0;
	FaceIterator fi;
	for(fi=m.face.begin();fi!=m.face.end();++fi)
		if( ! (*fi).IsD() )
			++m.fn;

	return 0;
}


	// Caricamento camera da un ply
int LoadCamera(const char * filename)
{
	PlyFile pf;
	if( pf.Open(filename,PlyFile::MODE_READ)==-1 )
	{
		pi.status = pf.GetError();
		return -1;
	}


	bool found = true;
	int i;
	for(i=0;i<23;++i)
	{
		if( pf.AddToRead(CameraDesc(i))==-1 )
		{
			found = false;
			break;
		}
	}

	if(!found)
		return -1;

	for(i=0;i<int(pf.elements.size());i++)
	{
		int n = pf.ElemNumber(i);

		if( !strcmp( pf.ElemName(i),"camera" ) )
		{
			pf.SetCurElement(i);

			LoadPly_Camera ca;

			for(int j=0;j<n;++j)
			{
				if( pf.Read( (void *)&(ca) )==-1 )
				{
					pi.status = PlyInfo::E_SHORTFILE;
					return -1;
				}	
				camera.valid     = true;
				camera.view_p[0] = ca.view_px;
				camera.view_p[1] = ca.view_py;
				camera.view_p[2] = ca.view_pz;
				camera.x_axis[0] = ca.x_axisx;
				camera.x_axis[1] = ca.x_axisy;
				camera.x_axis[2] = ca.x_axisz;
				camera.y_axis[0] = ca.y_axisx;
				camera.y_axis[1] = ca.y_axisy;
				camera.y_axis[2] = ca.y_axisz;
				camera.z_axis[0] = ca.z_axisx;
				camera.z_axis[1] = ca.z_axisy;
				camera.z_axis[2] = ca.z_axisz;
				camera.f         = ca.focal;
				camera.s[0]      = ca.scalex;
				camera.s[1]      = ca.scaley;
				camera.c[0]      = ca.centerx;
				camera.c[1]      = ca.centery;
				camera.viewport[0] = ca.viewportx;
				camera.viewport[1] = ca.viewporty;
				camera.k[0]      = ca.k1;
				camera.k[1]      = ca.k2;
				camera.k[2]      = ca.k3;
				camera.k[3]      = ca.k4;
			}
			break;
		}
	}

	return 0;
}


bool LoadMask(const char * filename, int &mask)
{
	mask=0;
	PlyFile pf;
	if( pf.Open(filename,PlyFile::MODE_READ)==-1 )
	{
		pi.status = pf.GetError();
		return false;
	}

	if( pf.AddToRead(VertDesc(0))!=-1 && 
	    pf.AddToRead(VertDesc(1))!=-1 && 
	    pf.AddToRead(VertDesc(2))!=-1 )   mask |= PLYMask::PM_VERTCOORD;

	if( pf.AddToRead(VertDesc(3))!=-1 )		mask |= PLYMask::PM_VERTFLAGS;
	if( pf.AddToRead(VertDesc(4))!=-1 )		mask |= PLYMask::PM_VERTQUALITY;
	if( pf.AddToRead(VertDesc(8))!=-1 )		mask |= PLYMask::PM_VERTQUALITY;
	if( ( pf.AddToRead(VertDesc(5))!=-1 ) && 
		  ( pf.AddToRead(VertDesc(6))!=-1 ) &&
			( pf.AddToRead(VertDesc(7))!=-1 )  )  mask |= PLYMask::PM_VERTCOLOR;
	
	if( pf.AddToRead(FaceDesc(0))!=-1 ) mask |= PLYMask::PM_FACEINDEX;
	if( pf.AddToRead(FaceDesc(1))!=-1 ) mask |= PLYMask::PM_FACEFLAGS;

	if( pf.AddToRead(FaceDesc(2))!=-1 ) mask |= PLYMask::PM_FACEQUALITY;
	if( pf.AddToRead(FaceDesc(3))!=-1 ) mask |= PLYMask::PM_WEDGTEXCOORD;
	if( pf.AddToRead(FaceDesc(5))!=-1 ) mask |= PLYMask::PM_WEDGTEXMULTI;
	if( pf.AddToRead(FaceDesc(4))!=-1 ) mask |= PLYMask::PM_WEDGCOLOR;
	if( ( pf.AddToRead(FaceDesc(6))!=-1 ) && 
		  ( pf.AddToRead(FaceDesc(7))!=-1 ) &&
			( pf.AddToRead(FaceDesc(8))!=-1 )  )  mask |= PLYMask::PM_FACECOLOR;


 return true;
}


}; // end class



} // end namespace tri
} // end namespace io
} // end namespace vcg
