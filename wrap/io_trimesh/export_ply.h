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
Revision 1.1  2004/03/03 15:00:51  cignoni
Initial commit

****************************************************************************/

/**
@name Load and Save in Ply format
*/
//@{

#ifndef __VCGLIB_EXPORT_PLY
#define __VCGLIB_EXPORT_PLY

#include<wrap/io_trimesh/io_mask.h>
#include<wrap/io_trimesh/io_ply.h>

#include <stdio.h>

namespace vcg {
namespace tri {
namespace io {


template <class SaveMeshType>
class ExporterPLY
{
public:
typedef ::vcg::ply::PropDescriptor PropDescriptor ;
typedef typename SaveMeshType::VertexPointer VertexPointer;
typedef typename SaveMeshType::ScalarType ScalarType;
typedef typename SaveMeshType::VertexType VertexType;
typedef typename SaveMeshType::FaceType FaceType;
typedef typename SaveMeshType::VertexIterator VertexIterator;
typedef typename SaveMeshType::FaceIterator FaceIterator;

static bool Save(SaveMeshType &m, const char * filename, bool binary=true)
{
  PlyInfo pi;
  return SavePly(m,filename,binary,pi);
}

static bool SavePly(SaveMeshType &m,  const char * filename, int savemask )
{
	PlyInfo pi;
  pi.mask=savemask;
  return SavePly(m,filename,true,pi);
}

static bool SavePly(SaveMeshType &m,  const char * filename, bool binary, PlyInfo &pi )	// V1.0
{
	FILE * fpout;
	int i;
  
	const char * hbin = "binary_little_endian";
	const char * hasc = "ascii";
	const char * h;

	bool multit = false;

	if(binary) h=hbin;
	else       h=hasc;

	fpout = fopen(filename,"wb");
  if(fpout==NULL)		{
    pi.status=::vcg::ply::E_CANTOPEN;
    return false;
  }
	fprintf(fpout,
		"ply\n"
		"format %s 1.0\n"
		"comment VCGLIB generated\n"
		,h
	);

	if( pi.mask & PLYMask::PM_WEDGTEXCOORD )
	{
		//const char * TFILE = "TextureFile";

		//for(i=0;i<textures.size();++i)
		//	fprintf(fpout,"comment %s %s\n", TFILE, (const char *)(textures[i]) );

		//if(textures.size()>1 && (HasPerWedgeTexture() || HasPerVertexTexture())) multit = true;
	}

	//if( (pi.mask & PLYMask::PM_CAMERA) && camera.IsValid() )
	/*{
		fprintf(fpout,
			"element camera 1\n"
			"property float view_px\n"
			"property float view_py\n"
			"property float view_pz\n"
			"property float x_axisx\n"
			"property float x_axisy\n"
			"property float x_axisz\n"
			"property float y_axisx\n"
			"property float y_axisy\n"
			"property float y_axisz\n"
			"property float z_axisx\n"
			"property float z_axisy\n"
			"property float z_axisz\n"
			"property float focal\n"
			"property float scalex\n"
			"property float scaley\n"
			"property float centerx\n"
			"property float centery\n"
			"property int viewportx\n"
			"property int viewporty\n"
			"property float k1\n"
			"property float k2\n"
			"property float k3\n"
			"property float k4\n"
		);
	}*/

	fprintf(fpout,
		"element vertex %d\n"
		"property float x\n"
		"property float y\n"
		"property float z\n"
		,m.vn
	);


	if( pi.mask & PLYMask::PM_VERTFLAGS )
	{
		fprintf(fpout,
			"property int flags\n"
		);
	}
	
	if( HasPerVertexColor()  && (pi.mask & PLYMask::PM_VERTCOLOR) )
	{
		fprintf(fpout,
			"property uchar red\n"
			"property uchar green\n"
			"property uchar blue\n"
			"property uchar alpha\n"
		);
	}

	if( HasPerVertexQuality() && (pi.mask & PLYMask::PM_VERTQUALITY) )
	{
		fprintf(fpout,
			"property float quality\n"
		);
	}

	for(i=0;i<pi.vdn;i++)
			fprintf(fpout,"property %s %s\n",pi.VertexData[i].stotypename(),pi.VertexData[i].propname);
	
	fprintf(fpout,
		"element face %d\n"
		"property list uchar int vertex_indices\n"
		,m.fn
	);

	if( pi.mask & PLYMask::PM_FACEFLAGS )
	{
		fprintf(fpout,
			"property int flags\n"
		);
	}

	if( ( (vertex_type::OBJ_TYPE & vertex_type::OBJ_TYPE_T ) && (pi.mask & PLYMask::PM_VERTTEXCOORD) ) ||
	    ( (face_type::  OBJ_TYPE & face_type::  OBJ_TYPE_WT) && (pi.mask & PLYMask::PM_WEDGTEXCOORD) )  )
	{
		fprintf(fpout,
			"property list uchar float texcoord\n"
		);

		if(multit)
			fprintf(fpout,
				"property int texnumber\n"
			);
	}

	if( HasPerFaceColor() && (pi.mask & PLYMask::PM_FACECOLOR) )
	{
		fprintf(fpout,
			"property uchar red\n"
			"property uchar green\n"
			"property uchar blue\n"
			"property uchar alpha\n"
		);
	}
	
	if ( HasPerWedgeColor() && (pi.mask & PLYMask::PM_WEDGCOLOR)  )
	{
		fprintf(fpout,
			"property list uchar float color\n"
		);
	}

	if( HasPerFaceQuality() && (pi.mask & PLYMask::PM_FACEQUALITY) )
	{
		fprintf(fpout,
			"property float quality\n"
		);
	}

	for(i=0;i<pi.fdn;i++)
			fprintf(fpout,"property %s %s\n",pi.FaceData[i].stotypename(),pi.FaceData[i].propname);

	fprintf(fpout, "end_header\n"	);

		// Salvataggio camera
	if( (pi.mask & PLYMask::PM_CAMERA) && camera.IsValid() )
	{
		//if(binary)
		//{
		//	float t[17];

		//	t[ 0] = camera.view_p[0];
		//	t[ 1] = camera.view_p[1];
		//	t[ 2] = camera.view_p[2];
		//	t[ 3] = camera.x_axis[0];
		//	t[ 4] = camera.x_axis[1];
		//	t[ 5] = camera.x_axis[2];
		//	t[ 6] = camera.y_axis[0];
		//	t[ 7] = camera.y_axis[1];
		//	t[ 8] = camera.y_axis[2];
		//	t[ 9] = camera.z_axis[0];
		//	t[10] = camera.z_axis[1];
		//	t[11] = camera.z_axis[2];
		//	t[12] = camera.f;
		//	t[13] = camera.s[0];
		//	t[14] = camera.s[1];
		//	t[15] = camera.c[0];
		//	t[16] = camera.c[1];
		//	fwrite(t,sizeof(float),17,fpout);

		//	fwrite( camera.viewport,sizeof(int),2,fpout );

		//	t[ 0] = camera.k[0];
		//	t[ 1] = camera.k[1];
		//	t[ 2] = camera.k[2];
		//	t[ 3] = camera.k[3];
		//	fwrite(t,sizeof(float),4,fpout);
		//}
		//else
		//{
		//	fprintf(fpout,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %d %d %g %g %g %g\n"
		//		,camera.view_p[0]
		//		,camera.view_p[1]
		//		,camera.view_p[2]
		//		,camera.x_axis[0]
		//		,camera.x_axis[1]
		//		,camera.x_axis[2]
		//		,camera.y_axis[0]
		//		,camera.y_axis[1]
		//		,camera.y_axis[2]
		//		,camera.z_axis[0]
		//		,camera.z_axis[1]
		//		,camera.z_axis[2]
		//		,camera.f
		//		,camera.s[0]
		//		,camera.s[1]
		//		,camera.c[0]
		//		,camera.c[1]
		//		,camera.viewport[0]
		//		,camera.viewport[1]
		//		,camera.k[0]
		//		,camera.k[1]
		//		,camera.k[2]
		//		,camera.k[3]
		//	);
		//}		
	}


	int j;
	vector<int> FlagV; 
	MVTYPE * vp;
	vertex_iterator vi;
	for(j=0,vi=m.vert.begin();vi!=m.vert.end();++vi)
	{
		vp=&(*vi);
		FlagV.push_back(vp->Supervisor_Flags()); // Salva in ogni caso flag del vertice
		if( ! vp->IsDeleted() )
		{
			if(binary)
			{
				float t;

				t = float(vp->Supervisor_P()[0]); fwrite(&t,sizeof(float),1,fpout);
				t = float(vp->Supervisor_P()[1]); fwrite(&t,sizeof(float),1,fpout);
				t = float(vp->Supervisor_P()[2]); fwrite(&t,sizeof(float),1,fpout);
				
				if( pi.mask & PLYMask::PM_VERTFLAGS )
					fwrite(&(vp->Supervisor_Flags()),sizeof(int),1,fpout);

				if( HasPerVertexColor() && (pi.mask & PLYMask::PM_VERTCOLOR) )
					fwrite(&( vp->C() ),sizeof(char),4,fpout);

				if( HasPerVertexQuality() && (pi.mask & PLYMask::PM_VERTQUALITY) )
					fwrite(&( vp->Q() ),sizeof(float),1,fpout);


				for(i=0;i<vdn;i++)
				{
					double td; float tf;int ti;short ts; char tc; unsigned char tuc;
					switch (VertexData[i].stotype1)
					{
					case T_FLOAT	:		PlyConv(VertexData[i].memtype1,  ((char *)vp)+VertexData[i].offset1, tf );	fwrite(&tf, sizeof(float),1,fpout); break;
					case T_DOUBLE :			PlyConv(VertexData[i].memtype1,  ((char *)vp)+VertexData[i].offset1, td );	fwrite(&td, sizeof(double),1,fpout); break;
					case T_INT		:		PlyConv(VertexData[i].memtype1,  ((char *)vp)+VertexData[i].offset1, ti );	fwrite(&ti, sizeof(int),1,fpout); break;
					case T_SHORT	:		PlyConv(VertexData[i].memtype1,  ((char *)vp)+VertexData[i].offset1, ts );	fwrite(&ts, sizeof(short),1,fpout); break;
					case T_CHAR		:		PlyConv(VertexData[i].memtype1,  ((char *)vp)+VertexData[i].offset1, tc );	fwrite(&tc, sizeof(char),1,fpout); break;
					case T_UCHAR	:		PlyConv(VertexData[i].memtype1,  ((char *)vp)+VertexData[i].offset1, tuc);	fwrite(&tuc,sizeof(unsigned char),1,fpout); break;
					default : assert(0);
					}
				}
			}
			else 	// ***** ASCII *****
			{
				fprintf(fpout,"%g %g %g " ,vp->P()[0],vp->P()[1],vp->P()[2]);

				if( pi.mask & PLYMask::PM_VERTFLAGS )
					fprintf(fpout,"%d ",vp->Supervisor_Flags());

				if( (vertex_type::OBJ_TYPE & vertex_type::OBJ_TYPE_C) && (pi.mask & PLYMask::PM_VERTCOLOR) )
					fprintf(fpout,"%d %d %d %d ",vp->C()[0],vp->C()[1],vp->C()[2],vp->C()[3] );

				if( (vertex_type::OBJ_TYPE & vertex_type::OBJ_TYPE_Q) && (pi.mask & PLYMask::PM_VERTQUALITY) )
					fprintf(fpout,"%g ",vp->Q());

				for(i=0;i<vdn;i++)
				{
					float tf;
					int ti;
					switch (VertexData[i].memtype1)
					{
					case T_FLOAT	:		tf=*( (float  *)        (((char *)vp)+VertexData[i].offset1));	fprintf(fpout,"%g ",tf); break;
					case T_DOUBLE :   tf=*( (double *)        (((char *)vp)+VertexData[i].offset1));	fprintf(fpout,"%g ",tf); break;
					case T_INT		:		ti=*( (int    *)        (((char *)vp)+VertexData[i].offset1));	fprintf(fpout,"%i ",ti); break;
					case T_SHORT	:		ti=*( (short  *)        (((char *)vp)+VertexData[i].offset1));  fprintf(fpout,"%i ",ti); break;
					case T_CHAR		:		ti=*( (char   *)        (((char *)vp)+VertexData[i].offset1));	fprintf(fpout,"%i ",ti); break;
					case T_UCHAR	:		ti=*( (unsigned char *) (((char *)vp)+VertexData[i].offset1));	fprintf(fpout,"%i ",ti); break;
					default : assert(0);
					}
				}

				fprintf(fpout,"\n");
			}

			vp->Supervisor_Flags()=j; // Trucco! Nascondi nei flags l'indice del vertice non deletato!
			j++;
		}
	}
	assert(j==vn);

	char c = 3;
	char k = 9;
	unsigned char b9 = 9;
	unsigned char b6 = 6;
	const MFTYPE * fp;
	int vv[3];
	face_iterator fi;
	int fcnt=0;
	for(j=0,fi=face.begin();fi!=face.end();++fi)
		{
			fp=&(*fi);
			if( ! fp->IsDeleted() )
			{ fcnt++;
				if(binary)
				{
					vv[0]=fp->cV(0)->Supervisor_Flags();
					vv[1]=fp->cV(1)->Supervisor_Flags();
					vv[2]=fp->cV(2)->Supervisor_Flags();
					fwrite(&c,1,1,fpout);
					fwrite(vv,sizeof(int),3,fpout);

					if( pi.mask & PLYMask::PM_FACEFLAGS )
						fwrite(&(fp->Flags()),sizeof(int),1,fpout);

					if( (vertex_type::OBJ_TYPE & vertex_type::OBJ_TYPE_T) && (pi.mask & PLYMask::PM_VERTTEXCOORD) )
					{
						fwrite(&b6,sizeof(char),1,fpout);
						float t[6];
						for(int k=0;k<3;++k)
						{
							t[k*2+0] = fp->V(k)->T().u();
							t[k*2+1] = fp->V(k)->T().v();
						}
						fwrite(t,sizeof(float),6,fpout);
					}
					else if( (face_type::  OBJ_TYPE & face_type::  OBJ_TYPE_WT) && (pi.mask & PLYMask::PM_WEDGTEXCOORD)  )
					{
						fwrite(&b6,sizeof(char),1,fpout);
						float t[6];
						for(int k=0;k<3;++k)
						{
							t[k*2+0] = fp->WT(k).u();
							t[k*2+1] = fp->WT(k).v();
						}
						fwrite(t,sizeof(float),6,fpout);
					}

					if(multit)
					{
						int t = fp->WT(0).n();
						fwrite(&t,sizeof(int),1,fpout);
					}
          
					if( HasPerFaceColor() && (pi.mask & PLYMask::PM_FACECOLOR) )
					   fwrite(&( fp->C() ),sizeof(char),4,fpout);


					if( HasPerWedgeColor() && (pi.mask & PLYMask::PM_WEDGCOLOR)  )
					{
						fwrite(&b9,sizeof(char),1,fpout);
						float t[3];
						for(int z=0;z<3;++z)
						{
							t[0] = float(fp->WC(z)[0])/255;
							t[1] = float(fp->WC(z)[1])/255;
							t[2] = float(fp->WC(z)[2])/255;
							fwrite( t,sizeof(float),3,fpout);
						}
					}

					if( (face_type::OBJ_TYPE & face_type::OBJ_TYPE_Q) && (pi.mask & PLYMask::PM_FACEQUALITY) )
						fwrite( &(fp->Q()),sizeof(float),1,fpout);


					for(i=0;i<fdn;i++)
						{
						double td; float tf;int ti;short ts; char tc; unsigned char tuc;
						switch (FaceData[i].stotype1){
								case T_FLOAT	:		PlyConv(FaceData[i].memtype1,  ((char *)fp)+FaceData[i].offset1, tf );	fwrite(&tf, sizeof(float),1,fpout); break;
								case T_DOUBLE :		PlyConv(FaceData[i].memtype1,  ((char *)fp)+FaceData[i].offset1, td );	fwrite(&td, sizeof(double),1,fpout); break;
								case T_INT		:		PlyConv(FaceData[i].memtype1,  ((char *)fp)+FaceData[i].offset1, ti );	fwrite(&ti, sizeof(int),1,fpout); break;
								case T_SHORT	:		PlyConv(FaceData[i].memtype1,  ((char *)fp)+FaceData[i].offset1, ts );	fwrite(&ts, sizeof(short),1,fpout); break;
								case T_CHAR		:		PlyConv(FaceData[i].memtype1,  ((char *)fp)+FaceData[i].offset1, tc );	fwrite(&tc, sizeof(char),1,fpout); break;
								case T_UCHAR	:		PlyConv(FaceData[i].memtype1,  ((char *)fp)+FaceData[i].offset1, tuc);	fwrite(&tuc,sizeof(unsigned char),1,fpout); break;
								default : assert(0);
						}
					}
				}
				else	// ***** ASCII *****
				{
					fprintf(fpout,"3 %d %d %d ",
						fp->cV(0)->Supervisor_Flags(),	fp->cV(1)->Supervisor_Flags(), fp->cV(2)->Supervisor_Flags() );

					if( pi.mask & PLYMask::PM_FACEFLAGS )
						fprintf(fpout,"%d ",fp->Flags());

					if( (vertex_type::OBJ_TYPE & vertex_type::OBJ_TYPE_T) && (pi.mask & PLYMask::PM_VERTTEXCOORD) )
					{
						fprintf(fpout,"6 ");
						for(int k=0;k<3;++k)
							fprintf(fpout,"%g %g "
								,fp->V(k)->T().u()
								,fp->V(k)->T().v()
							);
					}
					else if( (face_type::  OBJ_TYPE & face_type::  OBJ_TYPE_WT) && (pi.mask & PLYMask::PM_WEDGTEXCOORD)  )
					{
						fprintf(fpout,"6 ");
						for(int k=0;k<3;++k)
							fprintf(fpout,"%g %g "
								,fp->WT(k).u()
								,fp->WT(k).v()
							);
					}

					if(multit)
					{
						fprintf(fpout,"%d ",fp->WT(0).n());
					}

					if( (face_type::OBJ_TYPE & face_type::OBJ_TYPE_C) && (pi.mask & PLYMask::PM_FACECOLOR)  )
					{
						float t[3];
						t[0] = float(fp->C()[0])/255;
						t[1] = float(fp->C()[1])/255;
						t[2] = float(fp->C()[2])/255;
						fprintf(fpout,"9 ");
						fprintf(fpout,"%g %g %g ",t[0],t[1],t[2]);
						fprintf(fpout,"%g %g %g ",t[0],t[1],t[2]);
						fprintf(fpout,"%g %g %g ",t[0],t[1],t[2]);
					}
					else if( (face_type::OBJ_TYPE & face_type::OBJ_TYPE_WC) && (pi.mask & PLYMask::PM_WEDGCOLOR)  )
					{
						fprintf(fpout,"9 ");
						for(int z=0;z<3;++z)
							fprintf(fpout,"%g %g %g "
								,double(fp->WC(z)[0])/255
								,double(fp->WC(z)[1])/255
								,double(fp->WC(z)[2])/255
							);
					}

					if( (face_type::OBJ_TYPE & face_type::OBJ_TYPE_Q) && (pi.mask & PLYMask::PM_FACEQUALITY) )
						fprintf(fpout,"%g ",fp->Q());

					for(i=0;i<fdn;i++)
					{
						float tf;
						int ti;
						switch (FaceData[i].memtype1)
						{
						case T_FLOAT	:		tf=*( (float  *)        (((char *)fp)+FaceData[i].offset1));	fprintf(fpout,"%g ",tf); break;
						case T_DOUBLE   :		tf=*( (double *)        (((char *)fp)+FaceData[i].offset1));	fprintf(fpout,"%g ",tf); break;
						case T_INT		:		ti=*( (int    *)        (((char *)fp)+FaceData[i].offset1));	fprintf(fpout,"%i ",ti); break;
						case T_SHORT	:		ti=*( (short  *)        (((char *)fp)+FaceData[i].offset1));	fprintf(fpout,"%i ",ti); break;
						case T_CHAR		:		ti=*( (char   *)        (((char *)fp)+FaceData[i].offset1));	fprintf(fpout,"%i ",ti); break;
						case T_UCHAR	:		ti=*( (unsigned char *) (((char *)fp)+FaceData[i].offset1));	fprintf(fpout,"%i ",ti); break;
						default : assert(0);
						}
					}

					fprintf(fpout,"\n");
				}			  
			}
		}
	assert(fcnt==fn);
	fclose(fpout); 
	
	// Recupera i flag originali
	for(j=0,vi=vert.begin();vi!=vert.end();++vi)
		(*vi).Supervisor_Flags()=FlagV[j++]; 
	
	return 0;
}




}; // end class



} // end namespace tri
} // end namespace io
} // end namespace vcg

#endif