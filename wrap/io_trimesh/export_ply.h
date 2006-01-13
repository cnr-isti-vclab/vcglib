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
Revision 1.12  2006/01/10 13:20:42  cignoni
Changed ply::PlyMask to io::Mask

Revision 1.11  2005/11/23 15:48:25  pietroni
changed shot::similarity to shot::Similarity() and shot::camera to shot::Camera()

Revision 1.10  2004/10/28 00:52:45  cignoni
Better Doxygen documentation

Revision 1.9  2004/10/27 09:33:10  ganovelli
cast from scalar type to float added

Revision 1.8  2004/10/07 14:19:06  ganovelli
shot/camera io added

Revision 1.7  2004/07/15 10:54:48  ganovelli
std added

Revision 1.6  2004/05/28 14:11:13  ganovelli
changes to comply io_mask moving in vcg::ply namesp

Revision 1.5  2004/05/12 14:43:36  cignoni
removed warning of unused variables

Revision 1.4  2004/05/12 10:19:30  ganovelli
new line added at the end of file

Revision 1.3  2004/03/18 15:30:46  cignoni
Removed float/double warning

Revision 1.2  2004/03/09 21:26:47  cignoni
cr lf mismatch

Revision 1.1  2004/03/08 09:21:33  cignoni
Initial commit

Revision 1.1  2004/03/03 15:00:51  cignoni
Initial commit

****************************************************************************/

/**
@name Load and Save in Ply format
*/
//@{

#ifndef __VCGLIB_EXPORT_PLY
#define __VCGLIB_EXPORT_PLY

//#include<wrap/ply/io_mask.h>
#include<wrap/io_trimesh/io_mask.h>
#include<wrap/io_trimesh/io_ply.h>

#include <stdio.h>

namespace vcg {
namespace tri {
namespace io {


template <class SaveMeshType>
class ExporterPLY
{
  // Si occupa di convertire da un tipo all'altro.
// usata nella saveply per matchare i tipi tra stotype e memtype.
// Ad es se in memoria c'e' un int e voglio salvare un float
// src sara in effetti un puntatore a int il cui valore deve 
// essere convertito al tipo di ritorno desiderato (stotype)

template <class StoType> 
static void PlyConv(int mem_type, void *src, StoType &dest)
{
		switch (mem_type){
				case ply::T_FLOAT	:		dest = (StoType) (*  ((float  *) src)); break;
				case ply::T_DOUBLE:		dest = (StoType) (*  ((double *) src)); break;
				case ply::T_INT		:		dest = (StoType) (*  ((int    *) src)); break;
				case ply::T_SHORT	:		dest = (StoType) (*  ((short  *) src)); break;
				case ply::T_CHAR	:		dest = (StoType) (*  ((char   *) src)); break;
				case ply::T_UCHAR	:		dest = (StoType) (*  ((unsigned char *)src)); break;
			 	default : assert(0);
		}
}

public:
typedef ::vcg::ply::PropDescriptor PropDescriptor ;
typedef typename SaveMeshType::VertexPointer VertexPointer;
typedef typename SaveMeshType::ScalarType ScalarType;
typedef typename SaveMeshType::VertexType VertexType;
typedef typename SaveMeshType::FaceType FaceType;
typedef typename SaveMeshType::FacePointer FacePointer;
typedef typename SaveMeshType::VertexIterator VertexIterator;
typedef typename SaveMeshType::FaceIterator FaceIterator;

static int Save(SaveMeshType &m, const char * filename, bool binary=true)
{
  PlyInfo pi;
  return Save(m,filename,binary,pi);
}

static int Save(SaveMeshType &m,  const char * filename, int savemask )
{
	PlyInfo pi;
  pi.mask=savemask;
  return Save(m,filename,true,pi);
}

static int Save(SaveMeshType &m,  const char * filename, bool binary, PlyInfo &pi )	// V1.0
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

	if( pi.mask & Mask::IOM_WEDGTEXCOORD )
	{
		const char * TFILE = "TextureFile";

		for(i=0;i<m.textures.size();++i)
			fprintf(fpout,"comment %s %s\n", TFILE, (const char *)(m.textures[i].c_str()) );

		if(m.textures.size()>1 && (m.HasPerWedgeTexture() || m.HasPerVertexTexture())) multit = true;
	}

	if( (pi.mask & Mask::IOM_CAMERA) && m.shot.IsValid())
	 {
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
	} 

	fprintf(fpout,
		"element vertex %d\n"
		"property float x\n"
		"property float y\n"
		"property float z\n"
		,m.vn
	);


	if( pi.mask & Mask::IOM_VERTFLAGS )
	{
		fprintf(fpout,
			"property int flags\n"
		);
	}
	
	if( m.HasPerVertexColor()  && (pi.mask & Mask::IOM_VERTCOLOR) )
	{
		fprintf(fpout,
			"property uchar red\n"
			"property uchar green\n"
			"property uchar blue\n"
			"property uchar alpha\n"
		);
	}

	if( m.HasPerVertexQuality() && (pi.mask & Mask::IOM_VERTQUALITY) )
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

	if( pi.mask & Mask::IOM_FACEFLAGS )
	{
		fprintf(fpout,
			"property int flags\n"
		);
	}

	if( ( m.HasPerVertexTexture() && pi.mask & Mask::IOM_VERTTEXCOORD ) ||
	    ( m.HasPerWedgeTexture() && pi.mask & Mask::IOM_WEDGTEXCOORD ) )
	{
		fprintf(fpout,
			"property list uchar float texcoord\n"
		);

		if(multit)
			fprintf(fpout,
				"property int texnumber\n"
			);
	}

	if( m.HasPerFaceColor() && (pi.mask & Mask::IOM_FACECOLOR) )
	{
		fprintf(fpout,
			"property uchar red\n"
			"property uchar green\n"
			"property uchar blue\n"
			"property uchar alpha\n"
		);
	}
	
	if ( m.HasPerWedgeColor() && (pi.mask & Mask::IOM_WEDGCOLOR)  )
	{
		fprintf(fpout,
			"property list uchar float color\n"
		);
	}

	if( m.HasPerFaceQuality() && (pi.mask & Mask::IOM_FACEQUALITY) )
	{
		fprintf(fpout,
			"property float quality\n"
		);
	}

	for(i=0;i<pi.fdn;i++)
			fprintf(fpout,"property %s %s\n",pi.FaceData[i].stotypename(),pi.FaceData[i].propname);

	fprintf(fpout, "end_header\n"	);

		// Salvataggio camera
	 if( (pi.mask & Mask::IOM_CAMERA) && m.shot.IsValid() )
	 {
		 if(binary)
		 {
				float t[17];

				t[ 0] = (float) -m.shot.Similarity().tra[0];
				t[ 1] = (float)-m.shot.Similarity().tra[1];
				t[ 2] = (float)-m.shot.Similarity().tra[2];
				t[ 3] = (float)m.shot.Similarity().rot[0][0];
				t[ 4] = (float)m.shot.Similarity().rot[0][1];
				t[ 5] = (float)m.shot.Similarity().rot[0][2];
				t[ 6] = (float)m.shot.Similarity().rot[1][0];
				t[ 7] = (float)m.shot.Similarity().rot[1][1];
				t[ 8] = (float)m.shot.Similarity().rot[1][2];
				t[ 9] = (float)m.shot.Similarity().rot[2][0];
				t[10] = (float)m.shot.Similarity().rot[2][1];
				t[11] = (float)m.shot.Similarity().rot[2][2];
				t[12] = (float)m.shot.Camera().f;
				t[13] = (float)m.shot.Camera().s[0];
				t[14] = (float) m.shot.Camera().s[1];
				t[15] = (float)m.shot.Camera().c[0];
				t[16] = (float)m.shot.Camera().c[1];
				fwrite(t,sizeof(float),17,fpout);

				fwrite( &m.shot.Camera().viewport[0],sizeof(int),2,fpout );

				t[ 0] = (float)m.shot.Camera().k[0];
				t[ 1] = (float)m.shot.Camera().k[1];
				t[ 2] = (float)m.shot.Camera().k[2];
				t[ 3] = (float)m.shot.Camera().k[3];
				fwrite(t,sizeof(float),4,fpout);
		}
		else
		{
			fprintf(fpout,"%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %d %d %g %g %g %g\n"
				,-m.shot.Similarity().tra[0]
				,-m.shot.Similarity().tra[1]
				,-m.shot.Similarity().tra[2]
				,m.shot.Similarity().rot[0][0]
				,m.shot.Similarity().rot[0][1]
				,m.shot.Similarity().rot[0][2]
				,m.shot.Similarity().rot[1][0]
				,m.shot.Similarity().rot[1][1]
				,m.shot.Similarity().rot[1][2]
				,m.shot.Similarity().rot[2][0]
				,m.shot.Similarity().rot[2][1]
				,m.shot.Similarity().rot[2][2]
				,m.shot.Camera().f
				,m.shot.Camera().s[0]
				,m.shot.Camera().s[1]
				,m.shot.Camera().c[0]
				,m.shot.Camera().c[1]
				,m.shot.Camera().viewport[0]
				,m.shot.Camera().viewport[1]
				,m.shot.Camera().k[0]
				,m.shot.Camera().k[1]
				,m.shot.Camera().k[2]
				,m.shot.Camera().k[3]
			);
		}		
	}


	int j;
std::vector<int> FlagV; 
	VertexPointer  vp;
	VertexIterator vi;
	for(j=0,vi=m.vert.begin();vi!=m.vert.end();++vi)
	{
		vp=&(*vi);
		FlagV.push_back(vp->UberFlags()); // Salva in ogni caso flag del vertice
		if( ! vp->IsD() )
		{
			if(binary)
			{
				float t;

				t = float(vp->UberP()[0]); fwrite(&t,sizeof(float),1,fpout);
				t = float(vp->UberP()[1]); fwrite(&t,sizeof(float),1,fpout);
				t = float(vp->UberP()[2]); fwrite(&t,sizeof(float),1,fpout);
				
				if( pi.mask & Mask::IOM_VERTFLAGS )
					fwrite(&(vp->UberFlags()),sizeof(int),1,fpout);

				if( m.HasPerVertexColor() && (pi.mask & Mask::IOM_VERTCOLOR) )
					fwrite(&( vp->C() ),sizeof(char),4,fpout);

				if( m.HasPerVertexQuality() && (pi.mask & Mask::IOM_VERTQUALITY) )
					fwrite(&( vp->Q() ),sizeof(float),1,fpout);


				for(i=0;i<pi.vdn;i++)
				{
					double td(0); float tf(0);int ti;short ts; char tc; unsigned char tuc;
					switch (pi.VertexData[i].stotype1)
					{
          case ply::T_FLOAT	 :		PlyConv(pi.VertexData[i].memtype1,  ((char *)vp)+pi.VertexData[i].offset1, tf );	fwrite(&tf, sizeof(float),1,fpout); break;
					case ply::T_DOUBLE :		PlyConv(pi.VertexData[i].memtype1,  ((char *)vp)+pi.VertexData[i].offset1, td );	fwrite(&td, sizeof(double),1,fpout); break;
					case ply::T_INT		 :		PlyConv(pi.VertexData[i].memtype1,  ((char *)vp)+pi.VertexData[i].offset1, ti );	fwrite(&ti, sizeof(int),1,fpout); break;
					case ply::T_SHORT	 :		PlyConv(pi.VertexData[i].memtype1,  ((char *)vp)+pi.VertexData[i].offset1, ts );	fwrite(&ts, sizeof(short),1,fpout); break;
					case ply::T_CHAR	 :		PlyConv(pi.VertexData[i].memtype1,  ((char *)vp)+pi.VertexData[i].offset1, tc );	fwrite(&tc, sizeof(char),1,fpout); break;
					case ply::T_UCHAR	 :		PlyConv(pi.VertexData[i].memtype1,  ((char *)vp)+pi.VertexData[i].offset1, tuc);	fwrite(&tuc,sizeof(unsigned char),1,fpout); break;
					default : assert(0);
					}
				}
			}
			else 	// ***** ASCII *****
			{
				fprintf(fpout,"%g %g %g " ,vp->P()[0],vp->P()[1],vp->P()[2]);

				if( pi.mask & Mask::IOM_VERTFLAGS )
					fprintf(fpout,"%d ",vp->UberFlags());

				if( m.HasPerVertexColor() && (pi.mask & Mask::IOM_VERTCOLOR) )
					fprintf(fpout,"%d %d %d %d ",vp->C()[0],vp->C()[1],vp->C()[2],vp->C()[3] );

				if( m.HasPerVertexQuality() && (pi.mask & Mask::IOM_VERTQUALITY) )
					fprintf(fpout,"%g ",vp->Q());

				for(i=0;i<pi.vdn;i++)
				{
					float tf(0); double td(0);
					int ti;
					switch (pi.VertexData[i].memtype1)
					{
					case ply::T_FLOAT	 :		tf=*( (float  *)        (((char *)vp)+pi.VertexData[i].offset1));	fprintf(fpout,"%g ",tf); break;
					case ply::T_DOUBLE :    td=*( (double *)        (((char *)vp)+pi.VertexData[i].offset1));	fprintf(fpout,"%g ",tf); break;
					case ply::T_INT		 :		ti=*( (int    *)        (((char *)vp)+pi.VertexData[i].offset1));	fprintf(fpout,"%i ",ti); break;
					case ply::T_SHORT	 :		ti=*( (short  *)        (((char *)vp)+pi.VertexData[i].offset1)); fprintf(fpout,"%i ",ti); break;
					case ply::T_CHAR	 :		ti=*( (char   *)        (((char *)vp)+pi.VertexData[i].offset1));	fprintf(fpout,"%i ",ti); break;
					case ply::T_UCHAR	 :		ti=*( (unsigned char *) (((char *)vp)+pi.VertexData[i].offset1));	fprintf(fpout,"%i ",ti); break;
					default : assert(0);
					}
				}

				fprintf(fpout,"\n");
			}

			vp->UberFlags()=j; // Trucco! Nascondi nei flags l'indice del vertice non deletato!
			j++;
		}
	}
	assert(j==m.vn);

	char c = 3;
	unsigned char b9 = 9;
	unsigned char b6 = 6;
	FacePointer fp;
	int vv[3];
	FaceIterator fi;
	int fcnt=0;
	for(j=0,fi=m.face.begin();fi!=m.face.end();++fi)
		{
			fp=&(*fi);
			if( ! fp->IsD() )
			{ fcnt++;
				if(binary)
				{
					vv[0]=fp->cV(0)->UberFlags();
					vv[1]=fp->cV(1)->UberFlags();
					vv[2]=fp->cV(2)->UberFlags();
					fwrite(&c,1,1,fpout);
					fwrite(vv,sizeof(int),3,fpout);

					if( pi.mask & Mask::IOM_FACEFLAGS )
						fwrite(&(fp->Flags()),sizeof(int),1,fpout);

					if( m.HasPerVertexTexture() && (pi.mask & Mask::IOM_VERTTEXCOORD) )
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
					else if( m.HasPerWedgeTexture() && (pi.mask & Mask::IOM_WEDGTEXCOORD)  )
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
          
					if( m.HasPerFaceColor() && (pi.mask & Mask::IOM_FACECOLOR) )
					   fwrite(&( fp->C() ),sizeof(char),4,fpout);


					if( m.HasPerWedgeColor() && (pi.mask & Mask::IOM_WEDGCOLOR)  )
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

					if( m.HasPerFaceQuality() && (pi.mask & Mask::IOM_FACEQUALITY) )
						fwrite( &(fp->Q()),sizeof(float),1,fpout);


					for(i=0;i<pi.fdn;i++)
						{
						double td(0); float tf(0);int ti;short ts; char tc; unsigned char tuc;
						switch (pi.FaceData[i].stotype1){
								case ply::T_FLOAT	 :		PlyConv(pi.FaceData[i].memtype1,  ((char *)fp)+pi.FaceData[i].offset1, tf );	fwrite(&tf, sizeof(float),1,fpout); break;
								case ply::T_DOUBLE :		PlyConv(pi.FaceData[i].memtype1,  ((char *)fp)+pi.FaceData[i].offset1, td );	fwrite(&td, sizeof(double),1,fpout); break;
								case ply::T_INT		 :		PlyConv(pi.FaceData[i].memtype1,  ((char *)fp)+pi.FaceData[i].offset1, ti );	fwrite(&ti, sizeof(int),1,fpout); break;
								case ply::T_SHORT	 :		PlyConv(pi.FaceData[i].memtype1,  ((char *)fp)+pi.FaceData[i].offset1, ts );	fwrite(&ts, sizeof(short),1,fpout); break;
								case ply::T_CHAR	 :		PlyConv(pi.FaceData[i].memtype1,  ((char *)fp)+pi.FaceData[i].offset1, tc );	fwrite(&tc, sizeof(char),1,fpout); break;
								case ply::T_UCHAR	 :		PlyConv(pi.FaceData[i].memtype1,  ((char *)fp)+pi.FaceData[i].offset1, tuc);	fwrite(&tuc,sizeof(unsigned char),1,fpout); break;
								default : assert(0);
						}
					}
				}
				else	// ***** ASCII *****
				{
					fprintf(fpout,"3 %d %d %d ",
						fp->cV(0)->UberFlags(),	fp->cV(1)->UberFlags(), fp->cV(2)->UberFlags() );

					if( pi.mask & Mask::IOM_FACEFLAGS )
						fprintf(fpout,"%d ",fp->Flags());

					if( m.HasPerVertexTexture() && (pi.mask & Mask::IOM_VERTTEXCOORD) )
					{
						fprintf(fpout,"6 ");
						for(int k=0;k<3;++k)
							fprintf(fpout,"%g %g "
								,fp->V(k)->T().u()
								,fp->V(k)->T().v()
							);
					}
					else if( m.HasPerWedgeTexture() && (pi.mask & Mask::IOM_WEDGTEXCOORD)  )
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

					if( m.HasPerFaceColor() && (pi.mask & Mask::IOM_FACECOLOR)  )
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
					else if( m.HasPerWedgeColor() && (pi.mask & Mask::IOM_WEDGCOLOR)  )
					{
						fprintf(fpout,"9 ");
						for(int z=0;z<3;++z)
							fprintf(fpout,"%g %g %g "
								,double(fp->WC(z)[0])/255
								,double(fp->WC(z)[1])/255
								,double(fp->WC(z)[2])/255
							);
					}

					if( m.HasPerFaceQuality() && (pi.mask & Mask::IOM_FACEQUALITY) )
						fprintf(fpout,"%g ",fp->Q());

					for(i=0;i<pi.fdn;i++)
					{
						float tf(0); double td(0);
						int ti;
						switch (pi.FaceData[i].memtype1)
						{
						case  ply::T_FLOAT	:		tf=*( (float  *)        (((char *)fp)+pi.FaceData[i].offset1));	fprintf(fpout,"%g ",tf); break;
						case  ply::T_DOUBLE :		td=*( (double *)        (((char *)fp)+pi.FaceData[i].offset1));	fprintf(fpout,"%g ",tf); break;
						case  ply::T_INT		:		ti=*( (int    *)        (((char *)fp)+pi.FaceData[i].offset1));	fprintf(fpout,"%i ",ti); break;
						case  ply::T_SHORT	:		ti=*( (short  *)        (((char *)fp)+pi.FaceData[i].offset1));	fprintf(fpout,"%i ",ti); break;
						case  ply::T_CHAR		:		ti=*( (char   *)        (((char *)fp)+pi.FaceData[i].offset1));	fprintf(fpout,"%i ",ti); break;
						case  ply::T_UCHAR	:		ti=*( (unsigned char *) (((char *)fp)+pi.FaceData[i].offset1));	fprintf(fpout,"%i ",ti); break;
						default : assert(0);
						}
					}

					fprintf(fpout,"\n");
				}			  
			}
		}
	assert(fcnt==m.fn);
	fclose(fpout); 
	
	// Recupera i flag originali
	for(j=0,vi=m.vert.begin();vi!=m.vert.end();++vi)
		(*vi).UberFlags()=FlagV[j++]; 
	
	return 0;
}

static const char *ErrorMsg(int error)
{
  static std::vector<std::string> ply_error_msg;
  if(ply_error_msg.empty())
  {
    ply_error_msg.resize(PlyInfo::E_MAXPLYINFOERRORS );
    ply_error_msg[ply::E_NOERROR				]="No errors";
	  ply_error_msg[ply::E_CANTOPEN				]="Can't open file";
    ply_error_msg[ply::E_NOTHEADER ]="Header not found";
	  ply_error_msg[ply::E_UNESPECTEDEOF	]="Eof in header";
	  ply_error_msg[ply::E_NOFORMAT				]="Format not found";
	  ply_error_msg[ply::E_SYNTAX				]="Syntax error on header";
	  ply_error_msg[ply::E_PROPOUTOFELEMENT]="Property without element";
	  ply_error_msg[ply::E_BADTYPENAME		]="Bad type name";
	  ply_error_msg[ply::E_ELEMNOTFOUND		]="Element not found";
	  ply_error_msg[ply::E_PROPNOTFOUND		]="Property not found";
	  ply_error_msg[ply::E_BADTYPE				]="Bad type on addtoread";
	  ply_error_msg[ply::E_INCOMPATIBLETYPE]="Incompatible type";
	  ply_error_msg[ply::E_BADCAST				]="Bad cast";

    ply_error_msg[PlyInfo::E_NO_VERTEX      ]="No vertex field found";
    ply_error_msg[PlyInfo::E_NO_FACE        ]="No face field found";
	  ply_error_msg[PlyInfo::E_SHORTFILE      ]="Unespected eof";
	  ply_error_msg[PlyInfo::E_NO_3VERTINFACE ]="Face with more than 3 vertices";
	  ply_error_msg[PlyInfo::E_BAD_VERT_INDEX ]="Bad vertex index in face";
	  ply_error_msg[PlyInfo::E_NO_6TCOORD     ]="Face with no 6 texture coordinates";
	  ply_error_msg[PlyInfo::E_DIFFER_COLORS  ]="Number of color differ from vertices";
  }

  if(error>PlyInfo::E_MAXPLYINFOERRORS || error<0) return "Unknown error";
  else return ply_error_msg[error].c_str();
};



}; // end class



} // end namespace tri
} // end namespace io
} // end namespace vcg
//@}
#endif
