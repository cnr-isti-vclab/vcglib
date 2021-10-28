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

/**
@name Load and Save in Ply format
*/
//@{

#ifndef __VCGLIB_EXPORT_PLY
#define __VCGLIB_EXPORT_PLY

//#include<wrap/ply/io_mask.h>
#include<wrap/io_trimesh/io_mask.h>
#include<wrap/io_trimesh/io_ply.h>
#include<wrap/io_trimesh/precision.h>
#include<vcg/container/simple_temporary_data.h>
#include <vcg/complex/base.h>


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
	typedef typename SaveMeshType::ConstVertexPointer VertexPointer;
	typedef typename SaveMeshType::ScalarType ScalarType;
	typedef typename SaveMeshType::VertexType VertexType;
	typedef typename SaveMeshType::FaceType FaceType;
	typedef typename SaveMeshType::ConstFacePointer FacePointer;
	typedef typename SaveMeshType::ConstVertexIterator VertexIterator;
	typedef typename SaveMeshType::ConstFaceIterator FaceIterator;
	typedef typename SaveMeshType::ConstEdgeIterator EdgeIterator;
	typedef typename vcg::Shot<ScalarType>::ScalarType ShotScalarType;

	static int Save(const SaveMeshType &m, const char * filename, bool binary=true)
	{
		PlyInfo pi;
		return Save(m,filename,binary,pi);
	}

	static int Save(const SaveMeshType &m,  const char * filename, int savemask, bool binary = true, CallBackPos *cb=0 )
	{
		PlyInfo pi;
		pi.mask=savemask;
		return Save(m,filename,binary,pi,cb);
	}

	static int Save(const SaveMeshType &m,  const char * filename, bool binary, const PlyInfo &pi, CallBackPos *cb=0)	// V1.0
	{
		FILE * fpout;
		const char * hbin = "binary_little_endian";
		const char * hasc = "ascii";
		const char * h;
		//Coord ScalarType
		const int DGT = vcg::tri::io::Precision<ScalarType>::digits();
		const int DGTS = vcg::tri::io::Precision<ShotScalarType>::digits();
		const int DGTVQ = vcg::tri::io::Precision<typename VertexType::QualityType>::digits();
		const int DGTVR = vcg::tri::io::Precision<typename VertexType::RadiusType>::digits();
		const int DGTVT = vcg::tri::io::Precision<typename VertexType::TexCoordType::ScalarType>::digits();
		const int DGTFQ = vcg::tri::io::Precision<typename FaceType::QualityType>::digits();
		const int DGTFT = vcg::tri::io::Precision<typename FaceType::TexCoordType::ScalarType>::digits();
		bool saveTexIndexFlag = false;

		if(binary) h=hbin;
		else       h=hasc;

		fpout = fopen(filename,"wb");
		if(fpout==NULL)	{
			//pi.status=::vcg::ply::E_CANTOPEN;
			return ::vcg::ply::E_CANTOPEN;
		}
		fprintf(fpout,
				"ply\n"
				"format %s 1.0\n"
				"comment VCGLIB generated\n" ,
				h);

		if (((pi.mask & Mask::IOM_WEDGTEXCOORD) != 0) || ((pi.mask & Mask::IOM_VERTTEXCOORD) != 0))
		{
			const char * TFILE = "TextureFile";

			for(size_t i=0; i < m.textures.size(); ++i)
				fprintf(fpout,"comment %s %s\n", TFILE, (const char *)(m.textures[i].c_str()) );

			if(m.textures.size()>1 && (HasPerWedgeTexCoord(m) || HasPerVertexTexCoord(m))) saveTexIndexFlag = true;
		}

		if((pi.mask & Mask::IOM_CAMERA))
		{
			const char* cmtp = vcg::tri::io::Precision<ShotScalarType>::typeName();
			fprintf(fpout,"element camera 1\n");
			fprintf(fpout,"property %s view_px\n",cmtp);
			fprintf(fpout,"property %s view_py\n",cmtp);
			fprintf(fpout,"property %s view_pz\n",cmtp);
			fprintf(fpout,"property %s x_axisx\n",cmtp);
			fprintf(fpout,"property %s x_axisy\n",cmtp);
			fprintf(fpout,"property %s x_axisz\n",cmtp);
			fprintf(fpout,"property %s y_axisx\n",cmtp);
			fprintf(fpout,"property %s y_axisy\n",cmtp);
			fprintf(fpout,"property %s y_axisz\n",cmtp);
			fprintf(fpout,"property %s z_axisx\n",cmtp);
			fprintf(fpout,"property %s z_axisy\n",cmtp);
			fprintf(fpout,"property %s z_axisz\n",cmtp);
			fprintf(fpout,"property %s focal\n",cmtp);
			fprintf(fpout,"property %s scalex\n",cmtp);
			fprintf(fpout,"property %s scaley\n",cmtp);
			fprintf(fpout,"property %s centerx\n",cmtp);
			fprintf(fpout,"property %s centery\n",cmtp);
			fprintf(fpout,"property int viewportx\n");
			fprintf(fpout,"property int viewporty\n");
			fprintf(fpout,"property %s k1\n",cmtp);
			fprintf(fpout,"property %s k2\n",cmtp);
			fprintf(fpout,"property %s k3\n",cmtp);
			fprintf(fpout,"property %s k4\n",cmtp);
		}

		const char* vttp = vcg::tri::io::Precision<ScalarType>::typeName();
		fprintf(fpout,"element vertex %d\n",m.vn);
		fprintf(fpout,"property %s x\n",vttp);
		fprintf(fpout,"property %s y\n",vttp);
		fprintf(fpout,"property %s z\n",vttp);

		if( HasPerVertexNormal(m) &&( pi.mask & Mask::IOM_VERTNORMAL) )
		{
			fprintf(fpout,"property %s nx\n",vttp);
			fprintf(fpout,"property %s ny\n",vttp);
			fprintf(fpout,"property %s nz\n",vttp);
		}


		if( HasPerVertexFlags(m) &&( pi.mask & Mask::IOM_VERTFLAGS) )
		{
			fprintf(fpout,
					"property int flags\n");
		}

		if( HasPerVertexColor(m)  && (pi.mask & Mask::IOM_VERTCOLOR) )
		{
			fprintf(fpout,
					"property uchar red\n"
					"property uchar green\n"
					"property uchar blue\n"
					"property uchar alpha\n");
		}

		if( HasPerVertexQuality(m) && (pi.mask & Mask::IOM_VERTQUALITY) )
		{
			const char* vqtp = vcg::tri::io::Precision<typename VertexType::ScalarType>::typeName();
			fprintf(fpout,"property %s quality\n",vqtp);
		}

		if( tri::HasPerVertexRadius(m) && (pi.mask & Mask::IOM_VERTRADIUS) )
		{
			const char* rdtp = vcg::tri::io::Precision<typename VertexType::RadiusType>::typeName();
			fprintf(fpout,"property %s radius\n",rdtp);
		}
		if( ( HasPerVertexTexCoord(m) && pi.mask & Mask::IOM_VERTTEXCOORD ) )
		{
			const char* rdtp = vcg::tri::io::Precision<typename VertexType::TexCoordType::ScalarType>::typeName();
			fprintf(fpout,
					"property %s texture_u\n"
					"property %s texture_v\n",
					rdtp, rdtp);
		}
		for(size_t i=0;i<pi.VertDescriptorVec.size();i++){
			if (!pi.VertDescriptorVec[i].islist) {
				fprintf(
					fpout,
					"property %s %s\n",
					pi.VertDescriptorVec[i].stotypename(),
					pi.VertDescriptorVec[i].propname.c_str());
			}
			else {
				fprintf(
					fpout,
					"property list %s %s %s\n",
					pi.VertDescriptorVec[i].stotype2name(),
					pi.VertDescriptorVec[i].stotypename(),
					pi.VertDescriptorVec[i].propname.c_str());
			}
		}

		fprintf(fpout,
				"element face %d\n"
				"property list uchar int vertex_indices\n",
				m.fn );

		if(HasPerFaceFlags(m)   && (pi.mask & Mask::IOM_FACEFLAGS) )
		{
			fprintf(fpout,
					"property int flags\n");
		}

		if( (HasPerWedgeTexCoord(m)  && pi.mask & Mask::IOM_WEDGTEXCOORD ) ||
				(HasPerVertexTexCoord(m) && (!HasPerWedgeTexCoord(m)) && pi.mask & Mask::IOM_WEDGTEXCOORD ) )  // Note that you can save VT as WT if you really really want it...
		{
			const char* rdtp = vcg::tri::io::Precision<typename FaceType::TexCoordType::ScalarType>::typeName();
			fprintf(fpout,
					"property list uchar %s texcoord\n", rdtp );
		}
		// The texture index information has to be saved for each face (if necessary) both for PerVert and PerWedg
		if( saveTexIndexFlag &&
				( ( HasPerWedgeTexCoord(m)  && (pi.mask & Mask::IOM_WEDGTEXCOORD) ) ||
				  ( HasPerVertexTexCoord(m) && (pi.mask & Mask::IOM_VERTTEXCOORD) ) ||
				  ( HasPerVertexTexCoord(m) && (!HasPerWedgeTexCoord(m)) && (pi.mask & Mask::IOM_WEDGTEXCOORD) )
				  )
				)
		{
			fprintf(fpout,
					"property int texnumber\n");
		}

		if( HasPerFaceColor(m) && (pi.mask & Mask::IOM_FACECOLOR) )
		{
			fprintf(
						fpout,
						"property uchar red\n"
						"property uchar green\n"
						"property uchar blue\n"
						"property uchar alpha\n");
		}

		if ( HasPerWedgeColor(m) && (pi.mask & Mask::IOM_WEDGCOLOR)  )
		{
			fprintf(
						fpout,
						"property list uchar float color\n");
		}

		if (HasPerFaceNormal(m) && (pi.mask & Mask::IOM_FACENORMAL))
		{
			const char* fntp = vcg::tri::io::Precision<typename SaveMeshType::ScalarType>::typeName();
			fprintf(fpout, "property %s nx\n", fntp);
			fprintf(fpout, "property %s ny\n", fntp);
			fprintf(fpout, "property %s nz\n", fntp);
		}

		if( HasPerFaceQuality(m) && (pi.mask & Mask::IOM_FACEQUALITY) )
		{
			const char* fqtp = vcg::tri::io::Precision<typename SaveMeshType::FaceType::QualityType>::typeName();
			fprintf(fpout,"property %s quality\n",fqtp);
		}

		for(size_t i=0;i<pi.FaceDescriptorVec.size();i++) {
			if (!pi.FaceDescriptorVec[i].islist){
				fprintf(
					fpout,
					"property %s %s\n",
					pi.FaceDescriptorVec[i].stotypename(),
					pi.FaceDescriptorVec[i].propname.c_str());
			}
			else {
				fprintf(
					fpout,
					"property list %s %s %s\n",
					pi.FaceDescriptorVec[i].stotype2name(),
					pi.FaceDescriptorVec[i].stotypename(),
					pi.FaceDescriptorVec[i].propname.c_str());
			}
		}

		// Saving of edges is enabled if requested
		if( m.en>0 && (pi.mask & Mask::IOM_EDGEINDEX) )
			fprintf(
						fpout,
						"element edge %d\n" "property int vertex1\n""property int vertex2\n",m.en);
		fprintf(fpout, "end_header\n"	);

		// Salvataggio camera
		if((pi.mask & Mask::IOM_CAMERA))
		{
			if(binary)
			{
				ShotScalarType t[17];

				t[ 0] = (ShotScalarType)m.shot.Extrinsics.Tra()[0];
				t[ 1] = (ShotScalarType)m.shot.Extrinsics.Tra()[1];
				t[ 2] = (ShotScalarType)m.shot.Extrinsics.Tra()[2];
				t[ 3] = (ShotScalarType)m.shot.Extrinsics.Rot()[0][0];
				t[ 4] = (ShotScalarType)m.shot.Extrinsics.Rot()[0][1];
				t[ 5] = (ShotScalarType)m.shot.Extrinsics.Rot()[0][2];
				t[ 6] = (ShotScalarType)m.shot.Extrinsics.Rot()[1][0];
				t[ 7] = (ShotScalarType)m.shot.Extrinsics.Rot()[1][1];
				t[ 8] = (ShotScalarType)m.shot.Extrinsics.Rot()[1][2];
				t[ 9] = (ShotScalarType)m.shot.Extrinsics.Rot()[2][0];
				t[10] = (ShotScalarType)m.shot.Extrinsics.Rot()[2][1];
				t[11] = (ShotScalarType)m.shot.Extrinsics.Rot()[2][2];
				t[12] = (ShotScalarType)m.shot.Intrinsics.FocalMm;
				t[13] = (ShotScalarType)m.shot.Intrinsics.PixelSizeMm[0];
				t[14] = (ShotScalarType)m.shot.Intrinsics.PixelSizeMm[1];
				t[15] = (ShotScalarType)m.shot.Intrinsics.CenterPx[0];
				t[16] = (ShotScalarType)m.shot.Intrinsics.CenterPx[1];
				fwrite(t,sizeof(ShotScalarType),17,fpout);

				fwrite( &m.shot.Intrinsics.ViewportPx[0],sizeof(int),2,fpout );

				t[ 0] = (ShotScalarType)m.shot.Intrinsics.k[0];
				t[ 1] = (ShotScalarType)m.shot.Intrinsics.k[1];
				t[ 2] = (ShotScalarType)m.shot.Intrinsics.k[2];
				t[ 3] = (ShotScalarType)m.shot.Intrinsics.k[3];
				fwrite(t,sizeof(ShotScalarType),4,fpout);
			}
			else
			{
				fprintf(fpout,"%.*g %.*g %.*g %.*g %.*g %.*g %.*g %.*g %.*g %.*g %.*g %.*g %.*g %.*g %.*g %.*g %.*g %d %d %.*g %.*g %.*g %.*g\n"
							,DGTS,-m.shot.Extrinsics.Tra()[0]
						,DGTS,-m.shot.Extrinsics.Tra()[1]
						,DGTS,-m.shot.Extrinsics.Tra()[2]
						,DGTS,m.shot.Extrinsics.Rot()[0][0]
						,DGTS,m.shot.Extrinsics.Rot()[0][1]
						,DGTS,m.shot.Extrinsics.Rot()[0][2]
						,DGTS,m.shot.Extrinsics.Rot()[1][0]
						,DGTS,m.shot.Extrinsics.Rot()[1][1]
						,DGTS,m.shot.Extrinsics.Rot()[1][2]
						,DGTS,m.shot.Extrinsics.Rot()[2][0]
						,DGTS,m.shot.Extrinsics.Rot()[2][1]
						,DGTS,m.shot.Extrinsics.Rot()[2][2]
						,DGTS,m.shot.Intrinsics.FocalMm
						,DGTS,m.shot.Intrinsics.PixelSizeMm[0]
						,DGTS,m.shot.Intrinsics.PixelSizeMm[1]
						,DGTS,m.shot.Intrinsics.CenterPx[0]
						,DGTS,m.shot.Intrinsics.CenterPx[1]
						,m.shot.Intrinsics.ViewportPx[0]
						,m.shot.Intrinsics.ViewportPx[1]
						,DGTS,m.shot.Intrinsics.k[0]
						,DGTS,m.shot.Intrinsics.k[1]
						,DGTS,m.shot.Intrinsics.k[2]
						,DGTS,m.shot.Intrinsics.k[3]
						);
			}
		}


		int j;
		std::vector<int> FlagV;
		VertexPointer  vp;
		VertexIterator vi;
		SimpleTempData<typename SaveMeshType::VertContainer,int> indices(m.vert);

		std::vector<typename SaveMeshType:: template ConstPerVertexAttributeHandle<float > > thfv(pi.VertDescriptorVec.size());
		std::vector<typename SaveMeshType:: template ConstPerVertexAttributeHandle<double> > thdv(pi.VertDescriptorVec.size());
		std::vector<typename SaveMeshType:: template ConstPerVertexAttributeHandle<int   > > thiv(pi.VertDescriptorVec.size());
		std::vector<typename SaveMeshType:: template ConstPerVertexAttributeHandle<short > > thsv(pi.VertDescriptorVec.size());
		std::vector<typename SaveMeshType:: template ConstPerVertexAttributeHandle<char  > > thcv(pi.VertDescriptorVec.size());
		std::vector<typename SaveMeshType:: template ConstPerVertexAttributeHandle<unsigned char> >  thuv(pi.VertDescriptorVec.size());
		std::vector<typename SaveMeshType:: template ConstPerVertexAttributeHandle<vcg::Point3f> > thp3fv(pi.VertDescriptorVec.size());
		std::vector<typename SaveMeshType:: template ConstPerVertexAttributeHandle<vcg::Point3d> > thp3dv(pi.VertDescriptorVec.size());

		for(size_t i=0;i<pi.VertDescriptorVec.size();i++)
		{
			if(!pi.VertAttrNameVec.empty() && !pi.VertAttrNameVec[i].empty())
			{ // trying to use named attribute to retrieve the value to store
				assert(vcg::tri::HasPerVertexAttribute(m,pi.VertAttrNameVec[i]));
				if (!pi.VertDescriptorVec[i].islist){
					switch (pi.VertDescriptorVec[i].stotype1)
					{
					case ply::T_FLOAT  : thfv[i] = vcg::tri::Allocator<SaveMeshType>::template FindPerVertexAttribute<float>(m,pi.VertAttrNameVec[i]); break;
					case ply::T_DOUBLE : thdv[i] = vcg::tri::Allocator<SaveMeshType>::template FindPerVertexAttribute<double>(m,pi.VertAttrNameVec[i]); break;
					case ply::T_INT    : thiv[i] = vcg::tri::Allocator<SaveMeshType>::template FindPerVertexAttribute<int   >(m,pi.VertAttrNameVec[i]); break;
					case ply::T_SHORT  : thsv[i] = vcg::tri::Allocator<SaveMeshType>::template FindPerVertexAttribute<short >(m,pi.VertAttrNameVec[i]); break;
					case ply::T_CHAR   : thcv[i] = vcg::tri::Allocator<SaveMeshType>::template FindPerVertexAttribute<char>(m,pi.VertAttrNameVec[i]); break;
					case ply::T_UCHAR  : thuv[i] = vcg::tri::Allocator<SaveMeshType>::template FindPerVertexAttribute<unsigned char>(m,pi.VertAttrNameVec[i]); break;
					default : assert(0);
					}
				}
				else {
					switch (pi.VertDescriptorVec[i].stotype1)
					{
					case ply::T_FLOAT  : thp3fv[i] = vcg::tri::Allocator<SaveMeshType>::template FindPerVertexAttribute<vcg::Point3f>(m,pi.VertAttrNameVec[i]); break;
					case ply::T_DOUBLE : thp3dv[i] = vcg::tri::Allocator<SaveMeshType>::template FindPerVertexAttribute<vcg::Point3d>(m,pi.VertAttrNameVec[i]); break;
					default : assert(0);
					}
				}
			}
		}
		std::vector<typename SaveMeshType:: template ConstPerFaceAttributeHandle<float > > thff(pi.FaceDescriptorVec.size());
		std::vector<typename SaveMeshType:: template ConstPerFaceAttributeHandle<double> > thdf(pi.FaceDescriptorVec.size());
		std::vector<typename SaveMeshType:: template ConstPerFaceAttributeHandle<int   > > thif(pi.FaceDescriptorVec.size());
		std::vector<typename SaveMeshType:: template ConstPerFaceAttributeHandle<short > > thsf(pi.FaceDescriptorVec.size());
		std::vector<typename SaveMeshType:: template ConstPerFaceAttributeHandle<char  > > thcf(pi.FaceDescriptorVec.size());
		std::vector<typename SaveMeshType:: template ConstPerFaceAttributeHandle<unsigned char> >  thuf(pi.FaceDescriptorVec.size());
		std::vector<typename SaveMeshType:: template ConstPerFaceAttributeHandle<vcg::Point3f> > thp3ff(pi.FaceDescriptorVec.size());
		std::vector<typename SaveMeshType:: template ConstPerFaceAttributeHandle<vcg::Point3d> > thp3df(pi.FaceDescriptorVec.size());

		for(size_t i=0;i<pi.FaceDescriptorVec.size();i++)
		{
			if(!pi.FaceAttrNameVec.empty() && !pi.FaceAttrNameVec[i].empty())
			{ // trying to use named attribute to retrieve the value to store
				assert(vcg::tri::HasPerFaceAttribute(m,pi.FaceAttrNameVec[i]));
				if (!pi.FaceDescriptorVec[i].islist) {
					switch (pi.FaceDescriptorVec[i].stotype1)
					{
					case ply::T_FLOAT  : thff[i] = vcg::tri::Allocator<SaveMeshType>::template FindPerFaceAttribute<float>(m,pi.FaceAttrNameVec[i]); break;
					case ply::T_DOUBLE : thdf[i] = vcg::tri::Allocator<SaveMeshType>::template FindPerFaceAttribute<double>(m,pi.FaceAttrNameVec[i]); break;
					case ply::T_INT    : thif[i] = vcg::tri::Allocator<SaveMeshType>::template FindPerFaceAttribute<int   >(m,pi.FaceAttrNameVec[i]); break;
					case ply::T_SHORT  : thsf[i] = vcg::tri::Allocator<SaveMeshType>::template FindPerFaceAttribute<short >(m,pi.FaceAttrNameVec[i]); break;
					case ply::T_CHAR   : thcf[i] = vcg::tri::Allocator<SaveMeshType>::template FindPerFaceAttribute<char>(m,pi.FaceAttrNameVec[i]); break;
					case ply::T_UCHAR  : thuf[i] = vcg::tri::Allocator<SaveMeshType>::template FindPerFaceAttribute<unsigned char>(m,pi.FaceAttrNameVec[i]); break;
					default : assert(0);
					}
				}
				else {
					switch (pi.FaceDescriptorVec[i].stotype1)
					{
					case ply::T_FLOAT  : thp3ff[i] = vcg::tri::Allocator<SaveMeshType>::template FindPerFaceAttribute<vcg::Point3f>(m,pi.FaceAttrNameVec[i]); break;
					case ply::T_DOUBLE : thp3df[i] = vcg::tri::Allocator<SaveMeshType>::template FindPerFaceAttribute<vcg::Point3d>(m,pi.FaceAttrNameVec[i]); break;
					default : assert(0);
					}
				}
			}
		}

		for(j=0,vi=m.vert.begin();vi!=m.vert.end();++vi){
			vp=&(*vi);
			indices[vi] = j;
			//((m.vn+m.fn) != 0) all vertices and faces have been marked as deleted but the are still in the vert/face vectors
			if(cb && ((j%1000)==0) && ((m.vn+m.fn) != 0) )(*cb)( (100*j)/(m.vn+m.fn), "Saving Vertices");

			if( !HasPerVertexFlags(m) || !vp->IsD() )
			{
				if(binary)
				{
					ScalarType t;

					t = ScalarType(vp->P()[0]); fwrite(&t,sizeof(ScalarType),1,fpout);
					t = ScalarType(vp->P()[1]); fwrite(&t,sizeof(ScalarType),1,fpout);
					t = ScalarType(vp->P()[2]); fwrite(&t,sizeof(ScalarType),1,fpout);

					if( HasPerVertexNormal(m) && (pi.mask & Mask::IOM_VERTNORMAL) )
					{
						t = ScalarType(vp->N()[0]); fwrite(&t,sizeof(ScalarType),1,fpout);
						t = ScalarType(vp->N()[1]); fwrite(&t,sizeof(ScalarType),1,fpout);
						t = ScalarType(vp->N()[2]); fwrite(&t,sizeof(ScalarType),1,fpout);
					}
					if( HasPerVertexFlags(m) && (pi.mask & Mask::IOM_VERTFLAGS) )
						fwrite(&(vp->Flags()),sizeof(int),1,fpout);

					if( HasPerVertexColor(m) && (pi.mask & Mask::IOM_VERTCOLOR) ){
						auto c = vp->C();
						fwrite(&c,sizeof(char),4,fpout);
					}

					if( HasPerVertexQuality(m) && (pi.mask & Mask::IOM_VERTQUALITY) ){
						auto q = vp->Q();
						fwrite(&q, sizeof(typename VertexType::QualityType),1,fpout);
					}

					if( HasPerVertexRadius(m) && (pi.mask & Mask::IOM_VERTRADIUS) ){
						auto r = vp->R();
						fwrite(&r,sizeof(typename VertexType::RadiusType),1,fpout);
					}

					if( HasPerVertexTexCoord(m) && (pi.mask & Mask::IOM_VERTTEXCOORD) ){
						typename VertexType::TexCoordType::ScalarType t;
						t = ScalarType(vp->T().u()); fwrite(&t,sizeof(typename VertexType::TexCoordType::ScalarType),1,fpout);
						t = ScalarType(vp->T().v()); fwrite(&t,sizeof(typename VertexType::TexCoordType::ScalarType),1,fpout);
					}

					for(size_t i=0;i<pi.VertDescriptorVec.size();i++)
					{
						double td(0); float tf(0);int ti;short ts; char tc; unsigned char tu;
						if(!pi.VertAttrNameVec.empty() && !pi.VertAttrNameVec[i].empty())
						{ // trying to use named attribute to retrieve the value to store
							assert(vcg::tri::HasPerVertexAttribute(m,pi.VertAttrNameVec[i]));
							if (!pi.VertDescriptorVec[i].islist){
								switch (pi.VertDescriptorVec[i].stotype1)
								{
								case ply::T_FLOAT  : tf=thfv[i][vp]; fwrite(&tf, sizeof(float),1,fpout); break;
								case ply::T_DOUBLE : td=thdv[i][vp]; fwrite(&td, sizeof(double),1,fpout); break;
								case ply::T_INT    : ti=thiv[i][vp]; fwrite(&ti, sizeof(int),1,fpout); break;
								case ply::T_SHORT  : ts=thsv[i][vp]; fwrite(&ts, sizeof(short),1,fpout); break;
								case ply::T_CHAR   : tc=thcv[i][vp]; fwrite(&tc, sizeof(char),1,fpout); break;
								case ply::T_UCHAR  : tu=thuv[i][vp]; fwrite(&tu,sizeof(unsigned char),1,fpout); break;
								default : assert(0);
								}
							}
							else { //it is a Poin3f or a Point3d attribute. Saving it as a list
								static const unsigned char psize = 3;
								switch (pi.VertDescriptorVec[i].stotype1)
								{
								case ply::T_FLOAT  :
									fwrite(&psize, sizeof(unsigned char), 1,fpout);
									fwrite(&thp3fv[i][vp][0], sizeof(float), 1,fpout);
									fwrite(&thp3fv[i][vp][1], sizeof(float), 1,fpout);
									fwrite(&thp3fv[i][vp][2], sizeof(float), 1,fpout);
									break;
									//fprintf(fpout,"%d %f %f %f", 3, thp3fv[i][vp][0], thp3fv[i][vp][1], thp3fv[i][vp][2]); break;
								case ply::T_DOUBLE :
									fwrite(&psize, sizeof(unsigned char), 1,fpout);
									fwrite(&thp3dv[i][vp][0], sizeof(double), 1,fpout);
									fwrite(&thp3dv[i][vp][1], sizeof(double), 1,fpout);
									fwrite(&thp3dv[i][vp][2], sizeof(double), 1,fpout);
									break;
									//fprintf(fpout,"%d %lf %lf %lf", 3, thp3dv[i][vp][0], thp3dv[i][vp][1], thp3dv[i][vp][2]); break;
								default : assert(0);
								}
							}
						}
						else
						{
							switch (pi.VertDescriptorVec[i].stotype1)
							{
							case ply::T_FLOAT	 :		PlyConv(pi.VertDescriptorVec[i].memtype1,  ((char *)vp)+pi.VertDescriptorVec[i].offset1, tf );	fwrite(&tf, sizeof(float),1,fpout); break;
							case ply::T_DOUBLE :		PlyConv(pi.VertDescriptorVec[i].memtype1,  ((char *)vp)+pi.VertDescriptorVec[i].offset1, td );	fwrite(&td, sizeof(double),1,fpout); break;
							case ply::T_INT		 :		PlyConv(pi.VertDescriptorVec[i].memtype1,  ((char *)vp)+pi.VertDescriptorVec[i].offset1, ti );	fwrite(&ti, sizeof(int),1,fpout); break;
							case ply::T_SHORT	 :		PlyConv(pi.VertDescriptorVec[i].memtype1,  ((char *)vp)+pi.VertDescriptorVec[i].offset1, ts );	fwrite(&ts, sizeof(short),1,fpout); break;
							case ply::T_CHAR	 :		PlyConv(pi.VertDescriptorVec[i].memtype1,  ((char *)vp)+pi.VertDescriptorVec[i].offset1, tc );	fwrite(&tc, sizeof(char),1,fpout); break;
							case ply::T_UCHAR	 :		PlyConv(pi.VertDescriptorVec[i].memtype1,  ((char *)vp)+pi.VertDescriptorVec[i].offset1, tu );	fwrite(&tu,sizeof(unsigned char),1,fpout); break;
							default : assert(0);
							}
						}
					}
				}
				else 	// ***** ASCII *****
				{
					fprintf(fpout,"%.*g %.*g %.*g " ,DGT,vp->P()[0],DGT,vp->P()[1],DGT,vp->P()[2]);

					if( HasPerVertexNormal(m) && (pi.mask & Mask::IOM_VERTNORMAL) )
						fprintf(fpout,"%.*g %.*g %.*g " ,DGT,ScalarType(vp->N()[0]),DGT,ScalarType(vp->N()[1]),DGT,ScalarType(vp->N()[2]));

					if( HasPerVertexFlags(m) && (pi.mask & Mask::IOM_VERTFLAGS))
						fprintf(fpout,"%d ",vp->Flags());

					if( HasPerVertexColor(m) && (pi.mask & Mask::IOM_VERTCOLOR) )
						fprintf(fpout,"%d %d %d %d ",vp->C()[0],vp->C()[1],vp->C()[2],vp->C()[3] );

					if( HasPerVertexQuality(m) && (pi.mask & Mask::IOM_VERTQUALITY) )
						fprintf(fpout,"%.*g ",DGTVQ,vp->Q());

					if( HasPerVertexRadius(m) && (pi.mask & Mask::IOM_VERTRADIUS) )
						fprintf(fpout,"%.*g ",DGTVR,vp->R());

					if( HasPerVertexTexCoord(m) && (pi.mask & Mask::IOM_VERTTEXCOORD) )
						fprintf(fpout,"%.*g %.*g",DGTVT,vp->T().u(),vp->T().v());

					for(size_t i=0;i<pi.VertDescriptorVec.size();i++)
					{
						float tf(0); double td(0); int ti;
						if(!pi.VertAttrNameVec.empty() && !pi.VertAttrNameVec[i].empty())
						{ // trying to use named attribute to retrieve the value to store
							assert(vcg::tri::HasPerVertexAttribute(m,pi.VertAttrNameVec[i]));
							if (!pi.VertDescriptorVec[i].islist){
								switch (pi.VertDescriptorVec[i].stotype1)
								{
								case ply::T_FLOAT  : tf=thfv[i][vp]; fprintf(fpout,"%f ",tf); break;
								case ply::T_DOUBLE : td=thdv[i][vp]; fprintf(fpout,"%lf ",td); break;
								case ply::T_INT    : ti=thiv[i][vp]; fprintf(fpout,"%i ",ti); break;
								case ply::T_SHORT  : ti=thsv[i][vp]; fprintf(fpout,"%i ",ti); break;
								case ply::T_CHAR   : ti=thcv[i][vp]; fprintf(fpout,"%i ",ti); break;
								case ply::T_UCHAR  : ti=thuv[i][vp]; fprintf(fpout,"%i ",ti); break;
								default : assert(0);
								}
							}
							else { //it is a Poin3f or a Point3d attribute. Saving it as a list
								switch (pi.VertDescriptorVec[i].stotype1)
								{
								case ply::T_FLOAT  : fprintf(fpout,"%d %f %f %f", 3, thp3fv[i][vp][0], thp3fv[i][vp][1], thp3fv[i][vp][2]); break;
								case ply::T_DOUBLE : fprintf(fpout,"%d %lf %lf %lf", 3, thp3dv[i][vp][0], thp3dv[i][vp][1], thp3dv[i][vp][2]); break;
								default : assert(0);
								}
							}
						}
						else
						{
							switch (pi.VertDescriptorVec[i].memtype1)
							{
							case ply::T_FLOAT  : tf=*( (float  *)        (((char *)vp)+pi.VertDescriptorVec[i].offset1)); fprintf(fpout,"%f ",tf); break;
							case ply::T_DOUBLE : td=*( (double *)        (((char *)vp)+pi.VertDescriptorVec[i].offset1)); fprintf(fpout,"%lf ",tf); break;
							case ply::T_INT    : ti=*( (int	*)           (((char *)vp)+pi.VertDescriptorVec[i].offset1)); fprintf(fpout,"%i ",ti); break;
							case ply::T_SHORT  : ti=*( (short  *)        (((char *)vp)+pi.VertDescriptorVec[i].offset1)); fprintf(fpout,"%i ",ti); break;
							case ply::T_CHAR   : ti=*( (char   *)        (((char *)vp)+pi.VertDescriptorVec[i].offset1)); fprintf(fpout,"%i ",ti); break;
							case ply::T_UCHAR  : ti=*( (unsigned char *) (((char *)vp)+pi.VertDescriptorVec[i].offset1)); fprintf(fpout,"%i ",ti); break;
							default : assert(0);
							}
						}
					}

					fprintf(fpout,"\n");
				}
				j++;
			}
		}
		/*vcg::tri::*/
		// this assert triggers when the vn != number of vertexes in vert that are not deleted.
		assert(j==m.vn);

		unsigned char b3char = 3;
		unsigned char b9char = 9;
		unsigned char b6char = 6;
		FacePointer fp;
		int vv[3];
		FaceIterator fi;
		int fcnt=0;
		for(j=0,fi=m.face.begin();fi!=m.face.end();++fi)
		{
			//((m.vn+m.fn) != 0) all vertices and faces have been marked as deleted but the are still in the vert/face vectors
			if(cb && ((j%1000)==0) && ((m.vn+m.fn) != 0))
				(*cb)( 100*(m.vn+j)/(m.vn+m.fn), "Saving Vertices");

			fp=&(*fi);
			if( ! fp->IsD() )
			{ fcnt++;
				if(binary)
				{
					vv[0]=indices[fp->cV(0)];
					vv[1]=indices[fp->cV(1)];
					vv[2]=indices[fp->cV(2)];
					fwrite(&b3char,sizeof(char),1,fpout);
					fwrite(vv,sizeof(int),3,fpout);

					if(HasPerFaceFlags(m)&&( pi.mask & Mask::IOM_FACEFLAGS) ){
						auto fl = fp->Flags();
						fwrite(&fl,sizeof(int),1,fpout);
					}

					if( HasPerVertexTexCoord(m) && (!HasPerWedgeTexCoord(m)) && (pi.mask & Mask::IOM_WEDGTEXCOORD) )  // Note that you can save VT as WT if you really want it...
					{
						fwrite(&b6char,sizeof(char),1,fpout);
						typename FaceType::TexCoordType::ScalarType t[6];
						for(int k=0;k<3;++k)
						{
							t[k*2+0] = fp->V(k)->T().u();
							t[k*2+1] = fp->V(k)->T().v();
						}
						fwrite(t,sizeof(typename FaceType::TexCoordType::ScalarType),6,fpout);
					}
					else if( HasPerWedgeTexCoord(m) && (pi.mask & Mask::IOM_WEDGTEXCOORD)  )
					{
						fwrite(&b6char,sizeof(char),1,fpout);
						typename FaceType::TexCoordType::ScalarType t[6];
						for(int k=0;k<3;++k)
						{
							t[k*2+0] = fp->WT(k).u();
							t[k*2+1] = fp->WT(k).v();
						}
						fwrite(t,sizeof(typename FaceType::TexCoordType::ScalarType),6,fpout);
					}

					if(saveTexIndexFlag)
					{
						int t = fp->WT(0).n();
						fwrite(&t,sizeof(int),1,fpout);
					}

					if( HasPerFaceColor(m) && (pi.mask & Mask::IOM_FACECOLOR) )
						fwrite(&( fp->C() ),sizeof(char),4,fpout);


					if( HasPerWedgeColor(m) && (pi.mask & Mask::IOM_WEDGCOLOR)  )
					{
						fwrite(&b9char,sizeof(char),1,fpout);
						float t[3];
						for(int z=0;z<3;++z)
						{
							t[0] = float(fp->WC(z)[0])/255;
							t[1] = float(fp->WC(z)[1])/255;
							t[2] = float(fp->WC(z)[2])/255;
							fwrite( t,sizeof(float),3,fpout);
						}
					}

					if( HasPerFaceNormal(m) && (pi.mask & Mask::IOM_FACENORMAL) )
					{
						ScalarType t;
						t = ScalarType(fp->N()[0]); fwrite(&t,sizeof(ScalarType),1,fpout);
						t = ScalarType(fp->N()[1]); fwrite(&t,sizeof(ScalarType),1,fpout);
						t = ScalarType(fp->N()[2]); fwrite(&t,sizeof(ScalarType),1,fpout);
					}

					if( HasPerFaceQuality(m) && (pi.mask & Mask::IOM_FACEQUALITY) )
						fwrite( &(fp->Q()),sizeof(typename FaceType::ScalarType),1,fpout);


					for(size_t i=0;i<pi.FaceDescriptorVec.size();i++)
					{
						double td(0); float tf(0);int ti;short ts; char tc; unsigned char tu;
						if(!pi.FaceAttrNameVec.empty() && !pi.FaceAttrNameVec[i].empty())
						{ // trying to use named attribute to retrieve the value to store
							assert(vcg::tri::HasPerFaceAttribute(m,pi.FaceAttrNameVec[i]));
							if (!pi.FaceDescriptorVec[i].islist){
								switch (pi.FaceDescriptorVec[i].stotype1)
								{
								case ply::T_FLOAT  : tf=thff[i][fp]; fwrite(&tf, sizeof(float),1,fpout); break;
								case ply::T_DOUBLE : td=thdf[i][fp]; fwrite(&td, sizeof(double),1,fpout); break;
								case ply::T_INT    : ti=thif[i][fp]; fwrite(&ti, sizeof(int),1,fpout); break;
								case ply::T_SHORT  : ts=thsf[i][fp]; fwrite(&ts, sizeof(short),1,fpout); break;
								case ply::T_CHAR   : tc=thcf[i][fp]; fwrite(&tc, sizeof(char),1,fpout); break;
								case ply::T_UCHAR  : tu=thuf[i][fp]; fwrite(&tu,sizeof(unsigned char),1,fpout); break;
								default : assert(0);
								}
							}
							else {
								static const unsigned char psize = 3;
								switch (pi.FaceDescriptorVec[i].stotype1)
								{
								case ply::T_FLOAT  :
									fwrite(&psize, sizeof(unsigned char), 1,fpout);
									fwrite(&thp3ff[i][fp][0], sizeof(float), 1,fpout);
									fwrite(&thp3ff[i][fp][1], sizeof(float), 1,fpout);
									fwrite(&thp3ff[i][fp][2], sizeof(float), 1,fpout);
									break;
								case ply::T_DOUBLE :
									fwrite(&psize, sizeof(unsigned char), 1,fpout);
									fwrite(&thp3df[i][fp][0], sizeof(double), 1,fpout);
									fwrite(&thp3df[i][fp][1], sizeof(double), 1,fpout);
									fwrite(&thp3df[i][fp][2], sizeof(double), 1,fpout);
								default : assert(0);
								}
							}
						}
						else
						{
							switch (pi.FaceDescriptorVec[i].stotype1){
							case ply::T_FLOAT	 :		PlyConv(pi.FaceDescriptorVec[i].memtype1,  ((char *)fp)+pi.FaceDescriptorVec[i].offset1, tf );	fwrite(&tf, sizeof(float),1,fpout); break;
							case ply::T_DOUBLE :		PlyConv(pi.FaceDescriptorVec[i].memtype1,  ((char *)fp)+pi.FaceDescriptorVec[i].offset1, td );	fwrite(&td, sizeof(double),1,fpout); break;
							case ply::T_INT		 :		PlyConv(pi.FaceDescriptorVec[i].memtype1,  ((char *)fp)+pi.FaceDescriptorVec[i].offset1, ti );	fwrite(&ti, sizeof(int),1,fpout); break;
							case ply::T_SHORT	 :		PlyConv(pi.FaceDescriptorVec[i].memtype1,  ((char *)fp)+pi.FaceDescriptorVec[i].offset1, ts );	fwrite(&ts, sizeof(short),1,fpout); break;
							case ply::T_CHAR	 :		PlyConv(pi.FaceDescriptorVec[i].memtype1,  ((char *)fp)+pi.FaceDescriptorVec[i].offset1, tc );	fwrite(&tc, sizeof(char),1,fpout); break;
							case ply::T_UCHAR	 :		PlyConv(pi.FaceDescriptorVec[i].memtype1,  ((char *)fp)+pi.FaceDescriptorVec[i].offset1, tu );	fwrite(&tu, sizeof(unsigned char),1,fpout); break;
							default : assert(0);
							}
						}
					}
				}
				else	// ***** ASCII *****
				{
					fprintf(fpout,"%d " ,fp->VN());
					for(int k=0;k<fp->VN();++k)
						fprintf(fpout,"%d ",indices[fp->cV(k)]);

					if(HasPerFaceFlags(m)&&( pi.mask & Mask::IOM_FACEFLAGS ))
						fprintf(fpout,"%d ",fp->Flags());

					if( HasPerVertexTexCoord(m) && (pi.mask & Mask::IOM_WEDGTEXCOORD) ) // you can save VT as WT if you really want it...
					{
						fprintf(fpout,"%d ",fp->VN()*2);
						for(int k=0;k<fp->VN();++k)
							fprintf(fpout,"%.*g %.*g "
									,DGTFT
									,fp->V(k)->T().u()
									,fp->V(k)->T().v()
									);
					}
					else if( HasPerWedgeTexCoord(m) && (pi.mask & Mask::IOM_WEDGTEXCOORD)  )
					{
						fprintf(fpout,"%d ",fp->VN()*2);
						for(int k=0;k<fp->VN();++k)
							fprintf(fpout,"%f %f "
									,fp->WT(k).u()
									,fp->WT(k).v()
									);
					}

					if(saveTexIndexFlag)
					{
						fprintf(fpout,"%d ",fp->WT(0).n());
					}

					if( HasPerFaceColor(m) && (pi.mask & Mask::IOM_FACECOLOR)  )
					{
						fprintf(fpout, "%u %u %u %u ", fp->C()[0], fp->C()[1], fp->C()[2], fp->C()[3]);
					}
					else if( HasPerWedgeColor(m) && (pi.mask & Mask::IOM_WEDGCOLOR)  )
					{
						fprintf(fpout,"9 ");
						for(int z=0;z<3;++z)
							fprintf(fpout,"%g %g %g "
									,double(fp->WC(z)[0])/255
									,double(fp->WC(z)[1])/255
									,double(fp->WC(z)[2])/255
									);
					}

					if (HasPerFaceNormal(m) && (pi.mask & Mask::IOM_FACENORMAL))
						fprintf(fpout,"%.*g %.*g %.*g " ,DGT, ScalarType(fp->N()[0]),DGT,ScalarType(fp->N()[1]),DGT,ScalarType(fp->N()[2]));

					if( HasPerFaceQuality(m) && (pi.mask & Mask::IOM_FACEQUALITY) )
						fprintf(fpout,"%.*g ",DGTFQ,fp->Q());

					for(size_t i=0;i<pi.FaceDescriptorVec.size();i++)
					{
						float tf(0); double td(0); int ti;
						if(!pi.FaceAttrNameVec.empty() && !pi.FaceAttrNameVec[i].empty())
						{ // trying to use named attribute to retrieve the value to store
							assert(vcg::tri::HasPerFaceAttribute(m,pi.FaceAttrNameVec[i]));
							if(!pi.FaceDescriptorVec[i].islist) {
								switch (pi.FaceDescriptorVec[i].stotype1)
								{
								case ply::T_FLOAT  : tf=thff[i][fp]; fprintf(fpout,"%f ",tf); break;
								case ply::T_DOUBLE : td=thdf[i][fp]; fprintf(fpout,"%g ",td); break;
								case ply::T_INT	: ti=thif[i][fp]; fprintf(fpout,"%i ",ti); break;
								case ply::T_SHORT  : ti=thsf[i][fp]; fprintf(fpout,"%i ",ti); break;
								case ply::T_CHAR   : ti=thcf[i][fp]; fprintf(fpout,"%i ",ti); break;
								case ply::T_UCHAR  : ti=thuf[i][fp]; fprintf(fpout,"%i ",ti); break;
								default : assert(0);
								}
							}
							else {
								switch (pi.FaceDescriptorVec[i].stotype1)
								{
								case ply::T_FLOAT  : fprintf(fpout,"%d %f %f %f", 3, thp3ff[i][fp][0], thp3ff[i][fp][1], thp3ff[i][fp][2]); break;
								case ply::T_DOUBLE : fprintf(fpout,"%d %lf %lf %lf", 3, thp3df[i][fp][0], thp3df[i][fp][1], thp3df[i][fp][2]); break;
								default : assert(0);
								}
							}
						}
						else
						{
							switch (pi.FaceDescriptorVec[i].memtype1)
							{
							case  ply::T_FLOAT	:		tf=*( (float  *)		(((char *)fp)+pi.FaceDescriptorVec[i].offset1));	fprintf(fpout,"%g ",tf); break;
							case  ply::T_DOUBLE :		td=*( (double *)		(((char *)fp)+pi.FaceDescriptorVec[i].offset1));	fprintf(fpout,"%g ",tf); break;
							case  ply::T_INT		:		ti=*( (int	*)		(((char *)fp)+pi.FaceDescriptorVec[i].offset1));	fprintf(fpout,"%i ",ti); break;
							case  ply::T_SHORT	:		ti=*( (short  *)		(((char *)fp)+pi.FaceDescriptorVec[i].offset1));	fprintf(fpout,"%i ",ti); break;
							case  ply::T_CHAR		:		ti=*( (char   *)		(((char *)fp)+pi.FaceDescriptorVec[i].offset1));	fprintf(fpout,"%i ",ti); break;
							case  ply::T_UCHAR	:		ti=*( (unsigned char *) (((char *)fp)+pi.FaceDescriptorVec[i].offset1));	fprintf(fpout,"%i ",ti); break;
							default : assert(0);
							}
						}
					}

					fprintf(fpout,"\n");
				}
			}
		}
		assert(fcnt==m.fn);
		int eauxvv[2];
		if( pi.mask & Mask::IOM_EDGEINDEX )
		{
			int ecnt=0;
			for(EdgeIterator ei=m.edge.begin();ei!=m.edge.end();++ei)
			{
				if( ! ei->IsD() )
				{
					++ecnt;
					if(binary)
					{
						eauxvv[0]=indices[ei->cV(0)];
						eauxvv[1]=indices[ei->cV(1)];
						fwrite(eauxvv,sizeof(int),2,fpout);
					}
					else // ***** ASCII *****
						fprintf(fpout,"%d %d \n", indices[ei->cV(0)],	indices[ei->cV(1)]);
				}
			}
			assert(ecnt==m.en);
		}
		int result = 0;
		if (ferror(fpout)) result = ply::E_STREAMERROR;
		fclose(fpout);
		return result;
	}

	static const char *ErrorMsg(int error)
	{
		static std::vector<std::string> ply_error_msg;
		if(ply_error_msg.empty())
		{
			ply_error_msg.resize(PlyInfo::E_MAXPLYINFOERRORS );
			ply_error_msg[ply::E_NOERROR			]="No errors";
			ply_error_msg[ply::E_CANTOPEN		   ]="Can't open file";
			ply_error_msg[ply::E_NOTHEADER		  ]="Header not found";
			ply_error_msg[ply::E_UNESPECTEDEOF		]="Eof in header";
			ply_error_msg[ply::E_NOFORMAT		   ]="Format not found";
			ply_error_msg[ply::E_SYNTAX				]="Syntax error on header";
			ply_error_msg[ply::E_PROPOUTOFELEMENT   ]="Property without element";
			ply_error_msg[ply::E_BADTYPENAME		]="Bad type name";
			ply_error_msg[ply::E_ELEMNOTFOUND		]="Element not found";
			ply_error_msg[ply::E_PROPNOTFOUND		]="Property not found";
			ply_error_msg[ply::E_BADTYPE			]="Bad type on addtoread";
			ply_error_msg[ply::E_INCOMPATIBLETYPE   ]="Incompatible type";
			ply_error_msg[ply::E_BADCAST			]="Bad cast";

			ply_error_msg[ply::E_STREAMERROR		] = "Output Stream Error";

			ply_error_msg[PlyInfo::E_NO_VERTEX	  ]="No vertex field found";
			ply_error_msg[PlyInfo::E_NO_FACE		]="No face field found";
			ply_error_msg[PlyInfo::E_SHORTFILE	  ]="Unexpected EOF";
			ply_error_msg[PlyInfo::E_NO_3VERTINFACE ]="Face with more than 3 vertices";
			ply_error_msg[PlyInfo::E_BAD_VERT_INDEX ]="Bad vertex index in face";
			ply_error_msg[PlyInfo::E_NO_6TCOORD	 ]="Face with no 6 texture coordinates";
			ply_error_msg[PlyInfo::E_DIFFER_COLORS  ]="Number of color differ from vertices";
		}

		if(error>PlyInfo::E_MAXPLYINFOERRORS || error<0) return "Unknown error";
		else return ply_error_msg[error].c_str();
	};

	static int GetExportMaskCapability()
	{
		int capability = 0;
		capability |= vcg::tri::io::Mask::IOM_VERTCOORD	;
		capability |= vcg::tri::io::Mask::IOM_VERTFLAGS	;
		capability |= vcg::tri::io::Mask::IOM_VERTCOLOR	;
		capability |= vcg::tri::io::Mask::IOM_VERTQUALITY  ;
		capability |= vcg::tri::io::Mask::IOM_VERTNORMAL   ;
		capability |= vcg::tri::io::Mask::IOM_VERTRADIUS   ;
		capability |= vcg::tri::io::Mask::IOM_VERTTEXCOORD ;
		capability |= vcg::tri::io::Mask::IOM_FACEINDEX	;
		capability |= vcg::tri::io::Mask::IOM_FACEFLAGS	;
		capability |= vcg::tri::io::Mask::IOM_FACECOLOR	;
		capability |= vcg::tri::io::Mask::IOM_FACEQUALITY  ;
		// capability |= vcg::tri::io::Mask::IOM_FACENORMAL   ;
		capability |= vcg::tri::io::Mask::IOM_WEDGCOLOR	;
		capability |= vcg::tri::io::Mask::IOM_WEDGTEXCOORD ;
		capability |= vcg::tri::io::Mask::IOM_WEDGTEXMULTI ;
		capability |= vcg::tri::io::Mask::IOM_WEDGNORMAL   ;
		capability |= vcg::tri::io::Mask::IOM_CAMERA   ;
		//capability |= vcg::tri::io::Mask::IOM_BITPOLYGONAL;
		return capability;
	}


}; // end class



} // end namespace tri
} // end namespace io
} // end namespace vcg
//@}
#endif
