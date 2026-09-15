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

#include<wrap/callback.h>
#include<wrap/ply/plylib.h>
#include<wrap/io_trimesh/io_mask.h>
#include<wrap/io_trimesh/io_ply.h>
#include<wrap/io_trimesh/precision.h>
#include<vcg/container/simple_temporary_data.h>
#include <vcg/complex/base.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/polygon_support.h>
#include <vcg/complex/algorithms/update/topology.h>


#include <stdio.h>

namespace vcg {
namespace tri {
namespace io {



template <class SaveMeshType>
class ExporterPLY
{
	// It takes care of converting from one type to another.
	// It is used in saveply to match types between stotype and memtype (e.g. the type in the file and the type in memory).
	// For example if there is an int in memory and I want to save a float
	// src will actually be a pointer to int whose value must
	// be converted to the desired return type (stotype)

	template <class StoType>
	static void PlyConv(int mem_type, void *src, StoType &dest)
	{
		switch (mem_type){
		case ply::T_FLOAT	:		dest = (StoType) (*  ((float  *) src)); break;
		case ply::T_DOUBLE	:		dest = (StoType) (*  ((double *) src)); break;
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
	typedef typename SaveMeshType::FacePointer FacePointer;
	typedef typename SaveMeshType::ConstVertexIterator VertexIterator;
	typedef typename SaveMeshType::FaceIterator FaceIterator;
	typedef typename SaveMeshType::ConstEdgeIterator EdgeIterator;
	typedef typename vcg::Shot<ScalarType>::ScalarType ShotScalarType;

	// Preserve the historical const API for ordinary triangular export.
	// Polygon reconstruction needs mutable FF adjacency and visited scratch flags.
	static int Save(const SaveMeshType &m, const char *filename, bool binary=true)
	{
		return Save(const_cast<SaveMeshType &>(m),filename,binary);
	}

	static int Save(const SaveMeshType &m, const char *filename, int savemask, bool binary=true, CallBackPos *cb=0)
	{
		if(savemask & Mask::IOM_BITPOLYGONAL) return ply::E_STREAMERROR;
		return Save(const_cast<SaveMeshType &>(m),filename,savemask,binary,cb);
	}

	static int Save(const SaveMeshType &m, const char *filename, bool binary, const PlyInfo &pi, CallBackPos *cb=0)
	{
		if(pi.mask & Mask::IOM_BITPOLYGONAL) return ply::E_STREAMERROR;
		return Save(const_cast<SaveMeshType &>(m),filename,binary,pi,cb);
	}

	static int Save(SaveMeshType &m, const char * filename, bool binary=true)
	{
		PlyInfo pi;
		return Save(m,filename,binary,pi);
	}

	static int Save(SaveMeshType &m,  const char * filename, int savemask, bool binary = true, CallBackPos *cb=0 )
	{
		PlyInfo pi;
		pi.mask=savemask;
		return Save(m,filename,binary,pi,cb);
	}

	static int Save(SaveMeshType &m,  const char * filename, bool binary, const PlyInfo &pi, CallBackPos *cb=0)	// V1.0
	{
		const bool savePolygons = (pi.mask & Mask::IOM_BITPOLYGONAL) != 0;
		if(savePolygons)
		{
			// Faux edges encode polygon interiors. Recompute FF adjacency here so
			// traversal never relies on stale topology supplied by the caller.
			tri::RequireFFAdjacency(m);
			tri::UpdateTopology<SaveMeshType>::FaceFace(m);
		}

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

		// PLY declares the number of face records before their data, hence the
		// inexpensive first traversal when triangles must be grouped as polygons.
		const int faceCount = savePolygons
			? tri::Clean<SaveMeshType>::CountBitLargePolygons(m)
			: m.fn;
		fprintf(fpout,
				"element face %d\n"
				"property list uchar int vertex_indices\n",
				faceCount );

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
		{
			fprintf(
						fpout,
						"element edge %d\n" "property int vertex1\n""property int vertex2\n",m.en);
			if(HasPerEdgeColor(m) && (pi.mask & Mask::IOM_EDGECOLOR))
				fprintf(fpout,
						"property uchar red\n"
						"property uchar green\n"
						"property uchar blue\n"
						"property uchar alpha\n");
		}
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
        
        PlyAttributeHelper<SaveMeshType> PAH(pi,m);

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
								case ply::T_FLOAT  : tf=PAH.tchfv[i][vp]; fwrite(&tf, sizeof(float),1,fpout); break;
								case ply::T_DOUBLE : td=PAH.tchdv[i][vp]; fwrite(&td, sizeof(double),1,fpout); break;
								case ply::T_INT    : ti=PAH.tchiv[i][vp]; fwrite(&ti, sizeof(int),1,fpout); break;
								case ply::T_SHORT  : ts=PAH.tchsv[i][vp]; fwrite(&ts, sizeof(short),1,fpout); break;
								case ply::T_CHAR   : tc=PAH.tchcv[i][vp]; fwrite(&tc, sizeof(char),1,fpout); break;
								case ply::T_UCHAR  : tu=PAH.tchuv[i][vp]; fwrite(&tu,sizeof(unsigned char),1,fpout); break;
								default : assert(0);
								}
							}
							else { //it is a Poin3f or a Point3d attribute. Saving it as a list
								static const unsigned char psize = 3;
								switch (pi.VertDescriptorVec[i].stotype1)
								{
								case ply::T_FLOAT  :
									fwrite(&psize, sizeof(unsigned char), 1,fpout);
									fwrite(&PAH.tchp3fv[i][vp][0], sizeof(float), 1,fpout);
									fwrite(&PAH.tchp3fv[i][vp][1], sizeof(float), 1,fpout);
									fwrite(&PAH.tchp3fv[i][vp][2], sizeof(float), 1,fpout);
									break;
									//fprintf(fpout,"%d %f %f %f", 3, thp3fv[i][vp][0], thp3fv[i][vp][1], thp3fv[i][vp][2]); break;
								case ply::T_DOUBLE :
									fwrite(&psize, sizeof(unsigned char), 1,fpout);
									fwrite(&PAH.tchp3dv[i][vp][0], sizeof(double), 1,fpout);
									fwrite(&PAH.tchp3dv[i][vp][1], sizeof(double), 1,fpout);
									fwrite(&PAH.tchp3dv[i][vp][2], sizeof(double), 1,fpout);
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
						fprintf(fpout,"%.*g %.*g",DGTVT,vp->T().u(),DGTVT,vp->T().v());

					for(size_t i=0;i<pi.VertDescriptorVec.size();i++)
					{
						float tf(0); double td(0); int ti;
						if(!pi.VertAttrNameVec.empty() && !pi.VertAttrNameVec[i].empty())
						{ // trying to use named attribute to retrieve the value to store
							assert(vcg::tri::HasPerVertexAttribute(m,pi.VertAttrNameVec[i]));
							if (!pi.VertDescriptorVec[i].islist){
								switch (pi.VertDescriptorVec[i].stotype1)
								{
								case ply::T_FLOAT  : tf=PAH.tchfv[i][vp]; fprintf(fpout,"%f ",tf); break;
								case ply::T_DOUBLE : td=PAH.tchdv[i][vp]; fprintf(fpout,"%lf ",td); break;
								case ply::T_INT    : ti=PAH.tchiv[i][vp]; fprintf(fpout,"%i ",ti); break;
								case ply::T_SHORT  : ti=PAH.tchsv[i][vp]; fprintf(fpout,"%i ",ti); break;
								case ply::T_CHAR   : ti=PAH.tchcv[i][vp]; fprintf(fpout,"%i ",ti); break;
								case ply::T_UCHAR  : ti=PAH.tchuv[i][vp]; fprintf(fpout,"%i ",ti); break;
								default : assert(0);
								}
							}
							else { //it is a Poin3f or a Point3d attribute. Saving it as a list
								switch (pi.VertDescriptorVec[i].stotype1)
								{
								case ply::T_FLOAT  : fprintf(fpout,"%d %f %f %f", 3, PAH.tchp3fv[i][vp][0], PAH.tchp3fv[i][vp][1], PAH.tchp3fv[i][vp][2]); break;
								case ply::T_DOUBLE : fprintf(fpout,"%d %lf %lf %lf", 3, PAH.tchp3dv[i][vp][0], PAH.tchp3dv[i][vp][1], PAH.tchp3dv[i][vp][2]); break;
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

		FacePointer fp;
		FaceIterator fi;
		int fcnt=0;
		typedef std::pair<FacePointer, int> BoundaryCorner;
		std::vector<typename SaveMeshType::VertexPointer> polygonVertices;
		std::vector<FacePointer> polygonFaces;
		std::vector<BoundaryCorner> boundaryCorners;
		std::vector<int> vertexIndices;
		// ExtractPolygon uses VISITED as scratch state. The corner pairs retain the
		// source triangle and local corner needed for boundary wedge attributes.
		if(savePolygons)
			tri::UpdateFlags<SaveMeshType>::FaceClearV(m);
		for(j=0,fi=m.face.begin();fi!=m.face.end();++fi)
		{
			//((m.vn+m.fn) != 0) all vertices and faces have been marked as deleted but the are still in the vert/face vectors
			if(cb && ((j%1000)==0) && ((m.vn+m.fn) != 0))
				(*cb)( 100*(m.vn+j)/(m.vn+m.fn), "Saving Faces");

			fp=&(*fi);
			if( !fp->IsD() && !(savePolygons && fp->IsV()) )
			{ fcnt++;
				boundaryCorners.clear();
				polygonFaces.clear();
				if(savePolygons)
				{
					vcg::tri::PolygonSupport<SaveMeshType,SaveMeshType>::ExtractPolygon(
						fp, polygonVertices, polygonFaces, boundaryCorners);
					// Multi-triangle polygons are traversed clockwise; restore the
					// source winding. Plain triangles already have the correct order.
					if(boundaryCorners.size()>3)
						std::reverse(boundaryCorners.begin(), boundaryCorners.end());
				}
				else
				{
					polygonFaces.push_back(fp);
					for(int k=0;k<fp->VN();++k)
						boundaryCorners.emplace_back(fp,k);
				}
				const size_t cornerCount = boundaryCorners.size();
				// PLY list counts are declared as uchar; account for the number of
				// index/UV/color scalars written for each polygon corner.
				const size_t listScalarsPerCorner = (HasPerWedgeColor(m) && (pi.mask & Mask::IOM_WEDGCOLOR)) ? 3
					: ((pi.mask & Mask::IOM_WEDGTEXCOORD) ? 2 : 1);
				if(cornerCount<3 || cornerCount>255/listScalarsPerCorner)
				{
					if(savePolygons) tri::UpdateFlags<SaveMeshType>::FaceClearV(m);
					fclose(fpout);
					return ply::E_STREAMERROR;
				}
				// Progress is measured in source triangles. Scalar face attributes
				// below necessarily come from fp, the polygon's representative face.
				j += int(polygonFaces.size());
				if(binary)
				{
					const unsigned char listSize = static_cast<unsigned char>(cornerCount);
					vertexIndices.resize(cornerCount);
					for(size_t k=0;k<cornerCount;++k)
						vertexIndices[k]=indices[boundaryCorners[k].first->cV(boundaryCorners[k].second)];
					fwrite(&listSize,sizeof(char),1,fpout);
					fwrite(vertexIndices.data(),sizeof(int),cornerCount,fpout);

					if(HasPerFaceFlags(m)&&( pi.mask & Mask::IOM_FACEFLAGS) ){
						auto fl = fp->Flags();
						// VISITED is exporter scratch state and faux bits describe the
						// discarded triangulation, not the emitted polygon boundary.
						if(savePolygons)
							fl &= ~(FaceType::VISITED | FaceType::FAUX012);
						fwrite(&fl,sizeof(int),1,fpout);
					}

					if( HasPerVertexTexCoord(m) && (!HasPerWedgeTexCoord(m)) && (pi.mask & Mask::IOM_WEDGTEXCOORD) )  // Note that you can save VT as WT if you really want it...
					{
						const unsigned char listSize = static_cast<unsigned char>(cornerCount*2);
						fwrite(&listSize,sizeof(char),1,fpout);
						for(const BoundaryCorner &corner : boundaryCorners)
						{
							typename FaceType::TexCoordType::ScalarType t = corner.first->V(corner.second)->T().u();
							fwrite(&t,sizeof(t),1,fpout);
							t = corner.first->V(corner.second)->T().v();
							fwrite(&t,sizeof(t),1,fpout);
						}
					}
					else if( HasPerWedgeTexCoord(m) && (pi.mask & Mask::IOM_WEDGTEXCOORD)  )
					{
						const unsigned char listSize = static_cast<unsigned char>(cornerCount*2);
						fwrite(&listSize,sizeof(char),1,fpout);
						for(const BoundaryCorner &corner : boundaryCorners)
						{
							typename FaceType::TexCoordType::ScalarType t = corner.first->WT(corner.second).u();
							fwrite(&t,sizeof(t),1,fpout);
							t = corner.first->WT(corner.second).v();
							fwrite(&t,sizeof(t),1,fpout);
						}
					}

					if(saveTexIndexFlag)
					{
						const BoundaryCorner &corner = boundaryCorners.front();
						int t = corner.first->WT(corner.second).n();
						fwrite(&t,sizeof(int),1,fpout);
					}

					if( HasPerFaceColor(m) && (pi.mask & Mask::IOM_FACECOLOR) )
						fwrite(&( fp->C() ),sizeof(char),4,fpout);


					if( HasPerWedgeColor(m) && (pi.mask & Mask::IOM_WEDGCOLOR)  )
					{
						const unsigned char listSize = static_cast<unsigned char>(cornerCount*3);
						fwrite(&listSize,sizeof(char),1,fpout);
						float t[3];
						for(const BoundaryCorner &corner : boundaryCorners)
						{
							t[0] = float(corner.first->WC(corner.second)[0])/255;
							t[1] = float(corner.first->WC(corner.second)[1])/255;
							t[2] = float(corner.first->WC(corner.second)[2])/255;
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
								case ply::T_FLOAT  : tf=PAH.tchff[i][fp]; fwrite(&tf, sizeof(float),1,fpout); break;
								case ply::T_DOUBLE : td=PAH.tchdf[i][fp]; fwrite(&td, sizeof(double),1,fpout); break;
								case ply::T_INT    : ti=PAH.tchif[i][fp]; fwrite(&ti, sizeof(int),1,fpout); break;
								case ply::T_SHORT  : ts=PAH.tchsf[i][fp]; fwrite(&ts, sizeof(short),1,fpout); break;
								case ply::T_CHAR   : tc=PAH.tchcf[i][fp]; fwrite(&tc, sizeof(char),1,fpout); break;
								case ply::T_UCHAR  : tu=PAH.tchuf[i][fp]; fwrite(&tu,sizeof(unsigned char),1,fpout); break;
								default : assert(0);
								}
							}
							else {
								static const unsigned char psize = 3;
								switch (pi.FaceDescriptorVec[i].stotype1)
								{
								case ply::T_FLOAT  :
									fwrite(&psize, sizeof(unsigned char), 1,fpout);
									fwrite(&PAH.tchp3ff[i][fp][0], sizeof(float), 1,fpout);
									fwrite(&PAH.tchp3ff[i][fp][1], sizeof(float), 1,fpout);
									fwrite(&PAH.tchp3ff[i][fp][2], sizeof(float), 1,fpout);
									break;
								case ply::T_DOUBLE :
									fwrite(&psize, sizeof(unsigned char), 1,fpout);
									fwrite(&PAH.tchp3df[i][fp][0], sizeof(double), 1,fpout);
									fwrite(&PAH.tchp3df[i][fp][1], sizeof(double), 1,fpout);
									fwrite(&PAH.tchp3df[i][fp][2], sizeof(double), 1,fpout);
									break;
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
					fprintf(fpout,"%d " ,int(cornerCount));
					for(const BoundaryCorner &corner : boundaryCorners)
						fprintf(fpout,"%d ",indices[corner.first->cV(corner.second)]);

					if(HasPerFaceFlags(m)&&( pi.mask & Mask::IOM_FACEFLAGS ))
					{
						int flags = fp->Flags();
						// Do not leak traversal or internal-triangulation flags into PLY.
						if(savePolygons)
							flags &= ~(FaceType::VISITED | FaceType::FAUX012);
						fprintf(fpout,"%d ",flags);
					}

					// Match the binary path: genuine wedge UVs take precedence over
					// vertex UVs because only they can preserve texture seams.
					if( HasPerVertexTexCoord(m) && !HasPerWedgeTexCoord(m) && (pi.mask & Mask::IOM_WEDGTEXCOORD) ) // you can save VT as WT if you really want it...
					{
						fprintf(fpout,"%d ",int(cornerCount*2));
						for(const BoundaryCorner &corner : boundaryCorners)
							fprintf(fpout,"%.*g %.*g "
									,DGTFT,corner.first->V(corner.second)->T().u()
									,DGTFT,corner.first->V(corner.second)->T().v()
									);
					}
					else if( HasPerWedgeTexCoord(m) && (pi.mask & Mask::IOM_WEDGTEXCOORD)  )
					{
						fprintf(fpout,"%d ",int(cornerCount*2));
						for(const BoundaryCorner &corner : boundaryCorners)
							fprintf(fpout,"%f %f "
									,corner.first->WT(corner.second).u()
									,corner.first->WT(corner.second).v()
									);
					}

					if(saveTexIndexFlag)
					{
						const BoundaryCorner &corner = boundaryCorners.front();
						fprintf(fpout,"%d ",corner.first->WT(corner.second).n());
					}

					// Face and wedge colors are separate PLY properties and may coexist.
					if( HasPerFaceColor(m) && (pi.mask & Mask::IOM_FACECOLOR)  )
					{
						fprintf(fpout, "%u %u %u %u ", fp->C()[0], fp->C()[1], fp->C()[2], fp->C()[3]);
					}
					if( HasPerWedgeColor(m) && (pi.mask & Mask::IOM_WEDGCOLOR)  )
					{
						fprintf(fpout,"%d ",int(cornerCount*3));
						for(const BoundaryCorner &corner : boundaryCorners)
							fprintf(fpout,"%g %g %g "
									,double(corner.first->WC(corner.second)[0])/255
									,double(corner.first->WC(corner.second)[1])/255
									,double(corner.first->WC(corner.second)[2])/255
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
								case ply::T_FLOAT  : tf=PAH.tchff[i][fp]; fprintf(fpout,"%f ",tf); break;
								case ply::T_DOUBLE : td=PAH.tchdf[i][fp]; fprintf(fpout,"%g ",td); break;
								case ply::T_INT    : ti=PAH.tchif[i][fp]; fprintf(fpout,"%i ",ti); break;
								case ply::T_SHORT  : ti=PAH.tchsf[i][fp]; fprintf(fpout,"%i ",ti); break;
								case ply::T_CHAR   : ti=PAH.tchcf[i][fp]; fprintf(fpout,"%i ",ti); break;
								case ply::T_UCHAR  : ti=PAH.tchuf[i][fp]; fprintf(fpout,"%i ",ti); break;
								default : assert(0);
								}
							}
							else {
								switch (pi.FaceDescriptorVec[i].stotype1)
								{
								case ply::T_FLOAT  : fprintf(fpout,"%d %f %f %f", 3, PAH.tchp3ff[i][fp][0], PAH.tchp3ff[i][fp][1], PAH.tchp3ff[i][fp][2]); break;
								case ply::T_DOUBLE : fprintf(fpout,"%d %lf %lf %lf", 3, PAH.tchp3df[i][fp][0], PAH.tchp3df[i][fp][1], PAH.tchp3df[i][fp][2]); break;
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
		// Restore the only face state used as exporter scratch storage.
		if(savePolygons)
			tri::UpdateFlags<SaveMeshType>::FaceClearV(m);
		assert(fcnt==faceCount);
		(void)fcnt;
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
						if(HasPerEdgeColor(m) && (pi.mask & Mask::IOM_EDGECOLOR))
						{
							const vcg::Color4b edgeColor = ei->cC();
							fwrite(&edgeColor,sizeof(char),4,fpout);
						}
					}
					else // ***** ASCII *****
					{
						fprintf(fpout,"%d %d ", indices[ei->cV(0)], indices[ei->cV(1)]);
						if(HasPerEdgeColor(m) && (pi.mask & Mask::IOM_EDGECOLOR))
						{
							const vcg::Color4b edgeColor = ei->cC();
							fprintf(fpout,"%d %d %d %d ", edgeColor[0], edgeColor[1], edgeColor[2], edgeColor[3]);
						}
						fprintf(fpout,"\n");
					}
				}
			}
			assert(ecnt==m.en);
			(void)ecnt;
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
			ply_error_msg[PlyInfo::E_NO_6TCOORD	 ]="Texture coordinate count does not match face corners";
			ply_error_msg[PlyInfo::E_DIFFER_COLORS  ]="Wedge color count does not match face corners";
			ply_error_msg[PlyInfo::E_INVALID_POLYGON ]="Face is not a valid simple planar polygon";
		}

		if(error>=PlyInfo::E_MAXPLYINFOERRORS || error<0) return "Unknown error";
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
		capability |= vcg::tri::io::Mask::IOM_EDGEINDEX    ;
		capability |= vcg::tri::io::Mask::IOM_EDGECOLOR    ;
		capability |= vcg::tri::io::Mask::IOM_FACEINDEX	;
		capability |= vcg::tri::io::Mask::IOM_FACEFLAGS	;
		capability |= vcg::tri::io::Mask::IOM_FACECOLOR	;
		capability |= vcg::tri::io::Mask::IOM_FACEQUALITY  ;
		// Face normals have explicit nx/ny/nz PLY properties below.
		capability |= vcg::tri::io::Mask::IOM_FACENORMAL   ;
		capability |= vcg::tri::io::Mask::IOM_WEDGCOLOR	;
		capability |= vcg::tri::io::Mask::IOM_WEDGTEXCOORD ;
		capability |= vcg::tri::io::Mask::IOM_WEDGTEXMULTI ;
		// PLY has no standard per-corner normal property.
		capability |= vcg::tri::io::Mask::IOM_CAMERA   ;
		capability |= vcg::tri::io::Mask::IOM_BITPOLYGONAL;
		return capability;
	}


}; // end class



} // end namespace tri
} // end namespace io
} // end namespace vcg
//@}
#endif
