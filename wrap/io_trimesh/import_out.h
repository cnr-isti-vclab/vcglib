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
#ifndef __VCGLIB_IMPORTERBUNDLER
#define __VCGLIB_IMPORTERBUNDLER

#include <stddef.h>
#include <stdio.h>
#include <vcg/complex/complex.h>
//#include <vcg/complex/allocate.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <wrap/callback.h>
#include <wrap/io_trimesh/io_mask.h>
#include <QImageReader>
#include <QFileInfo>
#include "exif.h" //external easyexif lib

namespace vcg {
namespace tri {
namespace io {

struct Correspondence{
	Correspondence(unsigned int id_img_,unsigned int key_,float x_,float y_):id_img(id_img_),key(key_),x(x_),y(y_){}
	unsigned int id_img,key;
	float x;
	float y;
};

typedef std::vector<Correspondence> CorrVec;

/**
This class encapsulate a filter for opening bundler file
*/
template <class OpenMeshType>
class ImporterOUT
{
public:
	
	typedef typename OpenMeshType::VertexPointer VertexPointer;
	typedef typename OpenMeshType::ScalarType ScalarType;
	typedef typename OpenMeshType::CoordType CoordType;
	typedef typename OpenMeshType::VertexType VertexType;
	typedef typename OpenMeshType::FaceType FaceType;
	typedef typename OpenMeshType::VertexIterator VertexIterator;
	typedef typename OpenMeshType::FaceIterator FaceIterator;
	typedef typename OpenMeshType::EdgeIterator EdgeIterator;
	
	static void readline(FILE *fp, char *line, int max=100){
		fgets ( line, max, fp);
	}
	
	static bool ReadHeader(FILE *fp,unsigned int &num_cams, unsigned int &num_points){
		char line[100];
		readline(fp, line);
		if( line[0]=='\0' ) return false;
		line[18]='\0';
		if(0!=strcmp("# Bundle file v0.3", line))  return false;
		readline(fp, line);
		if(line[0]=='\0') return false;
		if (sscanf(line, "%d %d", &num_cams, &num_points) != 2) return false;
		return true;
	}
	
	static int Open(
			OpenMeshType &m, 
			std::vector<Shot<ScalarType> >  & shots,
			std::vector<std::string > & image_filenames,
			const char * filename,
			const char * filename_images, 
			CallBackPos *cb=0)
	{
		unsigned int   num_cams,num_points;
		typedef typename vcg::Matrix44<ScalarType> Matrix44x;
		typedef typename vcg::Matrix33<ScalarType> Matrix33x;
		FILE *fp = fopen(filename,"r");
		if(!fp) return false;
		ReadHeader(fp, num_cams,  num_points);
		char line[100];
		if(cb) cb(0,"Reading images");
		ReadImagesFilenames(filename_images, image_filenames);
		const QString path_im = QFileInfo(filename_images).absolutePath()+QString("/");
		
		if(cb) cb(50,"Reading cameras");
		shots.resize(num_cams);
		for(uint i = 0; i < num_cams;++i)
		{
			float f, k1, k2;
			float R[16]={0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,1};
			vcg::Point3f t;
			
			readline(fp, line); if(line[0]=='\0') return false; if (sscanf(line, "%f %f %f", &f, &k1, &k2) != 3) return false;
			
			readline(fp, line); if(line[0]=='\0') return false; if (sscanf(line, "%f %f %f", &(R[0]), &(R[1]), &(R[2])) != 3) return false;  R[3] = 0;
			readline(fp, line); if(line[0]=='\0') return false; if (sscanf(line, "%f %f %f", &(R[4]), &(R[5]), &(R[6])) != 3) return false;  R[7] = 0;
			readline(fp, line); if(line[0]=='\0') return false; if (sscanf(line, "%f %f %f", &(R[8]), &(R[9]), &(R[10])) != 3) return false; R[11] = 0;
			readline(fp, line); if(line[0]=='\0') return false; if (sscanf(line, "%f %f %f", &(t[0]), &(t[1]), &(t[2])) != 3) return false;
			
			Matrix44x mat = Matrix44x::Construct(Matrix44f(R));
			
			Matrix33x Rt = Matrix33x( Matrix44x(mat), 3);
			Rt.Transpose();
			
			CoordType pos = Rt * CoordType(t[0], t[1], t[2]);
			
			shots[i].Extrinsics.SetTra(CoordType(-pos[0],-pos[1],-pos[2]));
			shots[i].Extrinsics.SetRot(mat);
			
			shots[i].Intrinsics.FocalMm    = f;
			shots[i].Intrinsics.k[0] = 0.0;//k1; To be uncommented when distortion is taken into account reliably
			shots[i].Intrinsics.k[1] = 0.0;//k2;
			shots[i].Intrinsics.PixelSizeMm = vcg::Point2<ScalarType>(1,1);
			QSize size;
			QImageReader sizeImg(QString::fromStdString(image_filenames[i]));
			if(sizeImg.size()==QSize(-1,-1))
			{
				QImageReader sizeImg(QString::fromStdString(qPrintable(path_im)+image_filenames[i]));
				size=sizeImg.size();
			}
			else
				size=sizeImg.size();
			shots[i].Intrinsics.ViewportPx = vcg::Point2i(size.width(),size.height());
			shots[i].Intrinsics.CenterPx[0] = (int)((double)shots[i].Intrinsics.ViewportPx[0]/2.0f);
			shots[i].Intrinsics.CenterPx[1] = (int)((double)shots[i].Intrinsics.ViewportPx[1]/2.0f);
			//AddIntrinsics(shots[i], std::string(filename_images_path).append(image_filenames[i]).c_str());
		}
		
		// load all correspondences
		typename OpenMeshType::template PerVertexAttributeHandle<CorrVec> ch = vcg::tri::Allocator<OpenMeshType>::template GetPerVertexAttribute<CorrVec>(m,"correspondences");
		
		typename OpenMeshType::VertexIterator vi = vcg::tri::Allocator<OpenMeshType>::AddVertices(m,num_points);
		for(uint i = 0; i < num_points;++i,++vi){
			double x,y,z;
			unsigned int r,g,b,i_cam, key_sift,n_corr;
			if (fscanf(fp,"%lf %lf %lf ",&x,&y,&z) != 3) return false;
			(*vi).P() = vcg::Point3<typename OpenMeshType::ScalarType>(x,y,z);
			if (fscanf(fp,"%d %d %d ",&r,&g,&b) != 3) return false;
			(*vi).C() = vcg::Color4b(r,g,b,255);
			
			if (fscanf(fp,"%d ",&n_corr) != 1) return false;
			for(uint j = 0; j < n_corr; ++j){
				if (fscanf(fp,"%d %d %lf %lf ",&i_cam,&key_sift,&x,&y) != 4) return false;
				Correspondence corr(i_cam,key_sift,x,y);
				ch[i].push_back(corr);
			}
		}
		vcg::tri::UpdateBounding<OpenMeshType>::Box(m);
		fclose(fp);
		
		return (shots.size() == 0);
	}
	
	
	static bool ReadImagesFilenames(const char *  filename,std::vector<std::string> &image_filenames)
	{
		FILE * fp = fopen(filename,"r");
		if (!fp) return false;
		else
		{
			char line[1000], name[1000];
			while(!feof(fp)){
				readline(fp, line, 1000);
				if(line[0] == '\0') continue; //ignore empty lines (in theory, might happen only at end of file)
				if (sscanf(line, "%s", name) != 1) return false;
				std::string n(name);
				image_filenames.push_back(n);
			}
		}
		fclose(fp);
		return true;
	}
	
	static bool  AddIntrinsics(vcg::Shotf &shot, const char * image_file)
	{
		// Read the JPEG file into a buffer
		FILE *fp = fopen(qUtf8Printable(image_file), "rb");
		if (!fp) {
			std::cerr << "Exif Parsing: Unable to open file:\n\"%1\"\n\nError details: file %1 is not readable.";
			return false;
		}
		fseek(fp, 0, SEEK_END);
		unsigned long fsize = ftell(fp);
		rewind(fp);
		unsigned char *buf = new unsigned char[fsize];
		if (fread(buf, 1, fsize, fp) != fsize) {
			std::cerr << "Exif Parsing: Unable to read the content of the opened file:\n\"%1\"\n\nError details: file %1 is not readable.";
			delete[] buf;
			fclose(fp);
			return false;
		}
		fclose(fp);
		
		// Parse EXIF
		easyexif::EXIFInfo ImageInfo;
		int ret = ImageInfo.parseFrom(buf, fsize);
		delete[] buf;
		if (ret == 0) {
			std::cerr << "Warning unable to parse exif for file  %s" << qPrintable(image_file);
			return false;
		}
		
		shot.Intrinsics.ViewportPx = vcg::Point2i(ImageInfo.ImageWidth, ImageInfo.ImageHeight);
		shot.Intrinsics.CenterPx   = vcg::Point2f(float(ImageInfo.ImageWidth/2.0), float(ImageInfo.ImageHeight/2.0));
		
		return true;
	}
}; // end class



} // end namespace tri
} // end namespace io
} // end namespace vcg

#endif

