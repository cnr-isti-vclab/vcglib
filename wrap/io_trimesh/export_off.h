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
@name Save in OFF format
*/
//@{

#ifndef __VCGLIB_EXPORT_OFF
#define __VCGLIB_EXPORT_OFF

#include <stdio.h>
#include <wrap/io_trimesh/io_mask.h>
#include<wrap/io_trimesh/precision.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/polygon_support.h>


namespace vcg {
namespace tri {
namespace io {
template <class SaveMeshType>
class ExporterOFF
{

public:
  typedef typename SaveMeshType::VertexPointer VertexPointer;
  typedef typename SaveMeshType::ScalarType ScalarType;
  typedef typename SaveMeshType::VertexType VertexType;
  typedef typename SaveMeshType::FaceType FaceType;
  typedef typename SaveMeshType::FacePointer FacePointer;
  typedef typename SaveMeshType::VertexIterator VertexIterator;
  typedef typename SaveMeshType::FaceIterator FaceIterator;

  static int Save(SaveMeshType &m, const char * filename, int mask=0 )
  {
    vcg::face::Pos<FaceType> he;
    vcg::face::Pos<FaceType> hei;
    FILE * fpout = fopen(filename,"w");
    if(fpout==NULL)	return 1; // 1 is the error code for cant'open, see the ErrorMsg function


    if( tri::HasPerVertexColor(m)  && (mask & io::Mask::IOM_VERTNORMAL))
      fprintf(fpout,"N");
    if( tri::HasPerVertexColor(m)   && (mask & io::Mask::IOM_VERTCOLOR))
      fprintf(fpout,"C");
    if( tri::HasPerVertexTexCoord(m) && (mask & io::Mask::IOM_VERTTEXCOORD))
      fprintf(fpout,"ST");
    fprintf(fpout,"OFF\n");

    int polynumber;
    if (mask &io::Mask::IOM_BITPOLYGONAL)
      polynumber = tri::Clean<SaveMeshType>::CountBitLargePolygons(m);
    else
      polynumber = m.fn;

    fprintf(fpout,"%d %d 0\n", m.vn, polynumber); // note that as edge number we simply write zero

    //vertices
    const int DGT = vcg::tri::io::Precision<ScalarType>::digits();
    for(auto vp=m.vert.begin();vp!=m.vert.end();++vp)
    {
      if( ! vp->IsD() )
      {	// ***** ASCII *****

        fprintf(fpout,"%.*g %.*g %.*g " ,DGT,vp->P()[0],DGT,vp->P()[1],DGT,vp->P()[2]);
        if( tri::HasPerVertexColor(m)  && (mask & io::Mask::IOM_VERTCOLOR) )
          fprintf(fpout,"%d %d %d %d ",vp->C()[0],vp->C()[1],vp->C()[2],vp->C()[3] );

        if( tri::HasPerVertexNormal(m)  && (mask & io::Mask::IOM_VERTNORMAL) )
          fprintf(fpout,"%g %g %g ", double(vp->N()[0]),double(vp->N()[1]),double(vp->N()[2]));

        if( tri::HasPerVertexTexCoord(m)  && (mask & io::Mask::IOM_VERTTEXCOORD) )
          fprintf(fpout,"%g %g ",vp->T().u(),vp->T().v());

        fprintf(fpout,"\n");
      }
    }


    if (mask &io::Mask::IOM_BITPOLYGONAL) {
      tri::RequireFFAdjacency(m);
      std::vector<VertexPointer> polygon;
      tri::UpdateFlags<SaveMeshType>::FaceClearV(m);
      for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi) if (!fi->IsD()) if (!fi->IsV()) {
        vcg::tri::PolygonSupport<SaveMeshType,SaveMeshType>::ExtractPolygon(&*fi,polygon);
        if(!polygon.empty())
        {
          fprintf(fpout,"%d ", int(polygon.size()) );
          for (size_t i=0; i<polygon.size(); i++) fprintf(fpout,"%d ", polygon[i]->Flags() );
          fprintf(fpout,"\n");
        }
      }
    }
    else {
      for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
      {
        if( ! fi->IsD() )
        {
          fprintf(fpout,"%i ",fi->VN());
          for(int i=0;i<fi->VN();++i)
              fprintf(fpout,"%lu ",tri::Index(m,fi->V(i)));
          if( tri::HasPerFaceColor(m)  && (mask & io::Mask::IOM_FACECOLOR) )
            fprintf(fpout,"%i %i %i", fi->C()[0],fi->C()[1],fi->C()[2] );
          fprintf(fpout,"\n");
        }
      }
    }

	int result = 0;
	if (ferror(fpout)) result = 2;
	fclose(fpout);
	return result;
  }

  static const char *ErrorMsg(int error)
  {
	  static std::vector<std::string> off_error_msg;
	  if (off_error_msg.empty())
	  {
		  off_error_msg.resize(3);
		  off_error_msg[0] = "No errors";
		  off_error_msg[1] = "Can't open file";
		  off_error_msg[2] = "Output Stream error";
	  }

	  if (error>2 || error<0) return "Unknown error";
	  else return off_error_msg[error].c_str();
  }
  /*
            returns mask of capability one define with what are the saveable information of the format.
        */
  static int GetExportMaskCapability()
  {
    int capability = 0;
    capability |= vcg::tri::io::Mask::IOM_VERTCOORD;
    capability |= vcg::tri::io::Mask::IOM_VERTCOLOR;
    capability |= vcg::tri::io::Mask::IOM_VERTTEXCOORD;
    capability |= vcg::tri::io::Mask::IOM_FACEINDEX;
    capability |= vcg::tri::io::Mask::IOM_FACECOLOR;
    capability |= vcg::tri::io::Mask::IOM_BITPOLYGONAL;
    return capability;
  }

}; // end class
} // end namespace tri
} // end namespace io
} // end namespace vcg
//@}
#endif
