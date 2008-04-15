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
Revision 1.13  2007/11/06 10:58:25  cignoni
Changed the return value to the standard 0 in case of success and notzero for failures

Revision 1.12  2007/03/12 16:40:17  tarini
Texture coord name change!  "TCoord" and "Texture" are BAD. "TexCoord" is GOOD.

Revision 1.11  2006/12/07 00:37:58  cignoni
Corrected bug in the management of deleted vertices

****************************************************************************/

/**
@name Save in OFF format
*/
//@{

#ifndef __VCGLIB_EXPORT_OFF
#define __VCGLIB_EXPORT_OFF

#include <stdio.h>
#include <wrap/io_trimesh/io_mask.h>


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



					if( m.HasPerVertexNormal()  && (mask & io::Mask::IOM_VERTNORMAL)) 	fprintf(fpout,"N");
					if( tri::HasPerVertexColor(m)   && (mask & io::Mask::IOM_VERTCOLOR))		fprintf(fpout,"C");
					if( tri::HasPerVertexTexCoord(m) && (mask & io::Mask::IOM_VERTTEXCOORD))	fprintf(fpout,"ST");
					fprintf(fpout,"OFF\n");
					fprintf(fpout,"%d %d 0\n", m.vn, m.fn); // note that as edge number we simply write zero
          typename SaveMeshType::FaceIterator fi;
					
							//vertices
					int j;
					std::vector<int> FlagV; 
					VertexPointer  vp;
					VertexIterator vi;
					for(j=0,vi=m.vert.begin();vi!=m.vert.end();++vi)
					{
						vp=&(*vi);
            FlagV.push_back(vp->UberFlags()); // Salva in ogni caso flag del vertice
            if( ! vp->IsD() )
            {	// ***** ASCII *****

              fprintf(fpout,"%g %g %g " ,vp->P()[0],vp->P()[1],vp->P()[2]);
              if( tri::HasPerVertexColor(m)  && (mask & io::Mask::IOM_VERTCOLOR) )
                fprintf(fpout,"%d %d %d %d ",vp->C()[0],vp->C()[1],vp->C()[2],vp->C()[3] );

              if( m.HasPerVertexNormal()  && (mask & io::Mask::IOM_VERTNORMAL) )
                fprintf(fpout,"%g %g %g ", vp->N()[0],vp->N()[1],vp->N()[2]);

              if( m.HasPerVertexTexCoord()  && (mask & io::Mask::IOM_VERTTEXCOORD) )
                fprintf(fpout,"%g %g ",vp->T().u(),vp->T().v());
								
								fprintf(fpout,"\n");
								

              vp->UberFlags()=j; // Trucco! Nascondi nei flags l'indice del vertice non deletato!
              j++;
            }
					}

          assert(j==m.vn);
					FacePointer fp;
//					int vv[3];

					int fcnt=0;
					for(j=0,fi=m.face.begin();fi!=m.face.end();++fi)
					{
						fp=&(*fi);
						if( ! fp->IsD() )
						{ fcnt++;


						fprintf(fpout,"3 %d %d %d\n",
							fp->cV(0)->UberFlags(),	fp->cV(1)->UberFlags(), fp->cV(2)->UberFlags() );
						}
					}


					fclose(fpout);
					// Recupera i flag originali
					for(j=0,vi=m.vert.begin();vi!=m.vert.end();++vi)
						(*vi).UberFlags()=FlagV[j++]; 

					return 0;
				}

        static const char *ErrorMsg(int error)
        {
          static std::vector<std::string> off_error_msg;
          if(off_error_msg.empty())
          {
            off_error_msg.resize(2 );
            off_error_msg[0]="No errors";
	          off_error_msg[1]="Can't open file";
            }

          if(error>1 || error<0) return "Unknown error";
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
          capability |= vcg::tri::io::Mask::IOM_FACEINDEX;
	        return capability;
        }

			}; // end class
		} // end namespace tri
	} // end namespace io
} // end namespace vcg
//@}
#endif
