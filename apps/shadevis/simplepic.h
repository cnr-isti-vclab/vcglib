/*#***************************************************************************
 * Visibility.h                                                     o o      *
 *                                                                o     o    *
 * Visual Computing Group                                         _  O  _    *
 * IEI Institute, CNUCE Institute, CNR Pisa                        \/)\/     *
 *                                                                /\/|       *
 * Copyright(C) 1999 by Paolo Cignoni, Paolo Pingi, Claudio Rocchini |       *
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
/*#**************************************************************************
  History

 2003	Aug 26 First Working version
	    
****************************************************************************/
#ifndef __VCG_SIMPLE_PIC
#define __VCG_SIMPLE_PIC
#include <vcg/math/matrix44.h>

namespace vcg {
  template <class PixType> 
  class SimplePic
  {public:
    std::vector<PixType> img;
   int sx,sy;
   void Create(int tx,int ty)
   {
     sx=tx;sy=ty;
     img.resize(sx*sy);
   }
   PixType &Pix(int x, int y) {return img[sx*y+x];}

   void OpenGLSnap(GLenum format=0)
	{
		int vp[4];
		glGetIntegerv( GL_VIEWPORT,vp );		// Lettura viewport
		glPixelStorei( GL_PACK_ROW_LENGTH, 0);
		glPixelStorei( GL_PACK_ALIGNMENT, 1);
		int tx = vp[2];
		int ty = vp[3];

		Create(tx,ty);

		GLenum mtype  = 0;

		if(format==0) {
				format = GL_RGBA;
        mtype = GL_UNSIGNED_BYTE;
		}
		if(format==GL_DEPTH_COMPONENT) {
				format = GL_DEPTH_COMPONENT;
        mtype = GL_FLOAT;
		}
 		glReadPixels(vp[0],vp[1],vp[2],vp[3],format,mtype,(GLvoid *)&img[0]);
	}
	bool SavePPM( const char * filename )
	{
		FILE * fp = fopen(filename,"wb");
		if(fp==0) return false;


			fprintf(fp,"P6\n%d %d\n255\n",sx,sy);

			for(int i=0;i<sx*sy;++i)
			{
			 fwrite(&(img[i]),3,1,fp);
			}
	
		fclose(fp);
		return true;
	}

};

}
#endif // __VCG_MESH_VISIBILITY