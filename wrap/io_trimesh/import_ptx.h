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
#ifndef __VCGLIB_IMPORT_PTX
#define __VCGLIB_IMPORT_PTX

#include <io.h>
#include <stdio.h>
#include <wrap/callback.h>
#include <vcg/complex/trimesh/allocate.h>
#include <vcg/complex/trimesh/clean.h>

namespace vcg {
 namespace tri {
  namespace io {

   /** 
   This class encapsulate a filter for importing ptx meshes.
   */
   template <class OpenMeshType>
   class ImporterPTX
   {
   public:

    enum PTX_OPEN_MASK_ENUM
    {
     PTX_ONLY_POINTS   		= 0x08000000,  //BIT_27 no add faces (PTX_FLIPFACES and PTX_SWITCHSIDE are ignored!)
     PTX_COLOR         		= 0x10000000,  //BIT_28 must be VertexType::HasColor();
     PTX_COMPUTE_AABBOX 	= 0x20000000,  //BIT_29 compute axis aligned bbox
     PTX_FLIPFACES				= 0x40000000,  //BIT_30 flip all faces ( PTX_ONLY_POINTS must be false )
     PTX_SWITCHSIDE    		= 0x80000000   //BIT_31 inverse triangulation order (swaping row->cols) ( PTX_ONLY_POINTS must be false )
    };

    typedef typename OpenMeshType::VertexPointer VertexPointer;
    typedef typename OpenMeshType::ScalarType ScalarType;
    typedef typename OpenMeshType::VertexType VertexType;
    typedef typename OpenMeshType::FaceType FaceType;
    typedef typename OpenMeshType::VertexIterator VertexIterator;
    typedef typename OpenMeshType::FaceIterator FaceIterator;



    // skip a mesh 
    static bool skipmesh(FILE* fp)
    {
     int colnum;
     int rownum;
     int skiplines;
     char linebuf;

     if(feof(fp))	return false;

     // getting mesh size;
     fscanf(fp,"%i\n",&colnum);
     fscanf(fp,"%i\n",&rownum);

     if ( ( colnum <=0 ) || ( rownum <=0 ) ) return false;

     //printf("\n %i x %i \n", rownum, colnum);
     //printf(" expect V %i F %i\n",(rownum*colnum),((rownum-1)*(colnum-1)*2));

     if(feof(fp))	return false;

     skiplines = (colnum * rownum) + 8; // have to skip (col * row) lines plus 8 lines for the header
     for(int ii=0; ii<skiplines; ii++)
     {
      fread(&linebuf,1,1,fp);
      while(linebuf != '\n')  fread(&linebuf,1,1,fp);
     } 

     return true;
    }


    
    //if numMesh == -1 load all mesh
    static bool Open( OpenMeshType &m, const char * filename, int numMesh = -1, int mask = PTX_ONLY_POINTS, CallBackPos *cb=NULL)
    {
     FILE *fp;
     fp = fopen(filename, "rb");
     if(fp == NULL) return false;

     m.Clear();

     if ( numMesh>0 )
      for (int i=0; i!=numMesh; ++i)  if (!skipmesh(fp)) return false;

     int mn=0;
     if ( numMesh == -1 )
     {
      bool next = true;
      while ( next )
      {
       bool r = readPTX(m, fp, mask, mn, cb); 
       mn++;
       if ((r==false) && (m.vn==0) ) { fclose(fp); return false; }
       else if (r==false) next = false;
      }
     } else 
     {   
      bool r = readPTX(m, fp, mask, numMesh, cb); 
      if ((r==false) && (m.vn==0) ) { fclose(fp); return false; }
     }

     fclose(fp);

     // now i delete all points in (0,0,0) that are unsampled points
     for(VertexIterator vi = m.vert.begin(); vi != m.vert.end(); vi++)
     {
      if((*vi).P() == Point3f(0.0, 0.0, 0.0))
      {
       (*vi).SetD();
       m.vn--;
      }
     }

     bool onlypoints  =  ((mask & PTX_ONLY_POINTS) != 0);
     if(! onlypoints)
     {


      for(OpenMeshType::FaceIterator fi = m.face.begin(); fi != m.face.end(); fi++)
      {
       if( ((*fi).V(0)->IsD()) || ((*fi).V(1)->IsD()) || ((*fi).V(2)->IsD()) )
       {
        (*fi).SetD();
        m.fn--;
       }
      }

      // eliminate high angle triangles
      /*int angle =85;
      if(angle != 90)
      {
      printf(" culling by angle \n");
      float limit = cos( angle*3.14159265358979323846/180.0 );
      Point3f raggio;

      vcg::tri::UpdateNormals<OpenMeshType>::PerFaceNormalized(m);
      for(OpenMeshType::FaceIterator fi = m.face.begin(); fi != m.face.end(); fi++)
      if(!(*fi).IsD())
      {
      raggio = -((*fi).V(0)->P() + (*fi).V(1)->P() + (*fi).V(2)->P()) / 3.0;
      raggio.Normalize();
      if((raggio * (*fi).N()) < limit)
      {
      (*fi).SetD();
      m.fn--;
      }
      }

      }*/
     }
     if(cb) cb(60,"PTX Mesh Loading RemoveDuplicateVertex");	
     tri::Clean<OpenMeshType>::RemoveDuplicateVertex(m);

     if (!onlypoints) 
     { 
      if(cb) cb(60,"PTX Mesh Loading RemoveUnreferencedVertex");	
      tri::Clean<OpenMeshType>::RemoveUnreferencedVertex(m);
     }
     if(cb) cb(100,"PTX Mesh Loading finish!");	
     return true;
    }




    static bool readPTX( OpenMeshType &m, FILE *fp, int mask, int mn, CallBackPos *cb=NULL)
    {
     int colnum;
     int rownum;
     int numtokens;
     char linebuf[256];
     int ii;
     float xx,yy,zz;	 // position
     float rr,gg,bb;	 // color
     float rf;		     // reflectance
     Matrix44f		currtrasf;

     bool hascolor;

     bool savecolor   =  ((mask & PTX_COLOR) != 0) &&  VertexType::HasColor();
     bool computeBbox =  ((mask & PTX_COMPUTE_AABBOX) != 0);
     bool onlypoints  =  ((mask & PTX_ONLY_POINTS) != 0);
     bool switchside  =  ((mask & PTX_SWITCHSIDE) != 0);
     bool flipfaces   =  ((mask & PTX_FLIPFACES) != 0);
     int total = 50;
     if ( onlypoints ) total = 100;


     // getting mesh size;
     fscanf(fp,"%i\n",&colnum);
     fscanf(fp,"%i\n",&rownum);

     if ( ( colnum <=0 ) || ( rownum <=0 ) ) return false;

     // initial 4 lines [still don't know what is this :) :)]
     if ( !fscanf(fp,"%f %f %f\n", &xx, &yy, &zz) ) return false;
     if ( !fscanf(fp,"%f %f %f\n", &xx, &yy, &zz) ) return false;
     if ( !fscanf(fp,"%f %f %f\n", &xx, &yy, &zz) ) return false;
     if ( !fscanf(fp,"%f %f %f\n", &xx, &yy, &zz) ) return false;

     // now the transformation matrix
     if ( !fscanf(fp,"%f %f %f %f\n", &(currtrasf.ElementAt(0,0)), &(currtrasf.ElementAt(0,1)), &(currtrasf.ElementAt(0,2)), &(currtrasf.ElementAt(0,3))) )return false;
     if ( !fscanf(fp,"%f %f %f %f\n", &(currtrasf.ElementAt(1,0)), &(currtrasf.ElementAt(1,1)), &(currtrasf.ElementAt(1,2)), &(currtrasf.ElementAt(1,3))) )return false;
     if ( !fscanf(fp,"%f %f %f %f\n", &(currtrasf.ElementAt(2,0)), &(currtrasf.ElementAt(2,1)), &(currtrasf.ElementAt(2,2)), &(currtrasf.ElementAt(2,3))) )return false;
     if ( !fscanf(fp,"%f %f %f %f\n", &(currtrasf.ElementAt(3,0)), &(currtrasf.ElementAt(3,1)), &(currtrasf.ElementAt(3,2)), &(currtrasf.ElementAt(3,3))) )return false;

     // now the real data begins

     // first line, we should know if the format is
     // XX YY ZZ RF
     // or it is
     // XX YY ZZ RF RR GG BB

     // read the entire first line and then count the spaces. it's rude but it works :)
     ii=0;
     fread(&(linebuf[ii++]),1,1,fp);
     while(linebuf[ii-1] != '\n')  if ( fread(&(linebuf[ii++]),1,1,fp)==0 ) return false;
     linebuf[ii-1] = '\0'; // terminate the string
     numtokens=1;
     for(ii=0; ii<(int)strlen(linebuf); ii++) if(linebuf[ii] == ' ') numtokens++;
     if(numtokens == 4)  hascolor = false;
     else if(numtokens == 7)  hascolor = true;
     else  return false;

     Transpose(currtrasf);
     int vn = rownum*colnum;

     VertexIterator vi = Allocator<OpenMeshType>::AddVertices(m,vn); 

     // parse the first line....
     if(hascolor)
     {
      printf("\n hascolor ");
      sscanf(linebuf,"%f %f %f %f %f %f %f", &xx, &yy, &zz, &rf, &rr, &gg, &bb);
     }
     else
     {
      printf("\n no color ");
      sscanf(linebuf,"%f %f %f %f", &xx, &yy, &zz, &rf);
     }

     if (computeBbox) m.bbox.SetNull();


     //addthefirstpoint
     (*vi).P()[0]=xx;
     (*vi).P()[1]=yy;
     (*vi).P()[2]=zz;
     (*vi).P() = currtrasf * (*vi).P();
     if (computeBbox) m.bbox.Add( (*vi).P() );
     if(hascolor && savecolor)
     {
      (*vi).C()[0]=rr;
      (*vi).C()[1]=gg;
      (*vi).C()[2]=bb;
     }
     vi++;



     // now for each line until end of mesh (row*col)-1
     for(ii=0; ii<((rownum*colnum)-1); ii++)
     {

      char tmp[255];
      sprintf(tmp, "PTX Mesh Loading... mesh %i", mn);
      if(cb) cb((ii*total)/vn, tmp);	

      // read the stream
      if(hascolor)   fscanf(fp,"%f %f %f %f %f %f %f", &xx, &yy, &zz, &rf, &rr, &gg, &bb);
      else  fscanf(fp,"%f %f %f %f", &xx, &yy, &zz, &rf);


      // add the point
      (*vi).P()[0]=xx;
      (*vi).P()[1]=yy;
      (*vi).P()[2]=zz;
      (*vi).P() = currtrasf * (*vi).P();
      if (computeBbox) m.bbox.Add( (*vi).P() );
      if(hascolor && savecolor)
      {
       (*vi).C()[0]=rr;
       (*vi).C()[1]=gg;
       (*vi).C()[2]=bb;
      }
      vi++;



     }


     m.vn = m.vert.size();
     m.fn = 0;

     if(! onlypoints)
     {

      // now i can triangulate
      int trinum = (rownum-1) * (colnum-1) * 2;


      OpenMeshType::FaceIterator fi= Allocator<OpenMeshType>::AddFaces(m,trinum);


      m.fn = 0;

      int v0i,v1i,v2i, t;
      for(int rit=0; rit<rownum-1; rit++)
       for(int cit=0; cit<colnum-1; cit++)
       {

        if(cb) cb(50 + (t*50)/(rownum*colnum),"PTX Mesh Loading");	

        if(!switchside)
        {
         v0i = (rit  ) + ((cit  ) * rownum);
         v1i = (rit+1) + ((cit  ) * rownum);
         v2i = (rit  ) + ((cit+1) * rownum);
        }
        else
        {
         v0i = (cit  ) + ((rit  ) * colnum);
         v1i = (cit+1) + ((rit  ) * colnum);
         v2i = (cit  ) + ((rit+1) * colnum);
        }


        // upper tri
        (*fi).V(2) = &(m.vert[v0i]);
        (*fi).V(1) = &(m.vert[v1i]);
        (*fi).V(0) = &(m.vert[v2i]);

        if(flipfaces)
        {
         (*fi).V(2) = &(m.vert[v1i]);
         (*fi).V(1) = &(m.vert[v0i]);
        }

        m.fn++;
        fi++;

        if(!switchside)
        {
         v0i = (rit+1) + ((cit  ) * rownum);
         v1i = (rit+1) + ((cit+1) * rownum);
         v2i = (rit  ) + ((cit+1) * rownum);
        }
        else
        {
         v0i = (cit+1) + ((rit  ) * colnum);
         v1i = (cit+1) + ((rit+1) * colnum);
         v2i = (cit  ) + ((rit+1) * colnum);
        }

        // lower tri
        (*fi).V(2) = &(m.vert[v0i]);
        (*fi).V(1) = &(m.vert[v1i]);
        (*fi).V(0) = &(m.vert[v2i]);

        if(flipfaces)
        {
         (*fi).V(2) = &(m.vert[v1i]);
         (*fi).V(1) = &(m.vert[v0i]);
        }

        m.fn++;
        fi++;
       }



     }

     return true;




    }



   }; // end class


  } // end Namespace tri
 } // end Namespace io
} // end Namespace vcg

#endif
