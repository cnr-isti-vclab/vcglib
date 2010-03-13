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
 Revision 1.10  2007/10/17 09:49:50  cignoni
 correct management of point only files

 Revision 1.9  2007/03/12 16:40:16  tarini
 Texture coord name change!  "TCoord" and "Texture" are BAD. "TexCoord" is GOOD.

 Revision 1.8  2007/03/08 11:27:52  ganovelli
 added  include to tcoord2

 Revision 1.7  2007/02/16 21:12:13  m_di_benedetto
 Commented out strange abort in WriteMaterials()

 Revision 1.6  2006/12/07 00:37:58  cignoni
 Corrected bug in the management of deleted vertices

 Revision 1.5  2006/10/09 19:58:08  cignoni
 Added casts to remove warnings

 Revision 1.4  2006/09/18 12:14:38  cignoni
 Removed bug in the creation of the material filename

 Revision 1.3  2006/03/07 13:19:29  cignoni
 First Release with OBJ import support

 Revision 1.2  2006/02/28 14:38:09  corsini
 remove qt include

 Revision 1.1  2006/02/16 19:28:36  fmazzant
 transfer of Export_3ds.h, Export_obj.h, Io_3ds_obj_material.h from Meshlab to vcg

 Revision 1.7  2006/02/06 11:04:40  fmazzant
 added file material.h. it include struct Material, CreateNewMaterial(...) and MaterialsCompare(...)

 Revision 1.6  2006/02/04 10:18:46  fmazzant
 clean code

 Revision 1.5  2006/02/03 10:04:41  fmazzant
 no significant updated

 Revision 1.4  2006/01/30 14:02:05  fmazzant
 bug-fix

 Revision 1.3  2006/01/29 23:52:43  fmazzant
 correct a small bug

 Revision 1.2  2006/01/29 18:33:42  fmazzant
 added some comment to the code

 Revision 1.1  2006/01/29 16:33:03  fmazzant
 moved export_obj and export_3ds from test/io into meshio/

 Revision 1.34  2006/01/22 23:59:01  fmazzant
 changed default value of diffuse. 1.0 -> 0.8

 Revision 1.33  2006/01/19 09:36:29  fmazzant
 cleaned up history log

 Revision 1.32  2006/01/18 00:45:56  fmazzant
 added control on face's diffuse

 Revision 1.31  2006/01/17 13:48:54  fmazzant
 added capability mask on export file format

 Revision 1.30  2006/01/15 00:45:40  fmazzant
 extend mask exporter for all type file format +

 Revision 1.29  2006/01/14 00:03:26  fmazzant
 added more controls

****************************************************************************/

#ifndef __VCGLIB_EXPORT_OBJ
#define __VCGLIB_EXPORT_OBJ

#include <wrap/callback.h>
#include <vcg/complex/trimesh/allocate.h>
#include <vcg/space/texcoord2.h>
#include <wrap/io_trimesh/io_mask.h>
#include "io_material.h"
#include <iostream>
#include <fstream>
#include <map>

namespace vcg {
namespace tri {
namespace io {

	template <class SaveMeshType>
	class ExporterOBJ
	{
	public:	
		typedef typename SaveMeshType::FaceIterator FaceIterator;
		typedef typename SaveMeshType::VertexIterator VertexIterator;
		typedef typename SaveMeshType::VertexType VertexType;
	
		/*
			enum of all the types of error
		*/
		enum SaveError
		{
			E_NOERROR,					// 0
			E_CANTOPENFILE,				// 1
			E_CANTCLOSEFILE,			// 2
			E_UNESPECTEDEOF,			// 3
			E_ABORTED,					// 4
			E_NOTDEFINITION,			// 5
			E_NOTVEXTEXVALID,			// 6
			E_NOTFACESVALID				// 7
		};

		/*
			this function takes an index and the relative error message gets back
		*/
		static const char* ErrorMsg(int error)
		{
			static const char* obj_error_msg[] =
			{
					"No errors",							// 0
					"Can't open file",						// 1
					"can't close file",						// 2
					"Premature End of file",				// 3
					"File saving aborted",					// 4
					"Function not defined",					// 5
					"Vertices not valid",					// 6
					"Faces not valid"						// 7
				};

			if(error>7 || error<0) return "Unknown error";
			else return obj_error_msg[error];
		};

		/*
			returns mask of capability one define with what are the saveable information of the format.
		*/
		static int GetExportMaskCapability()
		{
			int capability = 0;
			
			//vert
			capability |= vcg::tri::io::Mask::IOM_VERTNORMAL;
      capability |= vcg::tri::io::Mask::IOM_VERTTEXCOORD;

			//face
			capability |= vcg::tri::io::Mask::IOM_FACECOLOR;

			//wedg
			capability |= vcg::tri::io::Mask::IOM_WEDGTEXCOORD;
			capability |= vcg::tri::io::Mask::IOM_WEDGNORMAL;

			return capability;
		}

		/*
			function which saves in OBJ file format
		*/
		static int Save(SaveMeshType &m, const char * filename, int mask, CallBackPos *cb=0)	
		{
      // texture coord and normal: cannot be saved BOTH per vertex and per wedge
      if (mask & vcg::tri::io::Mask::IOM_WEDGTEXCOORD &&
          mask & vcg::tri::io::Mask::IOM_VERTTEXCOORD ) {
        mask &= ~vcg::tri::io::Mask::IOM_VERTTEXCOORD;
      }
      if (mask & vcg::tri::io::Mask::IOM_WEDGCOLOR &&
          mask & vcg::tri::io::Mask::IOM_VERTCOLOR ) {
        mask &= ~vcg::tri::io::Mask::IOM_VERTCOLOR;
      }
      if(m.vn == 0)	return E_NOTVEXTEXVALID;
			// Commented out this control. You should be allowed to save a point cloud.
			// if(m.fn == 0)	return E_NOTFACESVALID;

			int current = 0;
			int totalPrimitives = m.vn+m.fn;
		
			std::vector<Material> materialVec;
			
			std::string fn(filename);
			int LastSlash=fn.size()-1;
			while(LastSlash>=0 && fn[LastSlash]!='/')
        --LastSlash;

			FILE *fp;
			fp = fopen(filename,"w");
			if(fp == NULL) return E_CANTOPENFILE;

			fprintf(fp,"####\n#\n# OBJ File Generated by Meshlab\n#\n####\n");
			fprintf(fp,"# Object %s\n#\n# Vertices: %d\n# Faces: %d\n#\n####\n",fn.substr(LastSlash+1).c_str(),m.vn,m.fn);
			
			//library materialVec
			if( (mask & vcg::tri::io::Mask::IOM_FACECOLOR)  || (mask & Mask::IOM_WEDGTEXCOORD) )
				fprintf(fp,"mtllib ./%s.mtl\n\n",fn.substr(LastSlash+1).c_str());
			
			//vertexs + normal
			VertexIterator vi;
			std::map<Point3f,int> NormalVertex;
      std::vector<int> VertexId(m.vert.size());
			int numvert = 0;
			int curNormalIndex = 1;
			for(vi=m.vert.begin(); vi!=m.vert.end(); ++vi) if( !(*vi).IsD() )
			{
        VertexId[vi-m.vert.begin()]=numvert;
				//saves normal per vertex
        if (mask & Mask::IOM_WEDGNORMAL )
				{
					if(AddNewNormalVertex(NormalVertex,(*vi).N(),curNormalIndex))
					{
						fprintf(fp,"vn %f %f %f\n",(*vi).N()[0],(*vi).N()[1],(*vi).N()[2]);
						curNormalIndex++;
					}
				}
        if (mask & Mask::IOM_VERTNORMAL ) {
          fprintf(fp,"vn %f %f %f\n",(*vi).N()[0],(*vi).N()[1],(*vi).N()[2]);
        }
        if (mask & Mask::IOM_VERTTEXCOORD ) {
          fprintf(fp,"vt %f %f\n",(*vi).T().P()[0],(*vi).T().P()[1]);
        }
        //if (mask & Mask::IOM_VERTCOLOR ) {
        //  fprintf(fp,"vc %f %f %f\n",(*vi).T().P()[0],(*vi).T().P()[1]);
        //}

				//saves vertex
				fprintf(fp,"v %f %f %f\n",(*vi).P()[0],(*vi).P()[1],(*vi).P()[2]);

				if (cb !=NULL) {
          if(!(*cb)((100*++current)/totalPrimitives, "writing vertices ")) 
            { fclose(fp); return E_ABORTED;} }
        numvert++;
			}
      assert(numvert == m.vn);

			fprintf(fp,"# %d vertices, %d vertices normals\n\n",m.vn,int(NormalVertex.size()));
			
			//faces + texture coords
			FaceIterator fi;
			std::map<vcg::TexCoord2<float>,int> CoordIndexTexture;
			unsigned int material_num = 0;
			int mem_index = 0; //var temporany
			int curTexCoordIndex = 1; 
			for(fi=m.face.begin(); fi!=m.face.end(); ++fi) if( !(*fi).IsD() )
			{
				if((mask & Mask::IOM_FACECOLOR) || (mask & Mask::IOM_WEDGTEXCOORD) )
				{
					int index = Materials<SaveMeshType>::CreateNewMaterial(m,materialVec,material_num,fi);
					
					if(index == (int)materialVec.size())//inserts a new element material
					{
						material_num++;
						fprintf(fp,"\nusemtl material_%d\n",materialVec[index-1].index);
						mem_index = index-1;
					}
					else
					{
						if(index != mem_index)//inserts old name elemente material
						{
							fprintf(fp,"\nusemtl material_%d\n",materialVec[index].index);
							mem_index=index;
						}
					}
				}

        //saves texture coord x wedge
        if(HasPerWedgeTexCoord(m) && (mask & Mask::IOM_WEDGTEXCOORD))
        for(unsigned int k=0;k<3;k++)
				{
					{
						if(AddNewTextureCoord(CoordIndexTexture,(*fi).WT(k),curTexCoordIndex))
						{
							fprintf(fp,"vt %f %f\n",(*fi).WT(k).u(),(*fi).WT(k).v());
							curTexCoordIndex++; //ncreases the value number to be associated to the Texture
						}
					}
				}


				fprintf(fp,"f ");
				for(unsigned int k=0;k<3;k++)
				{
				if(k!=0) fprintf(fp," ");	
					int vInd = -1; 
					// +1 because Obj file format begins from index = 1 but not from index = 0.
					vInd = VertexId[GetIndexVertex(m, (*fi).V(k))] + 1;//index of vertex per face
					
					int vt = -1;
					if(mask & Mask::IOM_WEDGTEXCOORD)
						vt = GetIndexVertexTexture(CoordIndexTexture,(*fi).WT(k));//index of vertex texture per face
          if (mask & Mask::IOM_VERTTEXCOORD)
            vt = vInd;

					int vn = -1;
          if(mask & Mask::IOM_WEDGNORMAL )
						vn = GetIndexVertexNormal(m, NormalVertex, (*fi).V(k)->cN());//index of vertex normal per face.
          if (mask & Mask::IOM_VERTNORMAL)
            vn = vInd;

					//writes elements on file obj
					WriteFacesElement(fp,vInd,vt,vn);
				}
				fprintf(fp,"\n");	
     		if (cb !=NULL) {
          if(!(*cb)((100*++current)/totalPrimitives, "writing vertices "))
              { fclose(fp); return E_ABORTED;} 
        }

			}//for
			fprintf(fp,"# %d faces, %d coords texture\n\n",m.fn,int(CoordIndexTexture.size()));
			
			fprintf(fp,"# End of File");
			fclose(fp);

			int r = 0;
			if((mask & Mask::IOM_WEDGTEXCOORD) || (mask & Mask::IOM_FACECOLOR) )
				r = WriteMaterials(materialVec, filename,cb);//write material 
			
			if(r!= E_NOERROR)
				return r;
			return E_NOERROR;
		}

		/*
			returns index of the vertex
		*/
		inline static int GetIndexVertex(SaveMeshType &m, VertexType *p)
		{
			return p-&*(m.vert.begin());
		}
		
		/*
			returns index of the texture coord
		*/
		inline static int GetIndexVertexTexture(std::map<TexCoord2<float>,int> &mapTexToInt, const vcg::TexCoord2<float> &wt)
		{
			std::map<vcg::TexCoord2<float>,int>::iterator iter= mapTexToInt.find(wt);
			if(iter != mapTexToInt.end()) return (*iter).second;
			else 		return -1;
			// Old wrong version.	
			// int index = mapTexToInt[wt];
			// if(index!=0){return index;}
		}

		/*
			returns index of the vertex normal
		*/
		inline static int GetIndexVertexNormal(SaveMeshType &/*m*/, std::map<Point3f,int> &mapNormToInt, const Point3f &norm )
		{
		  std::map<Point3f,int>::iterator iter= mapNormToInt.find(norm);
			if(iter != mapNormToInt.end()) return (*iter).second;
			else 		return -1;
			// Old wrong version.	
			// int index = mapNormToInt[m.vert[iv].N()];
			// if(index!=0){return index;}
		}
		
		/*
			write elements on file
	
				f v/tc/vn v/tc/vn v/tc/vn ...
				f v/tc v/tc v/tc ...
				f v//vn v//vn v//vn ...
				f v v v ...
				
		*/
		inline static void WriteFacesElement(FILE *fp,int v,int vt, int vn)
		{
			fprintf(fp,"%d",v);
			if(vt!=-1)
			{
				fprintf(fp,"/%d",vt);
				if(vn!=-1) 
					fprintf(fp,"/%d",vn);
			}
			else if(vn!=-1)
				fprintf(fp,"//%d",vn);
		}
		
		/*
			adds a new index to the coordinate of Texture if it is the first time 
			which is otherwise met not execute anything
		*/
		inline static bool AddNewTextureCoord(std::map<vcg::TexCoord2<float>,int> &m, const vcg::TexCoord2<float> &wt,int value)
		{
			int index = m[wt];
			if(index==0){m[wt]=value;return true;}
			return false;
		}

		/*
			adds a new index to the normal per vertex if it is the first time 
			which is otherwise met does not execute anything
		*/
		inline static bool AddNewNormalVertex(std::map<Point3f,int> &m, Point3f &n ,int value)
		{
			int index = m[n];
			if(index==0){m[n]=value;return true;}
			return false;
		}
		
		/*
			writes material into file
		*/
		inline static int WriteMaterials(std::vector<Material> &materialVec, const char * filename, CallBackPos *cb=0)
		{			
			std::string fileName = std::string(filename);
			fileName+=".mtl";
			
			if(materialVec.size() > 0)
			{
				FILE *fp;
				fp = fopen(fileName.c_str(),"w");
				if(fp==NULL)return E_ABORTED;
				
				fprintf(fp,"#\n# Wavefront material file\n# Converted by Meshlab Group\n#\n\n");
				
				int current = 0;

				for(unsigned int i=0;i<materialVec.size();i++)
				{
					if (cb !=NULL)
						(*cb)((100 * ++current)/materialVec.size(), "saving material file ");
					else
					{ /* fclose(fp); return E_ABORTED; */ }

					fprintf(fp,"newmtl material_%d\n",materialVec[i].index);
					fprintf(fp,"Ka %f %f %f\n",materialVec[i].Ka[0],materialVec[i].Ka[1],materialVec[i].Ka[2]);
					fprintf(fp,"Kd %f %f %f\n",materialVec[i].Kd[0],materialVec[i].Kd[1],materialVec[i].Kd[2]);
					fprintf(fp,"Ks %f %f %f\n",materialVec[i].Ks[0],materialVec[i].Ks[1],materialVec[i].Ks[2]);
					fprintf(fp,"Tr %f\n",materialVec[i].Tr);
					fprintf(fp,"illum %d\n",materialVec[i].illum);
					fprintf(fp,"Ns %f\n",materialVec[i].Ns);

					if(materialVec[i].map_Kd.size()>0)
						fprintf(fp,"map_Kd %s\n",materialVec[i].map_Kd.c_str());
					fprintf(fp,"\n");
				}
				fclose(fp);
			}
			return E_NOERROR;
		}

	}; // end class
} // end Namespace tri
} // end Namespace io
} // end Namespace vcg

#endif
