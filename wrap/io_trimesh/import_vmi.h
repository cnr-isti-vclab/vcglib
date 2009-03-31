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
 Revision 1.1  2007/02/14 01:20:37  ganovelli
 working draft of VCG Mesh Image importer and exporter. Does not consider optional attributes. The mesh atributes are only vn and fn (no bbox, texture coordiantes)

 

****************************************************************************/

#ifndef __VCGLIB_IMPORT_VMI
#define __VCGLIB_IMPORT_VMI

/*
	VMI VCG Mesh Image.
	The vmi image file consists of a header containing the description of the vertex and face type,
	the length of vectors containing vertices of faces and the memory image of the object mesh as it is when
	passed to the function Save(SaveMeshType m)
	NOTE: THIS IS NOT A FILE FORMAT. IT IS ONLY USEFUL FOR DUMPING MESH IMAGES FOR DEBUG PURPOSE.
	Example of use: say you are running a time consuming mesh processing and you want to save intermediate
	state, but no file format support all the attributes you need in your vertex/face type. 
	NOTE2: At the present if you add members to your TriMesh these will NOT be saved. More precisely, this file and
	import_vmi must be updated to reflect changes in vcg/complex/trimesh/base.h

*/

namespace vcg {
namespace tri {
namespace io {

	template <class OpenMeshType>
	class ImporterVMI
	{
	public:	
		typedef typename OpenMeshType::FaceIterator FaceIterator;
		typedef typename OpenMeshType::VertexIterator VertexIterator;
		typedef typename OpenMeshType::VertexType VertexType;

		static bool GetHeader(FILE * f,std::vector<std::string>& fnameV, std::vector<std::string>& fnameF, int & vertSize, int &faceSize){
			char name[100];
			char buf[4096];
			int nameFsize,nameVsize,i;
			fgets(buf,4096,f);
			sscanf(buf,"%s %d",&name[0],&nameFsize);
			for(i=0; i < nameFsize; ++i) {
				fgets(buf,4096,f);
				sscanf(buf,"%s ",&name[0]);fnameF.push_back(std::string(name));
			}
			fgets(buf,4096,f);
			sscanf(buf,"%s %d",&name[0],&faceSize);

			fgets(buf,4096,f);
			sscanf(buf,"%s %d",&name[0],&nameVsize);
			for(i=0; i < nameVsize; ++i) {
				fgets(buf,4096,f);
				sscanf(buf,"%s ",&name[0]);fnameV.push_back(std::string(name));}
			fgets(buf,4096,f);
			sscanf(buf,"%s %d",&name[0],&vertSize);
			fgets(buf,4096,f);
			assert(strstr(buf,"end_header")!=NULL);
			return true;
		}

		static bool GetHeader(char * filename,std::vector<std::string>& nameV, std::vector<std::string>& nameF, int & vertSize, int &faceSize){
				FILE * f = fopen(filename,"rb");
				return GetHeader(f,nameV, nameF, vertSize, faceSize);
				fclose(f);
	}
		static bool Open(OpenMeshType &m,char * filename){
			int i;
			typedef typename OpenMeshType::VertexType VertexType; 	
			typedef typename OpenMeshType::FaceType FaceType; 	
			typename OpenMeshType::FaceIterator fi;
			typename OpenMeshType::VertexIterator vi;
			FILE * f = fopen(filename,"rb");
			std::vector<std::string> nameF,nameV,fnameF,fnameV;
			int vertSize,faceSize;

			/* read the header */
			GetHeader(f,fnameV, fnameF, vertSize, faceSize);

			/* read the mesh type */
			OpenMeshType::FaceType::Name(nameF);	
			OpenMeshType::VertexType::Name(nameV);

			/* check if the type is the very same, otherwise return */
			if(fnameV != nameV) return false;
			if(fnameF != nameF) return false;

			 int offsetV,offsetF;
			/* read the address of the first vertex */
			fread(&offsetV,sizeof( int),1,f);

			/* read the address of the first face */
			fread(&offsetF,sizeof( int),1,f);

			/* read the object mesh */
			fread(&m.camera,sizeof(Camera<typename OpenMeshType::ScalarType>),1,f);
			fread(&m.shot,sizeof(Shot<typename OpenMeshType::ScalarType>),1,f);
			fread(&m.vn,sizeof(int),1,f);
			fread(&m.fn,sizeof(int),1,f);
			fread(&m.imark,sizeof(int),1,f);
			fread(&m.bbox,sizeof(Box3<typename OpenMeshType::ScalarType>),1,f);
			fread(&m.C(),sizeof(Color4b),1,f);

			/* resize the vector of vertices */
			m.vert.resize(vertSize);

			int read = 0;
			/* load the vertices */
			if(vertSize>0)
				read=fread((void*)& m.vert[0],sizeof(VertexType),vertSize,f);
			assert(ferror(f)==0);
			assert(read==vertSize);

			read = 0;
			m.face.resize(faceSize);
			if(faceSize>0)
				/* load the faces */
				read = fread((void*)& m.face[0],sizeof(FaceType),faceSize,f);
			assert(ferror(f)==0);
			assert(!feof(f));
			assert(read==faceSize);
			
			if(FaceType::HasVFAdjacency())
				for(vi = m.vert.begin(); vi != m.vert.end(); ++vi){
					(*vi).VFp() = (*vi).VFp()-(FaceType*)offsetF+ &m.face[0];
					(*vi).VFp() = (*vi).VFp()-(FaceType*)offsetF+ &m.face[0];
					(*vi).VFp() = (*vi).VFp()-(FaceType*)offsetF+ &m.face[0];
				}

			if(FaceType::HasVertexRef())
				for(fi = m.face.begin(); fi != m.face.end(); ++fi){
					(*fi).V(0) = (*fi).V(0)-(VertexType*)offsetV+ &m.vert[0];
					(*fi).V(1) = (*fi).V(1)-(VertexType*)offsetV+ &m.vert[0];
					(*fi).V(2) = (*fi).V(2)-(VertexType*)offsetV+ &m.vert[0];
				}

			if(FaceType::HasFFAdjacency())
				for(fi = m.face.begin(); fi != m.face.end(); ++fi){
					(*fi).FFp(0) = (*fi).FFp(0)-(FaceType*)offsetF+ &m.face[0];
					(*fi).FFp(1) = (*fi).FFp(1)-(FaceType*)offsetF+ &m.face[0];
					(*fi).FFp(2) = (*fi).FFp(2)-(FaceType*)offsetF+ &m.face[0];
				}

			if(FaceType::HasVFAdjacency())
				for(fi = m.face.begin(); fi != m.face.end(); ++fi){
					(*fi).VFp(0) = (*fi).VFp(0)-(FaceType*)offsetF+ &m.face[0];
					(*fi).VFp(1) = (*fi).VFp(1)-(FaceType*)offsetF+ &m.face[0];
					(*fi).VFp(2) = (*fi).VFp(2)-(FaceType*)offsetF+ &m.face[0];
				}

				fclose(f);
				return true;
		}

	}; // end class

} // end Namespace tri
} // end Namespace io
} // end Namespace vcg

#endif
