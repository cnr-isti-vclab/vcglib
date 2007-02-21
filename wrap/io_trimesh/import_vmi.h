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
			OpenMeshType::FaceIterator fi;
			OpenMeshType::VertexIterator vi;
			FILE * f = fopen(filename,"rb");
			std::vector<string> nameF,nameV,fnameF,fnameV;
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
			fread(&m,sizeof(OpenMeshType),1,f);

			/* overwrite che container because they could be inconsistent */
			OpenMeshType::VertContainer tvc;
			OpenMeshType::FaceContainer tfc;
			memcpy(&m.vert,&tvc,sizeof(OpenMeshType::VertContainer));
			memcpy(&m.face,&tfc,sizeof(OpenMeshType::FaceContainer));

			m.vert.resize(vertSize);

			int read;
			/* load the vertices */
			read=fread((void*)& m.vert[0],sizeof(OpenMeshType::VertexType),vertSize,f);
			assert(ferror(f)==0);
			assert(read==vertSize);

			m.face.resize(faceSize);
			/* load the faces */
			read = fread((void*)& m.face[0],sizeof(OpenMeshType::FaceType),faceSize,f);
			assert(ferror(f)==0);
			assert(!feof(f));
			assert(read==faceSize);
			
			if(OpenMeshType::FaceType::HasVFAdjacency())
				for(vi = m.vert.begin(); vi != m.vert.end(); ++vi){
					(*vi).VFp() = (*vi).VFp()-(OpenMeshType::FaceType*)offsetF+ &m.face[0];
					(*vi).VFp() = (*vi).VFp()-(OpenMeshType::FaceType*)offsetF+ &m.face[0];
					(*vi).VFp() = (*vi).VFp()-(OpenMeshType::FaceType*)offsetF+ &m.face[0];
				}

			if(OpenMeshType::FaceType::HasVertexRef())
				for(fi = m.face.begin(); fi != m.face.end(); ++fi){
					(*fi).V(0) = (*fi).V(0)-(OpenMeshType::VertexType*)offsetV+ &m.vert[0];
					(*fi).V(1) = (*fi).V(1)-(OpenMeshType::VertexType*)offsetV+ &m.vert[0];
					(*fi).V(2) = (*fi).V(2)-(OpenMeshType::VertexType*)offsetV+ &m.vert[0];
				}

			if(OpenMeshType::FaceType::HasFFAdjacency())
				for(fi = m.face.begin(); fi != m.face.end(); ++fi){
					(*fi).FFp(0) = (*fi).FFp(0)-(OpenMeshType::FaceType*)offsetF+ &m.face[0];
					(*fi).FFp(1) = (*fi).FFp(1)-(OpenMeshType::FaceType*)offsetF+ &m.face[0];
					(*fi).FFp(2) = (*fi).FFp(2)-(OpenMeshType::FaceType*)offsetF+ &m.face[0];
				}

			if(OpenMeshType::FaceType::HasVFAdjacency())
				for(fi = m.face.begin(); fi != m.face.end(); ++fi){
					(*fi).VFp(0) = (*fi).VFp(0)-(OpenMeshType::FaceType*)offsetF+ &m.face[0];
					(*fi).VFp(1) = (*fi).VFp(1)-(OpenMeshType::FaceType*)offsetF+ &m.face[0];
					(*fi).VFp(2) = (*fi).VFp(2)-(OpenMeshType::FaceType*)offsetF+ &m.face[0];
				}

				fclose(f);
		}

	}; // end class

} // end Namespace tri
} // end Namespace io
} // end Namespace vcg

#endif
