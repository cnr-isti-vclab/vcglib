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

#ifndef __VCGLIB_EXPORT_VMI
#define __VCGLIB_EXPORT_VMI

/*
	VMI VCG Mesh Image.
	The vmi image file consists of a header containing the description of the vertex and face type,
	the length of vectors containing vertices of faces and the memory image of the object mesh as it is when
	passed to the function Save(SaveMeshType m).
	NOTE: THIS IS NOT A FILE FORMAT. IT IS ONLY USEFUL FOR DUMPING MESH IMAGES FOR DEBUG PURPOSE.
	Example of use: say you are running a time consuming mesh processing and you want to save intermediate
	state, but no file format support all the attributes you need in your vertex/face type. 
	NOTE2: At the present if you add members to your TriMesh these will NOT be saved. More precisely, this file and
	import_vmi must be updated to reflect changes in vcg/complex/trimesh/base.h
	*/

namespace vcg {
namespace tri {
namespace io {

	template <class SaveMeshType>
	class ExporterVMI
	{
	public:	
		typedef typename SaveMeshType::FaceIterator FaceIterator;
		typedef typename SaveMeshType::VertexIterator VertexIterator;
		typedef typename SaveMeshType::VertexType VertexType;

		static void Save(const SaveMeshType &m,char * filename){
			unsigned int i;
			int vertSize,faceSize;
			FILE * f = fopen(filename,"wb");
			std::vector<std::string> nameF,nameV;
			SaveMeshType::FaceType::Name(nameF);
			SaveMeshType::VertexType::Name(nameV);
			vertSize = m.vert.size();
			faceSize = m.face.size();

			/* write header */
			fprintf(f,"FACE_TYPE %d\n",nameF.size());
			for(i=0; i < nameF.size(); ++i) fprintf(f,"%s\n",nameF[i].c_str());
			fprintf(f,"SIZE_VECTOR_FACES %d\n",faceSize);

			fprintf(f,"VERTEX_TYPE %d\n",nameV.size());
			for(i=0; i < nameV.size(); ++i) fprintf(f,"%s\n",nameV[i].c_str());
			fprintf(f,"SIZE_VECTOR_VERTS %d\n",vertSize);
			fprintf(f,"end_header\n");


			unsigned int offsetV = (unsigned int) &m.vert[0];
			/* write the address of the first vertex */
			fwrite(&offsetV,sizeof(unsigned int),1,f);

			 int offsetF= ( int) &m.face[0];
			/* write the address of the first face */
			fwrite(&offsetF,sizeof( int),1,f);

			/* save the object mesh */
			fwrite(&m.camera,sizeof(Camera<typename SaveMeshType::ScalarType>),1,f);
			fwrite(&m.shot,sizeof(Shot<typename SaveMeshType::ScalarType>),1,f);
			fwrite(&m.vn,sizeof(int),1,f);
			fwrite(&m.fn,sizeof(int),1,f);
			fwrite(&m.imark,sizeof(int),1,f);
			fwrite(&m.bbox,sizeof(Box3<typename SaveMeshType::ScalarType>),1,f);
			fwrite(&m.C(),sizeof(Color4b),1,f);

			int written;
			/* save the vertices */
			written = fwrite((void*)&m.vert[0],sizeof(typename SaveMeshType::VertexType),m.vert.size(),f);
			assert(written==m.vert.size());

			/* save the faces */
			written = fwrite((void*)&m.face[0],sizeof(typename SaveMeshType::FaceType),m.face.size(),f);
			assert(written==m.face.size());

		//	fflush(f);
			fclose(f);
		}

	}; // end class

} // end Namespace tri
} // end Namespace io
} // end Namespace vcg

#endif
