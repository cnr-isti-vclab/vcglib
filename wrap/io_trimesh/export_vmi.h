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

 	template <int N> struct PlaceHolderType{ char A[N];};

	template <class SaveMeshType>
	class ExporterVMI
	{
	public:	
		struct WriteString{	WriteString(FILE * f,const char * in){ unsigned int l = strlen(in); fwrite(&l,4,1,f); fwrite(in,1,l,f);}};
		struct WriteInt{ WriteInt(FILE *f, unsigned int i){ fwrite(&i,1,4,f);}};

		typedef typename SaveMeshType::FaceIterator FaceIterator;
		typedef typename SaveMeshType::VertexIterator VertexIterator;
		typedef typename SaveMeshType::VertexType VertexType;
		typedef SimpleTempDataBase<typename SaveMeshType::VertContainer> STDBv;
		typedef SimpleTempDataBase<typename SaveMeshType::FaceContainer> STDBf;
	//	typedef typename SaveMeshType::Attribute <SaveMeshType::FaceContainer> STDBm;

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
			WriteString(f,"FACE_TYPE");
			WriteInt(f,nameF.size());

			for(i=0; i < nameF.size(); ++i) WriteString(f,nameF[i].c_str());
			WriteString(f,"SIZE_VECTOR_FACES");
			WriteInt(f,faceSize);

			WriteString(f,"VERTEX_TYPE");
			WriteInt(f,nameV.size());

			for(i=0; i < nameV.size(); ++i) WriteString(f,nameV[i].c_str());

			WriteString(f,"SIZE_VECTOR_VERTS");
			WriteInt(f,vertSize);

			WriteString(f,"end_header");

			if(vertSize!=0){
				unsigned int offsetV = (unsigned int) &m.vert[0];
				/* write the address of the first vertex */
				fwrite(&offsetV,sizeof(unsigned int),1,f);
			}

			if(faceSize!=0){
				 int offsetF= ( int) &m.face[0];
				/* write the address of the first face */
				fwrite(&offsetF,sizeof( int),1,f);
			}
			/* save the object mesh */
			fwrite(&m.shot,sizeof(Shot<typename SaveMeshType::ScalarType>),1,f);
			fwrite(&m.vn,sizeof(int),1,f);
			fwrite(&m.fn,sizeof(int),1,f);
			fwrite(&m.imark,sizeof(int),1,f);
			fwrite(&m.bbox,sizeof(Box3<typename SaveMeshType::ScalarType>),1,f);
			fwrite(&m.C(),sizeof(Color4b),1,f);

			unsigned int written;

			if(vertSize!=0){
				/* save the vertices */
				written = fwrite((void*)&m.vert[0],sizeof(typename SaveMeshType::VertexType),m.vert.size(),f);
				assert(written==m.vert.size());
			}

			if(faceSize!=0){
				/* save the faces */
				written = fwrite((void*)&m.face[0],sizeof(typename SaveMeshType::FaceType),m.face.size(),f);
				assert(written==m.face.size());
			}

 
			/* save the attribtues */
            typename std::set< typename SaveMeshType::PointerToAttribute>::const_iterator ai;
 

			/* save the per vertex attributes */
			{
				typename std::set< typename SaveMeshType::PointerToAttribute>::const_iterator ai;
				unsigned int n_named_attr = 0;
				for(ai = m.vert_attr.begin(); ai != m.vert_attr.end(); ++ai) n_named_attr+=!(*ai)._name.empty();

				WriteString(f,"N_PER_VERTEX_ATTRIBUTES"); WriteInt(f,n_named_attr);
				for(ai = m.vert_attr.begin(); ai != m.vert_attr.end(); ++ai)
					if(!(*ai)._name.empty())
						{
							STDBv * stdb = (STDBv *) (*ai)._handle;

							WriteString(f,"PER_VERTEX_ATTR_NAME"); 
							WriteString(f,(*ai)._name.c_str() ); 

							WriteString(f,"PER_VERTEX_ATTR_SIZE");  
							WriteInt(f,stdb->SizeOf());

							fwrite(stdb->DataBegin(),m.vert.size(),stdb->SizeOf(),f);
 						}
			}

			/* save the per face attributes */
			{
				typename std::set< typename SaveMeshType::PointerToAttribute>::const_iterator ai;
				unsigned int n_named_attr = 0;
				for(ai = m.face_attr.begin(); ai != m.face_attr.end(); ++ai) n_named_attr+=!(*ai)._name.empty();

				WriteString(f,"N_PER_FACE_ATTRIBUTES");
				WriteInt(f,n_named_attr);

				for(ai = m.face_attr.begin(); ai != m.face_attr.end(); ++ai)
					if(!(*ai)._name.empty())
						{
							STDBf * stdb = (STDBf *) (*ai)._handle;

							WriteString(f,"PER_FACE_ATTR_NAME");
							WriteString(f,(*ai)._name.c_str());

							WriteString(f,"PER_FACE_ATTR_SIZE");
							WriteInt(f,stdb->SizeOf());

							fwrite(stdb->DataBegin(),m.face.size(),stdb->SizeOf(),f);
 						}
			}

			///* save the per mesh attributes */
			{
				typename std::set< typename SaveMeshType::PointerToAttribute>::const_iterator ai;
				unsigned int n_named_attr = 0;
				for(ai = m.mesh_attr.begin(); ai != m.mesh_attr.end(); ++ai) n_named_attr+=!(*ai)._name.empty();
				WriteString(f,"N_PER_MESH_ATTRIBUTES"); WriteInt(f,n_named_attr);			 
				for(ai = m.mesh_attr.begin(); ai != m.mesh_attr.end(); ++ai)
					if(!(*ai)._name.empty())
						{
							AttributeBase  *    handle =  (AttributeBase  *)   (*ai)._handle ;
	
							WriteString(f,"PER_MESH_ATTR_NAME");
							WriteString(f,(*ai)._name.c_str());

							WriteString(f,"PER_MESH_ATTR_SIZE");
							WriteInt(f,handle->SizeOf());

							fwrite(handle->DataBegin(),1,handle->SizeOf(),f);
 						}
			}

			//	fflush(f);
			fclose(f);
		}

	}; // end class

} // end Namespace tri
} // end Namespace io
} // end Namespace vcg

#endif
