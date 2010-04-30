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

#include <vcg/simplex/face/component.h>
#include <vcg/simplex/face/component_ocf.h>
#include <vcg/simplex/vertex/component.h>
#include <vcg/simplex/vertex/component_ocf.h>

namespace vcg {
namespace tri {
namespace io {

 	template <int N> struct PlaceHolderType{ char A[N];};


	template <class SaveMeshType>
	class ExporterVMI
	{
		
		static void WriteString(FILE *f,const char * in)		{ unsigned int l = strlen(in); fwrite(&l,4,1,f); fwrite(in,1,l,f);} 
		static void  WriteInt(FILE *f,const unsigned int i)	{ fwrite(&i,1,4,f);}  

		static void  WriteFloat(FILE *f,const float v)	{ fwrite(&v,1,sizeof(float),f);}

		/* save Ocf Vertex Components */
		template <typename OpenMeshType,typename CONT>
		struct SaveVertexOcf{
            SaveVertexOcf(FILE*f, const CONT & /*vert*/, bool only_header){
				// do nothing, it is a std::vector
				if(only_header){
					WriteString(f,"NOT_HAS_VERTEX_QUALITY_OCF");
					WriteString(f,"NOT_HAS_VERTEX_COLOR_OCF");
					WriteString(f,"NOT_HAS_VERTEX_NORMAL_OCF");
					WriteString(f,"NOT_HAS_VERTEX_MARK_OCF");
					WriteString(f,"NOT_HAS_VERTEX_TEXCOORD_OCF");
					WriteString(f,"NOT_HAS_VERTEX_VFADJACENCY_OCF");
					WriteString(f,"NOT_HAS_VERTEX_CURVATURE_OCF");
					WriteString(f,"NOT_HAS_VERTEX_CURVATUREDIR_OCF");
					WriteString(f,"NOT_HAS_VERTEX_RADIUS_OCF");
				}
			}
		};

		/* partial specialization for vector_ocf */
		template <typename MeshType>
		struct SaveVertexOcf<MeshType, vertex::vector_ocf<typename MeshType::VertexType> >{
			typedef typename MeshType::VertexType VertexType;
			SaveVertexOcf(FILE * f,const vertex::vector_ocf<VertexType> & vert, bool only_header){

                if( VertexType::HasQualityOcf() && vert.IsQualityEnabled()){
					WriteString(f,"HAS_VERTEX_QUALITY_OCF");
                    if(!only_header) fwrite(&vert.QV[0],sizeof(typename VertexType::QualityType),vert.size(),f);
				}else WriteString(f,"NOT_HAS_VERTEX_QUALITY_OCF");

                if( VertexType::HasColorOcf() && vert.IsColorEnabled()){
					WriteString(f,"HAS_VERTEX_COLOR_OCF");
                    if(!only_header) fwrite(&vert.CV[0],sizeof(typename VertexType::ColorType),vert.size(),f);
				}else WriteString(f,"NOT_HAS_VERTEX_COLOR_OCF");

                if( VertexType::HasNormalOcf() && vert.IsNormalEnabled()){
					WriteString(f,"HAS_VERTEX_NORMAL_OCF");
                    if(!only_header) fwrite(&vert.NV[0],sizeof(typename VertexType::NormalType),vert.size(),f);
				}else WriteString(f,"NOT_HAS_VERTEX_NORMAL_OCF");

                if( VertexType::HasMarkOcf() && vert.IsMarkEnabled()){
					WriteString(f,"HAS_VERTEX_MARK_OCF");
                    if(!only_header) fwrite(&vert.MV[0],sizeof(typename VertexType::MarkType),vert.size(),f);
				}else WriteString(f,"NOT_HAS_VERTEX_MARK_OCF");

				if( VertexType::HasTexCoordOcf() && vert.IsTexCoordEnabled()){
					WriteString(f,"HAS_VERTEX_TEXCOORD_OCF");
                    if(!only_header) fwrite(&vert.TV[0],sizeof(typename VertexType::TexCoordType),vert.size(),f);
				}else WriteString(f,"NOT_HAS_VERTEX_TEXCOORD_OCF");

				if( VertexType::HasVFAdjacencyOcf() && vert.IsVFAdjacencyEnabled()){
					WriteString(f,"HAS_VERTEX_VFADJACENCY_OCF");
                    if(!only_header) fwrite(&vert.AV[0],sizeof(typename vertex::vector_ocf<VertexType>::VFAdjType),vert.size(),f);
				}else WriteString(f,"NOT_HAS_VERTEX_VFADJACENCY_OCF");

				if( VertexType::HasCurvatureOcf() && vert.IsCurvatureEnabled()){
					WriteString(f,"HAS_VERTEX_CURVATURE_OCF");
                    if(!only_header) fwrite(&vert.CuV[0],sizeof(typename VertexType::CurvatureType),vert.size(),f);
				}else WriteString(f,"NOT_HAS_VERTEX_CURVATURE_OCF");

				if( VertexType::HasCurvatureDirOcf() && vert.IsCurvatureDirEnabled()){
					WriteString(f,"HAS_VERTEX_CURVATUREDIR_OCF");
                    if(!only_header) fwrite(&vert.CuDV[0],sizeof(typename VertexType::CurvatureDirType),vert.size(),f);
				}else WriteString(f,"NOT_HAS_VERTEX_CURVATUREDIR_OCF");

				if( VertexType::HasRadiusOcf() && vert.IsRadiusEnabled()){
					WriteString(f,"HAS_VERTEX_RADIUS_OCF");
                    if(!only_header) fwrite(&vert.RadiusV[0],sizeof(typename VertexType::RadiusType),vert.size(),f);
				}else WriteString(f,"NOT_HAS_VERTEX_RADIUS_OCF");

			}
		};


		/* save Ocf Face Components */
		template <typename MeshType,typename CONT>
		struct SaveFaceOcf{
            SaveFaceOcf(FILE * f,const CONT & /*face*/, bool only_header){
				// it is a std::vector
				if(only_header){
					WriteString(f,"NOT_HAS_FACE_QUALITY_OCF");
					WriteString(f,"NOT_HAS_FACE_COLOR_OCF");
					WriteString(f,"NOT_HAS_FACE_NORMAL_OCF");
					WriteString(f,"NOT_HAS_FACE_MARK_OCF");
					WriteString(f,"NOT_HAS_FACE_WEDGETEXCOORD_OCF");
					WriteString(f,"NOT_HAS_FACE_FFADJACENCY_OCF");
					WriteString(f,"NOT_HAS_FACE_VFADJACENCY_OCF");
					WriteString(f,"NOT_HAS_FACE_WEDGECOLOR_OCF");
					WriteString(f,"NOT_HAS_FACE_WEDGENORMAL_OCF");
				}
			}
		};

		/* partial specialization for vector_ocf */
		template <typename MeshType>
		struct SaveFaceOcf<  MeshType, face::vector_ocf<typename MeshType::FaceType> >{
			typedef typename MeshType::FaceType FaceType;
			SaveFaceOcf(FILE * f,const face::vector_ocf<FaceType> & face, bool only_header){

				if( FaceType::HasFaceQualityOcf() && face.IsQualityEnabled()){
					WriteString(f,"HAS_FACE_QUALITY_OCF");
                    if(!only_header) fwrite(&face.QV[0],sizeof(typename FaceType::QualityType),face.size(),f);
				}else WriteString(f,"NOT_HAS_FACE_QUALITY_OCF");

				if( FaceType::HasFaceColorOcf() && face.IsColorEnabled()){
					WriteString(f,"HAS_FACE_COLOR_OCF");
                    if(!only_header) fwrite(&face.CV[0],sizeof(typename FaceType::ColorType),face.size(),f);
				}else WriteString(f,"NOT_HAS_FACE_COLOR_OCF");

				if( FaceType::HasFaceNormalOcf() && face.IsNormalEnabled()){
					WriteString(f,"HAS_FACE_NORMAL_OCF");
                    if(!only_header) fwrite(&face.NV[0],sizeof(typename FaceType::NormalType),face.size(),f);
				}else WriteString(f,"NOT_HAS_FACE_NORMAL_OCF");

				if( FaceType::HasFaceMarkOcf() && face.IsMarkEnabled()){
					WriteString(f,"HAS_FACE_MARK_OCF");
                    if(!only_header) fwrite(&face.MV[0],sizeof(typename FaceType::MarkType),face.size(),f);
				}else WriteString(f,"NOT_HAS_FACE_MARK_OCF");

				if( FaceType::HasWedgeTexCoordOcf() && face.IsWedgeTexEnabled()){
					WriteString(f,"HAS_FACE_WEDGETEXCOORD_OCF");
                    if(!only_header) fwrite(&face.WTV[0],sizeof(typename FaceType::WedgeTexCoordType),face.size(),f);
				}else WriteString(f,"NOT_HAS_FACE_WEDGETEXCOORD_OCF");

				if( FaceType::HasFFAdjacencyOcf() && face.IsFFAdjacencyEnabled()){
					WriteString(f,"HAS_FACE_FFADJACENCY_OCF");
                    if(!only_header) fwrite(&face.AF[0],sizeof(typename face::vector_ocf<FaceType>::AdjTypePack),face.size(),f);
				}else WriteString(f,"NOT_HAS_FACE_FFADJACENCY_OCF");

				if( FaceType::HasVFAdjacencyOcf() && face.IsVFAdjacencyEnabled()){
					WriteString(f,"HAS_FACE_VFADJACENCY_OCF");
                    if(!only_header) fwrite(&face.AV[0],sizeof(typename face::vector_ocf<FaceType>::AdjTypePack),face.size(),f);
				}else WriteString(f,"NOT_HAS_FACE_VFADJACENCY_OCF");

				if( FaceType::HasWedgeColorOcf() && face.IsWedgeColorEnabled()){
					WriteString(f,"HAS_FACE_WEDGECOLOR_OCF");
                    if(!only_header) fwrite(&face.WCV[0],sizeof(typename face::vector_ocf<FaceType>::WedgeColorTypePack),face.size(),f);
				}else WriteString(f,"NOT_HAS_FACE_WEDGECOLOR_OCF");

				if( FaceType::HasWedgeNormalOcf() && face.IsWedgeNormalEnabled()){
					WriteString(f,"HAS_FACE_WEDGENORMAL_OCF");
                    if(!only_header) fwrite(&face.WNV[0],sizeof(typename face::vector_ocf<FaceType>::WedgeNormalTypePack),face.size(),f);
				}else WriteString(f,"NOT_HAS_FACE_WEDGENORMAL_OCF");
			}
		};



		static FILE *& F(){static FILE * f; return f;}

		typedef typename SaveMeshType::FaceContainer FaceContainer;
		typedef typename SaveMeshType::FaceIterator FaceIterator;
		typedef typename SaveMeshType::VertContainer VertContainer;
		typedef typename SaveMeshType::VertexIterator VertexIterator;
		typedef typename SaveMeshType::VertexType VertexType;
		typedef typename SaveMeshType::FaceType FaceType;
		typedef SimpleTempDataBase<typename SaveMeshType::VertContainer> STDBv;
		typedef SimpleTempDataBase<typename SaveMeshType::FaceContainer> STDBf;
	//	typedef typename SaveMeshType::Attribute <SaveMeshType::FaceContainer> STDBm;
		
		/* save Ocf Components */ 


	public:
        static int Save(const SaveMeshType &m,const char * filename){
			unsigned int i;
			unsigned int vertSize,faceSize;
			F() = fopen(filename,"wb");
            if(F()==NULL)	return 1; // 1 is the error code for cant'open, see the ErrorMsg function
            std::vector<std::string> nameF,nameV;
			SaveMeshType::FaceType::Name(nameF);
			SaveMeshType::VertexType::Name(nameV);
			vertSize = m.vert.size();
			faceSize = m.face.size();

			/* write header */
			WriteString(F(),"FACE_TYPE");
			WriteInt(F(),nameF.size());

			for(i=0; i < nameF.size(); ++i) WriteString(F(),nameF[i].c_str());
			SaveFaceOcf<SaveMeshType,FaceContainer>(F(),m.face,true);
			WriteString(F(),"SIZE_VECTOR_FACES");
			WriteInt(F(), faceSize );

			WriteString(F(),"VERTEX_TYPE");
			WriteInt(F(),nameV.size());

			for(i=0; i < nameV.size(); ++i) WriteString(F(),nameV[i].c_str());
			SaveVertexOcf<SaveMeshType,VertContainer>(F(),m.vert,true);

			WriteString(F(),"SIZE_VECTOR_VERTS");
			WriteInt(F(),vertSize);

			WriteString(F(),"BOUNDING_BOX");
			float float_value;
			for(unsigned int i =0; i < 2; ++i){float_value = m.bbox.min[i]; WriteFloat(F(),float_value);}
			for(unsigned int i =0; i < 2; ++i){float_value = m.bbox.max[i]; WriteFloat(F(),float_value);}

			WriteString(F(),"end_header");
			/* end header */

			if(vertSize!=0){
                                size_t offsetV = (size_t) &m.vert[0];
				/* write the address of the first vertex */
                                fwrite(&offsetV,sizeof(size_t),1,F());
			}

			if(faceSize!=0){
                                 size_t offsetF= ( size_t) &m.face[0];
				/* write the address of the first face */
                                fwrite(&offsetF,sizeof( size_t),1,F());
			}
			/* save the object mesh */
			fwrite(&m.shot,sizeof(Shot<typename SaveMeshType::ScalarType>),1,F());
			fwrite(&m.vn,sizeof(int),1,F());
			fwrite(&m.fn,sizeof(int),1,F());
			fwrite(&m.imark,sizeof(int),1,F());
			fwrite(&m.bbox,sizeof(Box3<typename SaveMeshType::ScalarType>),1,F());
			fwrite(&m.C(),sizeof(Color4b),1,F());

			unsigned int written;


			if(vertSize!=0){
				/* save the vertices */
				written = fwrite((void*)&m.vert[0],sizeof(typename SaveMeshType::VertexType),m.vert.size(),F());
				assert(written==m.vert.size());
				SaveVertexOcf<SaveMeshType,VertContainer>(F(),m.vert,false);
			}

			if(faceSize!=0){
				/* save the faces */
				written = fwrite((void*)&m.face[0],sizeof(typename SaveMeshType::FaceType),faceSize,F());
				assert(written==m.face.size());

				SaveFaceOcf<SaveMeshType,FaceContainer>(F(),m.face,false);

			}


 
			/* save the attributes */
            typename std::set< typename SaveMeshType::PointerToAttribute>::const_iterator ai;
 

			/* save the per vertex attributes */
			{
				typename std::set< typename SaveMeshType::PointerToAttribute>::const_iterator ai;
				unsigned int n_named_attr = 0;
				for(ai = m.vert_attr.begin(); ai != m.vert_attr.end(); ++ai) n_named_attr+=!(*ai)._name.empty();

				WriteString(F(),"N_PER_VERTEX_ATTRIBUTES"); WriteInt (F(),n_named_attr);
				for(ai = m.vert_attr.begin(); ai != m.vert_attr.end(); ++ai)
					if(!(*ai)._name.empty())
						{
							STDBv * stdb = (STDBv *) (*ai)._handle;

							WriteString(F(),"PER_VERTEX_ATTR_NAME"); 
							WriteString(F(),(*ai)._name.c_str() ); 

							WriteString(F(),"PER_VERTEX_ATTR_SIZE");  
							WriteInt(F(),stdb->SizeOf());

							fwrite(stdb->DataBegin(),m.vert.size(),stdb->SizeOf(),F());
 						}
			}

			/* save the per face attributes */
			{
				typename std::set< typename SaveMeshType::PointerToAttribute>::const_iterator ai;
				unsigned int n_named_attr = 0;
				for(ai = m.face_attr.begin(); ai != m.face_attr.end(); ++ai) n_named_attr+=!(*ai)._name.empty();

				WriteString(F(),"N_PER_FACE_ATTRIBUTES");
				WriteInt (F(),n_named_attr);

				for(ai = m.face_attr.begin(); ai != m.face_attr.end(); ++ai)
					if(!(*ai)._name.empty())
						{
							STDBf * stdb = (STDBf *) (*ai)._handle;

							WriteString(F(),"PER_FACE_ATTR_NAME");
							WriteString(F(),(*ai)._name.c_str());

							WriteString(F(),"PER_FACE_ATTR_SIZE");
							WriteInt(F(),stdb->SizeOf());

							fwrite(stdb->DataBegin(),m.face.size(),stdb->SizeOf(),F());
 						}
			}

			///* save the per mesh attributes */
			{
				typename std::set< typename SaveMeshType::PointerToAttribute>::const_iterator ai;
				unsigned int n_named_attr = 0;
				for(ai = m.mesh_attr.begin(); ai != m.mesh_attr.end(); ++ai) n_named_attr+=!(*ai)._name.empty();
				WriteString(F(),"N_PER_MESH_ATTRIBUTES"); WriteInt(F(),n_named_attr);			 
				for(ai = m.mesh_attr.begin(); ai != m.mesh_attr.end(); ++ai)
					if(!(*ai)._name.empty())
						{
							AttributeBase  *    handle =  (AttributeBase  *)   (*ai)._handle ;
	
							WriteString(F(),"PER_MESH_ATTR_NAME");
							WriteString(F(),(*ai)._name.c_str());

							WriteString(F(),"PER_MESH_ATTR_SIZE");
							WriteInt(F(),handle->SizeOf());

							fwrite(handle->DataBegin(),1,handle->SizeOf(),F());
 						}
			}

			//	fflush(F());
			fclose(F());
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
	}; // end class

} // end Namespace tri
} // end Namespace io
} // end Namespace vcg

#endif
