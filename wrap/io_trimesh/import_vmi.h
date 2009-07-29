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

	/* derivation chain */
	template <int N> struct DummyType{ char placeholder[N]; };

	template <class MeshType, class A, class T>
	struct Der:public T{
		typedef typename std::set<typename MeshType::PointerToAttribute >::iterator HWIte;
		static void AddAttrib(MeshType &m, char * name, int s, void * data){
			if(s == sizeof(A)){
				 typename MeshType::template PerVertexAttributeHandle<A> h = vcg::tri::Allocator<MeshType>:: template AddPerVertexAttribute<A>(m,name);
				for(int i  = 0; i < m.vert.size(); ++i)
					memcpy(&h[i], (void*) &((A*)data)[i],sizeof(A)); // we don't want the type conversion
			}
			else
				T::AddAttrib(m,name,s,data);
		}
	};
	template <class MeshType, class A, class T>
	struct DerK:public T{
		typedef typename std::set<typename MeshType::PointerToAttribute >::iterator HWIte;
		static void AddAttrib(MeshType &m, char * name, int s, void * data){
			if(s == sizeof(A)){
				typename MeshType::template PerVertexAttributeHandle<A> h = vcg::tri::Allocator<MeshType>::template AddPerVertexAttribute<A>(m,name);
				for(unsigned int i  = 0; i < m.vert.size(); ++i)
					memcpy((void*) &(h[i]), (void*) &((A*)data)[i],sizeof(A)); // we don't want the type conversion
			}
			else
				if(s < sizeof(A)){
					// padding
					int padd = sizeof(A) - s;
					typename MeshType::template PerVertexAttributeHandle<A> h = vcg::tri::Allocator<MeshType>::template AddPerVertexAttribute<A>(m,name);
					for(unsigned int i  = 0; i < m.vert.size(); ++i){
						char * dest =  &((char*)(&h[i]))[padd];
						memcpy( (void *)dest , (void*) &((A*)data)[i],s); // we don't want the type conversion
					}
					typename MeshType::PointerToAttribute pa;
					pa._name = std::string(name);
					HWIte res = m.vert_attr.find(pa);
					pa = *res;
					m.vert_attr.erase(res);
					pa._padding = padd;
					std::pair<HWIte,bool > new_pa = m.vert_attr.insert(pa);
					assert(new_pa.second);
				}
				else
 					T::AddAttrib(m,name,s,data);
		}
	};


	template <class MeshType>	struct K	{
		static void AddAttrib(MeshType &m, char * name, int s, void * data){

			assert(0);		
		}
	};

	template <class MeshType, class B0 >																												struct K0	: public DerK<  MeshType, B0,	K<MeshType> > {};
	template <class MeshType, class B0, class B1 >																										struct K1	: public DerK<  MeshType, B1,	K0<MeshType, B0> > {};
	template <class MeshType, class B0, class B1, class B2 >																							struct K2	: public DerK<  MeshType, B2,	K1<MeshType, B0, B1> > {};
	template <class MeshType, class B0, class B1, class B2,class B3>																					struct K3	: public DerK<  MeshType, B3,	K2<MeshType, B0, B1, B2> > {};
	template <class MeshType, class B0, class B1, class B2,class B3,class B4>																			struct K4	: public DerK<  MeshType, B4,	K3<MeshType, B0, B1, B2, B3> > {};
	template <class MeshType, class B0, class B1, class B2,class B3,class B4,class B5>																	struct K5	: public DerK<  MeshType, B5,	K4<MeshType, B0, B1, B2, B3, B4> > {};
	template <class MeshType, class B0, class B1, class B2,class B3,class B4,class B5,class B6>															struct K6	: public DerK<  MeshType, B6,	K5<MeshType, B0, B1, B2, B3, B4, B5> > {};
	template <class MeshType, class B0, class B1, class B2,class B3,class B4,class B5,class B6,class B7>												struct K7	: public DerK<  MeshType, B7,	K6<MeshType, B0, B1, B2, B3, B4, B5, B6> > {};
	template <class MeshType, class B0, class B1, class B2,class B3,class B4,class B5,class B6,class B7,class B8>										struct K8	: public DerK<  MeshType, B8,	K7<MeshType, B0, B1, B2, B3, B4, B5, B6, B7> > {};
	template <class MeshType, class B0, class B1, class B2,class B3,class B4,class B5,class B6,class B7,class B8,class B9>								struct K9	: public DerK<  MeshType, B9,	K8<MeshType, B0, B1, B2, B3, B4, B5, B6, B7, B8> > {};
	template <class MeshType, class B0, class B1, class B2,class B3,class B4,class B5,class B6,class B7,class B8,class B9,class B10>					struct K10	: public DerK<  MeshType, B10,	K9<MeshType, B0, B1, B2, B3, B4, B5, B6, B7, B8, B9> > {};
	template <class MeshType, class B0, class B1, class B2,class B3,class B4,class B5,class B6,class B7,class B8,class B9,class B10,class B11>			struct K11	: public DerK<  MeshType, B11,	K10<MeshType, B0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B11 > > {};
	template <class MeshType, class B0, class B1, class B2,class B3,class B4,class B5,class B6,class B7,class B8,class B9,class B10,class B11,class B12>struct K12	: public DerK<  MeshType, B12,	K11<MeshType, B0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B11, B12 > > {};

	template <class MeshType, class A0, 
		class B0  = DummyType<1048576>,
		class B1  = DummyType<2048>,
		class B2  = DummyType<1024>,
		class B3  = DummyType<512>,
		class B4  = DummyType<256>, 
		class B5  = DummyType<128>,
		class B6  = DummyType<64>,
		class B7  = DummyType<32>, 
		class B8  = DummyType<16>, 
		class B9  = DummyType<8>, 
		class B10 = DummyType<4>, 
		class B11 = DummyType<2>, 
		class B12 = DummyType<1> 
	>	struct C0		: public DerK<  MeshType, A0,    K12<MeshType, B0, B1, B2, B3, B4,B5,B6,B7,B8,B9,B10,B11,B12> > {};

	template <class MeshType, class A0, class A1>											struct C1		: public Der<  MeshType, A1,	C0<MeshType, A0> > {};
	template <class MeshType, class A0, class A1, class A2>	 								struct C2		: public Der<  MeshType, A2,	C1<MeshType, A0, A1> > {};
	template <class MeshType, class A0, class A1, class A2,class A3>	 					struct C3		: public Der<  MeshType, A3,	C2<MeshType, A0, A1, A2> > {};
	template <class MeshType, class A0, class A1, class A2,class A3,class A4>				struct AttrAll	: public Der<  MeshType, A4,	C3<MeshType, A0, A1, A2, A3> > {};


	/* end derivation chain */
 	 
	template <class OpenMeshType,class A0 = long, class A1 = double, class A2 = int,class A3 = short, class A4 = char > 
	class ImporterVMI: public AttrAll<OpenMeshType,A0,A1,A2,A3,A4>
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

			 if(vertSize!=0)
				/* read the address of the first vertex */
				fread(&offsetV,sizeof( int),1,f);

			 if(faceSize!=0)
				/* read the address of the first face */
				fread(&offsetF,sizeof( int),1,f);

			/* read the object mesh */
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
		

			/* load the per vertex attributes */
			char _string[65536],_trash[65536];
			int n,sz;
	
			fscanf(f,"%s %d",&_trash[0],&n);
			for(int ia = 0 ; ia < n; ++ia){
				fscanf(f,"%s %s",&_trash[0],&_string[0]);
				fscanf(f,"%s %d",&_trash[0],&sz);
				void * data = malloc(sz*m.vert.size());
				fread(data,sz,m.vert.size(),f);
				AddAttrib(m,_string,sz,data);
				free(data);
			}

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
