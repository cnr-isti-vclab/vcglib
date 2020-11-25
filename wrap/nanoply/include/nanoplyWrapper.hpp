/****************************************************************************
* NanoPLY                                                                   *
* NanoPLY is a C++11 header-only library to read and write PLY file         *
*                                                                           *
* Copyright(C) 2014-2015                                                    *
* Visual Computing Lab                                                      *
* ISTI - Italian National Research Council                                  *
*                                                                           *
* This Source Code Form is subject to the terms of the Mozilla Public       *
* License, v. 2.0. If a copy of the MPL was not distributed with this       *
* file, You can obtain one at http://mozilla.org/MPL/2.0/.                  *
*                                                                           *
****************************************************************************/

#ifndef NANOPLY_WRAPPER_VCG_H
#define NANOPLY_WRAPPER_VCG_H

#include <wrap/nanoply/include/nanoply.hpp>
#include <vcg/complex/complex.h>
#include <vcg/space/point.h>
#include <map>

#include  <typeinfo>

namespace nanoply
{

	template <class MeshType>
	class NanoPlyWrapper{

				
	private:

    typedef typename MeshType::PointerToAttribute                         PointerToAttribute;
    typedef typename MeshType::ScalarType                                 ScalarType;

    typedef typename MeshType::VertexType                                 VertexType;
    typedef typename MeshType::VertexType::ScalarType                     VertexCoordScalar;
    typedef typename MeshType::VertexType::NormalType                     VertexNormalType;
    typedef typename MeshType::VertexType::NormalType::ScalarType         VertexNormScalar;
    typedef typename MeshType::VertexType::ColorType                      VertexColorType;
    typedef typename MeshType::VertexType::ColorType::ScalarType          VertexColorScalar;
    typedef typename MeshType::VertexType::QualityType                    VertexQuality;
    typedef typename MeshType::VertexType::RadiusType                     VertexRadius;
    typedef typename MeshType::VertexType::FlagType                       VertexFlag;
    typedef typename MeshType::VertexType::TexCoordType                   VertexTexCoordType;
    typedef typename MeshType::VertexType::TexCoordType::ScalarType       VertexTexScalar;
    typedef typename MeshType::VertexType::CurvatureType                  VertexCurType;
    typedef typename MeshType::VertexType::CurvatureType::ScalarType      VertexCurScalar;
    typedef typename MeshType::VertexType::CurvatureDirType               VertexCurDirType;
    typedef typename MeshType::VertexType::CurScalarType                  VertexDirCurScalar;
    typedef typename MeshType::VertexType::CurVecType::ScalarType         VertexDirCurVecScalar;

    typedef typename MeshType::EdgeType                                   EdgeType;
    typedef typename MeshType::EdgeType::ColorType::ScalarType            EdgeColorScalar;
    typedef typename MeshType::EdgeType::QualityType                      EdgeQuality;
    typedef typename MeshType::EdgeType::FlagType                         EdgeFlag;

    typedef typename MeshType::FaceType                                   FaceType;
    typedef typename MeshType::FaceType::NormalType                       FaceNormalType;
    typedef typename MeshType::FaceType::NormalType::ScalarType           FaceNormScalar;
    typedef typename MeshType::FaceType::ColorType                        FaceColorType;
    typedef typename MeshType::FaceType::ColorType::ScalarType            FaceColorScalar;
    typedef typename MeshType::FaceType::QualityType                      FaceQuality;
    typedef typename MeshType::FaceType::FlagType                         FaceFlag;
    typedef typename vcg::face::vector_ocf<FaceType>::WedgeTexTypePack    FaceTexCoordType;
    typedef typename MeshType::FaceType::TexCoordType::ScalarType         FaceTexScalar;
    typedef typename MeshType::FaceType::CurvatureDirType                 FaceCurDirType;
    typedef typename MeshType::FaceType::CurScalarType                    FaceDirCurScalar;
    typedef typename MeshType::FaceType::CurVecType::ScalarType           FaceDirCurVecScalar;
    typedef typename vcg::face::vector_ocf<FaceType>::WedgeColorTypePack  WedgeColorType;
    typedef typename MeshType::FaceType::WedgeColorType::ScalarType       WedgeColorScalar;
    typedef typename vcg::face::vector_ocf<FaceType>::WedgeNormalTypePack WedgeNormalType;
    typedef typename MeshType::FaceType::WedgeNormalType::ScalarType      WedgeNormalScalar;

    typedef typename MeshType::FaceIterator                               FaceIterator;

		template<class T> static PlyType getEntity() { return NNP_UNKNOWN_TYPE; };
		template<> static PlyType getEntity<unsigned char>(){ return NNP_UINT8; };
		template<> static PlyType getEntity<char>(){ return NNP_INT8; };
		template<> static PlyType getEntity<unsigned short>(){ return NNP_UINT16; };
		template<> static PlyType getEntity<short>(){ return NNP_INT16; };
		template<> static PlyType getEntity<unsigned int>(){ return NNP_UINT32; };
		template<> static PlyType getEntity<int>(){ return NNP_INT32; };
		template<> static PlyType getEntity<float>(){ return NNP_FLOAT32; };
		template<> static PlyType getEntity<double>(){ return NNP_FLOAT64; };

    template<class T> static PlyType getEntityList() { return NNP_UNKNOWN_TYPE; };
		template<> static PlyType getEntityList<unsigned char>(){ return NNP_LIST_UINT8_UINT8; };
		template<> static PlyType getEntityList<char>(){ return NNP_LIST_UINT8_INT8; };
		template<> static PlyType getEntityList<unsigned short>(){ return NNP_LIST_UINT8_UINT16; };
		template<> static PlyType getEntityList<short>(){ return NNP_LIST_UINT8_INT16; };
		template<> static PlyType getEntityList<unsigned int>(){ return NNP_LIST_UINT8_UINT32; };
		template<> static PlyType getEntityList<int>(){ return NNP_LIST_UINT8_INT32; };
		template<> static PlyType getEntityList<float>(){ return NNP_LIST_UINT8_FLOAT32; };
		template<> static PlyType getEntityList<double>(){ return NNP_LIST_UINT8_FLOAT64; };


		template<class Container, class Type, int n>
		inline static void PushDescriport(std::vector<PlyProperty>& prop, ElementDescriptor& elem, PlyEntity entity, void* ptr)
		{
			prop.push_back(PlyProperty(getEntity<Type>(), entity));
			DescriptorInterface* di = new DataDescriptor<Container, n, Type>(entity, ptr);
			elem.dataDescriptor.push_back(di);
		}


		template<class Container, class Type, int n>
		inline static void PushDescriportList(std::vector<PlyProperty>& prop, ElementDescriptor& elem, PlyEntity entity, void* ptr)
		{
			prop.push_back(PlyProperty(getEntityList<Type>(), entity));
			DescriptorInterface* di = new DataDescriptor<Container, n, Type>(entity, ptr);
			elem.dataDescriptor.push_back(di);
		}


		template<class Container, class Type, int n>
		inline static void PushDescriport(std::vector<PlyProperty>& prop, ElementDescriptor& elem, const std::string& name, void* ptr)
		{
			prop.push_back(PlyProperty(getEntity<Type>(), name));
			DescriptorInterface* di = new DataDescriptor<Container, n, Type>(name, ptr);
			elem.dataDescriptor.push_back(di);
		}


		template<class Container, class Type, int n>
		inline static void PushDescriportList(std::vector<PlyProperty>& prop, ElementDescriptor& elem, const std::string& name, void* ptr)
		{
			prop.push_back(PlyProperty(getEntityList<Type>(), name));
			DescriptorInterface* di = new DataDescriptor<Container, n, Type>(name, ptr);
			elem.dataDescriptor.push_back(di);
		}
	
	public:


    typedef enum {
      IO_NONE = 0x00000000,

      IO_VERTCOORD = 0x00000001,
      IO_VERTFLAGS = 0x00000002,
      IO_VERTCOLOR = 0x00000004,
      IO_VERTQUALITY = 0x00000008,
      IO_VERTNORMAL = 0x00000010,
      IO_VERTTEXCOORD = 0x00000020,
      IO_VERTRADIUS = 0x00000040,
      IO_VERTCURV = 0x00000080,
      IO_VERTCURVDIR = 0x00000100,
      IO_VERTATTRIB = 0x00000200,

      IO_FACEINDEX = 0x00000400,
      IO_FACEFLAGS = 0x00000800,
      IO_FACECOLOR = 0x00001000,
      IO_FACEQUALITY = 0x00002000,
      IO_FACENORMAL = 0x00004000,
      IO_FACECURVDIR = 0x00008000,
      IO_FACEATTRIB = 0x00010000,

      IO_EDGEINDEX = 0x00020000,
      IO_EDGEQUALITY = 0x00040000,
      IO_EDGECOLOR = 0x00080000,
      IO_EDGEFLAGS = 0x00100000,
      IO_EDGEATTRIB = 0x00200000,

      IO_WEDGCOLOR = 0x00400000,
      IO_WEDGTEXCOORD = 0x00800000,
      IO_WEDGTEXMULTI = 0x01000000, // when texture index is explicit
      IO_WEDGNORMAL = 0x02000000,

      IO_ALL_ATTRIB = 0x03FFFFFF,

			IO_BITPOLYGONAL = 0x04000000, // loads explicit polygonal mesh

			IO_CAMERA = 0x08000000,
			IO_MESHATTRIB = 0x10000000,

			IO_FLAGS = IO_VERTFLAGS | IO_FACEFLAGS,

			IO_ALL = 0xFFFFFFFF
		}BitMask;



		class CustomAttributeDescriptor
		{
		public:

			typedef std::map<std::string, ElementDescriptor::PropertyDescriptor> MapMeshAttrib;
			typedef std::map<std::string, ElementDescriptor::PropertyDescriptor>::iterator MapMeshAttribIter;
			typedef std::map<std::string, std::vector<PlyProperty>> MapMeshAttribProp;
			typedef std::map<std::string, std::vector<PlyProperty>>::iterator MapMeshAttribPropIter;
			
			ElementDescriptor::PropertyDescriptor vertexAttrib;
			ElementDescriptor::PropertyDescriptor faceAttrib;
			ElementDescriptor::PropertyDescriptor edgeAttrib;
			std::vector<PlyProperty> vertexAttribProp;
			std::vector<PlyProperty> faceAttribProp;
			std::vector<PlyProperty> edgeAttribProp;

			MapMeshAttrib meshAttrib;
			MapMeshAttribProp meshAttribProp;
			std::map<std::string, int> meshAttribCnt;


			~CustomAttributeDescriptor()
			{
				for (size_t i = 0; i < vertexAttrib.size(); i++)
					delete vertexAttrib[i];
				for (size_t i = 0; i < edgeAttrib.size(); i++)
					delete edgeAttrib[i];
				for (size_t i = 0; i < faceAttrib.size(); i++)
					delete faceAttrib[i];
				CustomAttributeDescriptor::MapMeshAttribIter iter = meshAttrib.begin();
				for (; iter != meshAttrib.end(); iter++)
					for (size_t i = 0; i < (*iter).second.size(); i++)
						delete (*iter).second[i];
			}


			template<class Container, class Type, int n>
			void AddVertexAttribDescriptor(const std::string& name, PlyType type, void* ptr)
			{
        AddAttribDescriptor<Container, Type, n>(name, type, ptr, vertexAttrib, vertexAttribProp);
			}

      template<class Container, class Type, int n>
			void AddEdgeAttribDescriptor(const std::string& name, PlyType type, void* ptr)
			{
        AddAttribDescriptor<Container, Type, n>(name, type, ptr, edgeAttrib, edgeAttribProp);
			}

			template<class Container, class Type, int n>
			void AddFaceAttribDescriptor(const std::string& name, PlyType type, void* ptr)
			{
        AddAttribDescriptor<Container, Type, n>(name, type, ptr, faceAttrib, faceAttribProp);
			}

			template<class Container, class Type, int n>
			void AddMeshAttribDescriptor(const std::string& nameAttrib, const std::string& nameProp, PlyType type, void* ptr)
			{
				meshAttrib[nameAttrib].push_back(new DataDescriptor<Container, n, Type>(nameProp, ptr));
				meshAttribProp[nameAttrib].push_back(PlyProperty(type, nameProp));
			}

			void AddMeshAttrib(const std::string& name, int cnt)
			{
				meshAttribCnt[name] = cnt;
			}

			void GetMeshAttrib(std::string filename)
			{
				nanoply::Info info(filename);
				if (info.errInfo == nanoply::NNP_OK)
				{
					for (size_t i = 0; i < info.elemVec.size(); i++)
					{
						if (info.elemVec[i].plyElem == NNP_UNKNOWN_ELEM && info.elemVec[i].name != "camera")
							meshAttribCnt[info.elemVec[i].name] = info.elemVec[i].cnt;
					}
				}
			}
      
      bool CreateVertexAttribDescriptor(const PointerToAttribute* ptr)
      {
        return CreateAttribDescriptor(ptr, vertexAttrib, vertexAttribProp);
      }

      bool CreateEdgeAttribDescriptor(const PointerToAttribute* ptr)
      {
        return CreateAttribDescriptor(ptr, edgeAttrib, edgeAttribProp);
      }

      bool CreateFaceAttribDescriptor(const PointerToAttribute* ptr)
      {
        return CreateAttribDescriptor(ptr, faceAttrib, faceAttribProp);
      }

      bool AddVertexAttrib(MeshType& m, std::set<PointerToAttribute>& ptrAttrib, PlyElement* elem)
      {
        return AddCustomAttrib<0>(m, ptrAttrib, elem, vertexAttrib);
      }

      bool AddEdgeAttrib(MeshType& m, std::set<PointerToAttribute>& ptrAttrib, PlyElement* elem)
      {
        return AddCustomAttrib<1>(m, ptrAttrib, elem, edgeAttrib);
      }

      bool AddFaceAttrib(MeshType& m, std::set<PointerToAttribute>& ptrAttrib, PlyElement* elem)
      {
        return AddCustomAttrib<2>(m, ptrAttrib, elem, faceAttrib);
      }

     

     private:

			const std::vector<std::string> pointSuffix = { "x", "y", "z", "w" };
			const std::vector<std::string> colorSuffix = { "r", "g", "b", "a" };
      const std::map<std::string, int> pointSuffixMap = { { "x", 0 },{ "y", 1 },{ "z", 2 },{ "w", 3 } };
      const std::map<std::string, int> colorSuffixMap = { {"r", 0},{ "g", 1 },{ "b", 2 },{ "a", 3 } };
      const std::string separator = std::string("@_.#$>");

			template<class Container, class Type, int n>
			void AddAttribDescriptor(const std::string& name, PlyType type, void* ptr, ElementDescriptor::PropertyDescriptor& attrib, std::vector<PlyProperty>& attribProp)
			{
				attrib.push_back(new DataDescriptor<Container, n, Type>(name, ptr));
				attribProp.push_back(PlyProperty(type, name));
			}
			
      template<class Container>
      void AddScalarAttribDescriptor(const PointerToAttribute* ptr, ElementDescriptor::PropertyDescriptor& attrib, std::vector<PlyProperty>& attribProp)
      {
        AddAttribDescriptor<Container, Container, 1>(ptr->_name, getEntity<Container>(), ptr->_handle->DataBegin(), attrib, attribProp);
      }

			template<class Container>
			void AddPointAttribDescriptor(const PointerToAttribute* ptr, ElementDescriptor::PropertyDescriptor& attrib, std::vector<PlyProperty>& attribProp)
			{
				int size = int(Container::Dimension);
				std::string name(ptr->_name);
				Container* tmpPtr = (Container*)ptr->_handle->DataBegin();
				for (int i = 0 ; i < size; i++)
					AddAttribDescriptor<Container, typename Container::ScalarType, 1>((name + "." + pointSuffix[i]), getEntity<typename Container::ScalarType>(), &(*tmpPtr)[i], attrib, attribProp);
			}

			template<class Container>
			void AddColorAttribDescriptor(const PointerToAttribute* ptr, ElementDescriptor::PropertyDescriptor& attrib, std::vector<PlyProperty>& attribProp)
			{
				int size = int(Container::Dimension);
				std::string name(ptr->_name);
				Container* tmpPtr = (Container*)ptr->_handle->DataBegin();
				for (int i = 0; i < size; i++)
					AddAttribDescriptor<Container, typename Container::ScalarType, 1>((name + "." + colorSuffix[i]), getEntity<typename Container::ScalarType>(), &(*tmpPtr)[i], attrib, attribProp);
			}

      template<class Container, class Type>
      void AddListAttribDescriptor(const PointerToAttribute* ptr, ElementDescriptor::PropertyDescriptor& attrib, std::vector<PlyProperty>& attribProp)
      {
        AddAttribDescriptor<Container, Type, 0>(ptr->_name, getEntityList<Type>(), ptr->_handle->DataBegin(), attrib, attribProp);
      }

      bool CreateAttribDescriptor(const PointerToAttribute* ptr, ElementDescriptor::PropertyDescriptor& attrib, std::vector<PlyProperty>& attribProp)
			{
        if (ptr->_type == std::type_index(typeid(unsigned char)))
          AddScalarAttribDescriptor<unsigned char >(ptr, attrib, attribProp);
        else if (ptr->_type == std::type_index(typeid(char)))
          AddScalarAttribDescriptor<char>(ptr, attrib, attribProp);
        else if (ptr->_type == std::type_index(typeid(unsigned short)))
          AddScalarAttribDescriptor<unsigned short>(ptr, attrib, attribProp);
        else if (ptr->_type == std::type_index(typeid(short)))
          AddScalarAttribDescriptor<short>(ptr, attrib, attribProp);
        else if (ptr->_type == std::type_index(typeid(unsigned int)))
          AddScalarAttribDescriptor<unsigned int>(ptr, attrib, attribProp);
        else if (ptr->_type == std::type_index(typeid(int)))
          AddScalarAttribDescriptor<int>(ptr, attrib, attribProp);
        else if (ptr->_type == std::type_index(typeid(float)))
          AddScalarAttribDescriptor<float>(ptr, attrib, attribProp);
        else if (ptr->_type == std::type_index(typeid(double)))
          AddScalarAttribDescriptor<double>(ptr, attrib, attribProp);
				else if (ptr->_type == std::type_index(typeid(vcg::Point2s)))
					AddPointAttribDescriptor<vcg::Point2s>(ptr, attrib, attribProp);
				else if (ptr->_type == std::type_index(typeid(vcg::Point2i)))
					AddPointAttribDescriptor<vcg::Point2i>(ptr, attrib, attribProp);
				else if (ptr->_type == std::type_index(typeid(vcg::Point2f)))
					AddPointAttribDescriptor<vcg::Point2f>(ptr, attrib, attribProp);
				else if (ptr->_type == std::type_index(typeid(vcg::Point2d)))
					AddPointAttribDescriptor<vcg::Point2d>(ptr, attrib, attribProp);
				else if (ptr->_type == std::type_index(typeid(vcg::Point3s)))
					AddPointAttribDescriptor<vcg::Point3s>(ptr, attrib, attribProp);
				else if (ptr->_type == std::type_index(typeid(vcg::Point3i)))
					AddPointAttribDescriptor<vcg::Point3i>(ptr, attrib, attribProp);
				else if (ptr->_type == std::type_index(typeid(vcg::Point3f)))
					AddPointAttribDescriptor<vcg::Point3f>(ptr, attrib, attribProp);
				else if (ptr->_type == std::type_index(typeid(vcg::Point3d)))
					AddPointAttribDescriptor<vcg::Point3d>(ptr, attrib, attribProp);
				else if (ptr->_type == std::type_index(typeid(vcg::Point4s)))
					AddPointAttribDescriptor<vcg::Point4s>(ptr, attrib, attribProp);
				else if (ptr->_type == std::type_index(typeid(vcg::Point4i)))
					AddPointAttribDescriptor<vcg::Point4i>(ptr, attrib, attribProp);
				else if (ptr->_type == std::type_index(typeid(vcg::Point4f)))
					AddPointAttribDescriptor<vcg::Point4f>(ptr, attrib, attribProp);
				else if (ptr->_type == std::type_index(typeid(vcg::Point4d)))
					AddPointAttribDescriptor<vcg::Point4d>(ptr, attrib, attribProp);
				else if (ptr->_type == std::type_index(typeid(vcg::Color4b)))
					AddColorAttribDescriptor<vcg::Color4b>(ptr, attrib, attribProp);
				else if (ptr->_type == std::type_index(typeid(vcg::Color4f)))
					AddColorAttribDescriptor<vcg::Color4f>(ptr, attrib, attribProp);
				else if (ptr->_type == std::type_index(typeid(vcg::Color4d)))
					AddColorAttribDescriptor<vcg::Color4d>(ptr, attrib, attribProp);
        else if (ptr->_type == std::type_index(typeid(vcg::Color4d)))
          AddColorAttribDescriptor<vcg::Color4d>(ptr, attrib, attribProp);
        else if (ptr->_type == std::type_index(typeid(std::vector<unsigned char>)))
          AddListAttribDescriptor<std::vector<unsigned char>, unsigned char>(ptr, attrib, attribProp);
        else if (ptr->_type == std::type_index(typeid(std::vector<char>)))
          AddListAttribDescriptor<std::vector<char>,char>(ptr, attrib, attribProp);
        else if (ptr->_type == std::type_index(typeid(std::vector<unsigned short>)))
          AddListAttribDescriptor<std::vector<unsigned short>, unsigned short>(ptr, attrib, attribProp);
        else if (ptr->_type == std::type_index(typeid(std::vector<short>)))
          AddListAttribDescriptor<std::vector<short>, short>(ptr, attrib, attribProp);
        else if (ptr->_type == std::type_index(typeid(std::vector<unsigned int>)))
          AddListAttribDescriptor<std::vector<unsigned int>, unsigned int>(ptr, attrib, attribProp);
        else if (ptr->_type == std::type_index(typeid(std::vector<int>)))
          AddListAttribDescriptor<std::vector<int>, int>(ptr, attrib, attribProp);
        else if (ptr->_type == std::type_index(typeid(std::vector<float>)))
          AddListAttribDescriptor<std::vector<float>, float>(ptr, attrib, attribProp);
        else if (ptr->_type == std::type_index(typeid(std::vector<double>)))
          AddListAttribDescriptor<std::vector<double>, double>(ptr, attrib, attribProp);
				else
					return false;
				return true;
			}

      
      template <size_t ActionType, class Container>
      void AddScalarAttrib(MeshType& m, const std::string& name, PlyType& type)
      {
        if (ActionType == 0) //vertex
        {
          auto h = vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<Container>(m, name);
          AddVertexAttribDescriptor<Container, Container, 1>(name, type, h._handle->DataBegin());
        }
        else if (ActionType == 1) //Edge
        {
          auto h = vcg::tri::Allocator<MeshType>::template GetPerEdgeAttribute<Container>(m, name);
          AddEdgeAttribDescriptor<Container, Container, 1>(name, type, h._handle->DataBegin());
        }
        else if (ActionType == 2) //Face
        {
          auto h = vcg::tri::Allocator<MeshType>::template GetPerFaceAttribute<Container>(m, name);
          AddFaceAttribDescriptor<Container, Container, 1>(name, type, h._handle->DataBegin());
        }
      }

      template <size_t ActionType, class Container>
      void AddPointAttrib(MeshType& m, const std::string& name, std::vector<PlyProperty>& prop)
      {
        if (ActionType == 0) //vertex
        {
          auto h = vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<Container>(m, name);
          for (size_t i = 0; i < prop.size(); i++)
            AddVertexAttribDescriptor<Container, typename Container::ScalarType, 1>(prop[i].name, prop[i].type, &h[0][i]);
        }
        else if (ActionType == 1) //Edge
        {
          auto h = vcg::tri::Allocator<MeshType>::template GetPerEdgeAttribute<Container>(m, name);
          for (size_t i = 0; i < prop.size(); i++)
            AddEdgeAttribDescriptor<Container, typename Container::ScalarType, 1>(prop[i].name, prop[i].type, &h[0][i]);
        }
        else if (ActionType == 2) //Face
        {
          auto h = vcg::tri::Allocator<MeshType>::template GetPerFaceAttribute<Container>(m, name);
          for (size_t i = 0; i < prop.size(); i++)
            AddFaceAttribDescriptor<Container, typename Container::ScalarType, 1>(prop[i].name, prop[i].type, &h[0][i]);
        }
      }


      template <size_t ActionType, class Container, class Type>
      void AddListAttrib(MeshType& m, const std::string& name, PlyType& type)
      {
        if (ActionType == 0) //vertex
        {
          auto h = vcg::tri::Allocator<MeshType>::template GetPerVertexAttribute<Container>(m, name);
          AddVertexAttribDescriptor<Container, Type, 0>(name, type, h._handle->DataBegin());
        }
        else if (ActionType == 1) //Edge
        {
          auto h = vcg::tri::Allocator<MeshType>::template GetPerEdgeAttribute<Container>(m, name);
          AddEdgeAttribDescriptor<Container, Type, 0>(name, type, h._handle->DataBegin());
        }
        else if (ActionType == 2) //Face
        {
          auto h = vcg::tri::Allocator<MeshType>::template GetPerFaceAttribute<Container>(m, name);
          AddFaceAttribDescriptor<Container, Type, 0>(name, type, h._handle->DataBegin());
        }
      }


      template <size_t ActionType>
      bool AddCustomAttrib(MeshType& m, std::set<PointerToAttribute>& ptrAttrib, PlyElement* elem, ElementDescriptor::PropertyDescriptor& attrib)
      {
        //Custom attribute with already a data descriptor
        std::set<std::string> validCustomAttrib;
        for (size_t i = 0; i < attrib.size(); i++)
        {
          if (attrib[i]->base == NULL)
          {
            typename std::set<PointerToAttribute>::iterator ai;
            for (ai = ptrAttrib.begin(); ai != ptrAttrib.end(); ++ai)
            {
              if (attrib[i]->name == (*ai)._name)
              {
                attrib[i]->base = ai->_handle->DataBegin();
                validCustomAttrib.insert(attrib[i]->name);
                break;
              }
            }
          }
          else
            validCustomAttrib.insert(attrib[i]->name);
        }
        //Get custom properties without a data descriptor 
        std::map<std::string, std::pair<std::vector<PlyProperty>, int>> customAttribName;
        std::vector<std::string> attribName((*elem).propVec.size());
        for (size_t i = 0; i < (*elem).propVec.size(); i++)
        {
          if ((*elem).propVec[i].elem == NNP_UNKNOWN_ENTITY)
          {
            if (validCustomAttrib.find((*elem).propVec[i].name) == validCustomAttrib.end())
            {
              std::size_t found = (*elem).propVec[i].name.find_first_of(separator);
              if (found != std::string::npos)
              {
                attribName[i] = (*elem).propVec[i].name.substr(0, found);
                std::string suffix = (*elem).propVec[i].name.substr(found+1);
                std::map<std::string, int>::const_iterator it1, it2;
                it1 = pointSuffixMap.find(suffix);
                it2 = colorSuffixMap.find(suffix);
                if (it1 != pointSuffixMap.cend())
                  customAttribName[attribName[i]].second |= (1 << (*it1).second);
                else if (it2 != colorSuffixMap.cend())
                  customAttribName[attribName[i]].second |= (1 << ((*it2).second + 4));
              }
              else
                attribName[i] = (*elem).propVec[i].name;
              customAttribName[attribName[i]].first.push_back((*elem).propVec[i]);
              
            }
          }
        }
        //Create the attribute and the data descriptor
        std::map<std::string, std::pair<std::vector<PlyProperty>, int>>::iterator mapIter;
        for (mapIter = customAttribName.begin(); mapIter != customAttribName.end(); mapIter++)
        {
          unsigned int bitType = 0;
          std::vector<PlyProperty> tempProp((*mapIter).second.first.size());
          int checkMask = (1 << (*mapIter).second.first.size()) - 1;
          for (size_t i = 0; i < (*mapIter).second.first.size(); i++)
          {
            PlyProperty& p = (*mapIter).second.first[i];
            bitType |= p.type;
            std::size_t found = p.name.find_first_of(separator);
            if (found != std::string::npos)
            {
              std::string suffix = p.name.substr(found + 1);
              if ((*mapIter).second.second == checkMask)
                tempProp[(*pointSuffixMap.find(suffix)).second] = p;
              else if ((*mapIter).second.second == (checkMask << 4))
                tempProp[(*colorSuffixMap.find(suffix)).second] = p;
              else
                tempProp[i] = p;
            }
            else
              tempProp[i] = p;
          }
          unsigned int r = (bitType & (~tempProp[0].type)); (void)r;
          if (tempProp.size() > 1  && ((bitType & (~tempProp[0].type)) == 0))
          {
            if (tempProp.size() == 2)
            {
              switch (tempProp[0].type)
              {
              case NNP_FLOAT32: AddPointAttrib<ActionType, vcg::Point2f>(m, (*mapIter).first, tempProp); break;
              case NNP_FLOAT64: AddPointAttrib<ActionType, vcg::Point2d>(m, (*mapIter).first, tempProp); break;
              case NNP_INT8: AddPointAttrib<ActionType, vcg::Point2<char>>(m, (*mapIter).first, tempProp); break;
              case NNP_INT16: AddPointAttrib<ActionType, vcg::Point2s>(m, (*mapIter).first, tempProp); break;
              case NNP_INT32: AddPointAttrib<ActionType, vcg::Point2i>(m, (*mapIter).first, tempProp); break;
              case NNP_UINT8: AddPointAttrib<ActionType, vcg::Point2<unsigned char>>(m, (*mapIter).first, tempProp); break;
              case NNP_UINT16: AddPointAttrib<ActionType, vcg::Point2<unsigned short>>(m, (*mapIter).first, tempProp); break;
              case NNP_UINT32: AddPointAttrib<ActionType, vcg::Point2<unsigned int>>(m, (*mapIter).first, tempProp); break;
			  default: assert(0);
              }
            }
            else if (tempProp.size() == 3)
            {
              switch (tempProp[0].type)
              {
              case NNP_FLOAT32: AddPointAttrib<ActionType, vcg::Point3f>(m, (*mapIter).first, tempProp); break;
              case NNP_FLOAT64: AddPointAttrib<ActionType, vcg::Point3d>(m, (*mapIter).first, tempProp); break;
              case NNP_INT8: AddPointAttrib<ActionType, vcg::Point3<char>>(m, (*mapIter).first, tempProp); break;
              case NNP_INT16: AddPointAttrib<ActionType, vcg::Point3s>(m, (*mapIter).first, tempProp); break;
              case NNP_INT32: AddPointAttrib<ActionType, vcg::Point3i>(m, (*mapIter).first, tempProp); break;
              case NNP_UINT8: AddPointAttrib<ActionType, vcg::Point3<unsigned char>>(m, (*mapIter).first, tempProp); break;
              case NNP_UINT16: AddPointAttrib<ActionType, vcg::Point3<unsigned short>>(m, (*mapIter).first, tempProp); break;
              case NNP_UINT32: AddPointAttrib<ActionType, vcg::Point3<unsigned int>>(m, (*mapIter).first, tempProp); break;
			  default: assert(0);
              }
            }
            else if (tempProp.size() == 4 && (*mapIter).second.second != (checkMask << 4))
            {
              switch (tempProp[0].type)
              {
              case NNP_FLOAT32: AddPointAttrib<ActionType, vcg::Point4d>(m, (*mapIter).first, tempProp); break;
              case NNP_FLOAT64: AddPointAttrib<ActionType, vcg::Point4d>(m, (*mapIter).first, tempProp); break;
              case NNP_INT8: AddPointAttrib<ActionType, vcg::Point4<char>>(m, (*mapIter).first, tempProp); break;
              case NNP_INT16: AddPointAttrib<ActionType, vcg::Point4s>(m, (*mapIter).first, tempProp); break;
              case NNP_INT32: AddPointAttrib<ActionType, vcg::Point4i>(m, (*mapIter).first, tempProp); break;
              case NNP_UINT8: AddPointAttrib<ActionType, vcg::Point4<unsigned char>>(m, (*mapIter).first, tempProp); break;
              case NNP_UINT16: AddPointAttrib<ActionType, vcg::Point4<unsigned short>>(m, (*mapIter).first, tempProp); break;
              case NNP_UINT32: AddPointAttrib<ActionType, vcg::Point4<unsigned int>>(m, (*mapIter).first, tempProp); break;
			  default: assert(0);
              }
            }
            else if (tempProp.size() == 4)
            {
              switch (tempProp[0].type)
              {
              case NNP_INT8:
              case NNP_INT16:
              case NNP_INT32: 
              case NNP_UINT16:
              case NNP_UINT32:
              case NNP_FLOAT32: AddPointAttrib<ActionType, vcg::Color4f>(m, (*mapIter).first, tempProp); break;
              case NNP_FLOAT64: AddPointAttrib<ActionType, vcg::Color4d>(m, (*mapIter).first, tempProp); break;
              case NNP_UINT8: AddPointAttrib<ActionType, vcg::Color4b>(m, (*mapIter).first, tempProp); break;
			  default: assert(0);
              }
            }
          }
          else
          {
            for (size_t i = 0; i < tempProp.size(); i++)
            {
              switch (tempProp[i].type)
              {
              case NNP_FLOAT32: AddScalarAttrib<ActionType, float>(m, tempProp[i].name, tempProp[i].type); break;
              case NNP_FLOAT64: AddScalarAttrib<ActionType, double>(m, tempProp[i].name, tempProp[i].type); break;
              case NNP_INT8: AddScalarAttrib<ActionType, char>(m, tempProp[i].name, tempProp[i].type); break;
              case NNP_INT16: AddScalarAttrib<ActionType, short>(m, tempProp[i].name, tempProp[i].type); break;
              case NNP_INT32: AddScalarAttrib<ActionType, int>(m, tempProp[i].name, tempProp[i].type); break;
              case NNP_UINT8: AddScalarAttrib<ActionType, unsigned char>(m, tempProp[i].name, tempProp[i].type); break;
              case NNP_UINT16: AddScalarAttrib<ActionType, unsigned short>(m, tempProp[i].name, tempProp[i].type); break;
              case NNP_UINT32: AddScalarAttrib<ActionType, unsigned int>(m, tempProp[i].name, tempProp[i].type); break;
              case NNP_LIST_UINT8_UINT32:
              case NNP_LIST_INT8_UINT32: AddListAttrib<ActionType, std::vector<unsigned int>, unsigned int>(m, tempProp[i].name, tempProp[i].type); break;
              case NNP_LIST_UINT8_INT32:
              case NNP_LIST_INT8_INT32: AddListAttrib<ActionType, std::vector<int>, int>(m, tempProp[i].name, tempProp[i].type); break;
              case NNP_LIST_UINT8_FLOAT32:
              case NNP_LIST_INT8_FLOAT32: AddListAttrib<ActionType, std::vector<float>, float>(m, tempProp[i].name, tempProp[i].type); break;
              case NNP_LIST_UINT8_FLOAT64:
              case NNP_LIST_INT8_FLOAT64: AddListAttrib<ActionType, std::vector<double>, double>(m, tempProp[i].name, tempProp[i].type); break;
              case NNP_LIST_UINT8_UINT8:
              case NNP_LIST_INT8_UINT8: AddListAttrib<ActionType, std::vector<unsigned char>, unsigned char>(m, tempProp[i].name, tempProp[i].type); break;
              case NNP_LIST_UINT8_INT8:
              case NNP_LIST_INT8_INT8: AddListAttrib<ActionType, std::vector<char>, char>(m, tempProp[i].name, tempProp[i].type); break;
              case NNP_LIST_UINT8_UINT16:
              case NNP_LIST_INT8_UINT16: AddListAttrib<ActionType, std::vector<unsigned short>, unsigned short>(m, tempProp[i].name, tempProp[i].type); break;
              case NNP_LIST_UINT8_INT16:
              case NNP_LIST_INT8_INT16: AddListAttrib<ActionType, std::vector<short>, short>(m, tempProp[i].name, tempProp[i].type); break;
			  default: assert(0);
              }
            }
          }
        }
        return true;
      }
      
      

		};
		

    template <class T>
    struct OcfManager 
    {
      typedef typename T::VertexType VType;
      typedef typename T::FaceType FType;
      typedef typename T::EdgeType EType;
      

      template <bool f = std::is_same<typename T::VertContainer, vcg::vertex::vector_ocf<VType>>::value>
      static unsigned int EnableVertexOcf(typename T::VertContainer& cont, unsigned int mask) {
        (void)cont; (void)mask;
        return 0;
      }

      template <>
      static unsigned int EnableVertexOcf<true>(typename T::VertContainer& cont, unsigned int mask )
      { 
        unsigned int enabledMask = 0;
        if ((mask & BitMask::IO_VERTNORMAL) && VType::HasNormalOcf() && !cont.IsNormalEnabled())
        {
          cont.EnableNormal();
          enabledMask |= BitMask::IO_VERTNORMAL;
        }
        if ((mask & BitMask::IO_VERTCOLOR) && VType::HasColorOcf() && !cont.IsColorEnabled())
        {
          cont.EnableColor();
          enabledMask |= BitMask::IO_VERTCOLOR;
        }
        if ((mask & BitMask::IO_VERTQUALITY) && VType::HasQualityOcf() && !cont.IsQualityEnabled())
        {
          cont.EnableQuality();
          enabledMask |= BitMask::IO_VERTQUALITY;
        }
        if ((mask & BitMask::IO_VERTCURV) && VType::HasCurvatureOcf() && !cont.IsCurvatureEnabled())
        {
          cont.EnableCurvature();
          enabledMask |= BitMask::IO_VERTCURV;
        }
        if ((mask & BitMask::IO_VERTCURVDIR) && VType::HasCurvatureDirOcf() && !cont.IsCurvatureDirEnabled())
        {
          cont.EnableCurvatureDir();
          enabledMask |= BitMask::IO_VERTCURVDIR;
        }
        if ((mask & BitMask::IO_VERTRADIUS) && VType::HasRadiusOcf() && !cont.IsRadiusEnabled())
        {
          cont.EnableRadius();
          enabledMask |= BitMask::IO_VERTRADIUS;
        }
        if ((mask & BitMask::IO_VERTTEXCOORD) && VType::HasTexCoordOcf() && !cont.IsTexCoordEnabled())
        {
          cont.EnableTexCoord();
          enabledMask |= BitMask::IO_VERTTEXCOORD;
        }
        return enabledMask; 
      };


      template <bool f = std::is_same<typename T::FaceContainer, vcg::face::vector_ocf<FType>>::value>
      static unsigned int EnableFaceOcf(typename T::FaceContainer& cont, unsigned int mask) {
        (void)cont; (void)mask;
        return 0;
      }

      template <>
      static unsigned int EnableFaceOcf<true>(typename T::FaceContainer& cont, unsigned int mask)
      {
        unsigned int enabledMask = 0;
        if ((mask & BitMask::IO_FACENORMAL) && FType::HasNormalOcf() && !cont.IsNormalEnabled())
        {
          cont.EnableNormal();
          enabledMask |= BitMask::IO_FACENORMAL;
        }
        if ((mask & BitMask::IO_FACECOLOR) && FType::HasColorOcf() && !cont.IsColorEnabled())
        {
          cont.EnableColor();
          enabledMask |= BitMask::IO_FACECOLOR;
        }
        if ((mask & BitMask::IO_FACEQUALITY) && FType::HasQualityOcf() && !cont.IsQualityEnabled())
        {
          cont.EnableQuality();
          enabledMask |= BitMask::IO_FACEQUALITY;
        }
        if ((mask & BitMask::IO_FACECURVDIR) && FType::HasCurvatureDirOcf() && !cont.IsCurvatureDirEnabled())
        {
          cont.EnableCurvatureDir();
          enabledMask |= BitMask::IO_FACECURVDIR;
        }
        if ((mask & BitMask::IO_WEDGCOLOR) && FType::HasWedgeColorOcf() && !cont.IsWedgeColorEnabled())
        {
          cont.EnableWedgeColor();
          enabledMask |= BitMask::IO_WEDGCOLOR;
        }
        if ((mask & BitMask::IO_WEDGNORMAL) && FType::HasWedgeNormalOcf() && !cont.IsWedgeNormalEnabled())
        {
          cont.EnableWedgeNormal();
          enabledMask |= BitMask::IO_WEDGNORMAL;
        }
        if ((mask & BitMask::IO_WEDGTEXCOORD) && FType::HasWedgeTexCoordOcf() && !cont.IsWedgeTexCoordEnabled())
        {
          cont.EnableWedgeTexCoord();
          enabledMask |= BitMask::IO_WEDGTEXCOORD;
        }
        if ((mask & BitMask::IO_WEDGTEXMULTI) && FType::HasWedgeTexCoordOcf())
        {
          if (!cont.IsWedgeTexCoordEnabled())
            cont.EnableWedgeTexCoord();
          enabledMask |= BitMask::IO_WEDGTEXMULTI;
        }
        return enabledMask;
      };



      template <bool f = std::is_same<typename T::VertContainer, vcg::vertex::vector_ocf<VType>>::value>
      static unsigned int VertexOcfMask(typename T::VertContainer& cont) {
        (void)cont;
        return 0;
      }

      template <>
      static unsigned int VertexOcfMask<true>(typename T::VertContainer& cont)
      {
        unsigned int enabledMask = 0;
        if (VType::HasNormalOcf() && cont.IsNormalEnabled())
          enabledMask |= BitMask::IO_VERTNORMAL;
        if (VType::HasColorOcf() && cont.IsColorEnabled())
          enabledMask |= BitMask::IO_VERTCOLOR;
        if (VType::HasQualityOcf() && cont.IsQualityEnabled())
          enabledMask |= BitMask::IO_VERTQUALITY;
        if (VType::HasCurvatureOcf() && cont.IsCurvatureEnabled())
          enabledMask |= BitMask::IO_VERTCURV;
        if (VType::HasCurvatureDirOcf() && cont.IsCurvatureDirEnabled())
          enabledMask |= BitMask::IO_VERTCURVDIR;
        if (VType::HasRadiusOcf() && cont.IsRadiusEnabled())
          enabledMask |= BitMask::IO_VERTRADIUS;
        if (VType::HasTexCoordOcf() && cont.IsTexCoordEnabled())
          enabledMask |= BitMask::IO_VERTTEXCOORD;
        return enabledMask;
      };


      template <bool f = std::is_same<typename T::FaceContainer, vcg::face::vector_ocf<FType>>::value>
      static unsigned int FaceOcfMask(typename T::FaceContainer& cont) {
        (void)cont;
        return 0;
      }

      template <>
      static unsigned int FaceOcfMask<true>(typename T::FaceContainer& cont)
      {
        unsigned int enabledMask = 0;
        if (FType::HasNormalOcf() && cont.IsNormalEnabled())
          enabledMask |= BitMask::IO_FACENORMAL;
        if (FType::HasColorOcf() && cont.IsColorEnabled())
          enabledMask |= BitMask::IO_FACECOLOR;
        if (FType::HasQualityOcf() && cont.IsQualityEnabled())
          enabledMask |= BitMask::IO_FACEQUALITY;
        if (FType::HasCurvatureDirOcf() && cont.IsCurvatureDirEnabled())
          enabledMask |= BitMask::IO_FACECURVDIR;
        if (FType::HasWedgeColorOcf() && cont.IsWedgeColorEnabled())
          enabledMask |= BitMask::IO_WEDGCOLOR;
        if (FType::HasWedgeNormalOcf() && cont.IsWedgeNormalEnabled())
          enabledMask |= BitMask::IO_WEDGNORMAL;
        if (FType::HasWedgeTexCoordOcf() && cont.IsWedgeTexCoordEnabled())
        {
          enabledMask |= BitMask::IO_WEDGTEXCOORD;
          enabledMask |= BitMask::IO_WEDGTEXMULTI;
        }
        return enabledMask;
      };
            
    };

    
   
    static unsigned int GetFileBitMask(nanoply::Info& info)
    {
      unsigned int mask = 0;
      for (size_t j = 0; j < info.elemVec.size(); j++)
      {
        PlyElement& elem = info.elemVec[j];
        if (elem.plyElem == PlyElemEntity::NNP_VERTEX_ELEM)
        {
          for (size_t i = 0; i < elem.propVec.size(); i++)
          {
            switch (elem.propVec[i].elem)
            {
            case NNP_PXYZ: mask |= BitMask::IO_VERTCOORD; break;
            case NNP_NXYZ: mask |= BitMask::IO_VERTNORMAL; break;
            case NNP_CRGB:
            case NNP_CRGBA: mask |= BitMask::IO_VERTCOLOR; break;
            case NNP_BITFLAG: mask |= BitMask::IO_VERTFLAGS; break;
            case NNP_QUALITY: mask |= BitMask::IO_VERTQUALITY; break;
            case NNP_DENSITY:  mask |= BitMask::IO_VERTRADIUS; break;
            case NNP_TEXTURE2D:
            case NNP_TEXTURE3D: mask |= BitMask::IO_VERTTEXCOORD; break;
            case NNP_KH:
            case NNP_KG: mask |= BitMask::IO_VERTCURV; break;
            case NNP_K1:
            case NNP_K2:
            case NNP_K1DIR:
            case NNP_K2DIR: mask |= BitMask::IO_VERTCURVDIR; break;
            case NNP_UNKNOWN_ENTITY: mask |= BitMask::IO_VERTATTRIB; break;
			default: assert(0);
            }
          }
        }
        else if (elem.plyElem == PlyElemEntity::NNP_EDGE_ELEM)
        {
          for (size_t i = 0; i < elem.propVec.size(); i++)
          {
            switch (elem.propVec[i].elem)
            {
            case NNP_EDGE_V1:
            case NNP_EDGE_V2: mask |= BitMask::IO_EDGEINDEX; break;
            case NNP_CRGB:
            case NNP_CRGBA: mask |= BitMask::IO_EDGECOLOR; break;
            case NNP_BITFLAG: mask |= BitMask::IO_EDGEFLAGS; break;
            case NNP_QUALITY: mask |= BitMask::IO_EDGEQUALITY; break;
            case NNP_UNKNOWN_ENTITY: mask |= BitMask::IO_EDGEATTRIB; break;
			default: assert(0);
            }
          }
        }
        else if (elem.plyElem == PlyElemEntity::NNP_FACE_ELEM)
        {
          for (size_t i = 0; i < elem.propVec.size(); i++)
          {
            switch (elem.propVec[i].elem)
            {
            case NNP_FACE_VERTEX_LIST: mask |= BitMask::IO_FACEINDEX; break;
            case NNP_NXYZ: mask |= BitMask::IO_FACENORMAL; break;
            case NNP_CRGB:
            case NNP_CRGBA: mask |= BitMask::IO_FACECOLOR; break;
            case NNP_BITFLAG: mask |= BitMask::IO_FACEFLAGS; break;
            case NNP_QUALITY: mask |= BitMask::IO_FACEQUALITY; break;
            case NNP_K1:
            case NNP_K2:
            case NNP_K1DIR:
            case NNP_K2DIR: mask |= BitMask::IO_FACECURVDIR; break;
            case NNP_FACE_WEDGE_COLOR: mask |= BitMask::IO_WEDGCOLOR; break;
            case NNP_FACE_WEDGE_NORMAL: mask |= BitMask::IO_WEDGNORMAL; break;
            case NNP_FACE_WEDGE_TEX: mask |= BitMask::IO_WEDGTEXCOORD; break;
            case NNP_TEXTUREINDEX: mask |= BitMask::IO_WEDGTEXMULTI; break;
            case NNP_UNKNOWN_ENTITY: mask |= BitMask::IO_FACEATTRIB; break;
			default: assert(0);
            }
          }
        }
        else if (elem.plyElem == PlyElemEntity::NNP_UNKNOWN_ELEM)
        {
          if (elem.name == "camera")
            mask |= BitMask::IO_CAMERA;
          else
            mask |= BitMask::IO_MESHATTRIB;
        }
      }
      return mask;
    };



		static int LoadModel(const char* filename, MeshType& mesh, unsigned int bitMask, CustomAttributeDescriptor& custom)
		{
			nanoply::Info info(filename);
			if (info.errInfo != nanoply::NNP_OK)
				return info.errInfo;
      
      unsigned int headerMask = GetFileBitMask(info);
      bitMask &= headerMask;
      unsigned int ocfVertexMask = OcfManager<MeshType>::EnableVertexOcf(mesh.vert, bitMask);
      unsigned int ocfFaceMask = OcfManager<MeshType>::EnableFaceOcf(mesh.face, bitMask);
      
      //Camera
			ElementDescriptor cameraDescr(std::string("camera"));
			vcg::Point3<ScalarType> tra;
			vcg::Matrix44<ScalarType> rot;
			size_t count = info.GetElementCount(std::string("camera"));
			if (count > 0 && (bitMask & BitMask::IO_CAMERA))
			{
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("view_px"), &tra[0]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("view_py"), &tra[1]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("view_pz"), &tra[2]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("x_axisx"), &rot[0][0]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("x_axisy"), &rot[0][1]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("x_axisz"), &rot[0][2]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("y_axisx"), &rot[1][0]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("y_axisy"), &rot[1][1]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("y_axisz"), &rot[1][2]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("z_axisx"), &rot[2][0]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("z_axisy"), &rot[2][1]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("z_axisz"), &rot[2][2]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("focal"), &mesh.shot.Intrinsics.FocalMm));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("scalex"), &mesh.shot.Intrinsics.PixelSizeMm[0]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("scaley"), &mesh.shot.Intrinsics.PixelSizeMm[1]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("centerx"), &mesh.shot.Intrinsics.CenterPx[0]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("centery"), &mesh.shot.Intrinsics.CenterPx[1]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<int, 1, int>(std::string("viewportx"), &mesh.shot.Intrinsics.ViewportPx[0]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<int, 1, int>(std::string("viewporty"), &mesh.shot.Intrinsics.ViewportPx[1]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("k1"), &mesh.shot.Intrinsics.k[0]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("k2"), &mesh.shot.Intrinsics.k[1]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("k3"), &mesh.shot.Intrinsics.k[2]));
				cameraDescr.dataDescriptor.push_back(new DataDescriptor<ScalarType, 1, ScalarType>(std::string("k4"), &mesh.shot.Intrinsics.k[3]));
			}
			
			//Vertex
			std::vector<std::string> nameList;
			VertexType::Name(nameList);
			ElementDescriptor vertexDescr(NNP_VERTEX_ELEM);
			count = info.GetVertexCount();
			if (nameList.size() > 0 && count > 0)
			{
				vcg::tri::Allocator<MeshType>::AddVertices(mesh, count);
				if ((bitMask & BitMask::IO_VERTCOORD) && VertexType::HasCoord())
					vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 3, VertexCoordScalar>(NNP_PXYZ, (*mesh.vert.begin()).P().V()));
        if ((bitMask & BitMask::IO_VERTNORMAL) && vcg::tri::HasPerVertexNormal(mesh))
        {
          if (ocfVertexMask & BitMask::IO_VERTNORMAL)
            vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexNormalType, 3, VertexNormScalar>(NNP_NXYZ, (*mesh.vert.begin()).N().V()));
          else
            vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 3, VertexNormScalar>(NNP_NXYZ, (*mesh.vert.begin()).N().V()));
        }
        if ((bitMask & BitMask::IO_VERTCOLOR) && vcg::tri::HasPerVertexColor(mesh))
        {
          if (ocfVertexMask & BitMask::IO_VERTCOLOR)
            vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexColorType, 4, VertexColorScalar>(NNP_CRGBA, (*mesh.vert.begin()).C().V()));
          else
            vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 4, VertexColorScalar>(NNP_CRGBA, (*mesh.vert.begin()).C().V()));
        }
        if ((bitMask & BitMask::IO_VERTQUALITY) && vcg::tri::HasPerVertexQuality(mesh))
        {
          if (ocfVertexMask & BitMask::IO_VERTQUALITY)
            vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexQuality, 1, VertexQuality>(NNP_QUALITY, &(*mesh.vert.begin()).Q()));
          else
            vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 1, VertexQuality>(NNP_QUALITY, &(*mesh.vert.begin()).Q()));
        }
				if ((bitMask & BitMask::IO_VERTFLAGS) && vcg::tri::HasPerVertexFlags(mesh))
					vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 1, VertexFlag>(NNP_BITFLAG, &(*mesh.vert.begin()).Flags()));
        if ((bitMask & BitMask::IO_VERTRADIUS) && vcg::tri::HasPerVertexRadius(mesh))
        {
          if (ocfVertexMask & BitMask::IO_VERTRADIUS)
            vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexRadius, 1, VertexRadius>(NNP_DENSITY, &(*mesh.vert.begin()).R()));
          else
            vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 1, VertexRadius>(NNP_DENSITY, &(*mesh.vert.begin()).R()));
        }
        if ((bitMask & BitMask::IO_VERTTEXCOORD) && vcg::tri::HasPerVertexTexCoord(mesh))
        {
          if (ocfVertexMask & BitMask::IO_VERTTEXCOORD)
            vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexTexCoordType, 2, VertexTexScalar>(NNP_TEXTURE2D, (*mesh.vert.begin()).T().P().V()));
          else
            vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 2, VertexTexScalar>(NNP_TEXTURE2D, (*mesh.vert.begin()).T().P().V()));
        }
        if ((bitMask & BitMask::IO_VERTCURV) && vcg::tri::HasPerVertexCurvature(mesh))
				{
          if (ocfVertexMask & BitMask::IO_VERTCURV)
          {
            vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexCurType, 1, VertexCurScalar>(NNP_KG, &(*mesh.vert.begin()).Kg()));
            vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexCurType, 1, VertexCurScalar>(NNP_KH, &(*mesh.vert.begin()).Kh()));
          }
          else
          {
            vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 1, VertexCurScalar>(NNP_KG, &(*mesh.vert.begin()).Kg()));
            vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 1, VertexCurScalar>(NNP_KH, &(*mesh.vert.begin()).Kh()));
          }
				}
				if ((bitMask & BitMask::IO_VERTCURVDIR) && vcg::tri::HasPerVertexCurvatureDir(mesh))
				{
          if (ocfVertexMask & BitMask::IO_VERTCURVDIR)
          {
            vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexCurDirType, 1, VertexDirCurScalar>(NNP_K1, &(*mesh.vert.begin()).K1()));
            vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexCurDirType, 1, VertexDirCurScalar>(NNP_K2, &(*mesh.vert.begin()).K2()));
            vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexCurDirType, 3, VertexDirCurVecScalar>(NNP_K1DIR, (*mesh.vert.begin()).PD1().V()));
            vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexCurDirType, 3, VertexDirCurVecScalar>(NNP_K2DIR, (*mesh.vert.begin()).PD2().V()));
          }
          else
          {
            vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 1, VertexDirCurScalar>(NNP_K1, &(*mesh.vert.begin()).K1()));
            vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 1, VertexDirCurScalar>(NNP_K2, &(*mesh.vert.begin()).K2()));
            vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 3, VertexDirCurVecScalar>(NNP_K1DIR, (*mesh.vert.begin()).PD1().V()));
            vertexDescr.dataDescriptor.push_back(new DataDescriptor<VertexType, 3, VertexDirCurVecScalar>(NNP_K2DIR, (*mesh.vert.begin()).PD2().V()));
          }
        }
        if (bitMask & BitMask::IO_VERTATTRIB)
        {
          custom.AddVertexAttrib(mesh, mesh.vert_attr, info.GetVertexElement());
          for (size_t i = 0; i < custom.vertexAttrib.size(); i++)
          {
            if ((*custom.vertexAttrib[i]).base != NULL)
              vertexDescr.dataDescriptor.push_back(custom.vertexAttrib[i]);
          }
        }
			}

			//Edge
			nameList.clear();
			EdgeType::Name(nameList);
			ElementDescriptor edgeDescr(NNP_EDGE_ELEM);
			count = info.GetEdgeCount();
			std::vector<vcg::Point2i> edgeIndex;
			if (nameList.size() > 0 && count > 0)
			{
				vcg::tri::Allocator<MeshType>::AddEdges(mesh, count);
				if ((bitMask & BitMask::IO_EDGEINDEX) && MeshType::EdgeType::HasVertexRef())
				{
					edgeIndex.resize(count);
					edgeDescr.dataDescriptor.push_back(new DataDescriptor<vcg::Point2i, 1, int>(NNP_EDGE_V1, &(*edgeIndex.begin()).V()[0]));
					edgeDescr.dataDescriptor.push_back(new DataDescriptor<vcg::Point2i, 1, int>(NNP_EDGE_V2, &(*edgeIndex.begin()).V()[1]));
				}
				if ((bitMask & BitMask::IO_EDGEQUALITY) && vcg::tri::HasPerEdgeQuality(mesh))
					edgeDescr.dataDescriptor.push_back(new DataDescriptor<EdgeType, 1, EdgeQuality>(NNP_QUALITY, &(*mesh.edge.begin()).Q()));
				if ((bitMask & BitMask::IO_EDGECOLOR) && vcg::tri::HasPerEdgeColor(mesh))
					edgeDescr.dataDescriptor.push_back(new DataDescriptor<EdgeType, 4, EdgeColorScalar>(NNP_CRGBA, (*mesh.edge.begin()).C().V()));
				if ((bitMask & BitMask::IO_EDGEFLAGS) && vcg::tri::HasPerEdgeFlags(mesh))
					edgeDescr.dataDescriptor.push_back(new DataDescriptor<EdgeType, 1, EdgeFlag>(NNP_BITFLAG, &(*mesh.edge.begin()).Flags()));
        if (bitMask & BitMask::IO_EDGEATTRIB)
        {
          custom.AddEdgeAttrib(mesh, mesh.edge_attr, info.GetEdgeElement());
          for (size_t i = 0; i < custom.edgeAttrib.size(); i++)
          {
            if ((*custom.edgeAttrib[i]).base != NULL)
              edgeDescr.dataDescriptor.push_back(custom.edgeAttrib[i]);
          }
        }
			}

			//Face
			nameList.clear();
			FaceType::Name(nameList);
			ElementDescriptor faceDescr(NNP_FACE_ELEM);
			count = info.GetFaceCount();
      std::vector<std::vector<unsigned int>> faceIndex;
			std::vector<vcg::ndim::Point<6, FaceTexScalar>> wedgeTexCoord;
			if (nameList.size() > 0 && count > 0)
			{
				vcg::tri::Allocator<MeshType>::AddFaces(mesh, count);
				if ((bitMask & BitMask::IO_FACEINDEX) && FaceType::HasVertexRef())
				{
					faceIndex.resize(count);
          faceDescr.dataDescriptor.push_back(new DataDescriptor<std::vector<unsigned int>,0, unsigned int>(NNP_FACE_VERTEX_LIST, &faceIndex[0]));
				}
				if ((bitMask & BitMask::IO_FACEFLAGS) && vcg::tri::HasPerFaceFlags(mesh))
					faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceType, 1, FaceFlag>(NNP_BITFLAG, &(*mesh.face.begin()).Flags()));
        if ((bitMask & BitMask::IO_FACECOLOR) && vcg::tri::HasPerFaceColor(mesh))
        {
          if (ocfFaceMask & BitMask::IO_FACECOLOR)
            faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceColorType, 4, FaceColorScalar>(NNP_CRGBA, (*mesh.face.begin()).C().V()));
          else
            faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceType, 4, FaceColorScalar>(NNP_CRGBA, (*mesh.face.begin()).C().V()));
        }
        if ((bitMask & BitMask::IO_FACEQUALITY) && vcg::tri::HasPerFaceQuality(mesh))
        {
          if (ocfFaceMask & BitMask::IO_FACEQUALITY)
            faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceQuality, 1, FaceQuality>(NNP_QUALITY, &(*mesh.face.begin()).Q()));
          else
            faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceType, 1, FaceQuality>(NNP_QUALITY, &(*mesh.face.begin()).Q()));
        }
        if ((bitMask & BitMask::IO_FACENORMAL) && vcg::tri::HasPerFaceNormal(mesh))
        {
          if (ocfFaceMask & BitMask::IO_FACENORMAL)
            faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceNormalType, 3, FaceNormScalar>(NNP_NXYZ, (*mesh.face.begin()).N().V()));
          else
            faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceType, 3, FaceNormScalar>(NNP_NXYZ, (*mesh.face.begin()).N().V()));
        }
        if ((bitMask & BitMask::IO_FACECURVDIR) && vcg::tri::HasPerFaceCurvatureDir(mesh))
				{
          if (ocfFaceMask & BitMask::IO_FACECURVDIR)
          {
            faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceCurDirType, 1, FaceDirCurScalar>(NNP_K1, &(*mesh.face.begin()).K1()));
            faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceCurDirType, 1, FaceDirCurScalar>(NNP_K2, &(*mesh.face.begin()).K2()));
            faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceCurDirType, 3, FaceDirCurVecScalar>(NNP_K1DIR, (*mesh.face.begin()).PD1().V()));
            faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceCurDirType, 3, FaceDirCurVecScalar>(NNP_K2DIR, (*mesh.face.begin()).PD2().V()));
          }
          else
          {
            faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceType, 1, FaceDirCurScalar>(NNP_K1, &(*mesh.face.begin()).K1()));
            faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceType, 1, FaceDirCurScalar>(NNP_K2, &(*mesh.face.begin()).K2()));
            faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceType, 3, FaceDirCurVecScalar>(NNP_K1DIR, (*mesh.face.begin()).PD1().V()));
            faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceType, 3, FaceDirCurVecScalar>(NNP_K2DIR, (*mesh.face.begin()).PD2().V()));
          }
				}
				if (((bitMask & BitMask::IO_WEDGTEXCOORD) || (bitMask & BitMask::IO_WEDGTEXMULTI)) && vcg::tri::HasPerWedgeTexCoord(mesh))
				{
					wedgeTexCoord.resize(count);
          faceDescr.dataDescriptor.push_back(new DataDescriptor<vcg::ndim::Point<6, FaceTexScalar>, 6, FaceTexScalar>(NNP_FACE_WEDGE_TEX, (*wedgeTexCoord.begin()).V()));
				}
        if ((bitMask & BitMask::IO_WEDGTEXMULTI) && vcg::tri::HasPerWedgeTexCoord(mesh))
        {
          if (ocfFaceMask & BitMask::IO_WEDGTEXMULTI)
            faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceTexCoordType, 1, short>(NNP_TEXTUREINDEX, &(*mesh.face.begin()).WT(0).N()));
          else
            faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceType, 1, short>(NNP_TEXTUREINDEX, &(*mesh.face.begin()).WT(0).N()));
        }
        if ((bitMask & BitMask::IO_WEDGCOLOR) && vcg::tri::HasPerWedgeColor(mesh))
        {
          if (ocfFaceMask & BitMask::IO_WEDGCOLOR)
            faceDescr.dataDescriptor.push_back(new DataDescriptor<WedgeColorType, 12, WedgeColorScalar>(NNP_FACE_WEDGE_COLOR, (*mesh.face.begin()).WC(0).V()));
          else
            faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceType, 12, WedgeColorScalar>(NNP_FACE_WEDGE_COLOR, (*mesh.face.begin()).WC(0).V()));
        }
        if ((bitMask & BitMask::IO_WEDGNORMAL) && vcg::tri::HasPerWedgeNormal(mesh))
        {
          if (ocfFaceMask & BitMask::IO_WEDGNORMAL)
            faceDescr.dataDescriptor.push_back(new DataDescriptor<WedgeNormalType, 9, WedgeNormalScalar>(NNP_FACE_WEDGE_NORMAL, (*mesh.face.begin()).WN(0).V()));
          else
            faceDescr.dataDescriptor.push_back(new DataDescriptor<FaceType, 9, WedgeNormalScalar>(NNP_FACE_WEDGE_NORMAL, (*mesh.face.begin()).WN(0).V()));
        }
        if (bitMask & BitMask::IO_FACEATTRIB)
        {
          custom.AddFaceAttrib(mesh, mesh.face_attr, info.GetFaceElement());
          for (size_t i = 0; i < custom.faceAttrib.size(); i++)
          {
            if ((*custom.faceAttrib[i]).base != NULL)
              faceDescr.dataDescriptor.push_back(custom.faceAttrib[i]);
          }
        }
			}

			std::vector<ElementDescriptor*> meshDescr;
			meshDescr.push_back(&cameraDescr);
			meshDescr.push_back(&vertexDescr);
			meshDescr.push_back(&edgeDescr);
			meshDescr.push_back(&faceDescr);

			//Mesh attribute
			if ((bitMask & BitMask::IO_MESHATTRIB))
			{
				typename CustomAttributeDescriptor::MapMeshAttribIter iter = custom.meshAttrib.begin();
				for (; iter != custom.meshAttrib.end(); iter++)
				{
					std::string name((*iter).first);
					meshDescr.push_back(new ElementDescriptor(name));
					count = info.GetElementCount(name);
					if (count > 1)
					{
						meshDescr.back()->dataDescriptor = (*iter).second;
					}
				}

			}
			if (!OpenModel(info, meshDescr))
				return info.errInfo;

			mesh.shot.SetViewPoint(tra);
			mesh.shot.Extrinsics.SetRot(rot);
      bool triangleMesh = true;
      for (size_t i = 0; i < faceIndex.size(); i++)
        if (faceIndex[i].size() > 3)
          triangleMesh = false;

      if (!triangleMesh && !vcg::tri::HasPolyInfo(mesh))
      {
        bool hasFaceFlags = (bitMask & BitMask::IO_FACEFLAGS) && vcg::tri::HasPerFaceFlags(mesh);
        bool hasFaceColor = (bitMask & BitMask::IO_FACECOLOR) && vcg::tri::HasPerFaceColor(mesh);
        bool hasFaceQuality = (bitMask & BitMask::IO_FACEQUALITY) && vcg::tri::HasPerFaceQuality(mesh);
        bool hasFaceNormal = (bitMask & BitMask::IO_FACENORMAL) && vcg::tri::HasPerFaceNormal(mesh);
        for (size_t i = 0; i < faceIndex.size(); i++)
        {
          if (faceIndex[i].size() >= 3)
          {
            for (int j = 0; j < 3; j++)
              mesh.face[i].V(j) = &mesh.vert[faceIndex[i][j]];
            if (hasFaceFlags)
              mesh.face[i].SetF(2);
            FaceIterator fi = vcg::tri::Allocator<MeshType>::AddFaces(mesh, faceIndex[i].size() - 3);
            for (size_t j = 3; j < faceIndex[i].size(); j++)
            {
              (*fi).V(0) = &mesh.vert[faceIndex[i][0]];
              (*fi).V(1) = &mesh.vert[faceIndex[i][j - 1]];
              (*fi).V(2) = &mesh.vert[faceIndex[i][j]];
              if (hasFaceFlags)
              {
                (*fi).SetFlags(mesh.face[i].Flags());
                (*fi).SetF(0);
                if (j == faceIndex[i].size() - 1)
                  (*fi).ClearF(2);
              }
              if (hasFaceColor)
                (*fi).C() = mesh.face[i].C();
              if (hasFaceQuality)
                (*fi).Q() = mesh.face[i].Q();
              if (hasFaceNormal)
                (*fi).N() = mesh.face[i].N();
              fi++;
            }
          }
          else
            mesh.face[i].V(0) = mesh.face[i].V(1) = mesh.face[i].V(2) = &mesh.vert[0];
        }
      }
      else
      {
        for (size_t i = 0; i < faceIndex.size(); i++)
        {
          mesh.face[i].Alloc(faceIndex[i].size());
          for (size_t j = 0; j < faceIndex[i].size(); j++)
            mesh.face[i].V(j) = &mesh.vert[faceIndex[i][j]];
        }
        for (size_t i = 0; i < wedgeTexCoord.size(); i++)
        {
          for (int j = 0; j < 3; j++)
          {
            mesh.face[i].WT(j).U() = wedgeTexCoord[i][j * 2];
            mesh.face[i].WT(j).V() = wedgeTexCoord[i][j * 2 + 1];
            mesh.face[i].WT(j).N() = mesh.face[i].WT(0).N();
          }
        }
      }
      
      for (size_t i = 0; i < edgeIndex.size(); i++)
			{
				mesh.edge[i].V(0) = &mesh.vert[edgeIndex[i].X()];
				mesh.edge[i].V(1) = &mesh.vert[edgeIndex[i].Y()];
			}

			for (size_t i = 0; i < cameraDescr.dataDescriptor.size(); i++)
				delete cameraDescr.dataDescriptor[i];
			for (size_t i = 0; i < vertexDescr.dataDescriptor.size(); i++)
				if (vertexDescr.dataDescriptor[i]->elem != NNP_UNKNOWN_ENTITY)
					delete vertexDescr.dataDescriptor[i];
			for (size_t i = 0; i < edgeDescr.dataDescriptor.size(); i++)
				if (edgeDescr.dataDescriptor[i]->elem != NNP_UNKNOWN_ENTITY)
					delete edgeDescr.dataDescriptor[i];
			for (size_t i = 0; i < faceDescr.dataDescriptor.size(); i++)
				if (faceDescr.dataDescriptor[i]->elem != NNP_UNKNOWN_ENTITY)
					delete faceDescr.dataDescriptor[i];
			mesh.textures = info.textureFile;
			return info.errInfo;
		}


		static int LoadModel(const char* filename, MeshType& mesh, unsigned int bitMask)
		{
			CustomAttributeDescriptor custom;
			return LoadModel(filename, mesh, bitMask, custom);
		}
	

		static bool SaveModel(const char* filename, MeshType& mesh, unsigned int bitMask, CustomAttributeDescriptor& custom, bool binary)
		{
      unsigned int ocfVertexMask = OcfManager<MeshType>::VertexOcfMask(mesh.vert);
      unsigned int ocfFaceMask = OcfManager<MeshType>::FaceOcfMask(mesh.face);
      
			//Camera
			std::vector<PlyProperty> cameraProp;
			ElementDescriptor cameraDescr(std::string("camera"));
			vcg::Point3<ScalarType> tra = mesh.shot.Extrinsics.Tra();
			vcg::Matrix44<ScalarType> rot = mesh.shot.Extrinsics.Rot();
			if (bitMask & BitMask::IO_CAMERA)
			{
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("view_px"), &tra[0]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("view_py"), &tra[1]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("view_pz"), &tra[2]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("x_axisx"), &rot[0][0]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("x_axisy"), &rot[0][1]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("x_axisz"), &rot[0][2]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("y_axisx"), &rot[1][0]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("y_axisy"), &rot[1][1]);
				PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("y_axisz"), &rot[1][2]);
        PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("z_axisx"), &rot[2][0]);
        PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("z_axisy"), &rot[2][1]);
        PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("z_axisz"), &rot[2][2]);
        PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("focal"), &mesh.shot.Intrinsics.FocalMm);
        PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("scalex"), &mesh.shot.Intrinsics.PixelSizeMm[0]);
        PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("scaley"), &mesh.shot.Intrinsics.PixelSizeMm[1]);
        PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("centerx"), &mesh.shot.Intrinsics.CenterPx[0]);
        PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("centery"), &mesh.shot.Intrinsics.CenterPx[1]);
        PushDescriport<int, int, 1>(cameraProp, cameraDescr, std::string("viewportx"), &mesh.shot.Intrinsics.ViewportPx[0]);
        PushDescriport<int, int, 1>(cameraProp, cameraDescr, std::string("viewporty"), &mesh.shot.Intrinsics.ViewportPx[1]);
        PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("k1"), &mesh.shot.Intrinsics.k[0]);
        PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("k2"), &mesh.shot.Intrinsics.k[1]);
        PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("k3"), &mesh.shot.Intrinsics.k[2]);
        PushDescriport<ScalarType, ScalarType, 1>(cameraProp, cameraDescr, std::string("k4"), &mesh.shot.Intrinsics.k[3]);
      }

      //Vertex
      std::vector<std::string> nameList;
      VertexType::Name(nameList);
      std::vector<PlyProperty> vertexProp;
      ElementDescriptor vertexDescr(NNP_VERTEX_ELEM);
      if (nameList.size() > 0 && mesh.vert.size() > 0)
      {
        if ((bitMask & BitMask::IO_VERTCOORD) && VertexType::HasCoord())
          PushDescriport<VertexType, VertexCoordScalar, 3>(vertexProp, vertexDescr, NNP_PXYZ, (*mesh.vert.begin()).P().V());
        if ((bitMask & BitMask::IO_VERTNORMAL) && vcg::tri::HasPerVertexNormal(mesh))
        {
          if (ocfVertexMask & BitMask::IO_VERTNORMAL)
            PushDescriport<VertexNormalType, VertexNormScalar, 3>(vertexProp, vertexDescr, NNP_NXYZ, (*mesh.vert.begin()).N().V());
          else
            PushDescriport<VertexType, VertexNormScalar, 3>(vertexProp, vertexDescr, NNP_NXYZ, (*mesh.vert.begin()).N().V());
        }
        if ((bitMask & BitMask::IO_VERTCOLOR) && vcg::tri::HasPerVertexColor(mesh))
        {
          if (ocfVertexMask & BitMask::IO_VERTCOLOR)
            PushDescriport<VertexColorType, VertexColorScalar, 4>(vertexProp, vertexDescr, NNP_CRGBA, (*mesh.vert.begin()).C().V());
          else
            PushDescriport<VertexType, VertexColorScalar, 4>(vertexProp, vertexDescr, NNP_CRGBA, (*mesh.vert.begin()).C().V());
        }
        if ((bitMask & BitMask::IO_VERTQUALITY) && vcg::tri::HasPerVertexQuality(mesh))
        {
          if (ocfVertexMask & BitMask::IO_VERTQUALITY)
            PushDescriport<VertexQuality, VertexQuality, 1>(vertexProp, vertexDescr, NNP_QUALITY, &(*mesh.vert.begin()).Q());
          else
            PushDescriport<VertexType, VertexQuality, 1>(vertexProp, vertexDescr, NNP_QUALITY, &(*mesh.vert.begin()).Q());
        }
        if ((bitMask & BitMask::IO_VERTFLAGS) && vcg::tri::HasPerVertexFlags(mesh))
          PushDescriport<VertexType, VertexFlag, 1>(vertexProp, vertexDescr, NNP_BITFLAG, &(*mesh.vert.begin()).Flags());
        if ((bitMask & BitMask::IO_VERTRADIUS) && vcg::tri::HasPerVertexRadius(mesh))
        {
          if (ocfVertexMask & BitMask::IO_VERTRADIUS)
            PushDescriport<VertexRadius, VertexRadius, 1>(vertexProp, vertexDescr, NNP_DENSITY, &(*mesh.vert.begin()).R());
          else
            PushDescriport<VertexType, VertexRadius, 1>(vertexProp, vertexDescr, NNP_DENSITY, &(*mesh.vert.begin()).R());
        }
        if ((bitMask & BitMask::IO_VERTTEXCOORD) && vcg::tri::HasPerVertexTexCoord(mesh))
        {
          if (ocfVertexMask & BitMask::IO_VERTTEXCOORD)
            PushDescriport<VertexTexCoordType, VertexTexScalar, 2>(vertexProp, vertexDescr, NNP_TEXTURE2D, (*mesh.vert.begin()).T().P().V());
          else
            PushDescriport<VertexType, VertexTexScalar, 2>(vertexProp, vertexDescr, NNP_TEXTURE2D, (*mesh.vert.begin()).T().P().V());
        }
        if ((bitMask & BitMask::IO_VERTCURV) && vcg::tri::HasPerVertexCurvature(mesh))
        {
          if (ocfVertexMask & BitMask::IO_VERTTEXCOORD)
          {
            PushDescriport<VertexCurType, VertexCurScalar, 1>(vertexProp, vertexDescr, NNP_KG, &(*mesh.vert.begin()).Kg());
            PushDescriport<VertexCurType, VertexCurScalar, 1>(vertexProp, vertexDescr, NNP_KH, &(*mesh.vert.begin()).Kh());
          }
          else
          {
            PushDescriport<VertexType, VertexCurScalar, 1>(vertexProp, vertexDescr, NNP_KG, &(*mesh.vert.begin()).Kg());
            PushDescriport<VertexType, VertexCurScalar, 1>(vertexProp, vertexDescr, NNP_KH, &(*mesh.vert.begin()).Kh());
          }
        }
        if ((bitMask & BitMask::IO_VERTCURVDIR) && vcg::tri::HasPerVertexCurvatureDir(mesh))
        {
          if (ocfVertexMask & BitMask::IO_VERTCURVDIR)
          {
            PushDescriport<VertexCurDirType, VertexDirCurScalar, 1>(vertexProp, vertexDescr, NNP_K1, &(*mesh.vert.begin()).K1());
            PushDescriport<VertexCurDirType, VertexDirCurScalar, 1>(vertexProp, vertexDescr, NNP_K2, &(*mesh.vert.begin()).K2());
            PushDescriportList<VertexCurDirType, VertexDirCurVecScalar, 3>(vertexProp, vertexDescr, NNP_K1DIR, (*mesh.vert.begin()).PD1().V());
            PushDescriportList<VertexCurDirType, VertexDirCurVecScalar, 3>(vertexProp, vertexDescr, NNP_K2DIR, (*mesh.vert.begin()).PD2().V());
          }
          else
          {
            PushDescriport<VertexType, VertexDirCurScalar, 1>(vertexProp, vertexDescr, NNP_K1, &(*mesh.vert.begin()).K1());
            PushDescriport<VertexType, VertexDirCurScalar, 1>(vertexProp, vertexDescr, NNP_K2, &(*mesh.vert.begin()).K2());
            PushDescriportList<VertexType, VertexDirCurVecScalar, 3>(vertexProp, vertexDescr, NNP_K1DIR, (*mesh.vert.begin()).PD1().V());
            PushDescriportList<VertexType, VertexDirCurVecScalar, 3>(vertexProp, vertexDescr, NNP_K2DIR, (*mesh.vert.begin()).PD2().V());
          }
        }
        if ((bitMask & BitMask::IO_VERTATTRIB))
        {
          typename std::set<PointerToAttribute>::iterator ai;
          int userSize = custom.vertexAttrib.size();
          for (ai = mesh.vert_attr.begin(); ai != mesh.vert_attr.end(); ++ai)
          {
            bool userDescr = false;
            for (int i = 0; i < userSize; i++)
            {
              if ((*custom.vertexAttrib[i]).name == (*ai)._name)
              {
                userDescr = true;
                break;               
              }
            }
            if (!userDescr)
              custom.CreateVertexAttribDescriptor(&(*ai));
          }
          for (size_t i = 0; i < custom.vertexAttrib.size(); i++)
          {
            vertexProp.push_back(custom.vertexAttribProp[i]);
            vertexDescr.dataDescriptor.push_back(custom.vertexAttrib[i]);
          }
				}
			}

			//Edge
			nameList.clear();
			EdgeType::Name(nameList);
			std::vector<PlyProperty> edgeProp;
			ElementDescriptor edgeDescr(NNP_EDGE_ELEM);
			std::vector<vcg::Point2i> edgeIndex;
			for (size_t i = 0; i < mesh.edge.size(); i++)
				edgeIndex.push_back(vcg::Point2i(vcg::tri::Index(mesh, mesh.edge[i].V(0)), vcg::tri::Index(mesh, mesh.edge[i].V(1))));
			if (nameList.size() > 0 && mesh.edge.size() > 0)
			{
				if ((bitMask & BitMask::IO_EDGEINDEX) && EdgeType::HasVertexRef())
				{
        	PushDescriport<vcg::Point2i, int, 1>(edgeProp, edgeDescr, NNP_EDGE_V1, &(*edgeIndex.begin()).V()[0]);
					PushDescriport<vcg::Point2i, int, 1>(edgeProp, edgeDescr, NNP_EDGE_V2, &(*edgeIndex.begin()).V()[1]);
				}
				if ((bitMask & BitMask::IO_EDGEQUALITY) && vcg::tri::HasPerEdgeQuality(mesh))
					PushDescriport<EdgeType, EdgeQuality, 1>(edgeProp, edgeDescr, NNP_QUALITY, &(*mesh.edge.begin()).Q());
				if ((bitMask & BitMask::IO_EDGECOLOR) && vcg::tri::HasPerEdgeColor(mesh))
					PushDescriport<EdgeType, EdgeColorScalar, 4>(edgeProp, edgeDescr, NNP_CRGBA, (*mesh.edge.begin()).C().V());
				if ((bitMask & BitMask::IO_EDGEFLAGS) && vcg::tri::HasPerEdgeFlags(mesh))
					PushDescriport<EdgeType, EdgeFlag, 1>(edgeProp, edgeDescr, NNP_BITFLAG, &(*mesh.edge.begin()).Flags());
				
        if ((bitMask & BitMask::IO_EDGEATTRIB))
        {
          typename std::set<PointerToAttribute>::iterator ai;
          int userSize = custom.edgeAttrib.size();
          for (ai = mesh.edge_attr.begin(); ai != mesh.edge_attr.end(); ++ai)
          {
            bool userDescr = false;
            for (int i = 0; i < userSize; i++)
            {
              if ((*custom.edgeAttrib[i]).name == (*ai)._name)
              {
                userDescr = true;
                break;
              }
            }
            if (!userDescr)
              custom.CreateEdgeAttribDescriptor(&(*ai));
          }
          for (size_t i = 0; i < custom.edgeAttrib.size(); i++)
          {
            edgeProp.push_back(custom.edgeAttribProp[i]);
            edgeDescr.dataDescriptor.push_back(custom.edgeAttrib[i]);
          }
        }
			}
			
			//Face
			nameList.clear();
			FaceType::Name(nameList);
			std::vector<PlyProperty> faceProp;
			ElementDescriptor faceDescr(NNP_FACE_ELEM);
			std::vector<std::vector<unsigned int>> faceIndex;
			std::vector<vcg::ndim::Point<6, FaceTexScalar>> wedgeTexCoord;
      for (size_t i = 0; i < mesh.face.size(); i++)
      {
        faceIndex.push_back(std::vector<unsigned int>());
        for (int j = 0; j < mesh.face[i].VN(); j++)
          faceIndex.back().push_back(vcg::tri::Index(mesh, mesh.face[i].V(j)));
      }
			if (((bitMask & BitMask::IO_WEDGTEXCOORD) || (bitMask & BitMask::IO_WEDGTEXMULTI)) && vcg::tri::HasPerWedgeTexCoord(mesh))
			{
				for (size_t i = 0; i < mesh.face.size(); i++)
				{
          wedgeTexCoord.push_back(vcg::ndim::Point<6, FaceTexScalar>());
					for (int j = 0; j < 3; j++)
					{
						wedgeTexCoord.back()[j * 2] = mesh.face[i].WT(j).U();
						wedgeTexCoord.back()[j * 2 + 1] = mesh.face[i].WT(j).V();
					}
				}
			}
			if (nameList.size() > 0 && mesh.face.size() > 0)
			{
				if ((bitMask & BitMask::IO_FACEINDEX) && FaceType::HasVertexRef())
          PushDescriportList<std::vector<unsigned int>, unsigned int, 0>(faceProp, faceDescr, NNP_FACE_VERTEX_LIST, &faceIndex[0]);
				if ((bitMask & BitMask::IO_FACEFLAGS) && vcg::tri::HasPerFaceFlags(mesh))
					PushDescriport<FaceType, FaceFlag, 1>(faceProp, faceDescr, NNP_BITFLAG, &(*mesh.face.begin()).Flags());
        if ((bitMask & BitMask::IO_FACECOLOR) && vcg::tri::HasPerFaceColor(mesh))
        {
          if (ocfFaceMask & BitMask::IO_FACECOLOR)
            PushDescriport<FaceColorType, FaceColorScalar, 4>(faceProp, faceDescr, NNP_CRGBA, (*mesh.face.begin()).C().V());
          else
            PushDescriport<FaceType, FaceColorScalar, 4>(faceProp, faceDescr, NNP_CRGBA, (*mesh.face.begin()).C().V());
        }
        if ((bitMask & BitMask::IO_FACEQUALITY) && vcg::tri::HasPerFaceQuality(mesh))
        {
          if (ocfFaceMask & BitMask::IO_FACEQUALITY)
            PushDescriport<FaceQuality, FaceQuality, 1>(faceProp, faceDescr, NNP_QUALITY, &(*mesh.face.begin()).Q());
          else
            PushDescriport<FaceType, FaceQuality, 1>(faceProp, faceDescr, NNP_QUALITY, &(*mesh.face.begin()).Q());
        }
        if ((bitMask & BitMask::IO_FACENORMAL) && vcg::tri::HasPerFaceNormal(mesh))
        {
          if (ocfFaceMask & BitMask::IO_FACENORMAL)
            PushDescriport<FaceNormalType, FaceNormScalar, 3>(faceProp, faceDescr, NNP_NXYZ, (*mesh.face.begin()).N().V());
          else
            PushDescriport<FaceType, FaceNormScalar, 3>(faceProp, faceDescr, NNP_NXYZ, (*mesh.face.begin()).N().V());
        }
				if ((bitMask & BitMask::IO_VERTCURVDIR) && vcg::tri::HasPerFaceCurvatureDir(mesh))
				{
          if (ocfFaceMask & BitMask::IO_VERTCURVDIR)
          {
            PushDescriport<FaceCurDirType, FaceDirCurScalar, 1>(faceProp, faceDescr, NNP_K1, &(*mesh.face.begin()).K1());
            PushDescriport<FaceCurDirType, FaceDirCurScalar, 1>(faceProp, faceDescr, NNP_K2, &(*mesh.face.begin()).K2());
            PushDescriportList<FaceCurDirType, FaceDirCurVecScalar, 3>(faceProp, faceDescr, NNP_K1DIR, (*mesh.face.begin()).PD1().V());
            PushDescriportList<FaceCurDirType, FaceDirCurVecScalar, 3>(faceProp, faceDescr, NNP_K2DIR, (*mesh.face.begin()).PD2().V());
          }
          else
          {
            PushDescriport<FaceType, FaceDirCurScalar, 1>(faceProp, faceDescr, NNP_K1, &(*mesh.face.begin()).K1());
            PushDescriport<FaceType, FaceDirCurScalar, 1>(faceProp, faceDescr, NNP_K2, &(*mesh.face.begin()).K2());
            PushDescriportList<FaceType, FaceDirCurVecScalar, 3>(faceProp, faceDescr, NNP_K1DIR, (*mesh.face.begin()).PD1().V());
            PushDescriportList<FaceType, FaceDirCurVecScalar, 3>(faceProp, faceDescr, NNP_K2DIR, (*mesh.face.begin()).PD2().V());
          }
				}
				if (((bitMask & BitMask::IO_WEDGTEXCOORD) || (bitMask & BitMask::IO_WEDGTEXMULTI)) && vcg::tri::HasPerWedgeTexCoord(mesh))
          PushDescriportList<vcg::ndim::Point<6, FaceTexScalar>, FaceTexScalar, 6>(faceProp, faceDescr, NNP_FACE_WEDGE_TEX, (*wedgeTexCoord.begin()).V());
        if ((bitMask & BitMask::IO_WEDGTEXMULTI) && vcg::tri::HasPerWedgeTexCoord(mesh))
        {
          if (ocfFaceMask & BitMask::IO_VERTCURVDIR)
            PushDescriport<FaceTexCoordType, short, 1>(faceProp, faceDescr, NNP_TEXTUREINDEX, &(*mesh.face.begin()).WT(0).N());
          else
            PushDescriport<FaceType, short, 1>(faceProp, faceDescr, NNP_TEXTUREINDEX, &(*mesh.face.begin()).WT(0).N());
        }
        if ((bitMask & BitMask::IO_WEDGCOLOR) && vcg::tri::HasPerWedgeColor(mesh))
        {
          if (ocfFaceMask & BitMask::IO_WEDGCOLOR)
            PushDescriportList<WedgeColorType, WedgeColorScalar, 12>(faceProp, faceDescr, NNP_FACE_WEDGE_COLOR, (*mesh.face.begin()).WC(0).V());
          else
            PushDescriportList<FaceType, WedgeColorScalar, 12>(faceProp, faceDescr, NNP_FACE_WEDGE_COLOR, (*mesh.face.begin()).WC(0).V());
        }
        if ((bitMask & BitMask::IO_WEDGNORMAL) && vcg::tri::HasPerWedgeNormal(mesh))
        {
          if (ocfFaceMask & BitMask::IO_WEDGNORMAL)
            PushDescriportList<WedgeNormalType, WedgeNormalScalar, 9>(faceProp, faceDescr, NNP_FACE_WEDGE_NORMAL, (*mesh.face.begin()).WN(0).V());
          else
            PushDescriportList<FaceType, WedgeNormalScalar, 9>(faceProp, faceDescr, NNP_FACE_WEDGE_NORMAL, (*mesh.face.begin()).WN(0).V());
        }
        if ((bitMask & BitMask::IO_FACEATTRIB))
        {
          typename std::set<PointerToAttribute>::iterator ai;
          int userSize = custom.faceAttrib.size();
          for (ai = mesh.face_attr.begin(); ai != mesh.face_attr.end(); ++ai)
          {
            bool userDescr = false;
            for (int i = 0; i < userSize; i++)
            {
              if ((*custom.faceAttrib[i]).name == (*ai)._name)
              {
                userDescr = true;
                break;
              }
            }
            if (!userDescr)
              custom.CreateFaceAttribDescriptor(&(*ai));
          }
          for (size_t i = 0; i < custom.faceAttrib.size(); i++)
          {
            faceProp.push_back(custom.faceAttribProp[i]);
            faceDescr.dataDescriptor.push_back(custom.faceAttrib[i]);
          }
        }
			}
			
			Info infoSave;
			infoSave.filename = filename;
			infoSave.binary = binary;
			PlyElement cameraElem(std::string("camera"), cameraProp, 1);
			PlyElement vertexElem(NNP_VERTEX_ELEM, vertexProp, mesh.vert.size());
			PlyElement edgeElem(NNP_EDGE_ELEM, edgeProp, mesh.edge.size());
			PlyElement faceElem(NNP_FACE_ELEM, faceProp, mesh.face.size());
			infoSave.AddPlyElement(cameraElem);
			infoSave.AddPlyElement(vertexElem);
			infoSave.AddPlyElement(edgeElem);
			infoSave.AddPlyElement(faceElem);
			infoSave.textureFile = mesh.textures;
			std::vector<ElementDescriptor*> meshDescr;
			meshDescr.push_back(&cameraDescr);
			meshDescr.push_back(&vertexDescr);
			meshDescr.push_back(&edgeDescr);
			meshDescr.push_back(&faceDescr);

			//Mesh attribute
			if ((bitMask & BitMask::IO_MESHATTRIB))
			{
				typename CustomAttributeDescriptor::MapMeshAttribIter iter = custom.meshAttrib.begin();
				typename CustomAttributeDescriptor::MapMeshAttribPropIter iterProp = custom.meshAttribProp.begin();
				for (; iter != custom.meshAttrib.end(); iter++, iterProp++)
				{
					std::string name((*iter).first);
					PlyElement customElem(name, (*iterProp).second, custom.meshAttribCnt[(*iter).first]);
					infoSave.AddPlyElement(customElem);
					meshDescr.push_back(new ElementDescriptor(name));
					meshDescr.back()->dataDescriptor = (*iter).second;
				}
			}
			
			bool flag = nanoply::SaveModel(infoSave.filename, meshDescr, infoSave);

			for (size_t i = 0; i < cameraDescr.dataDescriptor.size(); i++)
				delete cameraDescr.dataDescriptor[i];
			for (size_t i = 0; i < vertexDescr.dataDescriptor.size(); i++)
				if (vertexDescr.dataDescriptor[i]->elem != NNP_UNKNOWN_ENTITY)
					delete vertexDescr.dataDescriptor[i];
			for (size_t i = 0; i < edgeDescr.dataDescriptor.size(); i++)
				if (edgeDescr.dataDescriptor[i]->elem != NNP_UNKNOWN_ENTITY)
					delete edgeDescr.dataDescriptor[i];
			for (size_t i = 0; i < faceDescr.dataDescriptor.size(); i++)
				if (faceDescr.dataDescriptor[i]->elem != NNP_UNKNOWN_ENTITY)
					delete faceDescr.dataDescriptor[i];
			
			return flag;
		}
		

		static bool SaveModel(const char* filename, MeshType& mesh, unsigned int bitMask, bool binary)
		{
			CustomAttributeDescriptor custom;
			return SaveModel(filename, mesh, bitMask, custom, binary);
		}


	};

}

#endif
