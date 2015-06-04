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

#ifndef __VCG_GL_MESH_ATTRIBUTES_FEEDER
#define __VCG_GL_MESH_ATTRIBUTES_FEEDER

#include <queue>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>

//#include <GL/glew.h>
#include <wrap/gl/space.h>
#include <wrap/gl/math.h>
#include <vcg/space/color4.h>


namespace vcg 
{
//WARNING: All the classes derived from MemoryInfo has been intended to be instantiated as a singleton in the host application 
	//(i.e. in every application using it just an instance of a class derived from MemoryInfo should be declared).

	class MemoryInfo
	{
	public:
		class MemoryInfoException : public std::exception
		{
		public:
			MemoryInfoException(const char* text)
				:std::exception(),_exctext(text){}

			~MemoryInfoException() throw() {}
			inline const char* what() const throw() {return _exctext;}
		private:
			const char* _exctext;
		};

		MemoryInfo(long long unsigned int originalmem)
			:_originaltotalmemory(originalmem),_currentfreememory(_originaltotalmemory)
		{       
		}

		virtual ~MemoryInfo() {}
		virtual void acquiredMemory(long long unsigned int mem) = 0;
		virtual long long unsigned int usedMemory() const = 0;
		virtual long long unsigned int currentFreeMemory() const = 0;
		virtual void releasedMemory(long long unsigned int mem = 0) = 0;
		virtual bool isAdditionalMemoryAvailable(long long unsigned int mem) = 0;

	protected:
		const long long unsigned int _originaltotalmemory;
		long long unsigned int _currentfreememory;
	};

	//WARNING: this is not a thread safe class. The object derived from MemoryInfo are intended to be used inside GLMeshAttributeFeeder as static variable in order to manage the available GPUMemory.
	//We strongly recommend you to define in your code a thread safe version of the class, defining mutexed access member functions. 
	//This class should be consider just as a basic example for the implementations of the required functionalities. 
	//It is safe to use it just when the user has only one mesh to pass to the GPU.

	class NotThreadSafeMemoryInfo : public MemoryInfo
	{
	public:
		NotThreadSafeMemoryInfo(long long unsigned int originalmem)
			:MemoryInfo(originalmem)
		{
		}

		~NotThreadSafeMemoryInfo() {}

		void acquiredMemory(long long unsigned int mem)
		{
			if (mem > _originaltotalmemory)
				throw MemoryInfo::MemoryInfoException("It has been requested more memory than the total one.\\n");
			else 
				if (mem > _currentfreememory)
					throw MemoryInfo::MemoryInfoException("It has been requested more memory than the free available one.\\n");
				else
					_currentfreememory -= mem;
		}

		long long unsigned int usedMemory() const
		{
			return _originaltotalmemory - _currentfreememory;
		}

		long long unsigned int currentFreeMemory() const
		{
			return _currentfreememory;
		}

		void releasedMemory(long long unsigned int mem = 0)
		{
			if (mem > _originaltotalmemory)
				throw MemoryInfo::MemoryInfoException("It has been released more memory than the total one. Something strange happened!\\n");
			else
				_currentfreememory += mem;
		}

		bool isAdditionalMemoryAvailable(long long unsigned int mem)
		{
			return (_currentfreememory >= mem);
		}
	};

	struct GLFeedEnum
	{
		enum MESH_ATTRIBUTE_UPDATED
		{
			ATT_VERTPOSITION   = 0x00000000,
			ATT_VERTNORMAL     = 0x00000001,
			ATT_VERTCOLOR      = 0x00000002,
			ATT_VERTTEXTURE    = 0x00000004,
			ATT_FACENORMAL     = 0x00000008,
			ATT_FACECOLOR      = 0x00000010,
			ATT_WEDGETEXTURE   = 0x00000020,
			ATT_ALL            = 0xffffffff 
		};

		enum COLOR_MODALITY
		{
			CL_NONE = 0,
			CL_PERVERT,
			CL_PERFACE,
			CL_PERMESH
		};

		enum NORMAL_MODALITY
		{
			NR_NONE = 0,
			NR_PERVERT,
			NR_PERFACE
		};

		enum TEXTURE_MODALITY
		{
			TX_NONE = 0,
			TX_PERVERT,
			TX_PERWEDGE
		};
	};

	//WARNING! member functions of this class should be called by the host application using concurrency  
	template <typename MESHTYPE>
	class GLMeshAttributesFeeder : public GLFeedEnum
	{
	public:
		GLMeshAttributesFeeder(/*const*/ MESHTYPE& mesh,MemoryInfo& meminfo, size_t perbatchtriangles = 100)
            :_mesh(mesh),_gpumeminfo(meminfo),_bo(8,NULL),_vaohandle(0),_lastfeedingusedreplicatedpipeline(false),_perbatchsimplex(perbatchtriangles),_chunkmap(),_borendering(false)
		{
			_bo[VERTPOSITIONBO] = new GLBufferObject(3,GL_FLOAT);
			_bo[VERTNORMALBO] = new GLBufferObject(3,GL_FLOAT);    
			_bo[FACENORMALBO] = new GLBufferObject(3,GL_FLOAT);
			_bo[VERTCOLORBO] = new GLBufferObject(4,GL_UNSIGNED_BYTE);
			_bo[FACECOLORBO] = new GLBufferObject(4,GL_UNSIGNED_BYTE);
			_bo[VERTTEXTUREBO] = new GLBufferObject(2,GL_FLOAT);
			_bo[WEDGETEXTUREBO] = new GLBufferObject(2,GL_FLOAT);
			_bo[VERTINDEXBO] = new GLBufferObject(3,GL_INT);
		}

		~GLMeshAttributesFeeder() 
		{
			size_t ii = 0;
			for(typename std::vector<GLBufferObject*>::iterator it = _bo.begin();it != _bo.end();++it)
			{
				if (*it != NULL)
				{
					glDeleteBuffers(1,&((*it)->_bohandle));
					(*it)->_bohandle = 0;
				}
				delete _bo[BO_NAMES(ii)];
				++ii;
			}
			_bo.clear();
			glDeleteVertexArrays(1,&_vaohandle);
		}

		void setPerBatchTriangles(size_t perbatchtriangles)
		{
			_perbatchsimplex = perbatchtriangles;
		}

		size_t perBatchTriangles() const
		{
			return _perbatchsimplex;
		}

		bool renderedWithBO() const
		{
			return _borendering;
		}

		void update(int mask)
		{  
			if ((mask & ATT_VERTPOSITION) || (mask & ATT_ALL))
				_bo[VERTPOSITIONBO]->_isvalid = false;

			if ((mask & ATT_VERTNORMAL) || (mask & ATT_ALL))
				_bo[VERTNORMALBO]->_isvalid = false;

			if ((mask & ATT_FACENORMAL) || (mask & ATT_ALL))
				_bo[FACENORMALBO]->_isvalid = false;

			if ((mask & ATT_VERTCOLOR) || (mask & ATT_ALL))
				_bo[VERTCOLORBO]->_isvalid = false;

			if ((mask & ATT_FACECOLOR) || (mask & ATT_ALL))
				_bo[FACECOLORBO]->_isvalid = false;

			if ((mask & ATT_VERTTEXTURE) || (mask & ATT_ALL))
				_bo[VERTTEXTUREBO]->_isvalid = false;

			if ((mask & ATT_WEDGETEXTURE) || (mask & ATT_ALL))
				_bo[WEDGETEXTUREBO]->_isvalid = false;

			if (mask & ATT_ALL)
				_bo[VERTINDEXBO]->_isvalid = false;
		}

		/*WARNING: the passTrianglesToOpenGL & the passPointsToOpenGL functions should be invoked with a reading mutex, 
		in order to be sure that the referenced mesh would not been changed when the mesh attributes are updated into the buffer objects.
		textureindex contained the texture OpenGL ids loaded in the gpu memory by the main application*/ 
		void passTrianglesToOpenGL(NORMAL_MODALITY nm,COLOR_MODALITY cm,TEXTURE_MODALITY tm,const std::vector<GLuint>& textureindex = std::vector<GLuint>())
		{
			std::vector<bool> importattribute(_bo.size());
			std::vector<bool> attributestobeupdated;
			attributesToBeImportedInTriangleBasedPipeline(importattribute, nm, cm, tm);
			bool replicated = !importattribute[VERTINDEXBO];
			bool immediatemode = !(buffersAllocationFunction(nm,cm,tm,importattribute,attributestobeupdated));
			if (immediatemode)
				immediateMode(nm,cm,tm,textureindex);
			else
			{
				bool somethingtoupdate = false;
				for(size_t hh = 0;hh < attributestobeupdated.size();++hh)
					somethingtoupdate = somethingtoupdate || attributestobeupdated[hh];
				if (somethingtoupdate) 
				{
					if (replicated)
					{
						//WARNING!In case we have to update the wedgetexture bo maybe (not always!) we must update also the other buffer already in memory
						//cause the wedgetexture pipeline force a change in the order of the triangles in GPU. 
						//they are now ordered by the texture seam and not more by the triangle index!
						if (attributestobeupdated[WEDGETEXTUREBO])
						{
							for(size_t jj = 0;jj < attributestobeupdated.size();++jj)
								attributestobeupdated[jj] = importattribute[jj] || attributestobeupdated[jj];
						}
						updateBuffersReplicatedPipeline(attributestobeupdated,tm);
					}
					else
						updateBuffersIndexedPipeline(attributestobeupdated);
				}
				drawTriangles(nm,cm,tm,textureindex);
			}

		}

		void passPointsToOpenGL(NORMAL_MODALITY nm,COLOR_MODALITY cm)
		{
			std::vector<bool> importattribute(_bo.size());
			std::vector<bool> attributestobeupdated;
			attributesToBeImportedInPointBasedPipeline(importattribute, nm, cm);
			GLenum err = glGetError();
			assert(err == GL_NO_ERROR);
			bool immediatemode = !(buffersAllocationFunction(nm,cm,TX_NONE,importattribute,attributestobeupdated));
			err = glGetError();
			assert(err == GL_NO_ERROR);
			std::vector<GLuint> textureindex;
			if (immediatemode)
				immediateMode(nm,cm,TX_NONE,textureindex);
			else
			{
				bool somethingtoupdate = false;
				for(size_t hh = 0;hh < attributestobeupdated.size();++hh)
					somethingtoupdate = somethingtoupdate || attributestobeupdated[hh];
				if (somethingtoupdate) 
				{
					GLenum err = glGetError();
					assert(err == GL_NO_ERROR);
					updateBuffersIndexedPipeline(attributestobeupdated);
					err = glGetError();
					assert(err == GL_NO_ERROR);
				}
				err = glGetError();
				assert(err == GL_NO_ERROR);
				drawPoints(nm,cm);
				err = glGetError();
				assert(err == GL_NO_ERROR);
			}
		}

	private:
		struct GLBufferObject
		{
			GLBufferObject(size_t components,GLenum gltype)
				:_size(0),_components(components),_isvalid(false),_gltype(gltype),_bohandle(0)
			{
			}

			size_t getSizeOfGLType() const
			{
				switch(_gltype)
				{
				case(GL_FLOAT):
					return sizeof(GLfloat);
				case(GL_INT):
					return sizeof(GLint);
				case(GL_UNSIGNED_INT):
					return sizeof(GLuint);
				case(GL_UNSIGNED_BYTE):
					return sizeof(GLubyte);
				}
				return 0;
			}

			size_t _size;
			const size_t _components;
			bool _isvalid;
			const GLenum _gltype;
			GLuint _bohandle;
		};

		enum BO_NAMES
		{
			VERTPOSITIONBO = 0,
			VERTNORMALBO = 1,
			FACENORMALBO = 2,
			VERTCOLORBO = 3,
			FACECOLORBO = 4,
			VERTTEXTUREBO = 5,
			WEDGETEXTUREBO = 6,
			VERTINDEXBO = 7
		};

		long long unsigned int bufferObjectsMemoryRequired(NORMAL_MODALITY nm,COLOR_MODALITY cm,TEXTURE_MODALITY tm,bool generateindex) const
		{
			bool replicated = (nm == NR_PERFACE) || (cm == CL_PERFACE) || (tm == TX_PERWEDGE);
			long long unsigned int result(0);

			result += (long long unsigned int) boExpectedDimension(VERTPOSITIONBO,replicated,generateindex);
			if (nm == NR_PERVERT)
				result += (long long unsigned int) boExpectedDimension(VERTNORMALBO,replicated,generateindex);
			if (nm == NR_PERFACE)
				result += (long long unsigned int) boExpectedDimension(FACENORMALBO,replicated,generateindex);
			if (cm == CL_PERVERT)
				result += (long long unsigned int) boExpectedDimension(VERTCOLORBO,replicated,generateindex);
			if (cm == CL_PERFACE)
				result += (long long unsigned int) boExpectedDimension(FACECOLORBO,replicated,generateindex);
			if (tm == TX_PERVERT)
				result += (long long unsigned int) boExpectedDimension(VERTTEXTUREBO,replicated,generateindex);
			if (tm == TX_PERWEDGE)
				result += (long long unsigned int) boExpectedDimension(WEDGETEXTUREBO,replicated,generateindex);
			result += (long long unsigned int) boExpectedDimension(VERTINDEXBO,replicated,generateindex);
			return result;
		}

		void attributesToBeImportedInTriangleBasedPipeline( std::vector<bool> &importattribute, NORMAL_MODALITY nm, COLOR_MODALITY cm, TEXTURE_MODALITY tm )
		{
			importattribute[VERTPOSITIONBO] = true;
			importattribute[VERTNORMALBO] = vcg::tri::HasPerVertexNormal(_mesh) && (nm == NR_PERVERT);
			importattribute[FACENORMALBO]  = vcg::tri::HasPerFaceNormal(_mesh) && (nm == NR_PERFACE);
			importattribute[VERTCOLORBO]  = vcg::tri::HasPerVertexColor(_mesh) && (cm == CL_PERVERT);
			importattribute[FACECOLORBO]  = vcg::tri::HasPerFaceColor(_mesh) && (cm == CL_PERFACE);
			importattribute[VERTTEXTUREBO]  = vcg::tri::HasPerVertexTexCoord(_mesh) && (tm == TX_PERVERT);
			importattribute[WEDGETEXTUREBO] = vcg::tri::HasPerWedgeTexCoord(_mesh) && (tm == TX_PERWEDGE);

			/*check if the requested attributes are available inside the mesh*/
			if (!importattribute[VERTNORMALBO] && !importattribute[FACENORMALBO]) 
				nm = NR_NONE;

			if (!importattribute[VERTCOLORBO] && !importattribute[FACECOLORBO])
				cm = CL_NONE;

			if (!importattribute[VERTTEXTUREBO] && !importattribute[WEDGETEXTUREBO])
				tm = TX_NONE;

			bool replicated = (nm == NR_PERFACE) || (cm == CL_PERFACE) || (tm == TX_PERWEDGE);

			importattribute[VERTINDEXBO] = !replicated; 
		}

		void attributesToBeImportedInPointBasedPipeline( std::vector<bool> &importattribute, NORMAL_MODALITY nm, COLOR_MODALITY cm)
		{
			importattribute[VERTPOSITIONBO] = true;
			importattribute[VERTNORMALBO] = vcg::tri::HasPerVertexNormal(_mesh) && (nm == NR_PERVERT);
			importattribute[FACENORMALBO]  = false;
			importattribute[VERTCOLORBO]  = vcg::tri::HasPerVertexColor(_mesh) && (cm == CL_PERVERT);
			importattribute[FACECOLORBO]  = false;
			importattribute[VERTTEXTUREBO]  = false;
			importattribute[WEDGETEXTUREBO] = false;
			importattribute[VERTINDEXBO] = false; 
		}

		bool buffersAllocationFunction(NORMAL_MODALITY nm,COLOR_MODALITY cm,TEXTURE_MODALITY tm,const std::vector<bool>& importattribute,std::vector<bool>& attributestobeupdated)
		{
            if (_vaohandle == 0)
				glGenVertexArrays(1,&_vaohandle);
			bool replicated = (importattribute[FACENORMALBO] || importattribute[FACECOLORBO] || importattribute[WEDGETEXTUREBO]);
			attributestobeupdated.clear();
			attributestobeupdated.resize(importattribute.size());
			long long unsigned int bomemoryrequiredbymesh = this->bufferObjectsMemoryRequired(nm,cm,tm,importattribute[VERTINDEXBO]);
			unsigned int ii = 0;
			for(typename std::vector<GLBufferObject*>::iterator it = _bo.begin();it != _bo.end();++it)
			{
				BO_NAMES boname = static_cast<BO_NAMES>(ii);
				size_t sz = boExpectedSize(boname,replicated,importattribute[VERTINDEXBO]);
				size_t dim = boExpectedDimension(boname,replicated,importattribute[VERTINDEXBO]); 
				if ((*it != NULL) && ((!(*it)->_isvalid) || 
					/*if _lastfeedingusedreplicatedpipeline == false means that maybe there are valid per vertex attribute buffer objects that MUST be reallocated anyway cause we have to switch to the replicated attributes pipeline*/  
					(replicated && !_lastfeedingusedreplicatedpipeline && (isPerVertexAttribute(boname) || (boname == VERTINDEXBO))) ||
					/*we switched back from the replicated pipeline to the normal one. All the bos have to be regenerated*/
					(!replicated && _lastfeedingusedreplicatedpipeline) ||
					/*the buffer object is valid but for same reason the number of cells of the bo don't suit anymore the required size. we have to reallocate the buffer object*/
					(((*it)->_isvalid) && (sz != (*it)->_size))) ||
					//the buffer is valid, but the attribute is not required to be displayed
					(((*it)->_isvalid) && !importattribute[ii]))
				{

					//disableClientState(boname,importattribute);

					/*WARNING! THIS CODE MUST BE INCAPSULATED INTO A DEALLOCATE FUNCTION IN A PROPER MADE BUFFER OBJECT CLASS
					I DON'T INSERT IT INTO THE GLBufferObjectInfo CLASS CAUSE I CONSIDER IT A TEMPORARY PRIVATE STRUCT*/ 
					glDeleteBuffers(1,&((*it)->_bohandle));
					(*it)->_bohandle = 0;
					if ((*it)->_size > 0)
						//we don't use dim cause dim is the value that is going to be allocated, instead use (*it)->_size * (*it)->getSizeOfGLType() is the value already in the buffer 
						_gpumeminfo.releasedMemory((*it)->_size * (*it)->getSizeOfGLType());
					(*it)->_isvalid = false;
					(*it)->_size = 0;
					/**********************************************************************************************************/
				}
				//there are already valid mesh attributes properly allocated in memory, we don't need to allocate them again
				if ((*it != NULL) && (((*it)->_isvalid) && (sz == (*it)->_size)))
					bomemoryrequiredbymesh -= dim;
				++ii;
			}

			if (!_gpumeminfo.isAdditionalMemoryAvailable(bomemoryrequiredbymesh))
			{
				std::cout << "no additional memory available!!! memory required: " << bomemoryrequiredbymesh << std::endl; 
				ii = 0;
				for(typename std::vector<GLBufferObject*>::iterator it = _bo.begin();it != _bo.end();++it)
				{  
					BO_NAMES boname = static_cast<BO_NAMES>(ii);
					size_t sz(boExpectedSize(boname,replicated,importattribute[VERTINDEXBO]));
					//there are already valid mesh attributes properly allocated in memory but there is not enough gpu memory for the remaining mesh.
					//we have to deallocate the previously allocated mesh attributes
					if ((*it != NULL) && (((*it)->_isvalid) && (sz == (*it)->_size)))
					{
						long long unsigned int dim(boExpectedDimension(boname,replicated,importattribute[VERTINDEXBO]));
						//disableClientState(boname,importattribute);
						if ((*it)->_size > 0)
						{


							/*WARNING! THIS CODE MUST BE INCAPSULATED INTO A DEALLOCATE FUNCTION IN A PROPER MADE BUFFER OBJECT CLASS
							I DON'T INSERT IT INTO THE GLBufferObjectInfo CLASS CAUSE I CONSIDER IT A TEMPORARY PRIVATE STRUCT*/  
							glDeleteBuffers(1,&(*it)->_bohandle);
							(*it)->_bohandle = 0;
							_gpumeminfo.releasedMemory(dim);
						}
						(*it)->_isvalid = false;
						(*it)->_size = 0;
						/*********************************************************************************************************/
					}
					++ii;
				}
				_borendering = false;
				_lastfeedingusedreplicatedpipeline = false;
				return false;
			}
			else
			{
				unsigned int ii = 0;
				//I have to update the invalid buffers requested to be imported
				attributestobeupdated = importattribute;
				for(typename std::vector<GLBufferObject*>::iterator it = _bo.begin();it != _bo.end();++it)
				{  
					BO_NAMES boname = static_cast<BO_NAMES>(ii);
					GLBufferObject* cbo = _bo.at(boname);
					bool importatt = importattribute.at(boname);
					GLenum target = GL_ARRAY_BUFFER;
					if (boname == VERTINDEXBO)
						target = GL_ELEMENT_ARRAY_BUFFER;

					glBindVertexArray(_vaohandle);
					bool notvalidbuttoberegenerated = (cbo != NULL) && (!cbo->_isvalid) && (importatt);
					if (notvalidbuttoberegenerated)
					{
						cbo->_size = boExpectedSize(boname,replicated,importattribute[VERTINDEXBO]);
						long long unsigned int dim = boExpectedDimension(boname,replicated,importattribute[VERTINDEXBO]);
						glGenBuffers(1, &cbo->_bohandle);
						glBindBuffer(target, cbo->_bohandle);
						glBufferData(target, dim, NULL, GL_STATIC_DRAW);
						setBufferPointerEnableClientState(boname);
						glBindBuffer(target, 0);
						_gpumeminfo.acquiredMemory(dim);
						attributestobeupdated[boname] = true;
						cbo->_isvalid = true;
					}
					else
					{
						attributestobeupdated[boname] = false;
						if (!importatt)
							cbo->_isvalid = false;
					}

					if ((cbo == NULL) || (!cbo->_isvalid))
					{
						//glBindVertexArray(_vaohandle);
						disableClientState(boname,importattribute);
						//glBindVertexArray(0);
					}
					++ii;
					glBindVertexArray(0);
				}
				_borendering = true;
				_lastfeedingusedreplicatedpipeline = replicated;
				return true;	
			}
		}

		bool updateBuffersIndexedPipeline(const std::vector<bool>& attributestobeupdated)
		{
			_chunkmap.clear();
			size_t vn = _mesh.vn;
			size_t tn = _mesh.fn;

			size_t facechunk = std::min(size_t(tn),_perbatchsimplex);
			size_t vertexchunk = std::min(size_t(vn),_perbatchsimplex);

			std::vector<vcg::Point3f> pv(vertexchunk);
			std::vector<vcg::Point3f> nv(vertexchunk);
			std::vector<vcg::Color4b> cv(vertexchunk); // Per vertex Colors
			std::vector<float> tv(vertexchunk * 2);

			size_t chunckingpu = 0;

			for(size_t i=0;i<vn;++i)
			{
				size_t chunckindex = i % vertexchunk;
				if (attributestobeupdated[VERTPOSITIONBO])
					pv[chunckindex].Import(_mesh.vert[i].cP());

				if (attributestobeupdated[VERTNORMALBO])
				{
					nv[chunckindex].Import(_mesh.vert[i].cN());
					nv[chunckindex].Normalize();
				}

				if (attributestobeupdated[VERTCOLORBO])
					cv[chunckindex] = _mesh.vert[i].cC();
				if (attributestobeupdated[VERTTEXTUREBO])
				{
					tv[chunckindex*2+0] = _mesh.vert[i].cT().U();
					tv[chunckindex*2+1] = _mesh.vert[i].cT().V();
				}

				if((i == vn - 1) || (chunckindex == vertexchunk - 1))
				{
					size_t chuncksize = vertexchunk;
					if (i == vn - 1)
						chuncksize = chunckindex + 1;

					if (attributestobeupdated[VERTPOSITIONBO])
					{
						GLBufferObject* buffobj = _bo[VERTPOSITIONBO];
						glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
						glBufferSubData(GL_ARRAY_BUFFER,chunckingpu * vertexchunk * buffobj->_components * buffobj->getSizeOfGLType(),buffobj->_components * buffobj->getSizeOfGLType() * chuncksize,&pv[0]);
						glBindBuffer(GL_ARRAY_BUFFER, 0);
					}
					if (attributestobeupdated[VERTNORMALBO])
					{
						GLBufferObject* buffobj = _bo[VERTNORMALBO];
						glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
						glBufferSubData(GL_ARRAY_BUFFER,chunckingpu * vertexchunk * buffobj->_components * buffobj->getSizeOfGLType(),buffobj->_components * buffobj->getSizeOfGLType() * chuncksize,&nv[0]);
						glBindBuffer(GL_ARRAY_BUFFER, 0);
					}
					if (attributestobeupdated[VERTCOLORBO])
					{
						GLBufferObject* buffobj = _bo[VERTCOLORBO];
						glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
						glBufferSubData(GL_ARRAY_BUFFER,chunckingpu * vertexchunk * buffobj->_components * buffobj->getSizeOfGLType(),buffobj->_components * buffobj->getSizeOfGLType() * chuncksize,&cv[0]);
						glBindBuffer(GL_ARRAY_BUFFER, 0);
					}
					if (attributestobeupdated[VERTTEXTUREBO])
					{
						GLBufferObject* buffobj = _bo[VERTTEXTUREBO];
						glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
						glBufferSubData(GL_ARRAY_BUFFER,chunckingpu * vertexchunk * buffobj->_components * buffobj->getSizeOfGLType(),buffobj->_components * buffobj->getSizeOfGLType() * chuncksize,&tv[0]);
						glBindBuffer(GL_ARRAY_BUFFER, 0);
					}

					++chunckingpu;
				}
			}

			pv.clear();
			nv.clear();
			cv.clear();
			tv.clear();

			chunckingpu = 0;
			std::vector<GLuint> ti(facechunk * 3);
			for(size_t i=0;i<tn;++i)
			{
				size_t chunckindex = i % facechunk;

				ti[chunckindex * 3 + 0] = GLuint(vcg::tri::Index(_mesh,_mesh.face[i].V(0)));
				ti[chunckindex * 3 + 1] = GLuint(vcg::tri::Index(_mesh,_mesh.face[i].V(1)));
				ti[chunckindex * 3 + 2] = GLuint(vcg::tri::Index(_mesh,_mesh.face[i].V(2)));

				if((i == tn - 1) || (chunckindex == facechunk - 1))
				{
					size_t chunksize = facechunk;
					if (i == tn - 1)
						chunksize = chunckindex + 1;

					if (attributestobeupdated[VERTINDEXBO])
					{
						glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _bo[VERTINDEXBO]->_bohandle);
						glBufferSubData(GL_ELEMENT_ARRAY_BUFFER,chunckingpu * facechunk *  _bo[VERTINDEXBO]->_components *  _bo[VERTINDEXBO]->getSizeOfGLType(),_bo[VERTINDEXBO]->_components *  _bo[VERTINDEXBO]->getSizeOfGLType() * chunksize,&ti[0]);
						glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
					}
					++chunckingpu;
				}
			}
			return true;
		}

		bool updateBuffersReplicatedPipeline(const std::vector<bool>& attributestobeupdated,GLFeedEnum::TEXTURE_MODALITY currtextmod)
		{
			size_t vn = _mesh.vn;
			size_t tn = _mesh.fn;

			size_t facechunk = std::min(size_t(tn),_perbatchsimplex);

			std::vector<vcg::Point3f> rpv(facechunk * 3);
			std::vector<vcg::Point3f> rnv(facechunk * 3);
			std::vector<vcg::Color4b> rcv(facechunk * 3); 
			std::vector<float> rtv(facechunk * 3 * 2);

			size_t chunckingpu = 0;

			//it's a map containing for each texture seams n a vector of all the triangle index ranges having n has texture seam
			//Suppose that in a mesh we have
			//TXS_0{t0,t1,t2,t3}, TXS_4{t4,t5},TXS_0{t6},TXS_-1{t7,t8,t9},TXS_4{t10,t11}
			//so chunkMap will contain
			// -1 -> [<t7,t9>]
			//  0 -> [<t0,t3>,<t6,t6>]           
			//  4 -> [<t4,t5>,<t10,t11>]
			//
			//if the map has no-texture coords at all in order to unify the code we fill the ChunkMap with texture seam -1 and a single triangle range going from face_0 to face_n-1 


			if (attributestobeupdated[WEDGETEXTUREBO] || attributestobeupdated[VERTTEXTUREBO])
			{
				_chunkmap.clear();
				if (attributestobeupdated[WEDGETEXTUREBO])
					fillChunckMap(_chunkmap);
				else
					if(attributestobeupdated[VERTTEXTUREBO])
						_chunkmap[0].push_back(std::make_pair(0,tn-1));
			}

			//default case: no texture is required to be rendered but a non texture attribute has to be updated
			//we have to init the _chunkmap with just one entry (-1...that means no texture) referring all the triangles in the mesh
			if ((currtextmod == GLFeedEnum::TX_NONE) &&
				(attributestobeupdated[VERTPOSITIONBO] ||
				attributestobeupdated[VERTNORMALBO] || attributestobeupdated[FACENORMALBO] ||
				attributestobeupdated[VERTCOLORBO] || attributestobeupdated[FACECOLORBO]))
			{
				_chunkmap.clear();
				_chunkmap[-1].push_back(std::make_pair(0,tn-1));
			}

			int t = 0;
			if (attributestobeupdated[WEDGETEXTUREBO] || attributestobeupdated[VERTTEXTUREBO])
			{
				_texindnumtriangles.clear();
				_texindnumtriangles.resize(_chunkmap.size());
			}
			int i = 0;
			size_t chunkindex = i;
			GLuint triangles = 0;
			for(ChunkMap::const_iterator mit = _chunkmap.begin();mit != _chunkmap.end();++mit)
			{
				for (ChunkVector::const_iterator cit = mit->second.begin();cit != mit->second.end();++cit)
				{
					for(size_t indf = cit->first;indf<=cit->second;++indf)
					{
						chunkindex = i % facechunk;
						if (attributestobeupdated[VERTPOSITIONBO])
						{
							rpv[chunkindex*3+0].Import(_mesh.face[indf].V(0)->P());
							rpv[chunkindex*3+1].Import(_mesh.face[indf].V(1)->P());
							rpv[chunkindex*3+2].Import(_mesh.face[indf].V(2)->P());
						}
						if (attributestobeupdated[VERTNORMALBO])
						{
							rnv[chunkindex*3+0].Import(_mesh.face[indf].V(0)->N().Normalize());
							rnv[chunkindex*3+1].Import(_mesh.face[indf].V(1)->N().Normalize());
							rnv[chunkindex*3+2].Import(_mesh.face[indf].V(2)->N().Normalize());
						}
						else if (attributestobeupdated[FACENORMALBO])
						{
							rnv[chunkindex*3+0].Import(_mesh.face[indf].N().Normalize());
							rnv[chunkindex*3+1].Import(_mesh.face[indf].N().Normalize());
							rnv[chunkindex*3+2].Import(_mesh.face[indf].N().Normalize());
						}

						if ((attributestobeupdated[VERTCOLORBO]))
						{
							rcv[chunkindex*3+0] = _mesh.face[indf].V(0)->C();
							rcv[chunkindex*3+1] = _mesh.face[indf].V(1)->C();
							rcv[chunkindex*3+2] = _mesh.face[indf].V(2)->C();
						}
						else if ((attributestobeupdated[FACECOLORBO]))
						{
							rcv[chunkindex*3+0] = _mesh.face[indf].C();
							rcv[chunkindex*3+1] = _mesh.face[indf].C();
							rcv[chunkindex*3+2] = _mesh.face[indf].C();
						}

						if (attributestobeupdated[WEDGETEXTUREBO])
						{
							rtv[chunkindex*6+0]=float(_mesh.face[indf].WT(0).U());
							rtv[chunkindex*6+1]=float(_mesh.face[indf].WT(0).V());
							rtv[chunkindex*6+2]=float(_mesh.face[indf].WT(1).U());
							rtv[chunkindex*6+3]=float(_mesh.face[indf].WT(1).V());
							rtv[chunkindex*6+4]=float(_mesh.face[indf].WT(2).U());
							rtv[chunkindex*6+5]=float(_mesh.face[indf].WT(2).V());
						}
						else if (attributestobeupdated[VERTTEXTUREBO])
						{
							rtv[chunkindex*6+0]=float(_mesh.face[indf].V(0)->T().U());
							rtv[chunkindex*6+1]=float(_mesh.face[indf].V(0)->T().V());
							rtv[chunkindex*6+2]=float(_mesh.face[indf].V(1)->T().U());
							rtv[chunkindex*6+3]=float(_mesh.face[indf].V(1)->T().V());
							rtv[chunkindex*6+4]=float(_mesh.face[indf].V(2)->T().U());
							rtv[chunkindex*6+5]=float(_mesh.face[indf].V(2)->T().V());
						}

						if((i == tn - 1) || (chunkindex == facechunk - 1))
						{
							size_t chunksize = facechunk;
							if (i == tn - 1)
								chunksize = chunkindex + 1;

							if (attributestobeupdated[VERTPOSITIONBO])
							{
								GLBufferObject* buffobj = _bo[VERTPOSITIONBO];
								glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
								glBufferSubData(GL_ARRAY_BUFFER,chunckingpu * facechunk * 3 * buffobj->_components * buffobj->getSizeOfGLType(),3 * buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&rpv[0]);
								glBindBuffer(GL_ARRAY_BUFFER, 0);
							}
							if (attributestobeupdated[VERTNORMALBO] || attributestobeupdated[FACENORMALBO])
							{
								GLBufferObject* buffobj;
								if (attributestobeupdated[VERTNORMALBO])
									buffobj = _bo[VERTNORMALBO];
								else
									buffobj = _bo[FACENORMALBO];
								glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
								glBufferSubData(GL_ARRAY_BUFFER,chunckingpu * facechunk * 3 * buffobj->_components * buffobj->getSizeOfGLType(),3 * buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&rnv[0]);
								glBindBuffer(GL_ARRAY_BUFFER, 0);
							}
							if (attributestobeupdated[VERTCOLORBO] || attributestobeupdated[FACECOLORBO])
							{
								GLBufferObject* buffobj;
								if (attributestobeupdated[VERTCOLORBO])
									buffobj = _bo[VERTCOLORBO];
								else
									buffobj = _bo[FACECOLORBO];
								glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
								glBufferSubData(GL_ARRAY_BUFFER,chunckingpu * facechunk * 3 *buffobj->_components * buffobj->getSizeOfGLType(),3 * buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&rcv[0]);
								glBindBuffer(GL_ARRAY_BUFFER, 0);
							}
							if (attributestobeupdated[VERTTEXTUREBO] || attributestobeupdated[WEDGETEXTUREBO])
							{
								GLBufferObject* buffobj;
								if (attributestobeupdated[VERTTEXTUREBO])
									buffobj = _bo[VERTTEXTUREBO];
								else
									buffobj = _bo[WEDGETEXTUREBO];
								glBindBuffer(GL_ARRAY_BUFFER, buffobj->_bohandle);
								glBufferSubData(GL_ARRAY_BUFFER,chunckingpu * facechunk * 3 *buffobj->_components * buffobj->getSizeOfGLType(),3 * buffobj->_components * buffobj->getSizeOfGLType() * chunksize,&rtv[0]);
								glBindBuffer(GL_ARRAY_BUFFER, 0);
							}

							++chunckingpu;
						}
						++i;
					}
					triangles += cit->second - cit->first + 1;
				}

				if (attributestobeupdated[WEDGETEXTUREBO] || attributestobeupdated[VERTTEXTUREBO])
					_texindnumtriangles[t] = std::make_pair(mit->first,triangles);
				++t;
			}

			//return (k != tn)
			//    throw MeshLabException("Mesh has not been properly partitioned");
			return true;
		}

		bool immediateMode(NORMAL_MODALITY nm,COLOR_MODALITY cm,TEXTURE_MODALITY tm,const std::vector<GLuint>& textureindex)
		{
			glPushAttrib(GL_ALL_ATTRIB_BITS);
			if(_mesh.fn==0) 
				return false;

			if(cm == CL_PERMESH)
				glColor(_mesh.C());

			//typename MESHTYPE::FaceContainer::iterator fp;
			typename MESHTYPE::FaceIterator fi = _mesh.face.begin();

			short curtexname=-1;
			if(tm == TX_PERWEDGE)
			{
				curtexname=(*fi).WT(0).n();
				if ((curtexname >= 0) && (curtexname < (int)textureindex.size()))
				{
					glEnable(GL_TEXTURE_2D);
					glBindTexture(GL_TEXTURE_2D,textureindex[curtexname]);
				}
				else
				{
					glDisable(GL_TEXTURE_2D);
				}
			}

			if(tm==TX_PERVERT && !textureindex.empty()) // in the case of per vertex tex coord we assume that we have a SINGLE texture.
			{
				curtexname = 0;
				glEnable(GL_TEXTURE_2D);
				glBindTexture(GL_TEXTURE_2D,textureindex[curtexname]);
			}

			glBegin(GL_TRIANGLES);

			while(fi!=_mesh.face.end())
			{
				typename MESHTYPE::FaceType & f = *fi;
				if(!f.IsD())
				{
					if(tm==TX_PERWEDGE)
						if(f.WT(0).n() != curtexname)
						{
							curtexname=(*fi).WT(0).n();
							glEnd();

							if (curtexname >= 0)
							{
								glEnable(GL_TEXTURE_2D);
								if(!textureindex.empty())
									glBindTexture(GL_TEXTURE_2D,textureindex[curtexname]);
							}
							else
							{
								glDisable(GL_TEXTURE_2D);
							}

							glBegin(GL_TRIANGLES);
						}

						if(nm == NR_PERFACE)	glNormal(f.cN());
						if(nm == NR_PERVERT)	glNormal(f.V(0)->cN());
						//if(nm == NMPerWedge)glNormal(f.WN(0));

						if(cm == CL_PERFACE)	glColor(f.C());
						if(cm == CL_PERVERT)	glColor(f.V(0)->C());
						if(tm == TX_PERVERT) glTexCoord(f.V(0)->T().P());
						if(tm == TX_PERWEDGE) glTexCoord(f.WT(0).t(0));
						glVertex(f.V(0)->P());

						if(nm == NR_PERVERT)	glNormal(f.V(1)->cN());
						//if(nm == NMPerWedge)glNormal(f.WN(1));
						if(cm == CL_PERVERT)	glColor(f.V(1)->C());
						if(tm==TX_PERVERT) glTexCoord(f.V(1)->T().P());
						if(tm==TX_PERWEDGE) glTexCoord(f.WT(1).t(0));
						glVertex(f.V(1)->P());

						if(nm == NR_PERVERT)	glNormal(f.V(2)->cN());
						//if(nm == NMPerWedge)glNormal(f.WN(2));
						if(cm == CL_PERVERT) glColor(f.V(2)->C());
						if(tm==TX_PERVERT) glTexCoord(f.V(2)->T().P());
						if(tm==TX_PERWEDGE) glTexCoord(f.WT(2).t(0));
						glVertex(f.V(2)->P());
				}
				++fi;
			}

			glEnd();
			glPopMatrix();
			glPopAttrib();
			return true;
		}

		bool drawTriangles(NORMAL_MODALITY nm,COLOR_MODALITY cm,TEXTURE_MODALITY tm,const std::vector<GLuint>& textureindex)
		{
			std::vector<bool> att(_bo.size(),false);
			attributesToBeImportedInTriangleBasedPipeline(att,nm,cm,tm);
			bool replicated = !att[VERTINDEXBO];
			glBindVertexArray(_vaohandle);

			if (replicated)
			{
				int firsttriangleoffset = 0;
				if(tm == GLFeedEnum::TX_NONE)
				{
					glDisable(GL_TEXTURE_2D);
					glDrawArrays(GL_TRIANGLES,0,_mesh.fn * 3);
				}
				else
				{
					glEnable(GL_TEXTURE_2D);
					for(std::vector< std::pair<short,GLuint> >::const_iterator it = _texindnumtriangles.begin();it != _texindnumtriangles.end();++it)
					{
						if ((it->first != -1) && (it->first < textureindex.size()))
							glBindTexture(GL_TEXTURE_2D,textureindex[it->first]);
						else
							glBindTexture(GL_TEXTURE_2D,0);
						glDrawArrays(GL_TRIANGLES,firsttriangleoffset,it->second * 3 - firsttriangleoffset);
						firsttriangleoffset = it->second * 3;
					}
					glBindTexture(GL_TEXTURE_2D,0);
					glDisable(GL_TEXTURE_2D);
				}

			}
			else
			{
				if(tm==TX_PERVERT)
				{
					if (textureindex.size() > 0)
					{
						glEnable(GL_TEXTURE_2D);
						glBindTexture(GL_TEXTURE_2D,textureindex[0]);

					}
				}
				else
					glDisable(GL_TEXTURE_2D);


				if  (_bo[VERTINDEXBO]->_isvalid)
				{
					glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,_bo[VERTINDEXBO]->_bohandle);
					glDrawElements( GL_TRIANGLES, _mesh.fn * _bo[VERTINDEXBO]->_components,GL_UNSIGNED_INT ,NULL);
					glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
				}

				glBindTexture(GL_TEXTURE_2D,0);
				glDisable(GL_TEXTURE_2D);
			}
			glBindVertexArray(0);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
			glBindBuffer(GL_ARRAY_BUFFER,0);

			return true;
		}

		void drawPoints(NORMAL_MODALITY nm,COLOR_MODALITY cm)
		{
			glDisable(GL_TEXTURE_2D);

			glBindVertexArray(_vaohandle);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _bo[VERTINDEXBO]->_bohandle);
			glDrawArrays(GL_POINTS,0,_mesh.vn);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);	
			glBindVertexArray(0);
		}

		void setBufferPointerEnableClientState( BO_NAMES boname)
		{
			if ((boname < VERTPOSITIONBO) || (boname > VERTINDEXBO))
				return; 
			GLBufferObject* cbo = _bo[boname];
			if (cbo == NULL)
				return;

			switch(boname)
			{
			case(VERTPOSITIONBO):
				{ 
					glVertexPointer(cbo->_components, cbo->_gltype, 0, 0);
					glEnableClientState(GL_VERTEX_ARRAY);
					break;
				}
			case(VERTNORMALBO):
			case(FACENORMALBO):
				{  
					glNormalPointer(cbo->_gltype, 0, 0);
					glEnableClientState(GL_NORMAL_ARRAY);
					break;
				}
			case(VERTCOLORBO):
			case(FACECOLORBO):
				{
					glColorPointer(cbo->_components, cbo->_gltype, 0, 0);
					glEnableClientState(GL_COLOR_ARRAY);
					break;
				}
			case(VERTTEXTUREBO):
			case(WEDGETEXTUREBO):
				{
					glTexCoordPointer(cbo->_components, cbo->_gltype, 0, 0);  
					glEnableClientState(GL_TEXTURE_COORD_ARRAY);
					break;
				}
			case(VERTINDEXBO):
				{
					glIndexPointer(cbo->_gltype, 0, 0);
					glEnableClientState(GL_INDEX_ARRAY);
					break;
				}
			}
		}

		void disableClientState( BO_NAMES boname,const std::vector<bool>& importatt)
		{
			if ((boname < VERTPOSITIONBO) || (boname > VERTINDEXBO))
				return; 

			switch(boname)
			{
			case(VERTPOSITIONBO):
				{ 
					glDisableClientState(GL_VERTEX_ARRAY);
					break;
				}
			case(VERTNORMALBO):
			case(FACENORMALBO):
				{  
					if (!importatt[VERTNORMALBO] && !importatt[FACENORMALBO])
						glDisableClientState(GL_NORMAL_ARRAY);
					break;
				}
			case(VERTCOLORBO):
			case(FACECOLORBO):
				{
					if (!importatt[VERTCOLORBO] && !importatt[FACECOLORBO])
						glDisableClientState(GL_COLOR_ARRAY);
					break;
				}
			case(VERTTEXTUREBO):
			case(WEDGETEXTUREBO):
				{  
					if (!importatt[VERTTEXTUREBO] && !importatt[WEDGETEXTUREBO])
						glDisableClientState(GL_TEXTURE_COORD_ARRAY);
					break;
				}
			case(VERTINDEXBO):
				{
					glDisableClientState(GL_INDEX_ARRAY);
					break;
				}
			}
		}

		//expected number of cells should have a the named bo
		//generateindex is true when i have a triangle based mesh
		//				is false when i have a point based mesh
		size_t boExpectedSize(BO_NAMES name,bool replicatedpipeline,bool generateindex) const
		{
			try
			{
				GLBufferObject& cbo = *_bo.at(name);
				size_t vertnum(_mesh.VN());
				size_t facenum(_mesh.FN());

				switch(name)
				{
				case(VERTPOSITIONBO):
				case(VERTNORMALBO):
				case(VERTCOLORBO):
				case(VERTTEXTUREBO):
					{
						if (replicatedpipeline)
							return facenum * 3 * cbo._components;
						else
							return vertnum * cbo._components;
					}

				case(FACENORMALBO):
				case(FACECOLORBO):
				case(WEDGETEXTUREBO):
					{
						if (replicatedpipeline)
							return facenum * 3 * cbo._components;
						else 
							return 0;
					}
				case(VERTINDEXBO):
					{
						if (replicatedpipeline || !generateindex)
							return 0;
						else
							return facenum * cbo._components;
					}
				}
			}
			catch(std::out_of_range& /*exc*/)
			{
				return 0;
			}
			return 0;
		}

		//generateindex is true when i have a triangle based mesh
		//				is false when i have a point based mesh		
		size_t boExpectedDimension(BO_NAMES name,bool replicatedpipeline,bool generateindex) const 
		{
			try
			{
				size_t sz = boExpectedSize(name,replicatedpipeline,generateindex);
				GLBufferObject& cbo = *_bo.at(name);
				return sz * cbo.getSizeOfGLType();
			}
			catch(std::out_of_range& /*exc*/)
			{
				return 0;
			}
			return 0;
		}

		static bool isPerVertexAttribute(BO_NAMES name)
		{
			return ((name == VERTPOSITIONBO) ||(name == VERTNORMALBO) || (name == VERTCOLORBO) || (name == VERTTEXTUREBO));
		}

		//it's a map containing for each texture seams n a vector of all the triangle index ranges having n has texture seam
		//Suppose that in a mesh we have
		//TXS_0{t0,t1,t2,t3}, TXS_4{t4,t5},TXS_0{t6},TXS_-1{t7,t8,t9},TXS_4{t10,t11}
		//so chunkMap will contain
		// -1 -> [<t7,t9>]
		//  0 -> [<t0,t3>,<t6,t6>]           
		//  4 -> [<t4,t5>,<t10,t11>]

		typedef std::vector< std::pair< GLuint,GLuint > > ChunkVector;
		typedef std::map< short, ChunkVector > ChunkMap;

		void fillChunckMap(ChunkMap& cmap)
		{
			if (!vcg::tri::HasPerWedgeTexCoord(_mesh))
				return;
			cmap.clear();
			typename MESHTYPE::FaceIterator infrange = _mesh.face.begin();
			short texind = std::numeric_limits<short>::max();
			int hh = 0;
			for(typename MESHTYPE::FaceIterator fit = _mesh.face.begin();fit != _mesh.face.end();++fit)
			{
				if (fit->WT(0).N() != texind)
				{
					if ((texind != std::numeric_limits<short>::max()) || (fit == _mesh.face.end() - 1))
					{
						GLuint lowind = std::distance(_mesh.face.begin(),infrange);
						GLuint topind = std::distance(_mesh.face.begin(),fit) - 1;
						cmap[texind].push_back(std::make_pair(lowind,topind));
						infrange = fit;
					}
					texind = fit->WT(0).N();
				}
				++hh;
			}
			cmap[texind].push_back(std::make_pair(std::distance(_mesh.face.begin(),infrange),std::distance(_mesh.face.begin(),_mesh.face.end() - 1)));
		}

		//ideally this should be const. I'm not yet sure if VCGLib will allow me to declare it as constant
		MESHTYPE& _mesh;
		MemoryInfo& _gpumeminfo;

		std::vector<GLBufferObject*> _bo;
		GLuint _vaohandle;
		std::vector< std::pair<short,GLuint> > _texindnumtriangles;

		bool _lastfeedingusedreplicatedpipeline;
		bool _borendering;
		size_t _perbatchsimplex;
		ChunkMap  _chunkmap;
	};

	template<typename MESHTYPE>
	class SceneToBeRendered
	{
	public:
		SceneToBeRendered(MemoryInfo& gpumeminfo,bool highprecision,size_t perbatchtriangles = 100)
			:_scene(),_gpumeminfo(gpumeminfo),_perbatchtriangles(perbatchtriangles),_globalscenecenter(),_highprecision(highprecision)
		{
		}

		~SceneToBeRendered()
		{
			_scene.clear();
		}

		void insert(MESHTYPE* mesh,const vcg::Matrix44<typename MESHTYPE::ScalarType>& transfmat)
		{
			MatrixedFeeder matfeed = std::make_pair(transfmat,GLMeshAttributesFeeder<MESHTYPE>(*mesh,_gpumeminfo,_perbatchtriangles));
			_scene[mesh] = matfeed;
			if (_highprecision)
				computeSceneGlobalCenter();
		}

		void remove(MESHTYPE* mesh)
		{
            typename std::map< MESHTYPE*, MatrixedFeeder >::iterator it = _scene.find(mesh);
			if (it != _scene.end())
				_scene.erase(it);
		}

		void update(MESHTYPE* mesh,int mask)
		{
            typename std::map< MESHTYPE*, MatrixedFeeder >::iterator it = _scene.find(mesh);
			if (it != _scene.end())
				it->second._feeder.update(mask);
		}

		void passTrianglesToOpenGL(MESHTYPE* mesh,GLFeedEnum::NORMAL_MODALITY nm,GLFeedEnum::COLOR_MODALITY cm,GLFeedEnum::TEXTURE_MODALITY tm,const std::vector<GLuint>& textureindex = std::vector<GLuint>())
		{
            typename std::map<MESHTYPE*,MatrixedFeeder>::iterator it = _scene.find(mesh);
			if(it == _scene.end())
				return;
			glPushAttrib(GL_TRANSFORM_BIT);
			glMatrixMode(GL_MODELVIEW);
			glPushMatrix();
			if (_highprecision)
			{
				vcg::glMultMatrix(it->second._localmeshmatrix);
				vcg::glTranslate(_globalscenecenter);
			}
			it->second._feeder.passTrianglesToOpenGL(nm,cm,tm,textureindex);
			glPopMatrix();
            glPopAttrib();
		}

		void passPointsToOpenGL(MESHTYPE* mesh,GLFeedEnum::NORMAL_MODALITY nm,GLFeedEnum::COLOR_MODALITY cm)
		{
            typename std::map<MESHTYPE*,MatrixedFeeder>::iterator it = _scene.find(mesh);
			if(it == _scene.end())
				return;
			glPushAttrib(GL_TRANSFORM_BIT);
			glMatrixMode(GL_MODELVIEW);
			glPushMatrix();
			if (_highprecision)
			{
				vcg::glMultMatrix(it->second._localmeshmatrix);
				vcg::glTranslate(_globalscenecenter);
			}
			it->second._feeder.passPointsToOpenGL(nm,cm);
			glPopMatrix();
            glPopAttrib();
		}

	private:
		struct MatrixedFeeder
		{
			vcg::Matrix44<typename MESHTYPE::ScalarType> _localmeshmatrix;
			GLMeshAttributesFeeder<MESHTYPE> _feeder;

            MatrixedFeeder(const vcg::Matrix44<typename MESHTYPE::ScalarType>& localmeshmatrix,GLMeshAttributesFeeder<MESHTYPE>& feeder)
				:_localmeshmatrix(localmeshmatrix),_feeder(feeder)
			{
			}
		};

		void computeSceneGlobalCenter()
		{
			vcg::Box3<typename MESHTYPE::ScalarType> scenebbox;
            for(typename std::map<MESHTYPE*,MatrixedFeeder >::const_iterator it = _scene.begin();it != _scene.end();++it)
				scenebbox.Add(it->first->bbox,it->second._localmeshmatrix);

			_globalscenecenter = -scenebbox.Center();
		}

		std::map< MESHTYPE*, MatrixedFeeder > _scene;
		MemoryInfo& _gpumeminfo;
		size_t _perbatchtriangles;
		vcg::Point3<typename MESHTYPE::ScalarType> _globalscenecenter;
		bool _highprecision;
	}; 
}

#endif
