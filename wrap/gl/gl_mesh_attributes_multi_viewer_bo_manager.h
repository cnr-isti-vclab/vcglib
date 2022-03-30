/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2022                                           \/)\/    *
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

#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <climits>
#include <string>

#include <wrap/gl/space.h>
#include <wrap/gl/math.h>
#include <vcg/space/color4.h>
#include <vcg/math/matrix44.h>
#include <wrap/system/memory_info.h>
#include <wrap/gl/gl_mesh_attributes_info.h>


namespace vcg
{

struct RenderingModalityGLOptions
{
	bool _perbbox_enabled = false;

	bool _perbbox_fixed_color_enabled = true;
	bool _perpoint_fixed_color_enabled = false;
	bool _perwire_fixed_color_enabled = true;
	bool _persolid_fixed_color_enabled = true;

	Color4b _perbbox_fixed_color = vcg::Color4b(Color4b::White);
	Color4b _perpoint_fixed_color = vcg::Color4b(Color4b::White);
	Color4b _perwire_fixed_color = Color4b(Color4b::DarkGray);
	Color4b _persolid_fixed_color = vcg::Color4b(Color4b::White);

	bool _perbbox_mesh_color_enabled = false;
	bool _perpoint_mesh_color_enabled = false;
	bool _perwire_mesh_color_enabled = false;
	bool _persolid_mesh_color_enabled = false;

	bool _perpoint_noshading = false;
	bool _perwire_noshading = true;
	bool _persolid_noshading = false;

	bool _perpoint_dot_enabled = false;

	float _perpoint_pointsize = 3.f;
	bool _perpoint_pointsmooth_enabled = false;
	bool _perpoint_pointattenuation_enabled = true;

	float _perwire_wirewidth = 1.f;
};

template<typename GL_OPTIONS_DERIVED_TYPE = RenderingModalityGLOptions>
class PerViewData : public GLMeshAttributesInfo
{
public:

	typedef GL_OPTIONS_DERIVED_TYPE GLOptionsType;

	PerViewData();

	virtual ~PerViewData() {};

	bool set(PRIMITIVE_MODALITY pm, const RendAtts& atts);

	bool set(PRIMITIVE_MODALITY pm, ATT_NAMES att, bool onoff);

	bool set(PRIMITIVE_MODALITY pm, bool onoff);

	void set(const GL_OPTIONS_DERIVED_TYPE& opts);

	bool isPrimitiveActive(PRIMITIVE_MODALITY pm) const;

	PRIMITIVE_MODALITY_MASK getPrimitiveModalityMask() const;

	bool get(PRIMITIVE_MODALITY pm, RendAtts& atts) const;

	bool get(GL_OPTIONS_DERIVED_TYPE& opts) const;

	void reset(bool deleteglopts = true);


	void serialize(std::string& str) const;

	bool deserialize(const std::string& str);

protected:
	template<typename MESH_TYPE, typename UNIQUE_VIEW_ID_TYPE, typename XX_GL_OPTIONS_DERIVED_TYPE> friend class NotThreadSafeGLMeshAttributesMultiViewerBOManager;

	typedef std::vector<InternalRendAtts> PerRendModData;

	PRIMITIVE_MODALITY_MASK _pmmask;
	PerRendModData _intatts;

	GL_OPTIONS_DERIVED_TYPE _glopts;
};

/*******************************WARNING!!!!!!!!!!!!!!!!!******************************************
 * You must inherit from NotThreadSafeGLMeshAttributesMultiViewerBOManager, providing thread safe
 * mechanisms, in order to use the bo facilities exposed by the class. In
 * wrap/qt/qt_thread_safe_memory_rendering.h you will find a ready to use class based on
 * QtConcurrency module.
 *************************************************************************************************/
template<typename MESH_TYPE, typename UNIQUE_VIEW_ID_TYPE = unsigned int, typename GL_OPTIONS_DERIVED_TYPE = RenderingModalityGLOptions>
class NotThreadSafeGLMeshAttributesMultiViewerBOManager : public GLMeshAttributesInfo
{
public:
	typedef PerViewData<GL_OPTIONS_DERIVED_TYPE> PVData;

protected:
	NotThreadSafeGLMeshAttributesMultiViewerBOManager(const MESH_TYPE& mesh, MemoryInfo& meminfo, size_t perbatchprimitives);
	~NotThreadSafeGLMeshAttributesMultiViewerBOManager();


	void meshAttributesUpdated(bool hasmeshconnectivitychanged, const RendAtts& changedrendatts);

	bool getPerViewInfo(UNIQUE_VIEW_ID_TYPE viewid, PVData& data) const;

	void setPerViewInfo(UNIQUE_VIEW_ID_TYPE viewid, const PVData& data);

	void setPerAllViewsInfo(const PVData& data);

	bool removeView(UNIQUE_VIEW_ID_TYPE viewid);

	void removeAllViews();

	void draw(UNIQUE_VIEW_ID_TYPE viewid, const std::vector<GLuint>& textid = std::vector<GLuint>()) const;

	void drawAllocatedAttributesSubset(UNIQUE_VIEW_ID_TYPE viewid, const PVData& dt, const std::vector<GLuint>& textid = std::vector<GLuint>()) const;

	bool isBORenderingAvailable() const;

	bool manageBuffers();

	void setGLOptions(UNIQUE_VIEW_ID_TYPE viewid, const GL_OPTIONS_DERIVED_TYPE& opts);

	void setTrMatrix(const vcg::Matrix44<typename MESH_TYPE::ScalarType>& tr);

	void setDebugMode(bool isdebug);

	void getLog(DebugInfo& info);

private:

	class EdgeVertInd
	{
	public:

		GLuint  _v[2];  // the two Vertex indices are ordered!

		EdgeVertInd() {}
		EdgeVertInd(const MESH_TYPE& m, const typename MESH_TYPE::FaceType* pf, const int nz) { this->set(m, pf, nz); }
		EdgeVertInd(const MESH_TYPE& m, const typename MESH_TYPE::EdgeType* pe, const int nz) { this->set(m, pe, nz); }

		void set(const MESH_TYPE& m, const typename MESH_TYPE::FaceType* pf, const int nz)
		{
			assert(pf != 0);
			assert(nz >= 0);
			assert(nz < pf->VN());

			_v[0] = GLuint(vcg::tri::Index(m, pf->V(nz)));;
			_v[1] = GLuint(vcg::tri::Index(m, pf->V(pf->Next(nz))));
			assert(_v[0] != _v[1]); // The face pointed by 'f' is Degenerate (two coincident vertexes)

			if (_v[0] > _v[1])
				std::swap(_v[0], _v[1]);
		}

		void set(const MESH_TYPE& m, const typename MESH_TYPE::EdgeType* pe, const int nz)
		{
			assert(pe != 0);
			assert(nz >= 0);
			assert(nz < 2);

			_v[0] = GLuint(vcg::tri::Index(m, pe->V(nz)));;
			_v[1] = GLuint(vcg::tri::Index(m, pe->V((nz + 1) % 2)));
			assert(_v[0] != _v[1]); // The face pointed by 'f' is Degenerate (two coincident vertexes)

			if (_v[0] > _v[1])
				std::swap(_v[0], _v[1]);
		}

		inline bool operator<(const EdgeVertInd& pe) const
		{
			if (_v[0] < pe._v[0])
				return true;
			else if (_v[0] > pe._v[0])
				return false;
			else
				return _v[1] < pe._v[1];
		}

		inline bool operator==(const EdgeVertInd & pe) const
		{
			return _v[0] == pe._v[0] && _v[1] == pe._v[1];
		}
	};

	struct GLBufferObject
	{
		GLBufferObject(size_t components, GLenum gltype, GLenum clientstatetag, GLenum target)
				:_size(0), _components(components), _isvalid(false), _gltype(gltype), _target(target), _clientstatetag(clientstatetag), _bohandle(0)
		{
		}

		GLBufferObject(size_t components, GLenum gltype, GLenum target)
				:_size(0), _components(components), _isvalid(false), _gltype(gltype), _target(target), _clientstatetag(), _bohandle(0)
		{
		}

		size_t getSizeOfGLType() const
		{
			switch (_gltype)
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
		const GLenum _target;

		/*WARNING!!!!!!!!!!!!!!!!! In openGL INDEX BO doesn't require to be enabled/disabled so has NOT a valid tag associated.
		In this case the client state tag remains not initialized and it's not meaningful */
		const GLenum _clientstatetag;
		/**********************************************************************************/

		GLuint _bohandle;
	};

	const MESH_TYPE& _mesh;

	MemoryInfo& _gpumeminfo;

	/*The buffer objects used for the rendering operation. They are shared among all the views*/
	std::vector<GLBufferObject*> _bo;

	typedef std::map< UNIQUE_VIEW_ID_TYPE, PVData > ViewsMap;

		 ///*_perviewreqatts contains a map of the requested atts by each single view. it's maintained for the actual rendering step*/
	ViewsMap _perviewreqatts;

	/*_currboatt contains the union of all the requested attributes by each single view on the scene. At the end it represents the BOs allocated in the GPU memory*/
	/* WARNING!!!! The currently allocated BOs are the union of all the BOs requested to be visualized in the _perviewreqatts plus, possibly, the edgeindex bo (depending by the kind of per view render primitive modality that is requested) and the vertexindex bo (depending of the solid rendering modality, per-face/per-vertex)*/
	/*             The EdgeIndex bo is allocated only if one of the requested rendering modality is PR_WIREFRAME_EDGES or PR_WIREFRAME_EDGES. If we have PR_SOLID the glPolygonMode function is used for rendering the triangle wireframe view*/
	InternalRendAtts _currallocatedboatt;

	bool _borendering;
	size_t _perbatchprim;

	/*Additional structures used for per wedge texturing modality*/
	typedef std::vector< std::pair< GLuint, GLuint > > ChunkVector;
	typedef std::map< short, ChunkVector > ChunkMap;

	std::vector< std::pair<short, GLuint> > _texindnumtriangles;
	ChunkMap  _chunkmap;

		 //Horrible waste of memory space...but computing the list of edges is too much expensive...we must minimize it!
	std::vector<EdgeVertInd> _edge;
	size_t _meshverticeswhenedgeindiceswerecomputed;
	size_t _meshtriangleswhenedgeindiceswerecomputed;

		 //vcg::GLOptions _glopts;
	vcg::Matrix44<typename MESH_TYPE::ScalarType> _tr;

	bool _debugmode;
	DebugInfo _loginfo;

	std::vector<InternalRendAtts> _meaningfulattsperprimitive;

	void initMeaningfulAttsMask();

	bool hasMeshAttribute(INT_ATT_NAMES attname) const;

	bool checkBuffersAllocationStatus(InternalRendAtts& tobeallocated, InternalRendAtts& tobedeallocated, InternalRendAtts& tobeupdated) const;

	bool manageAndFeedBuffersIfNeeded(const InternalRendAtts& tobeallocated, const InternalRendAtts& tobedeallocated, const InternalRendAtts& tobeupdated);

	void updateEdgeVertIndVector();

	bool buffersMemoryManagementFunction(const InternalRendAtts& tobeallocated, const InternalRendAtts& tobedeallocated, const InternalRendAtts& tobeupdated);

	bool updateBuffersIndexedPipeline(const InternalRendAtts& attributestobeupdated);

	bool isThereAReplicatedPipelineView() const;

	bool isThereAnEdgesView() const;

	bool updateBuffersReplicatedPipeline(const InternalRendAtts& attributestobeupdated);

	void buffersDeAllocationRequested(const InternalRendAtts& rq);

	void bufferDeAllocationRequested(INT_ATT_NAMES att);

	std::ptrdiff_t bufferObjectsMemoryRequired(const InternalRendAtts& rqatt) const;

	size_t boExpectedSize(INT_ATT_NAMES name, bool replicatedpipeline) const;

	size_t boExpectedDimension(INT_ATT_NAMES name, bool replicatedpipeline) const;

	void drawFun(const PVData& dt, const std::vector<GLuint>& textid = std::vector<GLuint>()) const;

	void drawFilledTriangles(const InternalRendAtts& req, const GL_OPTIONS_DERIVED_TYPE* glopts, const std::vector<GLuint>& textureindex = std::vector<GLuint>()) const;

	void drawWiredTriangles(const InternalRendAtts& req, const GL_OPTIONS_DERIVED_TYPE* glopts, const std::vector<GLuint>& textureindex = std::vector<GLuint>()) const;

	void drawTrianglesBO(const InternalRendAtts& req, const std::vector<GLuint>& textureindex = std::vector<GLuint>()) const;

	void drawTrianglesIM(const InternalRendAtts& req, const std::vector<GLuint>& textureindex = std::vector<GLuint>()) const;

	void drawPoints(const InternalRendAtts& req, const GL_OPTIONS_DERIVED_TYPE* glopts, const std::vector<GLuint>& textureindex = std::vector<GLuint>()) const;

	void drawPointsBO(const InternalRendAtts& req) const;

	void drawPointsIM(const InternalRendAtts& req) const;

	void drawEdges(const InternalRendAtts& req, const GL_OPTIONS_DERIVED_TYPE* glopts) const;

	void drawEdgesBO(const InternalRendAtts& req) const;

	void drawEdgesIM(const InternalRendAtts& req) const;

	void drawBBox(const GL_OPTIONS_DERIVED_TYPE* glopts) const;

	void drawBBoxBO() const;

	void drawBBoxIM() const;

	void updateClientState(const InternalRendAtts& req) const;

	void setBufferPointer(INT_ATT_NAMES boname) const;

	void disableClientState(INT_ATT_NAMES boname, const RendAtts& req) const;

	void fillchunkMap();

	void debug(const InternalRendAtts& tobeallocated, const InternalRendAtts& tobedeallocated, const InternalRendAtts& tobeupdated);

	static void fillEdgeVector(const MESH_TYPE &m, std::vector<EdgeVertInd> &edgeVec, bool includeFauxEdge = true);

	static void fillUniqueEdgeVector(const MESH_TYPE &m, std::vector<EdgeVertInd> &edgeVec);
};

} // namespace vcg

#include "gl_mesh_attributes_multi_viewer_bo_manager.ipp"

#endif // __VCG_GL_MESH_ATTRIBUTES_FEEDER
