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

#include "gl_mesh_attributes_multi_viewer_bo_manager.h"

namespace vcg {

/* PerViewData */

template<typename GL_OPTIONS_DERIVED_TYPE>
PerViewData<GL_OPTIONS_DERIVED_TYPE>::PerViewData() :
		_pmmask(), _intatts(PR_ARITY)
{
	reset();
}

template<typename GL_OPTIONS_DERIVED_TYPE>
bool PerViewData<GL_OPTIONS_DERIVED_TYPE>::set(PRIMITIVE_MODALITY pm, const RendAtts& atts)
{
	size_t pmind(pm);
	if (pm >= _intatts.size())
		return false;
	//_pmmask.set(pm);
	_intatts[pmind] = InternalRendAtts(atts, pm);
	_pmmask.set(size_t(pm), _intatts[pmind][INT_ATT_NAMES::ATT_VERTPOSITION]);
	return true;
}

template<typename GL_OPTIONS_DERIVED_TYPE>
bool PerViewData<GL_OPTIONS_DERIVED_TYPE>::set(PRIMITIVE_MODALITY pm, ATT_NAMES att, bool onoff)
{
	size_t pmind(pm);
	if (pm >= _intatts.size())
		return false;
	_intatts[pmind][att] = onoff;
	_pmmask.set(size_t(pm), _intatts[pmind][INT_ATT_NAMES::ATT_VERTPOSITION]);
	if (_pmmask.test(size_t(pm)))
		_intatts[pmind].setIndexingIfNeeded(pm);
	return true;
}

template<typename GL_OPTIONS_DERIVED_TYPE>
bool PerViewData<GL_OPTIONS_DERIVED_TYPE>::set(PRIMITIVE_MODALITY pm, bool onoff)
{
	return set(pm, INT_ATT_NAMES::ATT_VERTPOSITION, onoff);
}

template<typename GL_OPTIONS_DERIVED_TYPE>
void PerViewData<GL_OPTIONS_DERIVED_TYPE>::set(const GL_OPTIONS_DERIVED_TYPE& opts)
{
	_glopts = opts;
}

template<typename GL_OPTIONS_DERIVED_TYPE>
bool PerViewData<GL_OPTIONS_DERIVED_TYPE>::isPrimitiveActive(PRIMITIVE_MODALITY pm) const
{
	if (pm == PR_ARITY)
		return false;
	return (_pmmask.test(pm) && _intatts[size_t(pm)][INT_ATT_NAMES::ATT_VERTPOSITION]);
}

template<typename GL_OPTIONS_DERIVED_TYPE>
GLMeshAttributesInfo::PRIMITIVE_MODALITY_MASK PerViewData<GL_OPTIONS_DERIVED_TYPE>::getPrimitiveModalityMask() const
{
	return _pmmask;
}

template<typename GL_OPTIONS_DERIVED_TYPE>
bool PerViewData<GL_OPTIONS_DERIVED_TYPE>::get(PRIMITIVE_MODALITY pm, RendAtts& atts) const
{
	size_t pmind(pm);
	if (pm >= _intatts.size())
		return false;
	atts = _intatts[pmind];
	return true;
}

template<typename GL_OPTIONS_DERIVED_TYPE>
bool PerViewData<GL_OPTIONS_DERIVED_TYPE>::get(GL_OPTIONS_DERIVED_TYPE& opts) const
{
	opts = _glopts;
	return true;
}

template<typename GL_OPTIONS_DERIVED_TYPE>
void PerViewData<GL_OPTIONS_DERIVED_TYPE>::reset(bool deleteglopts)
{
	_pmmask.reset();
	for (typename PerRendModData::iterator it = _intatts.begin(); it != _intatts.end(); ++it)
		it->reset();
	if (deleteglopts)
	{
		_glopts = GL_OPTIONS_DERIVED_TYPE();
	}
}

template<typename GL_OPTIONS_DERIVED_TYPE>
void PerViewData<GL_OPTIONS_DERIVED_TYPE>::serialize(std::string& str) const
{
	str.append(_pmmask.to_string());
	for (typename PerRendModData::const_iterator it = _intatts.begin(); it != _intatts.end(); ++it)
	{
		std::string s = it->serialize();
		str.append(s);
	}
	std::string s;
	_glopts.serialize(s);
	str.append(s);
}

template<typename GL_OPTIONS_DERIVED_TYPE>
bool PerViewData<GL_OPTIONS_DERIVED_TYPE>::deserialize(const std::string& str)
{
	std::string::size_type pos = 0;
	std::string token[6];
	token[0] = str.substr(pos, _pmmask.size());
	if (token[0].length() < _pmmask.size())
		return false;
	int i = 1;
	pos = _pmmask.size();
	for (typename PerRendModData::iterator it = _intatts.begin(); it != _intatts.end(); ++it, i++)
	{
		token[i] = str.substr(pos, InternalRendAtts::AttName::enumArity());
		if (token[i].length() < InternalRendAtts::AttName::enumArity())
			return false;
		pos = pos + InternalRendAtts::AttName::enumArity();
	}

	std::string tmp;
	size_t size = _glopts.serialize(tmp);
	token[i] = str.substr(pos, size);
	if (token[i].length() < size)
		return false;
	_pmmask = PRIMITIVE_MODALITY_MASK(token[0]);
	i = 1;
	for (typename PerRendModData::iterator it = _intatts.begin(); it != _intatts.end(); ++it, i++)
		it->deserialize(token[i]);
	_glopts.deserialize(token[i]);
	return true;
}

/* NotThreadSafeGLMeshAttributesMultiViewerBOManager */

template<typename MESH_TYPE, typename UNIQUE_VIEW_ID_TYPE, typename GL_OPTIONS_DERIVED_TYPE>
NotThreadSafeGLMeshAttributesMultiViewerBOManager<
	MESH_TYPE,
	UNIQUE_VIEW_ID_TYPE,
	GL_OPTIONS_DERIVED_TYPE>::
	NotThreadSafeGLMeshAttributesMultiViewerBOManager(
		MESH_TYPE&  mesh,
		MemoryInfo& meminfo,
		size_t      perbatchprimitives) :
		_mesh(mesh),
		_gpumeminfo(meminfo),
		_bo(INT_ATT_NAMES::enumArity(), NULL),
		_currallocatedboatt(),
		_borendering(false),
		_perbatchprim(perbatchprimitives),
		_chunkmap(),
		_edge(),
		_meshverticeswhenedgeindiceswerecomputed(0),
		_meshtriangleswhenedgeindiceswerecomputed(0),
		_tr(),
		_debugmode(false),
		_loginfo(),
		_meaningfulattsperprimitive(PR_ARITY, InternalRendAtts())
{
	_tr.SetIdentity();
	_bo[INT_ATT_NAMES::ATT_VERTPOSITION] = new GLBufferObject(3, GL_FLOAT, GL_VERTEX_ARRAY, GL_ARRAY_BUFFER);
	_bo[INT_ATT_NAMES::ATT_VERTNORMAL] = new GLBufferObject(3, GL_FLOAT, GL_NORMAL_ARRAY, GL_ARRAY_BUFFER);
	_bo[INT_ATT_NAMES::ATT_FACENORMAL] = new GLBufferObject(3, GL_FLOAT, GL_NORMAL_ARRAY, GL_ARRAY_BUFFER);
	_bo[INT_ATT_NAMES::ATT_VERTCOLOR] = new GLBufferObject(4, GL_UNSIGNED_BYTE, GL_COLOR_ARRAY, GL_ARRAY_BUFFER);
	_bo[INT_ATT_NAMES::ATT_FACECOLOR] = new GLBufferObject(4, GL_UNSIGNED_BYTE, GL_COLOR_ARRAY, GL_ARRAY_BUFFER);
	/*MESHCOLOR has not a buffer object associated with it. It's just a call to glColor3f. it's anyway added to the _bo arrays for sake of coherence*/
	//_bo[INT_ATT_NAMES::ATT_FIXEDCOLOR] = NULL;
	_bo[INT_ATT_NAMES::ATT_VERTTEXTURE] = new GLBufferObject(2, GL_FLOAT, GL_TEXTURE_COORD_ARRAY, GL_ARRAY_BUFFER);
	_bo[INT_ATT_NAMES::ATT_WEDGETEXTURE] = new GLBufferObject(2, GL_FLOAT, GL_TEXTURE_COORD_ARRAY, GL_ARRAY_BUFFER);
	_bo[INT_ATT_NAMES::ATT_VERTINDICES] = new GLBufferObject(3, GL_UNSIGNED_INT, GL_ELEMENT_ARRAY_BUFFER);
	_bo[INT_ATT_NAMES::ATT_EDGEINDICES] = new GLBufferObject(2, GL_UNSIGNED_INT, GL_ELEMENT_ARRAY_BUFFER);

	initMeaningfulAttsMask();
}

template<typename MESH_TYPE, typename UNIQUE_VIEW_ID_TYPE, typename GL_OPTIONS_DERIVED_TYPE>
NotThreadSafeGLMeshAttributesMultiViewerBOManager<
	MESH_TYPE,
	UNIQUE_VIEW_ID_TYPE,
	GL_OPTIONS_DERIVED_TYPE>::~NotThreadSafeGLMeshAttributesMultiViewerBOManager()
{
	for (size_t ii = 0; ii < _bo.size(); ++ii)
		delete _bo[ii];
}

/**
 * meshAttributesUpdate will force the buffer allocation only of the bo rendered at least by one
 * viewer.
 * If a filter add to a mesh, for instance, a per vertex color attribute that was not previously
 * rendered, the meshAttributesUpdate() will ignore the attribute until a viewer require explicitly
 * to render the per-vertex-color, too.
 * In order to do it, please, call the setPerViewRendAtts() setting up for at least one the existing
 * viewer (or adding a new viewer) in the RendAtts reqatts parameter the reqatts[ATT_VERTCOLOR] to
 * true value.
 */
template<typename MESH_TYPE, typename UNIQUE_VIEW_ID_TYPE, typename GL_OPTIONS_DERIVED_TYPE>
void NotThreadSafeGLMeshAttributesMultiViewerBOManager<
	MESH_TYPE,
	UNIQUE_VIEW_ID_TYPE,
	GL_OPTIONS_DERIVED_TYPE>::
	meshAttributesUpdated(bool hasmeshconnectivitychanged, const RendAtts& changedrendatts)
{
	InternalRendAtts tobeupdated(changedrendatts);
	tobeupdated[INT_ATT_NAMES::ATT_VERTINDICES] = hasmeshconnectivitychanged;
	tobeupdated[INT_ATT_NAMES::ATT_EDGEINDICES] = hasmeshconnectivitychanged;
	for (unsigned int ii = 0; ii < INT_ATT_NAMES::enumArity(); ++ii) {
		INT_ATT_NAMES boname(ii);
		if (_bo[boname] != NULL)
			_bo[boname]->_isvalid = (_bo[boname]->_isvalid) && !(tobeupdated[boname]);
	}
}

template<typename MESH_TYPE, typename UNIQUE_VIEW_ID_TYPE, typename GL_OPTIONS_DERIVED_TYPE>
bool NotThreadSafeGLMeshAttributesMultiViewerBOManager<
	MESH_TYPE,
	UNIQUE_VIEW_ID_TYPE,
	GL_OPTIONS_DERIVED_TYPE>::getPerViewInfo(UNIQUE_VIEW_ID_TYPE viewid, PVData& data) const
{
	typename ViewsMap::const_iterator it = _perviewreqatts.find(viewid);
	if (it == _perviewreqatts.end())
		return false;
	data = it->second;
	return true;
}

template<typename MESH_TYPE, typename UNIQUE_VIEW_ID_TYPE, typename GL_OPTIONS_DERIVED_TYPE>
void NotThreadSafeGLMeshAttributesMultiViewerBOManager<
	MESH_TYPE,
	UNIQUE_VIEW_ID_TYPE,
	GL_OPTIONS_DERIVED_TYPE>::setPerViewInfo(UNIQUE_VIEW_ID_TYPE viewid, const PVData& data)
{
	/// cleanup stage...if an attribute impossible for a primitive modality is still here (it should
	/// not be...) we change the required atts into the view
	PVData copydt(data);
	for (PRIMITIVE_MODALITY pm = PRIMITIVE_MODALITY(0); pm < PR_ARITY; pm = next(pm))
		copydt._intatts[pm] = InternalRendAtts::intersectionSet(copydt._intatts[size_t(pm)], _meaningfulattsperprimitive[size_t(pm)]);
	_perviewreqatts[viewid] = copydt;
}

template<typename MESH_TYPE, typename UNIQUE_VIEW_ID_TYPE, typename GL_OPTIONS_DERIVED_TYPE>
void NotThreadSafeGLMeshAttributesMultiViewerBOManager<
	MESH_TYPE,
	UNIQUE_VIEW_ID_TYPE,
	GL_OPTIONS_DERIVED_TYPE>::setPerAllViewsInfo(const PVData& data)
{
	/// cleanup stage...if an attribute impossible for a primitive modality is still here (it should
	/// not be...) we change the required atts into the view
	PVData copydt(data);
	for (PRIMITIVE_MODALITY pm = PRIMITIVE_MODALITY(0); pm < PR_ARITY; pm = next(pm))
		copydt._intatts[pm] = InternalRendAtts::intersectionSet(copydt._intatts[size_t(pm)], _meaningfulattsperprimitive[size_t(pm)]);
	for (typename ViewsMap::iterator it = _perviewreqatts.begin(); it != _perviewreqatts.end(); ++it)
		it->second = copydt;
}

template<typename MESH_TYPE, typename UNIQUE_VIEW_ID_TYPE, typename GL_OPTIONS_DERIVED_TYPE>
bool NotThreadSafeGLMeshAttributesMultiViewerBOManager<
	MESH_TYPE,
	UNIQUE_VIEW_ID_TYPE,
	GL_OPTIONS_DERIVED_TYPE>::removeView(UNIQUE_VIEW_ID_TYPE viewid)
{
	typename ViewsMap::iterator it = _perviewreqatts.find(viewid);
	if (it == _perviewreqatts.end())
		return false;
	_perviewreqatts.erase(viewid);
	return true;
}

template<typename MESH_TYPE, typename UNIQUE_VIEW_ID_TYPE, typename GL_OPTIONS_DERIVED_TYPE>
void NotThreadSafeGLMeshAttributesMultiViewerBOManager<MESH_TYPE, UNIQUE_VIEW_ID_TYPE, GL_OPTIONS_DERIVED_TYPE>::removeAllViews()
{
	_perviewreqatts.clear();
}

template<typename MESH_TYPE, typename UNIQUE_VIEW_ID_TYPE, typename GL_OPTIONS_DERIVED_TYPE>
void NotThreadSafeGLMeshAttributesMultiViewerBOManager<
	MESH_TYPE,
	UNIQUE_VIEW_ID_TYPE,
	GL_OPTIONS_DERIVED_TYPE>::draw(UNIQUE_VIEW_ID_TYPE viewid, const std::vector<GLuint>& textid)
	const
{
	typename ViewsMap::const_iterator it = _perviewreqatts.find(viewid);
	if (it == _perviewreqatts.end())
		return;

	const PVData& dt = it->second;
	//const InternalRendAtts& atts = it->second._intatts;
	drawFun(dt, textid);
}

template<typename MESH_TYPE, typename UNIQUE_VIEW_ID_TYPE, typename GL_OPTIONS_DERIVED_TYPE>
void NotThreadSafeGLMeshAttributesMultiViewerBOManager<
	MESH_TYPE,
	UNIQUE_VIEW_ID_TYPE,
	GL_OPTIONS_DERIVED_TYPE>::
	drawAllocatedAttributesSubset(
		UNIQUE_VIEW_ID_TYPE        viewid,
		const PVData&              dt,
		const std::vector<GLuint>& textid) const
{
	typename ViewsMap::const_iterator it = _perviewreqatts.find(viewid);
	if (it == _perviewreqatts.end())
		return;

	PVData tmp = dt;

	if (!(_currallocatedboatt[INT_ATT_NAMES::ATT_VERTPOSITION])) {
		for (PRIMITIVE_MODALITY pm = PRIMITIVE_MODALITY(0); pm < PR_ARITY; pm = next(pm)) {
			tmp._pmmask[size_t(pm)] = 0;
			tmp._intatts[size_t(pm)] = InternalRendAtts();
		}
	}
	else {
		for (PRIMITIVE_MODALITY pm = PRIMITIVE_MODALITY(0); pm < PR_ARITY; pm = next(pm)) {
			tmp._intatts[size_t(pm)] = InternalRendAtts::intersectionSet(tmp._intatts[size_t(pm)], _meaningfulattsperprimitive[size_t(pm)]);
			tmp._intatts[size_t(pm)] = InternalRendAtts::intersectionSet(tmp._intatts[size_t(pm)], _currallocatedboatt);
		}
	}
	drawFun(dt, textid);
}

} // namespace vcg
