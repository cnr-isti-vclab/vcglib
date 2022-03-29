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

} // namespace vcg
