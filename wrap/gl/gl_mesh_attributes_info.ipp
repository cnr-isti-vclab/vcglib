/*****************************************************************************
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

#include "gl_mesh_attributes_info.h"

namespace vcg {

/* GLMeshAttributesInfo Class */

inline GLMeshAttributesInfo::PRIMITIVE_MODALITY GLMeshAttributesInfo::next(PRIMITIVE_MODALITY pm)
{
	int tmp = static_cast<int>(pm);
	if (tmp == PR_ARITY)
		throw Exception("PRIMITIVE_MODALITY iterator: PR_ARITY passed as parameter!");
	++tmp;
	return static_cast<PRIMITIVE_MODALITY>(tmp);
}

/* GLMeshAttributesInfo::Exception */

inline GLMeshAttributesInfo::Exception::Exception(const char* text) :
		std::exception(), _text(text)
{
}

inline const char* GLMeshAttributesInfo::Exception::what() const noexcept
{
	return _text.c_str();
}

/* GLMeshAttributesInfo::ATT_NAMES */

inline GLMeshAttributesInfo::ATT_NAMES::ATT_NAMES() : _val(ATT_VERTPOSITION)
{
}

inline GLMeshAttributesInfo::ATT_NAMES::ATT_NAMES(unsigned int att)
{
	if (att >= enumArity())
		throw Exception("Out of range value\n");
	else
		_val = att;
}

inline constexpr unsigned int GLMeshAttributesInfo::ATT_NAMES::enumArity()
{
	return ATT_NAMES::ATT_ARITY;
}

inline GLMeshAttributesInfo::ATT_NAMES::operator unsigned int() const
{
	return _val;
}

inline unsigned int GLMeshAttributesInfo::ATT_NAMES::value() const
{
	return _val;
}

template<typename ATT_NAMES_DERIVED_CLASS>
GLMeshAttributesInfo::RenderingAtts<ATT_NAMES_DERIVED_CLASS>::RenderingAtts(bool defaultvalue)
{
	reset(defaultvalue);
}

template<typename ATT_NAMES_DERIVED_CLASS>
bool &GLMeshAttributesInfo::RenderingAtts<ATT_NAMES_DERIVED_CLASS>::operator[](unsigned int ind)
{
	if (ind >= ATT_NAMES_DERIVED_CLASS::enumArity())
		throw Exception("Out of range value\n");
	return _atts[ind];
}

template<typename ATT_NAMES_DERIVED_CLASS>
bool GLMeshAttributesInfo::RenderingAtts<ATT_NAMES_DERIVED_CLASS>::operator[](unsigned int ind) const
{
	if (ind >= ATT_NAMES_DERIVED_CLASS::enumArity())
		throw Exception("Out of range value\n");
	return _atts[ind];
}

template<typename ATT_NAMES_DERIVED_CLASS>
void GLMeshAttributesInfo::RenderingAtts<ATT_NAMES_DERIVED_CLASS>::reset(bool defaultvalue)
{
	for (unsigned int ii = 0; ii < ATT_NAMES_DERIVED_CLASS::enumArity(); ++ii)
		_atts[ii] = defaultvalue;
}

template<typename ATT_NAMES_DERIVED_CLASS>
GLMeshAttributesInfo::RenderingAtts<ATT_NAMES_DERIVED_CLASS>
GLMeshAttributesInfo::RenderingAtts<ATT_NAMES_DERIVED_CLASS>::unionSet(
	const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& a,
	const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& b)
{
	RenderingAtts<ATT_NAMES_DERIVED_CLASS> res;
	for (unsigned int ii = 0; ii < ATT_NAMES_DERIVED_CLASS::enumArity(); ++ii)
		res[ii] = a[ii] || b[ii];
	return res;
}

template<typename ATT_NAMES_DERIVED_CLASS>
GLMeshAttributesInfo::RenderingAtts<ATT_NAMES_DERIVED_CLASS>
GLMeshAttributesInfo::RenderingAtts<ATT_NAMES_DERIVED_CLASS>::complementSet(
	const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& a,
	const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& b)
{
	/* TRUTH TABLE */
	// this[ATT_NAMES] | rq[ATT_NAMES] | res
	//     true        |     true      | false
	//     true        |     false     | true
	//     false       |     true      | false
	//     false       |     false     | false

	RenderingAtts<ATT_NAMES_DERIVED_CLASS> res = a;
	for (unsigned int ii = 0; ii < ATT_NAMES_DERIVED_CLASS::enumArity(); ++ii) {
		if (res[ii])
			res[ii] = !(b[ii]);
	}

	return res;
}

template<typename ATT_NAMES_DERIVED_CLASS>
GLMeshAttributesInfo::RenderingAtts<ATT_NAMES_DERIVED_CLASS>
GLMeshAttributesInfo::RenderingAtts<ATT_NAMES_DERIVED_CLASS>::intersectionSet(
	const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& a,
	const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& b)
{
	RenderingAtts<ATT_NAMES_DERIVED_CLASS> res;
	for (unsigned int ii = 0; ii < ATT_NAMES_DERIVED_CLASS::enumArity(); ++ii)
		res[ii] = a[ii] && b[ii];
	return res;
}

template<typename ATT_NAMES_DERIVED_CLASS>
std::string GLMeshAttributesInfo::RenderingAtts<ATT_NAMES_DERIVED_CLASS>::serialize() const
{
	std::string str;
	for (unsigned int ii = 0; ii < ATT_NAMES_DERIVED_CLASS::enumArity(); ++ii)
		str.append(((_atts[ii]) ? "1" : "0"));
	return str;
}

template<typename ATT_NAMES_DERIVED_CLASS>
void GLMeshAttributesInfo::RenderingAtts<ATT_NAMES_DERIVED_CLASS>::deserialize(
	const std::string& str)
{
	std::bitset<ATT_NAMES_DERIVED_CLASS::ATT_ARITY> bset(str);
	for (unsigned int ii = 0; ii < ATT_NAMES_DERIVED_CLASS::enumArity(); ++ii)
		_atts[ATT_NAMES_DERIVED_CLASS::enumArity() - ii - 1] = bset[ii];
}

/* GLMeshAttributesInfo::RenderingAtts */

} // namespace vcg
