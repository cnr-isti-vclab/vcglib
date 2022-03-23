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

} // namespace vcg
