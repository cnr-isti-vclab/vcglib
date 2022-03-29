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

#ifndef __VCG_GL_MESH_ATTRIBUTES_INFO
#define __VCG_GL_MESH_ATTRIBUTES_INFO

#include <bitset>
#include <string>
#include <vector>

#include <assert.h>

namespace vcg {

struct GLMeshAttributesInfo
{
	enum PRIMITIVE_MODALITY {
		PR_POINTS              = 0,
		PR_WIREFRAME_EDGES     = 1,
		PR_WIREFRAME_TRIANGLES = 2,
		PR_SOLID               = 3,
		PR_ARITY               = 4
	};

	class Exception;

	class ATT_NAMES;

	template<typename ATT_NAMES_DERIVED_CLASS>
	class RenderingAtts;

	typedef std::bitset<PR_ARITY> PRIMITIVE_MODALITY_MASK;
	typedef RenderingAtts<ATT_NAMES> RendAtts;

	class DebugInfo;

	static PRIMITIVE_MODALITY next(PRIMITIVE_MODALITY pm);

protected:
	class INT_ATT_NAMES;

	class InternalRendAtts;
};

class GLMeshAttributesInfo::Exception : public std::exception
{
public:
	Exception(const char* text);

	~Exception() {}
	const char* what() const noexcept;

private:
	std::string _text;
};

class GLMeshAttributesInfo::ATT_NAMES
{
public:
	enum {
		ATT_VERTPOSITION = 0,
		ATT_VERTNORMAL   = 1,
		ATT_FACENORMAL   = 2,
		ATT_VERTCOLOR    = 3,
		ATT_FACECOLOR    = 4,
		ATT_VERTTEXTURE  = 5,
		ATT_WEDGETEXTURE = 6,
		ATT_ARITY        = 7
	};

	ATT_NAMES();

	ATT_NAMES(unsigned int att);

	static constexpr unsigned int enumArity();

	operator unsigned int() const;

	unsigned int value() const;

protected:
	unsigned int _val;
};

template<typename ATT_NAMES_DERIVED_CLASS>
class GLMeshAttributesInfo::RenderingAtts
{
public:
	RenderingAtts(bool defaultvalue = false);

	bool operator[](unsigned int ind) const;

	bool& operator[](unsigned int ind);

	void reset(bool defaultvalue = false);

	static RenderingAtts<ATT_NAMES_DERIVED_CLASS> unionSet(
		const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& a,
		const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& b);

	static RenderingAtts<ATT_NAMES_DERIVED_CLASS> complementSet(
		const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& a,
		const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& b);

	static RenderingAtts<ATT_NAMES_DERIVED_CLASS> intersectionSet(
		const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& a,
		const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& b);

	std::string serialize() const;

	void deserialize(const std::string& str);

protected:
	bool _atts[ATT_NAMES_DERIVED_CLASS::ATT_ARITY]; // an array of enumArity() bool values
};

class GLMeshAttributesInfo::DebugInfo
{
public:
	std::string _tobeallocated;
	std::string _tobedeallocated;
	std::string _tobeupdated;

	std::string _currentlyallocated;

	std::vector<std::string> _perviewdata;

	DebugInfo();

	void reset();

	static const char* primitiveName(size_t ind);
};

class GLMeshAttributesInfo::INT_ATT_NAMES : public ATT_NAMES
{
public:
	/* WARNING!!!!!! the edges index bo it's just used only by the edges and quads meshes, NOT
	 * by the triangle meshes. Triangles meshes use just the vertex index array
	 * WHY? because quads meshes need both index arrays. One to render the "usual" mesh, the other
	 * one to render the wireframe quadrangulation on top of it A triangles meshes rendered in
	 * wireframe or in solid wireframe, in order to save GPU memory, use the glPolygonMode approach
	 * Edges meshes uses the edges index array just for a matter of coherence
	 */

	enum {
		ATT_VERTINDICES = 7,
		ATT_EDGEINDICES = 8,
		ATT_ARITY = 9
	};

	INT_ATT_NAMES();

	INT_ATT_NAMES(unsigned int att);

	static constexpr unsigned int enumArity();

	operator unsigned int() const;

	unsigned int value() const;
};

class GLMeshAttributesInfo::InternalRendAtts : public GLMeshAttributesInfo::RenderingAtts<INT_ATT_NAMES>
{
public:
	typedef INT_ATT_NAMES AttName;

	InternalRendAtts();

	InternalRendAtts(const RendAtts& reqatt);

	InternalRendAtts(const RendAtts& reqatt, PRIMITIVE_MODALITY pm);

	InternalRendAtts(const RenderingAtts<INT_ATT_NAMES>& r);

	// upcast from InternalRendAtts to RendAtts
	operator RendAtts() const;

	void setIndexingIfNeeded(PRIMITIVE_MODALITY pm);

	static bool isPerVertexAttribute(INT_ATT_NAMES name);

	static bool replicatedPipelineNeeded(const RendAtts& rqatt);

	static bool isVertexIndexingRequired(const RendAtts& rqatt, PRIMITIVE_MODALITY pm);

	static bool isEdgeIndexingRequired(PRIMITIVE_MODALITY pm);
};

} // namespace vcg

#include "gl_mesh_attributes_info.ipp"

#endif // __VCG_GL_MESH_ATTRIBUTES_INFO
