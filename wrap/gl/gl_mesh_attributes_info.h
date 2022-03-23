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
	enum Attribute {
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
	RenderingAtts(bool defaultvalue = false) { reset(defaultvalue); }

	RenderingAtts(const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& att)
	{
		reset();
		//_atts = new bool[ATT_NAMES_DERIVED_CLASS::enumArity()];
		for (unsigned int ii = 0; ii < ATT_NAMES_DERIVED_CLASS::enumArity(); ++ii)
			(*this)[ii] = att[ii];
	}

	~RenderingAtts() {}

	RenderingAtts<ATT_NAMES_DERIVED_CLASS>&
	operator=(const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& att)
	{
		reset();
		for (unsigned int ii = 0; ii < ATT_NAMES_DERIVED_CLASS::enumArity(); ++ii)
			(*this)[ii] = att[ii];
		return (*this);
	}

	bool operator[](unsigned int ind) const
	{
		if (ind >= ATT_NAMES_DERIVED_CLASS::enumArity())
			throw Exception("Out of range value\n");
		return _atts[ind];
	}

	bool& operator[](unsigned int ind)
	{
		if (ind >= ATT_NAMES_DERIVED_CLASS::enumArity())
			throw Exception("Out of range value\n");
		return _atts[ind];
	}

	void reset(bool defaultvalue = false)
	{
		for (unsigned int ii = 0; ii < ATT_NAMES_DERIVED_CLASS::enumArity(); ++ii)
			_atts[ii] = defaultvalue;
	}

	static RenderingAtts<ATT_NAMES_DERIVED_CLASS> unionSet(
		const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& a,
		const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& b)
	{
		RenderingAtts<ATT_NAMES_DERIVED_CLASS> res;
		for (unsigned int ii = 0; ii < ATT_NAMES_DERIVED_CLASS::enumArity(); ++ii)
			res[ii] = a[ii] || b[ii];
		return res;
	}

	static RenderingAtts<ATT_NAMES_DERIVED_CLASS> complementSet(
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

	static RenderingAtts<ATT_NAMES_DERIVED_CLASS> intersectionSet(
		const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& a,
		const RenderingAtts<ATT_NAMES_DERIVED_CLASS>& b)
	{
		RenderingAtts<ATT_NAMES_DERIVED_CLASS> res;
		for (unsigned int ii = 0; ii < ATT_NAMES_DERIVED_CLASS::enumArity(); ++ii)
			res[ii] = a[ii] && b[ii];
		return res;
	}

	size_t serialize(std::string& str) const
	{
		for (unsigned int ii = 0; ii < ATT_NAMES_DERIVED_CLASS::enumArity(); ++ii)
			str.append(((_atts[ii]) ? "1" : "0"));
		return ATT_NAMES_DERIVED_CLASS::enumArity();
	}

	void deserialize(const std::string& str)
	{
		std::bitset<ATT_NAMES_DERIVED_CLASS::ATT_ARITY> bset(str);
		for (unsigned int ii = 0; ii < ATT_NAMES_DERIVED_CLASS::enumArity(); ++ii)
			_atts[ATT_NAMES_DERIVED_CLASS::enumArity() - ii - 1] = bset[ii];
	}

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

	DebugInfo() :
			_tobeallocated(),
			_tobedeallocated(),
			_tobeupdated(),
			_currentlyallocated(),
			_perviewdata()
	{
	}

	void reset()
	{
		_tobeallocated.clear();
		_tobedeallocated.clear();
		_tobeupdated.clear();
		_currentlyallocated.clear();
		_perviewdata.clear();
	}

	static const char* primitiveName(size_t ind)
	{
		static std::string res;

		if (ind == size_t(PR_POINTS))
			res = std::string("PR_POINTS");

		if (ind == size_t(PR_WIREFRAME_EDGES))
			res = std::string("PR_WIREFRAME_EDGES");

		if (ind == size_t(PR_WIREFRAME_TRIANGLES))
			res = std::string("PR_WIREFRAME_TRIANGLES");

		if (ind == size_t(PR_SOLID))
			res = std::string("PR_SOLID");

		return res.c_str();
	}
};

class GLMeshAttributesInfo::INT_ATT_NAMES : public ATT_NAMES
{
public:
	/* WARNING!!!!!! the edges index bo it's just used only by the edges and quads meshes, NOT
	 * by the triangle meshes. Triangles meshes use just the vertex index array WHY? cause quads
	 * meshes need both index arrays. One to render the "usual" mesh, the other one to render
	 * the wireframe quadrangulation on top of it A triangles meshes rendered in wireframe or in
	 * solid wireframe, in order to save GPU memory, use the glPolygonMode approach Edges meshes
	 * uses the edges index array just for a matter of coherence
	 */

	/*WARNING!!!! to be changed whit ATT_NAMES::enumArity() (and so on...) as soon as constexpr
	 * will be supported by most of the old c++ compilers*/
	enum { ATT_VERTINDICES = 7, ATT_EDGEINDICES = 8, ATT_ARITY = 9 };

	INT_ATT_NAMES() : ATT_NAMES() {}

	INT_ATT_NAMES(unsigned int att) : ATT_NAMES()
	{
		if (att >= INT_ATT_NAMES::enumArity())
			throw Exception("Out of range value\n");
		else
			_val = att;
	}

	static unsigned int enumArity() { return INT_ATT_NAMES::ATT_ARITY; }

	operator unsigned int() const { return _val; }
};

class GLMeshAttributesInfo::InternalRendAtts : public GLMeshAttributesInfo::RenderingAtts<INT_ATT_NAMES>
{
public:
	typedef INT_ATT_NAMES AttName;

	InternalRendAtts() : RenderingAtts<INT_ATT_NAMES>() {}

	InternalRendAtts(const RendAtts& reqatt) : RenderingAtts<INT_ATT_NAMES>()
	{
		for (unsigned int ii = 0; ii < ATT_NAMES::enumArity(); ++ii) {
			(*this)[ii] = reqatt[ii];
		}

		(*this)[INT_ATT_NAMES::ATT_VERTINDICES] = false;
		(*this)[INT_ATT_NAMES::ATT_EDGEINDICES] = false;
	}

	InternalRendAtts(const RendAtts& reqatt, PRIMITIVE_MODALITY pm) :
			RenderingAtts<INT_ATT_NAMES>()
	{
		for (unsigned int ii = 0; ii < ATT_NAMES::enumArity(); ++ii)
			(*this)[ii] = reqatt[ii];

		(*this)[INT_ATT_NAMES::ATT_VERTINDICES] = isVertexIndexingRequired(reqatt, pm);
		(*this)[INT_ATT_NAMES::ATT_EDGEINDICES] = isEdgeIndexingRequired(pm);
	}

	InternalRendAtts(const RenderingAtts<INT_ATT_NAMES>& r) : RenderingAtts<INT_ATT_NAMES>(r) {}

	   // upcast from InternalRendAtts to RendAtts
	operator RendAtts() const
	{
		RendAtts rendatt;
		for (unsigned int ii = 0; ii < ATT_NAMES::enumArity(); ++ii)
			rendatt[ii] = _atts[ii];
		return rendatt;
	}

	InternalRendAtts& setIndexingIfNeeded(PRIMITIVE_MODALITY pm)
	{
		(*this)[INT_ATT_NAMES::ATT_VERTINDICES] = isVertexIndexingRequired((*this), pm);
		(*this)[INT_ATT_NAMES::ATT_EDGEINDICES] = isEdgeIndexingRequired(pm);
		return (*this);
	}

	static bool isPerVertexAttribute(INT_ATT_NAMES name)
	{
		return (
			(name.value() == INT_ATT_NAMES::ATT_VERTPOSITION) ||
			(name.value() == INT_ATT_NAMES::ATT_VERTNORMAL) ||
			(name.value() == INT_ATT_NAMES::ATT_VERTCOLOR) ||
			(name.value() == INT_ATT_NAMES::ATT_VERTTEXTURE));
	}

	static bool replicatedPipelineNeeded(const RendAtts& rqatt)
	{
		return (
			rqatt[INT_ATT_NAMES::ATT_FACENORMAL] || rqatt[INT_ATT_NAMES::ATT_FACECOLOR] ||
			rqatt[INT_ATT_NAMES::ATT_WEDGETEXTURE]);
	}

	static bool isVertexIndexingRequired(const RendAtts& rqatt, PRIMITIVE_MODALITY pm)
	{
		return (
			!replicatedPipelineNeeded(rqatt) &&
			((pm == PR_SOLID) || (pm == PR_WIREFRAME_TRIANGLES)));
	}

	static bool isEdgeIndexingRequired(PRIMITIVE_MODALITY pm)
	{
		return (pm == PR_WIREFRAME_EDGES);
	}
};

} // namespace vcg

#include "gl_mesh_attributes_info.ipp"

#endif
