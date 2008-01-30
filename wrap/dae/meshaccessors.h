#ifndef _MESH_ACCESSORS_H
#define _MESH_ACCESSORS_H


template<typename POINTTYPE>
struct CoordNumber
{
public:
	static unsigned int coord()
	{
		return 0;
	}
};

template<>
struct CoordNumber<vcg::Point2f>
{
public:
	static unsigned int coord() 
	{
		return 2;
	}
};

template<>
struct CoordNumber<vcg::Point3f>
{
public:
	static unsigned int coord() 
	{
		return 3;
	}
};

template<>
struct CoordNumber<vcg::Point4f>
{
public:
	static unsigned int coord() 
	{
		return 4;
	}
};




namespace MeshAccessors
{
template<typename MESHTYPE>
class VertexAccessor
{
public:
	typedef typename MESHTYPE::ConstVertexIterator ConstIterator;
	typedef typename MESHTYPE::VertexIterator Iterator;
	
	VertexAccessor(MESHTYPE& m)
		:_m(m)
	{

	}

	ConstIterator begin() const
	{
		return _m.vert.begin();
	}

	ConstIterator end() const
	{
		return _m.vert.end();
	}

	Iterator begin()
	{
		return _m.vert.begin();
	}

	Iterator end()
	{
		return _m.vert.end();
	}

	MESHTYPE& _m;
};
	
template<typename MESHTYPE>
class VertexPositionAccessor : public VertexAccessor<MESHTYPE>
{
public:
	VertexPositionAccessor(MESHTYPE& m)
		:VertexAccessor<MESHTYPE>(m)
	{
	}

	template<typename VERTEXTYPE>
	typename const VERTEXTYPE::CoordType operator()(VERTEXTYPE& v) const
	{
		return v.P();
	}
};

template<typename MESHTYPE>
class VertexNormalAccessor : public VertexAccessor<MESHTYPE>
{
public:
	VertexNormalAccessor(MESHTYPE& m)
		:VertexAccessor<MESHTYPE>(m)
	{
	}

	template<typename VERTEXTYPE>
	typename const VERTEXTYPE::NormalType operator()(VERTEXTYPE& v) const 
	{
		return v.N();
	}
};

template<typename MESHTYPE>
class VertexTextureCoordinateAccessor : public VertexAccessor<MESHTYPE>
{
public:
	VertexTextureCoordinateAccessor(MESHTYPE& m)
		:VertexAccessor<MESHTYPE>(m)
	{
	}

	template<typename VERTEXTYPE>
	typename const VERTEXTYPE::TexCoordType operator()(VERTEXTYPE& v) const 
	{
		return v.T();
	}
};

template<typename MESHTYPE>
class PolygonAccessor
{
public:
	typedef typename MESHTYPE::ConstFaceIterator ConstIterator;
	typedef typename MESHTYPE::FaceIterator Iterator;

	PolygonAccessor(MESHTYPE& m)
		:_m(m)
	{

	}

	ConstIterator begin() const
	{
		return _m.face.begin();
	}

	ConstIterator end() const
	{
		return _m.face.end();
	}

	Iterator begin()
	{
		return _m.face.begin();
	}

	Iterator end()
	{
		return _m.face.end();
	}

	MESHTYPE& _m;
};

template<typename MESHTYPE>
class PolygonPositionAccessor : public PolygonAccessor<MESHTYPE>
{
public:
	PolygonPositionAccessor<MESHTYPE>(MESHTYPE& m)
		:PolygonAccessor(m)
	{
	}

	template<typename POLYGONTYPE>
	typename const POLYGONTYPE::CoordType operator()(POLYGONTYPE& p,const unsigned int index) const 
	{
		return p.V(i)->P();
	}
};

template<typename MESHTYPE>
class PolygonNormalAccessor : public PolygonAccessor<MESHTYPE>
{
public:
	PolygonNormalAccessor(MESHTYPE& m)
		:PolygonAccessor<MESHTYPE>(m)
	{
	}

	template<typename POLYGONTYPE>
	typename const POLYGONTYPE::NormalType operator()(const POLYGONTYPE& p) const
	{
		POLYGONTYPE::NormalType r = p.cN();
		return r.Normalize();
	}
};

template<typename MESHTYPE>
class PolygonWedgeTextureCoordinateAccessor : public PolygonAccessor<MESHTYPE>
{
public:
	template<typename POLYGONTYPE>
	typename const POLYGONTYPE::TexCoordType operator()(POLYGONTYPE& p,const unsigned int index) const
	{
		return p.cWT(index);
	}

	PolygonWedgeTextureCoordinateAccessor(MESHTYPE& m)
		:PolygonAccessor<MESHTYPE>(m)
	{
	}
};
} //MeshAccessors

template<typename MESHTYPE,typename ACCESSOR>
struct AccessorInfo
{
	AccessorInfo(const MESHTYPE& m)
		:_a(m)
	{

	}

	ACCESSOR _a; 
};

template<typename MESHTYPE,typename ACCESSOR>
struct AccessorComponentNumberInfo  : public AccessorInfo<MESHTYPE,ACCESSOR>
{
	AccessorComponentNumberInfo(const MESHTYPE& m)
		:AccessorInfo(m)
	{
	}

	unsigned int _return_value_component_number;
};

template<typename MESHTYPE,typename ACCESSOR>
struct AccessorPolygonEdgeNumberInfo  : public AccessorInfo<MESHTYPE,ACCESSOR>
{
	AccessorPolygonEdgeNumberInfo(const MESHTYPE& m)
		:AccessorInfo(m)	
	{
	}

	unsigned int _edge_number;
};

#endif