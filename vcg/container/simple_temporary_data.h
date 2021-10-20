/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
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

#ifndef __VCGLIB_SIMPLE__
#define __VCGLIB_SIMPLE__

#include <string>
#include <cstring>
#include <limits>
#include <vector>
#include <cassert>

namespace vcg
{

class SimpleTempDataBase
{
public:
    virtual ~SimpleTempDataBase() {}
    SimpleTempDataBase() {}
    virtual void Resize(size_t sz) = 0;
    virtual void Reorder(std::vector<size_t> &newVertIndex) = 0;
    virtual size_t SizeOf() const = 0;
    virtual void *DataBegin() = 0;
    virtual const void* DataBegin() const = 0;

    virtual void       *At(size_t i) = 0;
    virtual const void *At(size_t i) const = 0;
    virtual void CopyValue(const size_t to, const size_t from, const SimpleTempDataBase *other) = 0;
};

template <class TYPE, class ...p>
class VectorNBW : public std::vector<TYPE, p...>
{
};

template <class ...p>
class VectorNBW<bool, p...>
{
public:
    VectorNBW() : booldata(nullptr), datasize(0), datareserve(0) {}

    ~VectorNBW()
    {
        if (booldata)
            delete[] booldata;
    }

	void reserve(size_t sz)
	{
		if (sz <= datareserve)
			return;
		bool* newdataLoc = new bool[sz];
		if (datasize != 0) {
			std::copy(booldata, booldata+datasize, newdataLoc);
			// memcpy(newdataLoc, booldata, sizeof(bool) * sizeof(datasize));
		}

		std::swap(booldata, newdataLoc);
		if (newdataLoc != 0)
			delete[] newdataLoc;
		datareserve = sz;
	}

	void resize(size_t sz)
	{
		int oldDatasize = datasize;
		if ((int) sz <= oldDatasize)
			return;
		if (sz > datareserve)
			reserve(sz);
		datasize = sz;
		for (unsigned int i = oldDatasize; i < datasize; ++i)
			booldata[i] = false;
	}
	void push_back(const bool &v)
    {
        resize(datasize + 1);
        booldata[datasize] = v;
    }

    void clear() { datasize = 0; }

    unsigned int size() const { return datasize; }

    bool empty() const { return datasize == 0; }

    bool* data() {return booldata;}
    const bool *data() const { return booldata; }

    bool &operator[](size_t i) { return booldata[i]; }
    const bool &operator[](size_t i) const { return booldata[i]; }

private:
    bool *booldata;
    size_t datasize;
    size_t datareserve;
};

template <class STL_CONT, class ATTR_TYPE>
class SimpleTempData : public SimpleTempDataBase
{

public:
    typedef SimpleTempData<STL_CONT, ATTR_TYPE> SimpTempDataType;
    typedef ATTR_TYPE AttrType;

    const STL_CONT &c;
    VectorNBW<ATTR_TYPE> data;
    int padding;

    SimpleTempData(const STL_CONT &_c) : c(_c), padding(0)
    {
        data.reserve(c.capacity());
        data.resize(c.size());
    };
    SimpleTempData(const STL_CONT &_c, const ATTR_TYPE &val) : c(_c)
    {
        data.reserve(c.capacity());
        data.resize(c.size());
        Init(val);
    };

    ~SimpleTempData()
    {
        data.clear();
    }

    void Init(const ATTR_TYPE &val)
    {
        std::fill(data.begin(), data.end(), val);
    }
    // access to data
    ATTR_TYPE &operator[](const typename STL_CONT::value_type &v)  { return data[&v - &*c.begin()]; }
    ATTR_TYPE &operator[](const typename STL_CONT::value_type *v)  { return data[v - &*c.begin()]; }
    ATTR_TYPE &operator[](const typename STL_CONT::const_iterator &cont) { return data[&(*cont) - &*c.begin()]; }
    ATTR_TYPE &operator[](const typename STL_CONT::iterator &cont) { return data[&(*cont) - &*c.begin()]; }
    ATTR_TYPE &operator[](size_t i) { return data[i]; }

    const ATTR_TYPE &operator[](const typename STL_CONT::value_type &v)  const { return data[&v - &*c.begin()]; }
    const ATTR_TYPE &operator[](const typename STL_CONT::value_type *v)  const { return data[v - &*c.begin()]; }
    const ATTR_TYPE &operator[](const typename STL_CONT::const_iterator &cont) const { return data[&(*cont) - &*c.begin()]; }
    const ATTR_TYPE &operator[](const typename STL_CONT::iterator &cont) const { return data[&(*cont) - &*c.begin()]; }
    const ATTR_TYPE &operator[](size_t i) const { return data[i]; }

    void       *At(size_t i) { return &(*this)[i]; }
    const void *At(size_t i) const { return &(*this)[i]; }

    void CopyValue(const size_t to, const size_t from, const SimpleTempDataBase *other)
    {
        assert(other != nullptr);
        data[to] = *(static_cast<const ATTR_TYPE *>(other->At(from)));
    }

    // update temporary data size
    bool UpdateSize()
    {
        if (data.size() != c.size())
        {
            data.resize(c.size());
            return false;
        }
        return true;
    }

    void Resize(size_t sz)
    {
        data.resize(sz);
    }

    void Reorder(std::vector<size_t> &newVertIndex)
    {
        for (size_t i = 0; i < data.size(); ++i)
        {
            if (newVertIndex[i] != (std::numeric_limits<size_t>::max)())
                data[newVertIndex[i]] = data[i];
        }
    }

    size_t SizeOf() const { return sizeof(ATTR_TYPE); }
    void *DataBegin() { return data.empty() ? nullptr : data.data(); }
    const void *DataBegin() const { return data.empty() ? nullptr : data.data(); }
};

template <class ATTR_TYPE>
class Attribute : public SimpleTempDataBase
{
public:
    typedef ATTR_TYPE AttrType;
    Attribute() { attribute = new ATTR_TYPE(); }
    ~Attribute() { delete attribute; }
    size_t SizeOf() const { return sizeof(ATTR_TYPE); }
    void *DataBegin() { return attribute; }
    const void* DataBegin() const {return attribute;}

    void Resize(size_t) { assert(0); }
    void Reorder(std::vector<size_t> &) { assert(0); }

    void *At(size_t)
    {
        assert(0);
        return (void *)0;
    }
    const void *At(size_t) const
    {
        assert(0);
        return (void *)0;
    }
    void CopyValue(const size_t, const size_t, const SimpleTempDataBase *) { assert(0); }
private:
    AttrType *attribute;
};

} // end namespace vcg

#endif
