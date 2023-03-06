/***************************************************************************
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
#ifndef __VCGLIB_TRIALLOCATOR
#define __VCGLIB_TRIALLOCATOR

#include <vector>
#include <set>

#include <vcg/container/simple_temporary_data.h>

#include "used_types.h"

namespace vcg {
namespace tri {
/** \addtogroup trimesh
@{
*/

template<class MeshType>
size_t Index(const MeshType &m, const typename MeshType::VertexType &v) {return &v-&*m.vert.begin();}
template<class MeshType>
size_t Index(const MeshType &m, const typename MeshType::FaceType &f) {return &f-&*m.face.begin();}
template<class MeshType>
size_t Index(const MeshType &m, const typename MeshType::EdgeType &e) {return &e-&*m.edge.begin();}
template<class MeshType>
size_t Index(const MeshType &m, const typename MeshType::HEdgeType &h) {return &h-&*m.hedge.begin();}
template <class MeshType>
size_t Index(const MeshType &m, const typename MeshType::TetraType &t) { return &t - &*m.tetra.begin(); }


template<class MeshType>
size_t Index(const MeshType &m, const typename MeshType::VertexType *vp) {return vp-&*m.vert.begin();}
template<class MeshType>
size_t Index(const MeshType &m, const typename MeshType::FaceType * fp) {return fp-&*m.face.begin();}
template<class MeshType>
size_t Index(const MeshType &m, const typename MeshType::EdgeType*  e) {return e-&*m.edge.begin();}
template<class MeshType>
size_t Index(const MeshType &m, const typename MeshType::HEdgeType*  h) {return h-&*m.hedge.begin();}
template <class MeshType>
size_t Index(const MeshType &m, const typename MeshType::TetraType *t) { return t - &*m.tetra.begin(); }


template<class MeshType>
bool IsValidPointer( MeshType & m, const typename MeshType::VertexType *vp) { return ( m.vert.size() > 0 && (vp >= &*m.vert.begin()) && (vp <= &m.vert.back()) ); }
template<class MeshType>
bool IsValidPointer(MeshType & m, const typename MeshType::EdgeType   *ep) { return ( m.edge.size() > 0 && (ep >= &*m.edge.begin()) && (ep <= &m.edge.back())); }
template<class MeshType>
bool IsValidPointer(MeshType & m, const typename MeshType::FaceType   *fp) { return ( m.face.size() > 0 && (fp >= &*m.face.begin()) && (fp <= &m.face.back())); }
template<class MeshType>
bool IsValidPointer(MeshType & m, const typename MeshType::HEdgeType  *hp) { return ( m.hedge.size() > 0 && (hp >= &*m.hedge.begin()) && (hp <= &m.hedge.back())); }
template <class MeshType>
bool IsValidPointer(MeshType &m, const typename MeshType::TetraType *tp) { return (m.tetra.size() > 0 && (tp >= &*m.tetra.begin()) && (tp <= &m.tetra.back())); }

template <class MeshType, class ATTR_CONT>
void ReorderAttribute(ATTR_CONT &c, std::vector<size_t> & newVertIndex, MeshType & /* m */){
  typename std::set<typename MeshType::PointerToAttribute>::iterator ai;
  for(ai = c.begin(); ai != c.end(); ++ai)
    ((typename MeshType::PointerToAttribute)(*ai)).Reorder(newVertIndex);
}

template <class MeshType, class ATTR_CONT>
void ResizeAttribute(ATTR_CONT &c, size_t sz, MeshType &/*m*/){
  typename std::set<typename MeshType::PointerToAttribute>::iterator ai;
  for(ai =c.begin(); ai != c.end(); ++ai)
    ((typename MeshType::PointerToAttribute)(*ai)).Resize(sz);
}

/*!
        \brief  Class to safely add and delete elements in a mesh.

        Adding elements to a mesh, like faces and vertices can involve the reallocation of the vectors of the involved elements.
        This class provide the only safe methods to add elements.
        It also provide an accessory class vcg::tri::PointerUpdater for updating pointers to mesh elements that are kept by the user.
        */
template <class MeshType>
class Allocator
{

public:
  typedef typename MeshType::VertexType     VertexType;
  typedef typename MeshType::VertexPointer  VertexPointer;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::VertContainer VertContainer;

  typedef typename MeshType::EdgeType     EdgeType;
  typedef typename MeshType::EdgePointer  EdgePointer;
  typedef typename MeshType::EdgeIterator EdgeIterator;
  typedef typename MeshType::EdgeContainer EdgeContainer;

  typedef typename MeshType::FaceType       FaceType;
  typedef typename MeshType::FacePointer    FacePointer;
  typedef typename MeshType::FaceIterator   FaceIterator;
  typedef typename MeshType::FaceContainer FaceContainer;

  typedef typename MeshType::HEdgeType     HEdgeType;
  typedef typename MeshType::HEdgePointer  HEdgePointer;
  typedef typename MeshType::HEdgeIterator HEdgeIterator;
  typedef typename MeshType::HEdgeContainer HEdgeContainer;

  typedef typename MeshType::TetraType TetraType;
  typedef typename MeshType::TetraPointer TetraPointer;
  typedef typename MeshType::TetraIterator TetraIterator;
  typedef typename MeshType::TetraContainer TetraContainer;

  typedef typename MeshType::CoordType     CoordType;

  typedef typename MeshType::PointerToAttribute PointerToAttribute;
  typedef typename std::set<PointerToAttribute>::iterator AttrIterator;
  typedef typename std::set<PointerToAttribute>::const_iterator AttrConstIterator;
  typedef typename std::set<PointerToAttribute >::iterator PAIte;

  /*!
            \brief Accessory class to update pointers after eventual reallocation caused by adding elements.

            This class is used whenever you trigger some allocation operation that can cause the invalidation of the pointers to mesh elements.
            Typical situations are when you are allocating new vertexes, edges, halfedges of faces or when you compact
            their containers to get rid of deleted elements.
            This object allows you to update an invalidate pointer immediately after an action that invalidate it.
            \note It can also be used to prevent any update of the various internal pointers caused by an invalidation.
            This can be useful in case you are building all the internal connections by hand as it happens in a importer;
            \sa \ref allocation
            */
  template<class SimplexPointerType>
  class PointerUpdater
  {
  public:
    PointerUpdater(void) : newBase(0), oldBase(0), newEnd(0), oldEnd(0), preventUpdateFlag(false) { ; }
    void Clear(){newBase=oldBase=newEnd=oldEnd=0; remap.clear();}
    /*! \brief Update a pointer to an element of a mesh after a reallocation

             The updating is correctly done only if this PointerUpdater have been passed to the corresponing allocation call. \sa \ref allocation
             */
    void Update(SimplexPointerType &vp)
    {
      //if(vp>=newBase && vp<newEnd) return;
      if(vp<oldBase || vp>oldEnd) return;
      assert(vp>=oldBase);
      assert(vp<oldEnd);
      vp=newBase+(vp-oldBase);
      if(!remap.empty())
        vp  = newBase + remap[vp-newBase];
    }
    /*!
  \brief return true if the allocation operation that initialized this PointerUpdater has caused a reallocation
  */
    bool NeedUpdate() {if((oldBase && newBase!=oldBase && !preventUpdateFlag) || !remap.empty()) return true; else return false;}

    SimplexPointerType newBase;
    SimplexPointerType oldBase;
    SimplexPointerType newEnd;
    SimplexPointerType oldEnd;
    std::vector<size_t> remap; // this vector keep the new position of an element. Uninitialized elements have max_int value to denote an element that has not to be remapped.

    bool preventUpdateFlag; /// when true no update is considered necessary.
  };

  /* +++++++++++++++ Add Vertices ++++++++++++++++ */

  /** \brief Add n vertices to the mesh.
            Function to add n vertices to the mesh.
            The elements are added always to the end of the vector.
            No attempt of reusing previously deleted element is done.
            \sa PointerUpdater
            \param m the mesh to be modified
            \param n the number of elements to be added
            \param pu  a PointerUpdater initialized so that it can be used to update pointers to vertices that could have become invalid after this adding.
            \retval the iterator to the first element added.
            */
  static VertexIterator AddVertices(MeshType &m, size_t n, PointerUpdater<VertexPointer> &pu)
  {
    VertexIterator last;
    if(n == 0) 
      return m.vert.end();
    pu.Clear();
    
    if(m.vert.empty())
      pu.oldBase=0;  // if the vector is empty we cannot find the last valid element
    else {
      pu.oldBase=&*m.vert.begin();
      pu.oldEnd=&m.vert.back()+1;
    }

    m.vert.resize(m.vert.size()+n);
    m.vn+=int(n);

    typename std::set<PointerToAttribute>::iterator ai;
    for(ai = m.vert_attr.begin(); ai != m.vert_attr.end(); ++ai)
      ((PointerToAttribute)(*ai)).Resize(m.vert.size());

    pu.newBase = &*m.vert.begin();
    pu.newEnd =  &m.vert.back()+1;
    if(pu.NeedUpdate())
    {
      for (FaceIterator fi=m.face.begin(); fi!=m.face.end(); ++fi)
        if(!(*fi).IsD())
          for(int i=0; i < (*fi).VN(); ++i)
            if ((*fi).cV(i)!=0) pu.Update((*fi).V(i));

      for (EdgeIterator ei=m.edge.begin(); ei!=m.edge.end(); ++ei)
        if(!(*ei).IsD())
        {
          // if(HasEVAdjacency (m)) 
          pu.Update((*ei).V(0));
          pu.Update((*ei).V(1));
          //							if(HasEVAdjacency(m))   pu.Update((*ei).EVp());
        }

      HEdgeIterator hi;
      for (hi=m.hedge.begin(); hi!=m.hedge.end(); ++hi)
        if(!(*hi).IsD())
        {
          if(HasHVAdjacency (m))
          {
            pu.Update((*hi).HVp());
          }
        }

      for (TetraIterator ti = m.tetra.begin(); ti != m.tetra.end(); ++ti)
        if (!(*ti).IsD())
          for (int i = 0; i < 4; ++i)
            if ((*ti).cV(i) != 0)
              pu.Update((*ti).V(i));

      // e poiche' lo spazio e' cambiato si ricalcola anche last da zero
    }
    size_t siz=(size_t)(m.vert.size()-n);

    last = m.vert.begin();
    advance(last,siz);

    return last;// deve restituire l'iteratore alla prima faccia aggiunta;
  }

  /** \brief Wrapper to AddVertices(); no PointerUpdater
            */
  static VertexIterator AddVertices(MeshType &m, size_t n)
  {
    PointerUpdater<VertexPointer> pu;
    return AddVertices(m, n,pu);
  }

  /** \brief Wrapper to AddVertices() no PointerUpdater but a vector of VertexPointer pointers to be updated
            */
  static VertexIterator AddVertices(MeshType &m, size_t n, std::vector<VertexPointer *> &local_vec)
  {
    PointerUpdater<VertexPointer> pu;
    VertexIterator v_ret =  AddVertices(m, n,pu);

    typename std::vector<VertexPointer *>::iterator vi;
    for(vi=local_vec.begin();vi!=local_vec.end();++vi)
      pu.Update(**vi);
    return v_ret;
  }
  
  /** \brief Wrapper to AddVertices() to add an eigen matrix of (vn,3)
   *  it returns the iterator to the first vertex added
            */
  static VertexIterator AddVertices(MeshType &m, const Eigen::MatrixXf &vm)
  {
      VertexIterator v_ret =  AddVertices(m, vm.rows());
      VertexIterator v_start = v_ret;
      for(int i=0;i<vm.rows();++i)
      {
          v_ret->P()[0]=vm(i,0);
          v_ret->P()[1]=vm(i,1);
          v_ret->P()[2]=vm(i,2);
          ++v_ret;
      }
      return v_start;
  }
  
   
  /** \brief Wrapper to AddVertices() to add a single vertex with given coords
            */
  static VertexIterator AddVertex(MeshType &m, const CoordType &p)
  {
    VertexIterator v_ret =  AddVertices(m, 1);
    v_ret->P()=p;
    return v_ret;
  }

  /** \brief Wrapper to AddVertices() to add a single vertex with given coords and normal
            */
  static VertexIterator AddVertex(MeshType &m, const CoordType &p,  const CoordType &n)
  {
    VertexIterator v_ret =  AddVertices(m, 1);
    v_ret->P()=p;
    v_ret->N()=n;
    return v_ret;
  }

  /** \brief Wrapper to AddVertices() to add a single vertex with given coords and color
            */
  static VertexIterator AddVertex(MeshType &m, const CoordType &p, const Color4b &c)
  {
    VertexIterator v_ret =  AddVertices(m, 1);
    v_ret->P()=p;
    v_ret->C()=c;
    return v_ret;
  }

  /* +++++++++++++++ Add Edges ++++++++++++++++ */

  /** \brief Add n edges to the mesh.
            Function to add n edges to the mesh.
            The elements are added always to the end of the vector. No attempt of reusing previously deleted element is done.
            \sa PointerUpdater
            \param m the mesh to be modified
            \param n the number of elements to be added
            \param pu  a PointerUpdater initialized so that it can be used to update pointers to edges that could have become invalid after this adding.
            \retval the iterator to the first element added.
            */
  static EdgeIterator AddEdges(MeshType &m, size_t n, PointerUpdater<EdgePointer> &pu)
  {
    if(n == 0) return m.edge.end();
    pu.Clear();
    if(m.edge.empty()) pu.oldBase=0;  // if the vector is empty we cannot find the last valid element
    else {
      pu.oldBase=&*m.edge.begin();
      pu.oldEnd=&m.edge.back()+1;
    }

    m.edge.resize(m.edge.size()+n);
    m.en+=int(n);
    size_t siz=(size_t)(m.edge.size()-n);
    EdgeIterator firstNewEdge = m.edge.begin();
    advance(firstNewEdge,siz);

    typename std::set<typename MeshType::PointerToAttribute>::iterator ai;
    for(ai = m.edge_attr.begin(); ai != m.edge_attr.end(); ++ai)
      ((typename MeshType::PointerToAttribute)(*ai)).Resize(m.edge.size());

    pu.newBase = &*m.edge.begin();
    pu.newEnd =  &m.edge.back()+1;
    if(pu.NeedUpdate())
    {
      if(HasFEAdjacency(m))
        for (FaceIterator fi=m.face.begin(); fi!=m.face.end(); ++fi){
          if(!(*fi).IsD())
            for(int i=0; i < (*fi).VN(); ++i)
              if ((*fi).cFEp(i)!=0) pu.Update((*fi).FEp(i));
        }

      if(HasVEAdjacency(m)){
        for (VertexIterator vi=m.vert.begin(); vi!=m.vert.end(); ++vi)
          if(!(*vi).IsD())
            if ((*vi).cVEp()!=0) pu.Update((*vi).VEp());
        for(EdgeIterator ei=m.edge.begin();ei!=firstNewEdge;++ei)
          if(!(*ei).IsD())
          {            
            if ((*ei).cVEp(0)!=0) pu.Update((*ei).VEp(0));
            if ((*ei).cVEp(1)!=0) pu.Update((*ei).VEp(1));
          }        
      }
      
      if(HasHEAdjacency(m))
        for (HEdgeIterator hi=m.hedge.begin(); hi!=m.hedge.end(); ++hi)
          if(!(*hi).IsD())
            if ((*hi).cHEp()!=0) pu.Update((*hi).HEp());
    }
    
    return firstNewEdge;// deve restituire l'iteratore alla prima faccia aggiunta;
  }

  /** Function to add a single edge to the mesh. and initializing it with two VertexPointer
            */
  static EdgeIterator AddEdge(MeshType &m, VertexPointer v0, VertexPointer v1)
  {
    EdgeIterator ei= AddEdges(m, 1);
    ei->V(0)=v0;
    ei->V(1)=v1;
    return ei;
  }
  
  /** Function to add a single edge to the mesh. and initializing it with two indexes to the vertexes
            */
  static EdgeIterator AddEdge(MeshType &m, size_t v0, size_t v1)
  {
    assert(v0!=v1);
    assert(v0>=0 && v0<m.vert.size());
    assert(v1>=0 && v1<m.vert.size());
    return AddEdge(m,&(m.vert[v0]),&(m.vert[v1]));
  }
  

  /** Function to add a face to the mesh and initializing it with the three given coords
            */
  static EdgeIterator AddEdge(MeshType &m, CoordType p0, CoordType p1)
  {
    VertexIterator vi = AddVertices(m,2);
    EdgeIterator ei = AddEdges(m,1);
    vi->P()=p0;
    ei->V(0)=&*vi++;
    vi->P()=p1;
    ei->V(1)=&*vi++;
    return ei;
  }

  /** Function to add n edges to the mesh.
            First wrapper, with no parameters
            */
  static EdgeIterator AddEdges(MeshType &m, size_t n)
  {
    PointerUpdater<EdgePointer> pu;
    return AddEdges(m, n,pu);
  }

  /** Function to add n edges to the mesh.
            Second Wrapper, with a vector of vertex pointers to be updated.
            */
  static EdgeIterator AddEdges(MeshType &m, size_t n, std::vector<EdgePointer*> &local_vec)
  {
    PointerUpdater<EdgePointer> pu;
    EdgeIterator v_ret =  AddEdges(m, n,pu);

    typename std::vector<EdgePointer *>::iterator ei;
    for(ei=local_vec.begin();ei!=local_vec.end();++ei)
      pu.Update(**ei);
    return v_ret;
  }

  /* +++++++++++++++ Add HalfEdges ++++++++++++++++ */

  /** Function to add n halfedges to the mesh. The second parameter hold a vector of
            pointers to pointer to elements of the mesh that should be updated after a
            possible vector realloc.
            \sa PointerUpdater
            \param m the mesh to be modified
            \param n the number of elements to be added
            \param pu  a PointerUpdater initialized so that it can be used to update pointers to edges that could have become invalid after this adding.
            \retval the iterator to the first element added.
            */
  static HEdgeIterator AddHEdges(MeshType &m, size_t n, PointerUpdater<HEdgePointer> &pu)
  {
    HEdgeIterator last;
    if(n == 0) return m.hedge.end();
    pu.Clear();
    if(m.hedge.empty()) pu.oldBase=0;  // if the vector is empty we cannot find the last valid element
    else {
      pu.oldBase=&*m.hedge.begin();
      pu.oldEnd=&m.hedge.back()+1;
    }

    m.hedge.resize(m.hedge.size()+n);
    m.hn+=int(n);

    pu.newBase = &*m.hedge.begin();
    pu.newEnd =  &m.hedge.back()+1;

    if(pu.NeedUpdate())
    {
      if(HasFHAdjacency(m)) {
        for (FaceIterator fi=m.face.begin(); fi!=m.face.end(); ++fi)
        {
          if(!(*fi).IsD() && (*fi).FHp())
            pu.Update((*fi).FHp());
        }
      }
      if(HasVHAdjacency(m)) {
        for (VertexIterator vi=m.vert.begin(); vi!=m.vert.end(); ++vi)
          if(!(*vi).IsD() && (*vi).cVHp()!=0)
              pu.Update((*vi).VHp());
      }
      if(HasEHAdjacency(m)) {
        for (EdgeIterator ei=m.edge.begin(); ei!=m.edge.end(); ++ei)
          if(!(*ei).IsD() && (*ei).cEHp()!=0)
              pu.Update((*ei).EHp());
      }

      int ii = 0;
      HEdgeIterator hi = m.hedge.begin();
      while(ii < m.hn - int(n))// cycle on all the faces except the new ones
      {
        if(!(*hi).IsD())
        {
          if(HasHNextAdjacency(m)) pu.Update((*hi).HNp());
          if(HasHPrevAdjacency(m)) pu.Update((*hi).HPp());
          if(HasHOppAdjacency(m)) pu.Update((*hi).HOp());
          ++ii;
        }
        ++hi;
      }
    }
    size_t siz = (size_t)(m.hedge.size()-n);

    last = m.hedge.begin();
    advance(last,siz);

    return last;// deve restituire l'iteratore alla prima faccia aggiunta;
  }

  /** Function to add n vertices to the mesh.
            First wrapper, with no parameters
            */
  static HEdgeIterator AddHEdges(MeshType &m, size_t n)
  {
    PointerUpdater<HEdgePointer> pu;
    return AddHEdges(m, n,pu);
  }

  /** Function to add n vertices to the mesh.
            Second Wrapper, with a vector of vertex pointers to be updated.
            */
  static HEdgeIterator AddHEdges(MeshType &m, size_t n, std::vector<HEdgePointer*> &local_vec)
  {
    PointerUpdater<HEdgePointer> pu;
    HEdgeIterator v_ret =  AddHEdges(m, n,pu);

    typename std::vector<HEdgePointer *>::iterator ei;
    for(ei=local_vec.begin();ei!=local_vec.end();++ei)
      pu.Update(**ei);
    return v_ret;
  }

  /* +++++++++++++++ Add Faces ++++++++++++++++ */

  /** Function to add a face to the mesh and initializing it with the three given VertexPointers
            */
  static FaceIterator AddFace(MeshType &m, VertexPointer v0, VertexPointer v1, VertexPointer v2)
  {
    assert(m.vert.size()>0);
    assert((v0!=v1) && (v1!=v2) && (v0!=v2));
    assert(v0>=&m.vert.front() && v0<=&m.vert.back());
    assert(v1>=&m.vert.front() && v1<=&m.vert.back());
    assert(v2>=&m.vert.front() && v2<=&m.vert.back());
    PointerUpdater<FacePointer> pu;
    FaceIterator fi = AddFaces(m,1,pu);
    fi->Alloc(3);
    fi->V(0)=v0;
    fi->V(1)=v1;
    fi->V(2)=v2;
    return fi;
  }

  /** Function to add a face to the mesh and initializing it with three indexes
            */
  static FaceIterator AddFace(MeshType &m, size_t v0, size_t v1, size_t v2)
  {
    assert((v0!=v1) && (v1!=v2) && (v0!=v2));
    assert(v0>=0 && v0<m.vert.size());
    assert(v1>=0 && v1<m.vert.size());
    assert(v2>=0 && v2<m.vert.size());
    return AddFace(m,&(m.vert[v0]),&(m.vert[v1]),&(m.vert[v2]));
  }
  /** Function to add a face to the mesh and initializing it with the three given coords
            */
  static FaceIterator AddFace(MeshType &m, CoordType p0, CoordType p1, CoordType p2)
  {
    VertexIterator vi = AddVertices(m,3);
    FaceIterator fi = AddFaces(m,1);
    fi->Alloc(3);
    vi->P()=p0;
    fi->V(0)=&*vi++;
    vi->P()=p1;
    fi->V(1)=&*vi++;
    vi->P()=p2;
    fi->V(2)=&*vi;
    return fi;
  }

  /** Function to add a quad face to the mesh and initializing it with the four given VertexPointers
             *
             * Note that this function add a single polygonal face if the mesh has polygonal info or two tris with the corresponding faux bit set in the standard common case of a triangular mesh.
            */
  static FaceIterator AddQuadFace(MeshType &m, VertexPointer v0, VertexPointer v1, VertexPointer v2, VertexPointer v3)
  {
    assert(m.vert.size()>0);
    assert(v0>=&m.vert.front() && v0<=&m.vert.back());
    assert(v1>=&m.vert.front() && v1<=&m.vert.back());
    assert(v2>=&m.vert.front() && v2<=&m.vert.back());
    assert(v3>=&m.vert.front() && v3<=&m.vert.back());
    PointerUpdater<FacePointer> pu;
    if(FaceType::HasPolyInfo())
    {
      FaceIterator fi = AddFaces(m,1,pu);
      fi->Alloc(4);
      fi->V(0)=v0; fi->V(1)=v1;
      fi->V(2)=v2; fi->V(3)=v3;
      return fi;
    }
    else
    {
      FaceIterator fi = AddFaces(m,2,pu);
      fi->Alloc(3); fi->V(0)=v0; fi->V(1)=v1; fi->V(2)=v2;
      fi->SetF(2);
      ++fi;
      fi->Alloc(3); fi->V(0)=v0; fi->V(1)=v2; fi->V(2)=v3;
      fi->SetF(0);
      return fi;
    }
  }
  /** \brief Function to add n faces to the mesh.
            First wrapper, with no parameters
            */
  static FaceIterator AddFaces(MeshType &m, size_t n)
  {
    PointerUpdater<FacePointer> pu;
    return AddFaces(m,n,pu);
  }

  /** \brief Function to add n faces to the mesh.
            Second Wrapper, with a vector of face pointer to be updated.
            */
  static FaceIterator AddFaces(MeshType &m, size_t n,std::vector<FacePointer *> &local_vec)
  {
    PointerUpdater<FacePointer> pu;
    FaceIterator f_ret= AddFaces(m,n,pu);

    typename std::vector<FacePointer *>::iterator fi;
    for(fi=local_vec.begin();fi!=local_vec.end();++fi)
      pu.Update(**fi);
    return f_ret;
  }
  
  /** \brief Function to add n faces to the mesh getting indexes
   *  from a (fn, 3) eigen matrix of int  .
            */
  static FaceIterator AddFaces(MeshType &m, const Eigen::MatrixXi &fm)
  {
      PointerUpdater<FacePointer> pu;
      FaceIterator f_start = AddFaces(m,fm.rows(),pu);
      FaceIterator f_ret=f_start;
      for(int i=0;i<fm.rows();++i)
      {
          f_start->V(0)= &m.vert[fm(i,0)];
          f_start->V(1)= &m.vert[fm(i,1)];
          f_start->V(2)= &m.vert[fm(i,2)];
          ++f_start;
      }
      return f_ret;
  }

  /** \brief Function to add n faces to the mesh.
            This is the only full featured function that is able to manage correctly
            all the official internal pointers of the mesh (like the VF and FF adjacency relations)
            \warning Calling this function can cause the invalidation of any not-managed FacePointer
            just because we resize the face vector.
            If you have such pointers you need to update them by mean of the PointerUpdater object.
            \sa PointerUpdater
            \param m the mesh to be modified
            \param n the number of elements to be added
            \param pu  a PointerUpdater initialized so that it can be used to update pointers to edges that could have become invalid after this adding.
            \retval the iterator to the first element added.
            */
  static FaceIterator AddFaces(MeshType &m, size_t n, PointerUpdater<FacePointer> &pu)
  {
    pu.Clear();
    if(n == 0) return m.face.end();
    if(!m.face.empty()) // if the vector is empty we cannot find the last valid element
    {
      pu.oldBase=&*m.face.begin();
      pu.oldEnd=&m.face.back()+1;
    }
    // The actual resize
    m.face.resize(m.face.size()+n);
    m.fn+=int(n);

    size_t siz=(size_t)(m.face.size()-n);
    FaceIterator firstNewFace = m.face.begin();
    advance(firstNewFace,siz);

    typename std::set<PointerToAttribute>::iterator ai;
    for(ai = m.face_attr.begin(); ai != m.face_attr.end(); ++ai)
      ((PointerToAttribute)(*ai)).Resize(m.face.size());

    pu.newBase = &*m.face.begin();
    pu.newEnd  = &m.face.back()+1;

    if(pu.NeedUpdate())
    {
      if(HasFFAdjacency(m))
      {  // cycle on all the faces except the new ones
        for(FaceIterator fi=m.face.begin();fi!=firstNewFace;++fi)
          if(!(*fi).IsD())
            for(int i  = 0; i < (*fi).VN(); ++i)
              if ((*fi).cFFp(i)!=0) pu.Update((*fi).FFp(i));
      }

      if(HasPerVertexVFAdjacency(m) && HasPerFaceVFAdjacency(m))
      {  // cycle on all the faces except the new ones
        for(FaceIterator fi=m.face.begin();fi!=firstNewFace;++fi)
          if(!(*fi).IsD())
            for(int i = 0; i < (*fi).VN(); ++i)
              if ((*fi).cVFp(i)!=0) pu.Update((*fi).VFp(i));

        for (VertexIterator vi=m.vert.begin(); vi!=m.vert.end(); ++vi)
          if(!(*vi).IsD() && (*vi).cVFp()!=0)
            pu.Update((*vi).VFp());
      }

      if(HasEFAdjacency(m))
      {
        for (EdgeIterator ei=m.edge.begin(); ei!=m.edge.end(); ++ei)
          if(!(*ei).IsD() && (*ei).cEFp()!=0)
            pu.Update((*ei).EFp());
      }

      if(HasHFAdjacency(m))
      {
        for (HEdgeIterator hi=m.hedge.begin(); hi!=m.hedge.end(); ++hi)
          if(!(*hi).IsD() && (*hi).cHFp()!=0)
            pu.Update((*hi).HFp());
      }
    }
    return firstNewFace;
  }

 //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  //:::::::::::::::::TETRAS ADDER FUNCTIONS:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  /** \brief Function to add n tetras to the mesh.
            This is the only full featured function that is able to manage correctly
            all the official internal pointers of the mesh (like the VT and TT adjacency relations)
            \warning Calling this function can cause the invalidation of any not-managed TetraPointer
            just because we resize the face vector.
            If you have such pointers you need to update them by mean of the PointerUpdater object.
            \sa PointerUpdater
            \param m the mesh to be modified
            \param n the number of elements to be added
            \param pu  a PointerUpdater initialized so that it can be used to update pointers to tetras that could have become invalid after this adding.
            \retval the iterator to the first element added.
            */
  static TetraIterator AddTetras(MeshType &m, size_t n, PointerUpdater<TetraPointer> &pu)
  {
    //nothing to do
    if (n == 0)
      return m.tetra.end();

    //prepare the pointerupdater info
    pu.Clear();
    if (m.tetra.empty())
      pu.oldBase = 0;
    else
    {
      pu.oldBase = &*m.tetra.begin();
      pu.oldEnd  = &m.tetra.back() + 1;
    }

    //resize the tetra list and update tetra count
    m.tetra.resize(m.tetra.size() + n);
    m.tn += n;

    //get the old size and advance to the first new tetrahedron position
    size_t oldSize = (size_t)(m.tetra.size() - n);

    TetraIterator firstNewTetra = m.tetra.begin();
    advance(firstNewTetra, oldSize);

    //for each attribute make adapt the list size
    typename std::set<typename MeshType::PointerToAttribute>::iterator ai;
    for (ai = m.tetra_attr.begin(); ai != m.tetra_attr.end(); ++ai)
      ((typename MeshType::PointerToAttribute)(*ai)).Resize(m.tetra.size());

    //do the update
    pu.newBase = &*m.tetra.begin();
    pu.newEnd  = &m.tetra.back() + 1;
    if (pu.NeedUpdate())
    {
      if (HasVTAdjacency(m))
      {
        for (VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
          if (!vi->IsD())
            pu.Update(vi->VTp());

        for (TetraIterator ti = m.tetra.begin(); ti != m.tetra.end(); ++ti)
          if (!ti->IsD())
          {
            pu.Update(ti->VTp(0));
            pu.Update(ti->VTp(1));
            pu.Update(ti->VTp(2));
            pu.Update(ti->VTp(3));
          } 
      }

      //do edge and face adjacency
      if (HasTTAdjacency(m))
        for (TetraIterator ti = m.tetra.begin(); ti != m.tetra.end(); ++ti)
          if (!ti->IsD())
          {
            pu.Update(ti->TTp(0));
            pu.Update(ti->TTp(1));
            pu.Update(ti->TTp(2));
            pu.Update(ti->TTp(3));
          }
    }

    return firstNewTetra;
  }

//TODO: ADD 4 FACES then add tetra
  /** Function to add a face to the mesh and initializing it with the three given VertexPointers
            */
  static TetraIterator AddTetra(MeshType &m, VertexPointer v0, VertexPointer v1, VertexPointer v2, VertexPointer v3)
  {
    assert(m.vert.size() > 0);
    assert((v0 != v1) && (v0 != v2) && (v0 != v3) && (v1 != v2) && (v1 != v3) && (v2 != v3));
    assert(v0 >= &m.vert.front() && v0 <= &m.vert.back());
    assert(v1 >= &m.vert.front() && v1 <= &m.vert.back());
    assert(v2 >= &m.vert.front() && v2 <= &m.vert.back());
    assert(v3 >= &m.vert.front() && v3 <= &m.vert.back());

//    AddFace(m, v0, v1, v2);
//    AddFace(m, v0, v3, v1);
//    AddFace(m, v0, v2, v3);
//    AddFace(m, v1, v3, v2);

    PointerUpdater<TetraPointer> pu;
    TetraIterator ti = AddTetras(m, 1, pu);
    ti->V(0) = v0;
    ti->V(1) = v1;
    ti->V(2) = v2;
    ti->V(3) = v3;
    return ti;
  }

  /** Function to add a face to the mesh and initializing it with three indexes
            */
  static TetraIterator AddTetra(MeshType &m, const size_t v0, const size_t v1, const size_t v2, const size_t v3)
  {
    assert(m.vert.size() > 0);
    assert((v0 != v1) && (v0 != v2) && (v0 != v3) && (v1 != v2) && (v1 != v3) && (v2 != v3));
    assert(v0 >= 0 && v0 < m.vert.size());
    assert(v1 >= 0 && v1 < m.vert.size());
    assert(v2 >= 0 && v2 < m.vert.size());
    assert(v3 >= 0 && v3 < m.vert.size());

    return AddTetra(m, &(m.vert[v0]), &(m.vert[v1]), &(m.vert[v2]), &(m.vert[v3]));
  }
  /** Function to add a face to the mesh and initializing it with the three given coords
            */
  static TetraIterator AddTetra(MeshType &m, const CoordType & p0, const CoordType & p1, const CoordType & p2, const CoordType & p3)
  {
    VertexIterator vi = AddVertices(m, 4);

    VertexPointer v0 = &*vi++;
    VertexPointer v1 = &*vi++;
    VertexPointer v2 = &*vi++;
    VertexPointer v3 = &*vi++;

    v0->P() = p0;
    v1->P() = p1;
    v2->P() = p2;
    v3->P() = p3;

    return AddTetra(m, v0, v1, v2, v3);
  }

// //requires no duplicate vertices on faces you use
//   static TetraIterator AddTetra(MeshType &m, const FaceType & f0, const FaceType & f1, const FaceType & f2, const FaceType & f3)
//   {
//     assert(m.face.size() > 0);
//     assert((f0 != f1) && (f0 != f2) && (f0 != f3) && (f1 != f2) && (f1 != f3) && (f2 != f3));
//     assert(f1 >= 0 && f1 < m.face.size());
//     assert(f2 >= 0 && f2 < m.face.size());
//     assert(f3 >= 0 && f3 < m.face.size());
//     assert(f0 >= 0 && f0 < m.face.size());
//     //TODO: decide if you want to address this like this
//     //ERROR: can't use position...so..could force to have no dup verts..and use pointers or avoid this kind of thing
//     assert(f0.V(0) == f1.V(0) && f0.V(0) == f2.V(0) &&   //v0
//             f0.V(1) == f1.V(2) && f0.V(1) == f3.V(0) &&  //v1
//             f0.V(2) == f2.V(1) && f0.V(2) == f3.V(2) &&  //v2
//             f1.V(1) == f2.V(2) && f1.V(1) == f3.V(1) )   //v3

//     //add a tetra...and set vertices correctly
//     PointerUpdater<TetraPointer> pu;
//     TetraIterator ti = AddTetras(m, 1, pu);
//     ti->V(0) = f0.V(0);
//     ti->V(1) = f0.V(1);
//     ti->V(2) = f0.V(2);
//     ti->V(3) = f1.V(1);

//     return ti;
//   }

  /** \brief Function to add n faces to the mesh.
            First wrapper, with no parameters
            */
  static TetraIterator AddTetras(MeshType &m, size_t n)
  {
    PointerUpdater<TetraPointer> pu;
    return AddTetras(m, n, pu);
  }

  /** \brief Function to add n faces to the mesh.
            Second Wrapper, with a vector of face pointer to be updated.
            */
  static TetraIterator AddTetras(MeshType &m, size_t n, std::vector<TetraPointer *> &local_vec)
  {
    PointerUpdater<TetraPointer> pu;
    TetraIterator t_ret = AddTetras(m, n, pu);

    typename std::vector<TetraPointer *>::iterator ti;
    for (ti = local_vec.begin(); ti != local_vec.end(); ++ti)
      pu.Update(**ti);
    return t_ret;
  }

  /* +++++++++++++++ Deleting  ++++++++++++++++ */

  /** Function to delete a face from the mesh.
            NOTE: THIS FUNCTION ALSO UPDATE FN
            */
  static void DeleteFace(MeshType &m, FaceType &f)
  {
    assert(&f >= &m.face.front() && &f <= &m.face.back());
    assert(!f.IsD());
    f.Dealloc();
    f.SetD();
    --m.fn;
  }

  /** Function to delete a vertex from the mesh.
            NOTE: THIS FUNCTION ALSO UPDATE vn
            */
  static void DeleteVertex(MeshType &m, VertexType &v)
  {
    assert(&v >= &m.vert.front() && &v <= &m.vert.back());
    assert(!v.IsD());
    v.SetD();
    --m.vn;
  }

  /** Function to delete an edge from the mesh.
            NOTE: THIS FUNCTION ALSO UPDATE en
            */
  static void DeleteEdge(MeshType &m, EdgeType &e)
  {
    assert(&e >= &m.edge.front() && &e <= &m.edge.back());
    assert(!e.IsD());
    e.SetD();
    --m.en;
  }

  /** Function to delete a hedge from the mesh.
            NOTE: THIS FUNCTION ALSO UPDATE en
           */
  static void DeleteHEdge(MeshType &m, HEdgeType &h)
  {
    assert(&h >= &m.hedge.front() && &h <= &m.hedge.back());
    assert(!h.IsD());
    h.SetD();
    --m.hn;
  }

  /** Function to delete a tetra from the mesh.
            NOTE: THIS FUNCTION ALSO UPDATE tn
           */
  static void DeleteTetra(MeshType &m, TetraType &t)
  {
    assert(&t >= &m.tetra.front() && &t <= &m.tetra.back());
    assert(!t.IsD());
    t.SetD();
    --m.tn;
  }

  /*
            Function to rearrange the vertex vector according to a given index permutation
            the permutation is vector such that after calling this function

                            m.vert[ newVertIndex[i] ] = m.vert[i];

            e.g. newVertIndex[i] is the new index of the vertex i

           */
  static void PermutateVertexVector(MeshType &m, PointerUpdater<VertexPointer> &pu)
  {
    if(m.vert.empty()) return;
    for(size_t i=0;i<m.vert.size();++i)
    {
      if(pu.remap[i]<size_t(m.vn))
      {
        assert(!m.vert[i].IsD());
        m.vert[ pu.remap [i] ].ImportData(m.vert[i]);
        if(HasVFAdjacency(m))
        {
          if (m.vert[i].IsVFInitialized())
          {
            m.vert[ pu.remap[i] ].VFp() = m.vert[i].cVFp();
            m.vert[ pu.remap[i] ].VFi() = m.vert[i].cVFi();
          }
          else m.vert [ pu.remap[i] ].VFClear();
        }
        if(HasVEAdjacency(m))
        {
          if (m.vert[i].IsVEInitialized())
          {
            m.vert[ pu.remap[i] ].VEp() = m.vert[i].cVEp();
            m.vert[ pu.remap[i] ].VEi() = m.vert[i].cVEi();
          }
          else m.vert [ pu.remap[i] ].VEClear();
        }
        if (HasVTAdjacency(m))
        {
          if (m.vert[i].IsVTInitialized())
          {
            m.vert[ pu.remap[i] ].VTp() = m.vert[i].cVTp();
            m.vert[ pu.remap[i] ].VTi() = m.vert[i].cVTi();
          }
          else m.vert[ pu.remap[i] ].VTClear();
        }
      }
    }

    // reorder the optional atttributes in m.vert_attr to reflect the changes
    ReorderAttribute(m.vert_attr,pu.remap,m);

    // setup the pointer updater
    pu.oldBase  = &m.vert[0];
    pu.oldEnd = &m.vert.back()+1;

    // resize
    m.vert.resize(m.vn);

    // setup the pointer updater
    pu.newBase  = (m.vert.empty())?0:&m.vert[0];
    pu.newEnd = (m.vert.empty())?0:&m.vert.back()+1;

    // resize the optional atttributes in m.vert_attr to reflect the changes
    ResizeAttribute(m.vert_attr,m.vn,m);

    // Loop on the face to update the pointers FV relation (vertex refs)
    for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
      if(!(*fi).IsD())
        for(int i=0;i<fi->VN();++i)
        {
          size_t oldIndex = (*fi).V(i) - pu.oldBase;
          assert(pu.oldBase <= (*fi).V(i) && oldIndex < pu.remap.size());
          (*fi).V(i) = pu.newBase+pu.remap[oldIndex];
        }
    // Loop on the tetras to update the pointers TV relation (vertex refs)
    for(TetraIterator ti = m.tetra.begin(); ti != m.tetra.end(); ++ti)
      if(!(*ti).IsD())
        for(int i = 0; i < 4; ++i)
        {
          size_t oldIndex = (*ti).V(i) - pu.oldBase;
          assert(pu.oldBase <= (*ti).V(i) && oldIndex < pu.remap.size());
          (*ti).V(i) = pu.newBase+pu.remap[oldIndex];
        }
    // Loop on the edges to update the pointers EV relation (vertex refs)
    // if(HasEVAdjacency(m))
      for(EdgeIterator ei=m.edge.begin();ei!=m.edge.end();++ei)
        if(!(*ei).IsD())
        {
          pu.Update((*ei).V(0));
          pu.Update((*ei).V(1));
        }
  }

  static void CompactEveryVector(MeshType &m)
  {
    CompactVertexVector(m);
    CompactEdgeVector(m);
    CompactFaceVector(m);
    CompactTetraVector(m);
  }

  /*!
        \brief Compact vector of vertices removing deleted elements.
        Deleted elements are put to the end of the vector and the vector is resized. Order between elements is preserved but not their position (hence the PointerUpdater)
        After calling this function the \c IsD() test in the scanning a vector, is no more necessary.

        \warning It should not be called when TemporaryData is active (but works correctly if attributes are present)
        */
  static void CompactVertexVector( MeshType &m,   PointerUpdater<VertexPointer> &pu   )
  {
    // If already compacted fast return please!
    if(m.vn==(int)m.vert.size()) return;

    // newVertIndex [ <old_vert_position> ] gives you the new position of the vertex in the vector;
    pu.remap.resize( m.vert.size(),std::numeric_limits<size_t>::max() );

    size_t pos=0;
    size_t i=0;

    for(i=0;i<m.vert.size();++i)
    {
      if(!m.vert[i].IsD())
      {
        pu.remap[i]=pos;
        ++pos;
      }
    }
    assert((int)pos==m.vn);

    PermutateVertexVector(m, pu);
  }

  /*! \brief Wrapper without the PointerUpdater. */
  static void CompactVertexVector( MeshType &m  ) {
    PointerUpdater<VertexPointer>  pu;
    CompactVertexVector(m,pu);
  }

  /*!
    \brief Compact vector of edges removing deleted elements.

    Deleted elements are put to the end of the vector and the vector is resized. Order between elements is preserved but not their position (hence the PointerUpdater)
    After calling this function the \c IsD() test in the scanning a vector, is no more necessary.

    \warning It should not be called when TemporaryData is active (but works correctly if attributes are present)
    */
  static void CompactEdgeVector( MeshType &m,   PointerUpdater<EdgePointer> &pu   )
  {
    // If already compacted fast return please!
    if(m.en==(int)m.edge.size()) return;

    // remap [ <old_edge_position> ] gives you the new position of the edge in the vector;
    pu.remap.resize( m.edge.size(),std::numeric_limits<size_t>::max() );

    size_t pos=0;
    size_t i=0;

    for(i=0;i<m.edge.size();++i)
    {
      if(!m.edge[i].IsD())
      {
        pu.remap[i]=pos;
        ++pos;
      }
    }
    assert((int)pos==m.en);

    // the actual copying of the data.
    for(size_t i=0;i<m.edge.size();++i)
    {
      if(pu.remap[i]<size_t(m.en))  // uninitialized entries in the remap vector has max_int value;
      {
        assert(!m.edge[i].IsD());
        m.edge[ pu.remap [i] ].ImportData(m.edge[i]);
        // copy the vertex reference (they are not data!)
        m.edge[ pu.remap[i] ].V(0) = m.edge[i].cV(0);
        m.edge[ pu.remap[i] ].V(1) = m.edge[i].cV(1);
        // Now just copy the adjacency pointers (without changing them, to be done later)
        if(HasVEAdjacency(m))
          //if (m.edge[i].cVEp(0)!=0)
          {
            m.edge[ pu.remap[i] ].VEp(0) = m.edge[i].cVEp(0);
            m.edge[ pu.remap[i] ].VEi(0) = m.edge[i].cVEi(0);
            m.edge[ pu.remap[i] ].VEp(1) = m.edge[i].cVEp(1);
            m.edge[ pu.remap[i] ].VEi(1) = m.edge[i].cVEi(1);
          }
        if(HasEEAdjacency(m))
//          if (m.edge[i].cEEp(0)!=0)
          {
            m.edge[ pu.remap[i] ].EEp(0) = m.edge[i].cEEp(0);
            m.edge[ pu.remap[i] ].EEi(0) = m.edge[i].cEEi(0);
            m.edge[ pu.remap[i] ].EEp(1) = m.edge[i].cEEp(1);
            m.edge[ pu.remap[i] ].EEi(1) = m.edge[i].cEEi(1);
          }
        if(HasEFAdjacency(m))
//          if (m.edge[i].cEEp(0)!=0)
          {
            m.edge[ pu.remap[i] ].EFp() = m.edge[i].cEFp();
            m.edge[ pu.remap[i] ].EFi() = m.edge[i].cEFi();
            m.edge[ pu.remap[i] ].EFp() = m.edge[i].cEFp();
            m.edge[ pu.remap[i] ].EFi() = m.edge[i].cEFi();
          }

      }
    }

    // reorder the optional attributes in m.vert_attr to reflect the changes
    ReorderAttribute(m.edge_attr, pu.remap,m);

    // setup the pointer updater
    pu.oldBase  = &m.edge[0];
    pu.oldEnd = &m.edge.back()+1;

    // THE resize
    m.edge.resize(m.en);

    // setup the pointer updater
    pu.newBase  = (m.edge.empty())?0:&m.edge[0];
    pu.newEnd = (m.edge.empty())?0:&m.edge.back()+1;

    // resize the optional atttributes in m.vert_attr to reflect the changes
    ResizeAttribute(m.edge_attr,m.en,m);

    // Loop on the vertices to update the pointers of VE relation
    if(HasVEAdjacency(m))
      for (VertexIterator vi=m.vert.begin(); vi!=m.vert.end(); ++vi)
        if(!(*vi).IsD())  pu.Update((*vi).VEp());

    // Loop on the edges to update the pointers EE VE relation
    for(EdgeIterator ei=m.edge.begin();ei!=m.edge.end();++ei)
      for(unsigned int i=0;i<2;++i)
      {
        if(HasVEAdjacency(m))
          pu.Update((*ei).VEp(i));
        if(HasEEAdjacency(m))
          pu.Update((*ei).EEp(i));
//        if(HasEFAdjacency(m))
//          pu.Update((*ei).EFp());
      }
  }

  /*! \brief Wrapper without the PointerUpdater. */
  static void CompactEdgeVector( MeshType &m  ) {
    PointerUpdater<EdgePointer>  pu;
    CompactEdgeVector(m,pu);
  }

  /*!
    \brief Compact face vector by removing deleted elements.

    Deleted elements are put to the end of the vector and the vector is resized.
    Order between elements is preserved, but not their position (hence the PointerUpdater)
    Immediately after calling this function the \c IsD() test during the scanning a vector, is no more necessary.
    \warning It should not be called when some TemporaryData is active (but works correctly if attributes are present)
    */
  static void CompactFaceVector( MeshType &m, PointerUpdater<FacePointer> &pu )
  {
    // If already compacted fast return please!
    if(m.fn==(int)m.face.size()) return;

    // newFaceIndex [ <old_face_position> ] gives you the new position of the face in the vector;
    pu.remap.resize( m.face.size(),std::numeric_limits<size_t>::max() );

    size_t pos=0;
    for(size_t i=0;i<m.face.size();++i)
    {
      if(!m.face[i].IsD())
      {
        if(pos!=i)
        {
          m.face[pos].ImportData(m.face[i]);
          if(FaceType::HasPolyInfo())
          {
            m.face[pos].Dealloc();
            m.face[pos].Alloc(m.face[i].VN());
          }
          for(int j=0;j<m.face[i].VN();++j)
            m.face[pos].V(j) = m.face[i].V(j);

          if(HasVFAdjacency(m))
            for(int j=0;j<m.face[i].VN();++j)
            {
              if (m.face[i].IsVFInitialized(j)) {
                m.face[pos].VFp(j) = m.face[i].cVFp(j);
                m.face[pos].VFi(j) = m.face[i].cVFi(j);
              }
              else m.face[pos].VFClear(j);
            }
          if(HasFFAdjacency(m))
            for(int j=0;j<m.face[i].VN();++j)
              {
                m.face[pos].FFp(j) = m.face[i].cFFp(j);
                m.face[pos].FFi(j) = m.face[i].cFFi(j);
              }
        }
        pu.remap[i]=pos;
        ++pos;
      }
    }
    assert((int)pos==m.fn);

    // reorder the optional atttributes in m.face_attr to reflect the changes
    ReorderAttribute(m.face_attr,pu.remap,m);

    FacePointer fbase=&m.face[0];

    // Loop on the vertices to correct VF relation
    if(HasVFAdjacency(m))
    {
      for (VertexIterator vi=m.vert.begin(); vi!=m.vert.end(); ++vi)
        if(!(*vi).IsD())
        {
          if ((*vi).IsVFInitialized() && (*vi).VFp()!=0 )
          {
            size_t oldIndex = (*vi).cVFp() - fbase;
            assert(fbase <= (*vi).cVFp() && oldIndex < pu.remap.size());
            (*vi).VFp() = fbase+pu.remap[oldIndex];
          }
        }
    }

    // Loop on the faces to correct VF and FF relations
    pu.oldBase  = &m.face[0];
    pu.oldEnd = &m.face.back()+1;
    for(size_t i=m.fn;i<m.face.size();++i)
      m.face[i].Dealloc();
    m.face.resize(m.fn);
    pu.newBase  = (m.face.empty())?0:&m.face[0];
    pu.newEnd = (m.face.empty())?0:&m.face.back()+1;


    // resize the optional atttributes in m.face_attr to reflect the changes
    ResizeAttribute(m.face_attr,m.fn,m);

    // now we update the various (not null) face pointers (inside VF and FF relations)
    for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
      if(!(*fi).IsD())
      {
        if(HasVFAdjacency(m))
          for(int i=0;i<(*fi).VN();++i)
            if ((*fi).IsVFInitialized(i) && (*fi).VFp(i)!=0 )
            {
              size_t oldIndex = (*fi).VFp(i) - fbase;
              assert(fbase <= (*fi).VFp(i) && oldIndex < pu.remap.size());
              (*fi).VFp(i) = fbase+pu.remap[oldIndex];
            }
        if(HasFFAdjacency(m))
          for(int i=0;i<(*fi).VN();++i)
            if ((*fi).cFFp(i)!=0)
            {
              size_t oldIndex = (*fi).FFp(i) - fbase;
			  assert(fbase <= (*fi).FFp(i) && oldIndex < pu.remap.size());
              (*fi).FFp(i) = fbase+pu.remap[oldIndex];
            }
      }



  }

  /*! \brief Wrapper without the PointerUpdater. */
  static void CompactFaceVector( MeshType &m  ) {
    PointerUpdater<FacePointer>  pu;
    CompactFaceVector(m,pu);
  }

/*!
    \brief Compact tetra vector by removing deleted elements.

    Deleted elements are put to the end of the vector and the vector is resized.
    Order between elements is preserved, but not their position (hence the PointerUpdater)
    Immediately after calling this function the \c IsD() test during the scanning a vector, is no more necessary.
    \warning It should not be called when some TemporaryData is active (but works correctly if attributes are present)
    */
  static void CompactTetraVector(MeshType & m, PointerUpdater<TetraPointer> & pu)
  {
    //nothing to do
    if (size_t(m.tn) == m.tetra.size())
      return;

    //init the remap 
    pu.remap.resize(m.tetra.size(), std::numeric_limits<size_t>::max());
    
    //cycle over all the tetras, pos is the last not D() position, I is the index
    //when pos != i and !tetra[i].IsD() => we need to compact and update adj
    size_t pos = 0;
    for (size_t i = 0; i < m.tetra.size(); ++i)
    {
      if (!m.tetra[i].IsD())
      {
        if (pos != i)
        {
          //import data 
          m.tetra[pos].ImportData(m.tetra[i]);
          //import vertex refs
          for (int j = 0; j < 4; ++j)
            m.tetra[pos].V(j) = m.tetra[i].V(j);
          //import VT adj
          if (HasVTAdjacency(m))
            for (int j = 0; j < 4; ++j)
            {
              if (m.tetra[i].IsVTInitialized(j))
              {
                m.tetra[pos].VTp(j) = m.tetra[i].VTp(j);
                m.tetra[pos].VTi(j) = m.tetra[i].VTi(j);
              }
              else
                m.tetra[pos].VTClear(j);
            }
          //import TT adj
          if (HasTTAdjacency(m))
            for (int j = 0; j < 4; ++j)
            {
              m.tetra[pos].TTp(j) = m.tetra[i].cTTp(j);
              m.tetra[pos].TTi(j) = m.tetra[i].cTTi(j);
            }
        }
        //update remapping and advance pos
        pu.remap[i] = pos;
        ++pos;
      }
    }

    assert(size_t(m.tn) == pos);
    //reorder the optional attributes in m.tetra_attr 
    ReorderAttribute(m.tetra_attr, pu.remap, m);
    // resize the optional atttributes in m.tetra_attr to reflect the changes
    ResizeAttribute(m.tetra_attr, m.tn, m);

    // Loop on the tetras to correct VT and TT relations
    pu.oldBase = &m.tetra[0];
    pu.oldEnd = &m.tetra.back() + 1;

    m.tetra.resize(m.tn);
    pu.newBase = (m.tetra.empty()) ? 0 : &m.tetra[0];
    pu.newEnd  = (m.tetra.empty()) ? 0 : &m.tetra.back() + 1;

    TetraPointer tbase = &m.tetra[0];

    //Loop on the vertices to correct VT relation (since we moved things around)
    if (HasVTAdjacency(m))
    {
      for (VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
        if (!(*vi).IsD())
        {
          if ((*vi).IsVTInitialized() && (*vi).VTp() != 0)
          {
            size_t oldIndex = (*vi).cVTp() - tbase;
            assert(tbase <= (*vi).cVTp() && oldIndex < pu.remap.size());
            (*vi).VTp() = tbase + pu.remap[oldIndex];
          }
        }
    }

    // Loop on the tetras to correct the VT and TT relations
    for (TetraIterator ti = m.tetra.begin(); ti != m.tetra.end(); ++ti)
      if (!(*ti).IsD())
      {
        //VT
        if (HasVTAdjacency(m))
          for (int i = 0; i < 4; ++i)
            if ((*ti).IsVTInitialized(i) && (*ti).VTp(i) != 0)
            {
              size_t oldIndex = (*ti).VTp(i) - tbase;
              assert(tbase <= (*ti).VTp(i) && oldIndex < pu.remap.size());
              (*ti).VTp(i) = tbase + pu.remap[oldIndex];
            }
        //TT
        if (HasTTAdjacency(m))
          for (int i = 0; i < 4; ++i)
            if ((*ti).cTTp(i) != 0)
            {
              size_t oldIndex = (*ti).TTp(i) - tbase;
              assert(tbase <= (*ti).TTp(i) && oldIndex < pu.remap.size());
              (*ti).TTp(i) = tbase + pu.remap[oldIndex];
            }
      }
  }

  /*! \brief Wrapper without the PointerUpdater. */
  static void CompactTetraVector(MeshType &m)
  {
    PointerUpdater<TetraPointer> pu;
    CompactTetraVector(m, pu);
  }


public:

  /**
   * @brief Checks if a handle to a Per-Vertex Attribute is valid
   */
  template <class ATTR_TYPE>
  static
  bool IsValidHandle( const MeshType & m,  const typename MeshType::template PerVertexAttributeHandle<ATTR_TYPE> & a){
    if(a._handle == nullptr) return false;
    for(AttrIterator i = m.vert_attr.begin(); i!=m.vert_attr.end();++i)
      if ( (*i).n_attr == a.n_attr ) return true;
    return false;
  }

  /**
   * @brief Checks if a const handle to a Per-Vertex Attribute is valid
   */
  template <class ATTR_TYPE>
  static
  bool IsValidHandle( const MeshType & m,  const typename MeshType::template ConstPerVertexAttributeHandle<ATTR_TYPE> & a){
    if(a._handle == nullptr) return false;
    for(AttrIterator i = m.vert_attr.begin(); i!=m.vert_attr.end();++i)
      if ( (*i).n_attr == a.n_attr ) return true;
    return false;
  }

  /*! \brief Add a Per-Vertex Attribute of the given ATTR_TYPE with the given name.

      No attribute with that name must exists (even of different type)
    */
  template <class ATTR_TYPE>
  static
  typename MeshType::template PerVertexAttributeHandle<ATTR_TYPE>
  AddPerVertexAttribute( MeshType & m, std::string name){
    PAIte i;
    PointerToAttribute h;
    h._name = name;
    if(!name.empty()){
      i = m.vert_attr.find(h);
      assert(i ==m.vert_attr.end() );// an attribute with this name exists
    }

    h._sizeof = sizeof(ATTR_TYPE);
    h._padding = 0;
    h._handle =   new SimpleTempData<VertContainer,ATTR_TYPE>(m.vert);
    h._type = typeid(ATTR_TYPE);
    m.attrn++;
    h.n_attr = m.attrn;
    std::pair < AttrIterator , bool> res =  m.vert_attr.insert(h);
    return typename MeshType::template PerVertexAttributeHandle<ATTR_TYPE>(res.first->_handle,res.first->n_attr );
  }

  template <class ATTR_TYPE>
  static typename MeshType::template PerVertexAttributeHandle<ATTR_TYPE>
  AddPerVertexAttribute( MeshType & m){
    return AddPerVertexAttribute<ATTR_TYPE>(m,std::string(""));
  }

  /*! \brief gives a handle to a per-vertex attribute with a given name and ATTR_TYPE
      \returns a valid handle. If the name is not empty and an attribute with that name and type exists returns a handle to it.
        Otherwise return a hanlde to a newly created.
      */
  template <class ATTR_TYPE>
  static
  typename MeshType::template PerVertexAttributeHandle<ATTR_TYPE>
  GetPerVertexAttribute( MeshType & m, std::string name = std::string("")){
    typename MeshType::template PerVertexAttributeHandle<ATTR_TYPE> h;
    if(!name.empty()){
      h =  FindPerVertexAttribute<ATTR_TYPE>(m,name);
      if(IsValidHandle(m,h))
        return h;
    }
    return AddPerVertexAttribute<ATTR_TYPE>(m,name);
  }

  /*! \brief gives a const handle to a per-vertex attribute with a given name and ATTR_TYPE
      \returns a valid handle. If the name is not empty and an attribute with that name and type exists returns a handle to it.
        Otherwise, returns an invalid handle (check it using IsValidHandle).
      */
  template <class ATTR_TYPE>
  static
  typename MeshType::template ConstPerVertexAttributeHandle<ATTR_TYPE>
  GetPerVertexAttribute( const MeshType & m, std::string name = std::string("")){
    return FindPerVertexAttribute<ATTR_TYPE>(m,name);
  }

  /*! \brief Try to retrieve an handle to an attribute with a given name and ATTR_TYPE
      \returns a invalid handle if no attribute with that name and type exists.
      */
  template <class ATTR_TYPE>
  static typename MeshType::template PerVertexAttributeHandle<ATTR_TYPE>
  FindPerVertexAttribute( MeshType & m, const std::string & name)
  {
    assert(!name.empty());
    PointerToAttribute h1; h1._name = name;
    typename std::set<PointerToAttribute > :: iterator i;

    i =m.vert_attr.find(h1);
    if(i!=m.vert_attr.end())
      if((*i)._sizeof == sizeof(ATTR_TYPE) ){
        if(	(*i)._padding != 0 ){
          PointerToAttribute attr = (*i);						// copy the PointerToAttribute
          m.vert_attr.erase(i);						// remove it from the set
          FixPaddedPerVertexAttribute<ATTR_TYPE>(m,attr);
          std::pair<AttrIterator,bool> new_i = m.vert_attr.insert(attr);	// insert the modified PointerToAttribute
          assert(new_i.second);
          i = new_i.first;
        }
        return typename MeshType::template PerVertexAttributeHandle<ATTR_TYPE>((*i)._handle,(*i).n_attr);
      }
    return typename MeshType:: template PerVertexAttributeHandle<ATTR_TYPE>(nullptr,0);
  }

  /**
   * @brief Try to retrieve a const handle to an attribute with a given name
   * and ATTR_TYPE, from the given const mesh.
   * If not found, an invalid handle will be returned.
   * Check it with the function IsValidHandle
   */
  template <class ATTR_TYPE>
  static typename MeshType::template ConstPerVertexAttributeHandle<ATTR_TYPE>
  FindPerVertexAttribute( const MeshType & m, const std::string & name)
  {
    if(!name.empty()){
      PointerToAttribute h1; h1._name = name;
      typename std::set<PointerToAttribute > :: iterator i;

      i =m.vert_attr.find(h1);
      if(i!=m.vert_attr.end()){
        if((*i)._sizeof == sizeof(ATTR_TYPE) ){
          return typename MeshType::template ConstPerVertexAttributeHandle<ATTR_TYPE>((*i)._handle,(*i).n_attr);
        }
      }
    }
    return typename MeshType:: template ConstPerVertexAttributeHandle<ATTR_TYPE>(nullptr,0);
  }

  /*! \brief query the mesh for all the attributes per vertex
      \returns the name of all attributes with a non-empy name.
      */
  template <class ATTR_TYPE>
  static void GetAllPerVertexAttribute(const MeshType & m, std::vector<std::string> &all){
    all.clear();
    typename std::set<PointerToAttribute > ::const_iterator i;
    for(i = m.vert_attr.begin(); i != m.vert_attr.end(); ++i )
      if(!(*i)._name.empty())
      {
        typename MeshType:: template ConstPerVertexAttributeHandle<ATTR_TYPE> hh;
        hh = Allocator<MeshType>:: template  FindPerVertexAttribute <ATTR_TYPE>(m,(*i)._name);
        if(IsValidHandle<ATTR_TYPE>(m,hh))
          all.push_back((*i)._name);
      }
  }

  template <class ATTR_TYPE>
  static
  void
  ClearPerVertexAttribute( MeshType & m,typename MeshType::template PerVertexAttributeHandle<ATTR_TYPE> & h, const  ATTR_TYPE & initVal = ATTR_TYPE()){
    typename std::set<PointerToAttribute > ::iterator i;
    for( i = m.vert_attr.begin(); i !=  m.vert_attr.end(); ++i)
      if( (*i)._handle == h._handle ){
        for(typename MeshType::VertexIterator vi = m.vert.begin(); vi != m.vert.end(); ++vi)
          h[vi] = initVal;
        return;}
    assert(0);
  }

  /*! \brief If  the per-vertex attribute exists, delete it.
    */
  template <class ATTR_TYPE>
  static
  void
  DeletePerVertexAttribute( MeshType & m,typename MeshType::template PerVertexAttributeHandle<ATTR_TYPE> & h){
    typename std::set<PointerToAttribute > ::iterator i;
    for( i = m.vert_attr.begin(); i !=  m.vert_attr.end(); ++i)
      if( (*i)._handle == h._handle ){
        delete ((SimpleTempData<VertContainer,ATTR_TYPE>*)(*i)._handle);
        m.vert_attr.erase(i);
        return;}
  }

  // Generic DeleteAttribute.
  // It must not crash if you try to delete a non existing attribute,
  // because you do not have a way of asking for a handle of an attribute for which you do not know the type.
  static
  bool DeletePerVertexAttribute( MeshType & m,  std::string name){
    AttrIterator i;
    PointerToAttribute h1; h1._name = name;
    i = m.vert_attr.find(h1);
    if(i==m.vert_attr.end()) return false;
    delete ((SimpleTempDataBase*)(*i)._handle);
    m.vert_attr.erase(i);
    return true;
  }



  /// Per Edge Attributes
  template <class ATTR_TYPE>
  static
  bool IsValidHandle( const MeshType & m,  const typename MeshType::template PerEdgeAttributeHandle<ATTR_TYPE> & a){
    if(a._handle == nullptr) return false;
    for(AttrIterator i = m.edge_attr.begin(); i!=m.edge_attr.end();++i)
      if ( (*i).n_attr == a.n_attr ) return true;
    return false;
  }

  template <class ATTR_TYPE>
  static
  bool IsValidHandle( const MeshType & m,  const typename MeshType::template ConstPerEdgeAttributeHandle<ATTR_TYPE> & a){
    if(a._handle == nullptr) return false;
    for(AttrIterator i = m.edge_attr.begin(); i!=m.edge_attr.end();++i)
      if ( (*i).n_attr == a.n_attr ) return true;
    return false;
  }

  template <class ATTR_TYPE>
  static
  typename MeshType::template PerEdgeAttributeHandle<ATTR_TYPE>
  AddPerEdgeAttribute( MeshType & m, std::string name){
    PAIte i;
    PointerToAttribute h;
    h._name = name;
    if(!name.empty()){
      i = m.edge_attr.find(h);
      assert(i ==m.edge_attr.end() );// an attribute with this name exists
    }
    h._sizeof = sizeof(ATTR_TYPE);
    h._padding = 0;
    //		h._typename = typeid(ATTR_TYPE).name();
    h._handle =  new SimpleTempData<EdgeContainer,ATTR_TYPE>(m.edge);
    h._type = typeid(ATTR_TYPE);
    m.attrn++;
    h.n_attr = m.attrn;
    std::pair < AttrIterator , bool> res =  m.edge_attr.insert(h);
    return typename MeshType::template PerEdgeAttributeHandle<ATTR_TYPE>(res.first->_handle,res.first->n_attr);
  }

  template <class ATTR_TYPE>
  static
  typename MeshType::template PerEdgeAttributeHandle<ATTR_TYPE>
  AddPerEdgeAttribute( MeshType & m){
    return AddPerEdgeAttribute<ATTR_TYPE>(m,std::string(""));
  }

  /*! \brief gives a handle to a per-edge attribute with a given name and ATTR_TYPE
      \returns a valid handle. If the name is not empty and an attribute with that name and type exists returns a handle to it.
        Otherwise return a hanlde to a newly created.
      */
  template <class ATTR_TYPE>
  static
  typename MeshType::template PerEdgeAttributeHandle<ATTR_TYPE>
  GetPerEdgeAttribute( MeshType & m, std::string name = std::string("")){
    typename MeshType::template PerEdgeAttributeHandle<ATTR_TYPE> h;
    if(!name.empty()){
      h =  FindPerEdgeAttribute<ATTR_TYPE>(m,name);
      if(IsValidHandle(m,h))
        return h;
    }
    return AddPerEdgeAttribute<ATTR_TYPE>(m,name);
  }

  template <class ATTR_TYPE>
  static
  typename MeshType::template ConstPerEdgeAttributeHandle<ATTR_TYPE>
  GetPerEdgeAttribute( const MeshType & m, std::string name = std::string("")){
    return FindPerEdgeAttribute<ATTR_TYPE>(m,name);
  }

  template <class ATTR_TYPE>
  static
  typename MeshType::template PerEdgeAttributeHandle<ATTR_TYPE>
  FindPerEdgeAttribute( MeshType & m, const std::string & name){
    assert(!name.empty());
    PointerToAttribute h1; h1._name = name;
    typename std::set<PointerToAttribute > ::const_iterator i;

    i =m.edge_attr.find(h1);
    if(i!=m.edge_attr.end())
      if((*i)._sizeof == sizeof(ATTR_TYPE) ){
        if(	(*i)._padding != 0 ){
          PointerToAttribute attr = (*i);						// copy the PointerToAttribute
          m.edge_attr.erase(i);						// remove it from the set
          FixPaddedPerEdgeAttribute<ATTR_TYPE>(m,attr);
          std::pair<AttrIterator,bool> new_i = m.edge_attr.insert(attr);	// insert the modified PointerToAttribute
          assert(new_i.second);
          i = new_i.first;
        }
        return typename MeshType::template PerEdgeAttributeHandle<ATTR_TYPE>((*i)._handle,(*i).n_attr);
      }

    return typename MeshType:: template PerEdgeAttributeHandle<ATTR_TYPE>(nullptr,0);
  }

  template <class ATTR_TYPE>
  static
  typename MeshType::template ConstPerEdgeAttributeHandle<ATTR_TYPE>
  FindPerEdgeAttribute( const MeshType & m, const std::string & name){
    if(!name.empty()){
      PointerToAttribute h1; h1._name = name;
      typename std::set<PointerToAttribute > ::const_iterator i;

      i =m.edge_attr.find(h1);
      if(i!=m.edge_attr.end()){
        if((*i)._sizeof == sizeof(ATTR_TYPE) ){
          return typename MeshType::template ConstPerEdgeAttributeHandle<ATTR_TYPE>((*i)._handle,(*i).n_attr);
        }
      }
    }
    return typename MeshType:: template ConstPerEdgeAttributeHandle<ATTR_TYPE>(nullptr,0);
  }

  template <class ATTR_TYPE>
  static void GetAllPerEdgeAttribute(const MeshType & m, std::vector<std::string> &all){
    all.clear();
    typename std::set<PointerToAttribute > :: const_iterator i;
    for(i = m.edge_attr.begin(); i != m.edge_attr.end(); ++i )
      if(!(*i)._name.empty())
      {
        typename MeshType:: template ConstPerEdgeAttributeHandle<ATTR_TYPE> hh;
        hh = Allocator<MeshType>:: template  FindPerEdgeAttribute <ATTR_TYPE>(m,(*i)._name);
        if(IsValidHandle<ATTR_TYPE>(m,hh))
          all.push_back((*i)._name);
      }
  }

  /*! \brief If  the per-edge attribute exists, delete it.
      */
  template <class ATTR_TYPE>
  static
  void
  DeletePerEdgeAttribute( MeshType & m,typename MeshType::template PerEdgeAttributeHandle<ATTR_TYPE> & h){
    typename std::set<PointerToAttribute > ::iterator i;
    for( i = m.edge_attr.begin(); i !=  m.edge_attr.end(); ++i)
      if( (*i)._handle == h._handle ){
        delete ((SimpleTempData<EdgeContainer,ATTR_TYPE>*)(*i)._handle);
        m.edge_attr.erase(i);
        return;}
  }

  // Generic DeleteAttribute.
  // It must not crash if you try to delete a non existing attribute,
  // because you do not have a way of asking for a handle of an attribute for which you do not know the type.
  static
  bool DeletePerEdgeAttribute( MeshType & m,  std::string name){
    AttrIterator i;
    PointerToAttribute h1; h1._name = name;
    i = m.edge_attr.find(h1);
    if(i==m.edge_attr.end()) return false;
    delete ((SimpleTempDataBase*)(*i)._handle);
    m.edge_attr.erase(i);
    return true;
  }

  /// Per Face Attributes
  /**
   * @brief Checks if a handle to a Per-Face attribute is valid
   */
  template <class ATTR_TYPE>
  static
  bool IsValidHandle( const MeshType & m,  const typename MeshType::template PerFaceAttributeHandle<ATTR_TYPE> & a){
    if(a._handle == nullptr) return false;
    for(AttrIterator i = m.face_attr.begin(); i!=m.face_attr.end();++i)
      if ( (*i).n_attr == a.n_attr ) return true;
    return false;
  }

  /**
   * @brief Checks if a const handle to a Per-Face attribute is valid
   */
  template <class ATTR_TYPE>
  static
  bool IsValidHandle( const MeshType & m,  const typename MeshType::template ConstPerFaceAttributeHandle<ATTR_TYPE> & a){
    if(a._handle == nullptr) return false;
    for(AttrIterator i = m.face_attr.begin(); i!=m.face_attr.end();++i)
      if ( (*i).n_attr == a.n_attr ) return true;
    return false;
  }

  template <class ATTR_TYPE>
  static
  typename MeshType::template PerFaceAttributeHandle<ATTR_TYPE>
  AddPerFaceAttribute( MeshType & m, std::string name){
    PAIte i;
    PointerToAttribute h;
    h._name = name;
    if(!name.empty()){
      i = m.face_attr.find(h);
      assert(i ==m.face_attr.end() );// an attribute with this name exists
    }

    h._sizeof = sizeof(ATTR_TYPE);
    h._padding = 0;
    h._handle =   new SimpleTempData<FaceContainer,ATTR_TYPE>(m.face);
    h._type = typeid(ATTR_TYPE);
    m.attrn++;
    h.n_attr = m.attrn;
    std::pair < AttrIterator , bool> res =  m.face_attr.insert(h);
    return typename MeshType::template PerFaceAttributeHandle<ATTR_TYPE>(res.first->_handle,res.first->n_attr);
  }

  template <class ATTR_TYPE>
  static
  typename MeshType::template PerFaceAttributeHandle<ATTR_TYPE>
  AddPerFaceAttribute( MeshType & m){
    return AddPerFaceAttribute<ATTR_TYPE>(m,std::string(""));
  }

  /*! \brief gives a handle to a per-face attribute with a given name and ATTR_TYPE
      \returns a valid handle. If the name is not empty and an attribute with that name and type exists returns a handle to it.
        Otherwise return a hanlde to a newly created.
      */
  template <class ATTR_TYPE>
  static
  typename MeshType::template PerFaceAttributeHandle<ATTR_TYPE>
  GetPerFaceAttribute( MeshType & m, std::string name = std::string("")){
    typename MeshType::template PerFaceAttributeHandle<ATTR_TYPE> h;
    if(!name.empty()){
      h =  FindPerFaceAttribute<ATTR_TYPE>(m,name);
      if(IsValidHandle(m,h))
        return h;
    }
    return AddPerFaceAttribute<ATTR_TYPE>(m,name);
  }
  
  /*! \brief gives a handle to a per-face attribute with a given name and ATTR_TYPE
      \returns a valid handle. If the name is not empty and an attribute with that name and type exists returns a handle to it.
        Otherwise, returns an invalid handle (check it using IsValidHandle).
      */
  template <class ATTR_TYPE>
  static
  typename MeshType::template ConstPerFaceAttributeHandle<ATTR_TYPE>
  GetPerFaceAttribute( const MeshType & m, std::string name = std::string("")){
    return FindPerFaceAttribute<ATTR_TYPE>(m,name);
  }

  template <class ATTR_TYPE>
  static
  typename MeshType::template PerFaceAttributeHandle<ATTR_TYPE>
  FindPerFaceAttribute( MeshType & m, const std::string & name){
    assert(!name.empty());
    PointerToAttribute h1; h1._name = name;
    typename std::set<PointerToAttribute > ::iterator i;

    i =m.face_attr.find(h1);
    if(i!=m.face_attr.end())
      if((*i)._sizeof == sizeof(ATTR_TYPE) ){
        if( (*i)._padding != 0 ){
          PointerToAttribute attr = (*i);											// copy the PointerToAttribute
          m.face_attr.erase(i);											// remove it from the set
          FixPaddedPerFaceAttribute<ATTR_TYPE>(m,attr);
          std::pair<AttrIterator,bool> new_i = m.face_attr.insert(attr);	// insert the modified PointerToAttribute
          assert(new_i.second);
          i = new_i.first;
        }
        return typename MeshType::template PerFaceAttributeHandle<ATTR_TYPE>((*i)._handle,(*i).n_attr);
      }
    return typename MeshType:: template PerFaceAttributeHandle<ATTR_TYPE>(nullptr,0);
  }
  
  /**
   * @brief Try to retrieve a const handle to an attribute with a given name
   * and ATTR_TYPE, from the given const mesh.
   * If not found, an invalid handle will be returned.
   * Check it with the function IsValidHandle
   */
  template <class ATTR_TYPE>
  static
  typename MeshType::template ConstPerFaceAttributeHandle<ATTR_TYPE>
  FindPerFaceAttribute( const MeshType & m, const std::string & name){
    if(!name.empty()){
      PointerToAttribute h1; h1._name = name;
      typename std::set<PointerToAttribute > ::iterator i;

      i =m.face_attr.find(h1);
      if(i!=m.face_attr.end()){
        if((*i)._sizeof == sizeof(ATTR_TYPE) ){
          return typename MeshType::template ConstPerFaceAttributeHandle<ATTR_TYPE>((*i)._handle,(*i).n_attr);
        }
      }
    }
    return typename MeshType:: template ConstPerFaceAttributeHandle<ATTR_TYPE>(nullptr,0);
  }

  template <class ATTR_TYPE>
  static void GetAllPerFaceAttribute(const MeshType & m, std::vector<std::string> &all){
    all.clear();
    typename std::set<PointerToAttribute > :: const_iterator i;
    for(i = m.face_attr.begin(); i != m.face_attr.end(); ++i )
      if(!(*i)._name.empty())
      {
        typename MeshType:: template ConstPerFaceAttributeHandle<ATTR_TYPE> hh;
        hh = Allocator<MeshType>:: template  FindPerFaceAttribute <ATTR_TYPE>(m,(*i)._name);
        if(IsValidHandle<ATTR_TYPE>(m,hh))
          all.push_back((*i)._name);
      }
  }

  /*! \brief If  the per-face attribute exists, delete it.
    */
  template <class ATTR_TYPE>
  static void DeletePerFaceAttribute( MeshType & m,typename MeshType::template PerFaceAttributeHandle<ATTR_TYPE> & h){
    typename std::set<PointerToAttribute > ::iterator i;
    for( i = m.face_attr.begin(); i !=  m.face_attr.end(); ++i)
      if( (*i)._handle == h._handle ){
        delete ((SimpleTempData<FaceContainer,ATTR_TYPE>*)(*i)._handle);
        m.face_attr.erase(i);
        return;}

  }

  // Generic DeleteAttribute.
  // It must not crash if you try to delete a non existing attribute,
  // because you do not have a way of asking for a handle of an attribute for which you do not know the type.
  static bool DeletePerFaceAttribute( MeshType & m,  std::string name){
    AttrIterator i;
    PointerToAttribute h1; h1._name = name;
    i = m.face_attr.find(h1);
    if(i==m.face_attr.end()) return false;
    delete ((SimpleTempDataBase*)(*i)._handle);
    m.face_attr.erase(i);
    return true;
  }

  /// Per Tetra Attributes
  template <class ATTR_TYPE>
  static bool IsValidHandle(const MeshType & m, const typename MeshType::template PerTetraAttributeHandle<ATTR_TYPE> & a)
  {
    if (a._handle == nullptr)
      return false;
    for (AttrIterator i = m.tetra_attr.begin(); i != m.tetra_attr.end(); ++i)
      if ((*i).n_attr == a.n_attr)
        return true;
    return false;
  }

  template <class ATTR_TYPE>
  static bool IsValidHandle(const MeshType & m, const typename MeshType::template ConstPerTetraAttributeHandle<ATTR_TYPE> & a)
  {
    if (a._handle == nullptr)
      return false;
    for (AttrIterator i = m.tetra_attr.begin(); i != m.tetra_attr.end(); ++i)
      if ((*i).n_attr == a.n_attr)
        return true;
    return false;
  }

  template <class ATTR_TYPE>
  static typename MeshType::template PerTetraAttributeHandle<ATTR_TYPE> AddPerTetraAttribute(MeshType & m, std::string name)
  {
    PAIte i;
    PointerToAttribute h;
    h._name = name;
    if (!name.empty())
    {
      i = m.tetra_attr.find(h);
      assert(i == m.tetra_attr.end());
    }

    h._sizeof = sizeof(ATTR_TYPE);
    h._padding = 0;
    h._handle = new SimpleTempData<TetraContainer, ATTR_TYPE>(m.tetra);
    h._type = typeid(ATTR_TYPE);
    m.attrn++;
    h.n_attr = m.attrn;
    std::pair<AttrIterator, bool> res = m.tetra_attr.insert(h);
    return typename MeshType::template PerTetraAttributeHandle<ATTR_TYPE>(res.first->_handle, res.first->n_attr);
  }

  template <class ATTR_TYPE>
  static typename MeshType::template PerTetraAttributeHandle<ATTR_TYPE> AddPerTetraAttribute(MeshType &m)
  {
    return AddPerTetraAttribute<ATTR_TYPE>(m, std::string(""));
  }

  /*! \brief gives a handle to a per-tetra attribute with a given name and ATTR_TYPE
      \returns a valid handle. If the name is not empty and an attribute with that name and type exists returns a handle to it.
        Otherwise return a hanlde to a newly created.
      */
  template <class ATTR_TYPE>
  static typename MeshType::template PerTetraAttributeHandle<ATTR_TYPE> GetPerTetraAttribute(MeshType &m, std::string name = std::string(""))
  {
    typename MeshType::template PerTetraAttributeHandle<ATTR_TYPE> h;
    if (!name.empty())
    {
      h = FindPerTetraAttribute<ATTR_TYPE>(m, name);
      if (IsValidHandle(m, h))
        return h;
    }
    return AddPerTetraAttribute<ATTR_TYPE>(m, name);
  }

  template <class ATTR_TYPE>
  static typename MeshType::template ConstPerTetraAttributeHandle<ATTR_TYPE> GetPerTetraAttribute(const MeshType &m, std::string name = std::string(""))
  {
    return FindPerTetraAttribute<ATTR_TYPE>(m, name);
  }

  template <class ATTR_TYPE>
  static typename MeshType::template PerTetraAttributeHandle<ATTR_TYPE> FindPerTetraAttribute(MeshType &m, const std::string &name)
  {
    assert(!name.empty());
    PointerToAttribute h1;
    h1._name = name;
    typename std::set<PointerToAttribute>::iterator i;

    i = m.tetra_attr.find(h1);
    if (i != m.tetra_attr.end())
      if ((*i)._sizeof == sizeof(ATTR_TYPE))
      {
        if ((*i)._padding != 0)
        {
          PointerToAttribute attr = (*i); // copy the PointerToAttribute
          m.tetra_attr.erase(i);           // remove it from the set
          FixPaddedPerTetraAttribute<ATTR_TYPE>(m, attr);
          std::pair<AttrIterator, bool> new_i = m.tetra_attr.insert(attr); // insert the modified PointerToAttribute
          assert(new_i.second);
          i = new_i.first;
        }
        return typename MeshType::template PerTetraAttributeHandle<ATTR_TYPE>((*i)._handle, (*i).n_attr);
      }
    return typename MeshType::template PerTetraAttributeHandle<ATTR_TYPE>(nullptr, 0);
  }

  template <class ATTR_TYPE>
  static typename MeshType::template ConstPerTetraAttributeHandle<ATTR_TYPE> FindPerTetraAttribute(const MeshType &m, const std::string &name)
  {
    if(!name.empty()){
      PointerToAttribute h1;
      h1._name = name;
      typename std::set<PointerToAttribute>::iterator i;

      i = m.tetra_attr.find(h1);
      if (i != m.tetra_attr.end()){
        if ((*i)._sizeof == sizeof(ATTR_TYPE))
        {
          return typename MeshType::template ConstPerTetraAttributeHandle<ATTR_TYPE>((*i)._handle, (*i).n_attr);
        }
      }
    }
    return typename MeshType::template ConstPerTetraAttributeHandle<ATTR_TYPE>(nullptr, 0);
  }

  template <class ATTR_TYPE>
  static void GetAllPerTetraAttribute(const MeshType &m, std::vector<std::string> &all)
  {
    all.clear();
    typename std::set<PointerToAttribute>::const_iterator i;
    for (i = m.tetra_attr.begin(); i != m.tetra_attr.end(); ++i)
      if (!(*i)._name.empty())
      {
        typename MeshType::template ConstPerTetraAttributeHandle<ATTR_TYPE> hh;
        hh = Allocator<MeshType>::template FindPerTetraAttribute<ATTR_TYPE>(m, (*i)._name);
        if (IsValidHandle<ATTR_TYPE>(m, hh))
          all.push_back((*i)._name);
      }
  }

  /*! \brief If  the per-face attribute exists, delete it.
    */
  template <class ATTR_TYPE>
  static void DeletePerTetraAttribute(MeshType &m, typename MeshType::template PerTetraAttributeHandle<ATTR_TYPE> &h)
  {
    typename std::set<PointerToAttribute>::iterator i;
    for (i = m.tetra_attr.begin(); i != m.tetra_attr.end(); ++i)
      if ((*i)._handle == h._handle)
      {
        delete ((SimpleTempData<TetraContainer, ATTR_TYPE> *)(*i)._handle);
        m.tetra_attr.erase(i);
        return;
      }
  }

  // Generic DeleteAttribute.
  // It must not crash if you try to delete a non existing attribute,
  // because you do not have a way of asking for a handle of an attribute for which you do not know the type.
  static bool DeletePerTetraAttribute(MeshType &m, std::string name)
  {
    AttrIterator i;
    PointerToAttribute h1;
    h1._name = name;
    i = m.tetra_attr.find(h1);
    if (i == m.tetra_attr.end())
      return false;
    delete ((SimpleTempDataBase *)(*i)._handle);
    m.tetra_attr.erase(i);
    return true;
  }


  /// Per Mesh Attributes
  template <class ATTR_TYPE>
  static
  bool IsValidHandle(const MeshType & m,  const typename MeshType::template PerMeshAttributeHandle<ATTR_TYPE> & a){
    if(a._handle == nullptr) return false;
    for(AttrIterator i = m.mesh_attr.begin(); i!=m.mesh_attr.end();++i)
      if ( (*i).n_attr == a.n_attr ) return true;
    return false;
  }

  template <class ATTR_TYPE>
  static
  bool IsValidHandle(const MeshType & m,  const typename MeshType::template ConstPerMeshAttributeHandle<ATTR_TYPE> & a){
    if(a._handle == nullptr) return false;
    for(AttrIterator i = m.mesh_attr.begin(); i!=m.mesh_attr.end();++i)
      if ( (*i).n_attr == a.n_attr ) return true;
    return false;
  }

  template <class ATTR_TYPE>
  static
  typename MeshType::template PerMeshAttributeHandle<ATTR_TYPE>
  AddPerMeshAttribute( MeshType & m, std::string name){
    PAIte i;
    PointerToAttribute h;
    h._name = name;
    if(!name.empty()){
      i = m.mesh_attr.find(h);
      assert(i ==m.mesh_attr.end() );// an attribute with this name exists
    }
    h._sizeof = sizeof(ATTR_TYPE);
    h._padding = 0;
    h._handle =  new Attribute<ATTR_TYPE>();
    h._type = typeid(ATTR_TYPE);
    m.attrn++;
    h.n_attr = m.attrn;
    std::pair < AttrIterator , bool> res =  m.mesh_attr.insert(h);
    return typename MeshType::template PerMeshAttributeHandle<ATTR_TYPE>(res.first->_handle,res.first->n_attr);
  }

  /*! \brief gives a handle to a per-edge attribute with a given name and ATTR_TYPE
      \returns a valid handle. If the name is not empty and an attribute with that name and type exists returns a handle to it.
        Otherwise return a hanlde to a newly created.
      */
  template <class ATTR_TYPE>
  static
  typename MeshType::template PerMeshAttributeHandle<ATTR_TYPE>
  GetPerMeshAttribute( MeshType & m, std::string name = std::string("")){
    typename MeshType::template PerMeshAttributeHandle<ATTR_TYPE> h;
    if(!name.empty()){
      h =  FindPerMeshAttribute<ATTR_TYPE>(m,name);
      if(IsValidHandle(m,h))
        return h;
    }
    return AddPerMeshAttribute<ATTR_TYPE>(m,name);
  }

  template <class ATTR_TYPE>
  static
  typename MeshType::template ConstPerMeshAttributeHandle<ATTR_TYPE>
  GetPerMeshAttribute(const MeshType & m, std::string name = std::string("")){
    return FindPerMeshAttribute<ATTR_TYPE>(m,name);
  }

  template <class ATTR_TYPE>
  static
  typename MeshType::template PerMeshAttributeHandle<ATTR_TYPE>
  FindPerMeshAttribute( MeshType & m, const std::string & name){
    assert(!name.empty());
    PointerToAttribute h1; h1._name = name;
    typename std::set<PointerToAttribute > ::iterator i;

    i =m.mesh_attr.find(h1);
    if(i!=m.mesh_attr.end())
      if((*i)._sizeof == sizeof(ATTR_TYPE)  ){
        if(	(*i)._padding != 0 ){
          PointerToAttribute attr = (*i);											// copy the PointerToAttribute
          m.mesh_attr.erase(i);											// remove it from the set
          FixPaddedPerMeshAttribute<ATTR_TYPE>(m,attr);
          std::pair<AttrIterator,bool> new_i = m.mesh_attr.insert(attr);	// insert the modified PointerToAttribute
          assert(new_i.second);
          i = new_i.first;
        }

        return typename MeshType::template PerMeshAttributeHandle<ATTR_TYPE>((*i)._handle,(*i).n_attr);
      }

    return typename MeshType:: template PerMeshAttributeHandle<ATTR_TYPE>(nullptr,0);
  }

  template <class ATTR_TYPE>
  static
  typename MeshType::template ConstPerMeshAttributeHandle<ATTR_TYPE>
  FindPerMeshAttribute( const MeshType & m, const std::string & name){
    if (!name.empty()){
      PointerToAttribute h1; h1._name = name;
      typename std::set<PointerToAttribute > ::iterator i;
      i =m.mesh_attr.find(h1);
      if(i!=m.mesh_attr.end()){
        if((*i)._sizeof == sizeof(ATTR_TYPE)  ){
          return typename MeshType::template ConstPerMeshAttributeHandle<ATTR_TYPE>((*i)._handle,(*i).n_attr);
        }
      }
    }

    return typename MeshType:: template ConstPerMeshAttributeHandle<ATTR_TYPE>(nullptr,0);
  }

  template <class ATTR_TYPE>
  static void GetAllPerMeshAttribute(const MeshType & m, std::vector<std::string> &all){
    typename std::set<PointerToAttribute > :: iterator i;
    for(i = m.mesh_attr.begin(); i != m.mesh_attr.end(); ++i )
      if((*i)._sizeof == sizeof(ATTR_TYPE))
        all.push_back((*i)._name);
  }

  /*! \brief If  the per-mesh attribute exists, delete it.
      */
  template <class ATTR_TYPE>
  static void DeletePerMeshAttribute( MeshType & m,typename MeshType::template PerMeshAttributeHandle<ATTR_TYPE> & h){
    typename std::set<PointerToAttribute > ::iterator i;
    for( i = m.mesh_attr.begin(); i !=  m.mesh_attr.end(); ++i)
      if( (*i)._handle == h._handle ){
        delete (( Attribute<ATTR_TYPE> *)(*i)._handle);
        m.mesh_attr.erase(i);
        return;}
  }

  // Generic DeleteAttribute.
  // It must not crash if you try to delete a non existing attribute,
  // because you do not have a way of asking for a handle of an attribute for which you do not know the type.
  static bool DeletePerMeshAttribute( MeshType & m,  std::string name){
    AttrIterator i;
    PointerToAttribute h1; h1._name = name;
    i = m.mesh_attr.find(h1);
    if (i==m.mesh_attr.end())
        return false;
    delete ((SimpleTempDataBase  *)(*i)._handle);
    m.mesh_attr.erase(i);
    return true;
  }

  template <class ATTR_TYPE>
  static void FixPaddedPerVertexAttribute (MeshType & m, PointerToAttribute & pa){

    // create the container of the right type
    SimpleTempData<VertContainer,ATTR_TYPE>* _handle =  new SimpleTempData<VertContainer,ATTR_TYPE>(m.vert);

    // copy the padded container in the new one
    _handle->Resize(m.vert.size());
    for(size_t i  = 0; i < m.vert.size(); ++i){
      ATTR_TYPE * dest = &(*_handle)[i];
      char * ptr = (char*)( ((SimpleTempDataBase *)pa._handle)->DataBegin());
      ATTR_TYPE* attrptr = (ATTR_TYPE*)ptr;
      *dest = attrptr[i * pa._sizeof ];
      //memcpy((void*)dest ,
            //(void*) &(ptr[i *  pa._sizeof ]) ,sizeof(ATTR_TYPE));
    }

    // remove the padded container
    delete ((SimpleTempDataBase*) pa._handle);

    // update the pointer to data
    pa._sizeof = sizeof(ATTR_TYPE);

    // update the pointer to data
    pa._handle = _handle;

    // zero the padding
    pa._padding = 0;
  }
  template <class ATTR_TYPE>
  static void FixPaddedPerEdgeAttribute (MeshType & m, PointerToAttribute & pa){

    // create the container of the right type
    SimpleTempData<EdgeContainer,ATTR_TYPE>* _handle =  new SimpleTempData<EdgeContainer,ATTR_TYPE>(m.edge);

    // copy the padded container in the new one
    _handle->Resize(m.edge.size());
    for(size_t i  = 0; i < m.edge.size(); ++i){
      ATTR_TYPE * dest = &(*_handle)[i];
      char * ptr = (char*)( ((SimpleTempDataBase *)pa._handle)->DataBegin());
      ATTR_TYPE* attrptr = (ATTR_TYPE*)ptr;
      *dest = attrptr[i * pa._sizeof ];
      //memcpy((void*)dest ,
             //(void*) &(ptr[i *  pa._sizeof ]) ,sizeof(ATTR_TYPE));
    }

    // remove the padded container
    delete ((SimpleTempDataBase*) pa._handle);

    // update the pointer to data
    pa._sizeof = sizeof(ATTR_TYPE);

    // update the pointer to data
    pa._handle = _handle;

    // zero the padding
    pa._padding = 0;
  }

  template <class ATTR_TYPE>
  static void FixPaddedPerFaceAttribute ( MeshType & m,PointerToAttribute & pa){

    // create the container of the right type
    SimpleTempData<FaceContainer,ATTR_TYPE>* _handle =  new SimpleTempData<FaceContainer,ATTR_TYPE>(m.face);

    // copy the padded container in the new one
    _handle->Resize(m.face.size());
    for(size_t i  = 0; i < m.face.size(); ++i){
      ATTR_TYPE * dest = &(*_handle)[i];
      char * ptr = (char*)( ((SimpleTempDataBase *)pa._handle)->DataBegin());
      ATTR_TYPE* attrptr = (ATTR_TYPE*)ptr;
      *dest = attrptr[i * pa._sizeof ];
      //memcpy((void*)dest ,
      //       (void*) &(ptr[i * pa._sizeof ]) ,sizeof(ATTR_TYPE));
    }

    // remove the padded container
    delete ((SimpleTempDataBase*) pa._handle);

    // update the pointer to data
    pa._sizeof = sizeof(ATTR_TYPE);

    // update the pointer to data
    pa._handle = _handle;

    // zero the padding
    pa._padding = 0;
  }

  template <class ATTR_TYPE>
  static void FixPaddedPerTetraAttribute(MeshType &m, PointerToAttribute &pa)
  {

    // create the container of the right type
    SimpleTempData<TetraContainer, ATTR_TYPE> *_handle = new SimpleTempData<TetraContainer, ATTR_TYPE>(m.tetra);

    // copy the padded container in the new one
    _handle->Resize(m.tetra.size());
    for (size_t i = 0; i < m.tetra.size(); ++i)
    {
      ATTR_TYPE *dest = &(*_handle)[i];
      char *ptr = (char *)(((SimpleTempDataBase *)pa._handle)->DataBegin());
      ATTR_TYPE* attrptr = (ATTR_TYPE*)ptr;
      *dest = attrptr[i * pa._sizeof ];
      //memcpy((void *)dest,
             //(void *)&(ptr[i * pa._sizeof]), sizeof(ATTR_TYPE));
    }

    // remove the padded container
    delete ((SimpleTempDataBase *)pa._handle);

    // update the pointer to data
    pa._sizeof = sizeof(ATTR_TYPE);

    // update the pointer to data
    pa._handle = _handle;

    // zero the padding
    pa._padding = 0;
  }
  
  template <class ATTR_TYPE>
  static void FixPaddedPerMeshAttribute ( MeshType & /* m */,PointerToAttribute & pa){

    // create the container of the right type
    Attribute<ATTR_TYPE> * _handle =  new Attribute<ATTR_TYPE>();

    // copy the padded container in the new one
    ATTR_TYPE* dest = (ATTR_TYPE*)_handle->DataBegin();
    char* ptr = (char*)( ((Attribute<ATTR_TYPE> *)pa._handle)->DataBegin());
    ATTR_TYPE* attrptr = (ATTR_TYPE*)ptr;
    *dest = *attrptr;
    //memcpy((void*)dest ,(void*) &(ptr[0]) ,sizeof(ATTR_TYPE));

    // remove the padded container
    delete ( (Attribute<ATTR_TYPE> *) pa._handle);

    // update the pointer to data
    pa._sizeof = sizeof(ATTR_TYPE);

    // update the pointer to data
    pa._handle = _handle;

    // zero the padding
    pa._padding = 0;
  }

}; // end Allocator class


/** @} */ // end doxygen group trimesh
} // end namespace tri
} // end namespace vcg

#endif
