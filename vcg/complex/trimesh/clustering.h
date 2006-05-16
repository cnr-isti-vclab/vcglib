/****************************************************************************
 * VCGLib                                                            o o     *
 * Visual and Computer Graphics Library                            o     o   *
 *                                                                _   O  _   *
 * Copyright(C) 2006                                                \/)\/    *
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
/****************************************************************************
  History

$Log: not supported by cvs2svn $
****************************************************************************/

#ifndef __VCGLIB_CLUSTERING
#define __VCGLIB_CLUSTERING

#include<vcg/complex/trimesh/base.h>
#include <vcg/complex/trimesh/clean.h>
#include<vcg/space/triangle3.h>
#include<vcg/complex/trimesh/update/topology.h>
#include<vcg/space/index/grid_util.h>

#include <iostream>
#include <math.h>
#include <limits>

// some stuff for portable hashes...
#ifdef WIN32
 #ifndef __MINGW32__
  #include <hash_map>
  #include <hash_set>
  #define STDEXT stdext
 #else
  #include <ext/hash_map>
  #include <ext/hash_set>
  #define STDEXT __gnu_cxx
 #endif
#else
 #include <ext/hash_map>
 #include <ext/hash_set>
 #define STDEXT __gnu_cxx
#endif



namespace vcg{
namespace tri{
#define HASH_P0 73856093
#define HASH_P1 19349663
#define HASH_P2 83492791

class HashedPoint3i : public Point3i
{
public:

  const size_t Hash() const 
  {
    return (V(0)*HASH_P0 ^ V(1)*HASH_P1 ^ V(2)*HASH_P2);
  }

  operator size_t () const
  {return Hash();}
};

template<class MeshType>
class  AverageCell
{
  typedef typename MeshType::CoordType CoordType;
  typedef typename MeshType::FaceType  FaceType;
  public:
    inline void Add(MeshType &m, FaceType &f, int i)
    {
      p+=f.cV(i)->cP();
      n+=f.cV(i)->cN();
      cnt++;
    }
    AverageCell(): p(0,0,0), n(0,0,0),cnt(0){}
   CoordType p;
   CoordType n;
   int cnt;
   int id;
  CoordType Pos() const 
  {
    return p/cnt;
  }
};

/*
  Metodo di clustering
*/
template<class MeshType, class CellType, bool Selected=true>
class Clustering
{	
 public:
  typedef typename MeshType::ScalarType  ScalarType;
  typedef typename MeshType::CoordType CoordType;
  typedef typename MeshType::VertexType  VertexType;
  typedef typename MeshType::FaceType  FaceType;
  typedef typename MeshType::VertexPointer  VertexPointer;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::FaceIterator   FaceIterator;

  class SimpleTri
  {
  public: 
    CellType *v[3];
    const int ii(int i) const {return *((int *)(&(v[i])));}
    bool operator < ( const SimpleTri &p) const {
      return	(v[2]!=p.v[2])?(v[2]<p.v[2]):
				(v[1]!=p.v[1])?(v[1]<p.v[1]):
				(v[0]<p.v[0]);
	  }

    void sort()
    {
      //std::sort(&(v[0]),&(v[3]));
      if(v[0] > v[1] ) swap(v[0],v[1]); // now v0 < v1
      if(v[0] > v[2] ) swap(v[0],v[2]); // now v0 is the minimum
      if(v[1] > v[2] ) swap(v[1],v[2]); // sorted!
    }
    // Hashing Function;
    operator size_t () const
    {
      return (ii(0)*HASH_P0 ^ ii(1)*HASH_P1 ^ ii(2)*HASH_P2);
    }
  };



  void Init(Box3<ScalarType> _mbb, int _size, ScalarType _cellsize=0)
  {
    Grid.bbox=_mbb;
	///inflate the bb calculated
	  ScalarType infl=Grid.bbox.Diag()/_size;
	  Grid.bbox.min-=CoordType(infl,infl,infl);
	  Grid.bbox.max+=CoordType(infl,infl,infl);
    Grid.dim  = Grid.bbox.max - Grid.bbox.min;
    if(		_cellsize==0)
        BestDim( _size, Grid.dim, Grid.siz );
    else
      Grid.siz = Point3i::Construct(Grid.dim / _cellsize);

				// find voxel size
		Grid.voxel[0] = Grid.dim[0]/Grid.siz[0];
		Grid.voxel[1] = Grid.dim[1]/Grid.siz[1];
		Grid.voxel[2] = Grid.dim[2]/Grid.siz[2];
  }
  
  
  BasicGrid<CellType,ScalarType> Grid;
  STDEXT::hash_set<SimpleTri> TriSet;
  STDEXT::hash_map<HashedPoint3i,CellType> GridCell;
  
  void Add(MeshType &m)
  {
    FaceIterator fi;
    for(fi=m.face.begin();fi!=m.face.end();++fi)
    {
      HashedPoint3i pi;
      SimpleTri st;
      for(int i=0;i<3;++i)
      {
        Grid.PToIP((*fi).cV(i)->cP(), pi );
        st.v[i]=&(GridCell[pi]);
        st.v[i]->Add(m,*(fi),i);
      }        
      if( (st.v[0]!=st.v[1]) && (st.v[0]!=st.v[2]) && (st.v[1]!=st.v[2]) )
      {
        st.sort();
        TriSet.insert(st);
      }
    //  printf("Inserted %8i triangles, clustered to %8i tri and %i cells\n",distance(m.face.begin(),fi),TriSet.size(),GridCell.size());
    }
  }
  
  void Extract(MeshType &m)
  {
    m.Clear();
    Allocator<MeshType>::AddVertices(m,GridCell.size());
    std::hash_map<HashedPoint3i,CellType>::iterator gi;
    int i=0;
    for(gi=GridCell.begin();gi!=GridCell.end();++gi)
    {
      m.vert[i].P()=(*gi).second.Pos();
      (*gi).second.id=i;
      ++i;
    }
    Allocator<MeshType>::AddFaces(m,TriSet.size());
    std::hash_set<SimpleTri>::iterator ti;
    i=0;
    for(ti=TriSet.begin();ti!=TriSet.end();++ti)
    {
      m.face[i].V(0)=&(m.vert[(*ti).v[0]->id]);
      m.face[i].V(1)=&(m.vert[(*ti).v[1]->id]);
      m.face[i].V(2)=&(m.vert[(*ti).v[2]->id]);
      CoordType N=Normal(m.face[i]);
      int badOrient=0;
      if( N*(*ti).v[0]->n <0) ++badOrient;
      if( N*(*ti).v[1]->n <0) ++badOrient;
      if( N*(*ti).v[2]->n <0) ++badOrient;
      if(badOrient>=2)
        swap(m.face[i].V(0),m.face[i].V(1));
      i++;
    }
    
  }
}; //end class clustering 
 } // namespace tri
} // namespace vcg

#endif
