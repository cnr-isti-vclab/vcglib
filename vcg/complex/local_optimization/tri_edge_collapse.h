/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    *
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
  $Log: not supported by cvs2svn $
  Revision 1.8  2004/09/08 14:33:31  ganovelli
  *** empty log message ***


****************************************************************************/

#ifndef __VCG_DECIMATION_TRICOLLAPSE
#define __VCG_DECIMATION_CTRIOLLAPSE

#include<vcg\complex\trimesh\edge_collapse.h>
#include<vcg\simplex\face\pos.h>
#include<vcg\complex\local_optimization.h>


namespace vcg{
namespace tri{

/** \addtogroup trimesh */
/*@{*/
/// This Class is specialization of LocalModification for the edge collapse
/// It wraps the atomic operation EdgeCollapse to be used in a optimizatin routine.
/// Note that it has knowledge of the heap of the class LocalOptimization because
/// it is responsible of updating it after a collapse has been performed

template<class TriMeshType,class MYTYPE>
class TriEdgeCollapse: public LocalOptimization<TriMeshType>::LocModType , public EdgeCollapse<TriMeshType> 
{
public:
 /// static data to gather statistical information about the reasons of collapse failures
 struct FailStat {
	static int &Volume()           {static int vol=0; return vol;}
	static int &LinkConditionFace(){static int lkf=0; return lkf;}
	static int &LinkConditionEdge(){static int lke=0; return lke;}
	static int &LinkConditionVert(){static int lkv=0; return lkv;}
	static int &OutOfDate()        {static int ofd=0; return ofd;}
	static int &Border()           {static int bor=0; return bor;}
  static void Init() 
  {
   Volume()           =0;
   LinkConditionFace()=0;
   LinkConditionEdge()=0;
   LinkConditionVert()=0;
   OutOfDate()        =0;
   Border()           =0;
  }
};

  typedef	typename TriMeshType::FaceType FaceType;
  typedef	typename FaceType::VertexType VertexType;
  typedef	typename FaceType::VertexType::CoordType CoordType;
  typedef	typename TriMeshType::VertexType::ScalarType ScalarType;
  typedef vcg::face::Pos<FaceType> PosType;
  typedef typename LocalOptimization<TriMeshType>::HeapElem HeapElem;
	

protected:
	///the pos of collapse 
	PosType pos;

	///mark for up_dating
	static int& _Imark(){ static int im=0; return im;}
	
	/// priority in the heap
	ScalarType _priority;

	public:
	/// Default Constructor
		TriEdgeCollapse()
			{}
	///Constructor with postype
	 TriEdgeCollapse(PosType p, int mark)
			{    
				_Imark() = mark;
				pos=p;
				_priority = ComputePriority();
			}

		~TriEdgeCollapse()
			{}

private:


public:


  ScalarType ComputePriority()
  { 
		_priority = Distance(pos.V()->P(),pos.VFlip()->P()); 
    return _priority;
  }

  virtual const char *Info(TriMeshType &m) {
    static char buf[60];
    sprintf(buf,"collapse %i -> %i %f\n", pos.V()-&m.vert[0], pos.VFlip()-&m.vert[0],_priority);
    return buf;
  }
 
  void Execute(TriMeshType &m)
  {	
    CoordType MidPoint=(pos.V()->P()+pos.VFlip()->P())/2.0;
	int FaceDel=DoCollapse(pos, MidPoint);
    m.fn-=FaceDel;
    --m.vn;
  }
  
  
  void UpdateHeap(typename LocalOptimization<TriMeshType>::HeapType & h_ret)
  {
		_Imark()++;
		vcg::face::VFIterator<FaceType> VFi(pos.V(1)->VFp(),pos.V(1)->VFi());
    while (!VFi.End())
    {
      for (int j=0;j<3;j++)
      {
				PosType p;
				p.Set(VFi.F(),VFi.I(),VFi.f->V(VFi.z));
				h_ret.push_back(HeapElem(new MYTYPE(p,_Imark())));
				std::push_heap(h_ret.begin(),h_ret.end());
				//// update the mark of the vertices
				VFi.f->V(VFi.z)->IMark() = _Imark();
				VFi.f->V( (VFi.z+1) % 3 )->IMark() = _Imark();
      }
      ++VFi;
    }
  }

  ModifierType IsOfType(){ return TriEdgeCollapseOp;}

  bool IsFeasible(){
		return LinkConditions(pos);
	}

  bool IsUpToDate(){
    if(pos.f->IsD()) {
				++FailStat::OutOfDate();
				return false;
			}
			
    if(pos.v->IsD()) {
				++FailStat::OutOfDate();
				return false;
			}

		  VertexType *v0=pos.V();
			VertexType *v1=pos.VFlip();
			
			if(! (( (!v0->IsD()) && (!v1->IsD())) &&
							 _Imark()>=v0->IMark() &&
							 _Imark()>=v1->IMark()))
			{
				++FailStat::OutOfDate();
				return false;
			}
				return true;
	}

	virtual ScalarType Priority() const {
	return _priority;
  }

	static void Init(TriMeshType&m,typename LocalOptimization<TriMeshType>::HeapType&h_ret){
		h_ret.clear();
		typename TriMeshType::FaceIterator fi;
		for(fi = m.face.begin(); fi != m.face.end();++fi)
		if(!(*fi).IsD()){
		   for (int j=0;j<3;j++)
      {
        PosType p=PosType(&*fi,j,(*fi).V(j));
        h_ret.push_back(HeapElem(new MYTYPE(p,m.IMark())));
        //printf("Inserting in heap coll %3i ->%3i %f\n",p.V()-&m.vert[0],p.VFlip()-&m.vert[0],h_ret.back().locModPtr->Priority());
      }
		}
	}

};
}//end namespace tri
}//end namespace vcg

#endif
