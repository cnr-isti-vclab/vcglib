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
  History
****************************************************************************/

#ifndef __VCG_TETRA_DECIMATION
#define __VCG_TETRA_DECIMATION
#include<vector>
#include<algorithm>
#include<time.h>
#include<math.h>
#include<vcg/complex/tetramesh/decimation/operationsdef.h>
#include<vcg/complex/tetramesh/decimation/collapse.h>

namespace vcg{
namespace tetra{


template<class TETRA_MESH_TYPE>
class Decimation
{
public:
 
  /// The tetrahedral mesh type
  typedef	typename TETRA_MESH_TYPE TetraMeshType;
  // ///the decimator's type
  //typedef	typename Decimation<TetraMeshType> DecimatorBase;
  /// The tetrahedron type
  typedef	typename TetraMeshType::TetraType TetraType;
	/// The vertex type
	typedef	typename TetraType::VertexType VertexType;
  /// The coordinate type
	typedef	typename TetraType::VertexType::CoordType CoordType;
  /// The scalar type
  typedef	typename TetraMeshType::VertexType::ScalarType ScalarType;
  //local modification type 
  typedef typename  vcg::tetra::LocalModification<TETRA_MESH_TYPE> LocalModification;
  /// The pos type
  //typedef typename vcg::tetra::Pos<TetraType> PosType;
  typedef typename vcg::tetra::LocalModification<TETRA_MESH_TYPE>::PosType PosType;

  ///the element of the heap
  struct HeapElem
  {
    ///the modifier's type
    vcg::tetra::ModifiersType Mt;
    ///the value of priority
    ScalarType priority;
    ///the pos where the modifier operate
    PosType pos;
    ///pointer to instance of local modifier
    LocalModification*LM;
    ///temporary mark for the opration
    char Imark;
   
    HeapElem(vcg::tetra::ModifiersType Mtype,PosType p,char mark)
    {
      Mt=Mtype;
      pos=p;
      Imark=mark;
      LM=NULL;
      Instanciate();
      priority=LM->ComputePriority();
      delete(LM);
    };

    const bool operator <(const HeapElem & h) const 
    {
		  return (priority<h.priority);
	  }
//this method instaciate each different type of operation
    void Instanciate()
    {
       if (Mt==MTEdgeCollapse)
        LM=new vcg::tetra::Collapse<TetraMeshType>(pos);
    }
 
    HeapElem & operator =( const HeapElem & h) 
    {
      Mt=h.Mt;
      priority=priority;
      pos=h.pos;
      LM=h.LM;
      Imark=h.Imark;
      return (*this);
	  }


    bool IsUpToDate()
    {
    	if (!pos.T()->IsD())
		  {
        VertexType *v0=pos.T()->V(Tetra::VofE(pos.E(),0));
			  VertexType *v1=pos.T()->V(Tetra::VofE(pos.E(),1));
			  
			return (( (!v0->IsD()) && (!v1->IsD())) &&
							 Imark>=v0->IMark() &&
							 Imark>=v1->IMark());
		  }
		else
		return false;
    }

  };

  /// The heap type
  typedef typename vector<HeapElem> HeapType;

  ///the pointer to tetramesh
  TetraMeshType &tm;
  ///the heap of operations
  HeapType h;

  /// Default Constructor
	Decimation(TETRA_MESH_TYPE &_tm):tm(_tm)
		{
		};

  ~Decimation()
		{
		};


 void AddOperation(vcg::tetra::ModifiersType Mt,Pos<TetraType> pos)
  { 
     h.push_back(HeapElem(Mt,pos,tm.GetMark()));
		 push_heap( h.begin(), h.end());
  }

  void DecimationStep(int nstep)
  {
    for (int i=0;i<nstep;i++)
    {
      if( ! h.back().IsUpToDate())	
		    h.pop_back();
	    else  
      {	
        h.back().Instanciate();
        h.back().LM->ComputeError();
        LocalModification* LastMod= h.back().LM;
        h.pop_back();
		    if (LastMod->PreserveTopology())
        {
			    LastMod->Execute();
          LocalModification::HeapRetType h_ret=LastMod->UpdateHeap();
          LocalModification::HeapRetType::iterator hi;
          for (hi=h_ret.begin();hi<h_ret.end();hi++)
            AddOperation((*hi).first,(*hi).second);
		    }
        delete (LastMod);     
      }	 
    }
  }

  void DoDecimation(int Time)
  {
    clock_t time0=clock();
    clock_t difftime=0;
    while ((difftime<(clock_t)Time)&&(HaveOperations()))
    {
      
      DecimationStep(1);
      clock_t time1=clock();
      difftime=time1-time0;
    }
  }

///initialize for all vertex the temporary mark must call only at the start of decimation
  void InitDecimation()
  {
    tm.InitIMark();
  }

///erase from the heap the operation that are not upto date
  void CleanHeap()
  {
    vector<vcg::tetra::HeapElem>::iterator hi;
	  for(hi=h.begin();hi!=h.end();++hi)
		if(!(*hi)->LM->IsUpToDate())
		{
			*hi=h.back();
			h.pop_back();
			if(hi==h.end()) break;
		}
		//printf("\nReduced heap from %i to %i",sz,h.size());
		make_heap(h.begin(),h.end());
  }

  ///return true if the decimation have some operation to perform
  bool HaveOperations()
  {
    return (h.size()>0);
  }
    
};//end class decimation

}//end namespace
}//end namespace
#endif