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
****************************************************************************/

#ifndef __VCGLIB_LOCALOPTIMIZATION
#define __VCGLIB_LOCALOPTIMIZATION
#include<vector>
#include<algorithm>
#include<time.h>
#include<math.h>

namespace vcg{
enum ModifierType{	TriEdgeCollapse, TriEdgeSwap, TriVertexSplit,
				TetraEdgeCollapse,TetraEdgeSplit,TetraEdgeSwap};

/*@{*/
/// This is an abstract class to define the interface of a local modification
template <class ScalarType, class HeapType>
class LocalModification
{

 public:
	LocalModification(){};
	~LocalModification(){};
  
	/// return the type of operation
	virtual ModifierType IsOfType() = 0 ;

	/// return true if the data have not changed since it was created
	virtual bool IsUpToDate() = 0 ;

	/// return true if no constraint disallow this operation to be performed (ex: change of topology in edge collapses)
	virtual bool IsFeasible() = 0;

	/// Compute the priority to be used in the heap
	virtual ScalarType ComputePriority()=0;

	/// Compute the error caused by this modification (can be the same as priority)
	virtual ScalarType ComputeError()=0;

	/// Update the heap as a consequence of this operation
	virtual void UpdateHeap(HeapType&)=0;

	/// Perform the operation
	virtual void Execute()=0;
};	//end class local modification



template<class MeshType>
class LocalOptimization
{
public:
	LocalOptimization(){}

	struct  HeapElem;
	// scalar type
	typedef typename MeshType::ScalarType ScalarType;
	// type of the heap
	typedef typename std::vector<HeapElem> HeapType;	
	// modification type	
	typedef  LocalModification <ScalarType, HeapType>  LocModType;
	// modification Pointer type	
	typedef  LocalModification <ScalarType, HeapType> * LocModPtrType;
	


	/// termination conditions	
	enum {	LOnSimplices	= 0x00,	// test number of simplicies	
			LOnVertices		= 0x01, // test number of verticies
			LOnOps			= 0x02, // test number of operations
			LOMetric		= 0x04, // test Metric (error, quality...instance dependent)
			LOTime			= 0x08  // test how much time is passed since the start
		} terminateFlags;
	int nPerfmormedOps,
		nTargetOps,
		nTargetSimplices,
		nTargetVertices;

	float	timeBudget,
			start, 
			currMetric,
			targetMetric;

	void SetTerminationFlag		(int v){tf |= v;}
	void ClearTerminationFlag	(int v){tf &= ~v;}
	bool IsTerminationFlag		(int v){return (tf & v);}

	void SetTargetSimplices	(int ts			){nTargetSimplices	= ts;	SetTerminationFlag(LOnSimplices);	}	 	
	void SetTargetVertices	(int tv			){nTargetSimplices	= tv;	SetTerminationFlag(LOnVertices);	} 
	void SetTargetOperations(int to			){nTargetOps		= to;	SetTerminationFlag(LOnOps);			} 
	void SetTargetMetric	(ScalarType tm	){targetMetric		= tm;	SetTerminationFlag(LOMetric);		} 
	void SetTimeBudget		(float tb		){timeBudget		= tb;	SetTerminationFlag(LOTime);			} 


	/// the mesh to optimize
	MeshType * m;



	///the heap of operations
	HeapType h;

  ///the element of the heap
  struct HeapElem
  {
	  ~HeapElem(){
		  delete locModPtr;
	  }

    ///pointer to instance of local modifier
    LocModPtrType locModPtr;

    ///temporary mark for the opration
    int imark;
   
    HeapElem(const LocModPtrType _locModPtr, int _imark)
    {
		locModPtr = _locModPtr;
		imark = _imark;
    };

    const bool operator <(const HeapElem & h) const 
    { 
		return (locModPtr->Priority()<h.locModPtr->Priority());
	}

    bool IsUpToDate()
    {
		return locModPtr->IsUpToDate(imark);

	}
/*
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
*/

  };



  /// Default Constructor
	LocalOptimization(MeshType &_m):m(_m){};
  /// Default distructor
  ~LocalOptimization(){};

  /// main cycle of optimization
  void DoOptimization()
  {
	int i=0;
	while( !GoalReached())
    {
      if( ! h.back().IsUpToDate())	
		    h.pop_back(); // if it is out of date simply discard it
	    else  
      {	
        h.back().locModPtr->ComputeError();
        LocModPtrType LastMod= h.back().locModPtr;
        h.pop_back();

		// check if it is feasible
		if (LastMod->IsFeasible())
		{
			LastMod->Execute();
			LastMod->UpdateHeap(h);
		 }
	  }
	}
  }
 
	///initialize for all vertex the temporary mark must call only at the start of decimation
	void Init()
	{
		m.InitIMark();
	}




	/// say if the process is to end or not: the process ends when any of the termination conditions is verified
	/// override this function to implemetn other tests
	bool GoalReached(){
		assert ( ( ( terminateFlags & LOnSimplices	)==0) ||  ( nTargetSimplices!= -1));
		assert ( ( ( terminateFlags & LOnVertices	)==0) ||  ( nTargetVertices	!= -1));
		assert ( ( ( terminateFlags & LOnOps		)==0) ||  ( nTargetOps		!= -1));
		assert ( ( ( terminateFlags & LOMetric		)==0) ||  ( targetMetric	!= -1));
		assert ( ( ( terminateFlags & LOTime		)==0) ||  ( timeBudget		!= -1));

		if ( ( terminateFlags & LOnSimplices)	&&	( m->SimplexNumber()< nTargetSimplices)) return true;
		if ( ( terminateFlags & LOnVertices)	&&  ( m->VertexNumber() < nTargetVertices)) return true;
		if ( ( terminateFlags & LOnOps)			&&  ( nPerfmormedOps	== nTargetOps)) return true;
		if ( ( terminateFlags & LOMetric)		&&  ( currMetric		> targetMetric)) return true;
		if ( ( terminateFlags & LOTime)			&&	( (clock()-start)/(float)CLOCKS_PER_SEC > timeBudget)) return true;
	}



///erase from the heap the operations that are out of date
  void ClearHeap()
  {
		typename HeapType::iterator hi;
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

};//end class decimation

}//end namespace
#endif
