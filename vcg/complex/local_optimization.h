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
  Revision 1.1  2004/07/08 08:25:15  ganovelli
  first draft

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
template <class ScalarType, class HeapType,class MESH_TYPE>
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

	/// Return the priority to be used in the heap (implement static priority)
	virtual ScalarType Priority()=0;

	/// Compute the error caused by this modification (can be the same as priority)
	virtual ScalarType ComputeError()=0;


	/// Perform the operation and return the variation in the number of simplicies (>0 is refinement, <0 is simplification)
	virtual int Execute()=0;

	/// perform initialization
	virtual void Init(MESH_TYPE&m,HeapType&)=0;


	/// Update the heap as a consequence of this operation
	virtual void UpdateHeap(HeapType&)=0;
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
	typedef  LocalModification <ScalarType, HeapType,MeshType>  LocModType;
	// modification Pointer type	
	typedef  LocalModification <ScalarType, HeapType,MeshType> * LocModPtrType;
	


	/// termination conditions	
	 enum {	LOnSimplices	= 0x00,	// test number of simplicies	
			LOnVertices		= 0x01, // test number of verticies
			LOnOps			= 0x02, // test number of operations
			LOMetric		= 0x04, // test Metric (error, quality...instance dependent)
			LOTime			= 0x08  // test how much time is passed since the start
		} ;

	 int tf;
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
		HeapElem(){locModPtr = NULL;}
	  ~HeapElem(){}

    ///pointer to instance of local modifier
    LocModPtrType locModPtr;

   
    HeapElem( LocModPtrType _locModPtr)
    {
		locModPtr = _locModPtr;
    };

    const bool operator <(const HeapElem & h) const 
    { 
		return (locModPtr->Priority()<h.locModPtr->Priority());
	}

    bool IsUpToDate()
    {
			return locModPtr->IsUpToDate();
		}
  };



  /// Default Constructor
	LocalOptimization(MeshType *_m):m(_m){};
  /// Default distructor
  ~LocalOptimization(){};

  /// main cycle of optimization
  void DoOptimization()
  {
		nPerfmormedOps =0;
	int i=0;
	while( !GoalReached())
    {int size = h.size();
			LocModPtrType  locMod   = h.back().locModPtr;
      if( ! h.back().IsUpToDate())	
			{
				h.pop_back(); // if it is out of date simply discard it
			}
	    else  
      {	
        locMod->ComputeError();
        h.pop_back();

				// check if it is feasible
			if (locMod->IsFeasible())
			{
				nPerfmormedOps++;
				int tmp = locMod->Execute();
				m->SimplexNumber()+= tmp;
				locMod->UpdateHeap(h);
				m->VertexNumber()--;
				}

			}
			delete locMod;
		}
  }
 
	///initialize for all vertex the temporary mark must call only at the start of decimation
	///by default it takes the first element in the heap and calls Init (static funcion) of that type
	///of local modification. 
	void Init()
	{
		m->InitIMark();
		if(!h.empty())
		{
			(*h.begin()).locModPtr->Init(*m,h);
		}
	}




	/// say if the process is to end or not: the process ends when any of the termination conditions is verified
	/// override this function to implemetn other tests
	bool GoalReached(){
		assert ( ( ( tf & LOnSimplices	)==0) ||  ( nTargetSimplices!= -1));
		assert ( ( ( tf & LOnVertices	)==0) ||  ( nTargetVertices	!= -1));
		assert ( ( ( tf & LOnOps		)==0) ||  ( nTargetOps		!= -1));
		assert ( ( ( tf & LOMetric		)==0) ||  ( targetMetric	!= -1));
		assert ( ( ( tf & LOTime		)==0) ||  ( timeBudget		!= -1));

		if(h.empty()) return true;

		if ( ( tf & LOnSimplices)	&&	( m->SimplexNumber()< nTargetSimplices)) return true;
		if ( ( tf & LOnVertices)	&&  ( m->VertexNumber() < nTargetVertices)) return true;
		if ( ( tf & LOnOps)			&&  ( nPerfmormedOps	== nTargetOps)) return true;
		if ( ( tf & LOMetric)		&&  ( currMetric		> targetMetric)) return true;
		if ( ( tf & LOTime)			&&	( (clock()-start)/(float)CLOCKS_PER_SEC > timeBudget)) return true;
		return false;
	}



///erase from the heap the operations that are out of date
  void ClearHeap()
  {
		typename HeapType::iterator hi;
		for(hi=h.begin();hi!=h.end();++hi)
			if(!(*hi).locModPtr->IsUpToDate())
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
