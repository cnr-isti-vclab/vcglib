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
  Revision 1.7  2004/09/29 17:08:39  ganovelli
  changed > to < in heapelem comparison

  Revision 1.6  2004/09/28 09:57:08  cignoni
  Better Doxygen docs

  Revision 1.5  2004/09/15 10:40:20  ponchio
  typedef LocalOptimization HeapType -> public:

  Revision 1.4  2004/09/08 15:10:59  ganovelli
  *** empty log message ***

  Revision 1.3  2004/07/27 09:46:15  cignoni
  First working version of the LocalOptimization/Simplification Framework

  Revision 1.1  2004/07/15 12:04:14  ganovelli
  minor changes

  Revision 1.2  2004/07/09 10:22:56  ganovelli
  working draft

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

template<class MeshType> 
class LocalOptimization;

enum ModifierType{	TetraEdgeCollapseOp, TriEdgeSwapOp, TriVertexSplitOp,
				TriEdgeCollapseOp,TetraEdgeSpliOpt,TetraEdgeSwapOp};
/** \addtogroup tetramesh */
/*@{*/
/// This abstract class define which functions  a local modification to be used in the LocalOptimization.
template <class MeshType>
class LocalModification
{
 public:
        typedef typename LocalOptimization<MeshType>::HeapType HeapType;
        typedef typename MeshType::ScalarType ScalarType;


	inline LocalModification(){};
	virtual ~LocalModification(){};
  
	/// return the type of operation
	virtual ModifierType IsOfType() = 0 ;

	/// return true if the data have not changed since it was created
	virtual bool IsUpToDate() = 0 ;

	/// return true if no constraint disallow this operation to be performed (ex: change of topology in edge collapses)
	virtual bool IsFeasible() = 0;

	/// Compute the priority to be used in the heap
	virtual ScalarType ComputePriority()=0;

	/// Return the priority to be used in the heap (implement static priority)
	virtual ScalarType Priority() const =0;

	/// Perform the operation and return the variation in the number of simplicies (>0 is refinement, <0 is simplification)
	virtual void Execute(MeshType &m)=0;

	/// perform initialization
	static void Init(MeshType &m, HeapType&);

  virtual const char *Info(MeshType &) {return 0;}
	/// Update the heap as a consequence of this operation
	virtual void UpdateHeap(HeapType&)=0;
};	//end class local modification


/// LocalOptimization:
/// This class implements the algorihms running on 0-1-2-3-simplicial complex that are based on local modification
/// The local modification can be and edge_collpase, or an edge_swap, a vertex plit...as far as they implement
/// the interface defined in LocalModification.
/// Implementation note: in order to keep the local modification itself indepented by its use in this class, they are not
/// really derived by LocalModification. Instead, a wrapper is done to this purpose (see vcg/complex/tetramesh/decimation/collapse.h)

template<class MeshType>
class LocalOptimization
{
public:
  LocalOptimization(MeshType &mm): m(mm){ ClearTermination();e=0.0;}

	struct  HeapElem;
	// scalar type
	typedef typename MeshType::ScalarType ScalarType;
	// type of the heap
	typedef typename std::vector<HeapElem> HeapType;	
	// modification type	
	typedef  LocalModification <MeshType>  LocModType;
	// modification Pointer type	
	typedef  LocalModification <MeshType> * LocModPtrType;
	


	/// termination conditions	
	 enum LOTermination {	
      LOnSimplices	= 0x01,	// test number of simplicies	
			LOnVertices		= 0x02, // test number of verticies
			LOnOps			= 0x04, // test number of operations
			LOMetric		= 0x08, // test Metric (error, quality...instance dependent)
			LOTime			= 0x10  // test how much time is passed since the start
		} ;

	int tf;
	
  int nPerfmormedOps,
		nTargetOps,
		nTargetSimplices,
		nTargetVertices;

	int	timeBudget,
			start;
	float
			currMetric,
			targetMetric;

	void SetTerminationFlag		(int v){tf |= v;}
	void ClearTerminationFlag	(int v){tf &= ~v;}
	bool IsTerminationFlag		(int v){return ((tf & v)!=0);}

	void SetTargetSimplices	(int ts			){nTargetSimplices	= ts;	SetTerminationFlag(LOnSimplices);	}	 	
	void SetTargetVertices	(int tv			){nTargetVertices	= tv;	SetTerminationFlag(LOnVertices);	} 
	void SetTargetOperations(int to			){nTargetOps		= to;	SetTerminationFlag(LOnOps);			} 

	void SetTargetMetric	(ScalarType tm	){targetMetric		= tm;	SetTerminationFlag(LOMetric);		} 
	void SetTimeBudget		(float tb		){timeBudget		= tb;	SetTerminationFlag(LOTime);			} 

  void ClearTermination()
  {
    tf=0;
    nTargetSimplices=0;
    nTargetOps=0;
    targetMetric=0;
    timeBudget=0;
    nTargetVertices=0;
  }
	/// the mesh to optimize
	MeshType & m;



	///the heap of operations
	HeapType h;

  ///the element of the heap
  // it is just a wrapper of the pointer to the localMod. 
  // std heap does not work for
  // pointers and we want pointers to have heterogenous heaps. 

  struct HeapElem
  {
		inline HeapElem(){locModPtr = NULL;}
	  ~HeapElem(){}

    ///pointer to instance of local modifier
    LocModPtrType locModPtr;

   
    inline HeapElem( LocModPtrType _locModPtr)
    {
		locModPtr = _locModPtr;
    };

    /// STL heap has the largest element as the first one.
    /// usually we mean priority as an error so we should invert the comparison
    const bool operator <(const HeapElem & h) const 
    { 
		  return (locModPtr->Priority() < h.locModPtr->Priority());
	  }

    bool IsUpToDate()
    {
			return locModPtr->IsUpToDate();
		}
  };



  /// Default distructor
  ~LocalOptimization(){};
	
	double e;

  /// main cycle of optimization
  bool DoOptimization()
  {
    start=clock();
		nPerfmormedOps =0;
#ifdef __SAVE__LOG__
 		FILE * fo=fopen("log.txt","w");
#endif __SAVE__LOG__    
		while( !GoalReached() && !h.empty())
			{
				std::pop_heap(h.begin(),h.end());
        LocModPtrType  locMod   = h.back().locModPtr;
				h.pop_back();
        				
				if( locMod->IsUpToDate() )	
				{	
          //printf("popped out: %s\n",locMod->Info(m));
         	// check if it is feasible
					if (locMod->IsFeasible())
					{
#ifdef __SAVE__LOG__
						fprintf(fo,"%s",locMod->Info(m));
#endif __SAVE__LOG__
	
						nPerfmormedOps++;
						locMod->Execute(m);
						locMod->UpdateHeap(h);
						}
				}
        //else printf("popped out unfeasible\n");
				delete locMod;
			}
#ifdef __SAVE__LOG__
			fclose(fo);
#endif __SAVE__LOG__
		return !(h.empty());
  }
 
	///initialize for all vertex the temporary mark must call only at the start of decimation
	///by default it takes the first element in the heap and calls Init (static funcion) of that type
	///of local modification. 
	template <class LocalModificationType> void Init()
	{
		m.InitVertexIMark();
		LocalModificationType::Init(m,h);
		std::make_heap(h.begin(),h.end());
	}




	/// say if the process is to end or not: the process ends when any of the termination conditions is verified
	/// override this function to implemetn other tests
	bool GoalReached(){
		assert ( ( ( tf & LOnSimplices	)==0) ||  ( nTargetSimplices!= -1));
		assert ( ( ( tf & LOnVertices	)==0) ||  ( nTargetVertices	!= -1));
		assert ( ( ( tf & LOnOps		)==0) ||  ( nTargetOps		!= -1));
		assert ( ( ( tf & LOMetric		)==0) ||  ( targetMetric	!= -1));
		assert ( ( ( tf & LOTime		)==0) ||  ( timeBudget		!= -1));

		if ( IsTerminationFlag(LOnSimplices) &&	( m.SimplexNumber()<= nTargetSimplices)) return true;
		if ( IsTerminationFlag(LOnVertices)  &&  ( m.VertexNumber() <= nTargetVertices)) return true;
		if ( IsTerminationFlag(LOnOps)		   && (nPerfmormedOps	== nTargetOps)) return true;
		if ( IsTerminationFlag(LOMetric)		 &&  ( currMetric		> targetMetric)) return true;
		if ( IsTerminationFlag(LOTime)			 &&	( (clock()-start)/(float)CLOCKS_PER_SEC > timeBudget)) return true;
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
