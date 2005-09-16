/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005                                                \/)\/    *
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
Revision 1.3  2005/09/16 10:04:15  m_di_benedetto
Modified interface, added GetKClosest().

Revision 1.2  2005/09/11 11:46:21  m_di_benedetto
First Commit




****************************************************************************/

#ifndef __VCGLIB_AABBBINARYTREESEARCH
#define __VCGLIB_AABBBINARYTREESEARCH

// stl headers
#include <limits>
#include <vector>
#include <queue>
#include <deque>

// vcg headers
#include <vcg/space/index/aabb_binary_tree.h>


/***************************************************************************************/

namespace vcg {

/*

Class AABBBinaryTreeSearch

SAMPLE USAGE:

NOTES:

*/

template <class OBJTYPE, class SCALARTYPE, class NODEUSERATATYPE>
class AABBBinaryTreeSearch {
	public:
		struct NodeSearchDataType {
			SCALARTYPE minDist;
		};
		struct NodeAuxDataType {
			NodeSearchDataType searchData;
			NODEUSERATATYPE userData;
		};
		typedef AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEUSERATATYPE> ClassType;
		typedef OBJTYPE ObjType;
		typedef SCALARTYPE ScalarType;
		typedef NODEUSERATATYPE NodeUserDataType;
		typedef AABBBinaryTree<ObjType, ScalarType, NodeAuxDataType> TreeType;
		typedef typename TreeType::ObjPtr ObjPtr;
		typedef typename TreeType::CoordType CoordType;

		inline AABBBinaryTreeSearch(void);
		inline ~AABBBinaryTreeSearch(void);

		inline TreeType & Tree(void);
		inline const TreeType & Tree(void) const;

		template <class OBJITERATOR, class OBJITERATORPTRFUNCT, class OBJBOXFUNCT, class OBJBARYCENTERFUNCT>
		inline bool Set(const OBJITERATOR & oBegin, const OBJITERATOR & oEnd, OBJITERATORPTRFUNCT & objPtr, OBJBOXFUNCT & objBox, OBJBARYCENTERFUNCT & objBarycenter, const unsigned int maxElemsPerLeaf = 1, const ScalarType & leafBoxMaxVolume = ((ScalarType)0), const bool useVariance = true);

		inline void Clear(void);

		template <class OBJPOINTDISTANCEFUNCT>
		ObjPtr GetClosest(OBJPOINTDISTANCEFUNCT & getPointDistance, const CoordType & p, ScalarType & minDist, CoordType & res) const;

		template <class OBJPOINTDISTANCEFUNCT, class OBJPTRCONTAINER, class DISTCONTAINER, class POINTCONTAINER>
		unsigned int GetKClosests(OBJPOINTDISTANCEFUNCT & getPointDistance, const unsigned int k, const CoordType & p, OBJPTRCONTAINER & closestObjs, DISTCONTAINER & distances, POINTCONTAINER & closestPts) const;

	protected:
		struct ClosestObjType {
			ObjPtr pObj;
			ScalarType minDist;
			CoordType closestPt;
		};

		class CompareClosest {
		public:
			bool operator () (const ClosestObjType & a, const ClosestObjType & b) {
				return (a.minDist < b.minDist);
			}
		};

		typedef std::priority_queue<typename ClassType::ClosestObjType, std::deque<typename ClassType::ClosestObjType>, typename ClassType::CompareClosest> PQueueType;

		TreeType tree;

		static inline CoordType Abs(const CoordType & p);
		static inline CoordType LowerClamp(const CoordType & p, const ScalarType & r);
		static inline ScalarType MinimumMaxSquareDistance(const ScalarType & currMinMaxDist, const std::vector<typename TreeType::NodeType *> & n, const CoordType & p);
		static inline void MinMaxSquareDistance(const typename TreeType::NodeType * n, const CoordType & p, ScalarType & dmin, ScalarType & dmax);

		template <class OBJPOINTDISTANCEFUNCT>
		static inline void DepthFirstCollect(const CoordType & p, typename TreeType::NodeType * node, ScalarType & mindmax, const unsigned int k, PQueueType & pq, OBJPOINTDISTANCEFUNCT & getPointDistance);

};

/***************************************************************************************/



template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::AABBBinaryTreeSearch(void) {

}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::~AABBBinaryTreeSearch(void) {
	this->Clear();
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
typename AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::TreeType & AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::Tree(void) {
	return (this->tree);
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
const typename AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::TreeType & AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::Tree(void) const {
	return (this->tree);
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
template <class OBJITERATOR, class OBJITERATORPTRFUNCT, class OBJBOXFUNCT, class OBJBARYCENTERFUNCT>
bool AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::Set(const OBJITERATOR & oBegin, const OBJITERATOR & oEnd, OBJITERATORPTRFUNCT & objPtr, OBJBOXFUNCT & objBox, OBJBARYCENTERFUNCT & objBarycenter, const unsigned int maxElemsPerLeaf, const ScalarType & leafBoxMaxVolume, const bool useVariance) {
	return (this->tree.Set(oBegin, oEnd, objPtr, objBox, objBarycenter, maxElemsPerLeaf, leafBoxMaxVolume, useVariance));
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
void AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::Clear(void) {
	this->tree.Clear();
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
template <class OBJPOINTDISTANCEFUNCT>
typename AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::ObjPtr AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::GetClosest(OBJPOINTDISTANCEFUNCT & getPointDistance, const CoordType & p, ScalarType & minDist, CoordType & res) const {
	typedef std::vector<typename TreeType::NodeType *> NodePtrVector;
	typedef typename NodePtrVector::const_iterator NodePtrVector_ci;

	const TreeType & t = this->tree;
	TreeType::NodeType * pRoot = t.pRoot;

	if (pRoot == 0) {
		return (0);
	}

	NodePtrVector clist1;
	NodePtrVector clist2;
	NodePtrVector leaves;

	NodePtrVector * candidates = &clist1;
	NodePtrVector * newCandidates = &clist2;

	clist1.reserve(t.pObjects.size());
	clist2.reserve(t.pObjects.size());
	leaves.reserve(t.pObjects.size());

	clist1.resize(0);
	clist2.resize(0);
	leaves.resize(0);

#ifdef max
#undef max
#endif
	ScalarType minMaxDist = std::numeric_limits<ScalarType>::max();

	candidates->push_back(t.pRoot);

	while (!candidates->empty()) {
		newCandidates->resize(0);
		minMaxDist = ClassType::MinimumMaxSquareDistance(minMaxDist, *candidates, p);
		for (NodePtrVector_ci ci=candidates->begin(); ci!=candidates->end(); ++ci) {
			if ((*ci)->auxData.searchData.minDist < minMaxDist) {
				if ((*ci)->IsLeaf()) {
					leaves.push_back(*ci);
				}
				else {
					if ((*ci)->children[0] != 0) {
						newCandidates->push_back((*ci)->children[0]);
					}
					if ((*ci)->children[1] != 0) {
						newCandidates->push_back((*ci)->children[1]);
					}
				}
			}
		}
		NodePtrVector * cSwap = candidates;
		candidates = newCandidates;
		newCandidates = cSwap;
	}

	clist1.clear();
	clist2.clear();

	ObjPtr closestObject = 0;
	CoordType closestPoint;
	ScalarType closestDist = math::Sqrt(minMaxDist) + std::numeric_limits<ScalarType>::epsilon();
	ScalarType closestDistSq = minMaxDist + std::numeric_limits<ScalarType>::epsilon();

	for (NodePtrVector_ci ci=leaves.begin(); ci!=leaves.end(); ++ci) {
		if ((*ci)->auxData.searchData.minDist < closestDistSq) {
			for (typename TreeType::ObjPtrVectorConstIterator si=(*ci)->oBegin; si!=(*ci)->oEnd; ++si) {
				if (getPointDistance(*(*si), p, closestDist, closestPoint)) {
					closestDistSq = closestDist * closestDist;
					closestObject = (*si);
				}
			}
		}
	}

	leaves.clear();

	res = closestPoint;
	minDist = closestDist;

	return (closestObject);
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
template <class OBJPOINTDISTANCEFUNCT>
void AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::DepthFirstCollect(const CoordType & p, typename TreeType::NodeType * node, ScalarType & mindmax, const unsigned int k, PQueueType & pq, OBJPOINTDISTANCEFUNCT & getPointDistance) {
	const CoordType dc = ClassType::Abs(p - node->boxCenter);

	if (pq.size() >= k) {
		const ScalarType dmin = ClassType::LowerClamp(dc - node->boxHalfDims, (ScalarType)0).SquaredNorm();
		if (dmin >= mindmax) {
			return;
		}
	}

	if (node->IsLeaf()) {
		bool someInserted = true;
		for (typename TreeType::ObjPtrVectorConstIterator si=node->oBegin; si!=node->oEnd; ++si) {
			ScalarType minDst = (pq.size() >= k) ? (pq.top().minDist) : (std::numeric_limits<ScalarType>::max());
			ClosestObjType cobj;
			if (getPointDistance(*(*si), p, minDst, cobj.closestPt)) {
				someInserted = true;
				cobj.pObj = (*si);
				cobj.minDist = minDst;
				if (pq.size() >= k) {
					pq.pop();
				}
				pq.push(cobj);
			}
		}
		if (someInserted) {
			if (pq.size() >= k) {
				const ScalarType dmax = pq.top().minDist;
				const ScalarType sqdmax = dmax * dmax;
				if (sqdmax < mindmax) {
					mindmax = sqdmax;
				}
			}
		}
	}
	else {
		if (node->children[0] != 0) {
			DepthFirstCollect(p, node->children[0], mindmax, k, pq, getPointDistance);
		}
		if (node->children[1] != 0) {
			DepthFirstCollect(p, node->children[1], mindmax, k, pq, getPointDistance);
		}
	}
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
template <class OBJPOINTDISTANCEFUNCT, class OBJPTRCONTAINER, class DISTCONTAINER, class POINTCONTAINER>
unsigned int AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::GetKClosests(OBJPOINTDISTANCEFUNCT & getPointDistance, const unsigned int k, const CoordType & p, OBJPTRCONTAINER & closestObjs, DISTCONTAINER & distances, POINTCONTAINER & closestPts) const {
	typedef std::vector<typename TreeType::NodeType *> NodePtrVector;
	typedef typename NodePtrVector::const_iterator NodePtrVector_ci;

	const TreeType & t = this->tree;
	TreeType::NodeType * pRoot = t.pRoot;

	if (pRoot == 0) {
		return (0);
	}

	PQueueType pq;
	ScalarType mindmax = std::numeric_limits<ScalarType>::max();

	ClassType::DepthFirstCollect(p, pRoot, mindmax, k, pq, getPointDistance);

	const unsigned int sz = pq.size();

	while (!pq.empty()) {
		ClosestObjType cobj = pq.top();
		pq.pop();
		closestObjs.push_back(cobj.pObj);
		distances.push_back(cobj.minDist);
		closestPts.push_back(cobj.closestPt);
	}

	return (sz);
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
typename AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::CoordType AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::Abs(const CoordType & p) {
	return (CoordType(math::Abs(p[0]), math::Abs(p[1]), math::Abs(p[2])));
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
typename AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::CoordType AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::LowerClamp(const CoordType & p, const ScalarType & r) {
	return (CoordType(math::Max<ScalarType>(p[0], r), math::Max<ScalarType>(p[1], r), math::Max<ScalarType>(p[2], r)));
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
typename AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::ScalarType AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::MinimumMaxSquareDistance(const ScalarType & currMinMaxDist, const std::vector<typename TreeType::NodeType *> & n, const CoordType & p) {
	typedef std::vector<typename TreeType::NodeType *> NodePtrVector;
	typedef typename NodePtrVector::const_iterator NodePtrVector_ci;

	ScalarType minMaxDist = currMinMaxDist;

	for (NodePtrVector_ci bv=n.begin(); bv!=n.end(); ++bv) {
		const CoordType dc = ClassType::Abs(p - (*bv)->boxCenter);
		const ScalarType maxDist = (dc + (*bv)->boxHalfDims).SquaredNorm();
		(*bv)->auxData.searchData.minDist = ClassType::LowerClamp(dc - (*bv)->boxHalfDims, (ScalarType)0).SquaredNorm();
		if (maxDist < minMaxDist) {
			minMaxDist = maxDist;
		}
	}

	return (minMaxDist);
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
void AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::MinMaxSquareDistance(const typename TreeType::NodeType * n, const CoordType & p, ScalarType & dmin, ScalarType & dmax) {
	const CoordType dc = ClassType::Abs(p - n->boxCenter);
	dmax = (dc + n->boxHalfDims).SquaredNorm();
	dmin = ClassType::LowerClamp(dc - n->boxHalfDims, (ScalarType)0).SquaredNorm();
}

}	// end namespace vcg

#endif // #ifndef __VCGLIB_AABBBINARYTREESEARCH
