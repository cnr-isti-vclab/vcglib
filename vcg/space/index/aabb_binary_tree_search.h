#ifndef __VCGLIB_AABBBINARYTREESEARCH
#define __VCGLIB_AABBBINARYTREESEARCH

// stl headers
#include <limits>
#include <vector>

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
			typename SCALARTYPE minDist;
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

		inline void Clear(void);

		template <class OBJITERATOR, class OBJITERATORPTRFUNCT, class OBJBOXFUNCT, class OBJBARYCENTERFUNCT>
		inline bool Set(const OBJITERATOR & oBegin, const OBJITERATOR & oEnd, const unsigned int size, const unsigned int maxElemsPerLeaf, const ScalarType & leafBoxMaxVolume, const bool useVariance, OBJITERATORPTRFUNCT & objPtr, OBJBOXFUNCT & objBox, OBJBARYCENTERFUNCT & objBarycenter);

		template <class OBJPOINTDISTANCEFUNCT>
		ObjPtr GetClosest(OBJPOINTDISTANCEFUNCT & getPointDistance, const CoordType & p, ScalarType & minDist, CoordType & res) const;

		inline TreeType & Tree(void);
		inline const TreeType & Tree(void) const;

	protected:
		TreeType tree;

		static inline CoordType Abs(const CoordType & p);
		static inline CoordType LowerClamp(const CoordType & p, const ScalarType & r);
		static inline ScalarType MinimumMaxDistance(const ScalarType & currMinMaxDist, const std::vector<typename TreeType::NodeType *> & n, const CoordType & p);

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
void AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::Clear(void) {
	this->tree.Clear();
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
template <class OBJITERATOR, class OBJITERATORPTRFUNCT, class OBJBOXFUNCT, class OBJBARYCENTERFUNCT>
bool AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::Set(const OBJITERATOR & oBegin, const OBJITERATOR & oEnd, const unsigned int size, const unsigned int maxElemsPerLeaf, const ScalarType & leafBoxMaxVolume, const bool useVariance, OBJITERATORPTRFUNCT & objPtr, OBJBOXFUNCT & objBox, OBJBARYCENTERFUNCT & objBarycenter) {
	return (this->tree.Set(oBegin, oEnd, size, maxElemsPerLeaf, leafBoxMaxVolume, useVariance, objPtr, objBox, objBarycenter));
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

	ScalarType minMaxDist = std::numeric_limits<ScalarType>::max();

	candidates->push_back(t.pRoot);

	while (!candidates->empty()) {
		newCandidates->resize(0);
		minMaxDist = ClassType::MinimumMaxDistance(minMaxDist, *candidates, p);
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
	ScalarType closestDist = std::numeric_limits<ScalarType>::max();
	ScalarType closestDistSq = std::numeric_limits<ScalarType>::max();

	for (NodePtrVector_ci ci=leaves.begin(); ci!=leaves.end(); ++ci) {
		if ((*ci)->auxData.searchData.minDist < closestDistSq) {
			for (TreeType::ObjPtrVectorConstIterator si=(*ci)->oBegin; si!=(*ci)->oEnd; ++si) {
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
typename AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::CoordType AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::Abs(const CoordType & p) {
	return (CoordType(math::Abs(p[0]), math::Abs(p[1]), math::Abs(p[2])));
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
typename AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::CoordType AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::LowerClamp(const CoordType & p, const ScalarType & r) {
	return (CoordType(math::Max<ScalarType>(p[0], r), math::Max<ScalarType>(p[1], r), math::Max<ScalarType>(p[2], r)));
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
typename AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::ScalarType AABBBinaryTreeSearch<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::MinimumMaxDistance(const ScalarType & currMinMaxDist, const std::vector<typename TreeType::NodeType *> & n, const CoordType & p) {
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

}	// end namespace vcg

#endif // #ifndef __VCGLIB_AABBBINARYTREESEARCH
