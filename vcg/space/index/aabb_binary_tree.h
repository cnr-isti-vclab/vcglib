#ifndef __VCGLIB_AABBBINARYTREE
#define __VCGLIB_AABBBINARYTREE

// stl headers
#include <vector>

// vcg headers
#include <vcg/space/point3.h>
#include <vcg/space/box3.h>


/***************************************************************************************/

namespace vcg {

/*

Class AABBBinaryTree

SAMPLE USAGE:

NOTES:

*/

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
class AABBBinaryTree {
	public:
		typedef AABBBinaryTree<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE> ClassType;
		typedef OBJTYPE ObjType;
		typedef ObjType * ObjPtr;
		typedef SCALARTYPE ScalarType;
		typedef NODEAUXDATATYPE NodeAuxDataType;
		typedef Point3<ScalarType> CoordType;

		typedef std::vector<ObjPtr> ObjPtrVector;
		typedef typename ObjPtrVector::iterator ObjPtrVectorIterator;
		typedef typename ObjPtrVector::const_iterator ObjPtrVectorConstIterator;

	public:
		class AABBBinaryTreeNode {
			public:
				AABBBinaryTreeNode * parent;
				AABBBinaryTreeNode * children[2];
				CoordType boxCenter;
				CoordType boxHalfDims;
				ObjPtrVectorIterator oBegin;
				ObjPtrVectorIterator oEnd;
				unsigned int numObjects;
				NodeAuxDataType auxData;

				inline AABBBinaryTreeNode(AABBBinaryTreeNode * pParent = 0);
				inline ~AABBBinaryTreeNode(void);

				inline void Clear(void);

				inline bool IsLeaf(void) const;
		};

		typedef AABBBinaryTreeNode NodeType;

		ObjPtrVector pObjects;
		NodeType * pRoot;

		inline AABBBinaryTree(void);
		inline ~AABBBinaryTree(void);

		inline void Clear(void);

		template <class OBJITERATOR, class OBJITERATORPTRFUNCT, class OBJBOXFUNCT, class OBJBARYCENTERFUNCT>
		inline bool Set(const OBJITERATOR & oBegin, const OBJITERATOR & oEnd, const unsigned int size, const unsigned int maxElemsPerLeaf, const ScalarType & leafBoxMaxVolume, const bool useVariance, OBJITERATORPTRFUNCT & objPtr, OBJBOXFUNCT & objBox, OBJBARYCENTERFUNCT & objBarycenter);

	protected:
		template <class OBJBOXFUNCT, class OBJBARYCENTERFUNCT>
		inline static NodeType * BoundObjects(NodeType * parent, const ObjPtrVectorIterator & oBegin, const ObjPtrVectorIterator & oEnd, const unsigned int size, const unsigned int maxElemsPerLeaf, const ScalarType & leafBoxMaxVolume, const bool useVariance, OBJBOXFUNCT & getBox, OBJBARYCENTERFUNCT & getBarycenter);

		template <class OBJBARYCENTERFUNCT>
		inline static int BalanceMedian(const ObjPtrVectorIterator & oBegin, const ObjPtrVectorIterator & oEnd, const int size, const int splitAxis, OBJBARYCENTERFUNCT & getBarycenter, ObjPtrVectorIterator & medianIter);
};

/***************************************************************************************/



template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
AABBBinaryTree<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::AABBBinaryTree(void) {
	this->pObjects.clear();
	this->pRoot = 0;
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
AABBBinaryTree<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::~AABBBinaryTree(void) {
	this->pObjects.clear();
	delete this->pRoot;
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
void AABBBinaryTree<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::Clear(void) {
	this->pObjects.clear();
	delete this->pRoot;
	this->pRoot = 0;
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
template <class OBJITERATOR, class OBJITERATORPTRFUNCT, class OBJBOXFUNCT, class OBJBARYCENTERFUNCT>
bool AABBBinaryTree<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::Set(const OBJITERATOR & oBegin, const OBJITERATOR & oEnd, const unsigned int size, const unsigned int maxElemsPerLeaf, const ScalarType & leafBoxMaxVolume, const bool useVariance, OBJITERATORPTRFUNCT & objPtr, OBJBOXFUNCT & objBox, OBJBARYCENTERFUNCT & objBarycenter) {
	this->Clear();

	if ((maxElemsPerLeaf == 0) && (leafBoxMaxVolume <= ((ScalarType)0))) {
		return (false);
	}

	this->pObjects.reserve(size);
	for (OBJITERATOR oi=oBegin; oi!=oEnd; ++oi) {
		this->pObjects.push_back(objPtr(*oi));
	}

	this->pRoot = ClassType::BoundObjects(0, this->pObjects.begin(), this->pObjects.end(), (unsigned int)(this->pObjects.size()), maxElemsPerLeaf, leafBoxMaxVolume, useVariance, objBox, objBarycenter);

	return (this->pRoot != 0);
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
template <class OBJBOXFUNCT, class OBJBARYCENTERFUNCT>
typename AABBBinaryTree<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::NodeType * AABBBinaryTree<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::BoundObjects(NodeType * parent, const ObjPtrVectorIterator & oBegin, const ObjPtrVectorIterator & oEnd, const unsigned int size, const unsigned int maxElemsPerLeaf, const ScalarType & leafBoxMaxVolume, const bool useVariance, OBJBOXFUNCT & getBox, OBJBARYCENTERFUNCT & getBarycenter) {
	if (size <= 0) {
		return (0);
	}

	NodeType * pNode = new NodeType(parent);
	if (pNode == 0) {
		return (0);
	}

	pNode->children[0] = 0;
	pNode->children[1] = 0;

	pNode->oBegin = oBegin;
	pNode->oEnd = oEnd;

	Box3<ScalarType> bbox;
	bbox.SetNull();
	for (ObjPtrVectorConstIterator oi=pNode->oBegin; oi!=pNode->oEnd; ++oi) {
		const Box3<ScalarType> tbox = getBox(*(*oi));
		bbox.Add(tbox);
	}

	pNode->boxCenter = bbox.Center();
	pNode->boxHalfDims = bbox.Dim() / ((ScalarType)2);

	const bool bMaxObjectsReached = (((maxElemsPerLeaf > 0) && (size <= maxElemsPerLeaf)) || (size == 1));
	const bool bMaxVolumeReached = ((leafBoxMaxVolume > ((ScalarType)0)) && (bbox.Volume() <= leafBoxMaxVolume));
	const bool isLeaf = bMaxObjectsReached || bMaxVolumeReached;

	if (isLeaf) {
		return (pNode);
	}

	pNode->numObjects = size;

	CoordType pSplit;

	if (useVariance) {
		CoordType mean((ScalarType)0, (ScalarType)0, (ScalarType)0);
		CoordType variance((ScalarType)0, (ScalarType)0, (ScalarType)0);
		for (ObjPtrVectorIterator oi=oBegin; oi!=oEnd; ++oi) {
			const CoordType bc = getBarycenter(*(*oi));
			mean += bc;
			variance[0] += bc[0] * bc[0];
			variance[1] += bc[1] * bc[1];
			variance[2] += bc[2] * bc[2];
		}
		variance[0] -= (mean[0] * mean[0]) / ((ScalarType)size);
		variance[1] -= (mean[1] * mean[1]) / ((ScalarType)size);
		variance[2] -= (mean[2] * mean[2]) / ((ScalarType)size);
		pSplit = variance;
	}
	else {
		pSplit = pNode->boxHalfDims;
	}

	ScalarType maxDim = pSplit[0];
	int splitAxis = 0;
	if (maxDim < pSplit[1]) {
		maxDim = pSplit[1];
		splitAxis = 1;
	}
	if (maxDim < pSplit[2]) {
		maxDim = pSplit[2];
		splitAxis = 2;
	}

	ObjPtrVectorIterator median;
	const int lSize = ClassType::BalanceMedian(pNode->oBegin, pNode->oEnd, size, splitAxis, getBarycenter, median);
	const int rSize = size - lSize;

	if (lSize > 0) {
		pNode->children[0] = ClassType::BoundObjects(pNode, pNode->oBegin, median, lSize, maxElemsPerLeaf, leafBoxMaxVolume, useVariance, getBox, getBarycenter);
		if (pNode->children[0] == 0) {
			delete pNode;
			return (0);
		}
	}

	if (rSize > 0) {
		pNode->children[1] = ClassType::BoundObjects(pNode, median, pNode->oEnd, rSize, maxElemsPerLeaf, leafBoxMaxVolume, useVariance, getBox, getBarycenter);
		if (pNode->children[1] == 0) {
			delete pNode;
			return (0);
		}
	}

	return (pNode);
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
template <class OBJBARYCENTERFUNCT>
int AABBBinaryTree<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::BalanceMedian(const ObjPtrVectorIterator & oBegin, const ObjPtrVectorIterator & oEnd, const int size, const int splitAxis, OBJBARYCENTERFUNCT & getBarycenter, ObjPtrVectorIterator & medianIter) {
	const int iMedian = (size + 1) / 2;

	ObjPtrVectorIterator l, r, i, j;
	ObjPtr iTmp;
	ScalarType pos;
	ObjPtrVectorIterator median = oBegin + iMedian;

	l = oBegin;
	r = oEnd - 1;

	while (l < r) {
		pos = getBarycenter(*(*r))[splitAxis];

		i = l;
		j = r - 1;

		while (true) {
			while ((getBarycenter(*(*i))[splitAxis] <= pos) && (i < r)) {
				i++;
			}
			while ((getBarycenter(*(*j))[splitAxis] > pos) && (j > l)) {
				j--;
			}
			if (i >= j) {
				break;
			}
			iTmp = (*i);
			(*i) = (*j);
			(*j) = iTmp;
		}

		iTmp = (*i);
		(*i) = (*r);
		(*r) = iTmp;

		if (i >= (median)) {
			r = i - 1;
		}
		if (i <= (median)) {
			l = i + 1;
		}
	}

	medianIter = median;

	return (iMedian);
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
AABBBinaryTree<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::AABBBinaryTreeNode::AABBBinaryTreeNode(NodeType * pParent) {
	this->parent = pParent;
	this->children[0] = 0;
	this->children[1] = 0;
	this->numObjects = 0;
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
AABBBinaryTree<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::AABBBinaryTreeNode::~AABBBinaryTreeNode(void) {
	delete this->children[0];
	delete this->children[1];
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
void AABBBinaryTree<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::AABBBinaryTreeNode::Clear(void) {
	delete this->children[0];
	this->children[0] = 0;

	delete this->children[1];
	this->children[1] = 0;

	this->numObjects = 0;
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
bool AABBBinaryTree<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::AABBBinaryTreeNode::IsLeaf(void) const {
	return ((this->children[0] == 0) && (this->children[1] == 0));
}

}	// end namespace vcg

#endif // #ifndef __VCGLIB_AABBBINARYTREE
