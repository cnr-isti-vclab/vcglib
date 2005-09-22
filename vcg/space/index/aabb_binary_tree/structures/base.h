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


****************************************************************************/

#ifndef TREE_TYPE
#pragma error message("\nThis file should never be directly included.\n")
#else

// standard headers
#include <assert.h>

// stl headers
#include <vector>

// vcg headers
#include <vcg/space/point3.h>
#include <vcg/space/box3.h>

/***************************************************************************************/

namespace vcg {

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
class TREE_TYPE {
	public:
		typedef TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE> ClassType;
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
			protected:
#ifdef __VCGLIB_AABBBINARYTREE_SA
				unsigned int splitAxis;
#endif
#ifdef __VCGLIB_AABBBINARYTREE_FL
				unsigned int flags;
#endif
#ifdef __VCGLIB_AABBBINARYTREE_IV
				ScalarType integralValue;
#endif
#ifdef __VCGLIB_AABBBINARYTREE_SV
				ScalarType scalarValue;
#endif
			public:
				AABBBinaryTreeNode * parent;
				AABBBinaryTreeNode * children[2];
				CoordType boxCenter;
				CoordType boxHalfDims;
				ObjPtrVectorIterator oBegin;
				ObjPtrVectorIterator oEnd;
				NodeAuxDataType auxData;

				inline AABBBinaryTreeNode(AABBBinaryTreeNode * pParent = 0);
				inline ~AABBBinaryTreeNode(void);

				inline void Clear(void);
				inline bool IsLeaf(void) const;
				inline unsigned int ObjectsCount(void) const;

				static inline bool HasSplitAxis(void);
				static inline bool HasFlags(void);
				static inline bool HasIntegralValue(void);
				static inline bool HasScalarValue(void);

				inline unsigned int & SplitAxis(void);
				inline const unsigned int & SplitAxis(void) const;

				inline unsigned int & Flags(void);
				inline const unsigned int & Flags(void) const;

				inline int & IntegralValue(void);
				inline const int & IntegralValue(void) const;

				inline ScalarType & ScalarValue(void);
				inline const ScalarType & ScalarValue(void) const;
		};

		typedef AABBBinaryTreeNode NodeType;

		ObjPtrVector pObjects;
		NodeType * pRoot;

		inline TREE_TYPE(void);
		inline ~TREE_TYPE(void);

		inline void Clear(void);

		template <class OBJITERATOR, class OBJITERATORPTRFUNCT, class OBJBOXFUNCT, class OBJBARYCENTERFUNCT>
		inline bool Set(const OBJITERATOR & oBegin, const OBJITERATOR & oEnd, OBJITERATORPTRFUNCT & objPtr, OBJBOXFUNCT & objBox, OBJBARYCENTERFUNCT & objBarycenter, const unsigned int maxElemsPerLeaf = 1, const ScalarType & leafBoxMaxVolume = ((ScalarType)0), const bool useVariance = true);

	protected:
		template <class OBJBOXFUNCT, class OBJBARYCENTERFUNCT>
		inline static NodeType * BoundObjects(NodeType * parent, const ObjPtrVectorIterator & oBegin, const ObjPtrVectorIterator & oEnd, const unsigned int size, const unsigned int maxElemsPerLeaf, const ScalarType & leafBoxMaxVolume, const bool useVariance, OBJBOXFUNCT & getBox, OBJBARYCENTERFUNCT & getBarycenter);

		template <class OBJBARYCENTERFUNCT>
		inline static int BalanceMedian(const ObjPtrVectorIterator & oBegin, const ObjPtrVectorIterator & oEnd, const int size, const int splitAxis, OBJBARYCENTERFUNCT & getBarycenter, ObjPtrVectorIterator & medianIter);
};



template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::TREE_TYPE(void) {
	this->pObjects.clear();
	this->pRoot = 0;
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::~TREE_TYPE(void) {
	this->pObjects.clear();
	delete this->pRoot;
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
void TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::Clear(void) {
	this->pObjects.clear();
	delete this->pRoot;
	this->pRoot = 0;
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
template <class OBJITERATOR, class OBJITERATORPTRFUNCT, class OBJBOXFUNCT, class OBJBARYCENTERFUNCT>
bool TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::Set(const OBJITERATOR & oBegin, const OBJITERATOR & oEnd, OBJITERATORPTRFUNCT & objPtr, OBJBOXFUNCT & objBox, OBJBARYCENTERFUNCT & objBarycenter, const unsigned int maxElemsPerLeaf, const ScalarType & leafBoxMaxVolume, const bool useVariance) {
	this->Clear();

	if ((maxElemsPerLeaf == 0) && (leafBoxMaxVolume <= ((ScalarType)0))) {
		return (false);
	}

	const unsigned int size = (unsigned int)std::distance(oBegin, oEnd);

	this->pObjects.reserve(size);
	for (OBJITERATOR oi=oBegin; oi!=oEnd; ++oi) {
		this->pObjects.push_back(objPtr(*oi));
	}

	this->pRoot = ClassType::BoundObjects(0, this->pObjects.begin(), this->pObjects.end(), size, maxElemsPerLeaf, leafBoxMaxVolume, useVariance, objBox, objBarycenter);

	return (this->pRoot != 0);
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
template <class OBJBOXFUNCT, class OBJBARYCENTERFUNCT>
typename TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::NodeType * TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::BoundObjects(NodeType * parent, const ObjPtrVectorIterator & oBegin, const ObjPtrVectorIterator & oEnd, const unsigned int size, const unsigned int maxElemsPerLeaf, const ScalarType & leafBoxMaxVolume, const bool useVariance, OBJBOXFUNCT & getBox, OBJBARYCENTERFUNCT & getBarycenter) {
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
		Box3<ScalarType> tbox;
		getBox(*(*oi), tbox);
		bbox.Add(tbox);
	}

	pNode->boxCenter = bbox.Center();
	pNode->boxHalfDims = bbox.Dim() / ((ScalarType)2);

	const bool bMaxObjectsReached = (((maxElemsPerLeaf > 0) && (size <= maxElemsPerLeaf)) || (size == 1));
	const bool bMaxVolumeReached = ((leafBoxMaxVolume > ((ScalarType)0)) && (bbox.Volume() <= leafBoxMaxVolume));
	const bool isLeaf = bMaxObjectsReached || bMaxVolumeReached;

	if (isLeaf) {
		if (NodeType::HasSplitAxis()) {
			pNode->SplitAxis() = 0;
		}
		return (pNode);
	}

	CoordType pSplit;

	if (useVariance) {
		CoordType mean((ScalarType)0, (ScalarType)0, (ScalarType)0);
		CoordType variance((ScalarType)0, (ScalarType)0, (ScalarType)0);
		for (ObjPtrVectorIterator oi=oBegin; oi!=oEnd; ++oi) {
			CoordType bc;
			getBarycenter(*(*oi), bc);
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

	if (NodeType::HasSplitAxis()) {
		pNode->SplitAxis() = splitAxis;
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
int TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::BalanceMedian(const ObjPtrVectorIterator & oBegin, const ObjPtrVectorIterator & oEnd, const int size, const int splitAxis, OBJBARYCENTERFUNCT & getBarycenter, ObjPtrVectorIterator & medianIter) {
	const int iMedian = (size + 1) / 2;

	ObjPtrVectorIterator l, r, i, j;
	ObjPtr iTmp;
	ScalarType pos;
	ObjPtrVectorIterator median = oBegin + iMedian;
	CoordType bc;

	l = oBegin;
	r = oEnd - 1;

	while (l < r) {
		getBarycenter(*(*r), bc);
		pos = bc[splitAxis];

		i = l;
		j = r - 1;

		while (true) {
			getBarycenter(*(*i), bc);
			while ((bc[splitAxis] <= pos) && (i < r)) {
				i++;
				getBarycenter(*(*i), bc);
			}
			getBarycenter(*(*j), bc);
			while ((bc[splitAxis] > pos) && (j > l)) {
				j--;
				getBarycenter(*(*j), bc);
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
TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::AABBBinaryTreeNode::AABBBinaryTreeNode(NodeType * pParent) {
	this->parent = pParent;
	this->children[0] = 0;
	this->children[1] = 0;
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::AABBBinaryTreeNode::~AABBBinaryTreeNode(void) {
	delete this->children[0];
	delete this->children[1];
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
void TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::AABBBinaryTreeNode::Clear(void) {
	delete this->children[0];
	this->children[0] = 0;

	delete this->children[1];
	this->children[1] = 0;
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
bool TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::AABBBinaryTreeNode::IsLeaf(void) const {
	return ((this->children[0] == 0) && (this->children[1] == 0));
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
unsigned int TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::AABBBinaryTreeNode::ObjectsCount(void) const {
	return ((unsigned int)(std::distance(this->oBegin, this->oEnd)));
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
bool TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::AABBBinaryTreeNode::HasSplitAxis(void) {
#ifdef __VCGLIB_AABBBINARYTREE_SA
	return (true);
#else
	return (false);
#endif
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
bool TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::AABBBinaryTreeNode::HasFlags(void) {
#ifdef __VCGLIB_AABBBINARYTREE_FL
	return (true);
#else
	return (false);
#endif
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
bool TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::AABBBinaryTreeNode::HasScalarValue(void) {
#ifdef __VCGLIB_AABBBINARYTREE_SV
	return (true);
#else
	return (false);
#endif
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
unsigned int & TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::AABBBinaryTreeNode::SplitAxis(void) {
#ifdef __VCGLIB_AABBBINARYTREE_FL
	return (this->splitAxis);
#else
	assert(0);
	return (*(unsigned int *)(this->children[0]));
#endif
}
	
template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
const unsigned int & TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::AABBBinaryTreeNode::SplitAxis(void) const {
#ifdef __VCGLIB_AABBBINARYTREE_FL
	return (this->splitAxis);
#else
	assert(0);
	return (*(const unsigned int *)(this->children[0]));
#endif
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
unsigned int & TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::AABBBinaryTreeNode::Flags(void) {
#ifdef __VCGLIB_AABBBINARYTREE_FL
	return (this->flags);
#else
	assert(0);
	return (*(const unsigned int *)(this->children[0]));
#endif
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
const unsigned int & TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::AABBBinaryTreeNode::Flags(void) const {
#ifdef __VCGLIB_AABBBINARYTREE_FL
	return (this->flags);
#else
	assert(0);
	return (*(const unsigned int *)(this->children[0]));
#endif
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
int & TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::AABBBinaryTreeNode::IntegralValue(void) {
#ifdef __VCGLIB_AABBBINARYTREE_IV
	return (this->integralValue);
#else
	assert(0);
	return ((unsigned int)(this->children[0]));
#endif
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
const int & TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::AABBBinaryTreeNode::IntegralValue(void) const {
#ifdef __VCGLIB_AABBBINARYTREE_IV
	return (this->integralValue);
#else
	assert(0);
	return ((unsigned int)(this->children[0]));
#endif
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
typename TREE_TYPE <OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::ScalarType & TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::AABBBinaryTreeNode::ScalarValue(void) {
#ifdef __VCGLIB_AABBBINARYTREE_SV
	return (this->scalarValue);
#else
	assert(0);
	return (this->boxCenter[0]);
#endif
}

template <class OBJTYPE, class SCALARTYPE, class NODEAUXDATATYPE>
const typename TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::ScalarType & TREE_TYPE<OBJTYPE, SCALARTYPE, NODEAUXDATATYPE>::AABBBinaryTreeNode::ScalarValue(void) const {
#ifdef __VCGLIB_AABBBINARYTREE_SV
	return (this->scalarValue);
#else
	assert(0);
	return (this->boxCenter[0]);
#endif
}

}	// end namespace vcg

#endif // #ifndef TREE_TYPE
