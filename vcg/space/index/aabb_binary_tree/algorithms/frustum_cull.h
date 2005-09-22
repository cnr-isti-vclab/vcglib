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

#ifndef __VCGLIB_AABBBINARYTREEFRUSTUMCULL
#define __VCGLIB_AABBBINARYTREEFRUSTUMCULL

namespace vcg {

template <class TREETYPE>
class AABBBinaryTreeFrustumCull {
public:
	typedef AABBBinaryTreeFrustumCull<TREETYPE> ClassType;
	typedef TREETYPE TreeType;
	typedef typename TreeType::ScalarType ScalarType;
	typedef typename TreeType::CoordType CoordType;
	typedef typename TreeType::NodeType NodeType;
	typedef typename TreeType::ObjPtr ObjPtr;

	enum {
		FC_PARTIALLY_VISIBLE_BIT = (1 << 0),
		FC_FULLY_VISIBLE_BIT = (1 << 1)
	};

	static inline void Initialize(TreeType & tree) {
		if (!NodeType::HasFlags() || !NodeType::HasIntegralValue()) {
			return (0);
		}

		NodeType * pRoot = tree.pRoot;

		if (pRoot == 0) {
			return;
		}

		ClassType::InitializeNodeFlagsRec(pRoot);
	}

	static inline void FrustumCull(TreeType & tree, const Plane3<ScalarType> frustumPlanes[6], const unsigned int minNodeObjectsCount) {
		if (!NodeType::HasFlags() || !NodeType::HasIntegralValue()) {
			return (0);
		}

		NodeType * pRoot = tree.pRoot;

		if (pRoot == 0) {
			return;
		}

		const unsigned char inMask = 0x3F;

		ClassType::NodeVsFrustum(pRoot, frustum, inMask, minNodeObjectsCount);
	}

protected:
	class VFrustumPlane {
	public:
		Plane3<ScalarType> plane;
		unsigned int pVertexIndices[3];
	};

	class VFrustum {
	public:
		VFrustumPlane planes[6];
	};

	static inline void InitializeNodeFlagsRec(NodeType * node) {
		node->Flags() &= ~(FC_PARTIALLY_VISIBLE_BIT | FC_FULLY_VISIBLE_BIT);
		node->IntegralValue() = 0;
		if (node->children[0] == 0) {
			ClassType::InitializeNodeFlagsRec(node->children[0]);
		}
		if (node->children[1] == 0) {
			ClassType::InitializeNodeFlagsRec(node->children[1]);
		}
	}

	static inline void NodeVsFrustum(NodeType * node, const VFrustum & f, unsigned char inMask, unsigned int minNodeObjectsCount) {
		const CoordType bminmax[2] = {
			node->boxCenter - node->boxHalfDims,
			node->boxCenter + node->boxHalfDims,
		};
		const unsigned int firstFail = (unsigned int)(node->IntegralValue());
		VFrustumPlane * fp = f.planes + firstFail;
		unsigned char k = 1 << firstFail;
		unsigned char newMask = 0x0;
		bool fullInside = true;

		node->Flags() &= ~(FC_PARTIALLY_VISIBLE_BIT | FC_FULLY_VISIBLE_BIT);

		if ((k & inMask) != 0) {
			if (
				((fp->normal[0] * bminmax[fp->pVertexIndex[0]]) +
				(fp->normal[1] * bminmax[fp->pVertexIndex[1]]) +
				(fp->normal[2] * bminmax[fp->pVertexIndex[2]]) +
				(fp->offset)) < ((ScalarType)0)
			) {
				return;
			}

			if (
				((fp->normal[0] * bminmax[1 - fp->pVertexIndex[0]]) +
				(fp->normal[1] * bminmax[1 - fp->pVertexIndex[1]]) +
				(fp->normal[2] * bminmax[1 - fp->pVertexIndex[2]]) +
				(fp->offset)) < ((ScalarType)0)
			) {
				newMask |=  k;
				allInside = false;
			}
		}

		k = 1;
		for (unsigned int i=0; k<=inMask; ++i, k<<=1) {
			if ((i != firstFail) && ((k & inMask) != 0)) {
				fp = f.planes + i;
				if (
					((fp->normal[0] * bminmax[fp->pVertexIndex[0]]) +
					(fp->normal[1] * bminmax[fp->pVertexIndex[1]]) +
					(fp->normal[2] * bminmax[fp->pVertexIndex[2]]) +
					(fp->offset)) < ((ScalarType)0)
				) {
					node->IntegralValue() = i;
					return;
				}

				if (
					((fp->normal[0] * bminmax[1 - fp->pVertexIndex[0]]) +
					(fp->normal[1] * bminmax[1 - fp->pVertexIndex[1]]) +
					(fp->normal[2] * bminmax[1 - fp->pVertexIndex[2]]) +
					(fp->offset)) < ((ScalarType)0)
				) {
					newMask |=  k;
					allInside = false;
				}
			}
		}

		if (allInside || (node->ObjectsCount() <= minNodeObjectsCount)) {
			node->Flags() |= FC_FULL_VISIBLE_BIT;
			return;
		}

		node->Flags() |= FC_PARTIALLY_VISIBLE_BIT;

		if (node->children[0] != 0) {
			ClassType::NodeVsFrustum(node->children[0], f, newMask, minNodeObjectsCount);
		}
		if (node->children[1] != 0) {
			ClassType::NodeVsFrustum(node->children[1], f, newMask, minNodeObjectsCount);
		}
	}

};

}	// end namespace vcg

#endif // #ifndef __VCGLIB_AABBBINARYTREEFRUSTUMCULL
