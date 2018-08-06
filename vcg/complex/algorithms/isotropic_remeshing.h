/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2017                                           \/)\/    *
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
#ifndef _VCG_ISOTROPICREMESHING_H
#define _VCG_ISOTROPICREMESHING_H

#include<vcg/complex/algorithms/update/quality.h>
#include<vcg/complex/algorithms/update/curvature.h>
#include<vcg/complex/algorithms/update/normal.h>
#include<vcg/complex/algorithms/refine.h>
#include<vcg/complex/algorithms/stat.h>
#include<vcg/complex/algorithms/smooth.h>
#include<vcg/complex/algorithms/local_optimization/tri_edge_collapse.h>
#include<vcg/space/index/spatial_hashing.h>

namespace vcg {
namespace tri {
template<class TRI_MESH_TYPE>
class IsotropicRemeshing
{
public:
	typedef TRI_MESH_TYPE MeshType;
	typedef typename MeshType::FaceType FaceType;
	typedef typename FaceType::VertexType VertexType;
	typedef	typename VertexType::ScalarType ScalarType;
	typedef typename face::Pos<FaceType> PosType;
	typedef BasicVertexPair<VertexType> VertexPair;
	typedef EdgeCollapser<MeshType, VertexPair> Collapser;
	typedef GridStaticPtr<FaceType, ScalarType> StaticGrid;


	typedef struct Params {
		typedef struct Stat {
			int splitNum;
			int collapseNum;
			int flipNum;

			void Reset() {
				splitNum=0;
				collapseNum=0;
				flipNum=0;
			}
		} Stat;


		ScalarType minLength; // minimal admitted length: no edge should be shorter than this value (used when collapsing)
		ScalarType maxLength; // maximal admitted length: no edge should be longer than this value  (used when refining)
		ScalarType lengthThr;
		ScalarType creaseAngleRadThr= math::ToRad(10.0) ;
		ScalarType creaseAngleCosThr= cos(math::ToRad(10.0)) ;
		bool splitFlag    = true;
		bool swapFlag     = true;
		bool collapseFlag = true;
		bool smoothFlag=true;
		bool projectFlag=true;
		bool selectedOnly = false;
		bool adapt=false;
		int iter=1;
		Stat stat;
		void SetTargetLen(ScalarType len)
		{
			minLength=len;
			maxLength=len*2.0;
			lengthThr=len*2.0;
		}
		void SetFeatureAngleDeg(ScalarType angle)
		{
			creaseAngleRadThr =  math::ToRad(angle);
			creaseAngleCosThr= cos(creaseAngleRadThr);
		}


	} Params;

	static void Do(MeshType &toRemesh, Params params, vcg::CallBackPos * cb=0)
	{
		MeshType toProjectCopy;
		tri::UpdateBounding<MeshType>::Box(toRemesh);
		tri::UpdateNormal<MeshType>::PerVertexNormalizedPerFaceNormalized(toRemesh);

		tri::Append<MeshType,MeshType>::MeshCopy(toProjectCopy, toRemesh);

		Do(toRemesh,toProjectCopy,params,cb);
	}
	static void Do(MeshType &toRemesh, MeshType &toProject, Params params, vcg::CallBackPos * cb=0)
	{
		assert(&toRemesh != &toProject);
		StaticGrid t;
		params.stat.Reset();
		t.Set(toProject.face.begin(), toProject.face.end());
		tri::FaceTmark<MeshType> mark;
		mark.SetMesh(&toProject);

		tri::UpdateTopology<MeshType>::FaceFace(toRemesh);
		tri::UpdateFlags<MeshType>::VertexBorderFromFaceAdj(toRemesh);

		tri::MeshAssert<MeshType>::FFTwoManifoldEdge(toRemesh);
		tri::UpdateTopology<MeshType>::VertexFace(toRemesh);
		computeQuality(toRemesh);
		tri::UpdateQuality<MeshType>::VertexSaturate(toRemesh);

		for(int i=0; i < params.iter; ++i)
		{
			params.stat.Reset();

			if(cb) cb(100*i/params.iter, "Remeshing");

			if(params.splitFlag)
				SplitLongEdges(toRemesh, params);

			if(params.swapFlag)
				ImproveValence(toRemesh, params);

			if(params.collapseFlag)
			{
				CollapseShortEdges(toRemesh, params);
				CollapseCrosses(toRemesh, params);
			}

			if(params.smoothFlag)
				ImproveByLaplacian(toRemesh, params);
			if(params.projectFlag)
				ProjectToSurface(toRemesh, t, mark);

		}
	}

private:
	IsotropicRemeshing() {}
	// this returns the value of cos(a) where a is the angle between n0 and n1. (scalar prod is cos(a))
	static inline ScalarType fastAngle(Point3<ScalarType> n0, Point3<ScalarType> n1)
	{
		return math::Clamp(n0*n1,(ScalarType)-1.0,(ScalarType)1.0);
	}
	// compare the value of the scalar prod with the cos of the crease threshold
	static inline bool testCreaseEdge(PosType &p, ScalarType creaseCosineThr)
	{
		ScalarType angle = fastAngle(NormalizedTriangleNormal(*(p.F())), NormalizedTriangleNormal(*(p.FFlip())));
		return angle <= creaseCosineThr;
		//        return (angle <= creaseCosineThr && angle >= -creaseCosineThr);
	}
	// this stores in minQ the value of the 10th percentile of the VertQuality distribution and in
	// maxQ the value of the 90th percentile.
	static inline void computeVQualityDistrMinMax(MeshType &m, ScalarType &minQ, ScalarType &maxQ)
	{
		Distribution<ScalarType> distr;
		tri::Stat<MeshType>::ComputePerVertexQualityDistribution(m,distr);

		maxQ = distr.Percentile(0.9f);
		minQ = distr.Percentile(0.1f);
	}

	//Computes PerVertexQuality as a function of the 'deviation' of the normals taken from
	//the faces incident to each vertex
	static void computeQuality(MeshType &m)
	{
		tri::RequirePerVertexQuality(m);
		tri::UpdateFlags<MeshType>::VertexClearV(m);

		for(auto vi=m.vert.begin(); vi!=m.vert.end(); ++vi)
			if(!(*vi).IsD())
			{
				vector<FaceType*> ff;
				face::VFExtendedStarVF(&*vi, 0, ff);

				ScalarType tot = 0.f;
				auto it = ff.begin();
				Point3<ScalarType> fNormal = NormalizedTriangleNormal(**it);
				++it;
				while(it != ff.end())
				{
					tot+= 1-math::Abs(fastAngle(fNormal, NormalizedTriangleNormal(**it)));
					++it;
				}
				vi->Q() = tot / (ScalarType)(std::max(1, ((int)ff.size()-1)));
				vi->SetV();
			}
	}

	/*
	 Computes the ideal valence for the vertex in pos p:
	 4 for border vertices
	 6 for internal vertices
	*/
	static inline int idealValence(PosType &p)
	{
		if(p.IsBorder()) return 4;
		return 6;
	}
	static inline int idealValence(VertexType &v)
	{
		if(v.IsB()) return 4;
		return 6;
	}
	static inline int idealValenceSlow(PosType &p)
	{
		std::vector<PosType> posVec;
		VFOrderedStarFF(p,posVec);
		float angleSumRad =0;
		for(PosType &ip : posVec)
		{
			angleSumRad += ip.AngleRad();
		}

		return (int)(std::ceil(angleSumRad / (M_PI/3.0f)));
	}
	/*
		Edge Swap Step:
		This method optimizes the valence of each vertex.
		oldDist is the sum of the absolute distance of each vertex from its ideal valence
		newDist is the sum of the absolute distance of each vertex from its ideal valence after
		the edge swap.
		If the swap decreases the total absolute distance, then it's applied, preserving the triangle
		quality.                        +1
				v1                    v1
			   / \                   /|\
			  /   \                 / | \
			 /     \               /  |  \
			/    _*p\           -1/   |   \ -1
		  v2--------v0 ========> v2   |   v0
		   \        /             \   |   /
			\      /               \  |  /
			 \    /                 \ | /
			  \  /                   \|/ +1
			   v3                     v3
			Before Swap             After Swap
	*/
	static bool testSwap(PosType p, ScalarType creaseAngleCosThr)
	{
		//if border or feature, do not swap
		if(p.IsBorder() || testCreaseEdge(p, creaseAngleCosThr)) return false;

		int oldDist = 0, newDist = 0, idealV, actualV;

		PosType tp=p;

		VertexType *v0=tp.V();
		idealV  = idealValence(tp); actualV = tp.NumberOfIncidentVertices();
		oldDist += abs(idealV - actualV); newDist += abs(idealV - (actualV - 1));

		tp.FlipF();tp.FlipE();tp.FlipV();
		VertexType *v1=tp.V();
		idealV  = idealValence(tp); actualV = tp.NumberOfIncidentVertices();
		oldDist += abs(idealV - actualV); newDist += abs(idealV - (actualV + 1));

		tp.FlipE();tp.FlipV();tp.FlipE();
		VertexType *v2=tp.V();
		idealV  = idealValence(tp); actualV = tp.NumberOfIncidentVertices();
		oldDist += abs(idealV - actualV); newDist += abs(idealV - (actualV - 1));

		tp.FlipF();tp.FlipE();tp.FlipV();
		VertexType *v3=tp.V();
		idealV  = idealValence(tp); actualV = tp.NumberOfIncidentVertices();
		oldDist += abs(idealV - actualV); newDist += abs(idealV - (actualV + 1));

		ScalarType qOld = std::min(Quality(v0->P(),v2->P(),v3->P()),Quality(v0->P(),v1->P(),v2->P()));
		ScalarType qNew = std::min(Quality(v0->P(),v1->P(),v3->P()),Quality(v2->P(),v3->P(),v1->P()));

		return (newDist < oldDist && qNew >= qOld * 0.50f) ||
		        (newDist == oldDist && qNew > qOld * 1.f) || qNew > 1.5f * qOld;
	}


	// Edge swap step: edges are flipped in order to optimize valence and triangle quality across the mesh
	static void ImproveValence(MeshType &m, Params &params)
	{
		tri::UpdateTopology<MeshType>::FaceFace(m); //collapser does not update FF
		ForEachFacePos(m, [&](PosType &p){
			if(p.FFlip() > p.F())
				if(((!params.selectedOnly) || (p.F()->IsS() && p.FFlip()->IsS())) &&
				   testSwap(p, params.creaseAngleCosThr) &&
				   face::CheckFlipEdgeNormal(*p.F(), p.E(), math::ToRad(10.f)) &&
				   face::CheckFlipEdge(*p.F(), p.E()) )
				{
					face::FlipEdge(*p.F(), p.E());
					++params.stat.flipNum;
				}
		});
	}

	// The predicate that defines which edges should be split
	class EdgeSplitAdaptPred
	{
	public:
		int count = 0;
		ScalarType length, lengthThr, minQ, maxQ;
		bool operator()(PosType &ep)
		{
			ScalarType mult = math::ClampedLerp((ScalarType)0.5,(ScalarType)1.5, (((math::Abs(ep.V()->Q())+math::Abs(ep.VFlip()->Q()))/(ScalarType)2.0)/(maxQ-minQ)));
			ScalarType dist = Distance(ep.V()->P(), ep.VFlip()->P());
			if(dist > std::max(mult*length,lengthThr*2))
			{
				++count;
				return true;
			}
			else
				return false;
		}
	};

	class EdgeSplitLenPred
	{
	public:
		int count = 0;
		ScalarType squaredlengthThr;
		bool operator()(PosType &ep)
		{
			if(SquaredDistance(ep.V()->P(), ep.VFlip()->P()) > squaredlengthThr)
			{
				++count;
				return true;
			}
			else
				return false;
		}
	};

	//Split pass: This pass uses the tri::RefineE from the vcglib to implement
	//the refinement step, using EdgeSplitPred as a predicate to decide whether to split or not
	static void SplitLongEdges(MeshType &m, Params &params)
	{
		tri::UpdateTopology<MeshType>::FaceFace(m);
		tri::MidPoint<MeshType> midFunctor(&m);

		ScalarType minQ,maxQ;
		if(params.adapt){
			computeVQualityDistrMinMax(m, minQ, maxQ);
			EdgeSplitAdaptPred ep;
			ep.minQ      = minQ;
			ep.maxQ      = maxQ;
			ep.length    = params.maxLength;
			ep.lengthThr = params.lengthThr;
			tri::RefineE(m,midFunctor,ep);
			params.stat.splitNum+=ep.count;
		}
		else {
			EdgeSplitLenPred ep;
			ep.squaredlengthThr = params.maxLength*params.maxLength;
			tri::RefineE(m,midFunctor,ep,params.selectedOnly);
			params.stat.splitNum+=ep.count;
		}
	}

	//Geometric check on feasibility of the collapse of the given pos
	//The check fails if:
	//  -new face has too bad quality.
	//  -new face normal changes too much after collapse.
	//  -new face has too long edges.
	// TRY: if the vertex has valence 4 (cross vertex) we relax the check on length
	static bool checkCollapseFacesAroundVert(PosType &p, Point3<ScalarType> &mp, Params & params, bool relaxed=false)
	{
		ScalarType minimalAdmittedArea = (params.minLength * params.minLength)/10000.0;
		vector<FaceType*> ff;
		vector<int> vi;

		face::VFStarVF<FaceType>(p.V(), ff, vi);
		bool allIncidentFaceSelected = true;
		for(FaceType *f: ff)
			if(!(*f).IsD() && f != p.F()) //i'm not a deleted face
			{
				allIncidentFaceSelected &= f->IsS();

				PosType pi(f, p.V()); //same vertex

				VertexType *v0 = pi.V();
				VertexType *v1 = pi.F()->V1(pi.VInd());
				VertexType *v2 = pi.F()->V2(pi.VInd());

				if( v1 == p.VFlip() || v2 == p.VFlip()) //i'm the other deleted face
					continue;

				float area = DoubleArea(*(pi.F()))/2.f;

				//quality and normal divergence checks
				ScalarType newQ = Quality(mp,      v1->P(), v2->P());
				ScalarType oldQ = Quality(v0->P(), v1->P(), v2->P());

				if(area > minimalAdmittedArea) // for triangles not too small
				{
					if( newQ <= 0.5*oldQ  )
						return false;

					Point3<ScalarType> oldN = NormalizedTriangleNormal(*(pi.F()));
					Point3<ScalarType> newN = Normal(mp, v1->P(), v2->P()).Normalize();
					float div = fastAngle(oldN, newN);
					//                if(AngleN(oldN,newN) > math::ToRad(1.0)) return false;
					if(div <= params.creaseAngleCosThr )
						return false;

					// we prevent collapse that makes edges too long (except for cross)
					if(!relaxed)
						if((Distance(mp, v1->P()) > params.maxLength || Distance(mp, v2->P()) > params.maxLength))
							return false;
				}
			}
		if(params.selectedOnly) return allIncidentFaceSelected;
		return true;
	}

	// Collapse test: Usual collapse test (check on target length) plus borders and crease handling
	// and adaptivity.
	static bool testCollapse(PosType &p, Point3<ScalarType> &mp, ScalarType minQ, ScalarType maxQ, Params &params, bool relaxed = false)
	{
		ScalarType mult = (params.adapt) ? math::ClampedLerp((ScalarType)0.5,(ScalarType)1.5, (((math::Abs(p.V()->Q())+math::Abs(p.VFlip()->Q()))/(ScalarType)2.0)/(maxQ-minQ))) : (ScalarType)1;
		ScalarType dist = Distance(p.V()->P(), p.VFlip()->P());
		ScalarType thr = mult*params.minLength;
		ScalarType area = DoubleArea(*(p.F()))/2.f;
		if(dist < thr || area < params.minLength*params.minLength/100.f)//if to collapse
		{
			PosType pp = p; p.FlipV();
			//check all faces around p() and p.vflip()
			return checkCollapseFacesAroundVert(p,  mp, params, relaxed) &&
			        checkCollapseFacesAroundVert(pp, mp, params, relaxed);
		}
		return false;
	}
	//This function is especially useful to enforce feature preservation during collapses
	//of boundary edges in planar or near planar section of the mesh
	static bool chooseBoundaryCollapse(PosType &p, VertexPair &pair)
	{
		Point3<ScalarType> collapseNV, collapsedNV0, collapsedNV1;
		collapseNV = (p.V()->P() - p.VFlip()->P()).normalized();

		vector<VertexType*> vv;
		face::VVStarVF<FaceType>(p.V(), vv);

		for(VertexType *v: vv)
			if(!(*v).IsD() && (*v).IsB()) //ignore non border
				collapsedNV0 = ((*v).P() - p.VFlip()->P()).normalized(); //edge vector after collapse

		face::VVStarVF<FaceType>(p.VFlip(), vv);

		for(VertexType *v: vv)
			if(!(*v).IsD() && (*v).IsB()) //ignore non border
				collapsedNV1 = ((*v).P() - p.V()->P()).normalized(); //edge vector after collapse

		float cosine = cos(math::ToRad(1.5f));
		float angle0 = fabs(fastAngle(collapseNV, collapsedNV0));
		float angle1 = fabs(fastAngle(collapseNV, collapsedNV1));
		//if on both sides we deviate too much after collapse => don't collapse
		if(angle0 <= cosine && angle1 <= cosine)
			return false;
		//choose the best collapse (the more parallel one to the previous edge..)
		pair = (angle0 >= angle1) ? VertexPair(p.V(), p.VFlip()) : VertexPair(p.VFlip(), p.V());
		return true;
	}

	//The actual collapse step: foreach edge it is collapse iff TestCollapse returns true AND
	// the linkConditions are preserved
	static void CollapseShortEdges(MeshType &m, Params &params)
	{
		ScalarType minQ, maxQ;
		int candidates = 0;

		if(params.adapt)
			computeVQualityDistrMinMax(m, minQ, maxQ);

		tri::UpdateTopology<MeshType>::VertexFace(m);
		tri::UpdateFlags<MeshType>::FaceBorderFromVF(m);
		tri::UpdateFlags<MeshType>::VertexBorderFromFaceBorder(m);

		for(auto fi=m.face.begin(); fi!=m.face.end(); ++fi)
			if(!(*fi).IsD() && (params.selectedOnly == false || fi->IsS()))
			{
				for(auto i=0; i<3; ++i)
				{
					PosType pi(&*fi, i);
					++candidates;
					VertexPair  bp = VertexPair(pi.V(), pi.VFlip());
					Point3<ScalarType> mp = (pi.V()->P()+pi.VFlip()->P())/2.f;

					bool boundary = false;
					//if both border or both internal use midpoint
					if(pi.V()->IsB() == pi.VFlip()->IsB())
					{
						//if both border but not on border edge..abort..if chooseBoundaryCollapse can't find good collapse..abort..
						if(pi.V()->IsB() && !fi->IsB(i) && !(boundary = chooseBoundaryCollapse(pi, bp)))
							continue;
						mp = (pi.V()->IsB()) ? bp.V(1)->P() : (pi.V()->P()+pi.VFlip()->P())/2.f;

					}
					else //if only one is border...collapse on it!
					{
						bp = (pi.V()->IsB()) ? VertexPair(pi.VFlip(), pi.V()) : VertexPair(pi.V(), pi.VFlip());
						mp = (pi.V()->IsB()) ? pi.V()->P() : pi.VFlip()->P();
					}

					if(testCollapse(pi, mp, minQ, maxQ, params, boundary) && Collapser::LinkConditions(bp))
					{
						Collapser::Do(m, bp, mp);
						++params.stat.collapseNum;
						break;
					}

				}
			}
		Allocator<MeshType>::CompactEveryVector(m);
	}


	//Here I just need to check the faces of the cross, since the other faces are not
	//affected by the collapse of the internal faces of the cross.
	static bool testCrossCollapse(PosType &p, Point3<ScalarType> &mp, Params &params)
	{
		if(!checkCollapseFacesAroundVert(p, mp, params, true))
			return false;
		return true;
	}

	//Choose the best way to collapse a cross based on the (external) cross vertices valence
	//and resulting face quality
	//                                      +0                   -1
	//             v1                    v1                    v1
	//            /| \                   /|\                  / \
	//           / |  \                 / | \                /   \
	//          /  |   \               /  |  \              /     \
	//         / *p|    \           -1/   |   \ -1       +0/       \+0
	//       v0-------- v2 ========> v0   |   v2    OR    v0-------v2
	//        \    |    /             \   |   /            \       /
	//         \   |   /               \  |  /              \     /
	//          \  |  /                 \ | /                \   /
	//           \ | /                   \|/ +0               \ / -1
	//             v3                     v3                   v3
	static VertexPair chooseBestCrossCollapse(PosType &p, vector<FaceType*> &ff)
	{
		vector<VertexType*> vv0, vv1, vv2, vv3;
		VertexType *v0, *v1, *v2, *v3;

		v0 = p.F()->V1(p.VInd());
		v1 = p.F()->V2(p.VInd());

		for(FaceType *f: ff)
			if(!(*f).IsD() && f != p.F())
			{
				PosType pi(f, p.V());
				VertexType *fv1 = pi.F()->V1(pi.VInd());
				VertexType *fv2 = pi.F()->V2(pi.VInd());

				if(fv1 == v0 || fv2 == v0)
					v3 = (fv1 == v0) ? fv2 : fv1;
				if(fv1 == v1 || fv2 == v1)
					v2 = (fv1 == v1) ? fv2 : fv1;
			}

		face::VVStarVF<FaceType>(v0, vv0);
		face::VVStarVF<FaceType>(v1, vv1);
		face::VVStarVF<FaceType>(v2, vv2);
		face::VVStarVF<FaceType>(v3, vv3);


		int nv0 = vv0.size(), nv1 = vv1.size();
		int nv2 = vv2.size(), nv3 = vv3.size();

		int delta1 = (idealValence(*v0) - nv0) + (idealValence(*v2) - nv2);
		int delta2 = (idealValence(*v1) - nv1) + (idealValence(*v3) - nv3);

		ScalarType Q1 = std::min(Quality(v0->P(), v1->P(), v3->P()), Quality(v1->P(), v2->P(), v3->P()));
		ScalarType Q2 = std::min(Quality(v0->P(), v1->P(), v2->P()), Quality(v2->P(), v3->P(), v0->P()));

		if(delta1 < delta2 && Q1 >= 0.6f*Q2)
			return VertexPair(p.V(), v1);
		else
			return VertexPair(p.V(), v0);
	}
	//Cross Collapse pass: This pass cleans the mesh from cross vertices, keeping in mind the link conditions
	//and feature preservations tests.
	static void CollapseCrosses(MeshType &m , Params &params)
	{
		tri::UpdateTopology<MeshType>::ClearFaceFace(m);
		tri::UpdateTopology<MeshType>::VertexFace(m);
		int count = 0;

		for(auto fi=m.face.begin(); fi!=m.face.end(); ++fi)
			if(!(*fi).IsD() && (params.selectedOnly == false || fi->IsS()))
			{
				for(auto i=0; i<3; ++i)
				{
					PosType pi(&*fi, i);
					if(!pi.V()->IsB())
					{
						vector<FaceType*> ff;
						vector<int> vi;
						face::VFStarVF<FaceType>(pi.V(), ff, vi);

						//removing crosses and tricuspidis only
						if(ff.size() == 4 || ff.size() == 3)
						{
							VertexPair bp  = (ff.size() == 4) ? chooseBestCrossCollapse(pi, ff) : VertexPair(pi.V(), pi.VFlip());
							Point3<ScalarType> mp = bp.V(1)->P();
							//todo: think about if you should try doing the other collapse if test or link fails for this one
							if(testCrossCollapse(pi, mp, params) && Collapser::LinkConditions(bp))
							{
								Collapser::Do(m, bp, mp);
								++count;
								break;
							}
						}
					}
				}
			}
		Allocator<MeshType>::CompactEveryVector(m);
	}

	// This function sets the selection bit on vertices that lie on creases
	static int selectVertexFromCrease(MeshType &m, ScalarType creaseThr)
	{
		int count = 0;

		ForEachFacePos(m, [&](PosType &p){
			if((p.FFlip() > p.F()) && testCreaseEdge(p, creaseThr))
			{
				p.V()->SetS();
				p.VFlip()->SetS();
				++count;
			}
		});

		return count;
	}
	/**
	  * Simple Laplacian Smoothing step
	  * Border and crease vertices are kept fixed.
	  * If there are selected faces and the param.onlySelected is true we compute
	  * the set of internal vertices to the selection and we combine it in and with
	  * the vertexes not on border or creases
	*/
	static void ImproveByLaplacian(MeshType &m, Params params)
	{
		SelectionStack<MeshType> ss(m);

		if(params.selectedOnly) {
			ss.push();
			tri::UpdateSelection<MeshType>::VertexFromFaceStrict(m);
			ss.push();
		}
		tri::UpdateTopology<MeshType>::FaceFace(m);
		tri::UpdateFlags<MeshType>::VertexBorderFromFaceAdj(m);
		tri::UpdateSelection<MeshType>::VertexFromBorderFlag(m);
		selectVertexFromCrease(m, params.creaseAngleCosThr);
		tri::UpdateSelection<MeshType>::VertexInvert(m);
		if(params.selectedOnly) {
			ss.popAnd();
		}
		tri::Smooth<MeshType>::VertexCoordPlanarLaplacian(m,1, params.creaseAngleRadThr,true);
		tri::UpdateSelection<MeshType>::VertexClear(m);
		if(params.selectedOnly) {
			ss.pop();
		}
	}
	/*
		Reprojection step, this method reprojects each vertex on the original surface
		sampling the nearest Point3 onto it using a uniform grid StaticGrid t
	*/
	static void ProjectToSurface(MeshType &m, StaticGrid t, FaceTmark<MeshType> mark)
	{
		face::PointDistanceBaseFunctor<ScalarType> distFunct;
		ScalarType maxDist = std::numeric_limits<ScalarType>::max(), minDist = 0.f;
		for(auto vi=m.vert.begin();vi!=m.vert.end();++vi)
			if(!(*vi).IsD())
			{
				Point3<ScalarType> newP;
				t.GetClosest(distFunct, mark, vi->P(), maxDist, minDist, newP);
				vi->P() = newP;
			}
	}
};
} // end namespace tri
} // end namespace vcg
#endif
