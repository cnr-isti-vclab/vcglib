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

#include <vcg/complex/algorithms/update/quality.h>
#include <vcg/complex/algorithms/update/curvature.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/refine.h>
#include <vcg/complex/algorithms/stat.h>
#include <vcg/complex/algorithms/smooth.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse.h>
#include <vcg/complex/algorithms/geodesic.h>
#include <vcg/space/index/spatial_hashing.h>
#include <vcg/complex/append.h>
#include <vcg/complex/allocate.h>
#include <wrap/io_trimesh/export.h>

namespace vcg {
namespace tri {
template<class TRI_MESH_TYPE>
class IsotropicRemeshing
{
public:
    typedef TRI_MESH_TYPE MeshType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::FacePointer FacePointer;
    typedef typename FaceType::VertexType VertexType;
    typedef typename FaceType::VertexPointer VertexPointer;
    typedef	typename VertexType::ScalarType ScalarType;
    typedef	typename VertexType::CoordType CoordType;
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

        ScalarType minAdaptiveMult = 1;
        ScalarType maxAdaptiveMult = 1;

        ScalarType minimalAdmittedArea;
        ScalarType maxSurfDist;

        ScalarType aspectRatioThr  = 0.05;                    //min aspect ratio: during relax bad triangles will be relaxed
        ScalarType foldAngleCosThr = cos(math::ToRad(140.));   //min angle to be considered folding: during relax folded triangles will be relaxed

        ScalarType creaseAngleRadThr = math::ToRad(10.0);
        ScalarType creaseAngleCosThr = cos(math::ToRad(10.0)); //min angle to be considered crease: two faces with normals diverging more than thr share a crease edge

        bool splitFlag    = true;
        bool swapFlag     = true;
        bool collapseFlag = true;
        bool smoothFlag=true;
        bool projectFlag=true;
        bool selectedOnly = false;
        bool cleanFlag = true;

        bool userSelectedCreases = false;
        bool surfDistCheck = true;

        bool adapt=false;
        int iter=1;
        Stat stat;

        void SetTargetLen(const ScalarType len)
        {
            minLength=len*4./5.;
            maxLength=len*4./3.;
            lengthThr=len*4./3.;
            minimalAdmittedArea = (minLength * minLength)/1000.0;
        }

        void SetFeatureAngleDeg(const ScalarType angle)
        {
            creaseAngleRadThr =  math::ToRad(angle);
            creaseAngleCosThr = cos(creaseAngleRadThr);
        }

        StaticGrid grid;
        MeshType* m;
        MeshType* mProject;

    } Params;

private:
    static void debug_crease (MeshType & toRemesh, std::string  prepend, int i)
    {
        ForEachVertex(toRemesh, [] (VertexType & v) {
            v.C() = Color4b::Gray;
            v.Q() = 0;
        });

        ForEachFacePos(toRemesh, [&](PosType &p){
            if (p.F()->IsFaceEdgeS(p.E()))
            {
                p.V()->Q() += 1;
                p.VFlip()->Q() += 1;
            }
        });

        ForEachVertex(toRemesh, [] (VertexType & v) {
            if (v.Q() >= 4)
                v.C() = Color4b::Green;
            else if (v.Q() >= 2)
                v.C() = Color4b::Red;
        });
        prepend += "_creases" + std::to_string(i) + ".ply";
        vcg::tri::io::Exporter<MeshType>::Save(toRemesh, prepend.c_str(), vcg::tri::io::Mask::IOM_ALL);
    }


    static void removeColinearFaces(MeshType & m, Params & params)
    {
        vcg::tri::UpdateTopology<MeshType>::FaceFace(m);

        int count = 0;
        int iter = 0;
        do
        {
            vcg::tri::UpdateTopology<MeshType>::FaceFace(m);
            vcg::tri::UnMarkAll(m);

            count = 0;
            for (size_t i = 0; i < size_t(m.FN()); ++i)
            {
                FaceType & f = m.face[i];

                ScalarType quality = vcg::QualityRadii(f.cP(0), f.cP(1), f.cP(2));

                if (quality <= 0.001)
                {
                    //find longest edge
                    double edges[3];
                    edges[0] = vcg::Distance(f.cP(0), f.cP(1));
                    edges[1] = vcg::Distance(f.cP(1), f.cP(2));
                    edges[2] = vcg::Distance(f.cP(2), f.cP(0));

                    ScalarType smallestEdge = std::min(edges[0], std::min(edges[1], edges[2]));
                    int longestIdx = std::find(edges, edges+3, std::max(std::max(edges[0], edges[1]), edges[2])) - (edges);

                    if (vcg::tri::IsMarked(m, f.V2(longestIdx)))
                        continue;


                    auto f1 = f.cFFp(longestIdx);
                    vcg::tri::Mark(m,f.V2(longestIdx));
                    if (!vcg::face::IsBorder(f, longestIdx) && vcg::face::IsManifold(f, longestIdx) && vcg::face::checkFlipEdgeNotManifold<FaceType>(f, longestIdx))  {

                        // Check if EdgeFlipping improves quality
                        FacePointer g = f.FFp(longestIdx); int k = f.FFi(longestIdx);
                        vcg::Triangle3<ScalarType> t1(f.P(longestIdx), f.P1(longestIdx), f.P2(longestIdx)), t2(g->P(k), g->P1(k), g->P2(k)),
                                t3(f.P(longestIdx), g->P2(k), f.P2(longestIdx)), t4(g->P(k), f.P2(longestIdx), g->P2(k));

                        auto n1 = vcg::TriangleNormal(t1);
                        auto n2 = vcg::TriangleNormal(t2);
                        auto n3 = vcg::TriangleNormal(t3);
                        auto n4 = vcg::TriangleNormal(t4);

                        auto biggestSmallest = vcg::DoubleArea(t1) > vcg::DoubleArea(t2) ? std::make_pair(t1, t2) : std::make_pair(t2, t1);
                        auto areaRatio = vcg::DoubleArea(biggestSmallest.first) / vcg::DoubleArea(biggestSmallest.second);

                        bool normalCheck = true;
                        //                        if (n1.Norm() > 0.001 && n2.Norm() > 0.001)
                        {
                            auto referenceNormal = vcg::NormalizedTriangleNormal(biggestSmallest.first);

                            normalCheck &= vcg::NormalizedTriangleNormal(t3) * referenceNormal >= 0.95;
                            normalCheck &= vcg::NormalizedTriangleNormal(t4) * referenceNormal >= 0.95;
                        }

                        bool areaCheck = false;
                        if (areaRatio > 1000)
                        {
                            areaCheck |= vcg::DoubleArea(t3) / vcg::DoubleArea(biggestSmallest.second) > 1000 && vcg::DoubleArea(t4) / vcg::DoubleArea(biggestSmallest.second) > 1000;
                        }

                        if ((normalCheck) && (areaCheck || std::min( QualityFace(t1), QualityFace(t2) ) <= std::min( QualityFace(t3), QualityFace(t4))))
                        {
                            ScalarType dist;
                            CoordType closest;
                            auto fp0 = vcg::tri::GetClosestFaceBase(*params.mProject, params.grid, vcg::Barycenter(t3), smallestEdge/4., dist, closest);
                            if (fp0 == NULL)
                                continue;

                            auto fp1 = vcg::tri::GetClosestFaceBase(*params.mProject, params.grid, vcg::Barycenter(t4), smallestEdge/4., dist, closest);
                            if (fp1 == NULL)
                                continue;

                            vcg::face::FlipEdgeNotManifold<FaceType>(f, longestIdx);
                            ++count;
                        }
                    }
                }
            }
        } while (count && ++iter < 75);

    }

    static void cleanMesh(MeshType & m, Params & params)
    {
        vcg::tri::Clean<MeshType>::RemoveDuplicateFace(m);
        vcg::tri::Clean<MeshType>::RemoveUnreferencedVertex(m);
        vcg::tri::Allocator<MeshType>::CompactEveryVector(m);

        vcg::tri::UpdateTopology<MeshType>::FaceFace(m);
        removeColinearFaces(m, params);
        vcg::tri::UpdateTopology<MeshType>::FaceFace(m);
    }

public:

    static void Do(MeshType &toRemesh, Params & params, vcg::CallBackPos * cb=0)
    {
        MeshType toProjectCopy;
        tri::UpdateBounding<MeshType>::Box(toRemesh);
        tri::UpdateNormal<MeshType>::PerVertexNormalizedPerFaceNormalized(toRemesh);

        tri::Append<MeshType,MeshType>::MeshCopy(toProjectCopy, toRemesh);

        Do(toRemesh,toProjectCopy,params,cb);
    }
    static void Do(MeshType &toRemesh, MeshType &toProject, Params & params, vcg::CallBackPos * cb=0)
    {
        assert(&toRemesh != &toProject);
        params.stat.Reset();


        tri::UpdateBounding<MeshType>::Box(toRemesh);

        {
            tri::UpdateBounding<MeshType>::Box(toProject);
            tri::UpdateNormal<MeshType>::PerFaceNormalized(toProject);
            params.m = &toRemesh;
            params.mProject = &toProject;
            params.grid.Set(toProject.face.begin(), toProject.face.end());
        }

        if (params.cleanFlag)
            cleanMesh(toRemesh, params);

        tri::UpdateTopology<MeshType>::FaceFace(toRemesh);
        tri::UpdateFlags<MeshType>::VertexBorderFromFaceAdj(toRemesh);
        tri::UpdateTopology<MeshType>::VertexFace(toRemesh);

        if (!params.userSelectedCreases)
            tagCreaseEdges(toRemesh, params);

        for(int i=0; i < params.iter; ++i)
        {
            if(cb) cb(100*i/params.iter, "Remeshing");


            if (params.adapt)
            {
                computeQualityDistFromRadii(toRemesh);
                tri::Smooth<MeshType>::VertexQualityLaplacian(toRemesh, 2);
            }

            if(params.splitFlag)
                SplitLongEdges(toRemesh, params);
#ifdef DEBUG_CREASE
            debug_crease(toRemesh, std::string("after_ref"), i);
#endif

            if(params.collapseFlag)
            {
                CollapseShortEdges(toRemesh, params);
                CollapseCrosses(toRemesh, params);
            }

            if(params.swapFlag)
                ImproveValence(toRemesh, params);

            if(params.smoothFlag)
                ImproveByLaplacian(toRemesh, params);

            if(params.projectFlag)
                ProjectToSurface(toRemesh, params);
        }
    }

    static int tagCreaseEdges(MeshType &m, Params & params, bool forceTag = false)
    {
        int count = 0;
        std::vector<char> creaseVerts(m.VN(), 0);

        vcg::tri::UpdateFlags<MeshType>::VertexClearV(m);
        std::queue<PosType> creaseQueue;

        //if we are forcing the crease taggin or we are not using user creases...reset the faceedgeS...
        if ((forceTag || !params.userSelectedCreases))
            vcg::tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(m);

        ForEachFacePos(m, [&](PosType &p){

            if (p.IsBorder())
                p.F()->SetFaceEdgeS(p.E());

            //			if((p.F1Flip() > p.F()))
            {
                FaceType *ff    = p.F();
                FaceType *ffAdj = p.FFlip();

                double quality    = vcg::QualityRadii(ff->cP(0), ff->cP(1), ff->cP(2));
                double qualityAdj = vcg::QualityRadii(ffAdj->cP(0), ffAdj->cP(1), ffAdj->cP(2));

                bool qualityCheck = quality > 0.00000001 && qualityAdj > 0.00000001;
                //				bool areaCheck    = vcg::DoubleArea(*ff) > 0.000001 && vcg::DoubleArea(*ffAdj) > 0.000001;

                if ((forceTag || !params.userSelectedCreases) && (testCreaseEdge(p, params.creaseAngleCosThr) /*&& areaCheck*//* && qualityCheck*/) || p.IsBorder())
                {
                    PosType pp = p;
                    std::vector<FacePointer> faces;
                    std::vector<int> edges;
                    bool allOk = true;

                    do {
                        faces.push_back(pp.F());
                        edges.push_back(pp.E());
                        //                        pp.F()->SetFaceEdgeS(pp.E());
                        if (vcg::QualityRadii(pp.F()->cP(0), pp.F()->cP(1), pp.F()->cP(2)) <= 0.0001)
                        {
                            allOk = false;
                            break;
                        }
                        pp.NextF();
                    } while (pp != p);

                    if (allOk)
                    {
                        for (int i = 0; i < faces.size(); ++i)
                        {
                            faces[i]->SetFaceEdgeS(edges[i]);
                        }
                    }

                    creaseQueue.push(p);
                }
            }
        });
        return count;
    }


private:
    /*
                TODO: Add better crease support: detect all creases at starting time, saving it on facedgesel flags
                          All operations must then preserve the faceedgesel flag accordingly:
                                Refinement -> Check that refiner propagates faceedgesel [should be doing it]
                                Collapse   -> Implement 1D edge collapse and better check on corners and creases
                                Swap       -> Totally avoid swapping crease edges [ok]
                                Smooth     -> Apply 1D smoothing to crease vertices + check on
                                                                (http://www.cs.ubc.ca/labs/imager/tr/2009/eltopo/sisc2009.pdf)
        */
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
        return angle <= creaseCosineThr && angle >= -0.98;
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

    static inline ScalarType computeLengthThrMult(const Params & params, const ScalarType & quality)
    {
        return (params.adapt) ? math::ClampedLerp(params.minAdaptiveMult, params.maxAdaptiveMult, quality) : (ScalarType) 1;
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
                face::VFExtendedStarVF(&*vi, 2, ff);

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
        tri::Smooth<MeshType>::VertexQualityLaplacian(m, 3);
    }

    static void computeQualityDistFromCrease(MeshType & m)
    {
        tri::RequirePerVertexQuality(m);
        tri::UpdateTopology<MeshType>::FaceFace(m);
//        tri::UpdateFlags<MeshType>::VertexClearV(m);
        for (size_t i=0;i<m.vert.size();i++)
            m.vert[i].IMark()=0;

        std::vector<VertexPointer> seeds;
        ForEachFace(m, [&] (FaceType & f) {
           for (int i = 0; i < 3; ++i)
           {
               if (f.IsFaceEdgeS(i))
               {
                   seeds.push_back(f.V0(i));
                   seeds.push_back(f.V1(i));
               }
           }
        });

        std::sort(seeds.begin(),seeds.end());
        auto last=std::unique(seeds.begin(),seeds.end());
        seeds.erase(last, seeds.end());

        tri::EuclideanDistance<MeshType> eu;
        tri::Geodesic<MeshType>::PerVertexDijkstraCompute(m, seeds, eu);
        tri::Smooth<MeshType>::VertexQualityLaplacian(m, 2);

        ForEachVertex(m, [] (VertexType & v) {
            v.Q()  = 1 / (v.Q() + 1);
        });
    }

    static void computeQualityDistFromRadii(MeshType & m)
    {
        tri::RequirePerVertexQuality(m);
        tri::RequirePerFaceQuality(m);

        ScalarType maxV = 0;
        ScalarType minV = 10;

        ForEachFace(m, [&] (FaceType & f) {
            f.Q() = 1. - vcg::QualityRadii(f.cP(0), f.cP(1), f.cP(2));
            maxV = std::max(maxV, f.Q());
            minV = std::min(minV, f.Q());
        });

        //normalize
        ForEachFace(m, [&] (FaceType & f) {
            f.Q() = std::pow((f.Q() - minV) / (maxV - minV), 2.);
        });

        std::vector<ScalarType> vertMax(m.VN(), 0);
        std::vector<ScalarType> vertMin(m.VN(), 10);

        ForEachFace(m, [&] (FaceType & f) {
            for (int i = 0; i < 3; ++i)
            {
                auto vidx = vcg::tri::Index(m, f.V(i));
                vertMax[vidx] = std::max(vertMax[vidx], f.Q());
                vertMin[vidx] = std::min(vertMin[vidx], f.Q());
            }
        });

        for (size_t v = 0; v < m.VN(); ++v)
        {
            m.vert[v].Q() = vertMax[v] - vertMin[v];
        }
       
//        tri::UpdateQuality<MeshType>::VertexFromFace(m);
    }

    static void computeQualityDistFromHeight(MeshType & m, const ScalarType lowerBound, const ScalarType higherBound)
    {
        tri::RequirePerVertexQuality(m);
        tri::RequirePerFaceQuality(m);

        ScalarType maxV = 0;
        ScalarType minV = 10;


        ForEachFace(m, [&] (FaceType & f) {
            ScalarType minH = std::numeric_limits<ScalarType>::max();
            for (int i = 0; i < 3; ++i)
            {
                CoordType closest;
                ScalarType dist = 0;
                vcg::Segment3<ScalarType> seg(f.cP1(i), f.cP2(i));
                vcg::SegmentPointDistance(seg, f.cP0(i), closest, dist);
                minH = std::min(minH, dist);
            }

            f.Q() = math::Clamp(minH, lowerBound, higherBound);
        });
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

    static bool testHausdorff (MeshType & m, StaticGrid & grid, const std::vector<CoordType> & verts, const ScalarType maxD, const CoordType checkOrientation = CoordType(0,0,0))
    {
        for (CoordType v : verts)
        {
            CoordType closest, normal, ip;
            ScalarType dist = 0;
            FaceType* fp = GetClosestFaceBase(m, grid, v, maxD, dist, closest);

            //you can't use this kind of orientation check, since when you stand on edges it fails
            if (fp == NULL || (checkOrientation != CoordType(0,0,0) && checkOrientation * fp->N() < 0.7))
            {
                return false;
            }
        }
        return true;
    }



    /*
                Edge Swap Step:
                This method optimizes the valence of each vertex.
                oldDist is the sum of the absolute distance of each vertex from its ideal valence
                newDist is the sum of the absolute distance of each vertex from its ideal valence after
                the edge swap.
                If the swap decreases the total absolute distance, then it's applied, preserving the triangle
                quality.                        +1
                           v1                     v1
                          /  \                   /|\
                         /    \                 / | \
                        /      \               /  |  \
                       /     _*p\           -1/   |   \ -1
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
        if (/*p.IsBorder() || */p.IsEdgeS()) return false;

        int oldDist = 0, newDist = 0, idealV, actualV;

        PosType tp=p;

        VertexType *v0=tp.V();

        std::vector<VertexType*> incident;

        vcg::face::VVStarVF<FaceType>(tp.V(), incident);
        idealV  = idealValence(tp); actualV = incident.size();
        oldDist += abs(idealV - actualV); newDist += abs(idealV - (actualV - 1));

        tp.NextF();tp.FlipE();tp.FlipV();
        VertexType *v1=tp.V();
        vcg::face::VVStarVF<FaceType>(tp.V(), incident);
        idealV  = idealValence(tp); actualV = incident.size();
        oldDist += abs(idealV - actualV); newDist += abs(idealV - (actualV + 1));

        tp.FlipE();tp.FlipV();tp.FlipE();
        VertexType *v2=tp.V();
        vcg::face::VVStarVF<FaceType>(tp.V(), incident);
        idealV  = idealValence(tp); actualV = incident.size();
        oldDist += abs(idealV - actualV); newDist += abs(idealV - (actualV - 1));

        tp.NextF();tp.FlipE();tp.FlipV();
        VertexType *v3=tp.V();
        vcg::face::VVStarVF<FaceType>(tp.V(), incident);
        idealV  = idealValence(tp); actualV = incident.size();
        oldDist += abs(idealV - actualV); newDist += abs(idealV - (actualV + 1));

        ScalarType qOld = std::min(Quality(v0->P(),v2->P(),v3->P()),Quality(v0->P(),v1->P(),v2->P()));
        ScalarType qNew = std::min(Quality(v0->P(),v1->P(),v3->P()),Quality(v2->P(),v3->P(),v1->P()));

        return (newDist < oldDist && qNew >= qOld * 0.50f) ||
                (newDist == oldDist && qNew > qOld * 1.f) || qNew > 1.5f * qOld;
    }

    static bool checkManifoldness(FaceType & f, int z)
    {
        PosType pos(&f, (z+2)%3, f.V2(z));
        PosType start = pos;

        do {
            pos.FlipE();
            if (!face::IsManifold(*pos.F(), pos.E()))
                break;
            pos.FlipF();
        } while (pos!=start);

        return pos == start;
    }

    // Edge swap step: edges are flipped in order to optimize valence and triangle quality across the mesh
    static void ImproveValence(MeshType &m, Params &params)
    {
        static ScalarType foldCheckRad = math::ToRad(5.);
        tri::UpdateTopology<MeshType>::FaceFace(m);
        tri::UpdateTopology<MeshType>::VertexFace(m);
        ForEachFace(m, [&] (FaceType & f) {
            //			if (face::IsManifold(f, 0) && face::IsManifold(f, 1) && face::IsManifold(f, 2))
            for (int i = 0; i < 3; ++i)
            {
                if (&f > f.cFFp(i))
                {
                    PosType pi(&f, i);
                    CoordType swapEdgeMidPoint = (f.cP2(i) + f.cFFp(i)->cP2(f.cFFi(i))) / 2.;
                    std::vector<CoordType> toCheck(1, swapEdgeMidPoint);


                    if(((!params.selectedOnly) || (f.IsS() && f.cFFp(i)->IsS())) &&
                            !face::IsBorder(f, i) &&
                            face::IsManifold(f, i) && /*checkManifoldness(f, i) &&*/
                            face::checkFlipEdgeNotManifold(f, i) &&
                            testSwap(pi, params.creaseAngleCosThr) &&
//                            face::CheckFlipEdge(f, i) &&
                            face::CheckFlipEdgeNormal(f, i, params.creaseAngleRadThr) && //vcg::math::ToRad(5.)) &&
                            (!params.surfDistCheck || testHausdorff(*params.mProject, params.grid, toCheck, params.maxSurfDist)))
                    {
                        //When doing the swap we need to preserve and update the crease info accordingly
                        FaceType* g = f.cFFp(i);
                        int w = f.FFi(i);

                        bool creaseF = g->IsFaceEdgeS((w + 1) % 3);
                        bool creaseG = f.IsFaceEdgeS((i + 1) % 3);

                        face::FlipEdgeNotManifold(f, i);

                        f.ClearFaceEdgeS((i + 1) % 3);
                        g->ClearFaceEdgeS((w + 1) % 3);

                        if (creaseF)
                            f.SetFaceEdgeS(i);
                        if (creaseG)
                            g->SetFaceEdgeS(w);

                        ++params.stat.flipNum;
                        break;
                    }
                }
            }
        });
    }

    // The predicate that defines which edges should be split
    class EdgeSplitAdaptPred
    {
    public:
        int count = 0;
        ScalarType length, lengthThr, minQ, maxQ;
        const Params & params;

        EdgeSplitAdaptPred(const Params & p) : params(p) {};

        bool operator()(PosType &ep)
        {
            ScalarType quality = (((math::Abs(ep.V()->Q())+math::Abs(ep.VFlip()->Q()))/(ScalarType)2.0)-minQ)/(maxQ-minQ);
            ScalarType mult = computeLengthThrMult(params, quality);
            ScalarType dist = Distance(ep.V()->P(), ep.VFlip()->P());
            if(dist > mult * length)
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
            EdgeSplitAdaptPred ep(params);
            ep.minQ      = minQ;
            ep.maxQ      = maxQ;
            ep.length    = params.maxLength;
            ep.lengthThr = params.lengthThr;
            tri::RefineMidpoint(m, ep, params.selectedOnly);
            params.stat.splitNum+=ep.count;
        }
        else {
            EdgeSplitLenPred ep;
            ep.squaredlengthThr = params.maxLength*params.maxLength;
            tri::RefineMidpoint(m, ep, params.selectedOnly);
            params.stat.splitNum+=ep.count;
        }
    }

    static int VtoE(const int v0, const int v1)
    {
        static /*constexpr*/ int Vmat[3][3] = { -1,  0,  2,
                                                0, -1,  1,
                                                2,  1, -1};
        return Vmat[v0][v1];
    }


    static bool checkCanMoveOnCollapse(PosType p, std::vector<FaceType*> & faces, std::vector<int> & vIdxes, Params &params)
    {
        bool allIncidentFaceSelected = true;

        PosType pi = p;

        CoordType dEdgeVector = (p.V()->cP() - p.VFlip()->cP()).Normalize();

        int incidentFeatures = 0;

        vcg::tri::UnMarkAll(*params.m);

        for (size_t i = 0; i < faces.size(); ++i)
        {
            if (faces[i]->IsFaceEdgeS(VtoE(vIdxes[i], (vIdxes[i]+1)%3)) && !vcg::tri::IsMarked(*params.m, faces[i]->cV1(vIdxes[i])))
            {
				vcg::tri::Mark(*params.m,faces[i]->V1(vIdxes[i]));
                incidentFeatures++;
                CoordType movingEdgeVector0 = (faces[i]->cP1(vIdxes[i]) - faces[i]->cP(vIdxes[i])).Normalize();
                if (std::fabs(movingEdgeVector0 * dEdgeVector) < .9f || !p.IsEdgeS())
                    return false;
            }
            if (faces[i]->IsFaceEdgeS(VtoE(vIdxes[i], (vIdxes[i]+2)%3)) && !vcg::tri::IsMarked(*params.m, faces[i]->cV2(vIdxes[i])))
            {
				vcg::tri::Mark(*params.m,faces[i]->V2(vIdxes[i]));
                incidentFeatures++;
                CoordType movingEdgeVector1 = (faces[i]->cP2(vIdxes[i]) - faces[i]->cP(vIdxes[i])).Normalize();
                if (std::fabs(movingEdgeVector1 * dEdgeVector) < .9f || !p.IsEdgeS())
                    return false;
            }
            allIncidentFaceSelected &= faces[i]->IsS();
        }

        if (incidentFeatures > 2)
            return false;

        return params.selectedOnly ? allIncidentFaceSelected : true;
    }

    static bool checkFacesAfterCollapse (std::vector<FaceType*> & faces, PosType p, const Point3<ScalarType> &mp, Params &params, bool relaxed)
    {
        for (FaceType* f : faces)
        {
            if(!(*f).IsD() && f != p.F()) //i'm not a deleted face
            {
                PosType pi(f, p.V()); //same vertex

                VertexType *v0 = pi.V();
                VertexType *v1 = pi.F()->V1(pi.VInd());
                VertexType *v2 = pi.F()->V2(pi.VInd());

                if( v1 == p.VFlip() || v2 == p.VFlip()) //i'm the other deleted face
                    continue;

                //check on new face quality
                {
                    ScalarType newQ = Quality(mp,      v1->P(), v2->P());
                    ScalarType oldQ = Quality(v0->P(), v1->P(), v2->P());

                    if(newQ <= 0.5*oldQ)
                        return false;
                }

                // we prevent collapse that makes edges too long (except for cross)
                if(!relaxed)
                    if((Distance(mp, v1->P()) > params.maxLength || Distance(mp, v2->P()) > params.maxLength))
                        return false;

                Point3<ScalarType> oldN = NormalizedTriangleNormal(*(pi.F()));
                Point3<ScalarType> newN = Normal(mp, v1->P(), v2->P()).Normalize();

//                if (oldN * newN < 0.5f)
//                    return false;

                std::vector<CoordType> baryP(1);
                baryP[0] = (v1->cP() + v2->cP() + mp) / 3.;

                if (!testHausdorff(*(params.mProject), params.grid, baryP, params.maxSurfDist, newN))
                    return false;

                //check on new face distance from original mesh
                if (params.surfDistCheck)
                {
                    std::vector<CoordType> points(3);
                    std::vector<CoordType> baryP(1);

                    baryP[0] = (v1->cP() + v2->cP() + mp) / 3.;

                    points[0] = (v1->cP() + mp) / 2.;
                    points[1] = (v2->cP() + mp) / 2.;
                    points[2] = mp;

                    if (!testHausdorff(*(params.mProject), params.grid, points, params.maxSurfDist))// ||
//                            !testHausdorff(*(params.mProject), params.grid, baryP, params.maxSurfDist, newN))
                        return false;
                }
            }
        }
        return true;
    }


    //TODO: Refactor code and implement the correct set up of crease info when collapsing towards a crease edge
    static bool checkCollapseFacesAroundVert1(PosType &p, VertexPair & pair, Point3<ScalarType> &mp, Params &params, bool relaxed)
    {
        PosType p0 = p, p1 = p;

        p1.FlipV();

        vector<int> vi0, vi1;
        vector<FaceType*> ff0, ff1;

        face::VFStarVF<FaceType>(p0.V(), ff0, vi0);
        face::VFStarVF<FaceType>(p1.V(), ff1, vi1);

        //check crease-moveability
        bool moveable0 = checkCanMoveOnCollapse(p0, ff0, vi0, params) && !p0.V()->IsS();
        bool moveable1 = checkCanMoveOnCollapse(p1, ff1, vi1, params) && !p1.V()->IsS();

        //if both moveable => go to midpoint
        // else collapse on movable one
        if (!moveable0 && !moveable1)
            return false;

        pair = moveable0 ? VertexPair(p0.V(), p1.V()) : VertexPair(p1.V(), p0.V());

        //casting int(true) is always 1 and int(false) = =0
        assert(int(true) == 1);
        assert(int(false) == 0);
        mp = (p0.V()->cP() * int(moveable1) + p1.V()->cP() * int(moveable0)) / (int(moveable0) + int(moveable1));

        if (checkFacesAfterCollapse(ff0, p0, mp, params, relaxed))
            return checkFacesAfterCollapse(ff1, p1, mp, params, relaxed);

        return false;
    }

    static bool testCollapse1(PosType &p, VertexPair & pair, Point3<ScalarType> &mp, ScalarType minQ, ScalarType maxQ, Params &params, bool relaxed = false)
    {
        ScalarType quality = (((math::Abs(p.V()->Q())+math::Abs(p.VFlip()->Q()))/(ScalarType)2.0)-minQ)/(maxQ-minQ);
        ScalarType mult = computeLengthThrMult(params, quality);
        ScalarType thr = mult*params.minLength;

        ScalarType dist = Distance(p.V()->P(), p.VFlip()->P());
        ScalarType area = DoubleArea(*(p.F()))/2.f;
        if(relaxed || (dist < thr || area < params.minLength*params.minLength/100.f))//if to collapse
        {
            return checkCollapseFacesAroundVert1(p, pair, mp, params, relaxed);
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
            if(!(*v).IsD() && (*v).IsB() && v != p.VFlip()) //ignore non border
                collapsedNV0 = ((*v).P() - p.VFlip()->P()).normalized(); //edge vector after collapse

        face::VVStarVF<FaceType>(p.VFlip(), vv);

        for(VertexType *v: vv)
            if(!(*v).IsD() && (*v).IsB() && v != p.V()) //ignore non border
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

        SelectionStack<MeshType> ss(m);
        ss.push();

        {
            tri::UpdateTopology<MeshType>::FaceFace(m);
            Clean<MeshType>::CountNonManifoldVertexFF(m,true);

            //FROM NOW ON VSelection is NotManifold

            for(auto fi=m.face.begin(); fi!=m.face.end(); ++fi)
                if(!(*fi).IsD() && (params.selectedOnly == false || fi->IsS()))
                {
                    for(auto i=0; i<3; ++i)
                    {
                        PosType pi(&*fi, i);
                        ++candidates;
                        VertexPair  bp = VertexPair(pi.V(), pi.VFlip());
                        Point3<ScalarType> mp = (pi.V()->P()+pi.VFlip()->P())/2.f;

                        if(testCollapse1(pi, bp, mp, minQ, maxQ, params) && Collapser::LinkConditions(bp))
                        {
                            Collapser::Do(m, bp, mp, true);
                            ++params.stat.collapseNum;
                            break;
                        }

                    }
                }
        }
        ss.pop();
    }


    //Here I just need to check the faces of the cross, since the other faces are not
    //affected by the collapse of the internal faces of the cross.
    static bool testCrossCollapse(PosType &p, std::vector<FaceType*> ff, std::vector<int> vi, Point3<ScalarType> &mp, Params &params)
    {
        if(!checkFacesAfterCollapse(ff, p, mp, params, true))
            return false;
        return true;
    }

    /*
         *Choose the best way to collapse a cross based on the (external) cross vertices valence
         *and resulting face quality
         *                                      +0                   -1
         *             v1                    v1                    v1
         *            /| \                   /|\                  / \
         *           / |  \                 / | \                /   \
         *          /  |   \               /  |  \              /     \
         *         / *p|    \           -1/   |   \ -1       +0/       \+0
         *       v0-------- v2 ========> v0   |   v2    OR    v0-------v2
         *        \    |    /             \   |   /            \       /
         *         \   |   /               \  |  /              \     /
         *          \  |  /                 \ | /                \   /
         *           \ | /                   \|/ +0               \ / -1
         *             v3                     v3                   v3
         */
    static bool chooseBestCrossCollapse(PosType &p, VertexPair& bp, vector<FaceType*> &ff)
    {
        vector<VertexType*> vv0, vv1, vv2, vv3;
        VertexType *v0, *v1, *v2, *v3;

        v0 = p.F()->V1(p.VInd());
        v1 = p.F()->V2(p.VInd());


        bool crease[4] = {false, false, false, false};

        crease[0] = p.F()->IsFaceEdgeS(VtoE(p.VInd(), (p.VInd()+1)%3));
        crease[1] = p.F()->IsFaceEdgeS(VtoE(p.VInd(), (p.VInd()+2)%3));

        for(FaceType *f: ff)
            if(!(*f).IsD() && f != p.F())
            {
                PosType pi(f, p.V());
                VertexType *fv1 = pi.F()->V1(pi.VInd());
                VertexType *fv2 = pi.F()->V2(pi.VInd());

                if(fv1 == v0 || fv2 == v0)
                {
                    if (fv1 == 0)
                    {
                        v3 = fv2;
                        crease[3] = f->IsFaceEdgeS(VtoE(pi.VInd(), (pi.VInd()+2)%3));
                    }
                    else
                    {
                        v3 = fv1;
                        crease[3] = f->IsFaceEdgeS(VtoE(pi.VInd(), (pi.VInd()+1)%3));
                    }
                    //					v3 = (fv1 == v0) ? fv2 : fv1;
                }

                if(fv1 == v1 || fv2 == v1)
                {
                    if (fv1 == v1)
                    {
                        v2 = fv2;
                        crease[2] = f->IsFaceEdgeS(VtoE(pi.VInd(), (pi.VInd()+2)%3));
                    }
                    else
                    {
                        v2 = fv1;
                        crease[2] = f->IsFaceEdgeS(VtoE(pi.VInd(), (pi.VInd()+1)%3));
                    }
                    //					v2 = (fv1 == v1) ? fv2 : fv1;
                }
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

        if (crease[0] || crease[1] || crease[2] || crease[3])
            return false;
        //		if (crease[0] && crease[1] && crease[2] && crease[3])
        //		{
        //			return false;
        //		}

        //		if (crease[0] || crease[2])
        //		{
        //			bp = VertexPair(p.V(), v0);
        //			return true;
        //		}

        //		if (crease[1] || crease[3])
        //		{
        //			bp = VertexPair(p.V(), v1);
        //			return true;
        //		}

        //no crease
        if(delta1 < delta2 && Q1 >= 0.6f*Q2)
        {
            bp = VertexPair(p.V(), v1);
            return true;
        }
        else
        {
            bp = VertexPair(p.V(), v0);
            return true;
        }
    }
    //Cross Collapse pass: This pass cleans the mesh from cross vertices, keeping in mind the link conditions
    //and feature preservations tests.
    static void CollapseCrosses(MeshType &m , Params &params)
    {
        tri::UpdateTopology<MeshType>::VertexFace(m);
        tri::UpdateFlags<MeshType>::VertexBorderFromNone(m);
        int count = 0;

        SelectionStack<MeshType> ss(m);
        ss.push();


        {
            tri::UpdateTopology<MeshType>::FaceFace(m);
            Clean<MeshType>::CountNonManifoldVertexFF(m,true);

            //From now on Selection on vertices is not manifoldness

            for(auto fi=m.face.begin(); fi!=m.face.end(); ++fi)
                if(!(*fi).IsD() && (!params.selectedOnly || fi->IsS()))
                {
                    for(auto i=0; i<3; ++i)
                    {
                        PosType pi(&*fi, i);
                        if(!pi.V()->IsB())
                        {
                            vector<FaceType*> ff;
                            vector<int> vi;
                            face::VFStarVF<FaceType>(pi.V(), ff, vi);

                            //if cross need to check what creases you have and decide where to collapse accordingly
                            //if tricuspidis need whenever you have at least one crease => can't collapse anywhere
                            if(ff.size() == 4 || ff.size() == 3)
                            {
                                //							VertexPair bp;
                                VertexPair  bp = VertexPair(pi.V(), pi.VFlip());
                                Point3<ScalarType> mp = (pi.V()->P()+pi.VFlip()->P())/2.f;

                                if(testCollapse1(pi, bp, mp, 0, 0, params, true) && Collapser::LinkConditions(bp))
                                {
                                    Collapser::Do(m, bp, mp, true);
                                    ++params.stat.collapseNum;
                                    ++count;
                                    break;
                                }
                            }
                        }
                    }
                }
        }

        ss.pop();
        Allocator<MeshType>::CompactEveryVector(m);
    }

    // This function sets the selection bit on vertices that lie on creases
    static int selectVertexFromCrease(MeshType &m, ScalarType creaseThr)
    {
        int count = 0;
        Clean<MeshType>::CountNonManifoldVertexFF(m, true, false);

        ForEachFacePos(m, [&](PosType &p){
            if(p.IsBorder() || p.IsEdgeS()/*testCreaseEdge(p, creaseThr)*/)
            {
                p.V()->SetS();
                p.VFlip()->SetS();
                ++count;
            }
        });
        return count;
    }

    static int selectVertexFromFold(MeshType &m, Params & params)
    {
        std::vector<char> creaseVerts(m.VN(), 0);
        ForEachFacePos(m, [&] (PosType & p) {
            if (p.IsEdgeS())
            {
                creaseVerts[vcg::tri::Index(m, p.V())] = 1;
                creaseVerts[vcg::tri::Index(m, p.VFlip())] = 1;
            }
        });


        ForEachFace(m, [&] (FaceType & f) {
            for (int i = 0; i < 3; ++i)
            {
                if (f.FFp(i) > &f)
                {
                    ScalarType angle = fastAngle(NormalizedTriangleNormal(f), NormalizedTriangleNormal(*(f.FFp(i))));
                    if (angle <= params.foldAngleCosThr)
                    {
                        if (creaseVerts[vcg::tri::Index(m, f.V0(i))] == 0)
                            f.V0(i)->SetS();
                        if (creaseVerts[vcg::tri::Index(m, f.V1(i))] == 0)
                            f.V1(i)->SetS();
                        if (creaseVerts[vcg::tri::Index(m, f.V2(i))] == 0)
                            f.V2(i)->SetS();
                        if (creaseVerts[vcg::tri::Index(m, f.FFp(i)->V2(f.FFi(i)))] == 0)
                            f.FFp(i)->V2(f.FFi(i))->SetS();
                    }
                }
            }
        });

        return 0;
    }




    static void FoldRelax(MeshType &m, Params params, const int step, const bool strict = true)
    {
        typename vcg::tri::Smooth<MeshType>::LaplacianInfo lpz(CoordType(0, 0, 0), 0);
        SimpleTempData<typename MeshType::VertContainer, typename vcg::tri::Smooth<MeshType>::LaplacianInfo> TD(m.vert, lpz);
        const ScalarType maxDist = (strict) ? params.maxSurfDist / 1000. : params.maxSurfDist;
        for (int i = 0; i < step; ++i)
        {
            TD.Init(lpz);
            vcg::tri::Smooth<MeshType>::AccumulateLaplacianInfo(m, TD, false);

            for (auto fi = m.face.begin(); fi != m.face.end(); ++fi)
            {
                std::vector<CoordType> newPos(4);
                bool moving = false;

                for (int j = 0; j < 3; ++j)
                {
                    newPos[j] = fi->cP(j);
                    if (!fi->V(j)->IsD() && TD[fi->V(j)].cnt > 0)
                    {
                        if (fi->V(j)->IsS())
                        {
                            newPos[j] = (fi->V(j)->P() + TD[fi->V(j)].sum) / (TD[fi->V(j)].cnt + 1);
                            moving = true;
                        }
                    }
                }

                if (moving)
                {
                    //					const CoordType oldN = vcg::NormalizedTriangleNormal(*fi);
                    //					const CoordType newN = vcg::Normal(newPos[0], newPos[1], newPos[2]).Normalize();

                    newPos[3] = (newPos[0] + newPos[1] + newPos[2]) / 3.;
                    if (/*(strict || oldN * newN > 0.99) &&*/ (!params.surfDistCheck || testHausdorff(*params.mProject, params.grid, newPos, maxDist)))
                    {
                        for (int j = 0; j < 3; ++j)
                            fi->V(j)->P() = newPos[j];
                    }
                }
            }
        }
    }

    static void VertexCoordPlanarLaplacian(MeshType &m, Params & params, int step, ScalarType delta = 0.2)
    {
        typename vcg::tri::Smooth<MeshType>::LaplacianInfo lpz(CoordType(0, 0, 0), 0);
        SimpleTempData<typename MeshType::VertContainer, typename vcg::tri::Smooth<MeshType>::LaplacianInfo> TD(m.vert, lpz);
        for (int i = 0; i < step; ++i)
        {
            TD.Init(lpz);
            vcg::tri::Smooth<MeshType>::AccumulateLaplacianInfo(m, TD, false);
            // First normalize the AccumulateLaplacianInfo
            for (auto vi = m.vert.begin(); vi != m.vert.end(); ++vi)
                if (!(*vi).IsD() && TD[*vi].cnt > 0)
                {
                    if ((*vi).IsS())
                        TD[*vi].sum = ((*vi).P() + TD[*vi].sum) / (TD[*vi].cnt + 1);
                }

//            for (auto fi = m.face.begin(); fi != m.face.end(); ++fi)
//            {
//                if (!(*fi).IsD())
//                {
//                    for (int j = 0; j < 3; ++j)
//                    {
//                        if (Angle(Normal(TD[(*fi).V0(j)].sum, (*fi).P1(j), (*fi).P2(j)),
//                                  Normal((*fi).P0(j), (*fi).P1(j), (*fi).P2(j))) > M_PI/2.)
//                            TD[(*fi).V0(j)].sum = (*fi).P0(j);
//                    }
//                }
//            }
//            for (auto fi = m.face.begin(); fi != m.face.end(); ++fi)
//            {
//                if (!(*fi).IsD())
//                {
//                    for (int j = 0; j < 3; ++j)
//                    {
//                        if (Angle(Normal(TD[(*fi).V0(j)].sum, TD[(*fi).V1(j)].sum, (*fi).P2(j)),
//                                  Normal((*fi).P0(j), (*fi).P1(j), (*fi).P2(j))) > M_PI/2.)
//                        {
//                            TD[(*fi).V0(j)].sum = (*fi).P0(j);
//                            TD[(*fi).V1(j)].sum = (*fi).P1(j);
//                        }
//                    }
//                }
//            }

            for (auto vi = m.vert.begin(); vi != m.vert.end(); ++vi)
                if (!(*vi).IsD() && TD[*vi].cnt > 0)
                {
                    std::vector<CoordType> newPos(1, TD[*vi].sum);
                    if ((*vi).IsS() && testHausdorff(*params.mProject, params.grid, newPos, params.maxSurfDist))
                        (*vi).P() = (*vi).P() * (1-delta) + TD[*vi].sum * (delta);
                }
        } // end step
    }

    //	static int
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

        VertexCoordPlanarLaplacian(m, params, 1);

        tri::UpdateSelection<MeshType>::VertexClear(m);

        selectVertexFromFold(m, params);
        FoldRelax(m, params, 2);

        tri::UpdateSelection<MeshType>::VertexClear(m);

        if(params.selectedOnly) {
            ss.pop();
        }
    }
    /*
                Reprojection step, this method reprojects each vertex on the original surface
                sampling the nearest Point3 onto it using a uniform grid StaticGrid t
        */
    //TODO: improve crease reprojection:
    //		crease verts should reproject only on creases.
    static void ProjectToSurface(MeshType &m, Params & params)
    {
        for(auto vi=m.vert.begin();vi!=m.vert.end();++vi)
            if(!(*vi).IsD())
            {
                Point3<ScalarType> newP, normP, barP;
                ScalarType maxDist = params.maxSurfDist * 2.5f, minDist = 0.f;
                FaceType* fp = GetClosestFaceBase(*params.mProject, params.grid, vi->cP(), maxDist, minDist, newP/*, normP, barP*/);

                if (fp != NULL)
                {
                    vi->P() = newP;
                }
            }
    }
};
} // end namespace tri
} // end namespace vcg
#endif
