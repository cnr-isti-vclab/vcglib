/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
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
#ifndef __VCGLIB_DUAL_MESH
#define __VCGLIB_DUAL_MESH

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/simplex/face/pos.h>
#include <vcg/complex/algorithms/clean.h>

namespace vcg {
namespace tri {

template <class PolyMeshType>
class DualMeshing
{

    typedef typename PolyMeshType::VertexType VertexType;
    typedef typename PolyMeshType::FaceType   FaceType;
    typedef typename PolyMeshType::CoordType  CoordType;
    typedef typename PolyMeshType::ScalarType ScalarType;

    static void ComposeFace(VertexType &startV,
                            FaceType &startF,
                            std::map<std::pair<int,int>, int> &EdgeMap,
                            const std::vector<int> &FaceMap,
                            const PolyMeshType &primal,
                            std::vector<int> &vertSeq)
    {
        vcg::face::Pos<FaceType> startP(&startF,&startV);

        //get the star of pos
        std::vector<vcg::face::Pos<FaceType> > posVec;
        vcg::face::VFOrderedStarFF(startP,posVec);

        for (size_t i=0;i<posVec.size();i++)
        {
            FaceType *f=posVec[i].F();
            int indexF=vcg::tri::Index(primal,f);
            int indexV=FaceMap[indexF];
            vertSeq.push_back(indexV);
        }

        if (startV.IsB())
        {
            vcg::face::Pos<FaceType> firstPos=posVec[0];
            firstPos.FlipE();
            assert(firstPos.IsBorder());

            int indexVt0=vcg::tri::Index(primal,firstPos.V());
            int indexVt1=vcg::tri::Index(primal,firstPos.VFlip());
            std::pair<int,int> key(std::min(indexVt0,indexVt1),
                                   std::max(indexVt0,indexVt1));
            assert(EdgeMap.count(key)>0);
            int indexV0=EdgeMap[key];

            vcg::face::Pos<FaceType> lastPos=posVec.back();
            assert(lastPos.IsBorder());

            indexVt0=vcg::tri::Index(primal,lastPos.V());
            indexVt1=vcg::tri::Index(primal,lastPos.VFlip());
            key=std::pair<int,int> (std::min(indexVt0,indexVt1),
                                    std::max(indexVt0,indexVt1));
            assert(EdgeMap.count(key)>0);
            int indexV1=EdgeMap[key];

            vertSeq.push_back(indexV1);
            vertSeq.push_back(indexV0);
        }
    }

    static void CreateBorderEdgeVert(PolyMeshType &primal,
                                     PolyMeshType &dual,
                                     std::map<std::pair<int,int>, int> &VertMap)
    {
        VertMap.clear();
        vcg::tri::UpdateFlags<PolyMeshType>::VertexClearB(primal);

        for (size_t i=0;i<primal.face.size();i++)
            for (int j=0;j<primal.face[i].VN();j++)
            {
                int edge_size=primal.face[i].VN();
                FaceType *nextF=primal.face[i].cFFp(j);

                if (nextF!=&primal.face[i])continue;

                VertexType *v0=primal.face[i].V(j);
                VertexType *v1=primal.face[i].V((j+1)%edge_size);

                v0->SetB();
                v1->SetB();

                int V0Index=vcg::tri::Index(primal,v0);
                int V1Index=vcg::tri::Index(primal,v1);
                CoordType pos=(v0->P()+v1->P())/2;
                vcg::tri::Allocator<PolyMeshType>::AddVertex(dual,pos);
                std::pair<int,int> key(std::min(V0Index,V1Index),
                                       std::max(V0Index,V1Index));

                VertMap[key]=dual.vert.size()-1;
            }
    }

    static void CreateFaceVert(PolyMeshType &primal,
                               PolyMeshType &dual,
                               std::vector<int> &VertMap,
                               std::vector<int> &VertFace,
                               bool snapBorder)
    {
        VertMap.clear();
        VertMap.resize(primal.face.size(),-1);
        VertFace.clear();
        VertFace.resize(primal.vert.size(),-1);
        for (size_t i=0;i<primal.face.size();i++)
        {
            CoordType pos(0,0,0);
            int num=0;
            if (snapBorder)//search for border edge
            {
                std::vector<CoordType> BorderPos;
                for (int j=0;j<primal.face[i].VN();j++)
                {
                   if (!primal.face[i].V(j)->IsB())continue;
                    pos+=primal.face[i].P(j);
                    num++;
                }
                if (num>0)
                pos/=num;
            }

            if (num==0)
            {
                for (int j=0;j<primal.face[i].VN();j++)
                {
                    pos+=primal.face[i].V(j)->P();
                    int indexV=vcg::tri::Index(primal,primal.face[i].V(j));
                    if (VertFace[indexV]!=-1)continue;
                    VertFace[indexV]=i;
                }
                pos/=(ScalarType)primal.face[i].VN();
            }
            vcg::tri::Allocator<PolyMeshType>::AddVertex(dual,pos);
            VertMap[i]=dual.vert.size()-1;
        }
    }

public:

    static void MakeDual(PolyMeshType &primal,
                         PolyMeshType &dual,
                         bool snapBorder=true)
    {
        dual.Clear();

        vcg::tri::RequirePolygonalMesh(primal);
        vcg::tri::RequirePolygonalMesh(dual);
        vcg::tri::RequireFFAdjacency(primal);

        vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(primal);

        std::map<std::pair<int,int>, int> VertEdgeMap;
        CreateBorderEdgeVert(primal,dual,VertEdgeMap);

        std::vector<int> VertFaceMap,VertFace;
        CreateFaceVert(primal,dual,VertFaceMap,VertFace,snapBorder);

        for (size_t i=0;i<primal.vert.size();i++)
        {
            if ((snapBorder)&&(primal.vert[i].IsB()))continue;

            FaceType *firstF=&primal.face[VertFace[i]];
            std::vector<int> VertSeq;
            ComposeFace(primal.vert[i],*firstF,VertEdgeMap,VertFaceMap,primal,VertSeq);
            vcg::tri::Allocator<PolyMeshType>::AddFaces(dual,1);
            dual.face.back().Alloc(VertSeq.size());
            for (size_t j=0;j<VertSeq.size();j++)
                dual.face.back().V(j)=&dual.vert[VertSeq[j]];
        }

        vcg::tri::Clean<PolyMeshType>::RemoveUnreferencedVertex(dual);
    }

};

}
}
#endif // __VCGLIB_DUAL_MESH
