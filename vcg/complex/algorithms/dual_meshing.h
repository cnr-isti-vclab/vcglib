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
#include <vcg/complex/algorithms/polygonal_algorithms.h>

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
        (void)EdgeMap;

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
                }
                pos/=(ScalarType)primal.face[i].VN();
            }
            vcg::tri::Allocator<PolyMeshType>::AddVertex(dual,pos);
            VertMap[i]=dual.vert.size()-1;
            dual.vert.back().Q()=i;
        }

        //then initialize VF first face
        for (size_t i=0;i<primal.face.size();i++)
            for (int j=0;j<primal.face[i].VN();j++)
            {
                int indexV=vcg::tri::Index(primal,primal.face[i].V(j));
                if (VertFace[indexV]!=-1)continue;
                VertFace[indexV]=i;
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
        vcg::tri::RequirePerFaceQuality(dual);

        vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(primal);
        vcg::tri::UpdateFlags<PolyMeshType>::VertexBorderFromFaceAdj(primal);

        std::map<std::pair<int,int>, int> VertEdgeMap;

        std::cout<<"Creating Dual Vertices"<<std::endl;
        std::vector<int> VertFaceMap,VertFace;
        CreateFaceVert(primal,dual,VertFaceMap,VertFace,snapBorder);

        std::cout<<"Creating Dual Faces"<<std::endl;
        for (size_t i=0;i<primal.vert.size();i++)
        {
            if (primal.vert[i].IsB())continue;

            FaceType *firstF=&primal.face[VertFace[i]];
            std::vector<int> VertSeq;
            ComposeFace(primal.vert[i],*firstF,VertEdgeMap,VertFaceMap,primal,VertSeq);
            vcg::tri::Allocator<PolyMeshType>::AddFaces(dual,1);
            dual.face.back().Alloc(VertSeq.size());
            for (size_t j=0;j<VertSeq.size();j++)
                dual.face.back().V(j)=&dual.vert[VertSeq[j]];
            dual.face.back().Q()=i;
        }
        vcg::tri::Clean<PolyMeshType>::RemoveUnreferencedVertex(dual);
        //finally remove valence 1 vertices on the border
        vcg::PolygonalAlgorithm<PolyMeshType>::RemoveValence2Vertices(dual);
    }

};

}
}
#endif // __VCGLIB_DUAL_MESH
