/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2022                                           \/)\/    *
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

#ifndef __VCGLIB_REFINE_DOOSABIN_H
#define __VCGLIB_REFINE_DOOSABIN_H
namespace vcg {
namespace tri {
/// \ingroup trimesh

/// \headerfile refine_doosabin.h vcg/complex/algorithms/refine_doosabin.h

/// \brief This class is used convert between polygonal meshes and triangular meshes

/**
    This class contains two members that allow to build a triangular mesh from a polygonal mesh
    and viceversa. In a trimesh, the generic polygons with n sides are codified represented by
    tagging the internal edge of the face as 'faux' with the SetF.
    */

template <class PolyMeshType>
class DooSabin {
    typedef typename PolyMeshType::FaceType FaceType;
    typedef typename PolyMeshType::FacePointer FacePointer;
    typedef typename PolyMeshType::FaceIterator FaceIterator;
    typedef typename PolyMeshType::VertexIterator VertexIterator;
    
public:
static void Refine(PolyMeshType &baseIn, PolyMeshType &refinedOut, int iterationNum=1)
{
	tri::RequirePolygonalMesh(baseIn);
	tri::RequirePolygonalMesh(refinedOut);
	tri::RequireFFAdjacency(baseIn);
//	fprintf(stdout,"Refining starting \n");fflush(stdout);
	
    PolyMeshType refined;
    PolyMeshType base;
    Append<PolyMeshType,PolyMeshType>::MeshCopy(base,baseIn);
    for(int step = 0; step<iterationNum;++step)
    {
        refined.Clear();
//		fprintf(stdout,"Refining iteration %i mesh has %i faces\n",step,base.FN() );fflush(stdout);
        UpdateTopology<PolyMeshType>::FaceFace(base);
        //    tri::UpdateTopology<PolyMeshType>::VertexFace(base);
        typedef typename PolyMeshType :: template PerVertexAttributeHandle<int> VertexIntHandleType;
        typedef typename PolyMeshType :: template PerFaceAttributeHandle<int> FaceIntHandleType;
        typedef typename PolyMeshType :: template PerVertexAttributeHandle<std::pair<FacePointer,int> > VertexPairHandleType;
        
        // This handle stores for each vertex v the index of the face that I have added to the refined mesh
        VertexIntHandleType FaceVertIndVH = vcg::tri::Allocator<PolyMeshType>:: template GetPerVertexAttribute<int>(base,"facevertInd");
        
        // This handle stores for each face f the index of the first of the f.VN() vertexes that I have added for that face       
        FaceIntHandleType FaceVertBaseFH = vcg::tri::Allocator<PolyMeshType>:: template GetPerFaceAttribute<int>(base,"FaceVertBase");
        
        // Computing the degree of each vertex and storing it in an attribute
        VertexIntHandleType degreeVH = vcg::tri::Allocator<PolyMeshType>:: template GetPerVertexAttribute<int>(base,"degree");
        VertexPairHandleType VFpH = vcg::tri::Allocator<PolyMeshType>:: template GetPerVertexAttribute<std::pair<FacePointer,int> >(base,"VFP");
        
        for(auto vi=base.vert.begin();vi!=base.vert.end();vi++)
            degreeVH[*vi]=0;
        for(auto fi=base.face.begin();fi!=base.face.end();fi++)
        {
            for(int j=0;j<fi->VN();++j)
            {
                degreeVH[fi->V(j)]= degreeVH[fi->V(j)]+1;
                VFpH[fi->V(j)] = std::make_pair(&*fi,j);
            }
        }
        
        // This map indicates for each corner of each face of the base mesh,
        // what is the index of the vertex I have created. 
        std::map<std::pair<int,int>,int> faceCornerToNewVertMap;
        
        // First create a new face for each face of the base mesh
        for(auto fi=base.face.begin();fi!=base.face.end();fi++)
        {
            Point3f b = PolyBarycenter(*fi);
            auto newf = tri::Allocator<PolyMeshType>::AddFaces(refined,1);
            newf->Alloc(fi->VN());
            FaceVertBaseFH[fi]= refined.vert.size();
            for(int j=0;j<fi->VN();++j)
            {
                auto newv = tri::Allocator<PolyMeshType>::AddVertex(refined,(fi->V(j)->P()+b)/2.0f);
                newf->V(j)=&*newv;
                faceCornerToNewVertMap[std::make_pair(tri::Index(base,*fi),j)] = tri::Index(refined,*newv);
            }
        }
        // second loop creating a face for each vertex 
        for(auto vi=base.vert.begin();vi!=base.vert.end();vi++)
        {
            auto newf = tri::Allocator<PolyMeshType>::AddFaces(refined,1);
            newf->Alloc(degreeVH[vi]);
            FaceVertIndVH[tri::Index(refined, *newf)];
            FacePointer curf = VFpH[&*vi].first;
            int curi = VFpH[&*vi].second;
            face::Pos<FaceType> startPos(curf,curi);
            assert(curf->V(curi) == &*vi);
            std::vector<face::Pos<FaceType> > starPosVec; 
            face::VFOrderedStarFF(startPos,starPosVec,false);
            assert(starPosVec.size() == (size_t)degreeVH[vi]);
            for(size_t i =0 ; i < starPosVec.size(); ++i)
            {
                int vind = starPosVec[i].VInd();
                int fpind = tri::Index(base, starPosVec[i].F());
                auto newvi = faceCornerToNewVertMap[std::make_pair(fpind,vind)];
                newf->V(i) = &refined.vert[newvi];
            }              
        }
        
        // Third loop creating the faces on the edges
        tri::UpdateFlags<PolyMeshType>::FaceClearV(base);
        for(auto fi=base.face.begin();fi!=base.face.end();fi++)
        {
            fi->SetV();
            for(int j=0;j<fi->VN();++j)
            {
                if(fi->FFp(j)->IsV())
                {
                    auto newf = tri::Allocator<PolyMeshType>::AddFaces(refined,1);
                    int newvi;
                    newf->Alloc(4);
                    face::Pos<FaceType> startPos(&*fi,j);
                    newvi = faceCornerToNewVertMap[std::make_pair((int)(tri::Index(base, startPos.F())),startPos.VInd())];
                    newf->V(3) = &refined.vert[newvi];
                    startPos.FlipV();
                    newvi = faceCornerToNewVertMap[std::make_pair((int)(tri::Index(base, startPos.F())),startPos.VInd())];
                    newf->V(2) = &refined.vert[newvi];
                    startPos.FlipF();
                    newvi = faceCornerToNewVertMap[std::make_pair((int)(tri::Index(base, startPos.F())),startPos.VInd())];
                    newf->V(1) = &refined.vert[newvi];
                    startPos.FlipV();
                    newvi = faceCornerToNewVertMap[std::make_pair((int)(tri::Index(base, startPos.F())),startPos.VInd())];
                    newf->V(0) = &refined.vert[newvi];
                }
            }
            
        }
        Append<PolyMeshType,PolyMeshType>::MeshCopy(base,refined);        
    }
    Append<PolyMeshType,PolyMeshType>::MeshCopy(refinedOut,refined);    
} // end refine Function

}; // end  DooSabin class 
} // end namespace tri
} // end namespace vcg


#endif // __VCGLIB_REFINE_DOOSABIN_H
