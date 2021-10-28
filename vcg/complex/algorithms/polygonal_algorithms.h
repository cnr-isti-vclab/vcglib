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
#ifndef __VCGLIB_POLY_MESH_ALGORITHM
#define __VCGLIB_POLY_MESH_ALGORITHM

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/space/polygon3.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/closest.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <wrap/io_trimesh/export_obj.h>

//define a temporary triangle mesh type
class TempFace;
class TempVertex;

struct TempUsedTypes: public vcg::UsedTypes<vcg::Use<TempVertex>::AsVertexType,
        vcg::Use<TempFace>::AsFaceType>{};

class TempVertex:public vcg::Vertex<TempUsedTypes,
        vcg::vertex::Coord3d,
        vcg::vertex::Normal3d,
        vcg::vertex::BitFlags>
{};

class TempFace:public vcg::Face<TempUsedTypes,
        vcg::face::VertexRef,
        vcg::face::BitFlags,
        vcg::face::FFAdj,
        vcg::face::Mark,
        vcg::face::Normal3d>
{};


class TempMesh: public vcg::tri::TriMesh< std::vector<TempVertex>,std::vector<TempFace > >
{};

namespace vcg{

/*!
\ingroup PolyMeshType

\headerfile color.h vcg/complex/algorithms/polygonal_algorithms.h

\brief processing and optimization of generic polygonal meshes.

This class is used to performs varisous kind of geometric optimization on generic polygonal mesh such as flattengin or imptove the shape of polygons.
*/

template <class PolyMeshType>
class PolygonalAlgorithm
{
    typedef typename PolyMeshType::FaceType      FaceType;
    typedef typename PolyMeshType::VertexType    VertexType;
    typedef typename PolyMeshType::VertexPointer VertexPointer;
    typedef typename PolyMeshType::CoordType     CoordType;
    typedef typename PolyMeshType::ScalarType    ScalarType;
    typedef typename vcg::face::Pos<FaceType>    PosType;

    static void SetFacePos(PolyMeshType &poly_m,
                           int IndexF,std::vector<CoordType> &Pos)
    {
        poly_m.face[IndexF].Dealloc();
        poly_m.face[IndexF].Alloc(Pos.size());
        //std::cout<<Pos.size()<<std::endl;
        int sizeV=poly_m.vert.size();
        for (size_t i=0;i<Pos.size();i++)
            vcg::tri::Allocator<PolyMeshType>::AddVertex(poly_m,Pos[i]);

        for (size_t i=0;i<Pos.size();i++)
            poly_m.face[IndexF].V(i)=&poly_m.vert[sizeV+i];
    }

public:

    static void SubdivideStep(PolyMeshType &poly_m)
    {
        //get the barycenters
        std::vector<CoordType> Bary;
        for (size_t i=0;i<poly_m.face.size();i++)
        {
            CoordType bary(0,0,0);
            for (size_t j=0;j<poly_m.face[i].VN();j++)
                bary+=poly_m.face[i].P(j);

            bary/=poly_m.face[i].VN();
            Bary.push_back(bary);
        }

        //get center of edge
        std::map<std::pair<CoordType,CoordType>, CoordType> EdgeVert;
        for (size_t i=0;i<poly_m.face.size();i++)
            for (size_t j=0;j<poly_m.face[i].VN();j++)
            {
                CoordType Pos0=poly_m.face[i].P0(j);
                CoordType Pos1=poly_m.face[i].P1(j);
                CoordType Avg=(Pos0+Pos1)/2;
                std::pair<CoordType,CoordType> Key(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
                EdgeVert[Key]=Avg;
            }

        int sizeF=poly_m.face.size();
        for (size_t i=0;i<sizeF;i++)
        {
            //retrieve the sequence of pos
            std::vector<CoordType> Pos;
            for (size_t j=0;j<poly_m.face[i].VN();j++)
            {
                CoordType Pos0=poly_m.face[i].P0(j);
                CoordType Pos1=poly_m.face[i].P1(j);
                std::pair<CoordType,CoordType> Key0(std::min(Pos0,Pos1),std::max(Pos0,Pos1));
                Pos0=EdgeVert[Key0];
                Pos.push_back(Pos0);
                Pos.push_back(Pos1);
            }
            //get also the barycenter
            CoordType BaryP=Bary[i];

            //then retrieve the face
            std::vector<CoordType> PosQ;
            PosQ.push_back(Pos[0]);
            PosQ.push_back(Pos[1]);
            PosQ.push_back(Pos[2]);
            PosQ.push_back(BaryP);
            SetFacePos(poly_m,i,PosQ);

            int sizeV=Pos.size();
            //int start=0;
            for (size_t j=2;j<sizeV;j+=2)
            {
                vcg::tri::Allocator<PolyMeshType>::AddFaces(poly_m,1);
                std::vector<CoordType> PosQ;
                PosQ.push_back(Pos[(j)%Pos.size()]);
                PosQ.push_back(Pos[(j+1)%Pos.size()]);
                PosQ.push_back(Pos[(j+2)%Pos.size()]);
                PosQ.push_back(BaryP);
                //start+=2;
                SetFacePos(poly_m,poly_m.face.size()-1,PosQ);
                //break;
            }
        }
        vcg::tri::Clean<PolyMeshType>::RemoveDuplicateVertex(poly_m);
        vcg::tri::Allocator<PolyMeshType>::CompactEveryVector(poly_m);
    }

    static bool CollapseEdges(PolyMeshType &poly_m,
                              const std::vector<PosType> &CollapsePos,
                              const std::vector<CoordType> &InterpPos)
    {

        //this set how to remap the vertices after deletion
        std::map<VertexType*,VertexType*> VertexRemap;

        vcg::tri::UpdateFlags<PolyMeshType>::VertexClearS(poly_m);

        bool collapsed=false;
        //go over all faces and check the ones needed to be deleted
        for (size_t i=0;i<CollapsePos.size();i++)
        {
            FaceType *currF=CollapsePos[i].F();
            int IndexE=CollapsePos[i].E();
            size_t NumV=currF->VN();
            VertexType *v0=currF->V(IndexE);
            VertexType *v1=currF->V((IndexE+1)%NumV);

            //safety check
            assert(v0!=v1);

            if (v0->IsS())continue;
            if (v1->IsS())continue;

            //put on the same position
            v0->P()=InterpPos[i];
            v1->P()=InterpPos[i];

            //select the the two vertices
            v0->SetS();
            v1->SetS();

            //set the remap
            VertexRemap[v1]=v0;

            collapsed=true;
        }

        //then remap vertices
        for (size_t i=0;i<poly_m.face.size();i++)
        {
            int NumV=poly_m.face[i].VN();
            for (int j=0;j<NumV;j++)
            {
                //get the two vertices of the edge
                VertexType *v0=poly_m.face[i].V(j);
                //see if it must substituted or not
                if (VertexRemap.count(v0)==0)continue;
                //in that case remap to the new one
                VertexType *newV=VertexRemap[v0];
                //assign new vertex
                poly_m.face[i].V(j)=newV;
            }
        }

        //then re-elaborate the face
        for (size_t i=0;i<poly_m.face.size();i++)
        {
            //get vertices of the face
            int NumV=poly_m.face[i].VN();
            std::vector<VertexType*> FaceV;
            for (int j=0;j<NumV;j++)
            {
                VertexType *v0=poly_m.face[i].V(j);
                VertexType *v1=poly_m.face[i].V((j+1)%NumV);
                if(v0==v1)continue;
                FaceV.push_back(v0);
            }

            //then deallocate face
            if ((int)FaceV.size()==NumV)continue;

            //otherwise deallocate and set new vertices
            poly_m.face[i].Dealloc();
            poly_m.face[i].Alloc(FaceV.size());
            for (size_t j=0;j<FaceV.size();j++)
                poly_m.face[i].V(j)=FaceV[j];
        }

        //remove unreferenced vertices
        vcg::tri::Clean<PolyMeshType>::RemoveUnreferencedVertex(poly_m);

        //and compact them
        vcg::tri::Allocator<PolyMeshType>::CompactEveryVector(poly_m);

        return collapsed;
    }
private:


    static bool CollapseBorderSmallEdgesStep(PolyMeshType &poly_m,
                                             const ScalarType edge_limit)
    {
        //update topology
        vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(poly_m);

        //update border vertices
        vcg::tri::UpdateFlags<PolyMeshType>::VertexBorderFromFaceAdj(poly_m);


        vcg::tri::UpdateSelection<PolyMeshType>::VertexCornerBorder(poly_m,math::ToRad(150.0));

        std::vector<PosType> CollapsePos;
        std::vector<CoordType> InterpPos;

        //go over all faces and check the ones needed to be deleted
        for (size_t i=0;i<poly_m.face.size();i++)
        {
            int NumV=poly_m.face[i].VN();
            for (int j=0;j<NumV;j++)
            {
                VertexType *v0=poly_m.face[i].V(j);
                VertexType *v1=poly_m.face[i].V((j+1)%NumV);
                assert(v0!=v1);

                bool IsBV0=v0->IsB();
                bool IsBV1=v1->IsB();

                bool IsS0=v0->IsS();
                bool IsS1=v1->IsS();

                if ((IsS0)&&(IsS1))continue;

                //in these cases is not possible to collapse
                if ((!IsBV0)&&(!IsBV1))continue;
                bool IsBorderE=(poly_m.face[i].FFp(j)==&poly_m.face[i]);
                if ((!IsBorderE)&&(IsBV0)&&(IsBV1))continue;

                assert((IsBV0)||(IsBV1));
                CoordType pos0=v0->P();
                CoordType pos1=v1->P();
                ScalarType currL=(pos0-pos1).Norm();
                if (currL>edge_limit)continue;

                //then collapse the point
                CoordType CurrInterpPos;
                if ((IsBV0)&&(!IsBV1))CurrInterpPos=pos0;
                if ((!IsBV0)&&(IsBV1))CurrInterpPos=pos1;
                if ((IsBV0)&&(IsBV1))
                {
                    if ((!IsS0)&&(!IsS1))
                        CurrInterpPos=(pos0+pos1)/2.0;
                    else
                    {
                        if ((!IsS0)&&(IsS1))
                            CurrInterpPos=pos1;
                        else
                        {
                            assert((IsS0)&&(!IsS1));
                            CurrInterpPos=pos0;
                        }
                    }
                }
                CollapsePos.push_back(PosType(&poly_m.face[i],j));
                InterpPos.push_back(CurrInterpPos);
            }
        }

        return CollapseEdges(poly_m,CollapsePos,InterpPos);
    }

    static void LaplacianPos(PolyMeshType &poly_m,std::vector<CoordType> &AvVert)
    {
        //cumulate step
        AvVert.clear();
        AvVert.resize(poly_m.vert.size(),CoordType(0,0,0));
        std::vector<ScalarType> AvSum(poly_m.vert.size(),0);
        for (size_t i=0;i<poly_m.face.size();i++) {
            if (poly_m.face[i].IsD())
                continue;

            for (size_t j=0;j<(size_t)poly_m.face[i].VN();j++)
            {
                //get current vertex
                VertexType *currV=poly_m.face[i].V(j);
                //and its position
                CoordType currP=currV->P();
                //cumulate over other positions
                ScalarType W=vcg::PolyArea(poly_m.face[i]);
                //assert(W!=0);
                for (size_t k=0;k<(size_t)poly_m.face[i].VN();k++)
                {
                    if (k==j) continue;
                    int IndexV=vcg::tri::Index(poly_m,poly_m.face[i].V(k));
                    AvVert[IndexV]+=currP*W;
                    AvSum[IndexV]+=W;
                }
            }
        }

        //average step
        for (size_t i=0;i<poly_m.vert.size();i++)
        {
            if (poly_m.vert[i].IsD())
                continue;

            if (AvSum[i]==0)continue;
            AvVert[i]/=AvSum[i];
        }
    }




    static void UpdateNormal(FaceType &F)
    {
        F.N()=vcg::PolygonNormal(F);
    }

    static void UpdateNormalByFitting(FaceType &F)
    {
        UpdateNormal(F);
        vcg::Plane3<ScalarType> PlF;
        PlF=PolyFittingPlane(F);
        if ((PlF.Direction()*F.N())<0)
            F.N()=-PlF.Direction();
        else
            F.N()=PlF.Direction();
    }

    static void DisplaceBySelected(FaceType &f,std::vector<CoordType> &TemplatePos,
                                   bool FixS,bool FixB)
    {
        CoordType AvPosF(0,0,0);
        CoordType AvPosT(0,0,0);
        size_t Num=0;
        for (size_t i=0;i<f.VN();i++)
        {
            bool AddVal=false;
            AddVal|=((FixS)&&(f.V(i)->IsS()));
            AddVal|=((FixB)&&(f.V(i)->IsB()));
            if (!AddVal)continue;
            Num++;
            AvPosF+=f.V(i)->P();
            AvPosT+=TemplatePos[i];
        }
        if (Num==0)return;
        AvPosF/=(ScalarType)Num;
        AvPosT/=(ScalarType)Num;
        CoordType Displ=AvPosF-AvPosT;
        for (size_t i=0;i<TemplatePos.size();i++)
            TemplatePos[i]+=Displ;
    }

public:

    static void SelectIrregularInternal(PolyMeshType &poly_m)
    {
        vcg::tri::UpdateQuality<PolyMeshType>::VertexValence(poly_m);
        vcg::tri::UpdateSelection<PolyMeshType>::VertexClear(poly_m);
        for (size_t i=0;i<poly_m.vert.size();i++)
        {
            if (poly_m.vert[i].IsB())continue;
            if (poly_m.vert[i].Q()==4)continue;
            poly_m.vert[i].SetS();
        }
    }

    static void SelectIrregularBorder(PolyMeshType &poly_m)
    {
        vcg::tri::UpdateQuality<PolyMeshType>::VertexValence(poly_m);
        for (size_t i=0;i<poly_m.vert.size();i++)
        {
            if (!poly_m.vert[i].IsB())continue;
            if (poly_m.vert[i].Q()==2)continue;
            poly_m.vert[i].SetS();
        }
    }

    static CoordType GetFaceGetBary(FaceType &F)
    {
        CoordType bary=PolyBarycenter(F);
        return bary;
    }

    /*! \brief update the face normal by averaging among vertex's
     * normals computed between adjacent edges
    */
    static void UpdateFaceNormals(PolyMeshType &poly_m)
    {
        for (size_t i=0;i<poly_m.face.size();i++)
            if (!poly_m.face[i].IsD())
            UpdateNormal(poly_m.face[i]);
    }

    /*! \brief update the face normal by fitting a plane
    */
    static void UpdateFaceNormalByFitting(PolyMeshType &poly_m)
    {
        for (size_t i=0;i<poly_m.face.size();i++)
            if (!poly_m.face[i].IsD())
            UpdateNormalByFitting(poly_m.face[i]);
    }


    enum PolyQualityType{QAngle,QPlanar,QTemplate};

    /*! \brief update the quality of the faces by considering different possibilities
     * QAngle   = consider the angle deviation from ideal one (ex 90° quad, 60° triangle...)
     * QPlanar  = consider the difference wrt interpolating plane
     * QTemplate= consider the difference wrt template polygon as in "Statics Aware Grid Shells"
    */
    static void UpdateQuality(PolyMeshType &poly_m,
                              const PolyQualityType &QType)
    {
        for (size_t i=0;i<poly_m.face.size();i++)
        {
            if (poly_m.face[i].IsD())continue;
            switch (QType)
            {
            case QAngle:
                ScalarType AvgDev,WorstDev;
                vcg::PolyAngleDeviation(poly_m.face[i],AvgDev,WorstDev);
                poly_m.face[i].Q()=WorstDev;
                break;
            case QPlanar:
                poly_m.face[i].Q()=vcg::PolyFlatness(poly_m.face[i]);
                break;
            default:
                poly_m.face[i].Q()=vcg::PolyAspectRatio(poly_m.face[i],true);
                break;
            }
        }
    }

    /*! \brief given a face this function returns the template positions as in "Statics Aware Grid Shells"
    */
    static void GetRotatedTemplatePos(FaceType &f,
                                      std::vector<CoordType> &TemplatePos)
    {
        vcg::GetPolyTemplatePos(f,TemplatePos,true);

        CoordType NormT=Normal(TemplatePos);

        //get the normal of vertices
        //CoordType AVN(0,0,0);
        //CoordType AVN0(0,0,0);
        CoordType Origin(0,0,0);
        //        for (int j=0;j<f.VN();j++)
        //            AVN0=AVN0+f.V(j)->N();

        CoordType AVN=vcg::PolygonNormal(f);
        //AVN0.Normalize();
        //        std::cout<<"AVN "<<AVN.X()<<","<<AVN.Y()<<","<<AVN.Z()<<std::endl;
        //        std::cout<<"AVN0 "<<AVN0.X()<<","<<AVN0.Y()<<","<<AVN0.Z()<<std::endl;
        //        std::cout<<"NormT "<<NormT.X()<<","<<NormT.Y()<<","<<NormT.Z()<<std::endl;

        for (size_t j=0;j<TemplatePos.size();j++)
            Origin+=TemplatePos[j];

        Origin/=(ScalarType)TemplatePos.size();
        AVN.Normalize();

        //find rotation matrix
        vcg::Matrix33<ScalarType> Rot=vcg::RotationMatrix(NormT,AVN);

        //apply transformation
        for (size_t j=0;j<TemplatePos.size();j++)
        {
            TemplatePos[j]=TemplatePos[j]-Origin;
            TemplatePos[j]=Rot*TemplatePos[j];
            TemplatePos[j]=TemplatePos[j]+Origin;
        }
    }

    /*! \brief This function performs the polygon regularization as in "Statics Aware Grid Shells"
    */
    static void SmoothPCA(PolyMeshType &poly_m,
                          int relax_step=10,
                          ScalarType Damp=0.5,
                          bool FixS=false,
                          bool isotropic=true,
                          ScalarType smoothTerm=0.1,
                          bool fixB=true,
                          bool WeightByQuality=false,
                          const std::vector<bool> *IgnoreF=NULL)
    {
        (void)isotropic;
        typedef typename PolyMeshType::FaceType PolygonType;
        //        // select irregular ones
        //        if (fixIrr)
        //            poly_m.NumIrregular(true);
        // compute the average edge
        ScalarType MeshArea=0;
        for (size_t i=0;i<poly_m.face.size();i++)
            MeshArea+=vcg::PolyArea(poly_m.face[i]);

        ScalarType AvgArea=MeshArea/(ScalarType)poly_m.face.size();

        if (WeightByQuality)
            UpdateQuality(poly_m,QTemplate);

        if (IgnoreF!=NULL){assert((*IgnoreF).size()==poly_m.face.size());}
        for (size_t s=0;s<(size_t)relax_step;s++)
        {
            //initialize the accumulation vector
            std::vector<CoordType> avgPos(poly_m.vert.size(),CoordType(0,0,0));
            std::vector<ScalarType> weightSum(poly_m.vert.size(),0);
            //then compute the templated positions

            for (size_t i=0;i<poly_m.face.size();i++)
            {
                if ((IgnoreF!=NULL)&&((*IgnoreF)[i]))continue;
                std::vector<typename PolygonType::CoordType> TemplatePos;
                GetRotatedTemplatePos(poly_m.face[i],TemplatePos);
                if ((FixS)||(fixB))
                    DisplaceBySelected(poly_m.face[i],TemplatePos,FixS,fixB);

                //then cumulate the position per vertex
                ScalarType val=vcg::PolyArea(poly_m.face[i]);
                if (val<(AvgArea*0.00001))
                    val=(AvgArea*0.00001);

                ScalarType W=1.0/val;
                if (WeightByQuality)
                    W=poly_m.face[i].Q()+0.00001;

                for (size_t j=0;j<TemplatePos.size();j++)
                {
                    int IndexV=vcg::tri::Index(poly_m,poly_m.face[i].V(j));
                    CoordType Pos=TemplatePos[j];
                    //sum up contributes
                    avgPos[IndexV]+=Pos*W;
                    weightSum[IndexV]+=W;
                }

            }

            //get the laplacian contribute
            std::vector<CoordType> AvVert;
            LaplacianPos(poly_m,AvVert);
            //then update the position
            for (size_t i=0;i<poly_m.vert.size();i++)
            {
                ScalarType alpha=smoothTerm;//PolyNormDeviation(poly_m.face[i]);
                //               if (alpha<0)alpha=0;
                //               if (alpha>1)alpha=1;
                //               if (isnan(alpha))alpha=1;

                CoordType newP=poly_m.vert[i].P();
                //safety checks
                if (weightSum[i]>0)
                    newP=avgPos[i]/weightSum[i];
                if (isnan(newP.X())||isnan(newP.Y())||isnan(newP.Z()))
                     newP=poly_m.vert[i].P();
                if ((newP-poly_m.vert[i].P()).Norm()>poly_m.bbox.Diag())
                    newP=poly_m.vert[i].P();
                //std::cout<<"W "<<weightSum[i]<<std::endl;
                newP=newP*(1-alpha)+AvVert[i]*alpha;
                //newP=AvVert[i];
                if ((fixB)&&(poly_m.vert[i].IsB()))continue;
                if ((FixS)&&(poly_m.vert[i].IsS()))continue;
                poly_m.vert[i].P()=poly_m.vert[i].P()*Damp+
                        newP*(1-Damp);
            }
        }
    }

    template <class TriMeshType>
    static void ReprojectBorder(PolyMeshType &poly_m,
                                TriMeshType &tri_mesh,
                                bool FixS=true)
    {
        //then reproject on border
        for (size_t i=0;i<poly_m.vert.size();i++)
        {
            if (!poly_m.vert[i].IsB())continue;
            if (FixS && poly_m.vert[i].IsS())continue;

            CoordType testPos=poly_m.vert[i].P();
            ScalarType minD=std::numeric_limits<ScalarType>::max();
            CoordType closPos;
            for (size_t j=0;j<tri_mesh.face.size();j++)
                for (size_t k=0;k<3;k++)
                {
                    //check if border edge
                    if (tri_mesh.face[j].FFp(k)!=(&tri_mesh.face[j]))continue;

                    CoordType P0,P1;
                    P0.Import(tri_mesh.face[j].cP0(k));
                    P1.Import(tri_mesh.face[j].cP1(k));
                    vcg::Segment3<ScalarType> Seg(P0,P1);
                    ScalarType testD;
                    CoordType closTest;
                    vcg::SegmentPointDistance(Seg,testPos,closTest,testD);
                    if (testD>minD)continue;
                    minD=testD;
                    closPos=closTest;
                }
            poly_m.vert[i].P()=closPos;
        }
    }

    /*! \brief This function smooth the borders of the polygonal mesh and reproject back to the triangolar one
     * except the vertices that are considered as corner wrt the angleDeg threshold
    */
    template <class TriMeshType>
    static void LaplacianReprojectBorder(PolyMeshType &poly_m,
                                         TriMeshType &tri_mesh,
                                         int nstep=100,
                                         ScalarType Damp=0.5,
                                         ScalarType angleDeg=100)
    {
        //first select corners
        vcg::tri::UpdateFlags<PolyMeshType>::VertexClearS(poly_m);

        //update topology
        vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(poly_m);

        //update border vertices
        vcg::tri::UpdateFlags<PolyMeshType>::VertexBorderFromFaceAdj(poly_m);

        //select corner vertices on the border
        ScalarType angleRad=angleDeg * M_PI / 180;
        vcg::tri::UpdateSelection<PolyMeshType>::VertexCornerBorder(poly_m,angleRad);

        for (int s=0;s<nstep;s++)
        {
            std::vector<CoordType> AvVert;
            LaplacianPos(poly_m,AvVert);

            for (size_t i=0;i<poly_m.vert.size();i++)
            {
                if (!poly_m.vert[i].IsB())continue;
                if (poly_m.vert[i].IsS())continue;
                poly_m.vert[i].P()=poly_m.vert[i].P()*Damp+
                        AvVert[i]*(1-Damp);
            }

            //            //then reproject on border
            //            for (size_t i=0;i<poly_m.vert.size();i++)
            //            {
            //                if (!poly_m.vert[i].IsB())continue;
            //                if (poly_m.vert[i].IsS())continue;

            //                CoordType testPos=poly_m.vert[i].P();
            //                ScalarType minD=std::numeric_limits<ScalarType>::max();
            //                CoordType closPos;
            //                for (size_t j=0;j<tri_mesh.face.size();j++)
            //                    for (size_t k=0;k<3;k++)
            //                    {
            //                        if (tri_mesh.face[j].FFp(k)!=(&tri_mesh.face[j]))continue;

            //                        CoordType P0,P1;
            //                        P0.Import(tri_mesh.face[j].cP0(k));
            //                        P1.Import(tri_mesh.face[j].cP1(k));
            //                        vcg::Segment3<ScalarType> Seg(P0,P1);
            //                        ScalarType testD;
            //                        CoordType closTest;
            //                        vcg::SegmentPointDistance(Seg,testPos,closTest,testD);
            //                        if (testD>minD)continue;
            //                        minD=testD;
            //                        closPos=closTest;
            //                    }
            //                poly_m.vert[i].P()=closPos;
            //            }
            ReprojectBorder(poly_m,tri_mesh);
        }
    }

    /*! \brief This function smooth the borders of the polygonal mesh and reproject back to its border
    */
    static void LaplacianReprojectBorder(PolyMeshType &poly_m,
                                         int nstep=100,
                                         ScalarType Damp=0.5,
                                         ScalarType Angle=100)
    {
        //transform into triangular
        TempMesh GuideSurf;
        vcg::tri::PolygonSupport<TempMesh,PolyMeshType>::ImportFromPolyMesh(GuideSurf,poly_m);
        vcg::tri::UpdateBounding<TempMesh>::Box(GuideSurf);
        vcg::tri::UpdateNormal<TempMesh>::PerVertexNormalizedPerFace(GuideSurf);
        vcg::tri::UpdateTopology<TempMesh>::FaceFace(GuideSurf);
        vcg::tri::UpdateFlags<TempMesh>::FaceBorderFromFF(GuideSurf);

        LaplacianReprojectBorder<TempMesh>(poly_m,GuideSurf,nstep,Damp,Angle);
    }

    /*! \brief This function performs the reprojection of the polygonal mesh onto a triangular one passed as input parameter
    */
    template <class TriMeshType>
    static void LaplacianReproject(PolyMeshType &poly_m,
                                   TriMeshType &tri_mesh,
                                   int nstep=100,
                                   ScalarType DampS=0.5,
                                   ScalarType DampR=0.5,
                                   bool OnlyOnSelected=false)
    {
        typedef typename TriMeshType::FaceType TriFaceType;
        typedef typename TriMeshType::ScalarType TriScalarType;
        typedef typename TriMeshType::CoordType TriCoordType;
        typedef vcg::GridStaticPtr<TriFaceType, TriScalarType> TriMeshGrid;
        TriMeshGrid grid;
        tri::MeshAssert<TriMeshType>::VertexNormalNormalized(tri_mesh);
        //initialize the grid
        grid.Set(tri_mesh.face.begin(),tri_mesh.face.end());

        TriScalarType MaxD=tri_mesh.bbox.Diag();

        for (int s=0;s<nstep;s++)
        {
            std::vector<CoordType> AvVert;
            LaplacianPos(poly_m,AvVert);

            for (size_t i=0;i<poly_m.vert.size();i++)
            {
                if (poly_m.vert[i].IsB()) continue;
                if(poly_m.vert[i].IsD() || (OnlyOnSelected && !poly_m.vert[i].IsS())) continue;
                poly_m.vert[i].P()=poly_m.vert[i].P()*DampS+
                        AvVert[i]*(1-DampS);
            }


            for (size_t i=0;i<poly_m.vert.size();i++)
            {
                if(poly_m.vert[i].IsD() || (OnlyOnSelected && !poly_m.vert[i].IsS())) continue;
                TriCoordType testPos;
                testPos.Import(poly_m.vert[i].P());
                TriCoordType closestPt;
                TriScalarType minDist;
                TriFaceType *f=NULL;
                TriCoordType norm,ip;
                f=vcg::tri::GetClosestFaceBase(tri_mesh,grid,testPos,MaxD,minDist,closestPt,norm,ip);
                CoordType closestImp;
                closestImp.Import(closestPt);
                poly_m.vert[i].P()=poly_m.vert[i].P()*DampR+
                        closestImp*(1-DampR);
                CoordType normalImp;
                normalImp.Import(norm);
                poly_m.vert[i].N()=normalImp;
            }
        }

    }

    static void LaplacianReproject(PolyMeshType &poly_m,
                                   int nstep=100,
                                   ScalarType Damp=0.5,
                                   bool OnlyOnSelected=false)
    {
        //transform into triangular
        TempMesh GuideSurf;
        //vcg::tri::PolygonSupport<TempMesh,PolyMeshType>:(GuideSurf,poly_m);
        TriangulateToTriMesh<TempMesh>(poly_m,GuideSurf);
        vcg::tri::UpdateBounding<TempMesh>::Box(GuideSurf);
        vcg::tri::UpdateNormal<TempMesh>::PerVertexNormalizedPerFace(GuideSurf);
        vcg::tri::UpdateTopology<TempMesh>::FaceFace(GuideSurf);
        vcg::tri::UpdateFlags<TempMesh>::FaceBorderFromFF(GuideSurf);
        LaplacianReproject<TempMesh>(poly_m,GuideSurf,nstep,Damp,0.5,OnlyOnSelected);
    }

    static void Laplacian(PolyMeshType &poly_m,
                          bool FixS=false,
                          int nstep=10,
                          ScalarType Damp=0.5)
    {
        for (int s=0;s<nstep;s++)
        {
            std::vector<CoordType> AvVert;
            LaplacianPos(poly_m,AvVert);

            for (size_t i=0;i<poly_m.vert.size();i++)
            {
                if ((FixS) && (poly_m.vert[i].IsS()))continue;
                poly_m.vert[i].P()=poly_m.vert[i].P()*Damp+
                        AvVert[i]*(1-Damp);
            }
        }

    }

    /*! \brief This function performs the polygon regularization as in "Statics Aware Grid Shells"
     * followed by a reprojection step on the triangle mesh passed as parameter
    */
    template <class TriMeshType>
    static void SmoothReprojectPCA(PolyMeshType &poly_m,
                                   TriMeshType &tri_mesh,
                                   int relaxStep=100,
                                   bool fixS=false,
                                   ScalarType Damp=0.5,
                                   ScalarType SharpDeg=0,
                                   bool WeightByQuality=false,
                                   bool FixB=true)
    {
        //vcg::tri::UpdateFlags<PolyMeshType>::VertexClearS(poly_m);

        vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(poly_m);

        //UpdateBorderVertexFromPFFAdj(poly_m);
        vcg::tri::UpdateFlags<PolyMeshType>::VertexBorderFromFaceAdj(poly_m);

        std::vector<std::vector<vcg::Line3<ScalarType> > > SharpEdge(poly_m.vert.size());
        //first select sharp features
        if (SharpDeg>0)
        {
            for (int i=0;i<(int)poly_m.face.size();i++)
                for (int j=0;j<(int)poly_m.face[i].VN();j++)
                {
                    //check only one side
                    if ((&poly_m.face[i])>=poly_m.face[i].FFp(j))continue;

                    CoordType N0=poly_m.face[i].N();
                    CoordType N1=poly_m.face[i].FFp(j)->N();

                    ScalarType Angle=vcg::Angle(N0,N1);
                    if (fabs(Angle)>(SharpDeg* (M_PI / 180.0)))
                    {
                        CoordType Pos0=poly_m.face[i].V0(j)->P();
                        CoordType Pos1=poly_m.face[i].V1(j)->P();
                        CoordType Ori=Pos0;
                        CoordType Dir=Pos1-Pos0;
                        Dir.Normalize();
                        vcg::Line3<ScalarType> L(Ori,Dir);
                        int Index0=vcg::tri::Index(poly_m,poly_m.face[i].V0(j));
                        int Index1=vcg::tri::Index(poly_m,poly_m.face[i].V1(j));
                        SharpEdge[Index0].push_back(L);
                        SharpEdge[Index1].push_back(L);
                    }
                }
            for (size_t i=0;i<poly_m.vert.size();i++)
            {
                if (SharpEdge[i].size()==0)continue;
                if (SharpEdge[i].size()>2)poly_m.vert[i].SetS();
            }
        }
        //        if (fixIrr)
        //        {
        //            vcg::tri::UpdateQuality<PolyMeshType>::VertexValence(poly_m);
        //            for (size_t i=0;i<poly_m.vert.size();i++)
        //            {
        //                if (poly_m.vert[i].IsB())continue;
        //                if (poly_m.vert[i].Q()==4)continue;
        //                poly_m.vert[i].SetS();
        //            }
        //        }


        typedef typename TriMeshType::FaceType FaceType;
        typedef vcg::GridStaticPtr<FaceType, typename TriMeshType::ScalarType> TriMeshGrid;
        TriMeshGrid grid;

        //initialize the grid
        grid.Set(tri_mesh.face.begin(),tri_mesh.face.end());

        ScalarType MaxD=tri_mesh.bbox.Diag();

        //        //update quality as area
        //        for (size_t i=0;i<poly_m.face.size();i++)
        //          poly_m.face[i].Q()=vcg::PolyArea(poly_m.face[i]);

        //        for (size_t i=0;i<poly_m.vert.size();i++)
        //        {
        //            typename TriMeshType::CoordType testPos;
        //            testPos.Import(poly_m.vert[i].P());
        //            typename TriMeshType::CoordType closestPt;
        //            typename TriMeshType::ScalarType minDist;
        //            typename TriMeshType::FaceType *f=NULL;
        //            typename TriMeshType::CoordType norm,ip;
        //            f=vcg::tri::GetClosestFaceBase(tri_mesh,grid,testPos,MaxD,minDist,closestPt,norm,ip);
        //            //poly_m.vert[i].N().Import(norm);
        //        }

        for(int k=0;k<relaxStep;k++)
        {
            //smooth PCA step
            SmoothPCA(poly_m,1,Damp,fixS,true,0.1,FixB,WeightByQuality);
            //reprojection step
            //laplacian smooth step
            //Laplacian(poly_m,Damp,1);

            for (size_t i=0;i<poly_m.vert.size();i++)
            {
                typename TriMeshType::CoordType testPos;
                testPos.Import(poly_m.vert[i].P());
                typename TriMeshType::CoordType closestPt;
                typename TriMeshType::ScalarType minDist;
                if ((FixB)&&(poly_m.vert[i].IsB()))
                {continue;}
                else
                    if (SharpEdge[i].size()==0)//reproject onto original mesh
                    {
                        FaceType *f=NULL;
                        typename TriMeshType::CoordType norm,ip;
                        f=vcg::tri::GetClosestFaceBase(tri_mesh,grid,testPos,MaxD,minDist,closestPt,norm,ip);
                        poly_m.vert[i].P().Import(testPos*Damp+closestPt*(1-Damp));
                        //poly_m.vert[i].N().Import(norm);
                    }
                    else //reproject onto segments
                    {
                        CoordType av_closest(0,0,0);
                        size_t sum=0;
                        for (size_t j=0;j<SharpEdge[i].size();j++)
                        {
                            CoordType currPos;
                            currPos.Import(testPos);
                            CoordType closest;
                            ScalarType dist;
                            vcg::LinePointDistance(SharpEdge[i][j],currPos,closest,dist);
                            av_closest+=closest;
                            sum++;
                        }
                        assert(sum>0);
                        poly_m.vert[i].P()=av_closest/sum;
                    }
            }
            if (!FixB)
                ReprojectBorder(poly_m,tri_mesh,true);
            UpdateFaceNormals(poly_m);
            vcg::tri::UpdateNormal<PolyMeshType>::PerVertexFromCurrentFaceNormal(poly_m);
        }

    }


    template <class TriMeshType>
    static void TriangulateToTriMesh(PolyMeshType &poly_m,TriMeshType &triangle_mesh, bool alsoTriangles = true)
    {
        triangle_mesh.Clear();

        PolyMeshType PolySwap;
        vcg::tri::Append<PolyMeshType,PolyMeshType>::Mesh(PolySwap,poly_m);
        Triangulate(PolySwap, alsoTriangles);

        //then copy onto the triangle mesh
        vcg::tri::Append<TriMeshType,PolyMeshType>::Mesh(triangle_mesh,PolySwap);
    }

    /*! \brief This function performs the polygon regularization as in "Statics Aware Grid Shells"
     * followed by a reprojection step on the original mesh
    */
    static void SmoothReprojectPCA(PolyMeshType &poly_m,
                                   int relaxStep=100,
                                   bool fixS=false,
                                   ScalarType Damp=0.5,
                                   ScalarType SharpDeg=0,
                                   bool WeightByQuality=false,
                                   bool FixB=true)
    {
        //transform into triangular
        TempMesh GuideSurf;

        //vcg::tri::PolygonSupport<TempMesh,PolyMeshType>:(GuideSurf,poly_m);
        TriangulateToTriMesh<TempMesh>(poly_m,GuideSurf);
        vcg::tri::UpdateBounding<TempMesh>::Box(GuideSurf);
        vcg::tri::UpdateNormal<TempMesh>::PerVertexNormalizedPerFace(GuideSurf);
        vcg::tri::UpdateTopology<TempMesh>::FaceFace(GuideSurf);
        vcg::tri::UpdateFlags<TempMesh>::FaceBorderFromFF(GuideSurf);

        //optimize it
        vcg::PolygonalAlgorithm<PolyMeshType>::SmoothReprojectPCA<TempMesh>(poly_m,GuideSurf,relaxStep,fixS,Damp,SharpDeg,WeightByQuality,FixB);
    }

    static void Reproject(PolyMeshType &poly_m,
                          PolyMeshType &target)
    {
        vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(poly_m);
        vcg::tri::UpdateFlags<PolyMeshType>::VertexBorderFromFaceAdj(poly_m);

        //transform into triangular
        TempMesh GuideSurf;
        //vcg::tri::PolygonSupport<TempMesh,PolyMeshType>:(GuideSurf,poly_m);
        TriangulateToTriMesh<TempMesh>(target,GuideSurf);
        vcg::tri::UpdateBounding<TempMesh>::Box(GuideSurf);
        vcg::tri::UpdateNormal<TempMesh>::PerVertexNormalizedPerFace(GuideSurf);
        vcg::tri::UpdateTopology<TempMesh>::FaceFace(GuideSurf);
        vcg::tri::UpdateFlags<TempMesh>::FaceBorderFromFF(GuideSurf);

        //initialize the grid
        typedef typename TempMesh::FaceType FaceType;
        typedef vcg::GridStaticPtr<FaceType, typename TempMesh::ScalarType> TriMeshGrid;
        TriMeshGrid grid;
        grid.Set(GuideSurf.face.begin(),GuideSurf.face.end());

        ScalarType MaxD=GuideSurf.bbox.Diag();

        for (size_t i=0;i<poly_m.vert.size();i++)
        {
            //reproject on border later
            if (poly_m.vert[i].IsB())continue;
            typename TempMesh::CoordType testPos;
            testPos.Import(poly_m.vert[i].P());
            typename TempMesh::CoordType closestPt;
            typename TempMesh::ScalarType minDist;
            typename TempMesh::FaceType *f=NULL;
            typename TempMesh::CoordType norm,ip;
            f=vcg::tri::GetClosestFaceBase(GuideSurf,grid,testPos,MaxD,minDist,closestPt,norm,ip);
            poly_m.vert[i].P()=closestPt;
        }
        //then reprojec the border
        ReprojectBorder(poly_m,GuideSurf);
    }

    template <class TriMesh>
    static void ReprojectonTriMesh(PolyMeshType &poly_m,
                                   TriMesh &target)
    {
        vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(poly_m);
        vcg::tri::UpdateFlags<PolyMeshType>::VertexBorderFromFaceAdj(poly_m);


        //initialize the grid
        typedef typename TriMesh::FaceType FaceType;
        typedef vcg::GridStaticPtr<FaceType, typename TriMesh::ScalarType> TriMeshGrid;
        TriMeshGrid grid;
        grid.Set(target.face.begin(),target.face.end());

        ScalarType MaxD=target.bbox.Diag();

        for (size_t i=0;i<poly_m.vert.size();i++)
        {
            //reproject on border later
            if (poly_m.vert[i].IsB())continue;
            typename TriMesh::CoordType testPos;
            testPos.Import(poly_m.vert[i].P());
            typename TriMesh::CoordType closestPt;
            typename TriMesh::ScalarType minDist;
            typename TriMesh::FaceType *f=NULL;
            typename TriMesh::CoordType norm,ip;
            f=vcg::tri::GetClosestFaceBase(target,grid,testPos,MaxD,minDist,closestPt,norm,ip);
            poly_m.vert[i].P()=closestPt;
        }
        //then reprojec the border
        ReprojectBorder(poly_m,target);
    }

    /*! \brief This function return average edge size
    */
    static ScalarType AverageEdge(const PolyMeshType &poly_m)
    {
        ScalarType AvL=0;
        size_t numE=0;
        for (size_t i=0;i<poly_m.face.size();i++)
        {
            int NumV=poly_m.face[i].VN();
            for (int j=0;j<NumV;j++)
            {
                CoordType pos0=poly_m.face[i].cV(j)->P();
                CoordType pos1=poly_m.face[i].cV((j+1)%NumV)->P();
                AvL+=(pos0-pos1).Norm();
                numE++;
            }
        }
        AvL/=numE;
        return AvL;
    }
    

    /*! \brief This function remove valence 2 faces from the mesh
    */
    static void RemoveValence2Faces(PolyMeshType &poly_m)
    {
        for (size_t i=0;i<poly_m.face.size();i++)
        {
            if (poly_m.face[i].VN()>=3)continue;
            vcg::tri::Allocator<PolyMeshType>::DeleteFace(poly_m,poly_m.face[i]);
        }

        //then remove unreferenced vertices
        vcg::tri::Clean<PolyMeshType>::RemoveUnreferencedVertex(poly_m);
        vcg::tri::Allocator<PolyMeshType>::CompactEveryVector(poly_m);

    }
    
    /*! \brief This function remove valence 2 vertices on the border by considering the degree threshold
     * bacause there could be eventually some corner that should be preserved
    */
    static void RemoveValence2Vertices(PolyMeshType &poly_m,
                                       ScalarType corner_degree=25)
    {
        //update topology
        vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(poly_m);

        //update border vertices
        //UpdateBorderVertexFromPFFAdj(poly_m);
        vcg::tri::UpdateFlags<PolyMeshType>::VertexBorderFromFaceAdj(poly_m);

        vcg::tri::UpdateFlags<PolyMeshType>::VertexClearS(poly_m);

        //select corners
        for (size_t i=0;i<poly_m.face.size();i++)
        {
            if (poly_m.face[i].IsD())continue;

            //get vertices of the face
            int NumV=poly_m.face[i].VN();

            for (int j=0;j<NumV;j++)
            {
                VertexType *v0=poly_m.face[i].V((j+NumV-1)%NumV);
                VertexType *v1=poly_m.face[i].V(j);
                VertexType *v2=poly_m.face[i].V((j+1)%NumV);
                //must be 3 borders
                bool IsB=((v0->IsB())&&(v1->IsB())&&(v2->IsB()));
                CoordType dir0=(v0->P()-v1->P());
                CoordType dir1=(v2->P()-v1->P());
                dir0.Normalize();
                dir1.Normalize();
                ScalarType testDot=(dir0*dir1);
                if ((IsB)&&(testDot>(-cos(corner_degree* (M_PI / 180.0)))))
                    v1->SetS();
            }
        }

        typename PolyMeshType::template PerVertexAttributeHandle<size_t> valenceVertH =
                vcg::tri::Allocator<PolyMeshType>:: template GetPerVertexAttribute<size_t> (poly_m);

        //initialize to zero
        for (size_t i=0;i<poly_m.vert.size();i++)
            valenceVertH[i]=0;

        //then sum up the valence
        for (size_t i=0;i<poly_m.face.size();i++)
            for (int j=0;j<poly_m.face[i].VN();j++)
                valenceVertH[poly_m.face[i].V(j)]++;

        //cannot collapse triangular vertices otherwise will collapse to a segment
        for (size_t i=0;i<poly_m.face.size();i++)
        {
            if (poly_m.face[i].VN()>3)continue;
            for (int j=0;j<poly_m.face[i].VN();j++)
                valenceVertH[poly_m.face[i].V(j)]=3;
        }

        //then re-elaborate the faces
        for (size_t i=0;i<poly_m.face.size();i++)
        {
            if (poly_m.face[i].IsD())continue;

            //get vertices of the face
            int NumV=poly_m.face[i].VN();

            std::vector<VertexType*> FaceV;
            for (int j=0;j<NumV;j++)
            {
                VertexType *v=poly_m.face[i].V(j);
                assert(!v->IsD());
                //if ((!v->IsS()) && (v->IsB()) && (valenceVertH[v]==1)) continue;
                if ((!v->IsS()) && (v->IsB()) && (valenceVertH[v]==1)) continue;
                if ((!v->IsB()) && (valenceVertH[v]<3)) continue;
                //if (!v->IsS()) continue;
                FaceV.push_back(v);
            }

            //then deallocate face
            if ((int)FaceV.size()==NumV)continue;

            //otherwise deallocate and set new vertices
            poly_m.face[i].Dealloc();
            poly_m.face[i].Alloc(FaceV.size());
            for (size_t j=0;j<FaceV.size();j++)
                poly_m.face[i].V(j)=FaceV[j];
        }

        //then remove unreferenced vertices
        vcg::tri::Clean<PolyMeshType>::RemoveUnreferencedVertex(poly_m);
        vcg::tri::Allocator<PolyMeshType>::CompactEveryVector(poly_m);

        vcg::tri::Allocator<PolyMeshType>::DeletePerVertexAttribute(poly_m,valenceVertH);
    }

    /*! \brief This function collapse small edges which are on the boundary of the mesh
     * this is sometimes useful to remove small edges coming out from a quadrangulation which is not
     * aligned to boundaries
    */
    static bool CollapseBorderSmallEdges(PolyMeshType &poly_m,
                                         const ScalarType perc_average=0.3)
    {
        //compute the average edge
        ScalarType AvEdge=AverageEdge(poly_m);
        ScalarType minLimit=AvEdge*perc_average;
        bool collapsed=false;
        while(CollapseBorderSmallEdgesStep(poly_m,minLimit)){collapsed=true;};

        RemoveValence2Faces(poly_m);

        //RemoveValence2BorderVertices(poly_m);
        RemoveValence2Vertices(poly_m);
        return collapsed;
    }

    /*! \brief This function use a local global approach to flatten polygonal faces
     * the approach is similar to "Shape-Up: Shaping Discrete Geometry with Projections"
    */
    static ScalarType FlattenFaces(PolyMeshType &poly_m, size_t steps=100,bool OnlySFaces=false)
    {
        ScalarType MaxDispl=0;
        for (size_t s=0;s<steps;s++)
        {
            std::vector<std::vector<CoordType> > VertPos(poly_m.vert.size());

            for (size_t i=0;i<poly_m.face.size();i++)
            {
                if (poly_m.face[i].IsD())continue;

                if (OnlySFaces && (!poly_m.face[i].IsS()))continue;
                //get vertices of the face
                int NumV=poly_m.face[i].VN();
                if (NumV<=3)continue;

                //save vertice's positions
                std::vector<CoordType> FacePos;
                for (int j=0;j<NumV;j++)
                {
                    VertexType *v=poly_m.face[i].V(j);
                    assert(!v->IsD());
                    FacePos.push_back(v->P());
                }

                //then fit the plane
                vcg::Plane3<ScalarType> FitPl;
                vcg::FitPlaneToPointSet(FacePos,FitPl);

                //project each point onto fitting plane
                for (int j=0;j<NumV;j++)
                {
                    VertexType *v=poly_m.face[i].V(j);
                    int IndexV=vcg::tri::Index(poly_m,v);
                    CoordType ProjP=FitPl.Projection(v->P());
                    VertPos[IndexV].push_back(ProjP);
                }
            }

            for (size_t i=0;i<poly_m.vert.size();i++)
            {
                CoordType AvgPos(0,0,0);

                for (size_t j=0;j<VertPos[i].size();j++)
                    AvgPos+=VertPos[i][j];

                if (VertPos[i].size()==0)continue;

                AvgPos/=(ScalarType)VertPos[i].size();

                MaxDispl=std::max(MaxDispl,(poly_m.vert[i].P()-AvgPos).Norm());
                poly_m.vert[i].P()=AvgPos;
            }
        }
        return MaxDispl;
    }

    static ScalarType Area(PolyMeshType &poly_m)
    {
        ScalarType MeshArea=0;
        for (size_t i=0;i<poly_m.face.size();i++)
            MeshArea+=vcg::PolyArea(poly_m.face[i]);
        return MeshArea;
    }

    static void InitQualityVertVoronoiArea(PolyMeshType &poly_m)
    {
        for (size_t i=0;i<poly_m.vert.size();i++)
            poly_m.vert[i].Q()=0;

        for (size_t i=0;i<poly_m.face.size();i++)
        {
            //            ScalarType AreaF=vcg::PolyArea(poly_m.face[i]);
            size_t sizeV=poly_m.face[i].VN()-1;
            CoordType baryF=vcg::PolyBarycenter(poly_m.face[i]);
            for (int j=0;j<poly_m.face[i].VN();j++)
            {
                CoordType P0=poly_m.face[i].P((j+sizeV-1)%sizeV);
                CoordType P1=poly_m.face[i].P(j);
                CoordType P2=poly_m.face[i].P1(j);
                vcg::Triangle3<ScalarType> T0(P1,(P0+P1)/2,baryF);
                vcg::Triangle3<ScalarType> T1(P1,(P1+P2)/2,baryF);

                poly_m.face[i].V(j)->Q()+=vcg::DoubleArea(T0)/2;
                poly_m.face[i].V(j)->Q()+=vcg::DoubleArea(T1)/2;
            }
        }
    }

    static ScalarType InitQualityFaceTorsion(PolyMeshType &poly_m)
    {
        UpdateFaceNormalByFitting(poly_m);
        vcg::tri::UpdateNormal<PolyMeshType>::PerVertexFromCurrentFaceNormal(poly_m);
        ScalarType MaxA=0;
        for (size_t i=0;i<poly_m.face.size();i++)
        {
            poly_m.face[i].Q()=PolygonTorsion(poly_m.face[i]);
            MaxA=std::max(MaxA,poly_m.face[i].Q());
        }
        return MaxA;
    }

    static ScalarType InitQualityFaceBending(PolyMeshType &poly_m)
    {
        UpdateFaceNormalByFitting(poly_m);
        vcg::tri::UpdateNormal<PolyMeshType>::PerVertexFromCurrentFaceNormal(poly_m);
        ScalarType MaxA=0;
        for (size_t i=0;i<poly_m.face.size();i++)
        {
            poly_m.face[i].Q()=PolygonBending(poly_m.face[i]);
            MaxA=std::max(MaxA,poly_m.face[i].Q());
        }
        return MaxA;
    }


    static void InitQualityVertEdgeLenght(PolyMeshType &poly_m)
    {
        for (size_t i=0;i<poly_m.vert.size();i++)
            poly_m.vert[i].Q()=0;

        for (size_t i=0;i<poly_m.face.size();i++)
        {
            for (int j=0;j<poly_m.face[i].VN();j++)
            {
                FaceType *f=&poly_m.face[i];
                FaceType *f1=f->FFp(j);
                if (f>f1)continue;
                ScalarType L=(poly_m.face[i].P0(j)-poly_m.face[i].P1(j)).Norm();
                poly_m.face[i].V0(j)->Q()+=L;
                poly_m.face[i].V1(j)->Q()+=L;
            }
        }
    }

    static void InterpolateQualityVertFormFaces(PolyMeshType &poly_m)
    {
        std::vector<ScalarType> SumW(poly_m.vert.size(),0);

        for (size_t i=0;i<poly_m.vert.size();i++)
            poly_m.vert[i].Q()=0;

        for (size_t i=0;i<poly_m.face.size();i++)
        {
            ScalarType AreaF=vcg::PolyArea(poly_m.face[i]);
            for (size_t j=0;j<poly_m.face[i].VN();j++)
            {
                poly_m.face[i].V(j)->Q()+=AreaF*(ScalarType)poly_m.face[i].Q();
                size_t IndexV=vcg::tri::Index(poly_m,poly_m.face[i].V(j));
                SumW[IndexV]+=AreaF;
            }
        }
        for (size_t i=0;i<poly_m.vert.size();i++)
        {
            if (SumW[i]>0)
                poly_m.vert[i].Q()/=SumW[i];
            else
                poly_m.vert[i].Q()=0;
        }
    }


    static void ClosestPoint(const PolyMeshType &poly_m,const CoordType &pos,
                             int &CloseF,CoordType &ClosePos)
    {
        ScalarType minD=std::numeric_limits<ScalarType>::max();
        CloseF=-1;
        for (size_t i=0;i<poly_m.face.size();i++)
        {
            CoordType closeTest;
            ScalarType currD=vcg::PolygonPointDistance(poly_m.face[i],pos,closeTest);
            if (currD>minD)continue;
            minD=currD;
            CloseF=i;
            ClosePos=closeTest;
        }
    }

    /*! \brief Triangulate a polygonal face with a triangle fan.
     * \returns pointer to the newly added vertex.
     */
    static VertexPointer Triangulate(PolyMeshType & poly_m, size_t IndexF)
    {

        const CoordType bary = vcg::PolyBarycenter(poly_m.face[IndexF]);
        size_t sizeV = poly_m.face[IndexF].VN();

        //add the new vertex
        VertexPointer newV = &(*vcg::tri::Allocator<PolyMeshType>::AddVertex(poly_m,bary));

        //then reupdate the faces
        for (size_t j=0;j<(sizeV-1);j++)
        {
            VertexType * v0=poly_m.face[IndexF].V0(j);
            VertexType * v1=poly_m.face[IndexF].V1(j);
            VertexType * v2=newV;

            vcg::tri::Allocator<PolyMeshType>::AddFaces(poly_m,1);

            poly_m.face.back().Alloc(3);
            poly_m.face.back().V(0)=v0;
            poly_m.face.back().V(1)=v1;
            poly_m.face.back().V(2)=v2;
            poly_m.face.back().Q()=poly_m.face[IndexF].Q();
        }

        VertexType * v0=poly_m.face[IndexF].V0((sizeV-1));
        VertexType * v1=poly_m.face[IndexF].V1((sizeV-1));
        poly_m.face[IndexF].Dealloc();
        poly_m.face[IndexF].Alloc(3);
        poly_m.face[IndexF].V(0)=v0;
        poly_m.face[IndexF].V(1)=v1;
        poly_m.face[IndexF].V(2)=newV;
        return newV;
    }

    static void ReorderFaceVert(FaceType &f,const size_t &StartI)
    {
        if (StartI==0)return;
        size_t sizeN=f.VN();
        assert(StartI>=0);
        assert(StartI<sizeN);
        std::vector<VertexType*> NewV;
        for (size_t i=0;i<sizeN;i++)
        {
            int IndexV=(i+StartI)%sizeN;
            NewV.push_back(f.V(IndexV));
        }
        //then reset all vertices
        for (size_t i=0;i<sizeN;i++)
            f.V(i)=NewV[i];
    }

    static void MergeAlongEdge(PolyMeshType &poly_m,
                               FaceType &f,
                               const size_t &EdgeI)
    {
        //cannot be a border
        assert(f.FFp(EdgeI)!=&f);
        FaceType *f1=f.FFp(EdgeI);
        int EdgeI1=f.FFi(EdgeI);

        //sort first face
        int FirstV0=(EdgeI+1) % f.VN();
        ReorderFaceVert(f,FirstV0);
        int FirstV1=(EdgeI1+1)%f1->VN();
        ReorderFaceVert(*f1,FirstV1);

        std::vector<VertexType*> NewV;
        for (size_t i=0;i<(f.VN()-1);i++)
            NewV.push_back(f.V(i));

        for (size_t i=0;i<(f1->VN()-1);i++)
            NewV.push_back(f1->V(i));

        f.Dealloc();
        f.Alloc(NewV.size());
        for (size_t i=0;i<NewV.size();i++)
            f.V(i)=NewV[i];

        vcg::tri::Allocator<PolyMeshType>::DeleteFace(poly_m,*f1);
    }

    static void MergeAlongEdges(PolyMeshType &poly_m,
                                const std::vector<FaceType*> &PolyF,
                                const std::vector<size_t> &EdgeI)
    {
        //create a table with all edges that have to be merged
        std::set<std::pair<CoordType,CoordType> > NeedMerge;
        for (size_t i=0;i<PolyF.size();i++)
        {
            CoordType P0=PolyF[i]->P0(EdgeI[i]);
            CoordType P1=PolyF[i]->P1(EdgeI[i]);
            std::pair<CoordType,CoordType> key(std::min(P0,P1),std::max(P0,P1));
            NeedMerge.insert(key);
        }

        //then cycle and collapse
        do{
            for (size_t i=0;i<poly_m.face.size();i++)
            {
                if (poly_m.face[i].IsD())continue;
                for (size_t j=0;j<poly_m.face[i].VN();j++)
                {
                    CoordType P0=poly_m.face[i].P0(j);
                    CoordType P1=poly_m.face[i].P1(j);
                    std::pair<CoordType,CoordType> key(std::min(P0,P1),std::max(P0,P1));
                    if (NeedMerge.count(key)==0)continue;

                    //do the merge
                    MergeAlongEdge(poly_m,poly_m.face[i],j);
                    //remove it
                    NeedMerge.erase(key);
                    break;
                }
            }
            vcg::tri::Allocator<PolyMeshType>::CompactEveryVector(poly_m);
        }while (!NeedMerge.empty());
    }

    static void Triangulate(PolyMeshType &poly_m,
                            bool alsoTriangles = true,
                            bool OnlyS=false)
    {
        size_t size0 = poly_m.face.size();
        if (alsoTriangles)
        {
            for (size_t i=0; i<size0; i++)
            {
                if ((OnlyS)&&(!poly_m.face[i].IsS()))continue;
                Triangulate(poly_m, i);
            }
        }
        else
        {
            for (size_t i=0; i<size0; i++)
            {
                if ((OnlyS)&&(!poly_m.face[i].IsS()))continue;
                if (poly_m.face[i].VN() > 3)
                {
                    Triangulate(poly_m, i);
                }
            }
        }
    }

};

}//end namespace vcg

#endif
