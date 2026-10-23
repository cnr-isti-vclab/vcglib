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

#ifndef VCG_COMPLEX_ALGORITHMS_PARAMETRIZATION_UV_UTILS_H
#define VCG_COMPLEX_ALGORITHMS_PARAMETRIZATION_UV_UTILS_H

#include <vcg/space/point2.h>
#include <vcg/space/box2.h>
#include <vcg/space/triangle2.h>
#include <vcg/space/triangle3.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/update/topology.h>

#include <utility>
#include <vector>

namespace vcg {
namespace tri {
template <class MeshType>
class UV_Utils
{
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::FacePointer FacePointer;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::VertexIterator VertexIterator;
    typedef typename MeshType::FaceIterator FaceIterator;
    typedef typename vcg::Point2<ScalarType> UVCoordType;

public:

    /// calculate the area in UV space
    static ScalarType PerVertUVArea(MeshType &m)
    {
        FaceIterator fi;
        ScalarType Area=0;
        for (fi=m.face.begin();fi!=m.face.end();fi++)
        {
            if ((*fi).IsD()) continue;
            UVCoordType E0= (*fi).V(1)->T().P()-(*fi).V(0)->T().P();
            UVCoordType E1= (*fi).V(2)->T().P()-(*fi).V(0)->T().P();
            ScalarType doubleA=fabs(E0^E1);
            Area+=doubleA/2;
        }
        return Area;
    }

	// uniform scale UV coordinates so that area in UV space is the same of the area of surface of the mesh
	static void PerWedgeRegularizeTexArea(MeshType &m)
	{
		float areaTex=0;
		float areaGeo=0;

		vcg::Box2f UVBox = tri::UV_Utils<MeshType>::PerWedgeUVBox(m);
		for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
		{
			areaTex+= fabs((fi->WT(1).P() - fi->WT(0).P()) ^ (fi->WT(2).P() - fi->WT(0).P())) ;
			areaGeo+= DoubleArea(*fi);
		}

		float ratio = sqrt(areaGeo/areaTex);

		for(FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
		{
			for(int j=0;j<3;++j)
				fi->WT(j).P() = (fi->WT(j).P()-UVBox.min) *ratio;
		}
	}

	static void PerWedgeScaleToUnitSpace(MeshType &m)
	{
		vcg::Box2f UVBox = tri::UV_Utils<MeshType>::PerWedgeUVBox(m);

		Point2f diff = UVBox.max - UVBox.min;
		for(auto &f : m.face) {
			for(int i = 0; i < 3; i++) {
				f.WT(i).P() -= UVBox.min;
				f.WT(i).P().X() /= diff.X();
				f.WT(i).P().Y() /= diff.Y();
			}
		}
	}

	static void PerVertScaleToUnitSpace(MeshType &m)
	{
		vcg::Box2f UVBox = tri::UV_Utils<MeshType>::PerVertUVBox(m);

		Point2f diff = UVBox.max - UVBox.min;

		for(auto &v : m.vert) {
			v.T().P() -= UVBox.min;
			v.T().P().X() /= diff.X();
			v.T().P().Y() /= diff.Y();
		}
	}


    ///scale vert UV to match 3D area
    static void ScaleVertUVToMatchArea(MeshType &m)
    {
        FaceIterator fi;
        ScalarType Area3D=0;
        for (fi=m.face.begin();fi!=m.face.end();fi++)
        {
            if ((*fi).IsD()) continue;
            Area3D+=vcg::DoubleArea((*fi))/2;
        }
        ScalarType Area2D=PerVertUVArea(m);
        ScalarType ScaleFact=sqrt( Area3D / Area2D );

        VertexIterator vi;
        for (vi=m.vert.begin();vi!=m.vert.end();vi++)
        {
            if ((*vi).IsD()) continue;
            (*vi).T().P()*=ScaleFact;
        }
    }

    ///calculate the BBox in UV space
    static vcg::Box2<ScalarType> PerWedgeUVBox(MeshType &m)
    {
        vcg::Box2<ScalarType> UVBox;
        FaceIterator fi;
        for (fi=m.face.begin();fi!=m.face.end();fi++)
        {
            if ((*fi).IsD()) continue;
            for (int i=0;i<3;i++)
                UVBox.Add((*fi).WT(i).P());
        }
        return UVBox;
    }

    ///calculate the BBox in UV space
    static vcg::Box2<ScalarType> PerVertUVBox(const MeshType &m)
    {
        vcg::Box2<ScalarType> UVBox;
        for (auto vi=m.vert.begin();vi!=m.vert.end();vi++)
        {
            if ((*vi).IsD()) continue;
            UVBox.Add((*vi).T().P());
        }
        return UVBox;
    }

    void PerWedgeMakeUnitaryUV(MeshType &m)
    {
        vcg::Box2<typename MeshType::ScalarType> UVBox = PerWedgeUVBox(m);

        typename MeshType::FaceIterator fi;
        Point2f boxSize(UVBox.max-UVBox.min);
        for (fi=m.face.begin();fi!=m.face.end();fi++)
        {
            if ((*fi).IsD()) continue;
            for (int i=0;i<3;i++)
            {
                (*fi).WT(i).U() = ((*fi).WT(i).U()-UVBox.min[0])/boxSize[0] ;
                (*fi).WT(i).V() = ((*fi).WT(i).V()-UVBox.min[1])/boxSize[1] ;
            }
        }
    }
    ///transform curvature to UV space
    static UVCoordType Coord3DtoUV(FaceType &f,const CoordType &dir)
    {
        ///then transform to UV
        CoordType bary3d=(f.P(0)+f.P(1)+f.P(2))/3.0;
        UVCoordType baryUV=(f.WT(0).P()+f.WT(1).P()+f.WT(2).P())/3.0;
        CoordType dir3d=bary3d+dir;
        CoordType baryCoordsUV;
        vcg::InterpolationParameters<FaceType,ScalarType>(f,dir3d,baryCoordsUV);
        UVCoordType dirUV=baryCoordsUV.X()*f.WT(0).P()+
                baryCoordsUV.Y()*f.WT(1).P()+
                baryCoordsUV.Z()*f.WT(2).P()-baryUV;
        dirUV.Normalize();
        return dirUV;
    }

    static void GloballyMirrorX(MeshType &m)
    {
        vcg::Box2<ScalarType> BBuv=PerVertUVBox(m);
        ScalarType Xmin=BBuv.min.X();
        ScalarType Xmax=BBuv.max.X();
        ScalarType XAv=(Xmax+Xmin)/2;
        VertexIterator vi;
        for (vi=m.vert.begin();vi!=m.vert.end();vi++)
        {
            ScalarType distAV=(*vi).T().P().X()-XAv;
            (*vi).T().P().X()=XAv-distAV;
        }
    }

    static void GloballyRotate(MeshType &m,ScalarType Angle)
    {
        vcg::Box2<ScalarType> BB=PerWedgeUVBox(m);
        UVCoordType Origin=BB.Center();
        typename MeshType::FaceIterator fi;
        for (fi=m.face.begin();fi!=m.face.end();fi++)
        {
            if ((*fi).IsD()) continue;
            for (int i=0;i<3;i++)
            {
                (*fi).WT(i).P()-=Origin;
                ScalarType X1=(*fi).WT(i).P().X()*cos(Angle)-(*fi).WT(i).P().Y()*sin(Angle);
                ScalarType Y1=(*fi).WT(i).P().X()*cos(Angle)+(*fi).WT(i).P().Y()*sin(Angle);
                (*fi).WT(i).P().X()=X1;
                (*fi).WT(i).P().Y()=Y1;
                (*fi).WT(i).P()+=Origin;
            }
        }
    }

    static void GloballyRotatePerVert(MeshType &m,ScalarType Angle)
    {
        vcg::Box2<ScalarType> BB=PerVertUVBox(m);
        UVCoordType Origin=BB.Center();
        typename MeshType::VertexIterator vi;
        for (vi=m.vert.begin();vi!=m.vert.end();vi++)
        {
            if ((*vi).IsD()) continue;
            (*vi).T().P()-=Origin;
            ScalarType X1=(*vi).T().P().X()*cos(Angle)-(*vi).T().P().Y()*sin(Angle);
            ScalarType Y1=(*vi).T().P().X()*cos(Angle)+(*vi).T().P().Y()*sin(Angle);
            (*vi).T().P().X()=X1;
            (*vi).T().P().Y()=Y1;
            (*vi).T().P()+=Origin;
        }
    }

    static void LaplacianUVVert(MeshType &m,bool fix_borders=false,int steps=3)
    {
        FaceIterator fi;
        for (int s=0;s<steps;s++)
        {
            std::vector<int> num(m.vert.size(),0);
            std::vector<UVCoordType> UVpos(m.vert.size(),UVCoordType(0,0));
            for (fi=m.face.begin();fi!=m.face.end();fi++)
            {
                for (int j=0;j<3;j++)
                {
                    VertexType *v0=(*fi).V(0);
                    VertexType *v1=(*fi).V1(0);
                    VertexType *v2=(*fi).V2(0);
                    assert(v0!=v1);
                    assert(v1!=v2);
                    assert(v0!=v2);
                    UVCoordType uv1=v1->T().P();
                    UVCoordType uv2=v2->T().P();
                    int index=v0-&(m.vert[0]);
                    num[index]+=2;
                    UVpos[index]+=uv1;
                    UVpos[index]+=uv2;
                }
            }
            VertexIterator vi;
            for (int i=0;i<m.vert.size();i++)
            {
                if ((fix_borders)&&(m.vert[i].IsB()))continue;
                if (num[i]==0)continue;
                m.vert[i].T().P()=UVpos[i]/(ScalarType)num[i];
            }
        }
    }

    static void CopyVertUVWedge(MeshType &m)
    {
        for (size_t i=0;i<m.face.size();i++)
            for (size_t j=0;j<3;j++)
                m.face[i].WT(j).P()=m.face[i].V(j)->T().P();
    }

    static void CopyWedgeVertUV(MeshType &m,bool onlyS=false)
    {
        for (size_t i=0;i<m.face.size();i++)
        {
            if ((onlyS)&&(!m.face[i].IsS()))continue;
            for (int j=0;j<3;j++)
                m.face[i].V(j)->T().P()=m.face[i].WT(j).P();
        }
    }

    /// Partition the faces into UV charts (also called UV islands): the maximal sets of faces
    /// connected across edges that the texture coordinates do not cut, that is the edges on
    /// whose two sides the faces agree on the UV of both ends. Fills faceChart with the chart
    /// of each face (-1 for a deleted face) and returns the number of charts.
    /// The face-face adjacency of m is rebuilt from the per-wedge texture coordinates
    /// (UpdateTopology::FaceFaceFromTexCoord) and left that way, so an edge cut by a seam
    /// reads as a border: CountUVEdgeNum needs it so, and a caller that wants the adjacency of
    /// the surface back must rebuild it with UpdateTopology::FaceFace.
    /// With per-vertex texture coordinates a seam already splits the vertices, so the plain
    /// face-face adjacency partitions the mesh as it is and no UV function is needed.
    static int PerWedgeCharts(MeshType &m, std::vector<int> &faceChart)
    {
        tri::RequireFFAdjacency(m);
        tri::RequirePerFaceWedgeTexCoord(m);
        UpdateTopology<MeshType>::FaceFaceFromTexCoord(m);
        std::vector< std::pair<int,FacePointer> > components;
        return Clean<MeshType>::ConnectedComponents(m, components, &faceChart);
    }

    /// Count the edges of the UV layout, and the ones on its boundary. An edge that a seam cuts
    /// is two edges of the layout, one for each side, and both are on its boundary, as are the
    /// edges on a border of the mesh. An edge shared by more than two faces counts once.
    /// Needs the face-face adjacency left by PerWedgeCharts.
    static void CountUVEdgeNum(MeshType &m, int &total_e, int &boundary_e)
    {
        tri::RequireFFAdjacency(m);
        total_e=0;
        boundary_e=0;
        for (FaceIterator fi=m.face.begin();fi!=m.face.end();++fi)
        {
            if ((*fi).IsD()) continue;
            for (int i=0;i<3;i++)
            {
                if (face::IsBorder(*fi,i))
                {
                    total_e++;
                    boundary_e++;
                    continue;
                }
                // The faces around an edge form a ring: the edge is counted from the first of them.
                const FacePointer first=&*fi;
                FacePointer g=first;
                int j=i;
                bool isFirst=true;
                do
                {
                    const FacePointer next=g->FFp(j);
                    j=g->FFi(j);
                    g=next;
                    if (g<first || (g==first && j<i)) isFirst=false;
                } while (isFirst && (g!=first || j!=i));
                if (isFirst) total_e++;
            }
        }
    }
};

} // end Namespace tri
} // end Namespace vcg
#endif // VCG_COMPLEX_ALGORITHMS_PARAMETRIZATION_UV_UTILS_H
