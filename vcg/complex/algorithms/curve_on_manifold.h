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
#ifndef __VCGLIB_CURVE_ON_SURF_H
#define __VCGLIB_CURVE_ON_SURF_H

#include<vcg/complex/complex.h>
#include<vcg/simplex/face/topology.h>
#include<vcg/complex/algorithms/update/topology.h>
#include<vcg/complex/algorithms/update/color.h>
#include<vcg/complex/algorithms/update/normal.h>
#include<vcg/complex/algorithms/update/quality.h>
#include<vcg/complex/algorithms/clean.h>
#include<vcg/complex/algorithms/refine.h>
#include<vcg/complex/algorithms/create/platonic.h>
#include<vcg/complex/algorithms/point_sampling.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/space/index/kdtree/kdtree.h>
#include <vcg/math/histogram.h>
#include<vcg/space/distance3.h>
#include <vcg/complex/algorithms/attribute_seam.h>
#include <wrap/io_trimesh/export_ply.h>

namespace vcg {
namespace tri {
/// \ingroup trimesh
/// \brief A class for managing curves on a 2manifold.
/**
  This class is used to project/simplify/smooth polylines over a triangulated surface. 
  
*/

template <class MeshType>
class CoM
{
public:
  typedef typename MeshType::ScalarType     ScalarType;
  typedef typename MeshType::CoordType      CoordType;
  typedef typename MeshType::VertexType     VertexType;
  typedef typename MeshType::VertexPointer  VertexPointer;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::EdgeIterator   EdgeIterator;
  typedef typename MeshType::EdgeType       EdgeType;
  typedef typename MeshType::FaceType       FaceType;
  typedef typename MeshType::FacePointer    FacePointer;
  typedef typename MeshType::FaceIterator   FaceIterator;
  typedef Box3<ScalarType>                  Box3Type;
  typedef Segment3<ScalarType>              Segment3Type;  
  typedef typename vcg::GridStaticPtr<FaceType, ScalarType> MeshGrid;  
  typedef typename vcg::GridStaticPtr<EdgeType, ScalarType> EdgeGrid;
  typedef typename face::Pos<FaceType> PosType;
  typedef typename tri::UpdateTopology<MeshType>::PEdge PEdge;
  class Param 
  {
  public:
    
    ScalarType surfDistThr; // Distance between surface and curve; used in simplify and refine
    ScalarType polyDistThr; // Distance between the 
    ScalarType minRefEdgeLen;  // Minimal admitted Edge Lenght (used in refine: never make edge shorther than this value) 
    ScalarType maxSimpEdgeLen; // Maximal admitted Edge Lenght (used in simplify: never make edges longer than this value) 
    ScalarType maxSmoothDelta; // The maximum movement that is admitted during smoothing.
    ScalarType maxSnapThr;     // The maximum distance allowed when snapping a vertex of the polyline onto a mesh vertex
    ScalarType gridBailout;    // The maximum distance bailout used in grid sampling
    ScalarType barycentricSnapThr;    // The maximum distance bailout used in grid sampling
    
    Param(MeshType &m) { Default(m);}
    
    void Default(MeshType &m)
    {
      surfDistThr   = m.bbox.Diag()/1000.0;
      polyDistThr   = m.bbox.Diag()/5000.0;
      minRefEdgeLen    = m.bbox.Diag()/16000.0;
      maxSimpEdgeLen    = m.bbox.Diag()/100.0;
      maxSmoothDelta =   m.bbox.Diag()/100.0;
      maxSnapThr =       m.bbox.Diag()/1000.0;
      gridBailout =      m.bbox.Diag()/20.0;
      barycentricSnapThr = 0.05;
    }
    void Dump() const
    {
      printf("surfDistThr    = %6.4f\n",surfDistThr   );
      printf("polyDistThr    = %6.4f\n",polyDistThr   );
      printf("minRefEdgeLen  = %6.4f\n",minRefEdgeLen    );
      printf("maxSimpEdgeLen = %6.4f\n",maxSimpEdgeLen    );
      printf("maxSmoothDelta = %6.4f\n",maxSmoothDelta);
    }
  };
  
  
  
  // The Data Members
  
  MeshType &base; 
  MeshGrid uniformGrid;
  
  Param par; 
  CoM(MeshType &_m) :base(_m),par(_m){}
 
  FaceType *GetClosestFace(const CoordType &p)
  {
    ScalarType closestDist;
    CoordType closestP;
    return vcg::tri::GetClosestFaceBase(base,uniformGrid,p, this->par.gridBailout, closestDist, closestP);
  }
  
  FaceType *GetClosestFaceIP(const CoordType &p, CoordType &ip)
    {
      ScalarType closestDist;
      CoordType closestP,closestN;
      return vcg::tri::GetClosestFaceBase(base,uniformGrid,p, this->par.gridBailout, closestDist, closestP,closestN,ip);
    }

  FaceType *GetClosestFaceIP(const CoordType &p, CoordType &ip, CoordType &in)
    {
      ScalarType closestDist;
      CoordType closestP;
      return vcg::tri::GetClosestFaceBase(base,uniformGrid,p, this->par.gridBailout, closestDist, closestP,in,ip);
    }

  FaceType *GetClosestFacePoint(const CoordType &p, CoordType &closestP)
  {
    ScalarType closestDist;
    return vcg::tri::GetClosestFaceBase(base,uniformGrid,p, this->par.gridBailout, closestDist, closestP);
  }
  
  bool IsSnappedEdge(CoordType &ip, int &ei)
  {
    for(int i=0;i<3;++i)
      if(ip[i]>0.0 && ip[(i+1)%3]>0.0 && ip[(i+2)%3]==0.0 ) {
        ei=i;
        return true; 
      }
    ei=-1;
    return false;
  }

  // Given a baricentric coordinate finds that we assume that snaps onto an edge, it finds the vertex on which it is snapping 
  bool IsSnappedVertex(CoordType &ip, int &vi)
  {
    for(int i=0;i<3;++i)
      if(ip[i]==1.0 && ip[(i+1)%3]==0.0 && ip[(i+2)%3]==0.0 ) {
        vi=i;
        return true; 
      }
    vi=-1;
    return false;
  }

  // Given a baricentric coordinate finds that we assume that snaps onto an edge, it finds the vertex on which it is snapping 
  VertexPointer FindVertexSnap(FacePointer fp, CoordType &ip)
  {
    for(int i=0;i<3;++i)
      if(ip[i]==1.0 && ip[(i+1)%3]==0.0 && ip[(i+2)%3]==0.0 ) return fp->V(i);
    return 0;
  }
  
  
  
  /**
   * @brief TagFaceEdgeSelWithPolyLine selects edges of basemesh when they coincide with the polyline ones  *
   * @param poly
   * @return true if all the edges of the polyline are snapped onto the mesh. 
   * 
   * Use this function together with the CutMeshAlongCrease function to actually cut the mesh with a snapped polyline.
   * 
   */
    
bool TagFaceEdgeSelWithPolyLine(MeshType &poly,bool markFlag=true)
{
	if (markFlag)
		tri::UpdateFlags<MeshType>::FaceClearFaceEdgeS(base);

	tri::UpdateTopology<MeshType>::VertexFace(base);
	tri::UpdateTopology<MeshType>::FaceFace(base);

	for(EdgeIterator ei=poly.edge.begin(); ei!=poly.edge.end();++ei)
	{
		CoordType ip0,ip1;
		FaceType *f0 = GetClosestFaceIP(ei->cP(0),ip0);
		FaceType *f1 = GetClosestFaceIP(ei->cP(1),ip1);

		if(BarycentricSnap(ip0) && BarycentricSnap(ip1))
		{
			VertexPointer v0 = FindVertexSnap(f0,ip0);
			VertexPointer v1 = FindVertexSnap(f1,ip1);

			if(v0==0 || v1==0)
				return false;
			if(v0==v1)
				return false;

			FacePointer ff0,ff1;
			int e0,e1;
			bool ret=face::FindSharedFaces<FaceType>(v0,v1,ff0,ff1,e0,e1);
			if(ret)
			{
				assert(ret);
				assert(ff0->V(e0)==v0 || ff0->V(e0)==v1);
				ff0->SetFaceEdgeS(e0);
				ff1->SetFaceEdgeS(e1);
			} else {
				return false;
			}
		}
		else {
			return false;
		}
	}
	return true;
}

  ScalarType MinDistOnEdge(CoordType samplePnt, EdgeGrid &edgeGrid, MeshType &poly, CoordType &closestPoint)
  {
      ScalarType polyDist;
      EdgeType *cep = vcg::tri::GetClosestEdgeBase(poly,edgeGrid,samplePnt,par.gridBailout,polyDist,closestPoint);        
      return polyDist;    
  }
  
  // Given an edge of a mesh, supposedly intersecting the polyline, 
  // we search on it the closest point to the polyline
  static ScalarType MinDistOnEdge(VertexType *v0,VertexType *v1, EdgeGrid &edgeGrid, MeshType &poly, CoordType &closestPoint)
  {
    ScalarType minPolyDist = std::numeric_limits<ScalarType>::max();
    const ScalarType sampleNum = 50;
    const ScalarType maxDist = poly.bbox.Diag()/10.0;
    for(ScalarType k = 0;k<sampleNum+1;++k)
    {
      ScalarType polyDist;
      CoordType closestPPoly;
      CoordType samplePnt = (v0->P()*k +v1->P()*(sampleNum-k))/sampleNum;          
      
      EdgeType *cep = vcg::tri::GetClosestEdgeBase(poly,edgeGrid,samplePnt,maxDist,polyDist,closestPPoly);        
      
      if(polyDist < minPolyDist)
      {
        minPolyDist = polyDist;
        closestPoint = samplePnt;
//        closestPoint = closestPPoly;
      }
    }
    return minPolyDist;    
  }
  
  
  /**
   * @brief ExtractVertex
   * must extract an unambiguous representation of a vertex 
   * to be used with attribute_seam.h
   * 
   */
  static inline void ExtractVertex(const MeshType & srcMesh, const FaceType & f, int whichWedge, const MeshType & dstMesh, VertexType & v)
  {
      (void)srcMesh;
      (void)dstMesh;
      // This is done to preserve every single perVertex property
      // perVextex Texture Coordinate is instead obtained from perWedge one.
      v.ImportData(*f.cV(whichWedge));
      v.C() = f.cC();
  }
  
  static inline bool CompareVertex(const MeshType & m, const VertexType & vA, const VertexType & vB)
  {
      (void)m;
      
      if(vA.C() == Color4b(Color4b::Red) && vB.C() == Color4b(Color4b::Blue) ) return false;
      if(vA.C() == Color4b(Color4b::Blue) && vB.C() == Color4b(Color4b::Red) ) return false;
      return true;      
  }
  
  
  
  static CoordType QLerp(VertexType *v0, VertexType *v1)
  {
    
    ScalarType qSum = fabs(v0->Q())+fabs(v1->Q());      
    ScalarType w0 = (qSum - fabs(v0->Q()))/qSum;
    ScalarType w1 = (qSum - fabs(v1->Q()))/qSum;      
    return v0->P()*w0 + v1->P()*w1;      
  }
  
  
  /**
   * @brief SnapPolyline snaps the vertexes of a polyline onto the base mesh
   * @param poly
   * @param newVertVec the vector of the indexes of the snapped vertices
   * 
   * Polyline vertices can be snapped either on vertexes or on edges. 
   * Usually the only points that we should allow to not be snapped are the endpoints and non manifold points.
   * Vertexes are colored according to their snapping state 
   */  
    
  void SnapPolyline(MeshType &poly)
  {
    tri::Allocator<MeshType>::CompactEveryVector(poly);     
    tri::UpdateTopology<MeshType>::VertexEdge(poly);
    int vertSnapCnt=0;
    int edgeSnapCnt=0;
    int borderCnt=0,midCnt=0,nonmanifCnt=0;
    for(VertexIterator vi=poly.vert.begin(); vi!=poly.vert.end();++vi)
    {
      CoordType ip;
      FaceType *f = GetClosestFaceIP(vi->cP(),ip);
      if(BarycentricSnap(ip))
      {
        if(ip[0]>0 && ip[1]>0) { vi->P() = f->P(0)*ip[0]+f->P(1)*ip[1]; edgeSnapCnt++; assert(ip[2]==0); vi->C()=Color4b::White;}
        if(ip[0]>0 && ip[2]>0) { vi->P() = f->P(0)*ip[0]+f->P(2)*ip[2]; edgeSnapCnt++; assert(ip[1]==0); vi->C()=Color4b::White;}
        if(ip[1]>0 && ip[2]>0) { vi->P() = f->P(1)*ip[1]+f->P(2)*ip[2]; edgeSnapCnt++; assert(ip[0]==0); vi->C()=Color4b::White;}
        
        if(ip[0]==1.0) { vi->P() = f->P(0); vertSnapCnt++; assert(ip[1]==0 && ip[2]==0); vi->C()=Color4b::Black;  }
        if(ip[1]==1.0) { vi->P() = f->P(1); vertSnapCnt++; assert(ip[0]==0 && ip[2]==0); vi->C()=Color4b::Black;}
        if(ip[2]==1.0) { vi->P() = f->P(2); vertSnapCnt++; assert(ip[0]==0 && ip[1]==0); vi->C()=Color4b::Black;}
      }
      else
      {
        int deg = edge::VEDegree<EdgeType>(&*vi);
        if (deg > 2) { nonmanifCnt++; vi->C()=Color4b::Magenta; }
        if (deg < 2) { borderCnt++;   vi->C()=Color4b::Green;} 
        if (deg== 2) { midCnt++;      vi->C()=Color4b::Blue;} 
      }
    }
    printf("SnapPolyline %i vertices:  snapped %i onto vert and %i onto edges %i nonmanif, %i border, %i mid\n",
           poly.vn, vertSnapCnt, edgeSnapCnt, nonmanifCnt,borderCnt,midCnt); fflush(stdout);
    int dupCnt=tri::Clean<MeshType>::RemoveDuplicateVertex(poly);
    tri::Allocator<MeshType>::CompactEveryVector(poly);     
    if(dupCnt) printf("SnapPolyline: Removed %i Duplicated vertices\n",dupCnt);
  }
  
   void SelectBoundaryVertex(MeshType &poly)
   {
     tri::UpdateSelection<MeshType>::VertexClear(poly);
     tri::UpdateTopology<MeshType>::VertexEdge(poly);
     ForEachVertex(poly, [&](VertexType &v){
       if(edge::VEDegree<EdgeType>(&v)==1) v.SetS();
     });
   }
  
   void SelectUniformlyDistributed(MeshType &poly, int k)
   {
     tri::TrivialPointerSampler<MeshType> tps;
     ScalarType samplingRadius = tri::Stat<MeshType>::ComputeEdgeLengthSum(poly)/ScalarType(k);
     tri::SurfaceSampling<MeshType, typename tri::TrivialPointerSampler<MeshType> >::EdgeMeshUniform(poly,tps,samplingRadius);     
     for(int i=0;i<tps.sampleVec.size();++i)
       tps.sampleVec[i]->SetS();
   }
   
   
    
  /*
   * Make an edge mesh 1-manifold by splitting all the
   * vertexes that have more than two incident edges
   * 
   * It performs the split in three steps. 
   * - First it collects and counts the vertices to be splitten. 
   * - Then it adds the vertices to the mesh and 
   * - lastly it updates the poly with the newly added vertices. 
   * 
   * singSplitFlag allows to ubersplit each singularity in a number of vertex of the same order of its degree. 
   * This is not really necessary but helps the management of sharp turns in the poly mesh.
   * \todo add corner detection and split.
   */
  
  void DecomposeNonManifoldPolyline(MeshType &poly, bool singSplitFlag = true)
  {
    tri::Allocator<MeshType>::CompactEveryVector(poly);
    std::vector<int> degreeVec(poly.vn, 0);
    tri::UpdateTopology<MeshType>::VertexEdge(poly);
    int neededVert=0;
    int delta;
    if(singSplitFlag) delta = 1;
                 else delta = 2;
      
    for(VertexIterator vi=poly.vert.begin(); vi!=poly.vert.end();++vi)
    {
      std::vector<EdgeType *> starVec;
      edge::VEStarVE(&*vi,starVec);
      degreeVec[tri::Index(poly, *vi)] = starVec.size();
      if(starVec.size()>2)
        neededVert += starVec.size()-delta;
    }
    printf("DecomposeNonManifold Adding %i vert to a polyline of %i vert\n",neededVert,poly.vn);
    VertexIterator firstVi = tri::Allocator<MeshType>::AddVertices(poly,neededVert);
    
    for(size_t i=0;i<degreeVec.size();++i)
    {
      if(degreeVec[i]>2)
      {
        std::vector<EdgeType *> edgeStarVec;
        edge::VEStarVE(&(poly.vert[i]),edgeStarVec);
        assert(edgeStarVec.size() == degreeVec[i]);
        for(size_t j=delta;j<edgeStarVec.size();++j)
        {
          EdgeType *ep = edgeStarVec[j];
          int ind; // index of the vertex to be changed
          if(tri::Index(poly,ep->V(0)) == i) ind = 0;
              else ind = 1;
  
          ep->V(ind) = &*firstVi;
          ep->V(ind)->P() = poly.vert[i].P();
          ep->V(ind)->N() = poly.vert[i].N();
          ++firstVi;
        }
      }
    }
    assert(firstVi == poly.vert.end());
  }
  

  
   
  /**
   * @brief SplitMeshWithPolyline
   * @param poly
   * 
   * First it splits the base mesh with all the non snapped points doing a standard 1 to 3 split;
   * 
   */
  void SplitMeshWithPolyline(MeshType &poly)
  {        
    std::vector< std::pair<int,VertexPointer> > toSplitVec;  // the index of the face to be split and the poly vertex to be used
    
    for(VertexIterator vi=poly.vert.begin(); vi!=poly.vert.end();++vi)
    {
      CoordType ip;
      FaceType *f = GetClosestFaceIP(vi->cP(),ip);
      if(!BarycentricSnap(ip))
        toSplitVec.push_back(std::make_pair(tri::Index(base,f),&*vi));
    }
    
    printf("SplitMeshWithPolyline found %lu non snapped points\n",toSplitVec.size());fflush(stdout);

    FaceIterator newFi = tri::Allocator<MeshType>::AddFaces(base,toSplitVec.size()*2);
    VertexIterator newVi = tri::Allocator<MeshType>::AddVertices(base,toSplitVec.size());    
    tri::UpdateColor<MeshType>::PerVertexConstant(base,Color4b::White);
    
    for(size_t i =0; i<toSplitVec.size();++i)
    {
        newVi->P() = toSplitVec[i].second->P(); 
        newVi->C()=Color4b::Green;      
        face::TriSplit(&base.face[toSplitVec[i].first],&*(newFi++),&*(newFi++),&*(newVi++));
    }
    Init(); //  need to reset everthing
    SnapPolyline(poly);
    
    // Second loop to perform the face-face Edge split **********************

    std::map<std::pair<CoordType,CoordType>, VertexPointer> edgeToPolyVertMap;
    for(VertexIterator vi=poly.vert.begin(); vi!=poly.vert.end();++vi)
    {
      CoordType ip;
      FaceType *f = GetClosestFaceIP(vi->cP(),ip);
      if(!BarycentricSnap(ip)) { assert(0); }            
      for(int i=0;i<3;++i)
      {
       if(ip[i]>0 && ip[(i+1)%3]>0 && ip[(i+2)%3]==0 )
       {
        CoordType p0=f->P0(i);
        CoordType p1=f->P1(i);
        if (p0>p1) std::swap(p0,p1);  
        if(edgeToPolyVertMap[make_pair(p0,p1)]) printf("Found an already used Edge %lu - %lu %lu!!!\n", tri::Index(base,f->V0(i)),tri::Index(base,f->V1(i)),tri::Index(poly,&*vi));
        edgeToPolyVertMap[make_pair(p0,p1)]=&*vi;
       }         
      }
    }
    printf("SplitMeshWithPolyline: Created a map of %lu edges to be split\n",edgeToPolyVertMap.size());
    EdgePointPred ePred(edgeToPolyVertMap);
    EdgePointSplit eSplit(edgeToPolyVertMap);
    tri::UpdateTopology<MeshType>::FaceFace(base);
    tri::RefineE(base,eSplit,ePred);     
    Init(); //  need to reset everthing    
 }
    
        
  
  void Init()
  {
    // Construction of the uniform grid
    UpdateNormal<MeshType>::PerFaceNormalized(base);
    UpdateTopology<MeshType>::FaceFace(base);    
    uniformGrid.Set(base.face.begin(), base.face.end());    
  }
  
  
  void SimplifyMidEdge(MeshType &poly)
  {
   int startVn;
   int midEdgeCollapseCnt=0;
   tri::Allocator<MeshType>::CompactEveryVector(poly); 
   do
   {
    startVn = poly.vn;
    for(int ei =0; ei<poly.en; ++ei)
    {
      VertexType *v0=poly.edge[ei].V(0);
      VertexType *v1=poly.edge[ei].V(1);
      CoordType ip0,ip1;    
      FaceType *f0=GetClosestFaceIP(v0->P(),ip0);
      FaceType *f1=GetClosestFaceIP(v1->P(),ip1);
      
      bool snap0=BarycentricSnap(ip0);
      bool snap1=BarycentricSnap(ip1);
      int e0i,e1i;
      bool e0 = IsSnappedEdge(ip0,e0i);
      bool e1 = IsSnappedEdge(ip1,e1i);
      if(e0 && e1)
        if( (          f0 == f1           &&          e0i == e1i) || 
            (          f0 == f1->FFp(e1i) &&          e0i == f1->FFi(e1i)) || 
            (f0->FFp(e0i) == f1           && f0->FFi(e0i) == e1i) || 
            (f0->FFp(e0i) == f1->FFp(e1i) && f0->FFi(e0i) == f1->FFi(e1i)) ) 
        {
          CoordType newp = (v0->P()+v1->P())/2.0;
          v0->P()=newp;
          v1->P()=newp;
          midEdgeCollapseCnt++;
        }
    }
    tri::Clean<MeshType>::RemoveDuplicateVertex(poly);
    tri::Allocator<MeshType>::CompactEveryVector(poly);     
//    printf("SimplifyMidEdge %5i -> %5i %i mid %i ve \n",startVn,poly.vn,midEdgeCollapseCnt);
   } while(startVn>poly.vn);
  } 
  
  /**
   * @brief SimplifyMidFace remove all the vertices that in the mid of a face 
   * and between two of the points snapped onto the edges of the same face
   * @param poly
   * 
   * It assumes that the mesh has been snapped and refined by the BaseMesh
   * 
   */
  void SimplifyMidFace(MeshType &poly)
  {
   int startVn= poly.vn;;
   int midFaceCollapseCnt=0;
   int vertexEdgeCollapseCnt=0;
   int curVn;
   do
   {
    tri::Allocator<MeshType>::CompactEveryVector(poly); 
    curVn = poly.vn;
    UpdateTopology<MeshType>::VertexEdge(poly);
    for(int i =0; i<poly.vn;++i)
    {
      std::vector<VertexPointer> starVecVp;
      edge::VVStarVE(&(poly.vert[i]),starVecVp);      
      if( (starVecVp.size()==2) )
      {
        CoordType ipP, ipN, ipI; 
        FacePointer fpP = GetClosestFaceIP(starVecVp[0]->P(),ipP);
        FacePointer fpN = GetClosestFaceIP(starVecVp[1]->P(),ipN);
        FacePointer fpI = GetClosestFaceIP(poly.vert[i].P(), ipI);
        
        bool snapP = (BarycentricSnap(ipP));
        bool snapN = (BarycentricSnap(ipN));
        bool snapI = (BarycentricSnap(ipI));
        VertexPointer vertexSnapP = 0;
        VertexPointer vertexSnapN = 0;
        VertexPointer vertexSnapI = 0;
        for(int j=0;j<3;++j)
        {
          if(ipP[j]==1.0) vertexSnapP=fpP->V(j);
          if(ipN[j]==1.0) vertexSnapN=fpN->V(j);
          if(ipI[j]==1.0) vertexSnapI=fpI->V(j);
        }
        
        bool collapseFlag=false;
        
        if((!snapI && snapP && snapN) ||              // First case a vertex that is not snapped between two snapped vertexes 
           (!snapI && !snapP && fpI==fpP) || // Or a two vertex not snapped but on the same face
           (!snapI && !snapN && fpI==fpN) )
        {
          collapseFlag=true;
          midFaceCollapseCnt++;
        } 
        
        else  // case 2) a vertex snap and edge snap we have to check that the edge do not share the same vertex of the vertex snap
          if(snapI && snapP && snapN && vertexSnapI==0 && (vertexSnapP!=0 || vertexSnapN!=0) )
          {
            for(int j=0;j<3;++j) {
              if(ipI[j]!=0 && (fpI->V(j)==vertexSnapP || fpI->V(j)==vertexSnapN)) {
                collapseFlag=true;                                          
                vertexEdgeCollapseCnt++;
              }
            }
          }            
        
        if(collapseFlag)  
          edge::VEEdgeCollapse(poly,&(poly.vert[i]));
      }
    }  
   } while(curVn>poly.vn);
   printf("SimplifyMidFace %5i -> %5i %i mid %i ve \n",startVn,poly.vn,midFaceCollapseCnt,vertexEdgeCollapseCnt);
  } 
  
  void Simplify(MeshType &poly)
  {
    int startEn = poly.en;
    Distribution<ScalarType> hist;
    for(int i =0; i<poly.en;++i) 
      hist.Add(edge::Length(poly.edge[i]));
        
    UpdateTopology<MeshType>::VertexEdge(poly);
    
    for(int i =0; i<poly.vn;++i)
    {
      std::vector<VertexPointer> starVecVp;
      edge::VVStarVE(&(poly.vert[i]),starVecVp);      
      if ((starVecVp.size()==2) && (!poly.vert[i].IsS()))
      {
        ScalarType newSegLen = Distance(starVecVp[0]->P(), starVecVp[1]->P());
        Segment3Type seg(starVecVp[0]->P(),starVecVp[1]->P());
        ScalarType segDist;
        CoordType closestPSeg;
        SegmentPointDistance(seg,poly.vert[i].cP(),closestPSeg,segDist);
        CoordType fp,fn;
        ScalarType maxSurfDist = MaxSegDist(starVecVp[0], starVecVp[1],fp,fn);
        
        if((maxSurfDist < par.surfDistThr) && (newSegLen < par.maxSimpEdgeLen) )
        {
          edge::VEEdgeCollapse(poly,&(poly.vert[i]));          
        }
      }
    }
    tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
    tri::Allocator<MeshType>::CompactEveryVector(poly);
    tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
//    printf("Simplify %5i -> %5i (total len %5.2f)\n",startEn,poly.en,hist.Sum());
  }
  
  void EvaluateHausdorffDistance(MeshType &poly, Distribution<ScalarType> &dist)
  {
    dist.Clear();
    tri::UpdateTopology<MeshType>::VertexEdge(poly);
    tri::UpdateQuality<MeshType>::VertexConstant(poly,0);
    for(int i =0; i<poly.edge.size();++i)
    {      
      CoordType farthestP, farthestN;      
      ScalarType maxDist = MaxSegDist(poly.edge[i].V(0),poly.edge[i].V(1), farthestP, farthestN, &dist);      
      poly.edge[i].V(0)->Q()+= maxDist;
      poly.edge[i].V(1)->Q()+= maxDist;
    }
    for(int i=0;i<poly.vn;++i)
    {
      ScalarType deg = edge::VEDegree<EdgeType>(&poly.vert[i]);
      poly.vert[i].Q()/=deg;
    }
    tri::UpdateColor<MeshType>::PerVertexQualityRamp(poly,0,dist.Max());    
  }
  

  
  /**
   * @brief BarycentricSnap
   * @param ip the baricentric coords to be snapped
   * @return true if they have been snapped.
   * 
   * This is the VERY important function that is used everywhere. 
   * Given a barycentric coord of a point inside a triangle it decides if it should be "snapped" either onto an edge or on a vertex. 
   * It relies on the  barycentricSnapThr parameter
   * 
   */
  bool BarycentricSnap(CoordType &ip)
  {
    for(int i=0;i<3;++i)
    {
      if(ip[i] <= par.barycentricSnapThr) ip[i]=0;
      if(ip[i] >= 1.0-par.barycentricSnapThr) ip[i]=1;
    }
    ScalarType sum = ip[0]+ip[1]+ip[2];
    
    for(int i=0;i<3;++i) 
      if(ip[i]!=1) ip[i]/=sum;
    
    if(ip[0]==0 || ip[1]==0 || ip[2]==0) return true;
    return false;
  }
  
  
  /**
   * @brief TestSplitSegWithMesh  Given a poly segment decide if it should be split along elements of base mesh. 
   * @param v0
   * @param v1
   * @param splitPt
   * @return true if it should be split
   * 
   * We make a few samples onto the edge and if some of them snaps onto a an edge we use it.
   * In case there are more than one candidate we choose the sample closeset to its snapping point.
   * We explicitly avoid snapping twice on the same edge by checking the starting and ending edges.
   * 
   * Two cases:
   * - poly edge pass near a vertex of the mesh
   * - poly edge cross one or more edges
   * 
   * Note that we have to check the case where 
   */
  bool TestSplitSegWithMesh(VertexType *v0, VertexType *v1, CoordType &splitPt)
  {
    Segment3Type segPoly(v0->P(),v1->P());
    const ScalarType sampleNum = 40;    
    CoordType ip0,ip1;
    
    FaceType *f0=GetClosestFaceIP(v0->P(),ip0);
    FaceType *f1=GetClosestFaceIP(v1->P(),ip1);
    if(f0==f1) return false;
    
    bool snap0=false,snap1=false; // true if the segment start/end on a edge/vert
    
    Segment3Type seg0; // The two segments to be avoided 
    Segment3Type seg1; // from which the current poly segment can start
    VertexPointer vertexSnap0 = 0;
    VertexPointer vertexSnap1 = 0;
    if(BarycentricSnap(ip0)) { 
      snap0=true; 
      for(int i=0;i<3;++i) {
        if(ip0[i]==1.0) vertexSnap0=f0->V(i);
        if(ip0[i]==0.0) seg0=Segment3Type(f0->P1(i),f0->P2(i)); 
      }        
    } 
    if(BarycentricSnap(ip1)) { 
      snap1=true; 
      for(int i=0;i<3;++i){
        if(ip1[i]==1.0) vertexSnap1=f1->V(i);
        if(ip1[i]==0.0) seg1=Segment3Type(f1->P1(i),f1->P2(i)); 
      }        
    } 
    
    CoordType bestSplitPt(0,0,0);
    ScalarType bestDist = std::numeric_limits<ScalarType>::max();
    for(ScalarType k = 1;k<sampleNum;++k)
    {
      CoordType samplePnt = segPoly.Lerp(k/sampleNum);    
      CoordType ip;
      FaceType *f=GetClosestFaceIP(samplePnt,ip);
//      BarycentricEdgeSnap(ip);
      if(BarycentricSnap(ip))
      {
        VertexPointer vertexSnapI = 0;        
        for(int i=0;i<3;++i)
          if(ip[i]==1.0) vertexSnapI=f->V(i);
        CoordType closestPt = f->P(0)*ip[0]+f->P(1)*ip[1]+f->P(2)*ip[2];
        if(Distance(samplePnt,closestPt) < bestDist )  
        {
          ScalarType dist0=std::numeric_limits<ScalarType>::max();
          ScalarType dist1=std::numeric_limits<ScalarType>::max();
          CoordType closestSegPt;
          if(snap0) SegmentPointDistance(seg0,closestPt,closestSegPt,dist0);
          if(snap1) SegmentPointDistance(seg1,closestPt,closestSegPt,dist1);
          if( (!vertexSnapI && (dist0 > par.surfDistThr/1000 && dist1>par.surfDistThr/1000) ) ||
              ( vertexSnapI!=vertexSnap0 && vertexSnapI!=vertexSnap1)  )
          {
            bestDist = Distance(samplePnt,closestPt);
            bestSplitPt = closestPt;            
          }
        }      
      }
    }
    if(bestDist < par.surfDistThr*100)
    {
      splitPt = bestSplitPt;
      return true;
    }
    
    return false;
  }
  /**
   * @brief SnappedOnSameFace Return true if the two points are snapped to a common face;
   * @param f0
   * @param i0
   * @param f1
   * @param i0
   * @return 
   * 
   * Require FFAdj. se assume that both SNAPPED. Three cases:
   * - Edge Edge - true iff the two edges belongs to a common face. 
   * - Vert Edge - true iff there is one of the two snapped edge faces has the vert as non-edge face;  
   * - Vert Vert 
   * 
   */
  bool SnappedOnSameFace(FacePointer f0, CoordType i0, FacePointer f1, CoordType i1)
  {
   if(f0==f1) return true;
   int e0,e1;
   int v0,v1;
   bool e0Snap = IsSnappedEdge(i0,e0);
   bool e1Snap = IsSnappedEdge(i1,e1);
   bool v0Snap = IsSnappedVertex(i0,v0);
   bool v1Snap = IsSnappedVertex(i1,v1);
   FacePointer f0p=0; int e0p=-1;  // When Edge snap the other face and the index of the snapped edge on the other face
   FacePointer f1p=0; int e1p=-1;
   assert((e0Snap != v0Snap) && (e1Snap != v1Snap));
   // For EdgeSnap compute the 'other' face stuff 
   if(e0Snap){
     f0p = f0->FFp(e0); e0p=f0->FFi(e0); assert(f0p->FFp(e0p)==f0);
   }
   if(e1Snap){
     f1p = f1->FFp(e1); e1p=f1->FFi(e1); assert(f1p->FFp(e1p)==f1);
   }
   
   if(e0Snap && e1Snap) {
    if(f0==f1p || f0p==f1p || f0p==f1 || f0==f1) return true;
   }
   
   if(e0Snap && v1Snap)  {
     assert(v1>=0 && v1<3 && v0==-1 && e1==-1);
     if(f0->V2(e0)  ==f1->V(v1)) return true;
     if(f0p->V2(e0p)==f1->V(v1)) return true;
   }
     
   if(e1Snap && v0Snap)  {
     assert(v0>=0 && v0<3 && v1==-1 && e0==-1);
     if(f1->V2(e1)  ==f0->V(v0)) return true;
     if(f1p->V2(e1p)==f0->V(v0)) return true;
   }
     
   if(v1Snap && v0Snap)  {
     PosType startPos(f0,f0->V(v0));
     PosType curPos=startPos;
     do
     {
       assert(curPos.V()==f0->V(v0));
       if(curPos.VFlip()==f1->V(v1)) return true;
       curPos.FlipE();
       curPos.FlipF();       
     }
     while(curPos!=startPos);   
   }   
   return false;    
  }
  
  /**
   * @brief TestSplitSegWithMesh  Given a poly segment decide if it should be split along elements of base mesh. 
   * @param v0
   * @param v1
   * @param splitPt
   * @return true if it should be split
   * 
   * We make a few samples onto the edge and if some of them snaps onto a an edge we use it.
   * In case there are more than one candidate we choose the sample closeset to its snapping point.
   * We explicitly avoid snapping twice on the same edge by checking the starting and ending edges.
   * 
   * Two cases:
   * - poly edge pass near a vertex of the mesh
   * - poly edge cross one or more edges
   * 
   * Note that we have to check the case where 
   */
  bool TestSplitSegWithMeshAdapt(VertexType *v0, VertexType *v1, CoordType &splitPt)
  {
    splitPt=(v0->P()+v1->P())/2.0;
      
    CoordType ip0,ip1,ipm;    
    FaceType *f0=GetClosestFaceIP(v0->P(),ip0);
    FaceType *f1=GetClosestFaceIP(v1->P(),ip1);
    FaceType *fm=GetClosestFaceIP(splitPt,ipm);
    
    if(f0==f1) return false;
    
    bool snap0=BarycentricSnap(ip0);
    bool snap1=BarycentricSnap(ip1);
    bool snapm=BarycentricSnap(ipm);
    
    splitPt = fm->P(0)*ipm[0]+fm->P(1)*ipm[1]+fm->P(2)*ipm[2];
    
    if(!snap0 && !snap1) {
      assert(f0!=f1);
      return true;
    }
    if(snap0 && snap1) 
    {
      if(SnappedOnSameFace(f0,ip0,f1,ip1)) 
        return false;            
    }
    
    if(snap0) {
      int e0,v0;
      if (IsSnappedEdge(ip0,e0)) {
        if(f0->FFp(e0) == f1) return false;
      }
      if(IsSnappedVertex(ip0,v0)) {
        for(int i=0;i<3;++i) 
          if(f1->V(i)==f0->V(v0)) return false;
      }
    }
    if(snap1) {
      int e1,v1;
      if (IsSnappedEdge(ip1,e1)) {
        if(f1->FFp(e1) == f0) return false;
      }
      if(IsSnappedVertex(ip1,v1)) {
        for(int i=0;i<3;++i) 
          if(f0->V(i)==f1->V(v1)) return false;
      }
    }
    
    return true;
  }
  
  
  bool TestSplitSegWithMeshAdaptOld(VertexType *v0, VertexType *v1, CoordType &splitPt)
  {
    Segment3Type segPoly(v0->P(),v1->P());
    const ScalarType sampleNum = 40;    
    CoordType ip0,ip1;    
    FaceType *f0=GetClosestFaceIP(v0->P(),ip0);
    FaceType *f1=GetClosestFaceIP(v1->P(),ip1);
    if(f0==f1) return false;
    
    bool snap0=BarycentricSnap(ip0);
    bool snap1=BarycentricSnap(ip1);
    
    if(!snap0 && !snap1) {
      assert(f0!=f1);
      splitPt=(v0->P()+v1->P())/2.0;
      return true;
    }
    if(snap0 && snap1) 
    {
      if(SnappedOnSameFace(f0,ip0,f1,ip1)) 
        return false;      
    }
    
    if(snap0) {
      int e0,v0;
      if (IsSnappedEdge(ip0,e0)) {
        if(f0->FFp(e0) == f1) return false;
      }
      if(IsSnappedVertex(ip0,v0)) {
        for(int i=0;i<3;++i) 
          if(f1->V(i)==f0->V(v0)) return false;
      }
    }
    splitPt=(v0->P()+v1->P())/2.0;
    return true;
  }
  
  // Given a segment find the maximum distance from it to the original surface. 
  // It is used to evaluate the Haustdorff distance of a Segment from the mesh.
  ScalarType MaxSegDist(VertexType *v0, VertexType *v1, CoordType &farthestPointOnSurf, CoordType &farthestN, Distribution<ScalarType> *dist=0)
  {
    ScalarType maxSurfDist = 0;
    const ScalarType sampleNum = 10;
    const ScalarType maxDist = base.bbox.Diag()/10.0;
    for(ScalarType k = 1;k<sampleNum;++k)
    {
      ScalarType surfDist;
      CoordType closestPSurf;
      CoordType samplePnt = (v0->P()*k +v1->P()*(sampleNum-k))/sampleNum;          
      FaceType *f = vcg::tri::GetClosestFaceBase(base,uniformGrid,samplePnt,maxDist, surfDist, closestPSurf);        
      if(dist)
        dist->Add(surfDist);
      assert(f);
      if(surfDist > maxSurfDist)
      {
        maxSurfDist = surfDist;
        farthestPointOnSurf = closestPSurf;
        farthestN = f->N();
      }
    }
    return maxSurfDist;
  }
  
  
  /**
   * @brief RefineCurve
   * @param poly the curve to be refined
   * @param uniformFlag
   * 
   * Make one pass of refinement for all the edges of the curve that are distant from the basemesh
   * uses two parameters:
   * - par.minRefEdgeLen 
   * - par.surfDistThr
   */
    
  void RefineCurveByDistance(MeshType &poly)
  {
    tri::Allocator<MeshType>::CompactEveryVector(poly);    
    int startEdgeSize = poly.en;
    for(int i =0; i<startEdgeSize;++i)
    {
      EdgeType &ei = poly.edge[i];
      if(edge::Length(ei)>par.minRefEdgeLen)  
      {      
        CoordType farthestP, farthestN;
        ScalarType maxDist = MaxSegDist(ei.V(0),ei.V(1),farthestP, farthestN);
        if(maxDist > par.surfDistThr)  
        {
          edge::VEEdgeSplit(poly, &ei, farthestP, farthestN); 
        }
      }
    }
//    tri::Allocator<MeshType>::CompactEveryVector(poly);
//    printf("Refine %i -> %i\n",startEdgeSize,poly.en);fflush(stdout);
  }
  
  /**
   * @brief RefineCurveByBaseMesh
   * @param poly
   */
  
  void RefineCurveByBaseMesh(MeshType &poly)
  {
    tri::Allocator<MeshType>::CompactEveryVector(poly);    
    std::vector<int> edgeToRefineVec;
    for(int i=0; i<poly.en;++i) 
      edgeToRefineVec.push_back(i);
    int startEn=poly.en;  
    int iterCnt=0;
    while (!edgeToRefineVec.empty() && iterCnt<100) {
      iterCnt++;
      std::vector<int> edgeToRefineVecNext;
      for(int i=0; i<edgeToRefineVec.size();++i)
      {
        EdgeType &e = poly.edge[edgeToRefineVec[i]];
        CoordType splitPt;
        if(TestSplitSegWithMeshAdapt(e.V(0),e.V(1),splitPt))  
        {
          edge::VEEdgeSplit(poly, &e, splitPt); 
          edgeToRefineVecNext.push_back(edgeToRefineVec[i]);
          edgeToRefineVecNext.push_back(poly.en-1);
        } 
      }
      tri::Allocator<MeshType>::CompactEveryVector(poly);
      swap(edgeToRefineVecNext,edgeToRefineVec);
       printf("RefineCurveByBaseMesh %i en -> %i en\n",startEn,poly.en); fflush(stdout);    
    }
//    
    SimplifyMidFace(poly);
    SimplifyMidEdge(poly);
    SnapPolyline(poly);    
    printf("RefineCurveByBaseMesh %i en -> %i en\n",startEn,poly.en); fflush(stdout);    
  }
  
  
  
  /**
   * @brief SmoothProject
   * @param poly
   * @param iterNum
   * @param smoothWeight  [0..1] range;  
   * @param projectWeight [0..1] range;
   * 
   * The very important function to adapt a polyline onto the base mesh
   * The projection process must be done slowly to guarantee some empirical convergence...
   * For each iteration it choose a new position of each vertex of the polyline. 
   * The new position is a blend between the smoothed position, the closest point on the surface and the original position. 
   * You need a good balance...
   * after each iteration the polyline is refined and simplified. 
   */
  void SmoothProject(MeshType &poly, int iterNum, ScalarType smoothWeight, ScalarType projectWeight)
  {
    tri::RequireCompactness(poly);
    tri::UpdateTopology<MeshType>::VertexEdge(poly);
//    printf("SmoothProject: Selected vert num %i\n",tri::UpdateSelection<MeshType>::VertexCount(poly));
    assert(poly.en>0 && base.fn>0);
    for(int k=0;k<iterNum;++k)
    {
      if(k==iterNum-1) projectWeight=1; 
      
      std::vector<CoordType> posVec(poly.vn,CoordType(0,0,0));
      std::vector<int>     cntVec(poly.vn,0);
  
      for(int i =0; i<poly.en;++i)
      {
        for(int j=0;j<2;++j)
        {
          int vertInd = tri::Index(poly,poly.edge[i].V0(j));
          posVec[vertInd] += poly.edge[i].V1(j)->P();
          cntVec[vertInd] += 1;
        }
      }
      
      const ScalarType maxDist = base.bbox.Diag()/10.0; 
      for(int i=0; i<poly.vn; ++i)
        if(!poly.vert[i].IsS())
        {
          CoordType smoothPos = (poly.vert[i].P() + posVec[i])/ScalarType(cntVec[i]+1);
          
          CoordType newP = poly.vert[i].P()*(1.0-smoothWeight) + smoothPos *smoothWeight;
          
//          CoordType delta =  newP - poly.vert[i].P();
//          if(delta.Norm() > par.maxSmoothDelta) 
//          {
//            newP =  poly.vert[i].P() + ( delta / delta.Norm()) * maxDist*0.5;
//          }
          
          ScalarType minDist;
          CoordType closestP;
          FaceType *f = vcg::tri::GetClosestFaceBase(base,uniformGrid,newP,maxDist, minDist, closestP);
          assert(f);
          poly.vert[i].P() = newP*(1.0-projectWeight) +closestP*projectWeight;
          poly.vert[i].N() = f->N();
        }
      
      //      Refine(poly);      
      tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
      RefineCurveByDistance(poly);      
      tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
      Simplify(poly);
      tri::UpdateTopology<MeshType>::TestVertexEdge(poly);
      int dupVertNum = Clean<MeshType>::RemoveDuplicateVertex(poly);
      if(dupVertNum) {
//        printf("****REMOVED %i Duplicated\n",dupVertNum);
        tri::Allocator<MeshType>::CompactEveryVector(poly);
        tri::UpdateTopology<MeshType>::VertexEdge(poly);
      }
    }
  }  
   
  
class EdgePointPred
{
public:
  std::map<std::pair<CoordType,CoordType>, VertexPointer> &edgeToPolyVertMap;
    
    EdgePointPred(std::map<std::pair<CoordType,CoordType>, VertexPointer> &_edgeToPolyVertMap):edgeToPolyVertMap(_edgeToPolyVertMap){};
    bool operator()(face::Pos<FaceType> ep) const
    {
      CoordType p0 = ep.V()->P();
      CoordType p1 = ep.VFlip()->P();
      if (p0>p1) std::swap(p0,p1);        
      VertexPointer vp=edgeToPolyVertMap[make_pair(p0,p1)];
      return vp!=0;
    }
};

struct EdgePointSplit
{
public:
  std::map<std::pair<CoordType,CoordType>, VertexPointer> &edgeToPolyVertMap;
  
  EdgePointSplit(std::map<std::pair<CoordType,CoordType>, VertexPointer> &_edgeToPolyVertMap):edgeToPolyVertMap(_edgeToPolyVertMap){};
    void operator()(VertexType &nv, face::Pos<FaceType> ep)
    {
      CoordType p0 = ep.V()->P();
      CoordType p1 = ep.VFlip()->P();
      if (p0>p1) std::swap(p0,p1);        
      VertexPointer vp=edgeToPolyVertMap[make_pair(p0,p1)];
      assert(vp);
      nv.P()=vp->P();
      return;
    }
    Color4b WedgeInterp(Color4b &c0, Color4b &c1)
    {
        Color4b cc;
        cc.lerp(c0,c1,0.5f);
        return Color4b::Red;
    }
    TexCoord2f WedgeInterp(TexCoord2f &t0, TexCoord2f &t1)
    {
        TexCoord2f tmp;
        assert(t0.n()== t1.n());
        tmp.n()=t0.n();
        tmp.t()=(t0.t()+t1.t())/2.0;
        return tmp;
    }
};

};

} // end namespace tri
} // end namespace vcg

#endif // __VCGLIB_CURVE_ON_SURF_H
