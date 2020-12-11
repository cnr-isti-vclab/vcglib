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

#ifndef __PLYMC_H__
#define __PLYMC_H__

#ifndef WIN32
#define _int64 long long
#define __int64 long long
#define __cdecl
#endif

#include <cstdio>
#include <time.h>
#include <float.h>
#include <math.h>

#include <vcg/complex/complex.h>

#include <vcg/math/histogram.h>
#include <vcg/complex/algorithms/geodesic.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export_ply.h>
//#include <wrap/ply/plystuff.h>

#include <vcg/complex/algorithms/create/marching_cubes.h>
#include <vcg/complex/algorithms/create/mc_trivial_walker.h>

// local optimization
#include <vcg/complex/algorithms/local_optimization.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h>

#include <stdarg.h>
#include "volume.h"
#include "tri_edge_collapse_mc.h"
namespace vcg {
namespace tri {

// Simple prototype for later use...
template<class MeshType>
int MCSimplify( MeshType &m, float perc, bool preserveBB=true, vcg::CallBackPos *cb=0);


/** Surface Reconstruction
 * 
 *  To allow the managment of a very large set of meshes to be merged,
 *  it is templated on a MeshProvider class that is able to provide the meshes to be merged. 
 *  IT is the surface reconstrction algorithm that have been used for a long time inside the ISTI-Visual Computer Lab.
 *  It is mostly a variant of the Curless et al. e.g. a volumetric approach with some original weighting schemes,"
 *  a different expansion rule, and another approach to hole filling through volume dilation/relaxations.
 */

template < class SMesh, class MeshProvider>
class PlyMC
{
public:

  typedef typename SMesh::FaceIterator     SFaceIterator;
  typedef typename SMesh::VertexIterator   SVertexIterator;

  class MCVertex;
  class MCEdge;
  class MCFace;

  class MCUsedTypes: public vcg::UsedTypes <  vcg::Use<MCVertex>::template AsVertexType,
                                              vcg::Use<MCEdge  >::template AsEdgeType,
                                              vcg::Use<MCFace  >::template AsFaceType >{};

  class MCVertex  : public Vertex< MCUsedTypes, vertex::Coord3f, vertex::Color4b, vertex::Mark, vertex::VFAdj, vertex::BitFlags, vertex::Qualityf>{};
  class MCEdge : public Edge<MCUsedTypes,edge::VertexRef>{};
  class MCFace    : public Face< MCUsedTypes, face::InfoOcf, face::VertexRef, face::FFAdjOcf, face::VFAdjOcf, face::BitFlags> {};
  class MCMesh	: public vcg::tri::TriMesh< std::vector< MCVertex>, face::vector_ocf< MCFace > > {};


  //******************************************
  //typedef Voxel<float> Voxelf;
  typedef Voxelfc Voxelf;
  //******************************************

  class Parameter
  {
  public:
    Parameter()
    {
      NCell=10000;
      WideNum= 3;
      WideSize=0;
      VoxSize=0;
      IPosS=Point3i(0,0,0);  // SubVolume Start
      IPosE=Point3i(0,0,0);  // SubVolume End
      IPosB=Point3i(0,0,0);  // SubVolume to restart from in lexicographic order (useful for crashes)
      IPos=Point3i(0,0,0);
      IDiv=Point3i(1,1,1);
      VerboseLevel=0;
      SliceNum=1;
      FillThr=12;
      ExpAngleDeg=30;
      SmoothNum=1;
      RefillNum=1;
      IntraSmoothFlag = false;
      QualitySmoothAbs = 0.0f; //  0 means un-setted value.
      QualitySmoothVox = 3.0f; // expressed in voxel
      OffsetFlag=false;
      OffsetThr=-3;
      GeodesicQualityFlag=true;
      PLYFileQualityFlag=false;
      FullyPreprocessedFlag=false;
      SaveVolumeFlag=false;
      SafeBorder=1;
      CleaningFlag=false;
      SimplificationFlag=false;
      VertSplatFlag=false;
      MergeColor=false;
      basename = "plymcout";
    }

    int NCell;
    int WideNum;
    float WideSize;
    float VoxSize;
    Point3i IPosS;  // SubVolume Start
    Point3i IPosE;  // SubVolume End
    Point3i IPosB;  // SubVolume to restart from in lexicographic order (useful for crashes)
    Point3i IPos;
    Point3i IDiv;
    int VerboseLevel;
    int SliceNum;
    int FillThr;
    float ExpAngleDeg;
    int SmoothNum;
    int RefillNum;
    bool IntraSmoothFlag;
    float QualitySmoothAbs; //  0 means un-setted value.
    float QualitySmoothVox; // expressed in voxel
    bool OffsetFlag;
    float OffsetThr;
    bool GeodesicQualityFlag;
    bool PLYFileQualityFlag;
    bool FullyPreprocessedFlag;
    bool CleaningFlag;
    bool SaveVolumeFlag;
    int SafeBorder;
    bool SimplificationFlag;
    bool VertSplatFlag;
    bool MergeColor;
    std::string basename;
    std::vector<std::string> OutNameVec;
    std::vector<std::string> OutNameSimpVec;
  }; //end Parameter class

  /// PLYMC Data
  MeshProvider MP;
  Parameter p;
  Volume<Voxelf> VV;
  char errorMessage[1024];

/// PLYMC Methods

  bool InitMesh(SMesh &m, const char *filename, Matrix44f Tr)
  {
    int loadmask;
    int ret = tri::io::Importer<SMesh>::Open(m,filename,loadmask);
    if(ret)
    {
      printf("Error: unabe to open mesh '%s'",filename);
      return false;
    }

    if(p.VertSplatFlag)
    {
      if(!(loadmask & tri::io::Mask::IOM_VERTNORMAL))
      {
        if(m.FN()==0)
        {
          sprintf(errorMessage,"%sError: mesh has not per vertex normals\n",errorMessage);
          return false;
        }
        else
        {          
          tri::Clean<SMesh>::RemoveUnreferencedVertex(m);
          tri::Allocator<SMesh>::CompactEveryVector(m);
          tri::UpdateNormal<SMesh>::PerVertexNormalizedPerFaceNormalized(m);          
        }
      }
      tri::UpdateNormal<SMesh>::NormalizePerVertex(m);    
      int badNormalCnt=0;
      for(SVertexIterator vi=m.vert.begin(); vi!=m.vert.end();++vi)
        if(math::Abs(SquaredNorm((*vi).N())-1.0)>0.0001)
        {
          badNormalCnt++;
          tri::Allocator<SMesh>::DeleteVertex(m,*vi);
        }
      tri::Allocator<SMesh>::CompactEveryVector(m);      
       if(badNormalCnt > m.VN()/10)
        {
          sprintf(errorMessage,"%sError: mesh has null normals\n",errorMessage);
          return false;
        }
      
      if(!(loadmask & tri::io::Mask::IOM_VERTQUALITY))
        tri::UpdateQuality<SMesh>::VertexConstant(m,0);
      tri::UpdateNormal<SMesh>::PerVertexMatrix(m,Tr);
    }
    else // processing for triangle meshes
    {
      if(!p.FullyPreprocessedFlag)
      {
        if(p.CleaningFlag){
          int dup = tri::Clean<SMesh>::RemoveDuplicateVertex(m);
          int unref =  tri::Clean<SMesh>::RemoveUnreferencedVertex(m);
          printf("Removed %i duplicates and %i unref",dup,unref);
        }

        tri::UpdateNormal<SMesh>::PerVertexNormalizedPerFaceNormalized(m);
        if(p.GeodesicQualityFlag) {
          tri::UpdateTopology<SMesh>::VertexFace(m);
          tri::UpdateFlags<SMesh>::FaceBorderFromVF(m);
          tri::Geodesic<SMesh>::DistanceFromBorder(m);
        }
      }

      tri::UpdatePosition<SMesh>::Matrix(m,Tr,true);
      tri::UpdateBounding<SMesh>::Box(m);
      printf("Init Mesh %s (%ivn,%ifn)\n",filename,m.vn,m.fn);
    }
    for(SVertexIterator vi=m.vert.begin(); vi!=m.vert.end();++vi)
      VV.Interize((*vi).P());
    return true;
  }

  // This function add a mesh (or a point cloud to the volume)
// the point cloud MUST have normalized vertex normals.
    bool AddMeshToVolumeM(SMesh &m, std::string meshname, const double w )
    {
      tri::RequireCompactness(m);
      if(!m.bbox.Collide(VV.SubBoxSafe)) return false;
      size_t found =meshname.find_last_of("/\\");
      std::string shortname = meshname.substr(found+1);

      Volume <Voxelf> B;
      B.Init(VV);

      bool res=false;
      double quality=0;

      // Now add the mesh to the volume
      if(!p.VertSplatFlag)
      {
        float minq=std::numeric_limits<float>::max(), maxq=-std::numeric_limits<float>::max();
            // Calcolo range qualita geodesica PER FACCIA come media di quelle per vertice
            for(SFaceIterator fi=m.face.begin(); fi!=m.face.end();++fi){
                (*fi).Q()=((*fi).V(0)->Q()+(*fi).V(1)->Q()+(*fi).V(2)->Q())/3.0f;
                minq=std::min((*fi).Q(),minq);
                maxq=std::max((*fi).Q(),maxq);
            }

            // La qualita' e' inizialmente espressa come distanza assoluta dal bordo della mesh
            printf("Q [%4.2f  %4.2f] \n",minq,maxq);
            bool closed=false;
            if(minq==maxq) closed=true;  // se la mesh e' chiusa la  ComputeGeodesicQuality mette la qualita a zero ovunque
            // Classical approach: scan each face
            int tt0=clock();
            printf("---- Face Rasterization");
            for(SFaceIterator fi=m.face.begin(); fi!=m.face.end();++fi)
                {
                    if(closed || (p.PLYFileQualityFlag==false && p.GeodesicQualityFlag==false)) quality=1.0;
                    else quality=w*(*fi).Q();
                    if(quality)
                            res |= B.ScanFace((*fi).V(0)->P(),(*fi).V(1)->P(),(*fi).V(2)->P(),quality,(*fi).N());
                }
            printf(" : %li\n",clock()-tt0);

    } else
    {	// Splat approach add only the vertices to the volume
        printf("Vertex Splatting\n");
        for(SVertexIterator vi=m.vert.begin();vi!=m.vert.end();++vi)
                {
                    if(p.PLYFileQualityFlag==false) quality=1.0;
                    else quality=w*(*vi).Q();
                    if(quality)
                        res |= B.SplatVert((*vi).P(),quality,(*vi).N(),(*vi).C());
                }
    }
    if(!res) return false;

    int vstp=0;
    if(p.VerboseLevel>0) {
      B.SlicedPPM(shortname.c_str(),std::string(SFormat("%02i",vstp)).c_str(),p.SliceNum	);
      B.SlicedPPMQ(shortname.c_str(),std::string(SFormat("%02i",vstp)).c_str(),p.SliceNum	);
      vstp++;
    }
    for(int i=0;i<p.WideNum;++i) {
    B.Expand(math::ToRad(p.ExpAngleDeg));
        if(p.VerboseLevel>1) B.SlicedPPM(shortname.c_str(),SFormat("%02ie",vstp++),p.SliceNum	);
        B.Refill(p.FillThr);
        if(p.VerboseLevel>1) B.SlicedPPM(shortname.c_str(),SFormat("%02if",vstp++),p.SliceNum	);
        if(p.IntraSmoothFlag)
        {
            Volume <Voxelf> SM;
            SM.Init(VV);
            SM.CopySmooth(B,1,p.QualitySmoothAbs);
            B=SM;
            if(p.VerboseLevel>1) B.SlicedPPM(shortname.c_str(),SFormat("%02is",vstp++),p.SliceNum	);
//			if(VerboseLevel>1) B.SlicedPPMQ(shortname,SFormat("%02is",vstp),SliceNum	);
        }
    }
    if(p.SmoothNum>0)
        {
            Volume <Voxelf> SM;
            SM.Init(VV);
            SM.CopySmooth(B,1,p.QualitySmoothAbs);
            B=SM;
            if(p.VerboseLevel>1) B.SlicedPPM(shortname.c_str(),SFormat("%02isf",vstp++),p.SliceNum	);
        }
    VV.Merge(B);
    if(p.VerboseLevel>0) VV.SlicedPPMQ(std::string("merge_").c_str(),shortname.c_str(),p.SliceNum	);
    return true;
}

bool Process(vcg::CallBackPos *cb=0)
{
  sprintf(errorMessage,"%s", "");
  printf("bbox scanning...\n"); fflush(stdout);
  Matrix44f Id; Id.SetIdentity();
  MP.InitBBox();
  printf("Completed BBox Scanning                   \n");
  Box3f fullb = MP.fullBB();
  assert (!fullb.IsNull());
  assert (!fullb.IsEmpty());
  // Calcolo gridsize
  Point3f voxdim;
  fullb.Offset(fullb.Diag() * 0.1 );

  int saveMask=0;
  saveMask|=tri::io::Mask::IOM_VERTQUALITY;
  if(p.MergeColor) saveMask |= tri::io::Mask::IOM_VERTCOLOR ;

  voxdim = fullb.max - fullb.min;

  // if kcell==0 the number of cells is computed starting from required voxel size;
  __int64 cells;
  if(p.NCell>0) cells = (__int64)(p.NCell)*(__int64)(1000);
  else cells = (__int64)(voxdim[0]/p.VoxSize) * (__int64)(voxdim[1]/p.VoxSize) *(__int64)(voxdim[2]/p.VoxSize) ;

  {
    Volume<Voxelf> B; // local to this small block

    Box3f fullbf; fullbf.Import(fullb);
    B.Init(cells,fullbf,p.IDiv,p.IPosS);
    B.Dump(stdout);
    if(p.WideSize>0) p.WideNum=p.WideSize/B.voxel.Norm();

    // Now the volume has been determined; the quality threshold in absolute units can be computed
    if(p.QualitySmoothAbs==0)
      p.QualitySmoothAbs= p.QualitySmoothVox * B.voxel.Norm();
  }


  int TotAdd=0,TotMC=0,TotSav=0; // partial timings counter

  for(p.IPos[0]=p.IPosS[0];p.IPos[0]<=p.IPosE[0];++p.IPos[0])
    for(p.IPos[1]=p.IPosS[1];p.IPos[1]<=p.IPosE[1];++p.IPos[1])
      for(p.IPos[2]=p.IPosS[2];p.IPos[2]<=p.IPosE[2];++p.IPos[2])
        if((p.IPos[2]+(p.IPos[1]*p.IDiv[2])+(p.IPos[0]*p.IDiv[2]*p.IDiv[1])) >=
           (p.IPosB[2]+(p.IPosB[1]*p.IDiv[2])+(p.IPosB[0]*p.IDiv[2]*p.IDiv[1]))) // skip until IPos >= IPosB
        {
          printf("----------- SubBlock %2i %2i %2i ----------\n",p.IPos[0],p.IPos[1],p.IPos[2]);
          //Volume<Voxelf> B;
          int t0=clock();

          Box3f fullbf; fullbf.Import(fullb);

          VV.Init(cells,fullbf,p.IDiv,p.IPos);
          printf("\n\n --------------- Allocated subcells. %i\n",VV.Allocated());

          std::string filename=p.basename;
          if(p.IDiv!=Point3i(1,1,1))
          {
            std::string subvoltag;
            VV.GetSubVolumeTag(subvoltag);
            filename+=subvoltag;
          }
          /********** Grande loop di scansione di tutte le mesh *********/
          bool res=false;
          if(!cb) printf("Step 1: Converting meshes into volume\n");
          for(int i=0;i<MP.size();++i)
          {
            Box3f bbb= MP.bb(i);
            /**********************/
            if(cb) cb((i+1)/MP.size(),"Step 1: Converting meshes into volume");
            /**********************/
            // if bbox of mesh #i is part of the subblock, then process it
            if(bbb.Collide(VV.SubBoxSafe))
            {
              SMesh *sm;
              if(!MP.Find(i,sm) )
              {
                res = InitMesh(*sm,MP.MeshName(i).c_str(),MP.Tr(i));
                if(!res)
                {
                  sprintf(errorMessage,"%sFailed Init of mesh %s\n",errorMessage,MP.MeshName(i).c_str());
                  return false ;
                }
              }
              res |= AddMeshToVolumeM(*sm, MP.MeshName(i),MP.W(i));
            }
          }

          //B.Normalize(1);
          printf("End Scanning\n");
          if(p.OffsetFlag)
          {
            VV.Offset(p.OffsetThr);
            if (p.VerboseLevel>0)
            {
              VV.SlicedPPM("finaloff","__",p.SliceNum);
              VV.SlicedPPMQ("finaloff","__",p.SliceNum);
            }
          }
          //if(p.VerboseLevel>1) VV.SlicedPPM(filename.c_str(),SFormat("_%02im",i),p.SliceNum	);

          for(int i=0;i<p.RefillNum;++i)
          {
            VV.Refill(3,6);
            if(p.VerboseLevel>1) VV.SlicedPPM(filename.c_str(),SFormat("_%02imsr",i),p.SliceNum	);
            //if(VerboseLevel>1) VV.SlicedPPMQ(filename,SFormat("_%02ips",i++),SliceNum	);
          }

          for(int i=0;i<p.SmoothNum;++i)
          {
            Volume <Voxelf> SM;
            SM.Init(VV);
            printf("%2i/%2i: ",i,p.SmoothNum);
            SM.CopySmooth(VV,1,p.QualitySmoothAbs);
            VV=SM;
            VV.Refill(3,6);
            if(p.VerboseLevel>1) VV.SlicedPPM(filename.c_str(),SFormat("_%02ims",i),p.SliceNum	);
          }

          int t1=clock();  //--------
          TotAdd+=t1-t0;
          printf("Extracting surface...\r");
          if (p.VerboseLevel>0)
          {
            VV.SlicedPPM("final","__",p.SliceNum);
            VV.SlicedPPMQ("final","__",p.SliceNum);
          }
          MCMesh me;
          if(res)
          {
            typedef vcg::tri::TrivialWalker<MCMesh, Volume <Voxelf> >	  Walker;
            typedef vcg::tri::MarchingCubes<MCMesh, Walker>             MarchingCubes;

            Walker walker;
            MarchingCubes	mc(me, walker);
            /**********************/
            if(cb) cb(50,"Step 2: Marching Cube...");
            else printf("Step 2: Marching Cube...\n");
            /**********************/
            walker.SetExtractionBox(VV.SubPartSafe);
            walker.BuildMesh(me,VV,mc,0);

            typename MCMesh::VertexIterator vi;
            Box3f bbb; bbb.Import(VV.SubPart);
            for(vi=me.vert.begin();vi!=me.vert.end();++vi)
            {
              if(!bbb.IsIn((*vi).P()))
                vcg::tri::Allocator< MCMesh >::DeleteVertex(me,*vi);
              VV.DeInterize((*vi).P());
            }
            for (typename MCMesh::FaceIterator fi = me.face.begin(); fi != me.face.end(); ++fi)
            {
              if((*fi).V(0)->IsD() || (*fi).V(1)->IsD() || (*fi).V(2)->IsD() )
                vcg::tri::Allocator< MCMesh >::DeleteFace(me,*fi);
              else std::swap((*fi).V1(0), (*fi).V2(0));
            }

            int t2=clock();  //--------
            TotMC+=t2-t1;
            if(me.vn >0 || me.fn >0)
            {
              p.OutNameVec.push_back(filename+std::string(".ply"));
              tri::io::ExporterPLY<MCMesh>::Save(me,p.OutNameVec.back().c_str(),saveMask);
              if(p.SimplificationFlag)
              {
                /**********************/
                if(cb) cb(50,"Step 3: Simplify mesh...");
                else printf("Step 3: Simplify mesh...\n");
                /**********************/
                p.OutNameSimpVec.push_back(filename+std::string(".d.ply"));
                me.face.EnableVFAdjacency();
                MCSimplify<MCMesh>(me, VV.voxel[0]/4.0);
                tri::Allocator<MCMesh>::CompactFaceVector(me);
                me.face.EnableFFAdjacency();
                tri::Clean<MCMesh>::RemoveTVertexByFlip(me,20,true);
                tri::Clean<MCMesh>::RemoveFaceFoldByFlip(me);
                tri::io::ExporterPLY<MCMesh>::Save(me,p.OutNameSimpVec.back().c_str(),saveMask);
              }
            }
            int t3=clock();  //--------
            TotSav+=t3-t2;

          }

          printf("Mesh Saved '%s':  %8d vertices, %8d faces                   \n",(filename+std::string(".ply")).c_str(),me.vn,me.fn);
          printf("Adding Meshes %8i\n",TotAdd);
          printf("MC            %8i\n",TotMC);
          printf("Saving        %8i\n",TotSav);
          printf("Total         %8i\n",TotAdd+TotMC+TotSav);
        }
        else
        {
          printf("----------- skipping SubBlock %2i %2i %2i ----------\n",p.IPos[0],p.IPos[1],p.IPos[2]);
        }
  return true;
}


}; //end PlyMC class


template < class MeshType, class VertexPair>
         class PlyMCTriEdgeCollapse: public MCTriEdgeCollapse< MeshType, VertexPair, PlyMCTriEdgeCollapse<MeshType,VertexPair> > {
                        public:
                        typedef  MCTriEdgeCollapse< MeshType, VertexPair, PlyMCTriEdgeCollapse  > MCTEC;
            inline PlyMCTriEdgeCollapse(  const VertexPair &p, int i, BaseParameterClass *pp) :MCTEC(p,i,pp){}
 };



template<   class MeshType>
int MCSimplify( MeshType &m, float absoluteError, bool preserveBB, vcg::CallBackPos *cb)
{

	typedef PlyMCTriEdgeCollapse<MeshType,BasicVertexPair<typename MeshType::VertexType> > MyColl;
	typedef typename MeshType::FaceIterator FaceIterator;
	typedef typename MeshType::CoordType CoordType;


	tri::UpdateBounding<MeshType>::Box(m);
	tri::UpdateTopology<MeshType>::VertexFace(m);
	TriEdgeCollapseMCParameter pp;
	pp.bb.Import(m.bbox);
	pp.preserveBBox=preserveBB;
	vcg::LocalOptimization<MeshType> DeciSession(m,&pp);
	if(absoluteError==0)
	{
      // guess the mc side.
      // In a MC mesh the vertices are on the egdes of the cells. and the edges are (mostly) on face of the cells.
      // If you have  2 vert over the same face xy they share z

		std::vector<float> ZSet;
		for(FaceIterator fi = m.face.begin();fi!=m.face.end();++fi)
		if(!(*fi).IsD())
		{
			CoordType v0=(*fi).V(0)->P();
			CoordType v1=(*fi).V(1)->P();
			CoordType v2=(*fi).V(2)->P();
			if(v0[2]==v1[2] && v0[1]!=v1[1] && v0[0]!=v1[0]) ZSet.push_back(v0[2]);
			if(v0[2]==v2[2] && v0[1]!=v2[1] && v0[0]!=v2[0]) ZSet.push_back(v0[2]);
			if(v1[2]==v2[2] && v1[1]!=v2[1] && v1[0]!=v2[0]) ZSet.push_back(v1[2]);
			if(ZSet.size()>100) break;
		}
		if (ZSet.size() == 0) return -1;	//no straight edges found. exit with error
		std::sort(ZSet.begin(),ZSet.end());
		std::vector<float>::iterator lastV = std::unique(ZSet.begin(),ZSet.end());
		ZSet.resize(lastV-ZSet.begin());
		float Delta=0;
		for(size_t i = 0; i< ZSet.size()-1;++i)
		{
			Delta = std::max(ZSet[i+1]-ZSet[i],Delta);
			//qDebug("%f",Delta);
		}
		absoluteError= Delta/4.0f;
	}
	//qDebug("Simplifying at absoluteError=%f",absoluteError);

	float TargetError = absoluteError;
	char buf[1024];
	DeciSession.template Init< MyColl > ();

	pp.areaThr=TargetError*TargetError;
	DeciSession.SetTimeBudget(1.0f);
	if(TargetError < std::numeric_limits<float>::max() ) DeciSession.SetTargetMetric(TargetError);
	while(DeciSession.DoOptimization() && DeciSession.currMetric < TargetError)
	{
		sprintf(buf,"Simplyfing %7i err %9g \r",m.fn,DeciSession.currMetric);
		if (cb) cb(int(100.0f*DeciSession.currMetric/TargetError),buf);
	}

	return 1; //success
}


} // end namespace tri
} // end namespace vcg

#endif
