#ifndef __MIQ__
#define __MIQ__


#define SIZEQUADS 512
#define V_SIZE 1
#define SIZEPARA 1024

#include <iostream>
#include <vector>
#include "mesh_type.h"
#include "quadrangulator.h"
#include "poisson_solver.h"
#include "param_stats.h"
#include "seams_initializer.h"
#include "vertex_indexing.h"

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

#define USECOMISO

template <class ScalarType>
ScalarType Gauss(ScalarType &value)
{
    const ScalarType E_NEPER=2.71828;
    const ScalarType den=sqrt(2.0*M_PI);
    ScalarType exponent=-0.5*pow(value,2);
    ScalarType res=((1.0/den)*pow(E_NEPER,exponent));
    return res;
}

template <class MeshType,class PolyMesh>
class MIQ{

public:
  typename MeshType::template PerFaceAttributeHandle<float> Handle_Stiffness;

  // Different stiffening mode
  enum StiffMode{NO_STIFF = 0,GAUSSIAN = 1,ITERATIVE = 2};

  // Init
  MIQ(MeshType &_mesh):mesh(_mesh),PSolver(mesh){};

  // Load a mesh from file
  bool LoadMesh(const std::string PathMesh)
  {
      int position=PathMesh.find(".ply");
      int err;
      if (position==-1)
      {
          position=PathMesh.find(".obj");
          //vcg::tri::io::ImporterOBJ<CMesh>::Info objInf;
          int mask;
          vcg::tri::io::ImporterOBJ<CMesh>::LoadMask(PathMesh.c_str(),mask);
          err=vcg::tri::io::ImporterOBJ<CMesh>::Open(mesh,PathMesh.c_str(),mask);
          assert(position!=-1);
          if (err!=0)return false;
      }
      else
      {
          err=vcg::tri::io::ImporterPLY<CMesh>::Open(mesh,PathMesh.c_str());
          if (err!=ply::E_NOERROR)return false;
      }
      ///UPDATE MESH STUFF
      vcg::tri::UpdateBounding<CMesh>::Box(mesh);
      vcg::tri::UpdateNormal<CMesh>::PerVertexNormalizedPerFace(mesh);
      vcg::tri::UpdateNormal<CMesh>::PerFaceNormalized(mesh);
      vcg::tri::UpdateTopology<CMesh>::FaceFace(mesh);
      vcg::tri::UpdateTopology<CMesh>::VertexFace(mesh);
      vcg::tri::UpdateFlags<CMesh>::FaceBorderFromFF(mesh);
      vcg::tri::UpdateFlags<CMesh>::VertexBorderFromFace(mesh);
      return true;
  }

  // Load a field from file
  bool LoadField(const std::string PathFField,const std::string PathMesh=NULL)
  {
      SInit.Init(&mesh);

      ///FIELD LOADING
      int position=PathFField.find(".ffield");
      if (position==-1)
      {
          position=PathFField.find(".grad");
          assert(position!=-1);
          SInit.InitFromGradOBJ(PathFField,PathMesh);
      }
      else
          SInit.InitFromFField(PathFField);

      VInd.Init(&mesh);
      //mesh.ScaleToUnitBox();
      VInd.InitMapping();
      VInd.InitFaceIntegerVal();
      VInd.InitSeamInfo();
      return true;
  }

  // Load a mesh and a field from file
  bool LoadMeshField(const std::string PathMesh, const std::string PathFField)
  {
      bool loaded=LoadMesh(PathMesh);
      if (!loaded)return false;
      LoadField(PathFField,PathMesh);

      return true;
  }

  // Parametrize the mesh
  void Solve(StiffMode stiffMode, double Stiffness = 5.0,
             double GradientSize = 30.0, bool DirectRound = false,
             int iter = 5, int localIter = 5, bool DoRound = true)
  {
      if (mesh.fn==0)return;
      InitDefaultStiffening();
      if (stiffMode==GAUSSIAN)
      {
        AddGaussStiffening(Stiffness);
        PSolver.SolvePoisson(GradientSize,1.f,DirectRound,localIter,DoRound);
      }
      else
      if (stiffMode==ITERATIVE)
      {
          for (int i=0;i<iter;i++)
          {
            PSolver.SolvePoisson(GradientSize,1.f,DirectRound,localIter,DoRound);
            int nflips=NumFlips(mesh);
            ScaleGLtoInt();
            bool folded=updateStiffening(GradientSize);
            ScaleInttoGL();
            printf("ITERATION %d FLIPS %d \n",i,nflips);
            if (!folded)break;
          }
      }
      else
      if (stiffMode==NO_STIFF)
      {
          PSolver.SolvePoisson(GradientSize,1.f,DirectRound,localIter,DoRound);
      }
      int nflips=NumFlips(mesh);
      printf("**** END OPTIMIZING #FLIPS %d  ****\n",nflips);
      fflush(stdout);
      SelectFlippedFaces(mesh);
      UpdateUVBox();
  }

  // Generate a quad mesh starting from the parametrization
  void Quadrangulate(PolyMesh &poly,double factor = 1)
  {
      Quad.Quadrangulate(mesh,poly,factor);
  }

//  void removeDuplicateVertices(std::vector<std::vector<double> >& V, std::vector< std::vector<int > >& F);

//  void exportQuadMesh(std::string out);

  //bool LoadData(std::string PathMesh,
   //             std::string PathFField);
  // Bounding box of the param domain
  vcg::Box2<double> UVBox;

  void ScaleGLtoInt()
  {
      double factor=(double)SIZEPARA/(double)SIZEQUADS;
      for (unsigned int i=0;i<mesh.face.size();i++)
      {
          if (mesh.face[i].IsD())continue;
          mesh.face[i].WT(0).P()*=factor;
          mesh.face[i].WT(1).P()*=factor;
          mesh.face[i].WT(2).P()*=factor;
      }
  }

  void ScaleInttoGL()
  {
      double factor=(double)SIZEQUADS/(double)SIZEPARA;
      for (unsigned int i=0;i<mesh.face.size();i++)
      {
          if (mesh.face[i].IsD())continue;
          mesh.face[i].WT(0).P()*=factor;
          mesh.face[i].WT(1).P()*=factor;
          mesh.face[i].WT(2).P()*=factor;
      }
  }

  void UpdateUVBox()
  {
      UVBox.SetNull();
      for (unsigned int i=0;i<mesh.face.size();i++)
      {
          UVBox.Add((mesh.face[i].WT(0).P()));
          UVBox.Add((mesh.face[i].WT(1).P()));
          UVBox.Add((mesh.face[i].WT(2).P()));
      }
  }

  // Quadrangulator
  Quadrangulator<CMesh,PolyMesh> Quad;

  void colorByStiffening(typename MeshType::ScalarType MaxVal=16)
  {
      bool hasStiffness = vcg::tri::HasPerFaceAttribute(mesh,std::string("Stiffness"));
      assert(hasStiffness);
      for (unsigned int i=0;i<mesh.face.size();i++)
      {
          //CMesh::ScalarType val=MaxVal-mesh.face[i].stiffening+1;
          CMesh::ScalarType val=MaxVal-Handle_Stiffness[i]+1;
          if (val<1)val=1;
          mesh.face[i].C()=vcg::Color4b::ColorRamp(1.0,MaxVal,val);
      }
  }

private:

  void AddStiffening(typename MeshType::ScalarType C,int radius=4)
  {
      bool hasStiffness = vcg::tri::HasPerFaceAttribute(mesh,std::string("Stiffness"));
      if(!hasStiffness)
          Handle_Stiffness=vcg::tri::Allocator<CMesh>::AddPerFaceAttribute<float>(mesh,std::string("Stiffness"));

      bool hasSingular = vcg::tri::HasPerVertexAttribute(mesh,std::string("Singular"));
      assert(hasSingular);
      CMesh::PerVertexAttributeHandle<bool> Handle_Singular;
      Handle_Singular=vcg::tri::Allocator<CMesh>::GetPerVertexAttribute<bool>(mesh,std::string("Singular"));

      std::vector<CMesh::VertexType*> to_stiff;
      for(unsigned int i=0;i<mesh.vert.size();i++)
      {
          CMesh::VertexType *v=&mesh.vert[i];
          if (v->IsD())continue;
          //if (!v->IsSingular())continue;
          if (!Handle_Singular[v])continue;
          to_stiff.push_back(v);
      }
      for(unsigned int i=0;i<mesh.face.size();i++)
      {

          CMesh::FaceType *f=&(mesh.face[i]);
          if (f->IsD())continue;
          if (!f->IsV())continue;
          to_stiff.push_back(f->V(0));
          to_stiff.push_back(f->V(1));
          to_stiff.push_back(f->V(2));
      }
      std::sort(to_stiff.begin(),to_stiff.end());
      std::vector<CMesh::VertexType*>::iterator new_end=std::unique(to_stiff.begin(),to_stiff.end());
      int dist=distance(to_stiff.begin(),new_end);
      to_stiff.resize(dist);
      for (unsigned int i=0;i<to_stiff.size();i++)
      {
          CMesh::VertexType *v=to_stiff[i];
          for (int r=0;r<radius;r++)
          {
              CMesh::ScalarType stiffVal=((CMesh::ScalarType)r)/(CMesh::ScalarType)radius;//((ScalarType)(radius-r))/(ScalarType)radius;
              stiffVal*=3.0;
              stiffVal=Gauss(stiffVal)/0.4;
              stiffVal=1+(stiffVal*C);
              std::vector<CMesh::FaceType*> ring;
              //mesh.GetConnectedFaces(v,r,ring);
              VFExtendedStarVF(v,r,ring);
              ///then set stiffening
              for (unsigned int k=0;k<ring.size();k++)
              {
                  CMesh::FaceType* f=ring[k];
                  //if (f->stiffening<stiffVal)
                  //    f->stiffening=stiffVal;
                  if (Handle_Stiffness[f]<stiffVal)
                      Handle_Stiffness[f]=stiffVal;
              }
          }
      }
  }


  bool updateStiffening(typename MeshType::ScalarType grad_size)
  {
      bool hasStiffness = vcg::tri::HasPerFaceAttribute(mesh,std::string("Stiffness"));
      if(!hasStiffness)
          Handle_Stiffness=vcg::tri::Allocator<CMesh>::AddPerFaceAttribute<float>(mesh,std::string("Stiffness"));

      bool flipped = NumFlips(mesh)>0;
      //if (h == 0.0)
      //    return flipped;
      //
      //assert(h != 0.0);

      if (!flipped)
          return false;
      CMesh::ScalarType maxL=0;
      CMesh::ScalarType maxD=0;
      if (flipped)
      {
          const double c = 1.0;
          const double d = 5.0;

          for (unsigned int i = 0; i < mesh.face.size(); ++i)
          {
              CMesh::ScalarType dist=Distortion(mesh.face[i],grad_size);
              if (dist>maxD)maxD=dist;
              CMesh::ScalarType absLap=fabs(LaplaceDistortion(mesh.face[i], grad_size));
              if (absLap>maxL)maxL=absLap;
              CMesh::ScalarType stiffDelta = std::min(c * absLap, d);
              //mesh.face[i].stiffening+=stiffDelta;
              Handle_Stiffness[i]+=stiffDelta;
          }
      }
      printf("Maximum Distorsion %4.4f \n",maxD);
      printf("Maximum Laplacian %4.4f \n",maxL);
      return flipped;
  }

  void InitDefaultStiffening()
  {
      bool hasStiffness = vcg::tri::HasPerFaceAttribute(mesh,std::string("Stiffness"));
      if(!hasStiffness)
          Handle_Stiffness=vcg::tri::Allocator<CMesh>::AddPerFaceAttribute<float>(mesh,std::string("Stiffness"));

      for(unsigned int i=0;i<mesh.face.size();i++)
      {
          CMesh::FaceType *f=&(mesh.face[i]);
          //f->stiffening=1;
          Handle_Stiffness[f]=1;
      }
  }

  void AddGaussStiffening(typename MeshType::ScalarType C)
  {
      int radius=floor(C);
      if (C<4)radius=4;
      AddStiffening(C,radius);
  }

  // Mesh class
  MeshType &mesh;

  // Quad mesh class
  //CMesh quadmesh;

  ///seams initializer
  SeamsInitializer<MeshType> SInit;

  // Vertex indexing class used for the solver
  VertexIndexing<MeshType> VInd;

  // Poisson solver
  PoissonSolver<MeshType> PSolver;

};

#endif
