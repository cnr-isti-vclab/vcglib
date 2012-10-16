#ifndef __MIQ__
#define __MIQ__

#include <iostream>
#include <vector>
#include "mesh_type.h"
#include "quadrangulator.h"
#include "poisson_solver.h"
#include "param_stats.h"
#include "seams_initializer.h"
#include "vertex_indexing.h"
#include "stiffening.h"
#include <vcg/complex/algorithms/clean.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

#define USECOMISO

template <class MeshType>
class MIQ_parametrization{

public:
  typename MeshType::template PerFaceAttributeHandle<float> Handle_Stiffness;

  // Different stiffening mode
  enum StiffMode{NO_STIFF = 0,GAUSSIAN = 1,ITERATIVE = 2};

  // Parametrize the mesh
  static void DoParameterize(MeshType &mesh,StiffMode stiffMode,
                             double Stiffness = 5.0,double GradientSize = 30.0,
                             bool DirectRound = false,int iter = 5,
                             int localIter = 5, bool DoRound = true)
  {
      PoissonSolver<MeshType> PSolver(mesh);
      if (mesh.fn==0)return;
      StiffeningInitializer<MeshType>::InitDefaultStiffening(mesh);
      if (stiffMode==GAUSSIAN)
      {
        StiffeningInitializer<MeshType>::AddGaussStiffening(mesh,Stiffness);
        PSolver.SolvePoisson(GradientSize,1.f,DirectRound,localIter,DoRound);
      }
      else
      if (stiffMode==ITERATIVE)
      {
          for (int i=0;i<iter;i++)
          {
            PSolver.SolvePoisson(GradientSize,1.f,DirectRound,localIter,DoRound);
            int nflips=NumFlips(mesh);
            bool folded=StiffeningInitializer<MeshType>::updateStiffeningJacobianDistorsion(mesh,GradientSize);
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

  }

public:

  static bool IsValid(MeshType &mesh)
  {
      int n_comp=vcg::tri::Clean<MeshType>::CountConnectedComponents(mesh);
      int non_manifE=vcg::tri::Clean<MeshType>::CountNonManifoldEdgeFF(mesh);
      int non_manifV=vcg::tri::Clean<MeshType>::CountNonManifoldVertexFF(mesh);
      return ((n_comp==1)&&(non_manifE==0)&&(non_manifV==0));
  }
  static void InitSeamsSing(MeshType &mesh,
                            bool orient_globally,
                            bool initMM,
                             bool initCuts)
  {
      SeamsInitializer<MeshType> SInit;
      SInit.Init(&mesh,orient_globally,initMM,initCuts);
  }

  static void Parametrize(MeshType &mesh,StiffMode stiffMode,
                          double Stiffness = 5.0,
                          double GradientSize = 30.0,
                          bool DirectRound = false,int iter = 5,
                          int localIter = 5, bool DoRound = true)
  {
      VertexIndexing<MeshType> VInd;

      VInd.Init(&mesh);
      VInd.InitMapping();
      VInd.InitFaceIntegerVal();
      VInd.InitSeamInfo();

      DoParameterize(mesh,stiffMode,Stiffness,GradientSize,DirectRound,iter,localIter , DoRound);
  }


};

#endif
