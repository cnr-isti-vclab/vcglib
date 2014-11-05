/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004                                                \/)\/    *
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
#ifndef __VCG_IMPLICIT_SMOOTHER
#define __VCG_IMPLICIT_SMOOTHER

#include <eigenlib/Eigen/Sparse>
#include <vcg/complex/algorithms/mesh_to_matrix.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <vcg/complex/algorithms/smooth.h>

#define PENALTY 10000

namespace vcg{


template <class MeshType>
class ImplicitSmoother
{
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> MatrixXm;

    MeshType &to_smooth_mesh;

public:

    struct FaceConstraint
    {
        int numF;
        std::vector<std::vector<ScalarType> > BarycentricW;
        CoordType TargetPos;
        ScalarType Weight;

        FaceConstraint()
        {
            numF=-1;
            Weight=0;
        }

        FaceConstraint(int _numF,const std::vector<std::vector<ScalarType> > &_BarycentricW,
                       const CoordType &_TargetPos,ScalarType _Weight)
        {
            numF=_numF;
            BarycentricW= std::vector<std::vector<ScalarType> > (_BarycentricW.begin(),_BarycentricW.end());
            TargetPos=_TargetPos;
            Weight=_Weight;
        }
    };

    struct SmoothParam
    {
        //the amount of smoothness, useful only if we set the mass matrix
        ScalarType lambda;
        //the use of mass matrix to keep the mesh close to its original position
        //(weighted per area distributed on vertices)
        bool useMassMatrix;
        //this bool is used to fix the border vertices of the mesh or not
        bool fixBorder;
        //this bool is used to set if cotangent weight is used, this flag to false means uniform laplacian
        bool useCotWeight;
        //the set of fixed vertices
        std::vector<int> FixedV;
        //the set of faces for barycentric constraints
        std::vector<FaceConstraint> ConstrainedF;

        SmoothParam()
        {
            lambda=0.2;
            useMassMatrix=true;
            fixBorder=false;
            useCotWeight=false;
        }
    };

private:

    //    void GetMassMatrix(MeshType &mesh,
    //                       std::vector<std::pair<int,int> > &index,
    //                       std::vector<ScalarType> &entry)
    //    {
    //        entry.clear();
    //        index.clear();

    //        //calculate area
    //        vcg::tri::UpdateQuality<MeshType>::FaceArea(mesh);
    //        //then distribute per vertex
    //        vcg::tri::UpdateQuality<MeshType>::VertexFromFace(mesh);

    //        //3, one per coordinate
    //        index.resize(mesh.vert.size()*3,-1);
    //        entry.resize(mesh.vert.size()*3,0);

    //        //store the index and the scalar for the sparse matrix
    //        for (size_t i=0;i<mesh.vert.size();i++)
    //        {
    //            for (size_t j=0;j<3;j++)
    //            {
    //                int currI=(i*3)+j;
    //                index[currI]=currI;
    //                entry[currI]=mesh.vert[i].Q();
    //            }
    //        }
    //    }

    void InitSparse(const std::vector<std::pair<int,int> > &Index,
                    const std::vector<ScalarType> &Values,
                    const size_t m,
                    const size_t n,
                    Eigen::SparseMatrix<ScalarType>& X)
    {

        assert(Index.size()==Values.size());

        std::vector<Eigen::Triplet<ScalarType> > IJV;
        IJV.reserve(Index.size());

        for(size_t i= 0;i<Index.size();i++)
        {
            int row=Index[i].first;
            int col=Index[i].second;
            ScalarType val=Values[i];

            assert(row<m);
            assert(col<n);

            IJV.push_back(Eigen::Triplet<ScalarType>(row,col,val));
        }
        X.resize(m,n);
        X.setFromTriplets(IJV.begin(),IJV.end());
    }

    int SystemSize(SmoothParam & SParam)
    {
        int basic_size=to_smooth_mesh.vert.size();
        int constr_size=SParam.ConstrainedF.size();
        return (basic_size+constr_size);
    }

    void CollectHardConstraints(const SmoothParam &SParam,
                                std::vector<std::pair<int,int> > &IndexC,
                                std::vector<ScalarType> &WeightC)
    {
        std::vector<int> To_Fix;

        //collect fixed vert
        if (SParam.fixBorder)
        {
            //add penalization constra
            for (int i=0;i<to_smooth_mesh.vert.size();i++)
            {
                if (!to_smooth_mesh.vert[i].IsB())continue;
                To_Fix.push_back(i);
            }
        }
        //add additional fixed vertices constraint
        To_Fix.insert(To_Fix.end(),SParam.FixedV.begin(),SParam.FixedV.end());

        //sort and make them unique
        std::sort(To_Fix.begin(),To_Fix.end());
        typename std::vector<int>::iterator it= std::unique (To_Fix.begin(), To_Fix.end());
        To_Fix.resize( std::distance(To_Fix.begin(),it) );

        for (size_t i=0;i<To_Fix.size();i++)
        {
            for (int j=0;j<3;j++)
            {
                int IndexV=(i*3)+j;
                IndexC.push_back(std::pair<int,int>(IndexV,IndexV));
                WeightC.push_back((ScalarType)PENALTY);
            }
        }
    }

public:



    //static void Smooth
    //    static void SmoothLIbIGL(MeshType &mesh,
    //                             ScalarType lambda=0.5)//0.0001)
    //    {
    //        Eigen::MatrixXi F;
    //        Eigen::MatrixXd V,U;
    //        vcg::tri::MeshToMatrix<MeshType>::GetTriMeshData(mesh,F,V);
    //        U=V;

    //        Eigen::SparseMatrix<ScalarType> L,M;
    //        igl::cotmatrix(V,F,L);

    //        //compute the mass matrix
    //        static bool computed=false;
    //        igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);
    //        M.setIdentity();
    //        computed=true;

    //        //const auto & S = (M - 0.001*L);
    //        Eigen::SparseMatrix<double> S = (M - lambda*L);


    //        Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
    //        assert(solver.info() == Eigen::Success);

    //        U = solver.solve(M*U).eval();
    //        // Normalize to unit surface area (important for numerics)
    //        U.array() /= sqrt(M.diagonal().array().sum());

    //        //then assing per vertex positions
    //        printf("smoothed \n");
    //        fflush(stdout);
    //        for (size_t i=0;i<mesh.vert.size();i++)
    //        {
    //            mesh.vert[i].P().X()=U(i,0);
    //            mesh.vert[i].P().Y()=U(i,1);
    //            mesh.vert[i].P().Z()=U(i,2);
    //        }
    //    }

    void Smooth(SmoothParam &SParam)
    {

        //the laplacian and the mass matrix
        Eigen::SparseMatrix<ScalarType> L,M;

        //initialize the mass matrix
        std::vector<std::pair<int,int> > IndexM;
        std::vector<ScalarType> ValuesM;

        //add the entries for mass matrix
        if (SParam.useMassMatrix) MeshToMatrix<MeshType>::MassMatrixEntry(to_smooth_mesh,IndexM,ValuesM);

        //then also collect hard constraints
        //CollectHardConstraints(SParam,IndexM,ValuesM);

        //initialize sparse mass matrix
        InitSparse(IndexM,ValuesM,to_smooth_mesh.vert.size()*3,to_smooth_mesh.vert.size()*3,M);

        //get the entries for laplacian matrix
        std::vector<std::pair<int,int> > IndexL;
        std::vector<ScalarType> ValuesL;
        MeshToMatrix<MeshType>::GetLaplacianMatrix(to_smooth_mesh,IndexL,ValuesL,false);//SParam.useCotWeight);

        //initialize sparse laplacian matrix
        InitSparse(IndexL,ValuesL,to_smooth_mesh.vert.size()*3,to_smooth_mesh.vert.size()*3,L);

        //then solve the system
        Eigen::SparseMatrix<ScalarType> S = (M + SParam.lambda*L);

        //SimplicialLDLT
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<ScalarType > > solver(S);
        printf("output %d \n",solver.info());
        fflush(stdout);
        assert(solver.info() == Eigen::Success);

        int size=to_smooth_mesh.vert.size()*3;
        MatrixXm V(size,1);

        for (size_t i=0;i<to_smooth_mesh.vert.size();i++)
        {
            int index=i*3;
            assert(index<size);
            V(index,0)=to_smooth_mesh.vert[i].P().X();
            V(index+1,0)=to_smooth_mesh.vert[i].P().Y();
            V(index+2,0)=to_smooth_mesh.vert[i].P().Z();
        }

        V = solver.solve(M*V).eval();

        for (size_t i=0;i<to_smooth_mesh.vert.size();i++)
        {
            int index=i*3;
            to_smooth_mesh.vert[i].P().X()=V(index,0);
            to_smooth_mesh.vert[i].P().Y()=V(index+1,0);
            to_smooth_mesh.vert[i].P().Z()=V(index+2,0);
        }
    }

    ImplicitSmoother(MeshType &_to_smooth_mesh):to_smooth_mesh(_to_smooth_mesh){}
};

}//end namespace vcg

#endif
