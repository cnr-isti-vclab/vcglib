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

#ifndef VCG_POISSON_SOLVER
#define VCG_POISSON_SOLVER

#include <Eigen/Sparse>

#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/parametrization/distortion.h>
#include <vcg/complex/algorithms/parametrization/uv_utils.h>
#include <wrap/callback.h>

namespace vcg {
namespace tri{

/** Parametrizes a disk-like mesh by solving one sparse linear system, in one of two modes.

  - Conformal (the default): least squares conformal maps, the cotangent Dirichlet energy minus
    the signed UV area. It needs at least two fixed vertices, which fix the similarity the
    energy is blind to.
  - Field-aligned (Init(true)): the gradients of U and V are fitted to the faces' normalized
    PD1 and PD2, so a unit of UV is a unit on the surface. It needs one fixed vertex, and up
    to date face normals.

  Vertices are fixed with FixVertex, SetBorderAsFixed, SetSelectedAsFixed or
  FixDefaultVertices; AddPointConstraint adds a soft constraint on one coordinate at a point
  inside a face. The result goes to the per-vertex and per-wedge texture coordinates.

  Typical use, on a mesh with face-face adjacency:
  \code
  tri::PoissonSolver<MyMesh> PS(m);
  if(PS.IsFeasible()) { PS.Init(); PS.FixDefaultVertices(); PS.SolvePoisson(); }
  \endcode
*/
template <class MeshType>
class PoissonSolver
{
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::CoordType CoordType;
    typedef Eigen::Triplet<double> Triplet;

    MeshType &mesh;

    /// variable of each vertex (by position in mesh.vert, -1 for deleted ones), and back
    std::vector<int> VertexToInd;
    std::vector<VertexType*> IndToVertex;

    std::vector<VertexType *> to_fix;
    struct PointConstraint { FaceType *f; CoordType bary; int coord; ScalarType value, weight; };
    std::vector<PointConstraint> pointConstraints;

    /// the system A x = b, in which x interleaves the U and V of each vertex
    std::vector<Triplet> A;
    Eigen::VectorXd b,x;

    bool use_direction_field,correct_fixed;

    int VertexIndex(const VertexType* v) const { return VertexToInd[tri::Index(mesh,v)]; }

    void SetValA(int Xindex,int Yindex,ScalarType val) { A.emplace_back(Xindex,Yindex,val); }

    ///add the area term, the antisymmetric coupling of U and V that turns the
    ///Dirichlet energy into the conformal one
    void AddAreaTerm(const int Vindexes[3])
    {
        const ScalarType entry=0.5;
        const ScalarType val[3][3]= { {0,  entry, -entry},
                                      {-entry,  0,  entry},
                                      {entry, -entry,  0} };
        for (int i=0;i<3;i++)
            for (int j=0;j<3;j++)
            {
                SetValA(Vindexes[i]*2+1,Vindexes[j]*2,-val[i][j]);
                SetValA(Vindexes[i]*2,Vindexes[j]*2+1,val[i][j]);
            }
    }

    ///add the cotangent stiffness matrix of a face, for both U and V
    void AddLaplacianTerm(const FaceType *f,const int Vindexes[3])
    {
        CoordType e[3];
        for (int k=0;k<3;k++)
            e[k]=f->cP((k+2)%3)-f->cP((k+1)%3);
        const ScalarType areaT=DoubleArea(*f)/2.0;
        for (int i=0;i<3;i++)
            for (int j=0;j<3;j++)
            {
                ScalarType val=(e[i]*e[j])/(4.0*areaT);
                SetValA(Vindexes[i]*2,Vindexes[j]*2,val);
                SetValA(Vindexes[i]*2+1,Vindexes[j]*2+1,val);
            }
    }

    ///add the right hand side of a face in the field-aligned mode: the divergence of PD1
    ///for U and of PD2 for V, area * K . grad(phi_i), with grad(phi_i) = (n ^ e_i) / (2 area)
    void AddFieldTerm(const FaceType *f,const int Vindexes[3],ScalarType vector_field_scale)
    {
        CoordType fNorm=f->cN();
        fNorm.Normalize();
        CoordType K1 = CoordType::Construct(f->cPD1());
        CoordType K2 = CoordType::Construct(f->cPD2());
        K1.Normalize();
        K2.Normalize();
        vector_field_scale/=2;
        for (int i=0;i<3;i++)
        {
            const CoordType neg_t = fNorm ^ (f->cP((i+2)%3) - f->cP((i+1)%3));
            b[Vindexes[i]*2]   += (K1*vector_field_scale) * neg_t;
            b[Vindexes[i]*2+1] += (K2*vector_field_scale) * neg_t;
        }
    }

    ///fixed vertices are a stiff quadratic penalty towards their current T()
    void AddFixedTerms()
    {
        const ScalarType penalization=1000000;
        for (VertexType *v : to_fix)
            for (int c=0;c<2;c++)
            {
                const int index=VertexIndex(v)*2+c;
                SetValA(index,index,penalization);
                b[index]+=v->T().P()[c]*penalization;
            }
    }

    ///each point constraint adds weight*(sum_k bary_k X_k - value)^2, with X the chosen
    ///coordinate of the face's corners
    void AddPointTerms()
    {
        for (const PointConstraint &pc : pointConstraints)
            for (int i=0;i<3;i++)
            {
                const int row=VertexIndex(pc.f->cV(i))*2+pc.coord;
                b[row]+=pc.weight*pc.bary[i]*pc.value;
                for (int j=0;j<3;j++)
                    SetValA(row,VertexIndex(pc.f->cV(j))*2+pc.coord,pc.weight*pc.bary[i]*pc.bary[j]);
            }
    }

    void BuildSystem(ScalarType vfscale)
    {
        VertexToInd.assign(mesh.vert.size(),-1);
        IndToVertex.clear();
        ForEachVertex(mesh, [&](VertexType &v){
            VertexToInd[tri::Index(mesh,&v)]=int(IndToVertex.size());
            IndToVertex.push_back(&v);
        });
        const int total_size=int(IndToVertex.size())*2;
        A.clear();
        A.reserve(size_t(mesh.fn)*36+to_fix.size()*2+pointConstraints.size()*9);
        b=Eigen::VectorXd::Zero(total_size);
        ForEachFace(mesh, [&](FaceType &f){
            const int Vindexes[3]={VertexIndex(f.V(0)),VertexIndex(f.V(1)),VertexIndex(f.V(2))};
            AddLaplacianTerm(&f,Vindexes);
            if (use_direction_field)
                AddFieldTerm(&f,Vindexes,vfscale);
            else
                AddAreaTerm(Vindexes);
        });
        AddFixedTerms();
        AddPointTerms();
    }

    bool Solve()
    {
        Eigen::SparseMatrix<double> As(b.size(),b.size());
        As.setFromTriplets(A.begin(),A.end());
        A.clear();
        A.shrink_to_fit();
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver(As);
        if (solver.info()!=Eigen::Success) return false;
        x = solver.solve(b);
        return (solver.info()==Eigen::Success);
    }

    ///map the solution back to the vertices and then to the wedges; fixed vertices keep
    ///their exact value if correct_fixed
    void MapCoords()
    {
        for (size_t i=0;i<IndToVertex.size();i++)
            IndToVertex[i]->T().P()=Point2<ScalarType>(x[i*2],x[i*2+1]);
        if (correct_fixed)
            for (size_t i=0;i<to_fix.size();i++)
                to_fix[i]->T().P()=fixedUV[i];
    }

    void CopyToWedges()
    {
        ForEachFace(mesh, [](FaceType &f){
            for (int j=0;j<3;j++)
            {
                f.WT(j).P()=f.V(j)->T().P();
                f.WT(j).N()=0; // a texture id, or the face reads as unmapped (see TexCoord2)
            }
        });
    }

    std::vector<Point2<ScalarType>> fixedUV;

public:

    ///return true if the mesh is a disk: manifold, genus 0, with at least one boundary.
    ///It updates the face-face adjacency.
    bool IsFeasible()
    {
        tri::UpdateTopology<MeshType>::FaceFace(mesh);
        if (tri::Clean<MeshType>::CountNonManifoldEdgeFF(mesh)!=0) return false;
        if (tri::Clean<MeshType>::CountNonManifoldVertexFF(mesh)!=0) return false;
        if (tri::Clean<MeshType>::CountHoles(mesh)==0) return false;
        return tri::Clean<MeshType>::MeshGenus(mesh)==0;
    }

    ///fix a vertex at the given texture coordinate
    void FixVertex(VertexType *v,const Point2<ScalarType> &uv)
    {
        v->T().P()=uv;
        to_fix.push_back(v);
    }

    ///fix the border vertices at their current T(); needs up to date border flags
    void SetBorderAsFixed()
    {
        ForEachVertex(mesh, [&](VertexType &v){ if(v.IsB()) to_fix.push_back(&v); });
    }

    ///fix the selected vertices at their current T()
    void SetSelectedAsFixed()
    {
        ForEachVertex(mesh, [&](VertexType &v){ if(v.IsS()) to_fix.push_back(&v); });
    }

    ///add a soft constraint pulling coordinate coord (0 for U, 1 for V) of the point of face f
    ///with barycentric coordinates bary towards value, with the given weight; the weight is
    ///relative to the Dirichlet energy, which is scale invariant
    void AddPointConstraint(FaceType *f,const CoordType &bary,int coord,ScalarType value,ScalarType weight)
    {
        if (coord<0 || coord>1) throw std::invalid_argument("PoissonSolver: coord must be 0 (U) or 1 (V)");
        pointConstraints.push_back({f,bary,coord,value,weight});
    }

    ///fix the vertices the solve needs when none are fixed yet: one at the origin in the
    ///field-aligned mode; otherwise the two extremes of the longest side of the bounding box,
    ///at (0,0) and (1,1)
    void FixDefaultVertices()
    {
        if (!to_fix.empty()) return;
        if (use_direction_field)
        {
            for (VertexType &v : mesh.vert)
                if (!v.IsD()) { FixVertex(&v,Point2<ScalarType>(0,0)); return; }
            return;
        }
        tri::UpdateBounding<MeshType>::Box(mesh);
        const int bestAxis = mesh.bbox.MaxDim();
        VertexType *v0=nullptr, *v1=nullptr;
        ForEachVertex(mesh, [&](VertexType &v){
            if(!v0 || v.P()[bestAxis] < v0->P()[bestAxis]) v0 = &v;
            if(!v1 || v.P()[bestAxis] > v1->P()[bestAxis]) v1 = &v;
        });
        if (v0==v1) throw std::invalid_argument("PoissonSolver: the mesh has no extent");
        FixVertex(v0,Point2<ScalarType>(0,0));
        FixVertex(v1,Point2<ScalarType>(1,1));
    }

    ///choose the mode and clear the fixed vertices and point constraints
    void Init(bool _use_direction_field=false,
              bool _correct_fixed=true)
    {
        use_direction_field=_use_direction_field;
        correct_fixed=_correct_fixed;
        to_fix.clear();
        pointConstraints.clear();
    }

    ///solve the system; it returns false if the matrix is singular. If solve_global_fold, a
    ///result with most faces folded is mirrored. Throws if too few vertices are fixed.
    bool SolvePoisson(CallBackPos *cb=nullptr,
                      ScalarType fieldScale=1.0,
                      bool solve_global_fold=true)
    {
        RequirePerVertexTexCoord(mesh);
        RequirePerFaceWedgeTexCoord(mesh);
        if (use_direction_field) RequirePerFaceCurvatureDir(mesh);
        std::sort(to_fix.begin(),to_fix.end());
        to_fix.erase(std::unique(to_fix.begin(),to_fix.end()),to_fix.end());
        if (to_fix.size() < (use_direction_field ? 1u : 2u))
            throw std::invalid_argument(use_direction_field ? "PoissonSolver: the field-aligned mode needs a fixed vertex"
                                                            : "PoissonSolver: the conformal mode needs two fixed vertices");
        fixedUV.clear();
        for (VertexType *v : to_fix) fixedUV.push_back(v->T().P());

        if (cb) cb(0,"Building the system");
        BuildSystem(fieldScale);
        if (cb) cb(30,"Solving the system");
        if (!Solve()) return false;
        if (cb) cb(90,"Assigning texture coordinates");
        MapCoords();

        if (solve_global_fold && tri::Distortion<MeshType,false>::GloballyUnFolded(mesh))
            tri::UV_Utils<MeshType>::GloballyMirrorX(mesh);
        CopyToWedges();
        return true;
    }

    PoissonSolver(MeshType &_mesh):mesh(_mesh),use_direction_field(false),correct_fixed(true) {}

}; // end class
} // End Namespace tri
} // End Namespace vcg
#endif
